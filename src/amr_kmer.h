#ifndef AMR_KMER_H
#define AMR_KMER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <ctype.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define DEFAULT_K 21
#define MAX_GENES 16
#define MAX_GENE_NAME 64
#define MAX_SEQ_LEN 16384
#define HASH_TABLE_CAPACITY 8192 // Power of 2, plenty of headroom for ~400 reference kmers

// 2-bit nucleotide encoding: A=00, C=01, G=10, T=11
static inline int encode_base(char c) {
    switch (toupper((unsigned char)c)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1; // non-canonical / ambiguous
    }
}

static inline char decode_base(int b) {
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    return (b >= 0 && b < 4) ? bases[b] : 'N';
}

// 64-bit integer hash function (Thomas Wang)
static inline uint64_t hash_u64(uint64_t key) {
    key = (~key) + (key << 21);
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8);
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4);
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

// Hash Table Entry for reference k-mer lookup
typedef struct {
    uint64_t kmer;
    int gene_id;     // 1-based gene index (1 to num_genes)
    int kmer_idx;    // 0-based index of distinct k-mer within this gene
    bool occupied;
} HashEntry;

// Reference Gene representation
typedef struct {
    char name[MAX_GENE_NAME];
    char *sequence;
    int seq_len;
    int total_distinct_kmers;
} ReferenceGene;

// Global Reference Database
typedef struct {
    int k;
    int num_genes;
    ReferenceGene genes[MAX_GENES];
    HashEntry table[HASH_TABLE_CAPACITY];
} ReferenceDB;

// Bacterial Genome Record loaded from CSV
typedef struct {
    char sample_id[32];
    char species[64];
    char animal_source[32];
    char ground_truth_status[16]; // "resistant" or "susceptible"
    char ground_truth_gene[32];   // "none", "blaTEM", etc.
    char *sequence;
    int seq_len;
} GenomeRecord;

// Screening Result for a Genome
typedef struct {
    char sample_id[32];
    char predicted_status[16];
    char predicted_gene[32];
    double best_coverage;
    int detected_gene_id; // 0 if none
    int hits_per_gene[MAX_GENES];
    bool is_correct;
} ScreeningResult;

// Initialize Reference Database
static inline void init_ref_db(ReferenceDB *db, int k) {
    memset(db, 0, sizeof(ReferenceDB));
    db->k = k;
    db->num_genes = 0;
    for (int i = 0; i < HASH_TABLE_CAPACITY; i++) {
        db->table[i].occupied = false;
    }
}

// Insert distinct k-mer into Reference Hash Table
static inline bool insert_ref_kmer(ReferenceDB *db, uint64_t kmer, int gene_id, int kmer_idx) {
    uint64_t mask = HASH_TABLE_CAPACITY - 1;
    uint64_t idx = hash_u64(kmer) & mask;
    
    for (int step = 0; step < HASH_TABLE_CAPACITY; step++) {
        uint64_t curr = (idx + step) & mask;
        if (!db->table[curr].occupied) {
            db->table[curr].kmer = kmer;
            db->table[curr].gene_id = gene_id;
            db->table[curr].kmer_idx = kmer_idx;
            db->table[curr].occupied = true;
            return true;
        } else if (db->table[curr].kmer == kmer && db->table[curr].gene_id == gene_id) {
            // Already present for this gene
            return false;
        }
    }
    return false; // table full
}

// Lookup k-mer in Reference Hash Table
static inline bool lookup_ref_kmer(const ReferenceDB *db, uint64_t kmer, int *out_gene_id, int *out_kmer_idx) {
    uint64_t mask = HASH_TABLE_CAPACITY - 1;
    uint64_t idx = hash_u64(kmer) & mask;
    
    for (int step = 0; step < HASH_TABLE_CAPACITY; step++) {
        uint64_t curr = (idx + step) & mask;
        if (!db->table[curr].occupied) {
            return false;
        }
        if (db->table[curr].kmer == kmer) {
            if (out_gene_id) *out_gene_id = db->table[curr].gene_id;
            if (out_kmer_idx) *out_kmer_idx = db->table[curr].kmer_idx;
            return true;
        }
    }
    return false;
}

// Process and insert k-mers from a gene sequence
static inline void process_gene_kmers(ReferenceDB *db, int gene_idx, const char *seq, int seq_len) {
    int k = db->k;
    int gene_id = gene_idx + 1; // 1-based
    uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
    uint64_t rolling_kmer = 0;
    int valid_bases = 0;
    int distinct_count = 0;
    
    for (int i = 0; i < seq_len; i++) {
        int b = encode_base(seq[i]);
        if (b >= 0) {
            rolling_kmer = ((rolling_kmer << 2) | (uint64_t)b) & mask;
            valid_bases++;
            if (valid_bases >= k) {
                // Try inserting with next distinct index
                if (insert_ref_kmer(db, rolling_kmer, gene_id, distinct_count)) {
                    distinct_count++;
                }
            }
        } else {
            valid_bases = 0;
            rolling_kmer = 0;
        }
    }
    db->genes[gene_idx].total_distinct_kmers = distinct_count;
}

// Resolve relative path if run from src/ directory or project root
static inline const char* resolve_input_path(const char *path, char *buf, size_t buf_size) {
    FILE *fp = fopen(path, "r");
    if (fp) {
        fclose(fp);
        return path;
    }
    snprintf(buf, buf_size, "../%s", path);
    fp = fopen(buf, "r");
    if (fp) {
        fclose(fp);
        return buf;
    }
    return path;
}

static inline const char* resolve_output_path(const char *path, char *buf, size_t buf_size) {
    FILE *fp = fopen(path, "w");
    if (fp) {
        fclose(fp);
        return path;
    }
    snprintf(buf, buf_size, "../%s", path);
    fp = fopen(buf, "w");
    if (fp) {
        fclose(fp);
        return buf;
    }
    return path;
}

// Load Reference Genes from FASTA file
static inline bool load_reference_genes(const char *fasta_path, ReferenceDB *db, int k) {
    init_ref_db(db, k);
    char path_buf[1024];
    const char *actual_path = resolve_input_path(fasta_path, path_buf, sizeof(path_buf));
    FILE *fp = fopen(actual_path, "r");
    if (!fp) {
        fprintf(stderr, "Error opening reference FASTA file: %s (tried %s)\n", fasta_path, actual_path);
        return false;
    }
    
    char line[1024];
    ReferenceGene *curr_gene = NULL;
    char seq_buffer[MAX_SEQ_LEN];
    seq_buffer[0] = '\0';
    
    while (fgets(line, sizeof(line), fp)) {
        size_t len = strlen(line);
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
            line[--len] = '\0';
        }
        if (len == 0) continue;
        
        if (line[0] == '>') {
            if (curr_gene && strlen(seq_buffer) > 0) {
                curr_gene->seq_len = (int)strlen(seq_buffer);
                curr_gene->sequence = strdup(seq_buffer);
                process_gene_kmers(db, db->num_genes - 1, curr_gene->sequence, curr_gene->seq_len);
            }
            
            if (db->num_genes >= MAX_GENES) {
                fprintf(stderr, "Exceeded MAX_GENES limit (%d)\n", MAX_GENES);
                break;
            }
            
            curr_gene = &db->genes[db->num_genes];
            db->num_genes++;
            
            char *token = strtok(line + 1, " \t\r\n");
            if (token) {
                strncpy(curr_gene->name, token, sizeof(curr_gene->name) - 1);
            }
            seq_buffer[0] = '\0';
        } else {
            strncat(seq_buffer, line, sizeof(seq_buffer) - strlen(seq_buffer) - 1);
        }
    }
    
    if (curr_gene && strlen(seq_buffer) > 0) {
        curr_gene->seq_len = (int)strlen(seq_buffer);
        curr_gene->sequence = strdup(seq_buffer);
        process_gene_kmers(db, db->num_genes - 1, curr_gene->sequence, curr_gene->seq_len);
    }
    
    fclose(fp);
    return true;
}

// Free reference DB allocated memory
static inline void free_ref_db(ReferenceDB *db) {
    for (int i = 0; i < db->num_genes; i++) {
        if (db->genes[i].sequence) free(db->genes[i].sequence);
    }
}

// Load Genome Records from CSV
static inline GenomeRecord* load_dataset_csv(const char *csv_path, int *out_count, int max_records) {
    char path_buf[1024];
    const char *actual_path = resolve_input_path(csv_path, path_buf, sizeof(path_buf));
    FILE *fp = fopen(actual_path, "r");
    if (!fp) {
        fprintf(stderr, "Error opening CSV: %s (tried %s)\n", csv_path, actual_path);
        return NULL;
    }
    
    int capacity = (max_records > 0) ? max_records : 1000;
    GenomeRecord *records = (GenomeRecord*)malloc(capacity * sizeof(GenomeRecord));
    if (!records) {
        fclose(fp);
        return NULL;
    }
    
    char line_buf[16384];
    
    // Read header line
    if (!fgets(line_buf, sizeof(line_buf), fp)) {
        free(records);
        fclose(fp);
        return NULL;
    }
    
    int count = 0;
    while (fgets(line_buf, sizeof(line_buf), fp)) {
        if (max_records > 0 && count >= max_records) {
            break;
        }
        
        size_t len = strlen(line_buf);
        while (len > 0 && (line_buf[len - 1] == '\n' || line_buf[len - 1] == '\r')) {
            line_buf[--len] = '\0';
        }
        if (len == 0) continue;
        
        char *tok;
        char *rest = line_buf;
        
        // sample_id
        tok = strtok_s(rest, ",", &rest);
        if (!tok) continue;
        strncpy(records[count].sample_id, tok, sizeof(records[count].sample_id) - 1);
        records[count].sample_id[sizeof(records[count].sample_id) - 1] = '\0';
        
        // species
        tok = strtok_s(rest, ",", &rest);
        if (!tok) continue;
        strncpy(records[count].species, tok, sizeof(records[count].species) - 1);
        records[count].species[sizeof(records[count].species) - 1] = '\0';
        
        // animal_source
        tok = strtok_s(rest, ",", &rest);
        if (!tok) continue;
        strncpy(records[count].animal_source, tok, sizeof(records[count].animal_source) - 1);
        records[count].animal_source[sizeof(records[count].animal_source) - 1] = '\0';
        
        // amr_status
        tok = strtok_s(rest, ",", &rest);
        if (!tok) continue;
        strncpy(records[count].ground_truth_status, tok, sizeof(records[count].ground_truth_status) - 1);
        records[count].ground_truth_status[sizeof(records[count].ground_truth_status) - 1] = '\0';
        
        // amr_gene
        tok = strtok_s(rest, ",", &rest);
        if (!tok) continue;
        strncpy(records[count].ground_truth_gene, tok, sizeof(records[count].ground_truth_gene) - 1);
        records[count].ground_truth_gene[sizeof(records[count].ground_truth_gene) - 1] = '\0';
        
        // sequence
        tok = strtok_s(rest, ",", &rest);
        if (!tok) continue;
        records[count].seq_len = (int)strlen(tok);
        records[count].sequence = strdup(tok);
        
        count++;
        if (count >= capacity && max_records <= 0) {
            capacity *= 2;
            GenomeRecord *tmp = (GenomeRecord*)realloc(records, capacity * sizeof(GenomeRecord));
            if (!tmp) {
                fprintf(stderr, "Out of memory resizing genome array\n");
                break;
            }
            records = tmp;
        }
    }
    
    fclose(fp);
    *out_count = count;
    return records;
}

// Free Genome dataset
static inline void free_dataset(GenomeRecord *records, int count) {
    if (!records) return;
    for (int i = 0; i < count; i++) {
        if (records[i].sequence) free(records[i].sequence);
    }
    free(records);
}

// Screen a single genome against the reference DB
static inline void screen_genome(const GenomeRecord *rec, const ReferenceDB *db, ScreeningResult *res, double threshold) {
    strncpy(res->sample_id, rec->sample_id, sizeof(res->sample_id) - 1);
    res->sample_id[sizeof(res->sample_id) - 1] = '\0';
    
    int k = db->k;
    int seq_len = rec->seq_len;
    
    uint8_t hit_kmers[MAX_GENES][1024];
    memset(hit_kmers, 0, sizeof(hit_kmers));
    
    for (int g = 0; g < MAX_GENES; g++) {
        res->hits_per_gene[g] = 0;
    }
    
    uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
    uint64_t rolling_kmer = 0;
    int valid_bases = 0;
    
    for (int i = 0; i < seq_len; i++) {
        int b = encode_base(rec->sequence[i]);
        if (b >= 0) {
            rolling_kmer = ((rolling_kmer << 2) | (uint64_t)b) & mask;
            valid_bases++;
            if (valid_bases >= k) {
                int gene_id = 0;
                int kmer_idx = 0;
                if (lookup_ref_kmer(db, rolling_kmer, &gene_id, &kmer_idx)) {
                    int g_idx = gene_id - 1;
                    if (g_idx >= 0 && g_idx < db->num_genes && kmer_idx >= 0 && kmer_idx < 1024) {
                        if (!hit_kmers[g_idx][kmer_idx]) {
                            hit_kmers[g_idx][kmer_idx] = 1;
                            res->hits_per_gene[g_idx]++;
                        }
                    }
                }
            }
        } else {
            valid_bases = 0;
            rolling_kmer = 0;
        }
    }
    
    double best_cov = 0.0;
    int best_gene_idx = -1;
    
    for (int g = 0; g < db->num_genes; g++) {
        int total = db->genes[g].total_distinct_kmers;
        if (total > 0) {
            double cov = (double)res->hits_per_gene[g] / (double)total;
            if (cov >= threshold && cov > best_cov) {
                best_cov = cov;
                best_gene_idx = g;
            }
        }
    }
    
    res->best_coverage = best_cov;
    if (best_gene_idx >= 0) {
        res->detected_gene_id = best_gene_idx + 1;
        strcpy(res->predicted_status, "resistant");
        strncpy(res->predicted_gene, db->genes[best_gene_idx].name, sizeof(res->predicted_gene) - 1);
    } else {
        res->detected_gene_id = 0;
        strcpy(res->predicted_status, "susceptible");
        strcpy(res->predicted_gene, "none");
    }
    
    bool status_match = (strcmp(res->predicted_status, rec->ground_truth_status) == 0);
    bool gene_match = (strcmp(res->predicted_gene, rec->ground_truth_gene) == 0);
    res->is_correct = (status_match && gene_match);
}

#endif // AMR_KMER_H
