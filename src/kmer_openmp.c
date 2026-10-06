#include "amr_common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

void print_usage(const char *prog) {
    printf("Usage: %s [options]\n", prog);
    printf("Options:\n");
    printf("  -i <path>        Input dataset CSV (default: data/veterinary_amr_kmer_dataset.csv)\n");
    printf("  -r <path>        Reference genes FASTA (default: data/reference_genes.fasta)\n");
    printf("  -o <path>        Output predictions CSV (default: results/predictions_omp.csv)\n");
    printf("  -n <count>       Number of genomes to screen (default: all)\n");
    printf("  -k <value>       k-mer length (default: 21)\n");
    printf("  -t <threshold>   Coverage threshold for match (default: 0.80)\n");
    printf("  -threads <count> Number of OpenMP threads to use (default: max available)\n");
    printf("  -chunk <size>    OpenMP chunk size (default: 16)\n");
    printf("  -b <runs>        Number of benchmark repeat runs (default: 1)\n");
    printf("  -h               Show this help message\n");
}

int main(int argc, char *argv[]) {
    const char *csv_path = "data/veterinary_amr_kmer_dataset.csv";
    const char *fasta_path = "data/reference_genes.fasta";
    const char *output_csv = "results/predictions_omp.csv";
    int max_records = 0;
    int k = DEFAULT_K;
    double threshold = 0.80;
    int num_threads = 0;
    int chunk_size = 16;
    int benchmark_runs = 1;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) csv_path = argv[++i];
        else if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) fasta_path = argv[++i];
        else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) output_csv = argv[++i];
        else if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) max_records = atoi(argv[++i]);
        else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc) k = atoi(argv[++i]);
        else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) threshold = atof(argv[++i]);
        else if (strcmp(argv[i], "-threads") == 0 && i + 1 < argc) num_threads = atoi(argv[++i]);
        else if (strcmp(argv[i], "-chunk") == 0 && i + 1 < argc) chunk_size = atoi(argv[++i]);
        else if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) benchmark_runs = atoi(argv[++i]);
        else if (strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        }
    }

    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    } else {
        num_threads = omp_get_max_threads();
    }

    printf("============================================================\n");
    printf("  PARALLEL OpenMP AMR K-MER SCREENING PIPELINE              \n");
    printf("============================================================\n");
    printf("Dataset CSV:       %s\n", csv_path);
    printf("Reference FASTA:   %s\n", fasta_path);
    printf("Output CSV:        %s\n", output_csv);
    printf("k-mer Size:        %d\n", k);
    printf("Coverage Thresh:   %.2f\n", threshold);
    printf("OpenMP Threads:    %d (System Max: %d)\n", num_threads, omp_get_max_threads());
    printf("Dynamic Chunk:     %d\n", chunk_size);
    printf("Benchmark Runs:    %d\n", benchmark_runs);

    // 1. Load Reference Database
    ReferenceDB ref_db;
    double t_ref_start = omp_get_wtime();
    if (!load_reference_genes(fasta_path, &ref_db, k)) {
        fprintf(stderr, "Failed to load reference genes.\n");
        return 1;
    }
    double t_ref_end = omp_get_wtime();
    printf("[Init] Loaded %d reference genes in %.4f ms:\n", ref_db.num_genes, (t_ref_end - t_ref_start) * 1000.0);
    for (int g = 0; g < ref_db.num_genes; g++) {
        printf("       Gene %d: %-8s (len=%d bp, kmers=%d)\n", 
               g + 1, ref_db.genes[g].name, ref_db.genes[g].seq_len, ref_db.genes[g].total_distinct_kmers);
    }

    // 2. Load Dataset Genomes
    int num_genomes = 0;
    double t_data_start = omp_get_wtime();
    GenomeRecord *genomes = load_dataset_csv(csv_path, &num_genomes, max_records);
    double t_data_end = omp_get_wtime();
    if (!genomes || num_genomes == 0) {
        fprintf(stderr, "Failed to load genome dataset.\n");
        free_ref_db(&ref_db);
        return 1;
    }
    printf("[Init] Loaded %d genomes in %.4f ms (sequence length = %d bp)\n", 
           num_genomes, (t_data_end - t_data_start) * 1000.0, genomes[0].seq_len);

    ScreeningResult *results = (ScreeningResult*)malloc(num_genomes * sizeof(ScreeningResult));
    if (!results) {
        fprintf(stderr, "Failed to allocate memory for results.\n");
        free_dataset(genomes, num_genomes);
        free_ref_db(&ref_db);
        return 1;
    }

    // 3. Warm-up parallel run
    #pragma omp parallel for schedule(dynamic, chunk_size)
    for (int i = 0; i < num_genomes; i++) {
        screen_genome(&genomes[i], &ref_db, &results[i], threshold);
    }

    // 4. Timed Parallel Benchmark Screening
    printf("\n[Execution] Running parallel OpenMP screening (%d threads, %d benchmark runs)...\n", 
           num_threads, benchmark_runs);
    double total_screening_time = 0.0;
    double min_screening_time = 1e9;
    double max_screening_time = 0.0;

    for (int run = 0; run < benchmark_runs; run++) {
        double t_start = omp_get_wtime();

        #pragma omp parallel for schedule(dynamic, chunk_size)
        for (int i = 0; i < num_genomes; i++) {
            screen_genome(&genomes[i], &ref_db, &results[i], threshold);
        }

        double t_end = omp_get_wtime();
        double elapsed = t_end - t_start;
        total_screening_time += elapsed;
        if (elapsed < min_screening_time) min_screening_time = elapsed;
        if (elapsed > max_screening_time) max_screening_time = elapsed;
    }

    double avg_screening_time = total_screening_time / benchmark_runs;

    // 5. Verification & Correctness Statistics
    int tp = 0, tn = 0, fp = 0, fn = 0;
    int exact_gene_matches = 0;
    int total_resistant_true = 0;
    int total_susceptible_true = 0;

    for (int i = 0; i < num_genomes; i++) {
        bool true_res = (strcmp(genomes[i].ground_truth_status, "resistant") == 0);
        bool pred_res = (strcmp(results[i].predicted_status, "resistant") == 0);

        if (true_res) total_resistant_true++;
        else total_susceptible_true++;

        if (true_res && pred_res) {
            tp++;
            if (strcmp(results[i].predicted_gene, genomes[i].ground_truth_gene) == 0) {
                exact_gene_matches++;
            }
        } else if (!true_res && !pred_res) {
            tn++;
            if (strcmp(results[i].predicted_gene, "none") == 0) {
                exact_gene_matches++;
            }
        } else if (!true_res && pred_res) {
            fp++;
        } else if (true_res && !pred_res) {
            fn++;
        }
    }

    double accuracy = (double)(tp + tn) / num_genomes * 100.0;
    double sensitivity = (total_resistant_true > 0) ? (double)tp / total_resistant_true * 100.0 : 0.0;
    double specificity = (total_susceptible_true > 0) ? (double)tn / total_susceptible_true * 100.0 : 0.0;
    double gene_accuracy = (double)exact_gene_matches / num_genomes * 100.0;

    printf("\n============================================================\n");
    printf("  CORRECTNESS EVALUATION (Ground Truth Comparison)          \n");
    printf("============================================================\n");
    printf("Total Genomes Screened:      %d\n", num_genomes);
    printf("Ground Truth Resistant:      %d\n", total_resistant_true);
    printf("Ground Truth Susceptible:    %d\n", total_susceptible_true);
    printf("True Positives (TP):         %d\n", tp);
    printf("True Negatives (TN):         %d\n", tn);
    printf("False Positives (FP):        %d\n", fp);
    printf("False Negatives (FN):        %d\n", fn);
    printf("Status Classification Acc:   %.2f %%\n", accuracy);
    printf("Sensitivity (Recall):        %.2f %%\n", sensitivity);
    printf("Specificity:                 %.2f %%\n", specificity);
    printf("Exact Gene Match Accuracy:   %.2f %%\n", gene_accuracy);

    printf("\n============================================================\n");
    printf("  PARALLEL OpenMP EXECUTION TIMING                          \n");
    printf("============================================================\n");
    printf("Threads:                     %d\n", num_threads);
    printf("Average Screening Time:      %.6f seconds (%.2f ms)\n", avg_screening_time, avg_screening_time * 1000.0);
    printf("Min Screening Time:          %.6f seconds\n", min_screening_time);
    printf("Max Screening Time:          %.6f seconds\n", max_screening_time);
    printf("Throughput:                  %.1f genomes/second\n", num_genomes / avg_screening_time);
    printf("K-mer Lookups / Sec:         %.2f million k-mers/sec\n", 
           ((double)num_genomes * (genomes[0].seq_len - k + 1) / avg_screening_time) / 1e6);

    // 6. Write Predictions CSV
    char out_buf[1024];
    const char *actual_output_csv = resolve_output_path(output_csv, out_buf, sizeof(out_buf));
    FILE *out_fp = fopen(actual_output_csv, "w");
    if (out_fp) {
        fprintf(out_fp, "sample_id,species,animal_source,true_status,true_gene,pred_status,pred_gene,coverage,is_correct\n");
        for (int i = 0; i < num_genomes; i++) {
            fprintf(out_fp, "%s,%s,%s,%s,%s,%s,%s,%.4f,%d\n",
                    genomes[i].sample_id, genomes[i].species, genomes[i].animal_source,
                    genomes[i].ground_truth_status, genomes[i].ground_truth_gene,
                    results[i].predicted_status, results[i].predicted_gene,
                    results[i].best_coverage, results[i].is_correct ? 1 : 0);
        }
        fclose(out_fp);
        printf("\nSaved predictions to:        %s\n", actual_output_csv);
    }

    // Output machine-readable timing for benchmarking harness: [TIME_SEC: <value>]
    printf("\n[BENCHMARK_RESULT] threads=%d genomes=%d time_sec=%.6f\n", 
           num_threads, num_genomes, avg_screening_time);

    // Free resources
    free(results);
    free_dataset(genomes, num_genomes);
    free_ref_db(&ref_db);
    return 0;
}
