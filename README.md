# Parallel $k$-mer-Based Screening of Synthetic Veterinary Bacterial Genomes for Antimicrobial-Resistance Genes using OpenMP

A high-performance, alignment-free genomic screening system implemented in **C** and parallelized with **OpenMP** to rapidly identify antimicrobial resistance (AMR) genes in veterinary bacterial genomes.
Antimicrobial resistance (AMR) makes bacterial infections harder to treat and creates a need for 
efficient genome-based surveillance. A bacterial genome is a DNA sequence made of the bases 
A, C, G and T, and an AMR gene is a sequence associated with resistance to an antimicrobial 
drug. A genome can be screened by comparing it with known AMR reference genes to find 
evidence that a particular resistance gene is present. 
The problem statement for this project is to develop a computational system for processing large 
genomic datasets to identify useful patterns, similarities or variations, and to investigate how 
computational performance changes as the dataset size increases. In this project, the useful 
pattern is the occurrence of k-mers belonging to known AMR reference genes in synthetic 
veterinary bacterial genomes. 
A k-mer is a substring of length k taken from a DNA sequence. The project uses k = 21 and three 
reference genes: blaTEM, tetA and sul1. The proportion of a reference gene’s k-mers that are 
found in a genome is the coverage of that gene. A genome is classified as resistant to a gene 
when the coverage reaches the threshold of 0.80. 
The major computational challenge is that every genome must be compared with every reference 
gene, and this work is repeated for each genome in the dataset. As the number of genomes 
increases, the execution time also increases. Since the screening of one genome does not depend 
on the result of any other genome, parallel processing can be used to reduce the computation 
time.
---

## 📂 Repository Structure

```
parallel_kmer_amr_project/
│
├── data/
│   ├── veterinary_amr_kmer_dataset.csv   ← Original 1,000-genome dataset (5,000 bp each)
│   ├── reference_genes.fasta             ← Isolated reference AMR genes (blaTEM, tetA, sul1)
│   └── scaled/
│       ├── amr_100.csv                  ← Scaled subset: 100 genomes
│       ├── amr_200.csv                  ← Scaled subset: 200 genomes
│       ├── amr_400.csv                  ← Scaled subset: 400 genomes
│       ├── amr_800.csv                  ← Scaled subset: 800 genomes
│       └── amr_1000.csv                 ← Scaled subset: 1,000 genomes
│
├── src/
│   ├── amr_common.h                     ← Common 2-bit encoder, hash table, and screening logic
│   ├── kmer_sequential.c                ← Sequential C baseline implementation
│   └── kmer_openmp.c                    ← Parallel OpenMP C implementation
│
├── scripts/
│   ├── generate_datasets.py             ← Generates scaled datasets in data/scaled/
│   ├── benchmark.py                     ← Runs thread and dataset scaling experiments
│   └── plot_results.py                  ← Generates publication-quality performance graphs
│
├── results/                             ← Benchmarks, predictions, and charts
│   ├── benchmark_thread_scaling.csv
│   ├── benchmark_dataset_scaling.csv
│   ├── chart1_execution_time_vs_threads.png
│   ├── chart2_speedup_vs_threads.png
│   ├── chart3_parallel_efficiency_vs_threads.png
│   ├── chart4_dataset_scaling.png
│   └── chart5_screening_throughput.png
│
├── Makefile                             ← Build and run automation
├── README.md                            ← Project documentation and user guide
├── LLM_USAGE_LOG.md                     ← Transparency log of AI assistance and auditing
└── TEAM_CONTRIBUTIONS.md                ← Project team workload breakdown
```

---

## 🚀 Quick Start Guide

### 1. Prerequisites
- **C Compiler:** GCC (MinGW-w64 on Windows or Linux GCC) with OpenMP support (`-fopenmp`)
- **Python:** Python 3.8+ with `matplotlib` for generating graphs

```powershell
# Install Matplotlib (if not already installed)
python -m pip install matplotlib
```

---

### 2. Compilation

Compile both implementations using the `-O3` optimization flag:

```powershell
# Compile Sequential Version
gcc -O3 src/kmer_sequential.c -o bin/kmer_sequential.exe

# Compile Parallel OpenMP Version
gcc -O3 -fopenmp src/kmer_openmp.c -o bin/kmer_openmp.exe
```

*Alternatively, if using `make`:*
```powershell
make
```

---

### 3. Execution

#### Sequential Baseline
```powershell
.\bin\kmer_sequential.exe -b 5
```

#### Parallel OpenMP Execution
Test across different thread counts ($1, 2, 4, 8, 16$):
```powershell
# 1 Thread
.\bin\kmer_openmp.exe -threads 1 -b 5

# 2 Threads
.\bin\kmer_openmp.exe -threads 2 -b 5

# 4 Threads
.\bin\kmer_openmp.exe -threads 4 -b 5

# 8 Threads
.\bin\kmer_openmp.exe -threads 8 -b 5

# 16 Threads
.\bin\kmer_openmp.exe -threads 16 -b 5
```

---

### 4. Running Benchmarks & Generating Plots

Generate the scaled datasets and run the full benchmarking experiment:

```powershell
# Step 1: Generate scaled datasets (100, 200, 400, 800, 1000 genomes)
python scripts/generate_datasets.py

# Step 2: Run thread & dataset scaling experiments
python scripts/benchmark.py

# Step 3: Generate the performance graphs
python scripts/plot_results.py
```

All figures are saved to `results/`:
- `chart1_execution_time_vs_threads.png`: Execution time drop across thread counts
- `chart2_speedup_vs_threads.png`: Speedup curve vs. ideal linear speedup ($y=x$)
- `chart3_parallel_efficiency_vs_threads.png`: Parallel efficiency (%) across thread counts
- `chart4_dataset_scaling.png`: Execution time scaling with dataset size
- `chart5_screening_throughput.png`: Genomes processed per second bar chart

---

## 🔬 Algorithmic Overview

1. **2-Bit Nucleotide Encoding (`uint64_t`):**
   Nucleotides are mapped to 2-bit values ($A=00_2, C=01_2, G=10_2, T=11_2$). A $21$-mer requires $42\text{ bits}$, which fits into a single $64$-bit unsigned integer.
2. **$O(1)$ Rolling Window:**
   Sliding to the next $21$-mer requires only bit-shifts and a bitmask, avoiding string parsing:
   $$\text{kmer}_{i+1} = \left((\text{kmer}_i \ll 2) \mid \text{base}_{i+21}\right) \;\&\; \text{MASK}_{42}$$
3. **Read-Only Hash Table:**
   Reference $k$-mers are indexed in an open-addressing hash table with $\sim 5\%$ load factor. During screening, lookups are completely read-only, eliminating mutex locks and data races.
4. **Data-Parallel OpenMP Architecture:**
   Each genome is screened independently using `#pragma omp parallel for schedule(dynamic)`. Output is written directly into isolated array indices (`results[i]`), preventing false sharing.

---

## 📊 Summary of Benchmark Findings

- **Diagnostic Accuracy:** **$100.00\%$** Sensitivity, Specificity, and Gene Classification ($TP=750, TN=250, FP=0, FN=0$).
- **Bit-for-Bit Consistency:** Parallel outputs match sequential outputs identically.
- **Speedup:** Reaches **$5.15\times$ speedup** at 12 threads on strong scaling, and up to **$7.07\times$ speedup** across dataset sizes.
- **Throughput:** Processes up to **$205,000\text{ genomes/second}$** ($>1.3\text{ billion } k\text{-mer lookups/second}$).
