# LLM Usage & Technical Audit Log

This document provides a transparent record of AI-assisted engineering, technical auditing, and algorithmic debugging conducted during the development of this project.

---

## 1. Project Inception & Dataset Audit
- **Objective:** Build a sequential and OpenMP parallel screening engine for antimicrobial resistance genes in synthetic veterinary genomes.
- **Assistance Provided:**
  - Programmatically inspected `data/veterinary_amr_kmer_dataset.csv` without opening the entire 5 MB file in text editors.
  - Determined column schema (`sample_id`, `species`, `animal_source`, `amr_status`, `amr_gene`, `sequence`, `k`).
  - Identified sequence length ($5,000\text{ bp}$ fixed) and confirmed that reference gene sequences were embedded inside resistant genomes rather than provided as a separate file.
  - Reconstructed reference sequences for `blaTEM` ($224\text{ bp}$), `tetA` ($154\text{ bp}$), and `sul1` ($119\text{ bp}$) and saved them to `data/reference_genes.fasta`.

---

## 2. Key Technical Challenges & Algorithmic Debugging

### Challenge 1: The `sul1` Repetitive $k$-mer Bug (FN = 250)
- **Problem:** During early sequential testing, `blaTEM` and `tetA` achieved 100% detection, but `sul1` yielded 250 False Negatives (75% overall accuracy).
- **Diagnosis:** The `sul1` sequence ($119\text{ bp}$) contains internal tandem repeats. While there are $119 - 21 + 1 = 99$ sliding window positions, only $69$ are *distinct* $21$-mers. When computing coverage against the window count ($69 / 99 = 69.7\%$), it fell below the $80\%$ detection threshold.
- **Solution:** Modified `src/amr_common.h` to calculate coverage against the count of *unique* distinct $k$-mers in the reference index ($69 / 69 = 100\%$). Sensitivity rose to **$100.00\%$** across all genes.

### Challenge 2: Relative Path Resolution across Working Directories
- **Problem:** Executing binaries from `src/` failed to locate `data/reference_genes.fasta` because the relative path expected the root directory.
- **Solution:** Implemented `resolve_input_path()` and `resolve_output_path()` in `src/amr_common.h` to automatically test `./data/` and fallback to `../data/`, making execution robust from any directory.

---

## 3. High-Performance Audit: Resolving the 38.42x Speedup Anomaly
- **Observation:** A user test reported sequential time $= 192.10\text{ ms}$ and OpenMP 8-thread time $= 5.00\text{ ms}$, claiming a $38.42\times$ speedup on 8 cores ($480\%$ efficiency).
- **Audit Findings:**
  1. Both sequential and parallel implementations perform the identical workload: $1,000$ genomes, $5,000\text{ bp}$ sequences, and $4,980,000$ $k$-mer lookups.
  2. The timing blocks, data structures, and hash lookups are strictly identical.
  3. **Root Cause:** The sequential program had been compiled without optimization (`gcc amr_screen_seq.c -o amr_screen_seq`, defaulting to `-O0`), whereas OpenMP was compiled with `gcc -O3 -fopenmp`.
- **Resolution:** Re-benchmarked with fair, identical `-O3` optimization for both programs. The genuine 8-core speedup was measured at **$4.63\times$ to $5.10\times$** ($58\%\text{--}64\%$ efficiency), adhering to Amdahl's Law and realistic multicore scaling.

---

## 4. Verification Checklist
- [x] Sequential and parallel C codes compiled with identical flags (`-O3`).
- [x] Correctness verified with $100.00\%$ precision, recall, and specificity ($TP=750, TN=250, FP=0, FN=0$).
- [x] Sequential and parallel output CSVs verified bit-for-bit identical.
- [x] Thread scaling evaluated across 1, 2, 4, 6, 8, 12, 16, 18 threads.
- [x] Dataset size scaling evaluated across 100, 200, 400, 800, 1,000 genomes.
