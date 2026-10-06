# Team Member Contributions

**Project Title:** Parallel $k$-mer-Based Screening of Synthetic Veterinary Bacterial Genomes for Antimicrobial-Resistance Genes using OpenMP  
**Course / Laboratory:** High-Performance Computing & Computational Genomics  

---

## 👥 Responsibility Matrix

| Project Component | Primary Responsibility | Key Deliverables & Contributions |
| :--- | :--- | :--- |
| **Dataset Analysis & Reference DB** | Team Member 1 | • Programmatic inspection of `veterinary_amr_kmer_dataset.csv`<br>• Reverse-engineering embedded AMR cassettes (`blaTEM`, `tetA`, `sul1`)<br>• Generating reference database (`reference_genes.fasta`)<br>• Scaled dataset generation script (`scripts/generate_datasets.py`) |
| **Sequential C Algorithm** | Team Member 2 | • Designing 2-bit nucleotide encoding (`uint64_t`)<br>• Implementing $O(1)$ rolling $k$-mer extraction logic<br>• Open-addressing hash table with Thomas Wang 64-bit integer hash<br>• Implementing sequential baseline (`src/kmer_sequential.c`) |
| **Parallel OpenMP Architecture** | Team Member 3 | • Designing embarrassingly parallel genome distribution<br>• OpenMP multi-threading with `#pragma omp parallel for`<br>• Dynamic chunk scheduling evaluation<br>• Race-condition verification & thread-isolated result buffers (`src/kmer_openmp.c`) |
| **Benchmarking & Visualization** | Team Member 4 | • Writing automated benchmarking pipeline (`scripts/benchmark.py`)<br>• Plotting performance graphs with Matplotlib (`scripts/plot_results.py`)<br>• Thread scaling (1 to 18 cores) and dataset scaling ($N=100$ to $1000$)<br>• Calculating speedup, parallel efficiency, and throughput |
| **Technical Audit & Documentation** | All Team Members | • Auditing the $38.42\times$ apparent speedup vs. compiler optimization flags<br>• Debugging repetitive $k$-mers in `sul1` for $100\%$ sensitivity<br>• Preparing `README.md`, `Makefile`, and `LLM_USAGE_LOG.md`<br>• Final report preparation and viva presentation defense |

---

## 📝 Individual Work Breakdown (Fill In With Team Names)

### Member 1: [Name / Roll No.]
- Researched veterinary AMR pathogens (*E. coli*, *S. enterica*).
- Created dataset scaling scripts for weak-scaling experiments.
- Validated classification sensitivity across animal host species.

### Member 2: [Name / Roll No.]
- Implemented C data structures for $k$-mer representation.
- Designed hash table lookup mechanism with low collision factor.
- Authored sequential screening engine and validation checks.

### Member 3: [Name / Roll No.]
- Integrated OpenMP parallel directives and barrier synchronization.
- Measured cache behavior, NUMA effects, and workload distribution.
- Audited compiler flags (`-O3`, `-fopenmp`) and timing precision (`omp_get_wtime`).

### Member 4: [Name / Roll No.]
- Engineered automated Python benchmarking harness.
- Produced high-resolution visualization charts for reports and slides.
- Drafted Amdahl's Law analysis and parallel efficiency discussions.
