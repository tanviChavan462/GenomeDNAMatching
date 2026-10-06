import os
import subprocess
import re
import csv

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN_DIR = os.path.join(BASE_DIR, "bin")

# Look for kmer_sequential or amr_screen_seq
SEQ_EXE = os.path.join(BIN_DIR, "kmer_sequential.exe")
if not os.path.exists(SEQ_EXE):
    SEQ_EXE = os.path.join(BIN_DIR, "amr_screen_seq.exe")

OMP_EXE = os.path.join(BIN_DIR, "kmer_openmp.exe")
if not os.path.exists(OMP_EXE):
    OMP_EXE = os.path.join(BIN_DIR, "amr_screen_omp.exe")

DATA_DIR = os.path.join(BASE_DIR, "data")
ORIGINAL_CSV = os.path.join(DATA_DIR, "veterinary_amr_kmer_dataset.csv")
REF_FASTA = os.path.join(DATA_DIR, "reference_genes.fasta")
SCALED_DIR = os.path.join(DATA_DIR, "scaled")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

os.makedirs(RESULTS_DIR, exist_ok=True)

def run_cmd(cmd):
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running command: {' '.join(cmd)}")
        print("STDERR:", result.stderr)
        return None
    return result.stdout

def parse_time(output):
    m = re.search(r"Average Screening Time:\s+([0-9\.]+)\s+seconds", output)
    if m:
        return float(m.group(1))
    m = re.search(r"time_sec=([0-9\.]+)", output)
    if m:
        return float(m.group(1))
    return None

def parse_throughput(output):
    m = re.search(r"Throughput:\s+([0-9\.]+)\s+genomes/second", output)
    if m:
        return float(m.group(1))
    return None

def parse_kmer_lookups(output):
    m = re.search(r"K-mer Lookups / Sec:\s+([0-9\.]+)\s+million", output)
    if m:
        return float(m.group(1))
    return None

def benchmark_threads(benchmark_runs=15):
    print("=" * 70)
    print(f" 1. BENCHMARKING THREAD SCALING (1000 Genomes, {benchmark_runs} Repetitions each)")
    print("=" * 70)
    
    # Run sequential baseline
    print(f"Running Sequential Baseline using {os.path.basename(SEQ_EXE)}...")
    seq_out = run_cmd([SEQ_EXE, "-i", ORIGINAL_CSV, "-r", REF_FASTA, "-b", str(benchmark_runs)])
    t_seq = parse_time(seq_out)
    tp_seq = parse_throughput(seq_out)
    kmer_seq = parse_kmer_lookups(seq_out)
    print(f"Sequential Time: {t_seq*1000:.3f} ms | Throughput: {tp_seq:.0f} genomes/s\n")
    
    threads_list = [1, 2, 4, 6, 8, 12, 16, 18]
    thread_results = []
    
    for th in threads_list:
        print(f"Running OpenMP with {th:2d} threads...")
        omp_out = run_cmd([OMP_EXE, "-i", ORIGINAL_CSV, "-r", REF_FASTA, "-threads", str(th), "-b", str(benchmark_runs)])
        t_omp = parse_time(omp_out)
        tp_omp = parse_throughput(omp_out)
        kmer_omp = parse_kmer_lookups(omp_out)
        
        speedup = t_seq / t_omp if t_omp > 0 else 0.0
        efficiency = (speedup / th) * 100.0 if th > 0 else 0.0
        
        thread_results.append({
            'threads': th,
            'time_sec': t_omp,
            'time_ms': t_omp * 1000.0,
            'speedup': speedup,
            'efficiency_pct': efficiency,
            'throughput_gps': tp_omp,
            'kmer_lookups_m_per_s': kmer_omp
        })
        print(f"  -> Time: {t_omp*1000:6.3f} ms | Speedup: {speedup:5.2f}x | Efficiency: {efficiency:5.1f}% | Throughput: {tp_omp:8.0f} genomes/s")
        
    csv_file = os.path.join(RESULTS_DIR, "benchmark_thread_scaling.csv")
    with open(csv_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            'threads', 'time_sec', 'time_ms', 'speedup', 'efficiency_pct', 'throughput_gps', 'kmer_lookups_m_per_s'
        ])
        writer.writeheader()
        writer.writerows(thread_results)
    print(f"\nSaved thread scaling data to: {csv_file}")
    return t_seq, thread_results

def benchmark_dataset_size(benchmark_runs=10):
    print("\n" + "=" * 70)
    print(f" 2. BENCHMARKING DATASET SIZE SCALING (data/scaled/ datasets)")
    print("=" * 70)
    
    target_sizes = [100, 200, 400, 800, 1000]
    config_threads = [1, 4, 8, 16]
    dataset_results = []
    
    for n in target_sizes:
        scaled_csv = os.path.join(SCALED_DIR, f"amr_{n}.csv")
        # Fallback to original with -n if scaled_csv doesn't exist
        csv_to_use = scaled_csv if os.path.exists(scaled_csv) else ORIGINAL_CSV
        n_arg = ["-n", str(n)] if csv_to_use == ORIGINAL_CSV else []
        
        print(f"\nTesting dataset size N = {n} genomes ({os.path.basename(csv_to_use)})...")
        seq_cmd = [SEQ_EXE, "-i", csv_to_use, "-r", REF_FASTA, "-b", str(benchmark_runs)] + n_arg
        seq_out = run_cmd(seq_cmd)
        t_seq = parse_time(seq_out)
        
        row = {'genomes': n, 'sequential_ms': t_seq * 1000.0}
        
        for th in config_threads:
            omp_cmd = [OMP_EXE, "-i", csv_to_use, "-r", REF_FASTA, "-threads", str(th), "-b", str(benchmark_runs)] + n_arg
            omp_out = run_cmd(omp_cmd)
            t_omp = parse_time(omp_out)
            row[f'omp_{th}th_ms'] = t_omp * 1000.0
            row[f'speedup_{th}th'] = t_seq / t_omp if t_omp > 0 else 0.0
            print(f"  N={n:4d} | Seq: {t_seq*1000:6.2f} ms | OpenMP ({th:2d} th): {t_omp*1000:6.2f} ms (Speedup: {t_seq/t_omp:.2f}x)")
            
        dataset_results.append(row)
        
    csv_file = os.path.join(RESULTS_DIR, "benchmark_dataset_scaling.csv")
    fieldnames = ['genomes', 'sequential_ms'] + [f'omp_{th}th_ms' for th in config_threads] + [f'speedup_{th}th' for th in config_threads]
    with open(csv_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(dataset_results)
    print(f"\nSaved dataset scaling data to: {csv_file}")
    return dataset_results

if __name__ == "__main__":
    benchmark_threads(benchmark_runs=10)
    benchmark_dataset_size(benchmark_runs=10)
    print("\nBenchmark experiment finished! Run `python scripts/plot_results.py` to plot charts.")
