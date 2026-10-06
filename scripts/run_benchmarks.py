import os
import subprocess
import re
import csv
import matplotlib.pyplot as plt

# Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEQ_EXE = os.path.join(BASE_DIR, "bin", "amr_screen_seq.exe")
OMP_EXE = os.path.join(BASE_DIR, "bin", "amr_screen_omp.exe")
DATA_CSV = os.path.join(BASE_DIR, "data", "veterinary_amr_kmer_dataset.csv")
REF_FASTA = os.path.join(BASE_DIR, "data", "reference_genes.fasta")
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
    # Match Average Screening Time: 0.018184 seconds
    m = re.search(r"Average Screening Time:\s+([0-9\.]+)\s+seconds", output)
    if m:
        return float(m.group(1))
    # Match benchmark result line
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

def benchmark_threads():
    print("=" * 70)
    print(" 1. BENCHMARKING THREAD SCALING (1000 Genomes, 15 Repetitions each)")
    print("=" * 70)
    
    # Run sequential baseline
    print("Running Sequential Baseline...")
    seq_out = run_cmd([SEQ_EXE, "-i", DATA_CSV, "-r", REF_FASTA, "-b", "15"])
    t_seq = parse_time(seq_out)
    tp_seq = parse_throughput(seq_out)
    kmer_seq = parse_kmer_lookups(seq_out)
    print(f"Sequential Time: {t_seq*1000:.3f} ms | Throughput: {tp_seq:.0f} genomes/s\n")
    
    threads_list = [1, 2, 4, 6, 8, 12, 16, 18]
    thread_results = []
    
    for th in threads_list:
        print(f"Running OpenMP with {th:2d} threads...")
        omp_out = run_cmd([OMP_EXE, "-i", DATA_CSV, "-r", REF_FASTA, "-threads", str(th), "-b", "15"])
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
        
    # Save CSV
    csv_file = os.path.join(RESULTS_DIR, "benchmark_thread_scaling.csv")
    with open(csv_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            'threads', 'time_sec', 'time_ms', 'speedup', 'efficiency_pct', 'throughput_gps', 'kmer_lookups_m_per_s'
        ])
        writer.writeheader()
        writer.writerows(thread_results)
    print(f"\nSaved thread scaling data to: {csv_file}")
    return t_seq, thread_results

def benchmark_dataset_size():
    print("\n" + "=" * 70)
    print(" 2. BENCHMARKING DATASET SIZE SCALING")
    print("=" * 70)
    
    sizes = [100, 250, 500, 750, 1000]
    config_threads = [1, 4, 8, 16]
    dataset_results = []
    
    for n in sizes:
        print(f"\nTesting dataset size N = {n} genomes...")
        # Sequential
        seq_out = run_cmd([SEQ_EXE, "-i", DATA_CSV, "-r", REF_FASTA, "-n", str(n), "-b", "10"])
        t_seq = parse_time(seq_out)
        
        row = {'genomes': n, 'sequential_ms': t_seq * 1000.0}
        
        # Multi-threaded OpenMP
        for th in config_threads:
            omp_out = run_cmd([OMP_EXE, "-i", DATA_CSV, "-r", REF_FASTA, "-n", str(n), "-threads", str(th), "-b", "10"])
            t_omp = parse_time(omp_out)
            row[f'omp_{th}th_ms'] = t_omp * 1000.0
            row[f'speedup_{th}th'] = t_seq / t_omp if t_omp > 0 else 0.0
            print(f"  N={n:4d} | Seq: {t_seq*1000:6.2f} ms | OpenMP ({th} th): {t_omp*1000:6.2f} ms (Speedup: {t_seq/t_omp:.2f}x)")
            
        dataset_results.append(row)
        
    # Save CSV
    csv_file = os.path.join(RESULTS_DIR, "benchmark_dataset_scaling.csv")
    fieldnames = ['genomes', 'sequential_ms'] + [f'omp_{th}th_ms' for th in config_threads] + [f'speedup_{th}th' for th in config_threads]
    with open(csv_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(dataset_results)
    print(f"\nSaved dataset scaling data to: {csv_file}")
    return dataset_results

def plot_performance_graphs(t_seq, thread_results, dataset_results):
    print("\n" + "=" * 70)
    print(" 3. GENERATING PERFORMANCE VISUALIZATION GRAPHS")
    print("=" * 70)
    
    threads = [r['threads'] for r in thread_results]
    time_ms = [r['time_ms'] for r in thread_results]
    speedups = [r['speedup'] for r in thread_results]
    efficiencies = [r['efficiency_pct'] for r in thread_results]
    throughputs = [r['throughput_gps'] for r in thread_results]
    
    # -------------------------------------------------------------
    # Plot 1: Execution Time vs Number of Threads
    # -------------------------------------------------------------
    plt.figure(figsize=(8, 5), dpi=300)
    plt.plot(threads, time_ms, marker='o', color='#1f77b4', linewidth=2.5, markersize=8, label='OpenMP Execution Time')
    plt.axhline(y=t_seq*1000.0, color='#d62728', linestyle='--', linewidth=2, label=f'Sequential Baseline ({t_seq*1000.0:.2f} ms)')
    plt.title('Execution Time vs Number of Threads\n(1,000 Bacterial Genomes, 5 Mbp total)', fontsize=13, fontweight='bold', pad=12)
    plt.xlabel('Number of OpenMP Threads', fontsize=11, fontweight='semibold')
    plt.ylabel('Execution Time (milliseconds)', fontsize=11, fontweight='semibold')
    plt.xticks(threads)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(frameon=True, facecolor='white', framealpha=0.9)
    for x, y in zip(threads, time_ms):
        plt.annotate(f'{y:.2f} ms', (x, y), textcoords="offset points", xytext=(0, 8), ha='center', fontsize=9, fontweight='bold')
    plt.tight_layout()
    chart1_path = os.path.join(RESULTS_DIR, "chart1_execution_time_vs_threads.png")
    plt.savefig(chart1_path)
    plt.close()
    print("Saved:", chart1_path)

    # -------------------------------------------------------------
    # Plot 2: Speedup vs Number of Threads
    # -------------------------------------------------------------
    plt.figure(figsize=(8, 5), dpi=300)
    plt.plot(threads, speedups, marker='s', color='#2ca02c', linewidth=2.5, markersize=8, label='Measured Speedup (OpenMP)')
    plt.plot(threads, threads, linestyle='--', color='#7f7f7f', linewidth=2, label='Ideal Linear Speedup (S = p)')
    plt.title('Parallel Speedup vs Number of Threads\n(Strong Scaling on 1,000 Genomes)', fontsize=13, fontweight='bold', pad=12)
    plt.xlabel('Number of OpenMP Threads (p)', fontsize=11, fontweight='semibold')
    plt.ylabel('Speedup (S = T_seq / T_par)', fontsize=11, fontweight='semibold')
    plt.xticks(threads)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(frameon=True, facecolor='white', framealpha=0.9)
    for x, y in zip(threads, speedups):
        plt.annotate(f'{y:.2f}x', (x, y), textcoords="offset points", xytext=(0, 8), ha='center', fontsize=9, fontweight='bold')
    plt.tight_layout()
    chart2_path = os.path.join(RESULTS_DIR, "chart2_speedup_vs_threads.png")
    plt.savefig(chart2_path)
    plt.close()
    print("Saved:", chart2_path)

    # -------------------------------------------------------------
    # Plot 3: Parallel Efficiency (%) vs Number of Threads
    # -------------------------------------------------------------
    plt.figure(figsize=(8, 5), dpi=300)
    plt.plot(threads, efficiencies, marker='^', color='#ff7f0e', linewidth=2.5, markersize=8, label='Parallel Efficiency')
    plt.axhline(y=100.0, color='#7f7f7f', linestyle='--', linewidth=2, label='Ideal 100% Efficiency')
    plt.title('Parallel Efficiency vs Number of Threads\n(E = Speedup / Threads * 100%)', fontsize=13, fontweight='bold', pad=12)
    plt.xlabel('Number of OpenMP Threads (p)', fontsize=11, fontweight='semibold')
    plt.ylabel('Parallel Efficiency (%)', fontsize=11, fontweight='semibold')
    plt.xticks(threads)
    plt.ylim(0, 115)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(frameon=True, facecolor='white', framealpha=0.9)
    for x, y in zip(threads, efficiencies):
        plt.annotate(f'{y:.1f}%', (x, y), textcoords="offset points", xytext=(0, 8), ha='center', fontsize=9, fontweight='bold')
    plt.tight_layout()
    chart3_path = os.path.join(RESULTS_DIR, "chart3_parallel_efficiency_vs_threads.png")
    plt.savefig(chart3_path)
    plt.close()
    print("Saved:", chart3_path)

    # -------------------------------------------------------------
    # Plot 4: Dataset Size Scaling (Execution Time vs Number of Genomes)
    # -------------------------------------------------------------
    plt.figure(figsize=(8, 5), dpi=300)
    genomes_n = [r['genomes'] for r in dataset_results]
    seq_time = [r['sequential_ms'] for r in dataset_results]
    omp4_time = [r['omp_4th_ms'] for r in dataset_results]
    omp8_time = [r['omp_8th_ms'] for r in dataset_results]
    omp16_time = [r['omp_16th_ms'] for r in dataset_results]
    
    plt.plot(genomes_n, seq_time, marker='o', linewidth=2.2, color='#d62728', label='Sequential (1 Core)')
    plt.plot(genomes_n, omp4_time, marker='s', linewidth=2.2, color='#1f77b4', label='OpenMP (4 Threads)')
    plt.plot(genomes_n, omp8_time, marker='^', linewidth=2.2, color='#2ca02c', label='OpenMP (8 Threads)')
    plt.plot(genomes_n, omp16_time, marker='D', linewidth=2.2, color='#9467bd', label='OpenMP (16 Threads)')
    
    plt.title('Dataset Size Scaling: Execution Time vs Number of Genomes', fontsize=13, fontweight='bold', pad=12)
    plt.xlabel('Dataset Size (Number of Genomes)', fontsize=11, fontweight='semibold')
    plt.ylabel('Execution Time (milliseconds)', fontsize=11, fontweight='semibold')
    plt.xticks(genomes_n)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(frameon=True, facecolor='white', framealpha=0.9)
    plt.tight_layout()
    chart4_path = os.path.join(RESULTS_DIR, "chart4_dataset_scaling.png")
    plt.savefig(chart4_path)
    plt.close()
    print("Saved:", chart4_path)

    # -------------------------------------------------------------
    # Plot 5: Screening Throughput Comparison (Genomes per Second)
    # -------------------------------------------------------------
    plt.figure(figsize=(8, 5), dpi=300)
    bars = plt.bar([str(t) for t in threads], throughputs, color='#3b528b', width=0.55, edgecolor='black', linewidth=0.8)
    plt.title('Screening Throughput vs Thread Count\n(Genomes Processed Per Second)', fontsize=13, fontweight='bold', pad=12)
    plt.xlabel('Number of OpenMP Threads', fontsize=11, fontweight='semibold')
    plt.ylabel('Throughput (Genomes / Second)', fontsize=11, fontweight='semibold')
    plt.grid(axis='y', linestyle=':', alpha=0.6)
    for bar in bars:
        height = bar.get_height()
        plt.annotate(f'{height:,.0f}',
                     xy=(bar.get_x() + bar.get_width() / 2, height),
                     xytext=(0, 4), textcoords="offset points",
                     ha='center', va='bottom', fontsize=8.5, fontweight='bold')
    plt.tight_layout()
    chart5_path = os.path.join(RESULTS_DIR, "chart5_screening_throughput.png")
    plt.savefig(chart5_path)
    plt.close()
    print("Saved:", chart5_path)

if __name__ == "__main__":
    t_seq, thread_results = benchmark_threads()
    dataset_results = benchmark_dataset_size()
    plot_performance_graphs(t_seq, thread_results, dataset_results)
    print("\nAll benchmarks and graphs successfully completed!")
