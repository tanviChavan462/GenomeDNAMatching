import os
import csv
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "results")

THREAD_CSV = os.path.join(RESULTS_DIR, "benchmark_thread_scaling.csv")
DATASET_CSV = os.path.join(RESULTS_DIR, "benchmark_dataset_scaling.csv")

def plot_all():
    if not os.path.exists(THREAD_CSV) or not os.path.exists(DATASET_CSV):
        print("Benchmark CSV files not found. Please run `python scripts/benchmark.py` first.")
        return

    # 1. Read Thread Scaling
    thread_results = []
    with open(THREAD_CSV, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            thread_results.append({
                'threads': int(row['threads']),
                'time_ms': float(row['time_ms']),
                'speedup': float(row['speedup']),
                'efficiency_pct': float(row['efficiency_pct']),
                'throughput_gps': float(row['throughput_gps']),
            })

    # 2. Read Dataset Scaling
    dataset_results = []
    with open(DATASET_CSV, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            dataset_results.append({
                'genomes': int(row['genomes']),
                'sequential_ms': float(row['sequential_ms']),
                'omp_4th_ms': float(row.get('omp_4th_ms', 0)),
                'omp_8th_ms': float(row.get('omp_8th_ms', 0)),
                'omp_16th_ms': float(row.get('omp_16th_ms', 0)),
            })

    threads = [r['threads'] for r in thread_results]
    time_ms = [r['time_ms'] for r in thread_results]
    speedups = [r['speedup'] for r in thread_results]
    efficiencies = [r['efficiency_pct'] for r in thread_results]
    throughputs = [r['throughput_gps'] for r in thread_results]

    # Chart 1: Execution Time vs Threads
    plt.figure(figsize=(8, 5), dpi=300)
    plt.plot(threads, time_ms, marker='o', color='#1f77b4', linewidth=2.5, markersize=8, label='OpenMP Execution Time')
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

    # Chart 2: Speedup vs Threads
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

    # Chart 3: Parallel Efficiency vs Threads
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

    # Chart 4: Dataset Size Scaling
    plt.figure(figsize=(8, 5), dpi=300)
    genomes_n = [r['genomes'] for r in dataset_results]
    seq_time = [r['sequential_ms'] for r in dataset_results]
    omp4_time = [r['omp_4th_ms'] for r in dataset_results]
    omp8_time = [r['omp_8th_ms'] for r in dataset_results]
    omp16_time = [r['omp_16th_ms'] for r in dataset_results]
    
    plt.plot(genomes_n, seq_time, marker='o', linewidth=2.2, color='#d62728', label='Sequential (1 Core)')
    if any(omp4_time):
        plt.plot(genomes_n, omp4_time, marker='s', linewidth=2.2, color='#1f77b4', label='OpenMP (4 Threads)')
    if any(omp8_time):
        plt.plot(genomes_n, omp8_time, marker='^', linewidth=2.2, color='#2ca02c', label='OpenMP (8 Threads)')
    if any(omp16_time):
        plt.plot(genomes_n, omp16_time, marker='D', linewidth=2.2, color='#9467bd', label='OpenMP (16 Threads)')
    
    plt.title('Dataset Size Scaling: Execution Time vs Number of Genomes\n(Scaled Datasets: amr_100.csv to amr_1000.csv)', fontsize=12, fontweight='bold', pad=12)
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

    # Chart 5: Throughput Bar Chart
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
    plot_all()
    print("All charts successfully generated in results/!")
