import base64, re, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(BASE, "results")
DASHBOARD = os.path.join(BASE, "dashboard.html")

charts = [
    ("chart1_execution_time_vs_threads.png",  "Time vs threads"),
    ("chart2_speedup_vs_threads.png",          "Speedup vs threads"),
    ("chart3_parallel_efficiency_vs_threads.png", "Efficiency vs threads"),
    ("chart4_dataset_scaling.png",             "Time vs number of genomes"),
    ("chart5_screening_throughput.png",        "Genomes screened per second"),
]

with open(DASHBOARD, "r", encoding="utf-8") as f:
    html = f.read()

for fname, alt in charts:
    path = os.path.join(RESULTS, fname)
    with open(path, "rb") as img:
        b64 = base64.b64encode(img.read()).decode("utf-8")
    pattern = r'src="data:image/png;base64,[^"]+" alt="' + re.escape(alt) + '"'
    replacement = 'src="data:image/png;base64,' + b64 + '" alt="' + alt + '"'
    new_html, n = re.subn(pattern, replacement, html)
    if n:
        html = new_html
        print(f"Updated: {fname}")
    else:
        print(f"WARNING: alt text not matched for {fname}")

with open(DASHBOARD, "w", encoding="utf-8") as f:
    f.write(html)

print("Done — dashboard.html updated with corrected charts.")
