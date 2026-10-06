import os
import csv

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
ORIGINAL_CSV = os.path.join(DATA_DIR, "veterinary_amr_kmer_dataset.csv")
SCALED_DIR = os.path.join(DATA_DIR, "scaled")

os.makedirs(SCALED_DIR, exist_ok=True)

def generate_scaled_datasets():
    print(f"Reading original dataset: {ORIGINAL_CSV}")
    with open(ORIGINAL_CSV, mode="r", encoding="utf-8") as f:
        reader = list(csv.reader(f))
        header = reader[0]
        rows = reader[1:]

    total_available = len(rows)
    print(f"Total available genomes: {total_available}")

    target_sizes = [100, 200, 400, 800, 1000]

    for size in target_sizes:
        if size <= total_available:
            subset = rows[:size]
        else:
            # Replicate if size > available
            subset = (rows * ((size // total_available) + 1))[:size]
            
        out_path = os.path.join(SCALED_DIR, f"amr_{size}.csv")
        with open(out_path, mode="w", newline="", encoding="utf-8") as out_f:
            writer = csv.writer(out_f)
            writer.writerow(header)
            writer.writerows(subset)
            
        print(f"  -> Generated {out_path} ({len(subset)} genomes)")

if __name__ == "__main__":
    generate_scaled_datasets()
    print("All scaled datasets generated successfully in data/scaled/!")
