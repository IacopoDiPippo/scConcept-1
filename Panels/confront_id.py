import pandas as pd
import os
from tabulate import tabulate

def confront(pd1, pd2):
    # assume first column contains gene IDs (ENSMUSG)
    col1 = pd1.iloc[:, 0].astype(str)
    col2 = pd2.iloc[:, 0].astype(str)

    return len(set(col1) & set(col2))


base_path = os.getcwd()
panels_path = os.path.join(base_path, "done_panels")

files = []
labels = []
row_counts = {}

# load CSV files
for file in sorted(os.listdir(panels_path)):
    if file.endswith(".csv"):
        full_path = os.path.join(panels_path, file)
        df = pd.read_csv(full_path)

        files.append(df)

        # use filename without .csv as label
        label = os.path.splitext(file)[0]
        labels.append(label)

        row_counts[label] = len(df)

# print number of rows
print("=== ROW COUNTS ===")
for name, count in row_counts.items():
    print(f"{name}: {count} rows")

# create overlap matrix
overlap_matrix = pd.DataFrame(
    index=labels,
    columns=labels,
    dtype=int
)

# fill overlap matrix
for i in range(len(files)):
    for j in range(len(files)):
        overlap_matrix.iloc[i, j] = confront(files[i], files[j])

# print nicely
print("\n=== OVERLAP MATRIX (shared gene IDs) ===\n")
print(
    tabulate(
        overlap_matrix,
        headers="keys",
        tablefmt="grid"
    )
)
