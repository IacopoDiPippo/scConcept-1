import pandas as pd
import os

def confront(pd1, pd2):
    # extract the only column from each df
    col1 = pd1.iloc[:, 0].astype(str)
    col2 = pd2.iloc[:, 0].astype(str)

    # convert to sets
    s1 = set(col1)
    s2 = set(col2)

    # intersection
    common = s1 & s2

    return len(common)   # number of shared IDs


base_path = os.getcwd()
panels_path = os.path.join(base_path, "done_panels")

files = []
file_names = []
row_counts = {}

# load CSV files properly
for file in os.listdir(panels_path):
    if file.endswith(".csv"):
        full_path = os.path.join(panels_path, file)
        df = pd.read_csv(full_path)

        files.append(df)
        file_names.append(file)
        row_counts[file] = len(df)

# print number of rows in each file
print("=== ROW COUNTS ===")
for name, count in row_counts.items():
    print(f"{name}: {count} rows")

# compare all pairs
for i in range(len(files)):
    for j in range(i + 1, len(files)):
        count = confront(files[i], files[j])
        print(f"{file_names[i]} <-> {file_names[j]}: {count} common IDs")

for i in range(len(files)):
    print("")
