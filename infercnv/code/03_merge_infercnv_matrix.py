import os
import pandas as pd

# ===================== Parameter Section =====================
input_dir = "/hanlab/xueqinliu/test/infercnv/02expr_z"
output_file = "/hanlab/xueqinliu/test/infercnv/result.csv"
# =============================================================

# Get full path for every csv file
csv_files = [
    os.path.join(input_dir, f)
    for f in os.listdir(input_dir)
    if f.endswith(".csv")
]
csv_files = sorted(csv_files)

if len(csv_files) == 0:
    raise FileNotFoundError(f"No CSV files found in {input_dir}")

# Initialize
df_merged = pd.read_csv(csv_files[0], index_col=0)
df_merged["key"] = df_merged.index
print(f"Start with: {csv_files[0]}, shape = {df_merged.shape}")

# Merge all samples
for file_path in csv_files[1:]:
    df_sub = pd.read_csv(file_path, index_col=0)
    df_sub["key"] = df_sub.index
    df_merged = pd.merge(df_merged, df_sub, how="outer", on="key")
    print(f"Merged {os.path.basename(file_path)}, current shape: {df_merged.shape}")

df_merged = df_merged.set_index("key")

# Make output folder
out_dir = os.path.dirname(output_file)
os.makedirs(out_dir, exist_ok=True)
df_merged.to_csv(output_file)
print(f"\nDone! Final matrix shape: {df_merged.shape}")
print(f"Saved to: {output_file}")

