
import pyarrow.parquet as pq
import pandas as pd
import os
import math

in_path = "TGWASmissing.parquet"
out_dir = "split_txt"
os.makedirs(out_dir, exist_ok=True)

# Load schema only (no data loaded)
pq_file = pq.ParquetFile(in_path)
all_cols = pq_file.schema.names

# Split into first 9 and the rest
meta_cols = all_cols[:9]
remaining_cols = all_cols[9:]

group_size = 3  # take 3 columns from the remaining ones
num_groups = math.ceil(len(remaining_cols) / group_size)

print(f"Total columns: {len(all_cols)}")
print(f"Meta columns: {meta_cols}")
print(f"Number of 3-column groups: {num_groups}")

for i in range(num_groups):
    start = i * group_size
    end = min(start + group_size, len(remaining_cols))
    cols = meta_cols + remaining_cols[start:end]
    out_path = os.path.join(out_dir, f"group_{i+1:04d}.txt")

    print(f"Writing columns {start+1}-{end} of remaining -> {out_path}")

    # Read only the selected 9 + 3 columns
    table = pq_file.read(columns=cols)
    df = table.to_pandas()
    df.to_csv(out_path, sep="\t", index=False)

    del table, df
