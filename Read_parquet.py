import pandas as pd
pd.set_option('display.max_columns', None)

# Show all rows (optional)
pd.set_option('display.max_rows', None)
# Read the parquet file
df = pd.read_parquet("TGWASmissing-x1-x10.parquet")
print("Rows:", df.shape[0])
print("Columns:", df.shape[1])
print(df.iloc[1:3, :])
# Loop through columns in pairs
# for i in range(0, 13, 2): #df.shape[1]
#     # Select two columns at a time
#     cols = df.iloc[:, [i, i+1]] if i+1 < df.shape[1] else df.iloc[:, [i]]
    
#     base_name = df.columns[i]

#     # Create filename using that base name
#     filename = f"TGWAS_{base_name}.tsv"
    
#     # Save to TSV
#     cols.to_csv(filename, sep="\t", index=False)

#     print(f"Saved {filename} with columns {df.columns[i:i+2].tolist()}")
