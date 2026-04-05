import pandas as pd
import os
import sys

base_dir = os.environ.get("PROJ_BASE", os.path.dirname(os.path.dirname(__file__)))
path_to_pg_file = os.path.join(base_dir, "data", "reference", "PG.txt")
npo_homo_input = os.path.join(base_dir, "data", "reference", "nPO_homo.xlsx")
final_output = os.path.join(base_dir, "output", "nPO_pgid.xlsx")

print("Reading PG.txt to build dictionary using robust parser...")
pg_id_map = {}

if not os.path.exists(path_to_pg_file):
    print(f"Error: Dictionary file not found at {path_to_pg_file}")
    sys.exit(1)

try:
    pg_df = pd.read_csv(path_to_pg_file, sep=r'\s+', dtype=str, engine='python')
except UnicodeDecodeError:
    pg_df = pd.read_csv(path_to_pg_file, sep=r'\s+', dtype=str, encoding='gbk', engine='python')

for index, row in pg_df.iterrows():
    pg_id = str(row.iloc[0]).strip()
    if pg_id == 'nan' or not pg_id:
        continue
    
    for localtag in row.iloc[1:].dropna():
        localtag_str = str(localtag).strip()
        if localtag_str == '-' or not localtag_str:
            continue
            
        if '|' in localtag_str:
            for tag in localtag_str.split('|'):
                pg_id_map[tag.strip()] = pg_id
        else:
            pg_id_map[localtag_str] = pg_id

if "-" in pg_id_map:  
    del pg_id_map["-"]

print(f"Replacing tags in {os.path.basename(npo_homo_input)}...")
df_homo = pd.read_excel(npo_homo_input)

columns_to_drop = ['Single/Mul', 'single', 'multi', 'sum'] 
existing_cols_to_drop = [col for col in columns_to_drop if col in df_homo.columns]

if existing_cols_to_drop:
    df_homo = df_homo.drop(columns=existing_cols_to_drop)

# CRITICAL FIX: Sort keys by length descending to prevent substring overlap
sorted_keys = sorted(pg_id_map.keys(), key=len, reverse=True)

def partial_replace(cell_value):
    cell_str = str(cell_value)
    if cell_str == 'nan' or cell_str == '-':
        return cell_value
        
    for key in sorted_keys:   
        if key in cell_str:  
            cell_str = cell_str.replace(key, pg_id_map[key])  
    return cell_str 

if hasattr(df_homo, 'map'):
    df_homo = df_homo.map(partial_replace)
else:
    df_homo = df_homo.applymap(partial_replace)

df_homo.to_excel(final_output, index=False)
print(f"Stage 2: File generated successfully at: {final_output}")
