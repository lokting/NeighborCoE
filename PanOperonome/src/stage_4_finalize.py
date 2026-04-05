import os
import sys
import pandas as pd
import re
import shutil

base_dir = os.environ.get("PROJ_BASE", os.path.dirname(os.path.dirname(__file__)))
output_dir = os.path.join(base_dir, "output")
temp_dir = os.path.join(output_dir, "temp_po_out")
network_dir = os.path.join(output_dir, "network_results")
network_path = os.path.join(network_dir, "network_final.csv")
pg_file = os.path.join(base_dir, "data", "reference", "PG.txt")

print("=== Stage 4: Absolute Operon Correction ===")

if not os.path.exists(network_path):
    print(f"Error: Cannot find {network_path}.")
    sys.exit(1)

df = pd.read_csv(network_path)
strain_cols = [col for col in df.columns if col != 'poid']

# 1. Build the Absolute Source of Truth from Stage 0
print("Loading absolute truth from original operon files...")
operon_truth = {strain: {} for strain in strain_cols}
pg_to_loc_truth = {strain: {} for strain in strain_cols}

for strain in strain_cols:
    temp_file = os.path.join(temp_dir, f"{strain}_PanOperon_new.xlsx")
    if os.path.exists(temp_file):
        temp_df = pd.read_excel(temp_file)
        for _, row in temp_df.iterrows():
            # Extract Operon ID (e.g., operon_1789)
            op_raw = str(row['operon']).strip()
            op_id = op_raw if op_raw.startswith('operon_') else f"operon_{op_raw}"
            
            # This is the pure, unaltered original local tag (e.g., STM14_RS24910)
            loc_tags = str(row['localtag']).strip()
            operon_truth[strain][op_id] = loc_tags
            
            # Also map individual genes for any standalone PG IDs
            pg_genes = str(row['gene']).strip().split(',')
            loc_list = loc_tags.split(',')
            if len(pg_genes) == len(loc_list):
                for p, l in zip(pg_genes, loc_list):
                    pg_to_loc_truth[strain][p.strip().replace('*', '')] = l.strip()

# 2. Build PG.txt Fallback (for merged genes from other strains)
global_pg_map = {}
try:
    pg_df = pd.read_csv(pg_file, sep=r'\s+', dtype=str, engine='python')
    for _, row in pg_df.iterrows():
        pid = str(row.iloc[0]).strip()
        tags = [str(v).strip() for v in row.iloc[1:].dropna() if str(v).strip() != '-']
        if tags: 
            # Drop the | to prevent clutter in non-operon cells
            global_pg_map[pid] = tags[0].split('|')[0] 
except: pass

# 3. The Core Overwrite Function
def final_revert(cell_val, strain_name):
    if pd.isna(cell_val) or str(cell_val) == '-': return cell_val
    cell_str = str(cell_val)
    
    # --- Fix 1: Absolute Operon Overwrite (The User's Logic) ---
    def fix_operon_match(match):
        op_id = match.group(1)   # e.g., operon_1789
        content = match.group(2) # e.g., STM14_RS24910|STM14_RS06995
        
        # If we know the exact composition of this operon, FORCE overwrite it!
        if op_id in operon_truth[strain_name]:
            return f"{op_id}({operon_truth[strain_name][op_id]})"
        
        # If the operon is foreign (merged), at least clean up the '|'
        if '|' in content:
            cleaned = [item.split('|')[0].strip() for item in content.split(',')]
            return f"{op_id}({','.join(cleaned)})"
            
        return match.group(0)

    # Apply the operon fix to all operons in the cell
    cell_str = re.sub(r'(operon_[A-Za-z0-9_-]+)\(([^)]+)\)', fix_operon_match, cell_str)
    
    # --- Fix 2: Translate any remaining raw PG IDs (e.g., 3PG0996) ---
    def fix_standalone_pg(match):
        pg_id = match.group(0).replace('*', '')
        if pg_id in pg_to_loc_truth[strain_name]:
            return pg_to_loc_truth[strain_name][pg_id]
        if pg_id in global_pg_map:
            return global_pg_map[pg_id]
        return pg_id
        
    cell_str = re.sub(r'\d+PG\d+', fix_standalone_pg, cell_str)
    
    return cell_str

print("Applying absolute replacements...")
for strain in strain_cols:
    df[strain] = df[strain].apply(lambda x: final_revert(x, strain))

# 4. Final Recount and Cleanup
print("Finalizing file structure...")
df['lvl'] = df[strain_cols].apply(lambda r: sum(1 for x in r if pd.notna(x) and 'operon' in str(x)), axis=1)
df = df.sort_values('lvl').reset_index(drop=True)

counts = {}
poids = []
for l in df['lvl']:
    counts[l] = counts.get(l, 0) + 1
    poids.append(f"{l}PO{counts[l]}")
df['poid'] = poids
df[['poid'] + strain_cols].to_csv(network_path, index=False)

# Cleanup
safe_file = os.path.join(base_dir, "final_backup.csv")
shutil.copy2(network_path, safe_file)
if os.path.exists(output_dir): shutil.rmtree(output_dir)
os.makedirs(network_dir)
shutil.move(safe_file, network_path)

print("Done! Check network_results/network_final.csv")