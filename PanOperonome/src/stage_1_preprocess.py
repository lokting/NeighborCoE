import os
import pandas as pd


# 1. Set paths dynamically

base_dir = os.environ.get("PROJ_BASE", os.path.dirname(os.path.dirname(__file__)))

input_folder = os.path.join(base_dir, "output", "temp_po_out")
output_folder = os.path.join(base_dir, "output")

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Automatically detect all strains
strains = []
file_list = []
for file_name in os.listdir(input_folder):
    if file_name.endswith('_PanOperon_new.xlsx'):
        strain = file_name.split('_PanOperon_new.xlsx')[0]
        strains.append(strain)
        file_list.append(file_name)

print(f"Detected {len(strains)} strains for preprocessing: {strains}")

all_records = []


# 2. Extract data and assemble operon(localtag) format

for strain, file_name in zip(strains, file_list):
    file_path = os.path.join(input_folder, file_name)
    print(f"Reading and processing: {strain} ...")
    
    df = pd.read_excel(file_path, usecols=['operon', 'localtag', 'gene'])
    
    def format_val(row):
        if pd.notna(row['localtag']) and str(row['localtag']).strip() != '':
            return f"{row['operon']}({row['localtag']})"
        return str(row['operon'])
        
    df['formatted'] = df.apply(format_val, axis=1)
    
    grouped = df.groupby('gene')['formatted'].apply(lambda x: ', '.join(list(x))).reset_index()
    grouped['strain'] = strain
    
    all_records.append(grouped)


# 3. Merge, calculate, and generate matrix

combined_df = pd.concat(all_records, ignore_index=True)

gene_counts = combined_df['gene'].value_counts()
combined_df['PO_level'] = combined_df['gene'].map(gene_counts)

pivot_df = combined_df.pivot(index='gene', columns='strain', values='formatted')
pivot_df.fillna('-', inplace=True)
pivot_df = pivot_df.reset_index()

pivot_df['PO_level'] = pivot_df['gene'].map(gene_counts)
pivot_df.sort_values(by=['PO_level', 'gene'], inplace=True)


# 4. Generate continuous PO_IDs

po_ids = []
counter = 1
for _, row in pivot_df.iterrows():
    level = row['PO_level']
    po_ids.append(f"{level}PO{counter:04d}")
    counter += 1

pivot_df.insert(0, 'PO_ID', po_ids)

final_cols = ['PO_ID'] + strains
final_df = pivot_df[final_cols]


# 5. Save the final file

output_path = os.path.join(output_folder, 'allPO.xlsx')
final_df.to_excel(output_path, index=False)

print(f"Stage 1: Data processing complete! File saved to: {output_path}")

