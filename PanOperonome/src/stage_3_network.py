import pandas as pd
import os


# 1. Set input and output paths dynamically

base_dir = os.environ.get("PROJ_BASE", os.path.dirname(os.path.dirname(__file__)))

input_excel = os.path.join(base_dir, "output", "nPO_pgid.xlsx")
output_dir = os.path.join(base_dir, "output", "network_results")

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

print(f"Reading data: {input_excel}")
POdata = pd.read_excel(input_excel)

columns = POdata.columns
strains = [col for col in columns if col != columns[0]]
print(f"Automatically detected {len(strains)} strains: {strains}\n")


# Step 1: Automatically iterate through strains

print("=== Step 1: Generating initial networks for each strain ===")
for target_column in strains:
    print(f"Processing network starting from {target_column} ...")
    output_data = []
    
    for index, row in POdata.iterrows():
        column_value = row[target_column]
        if pd.notna(column_value) and 'operon' in str(column_value):
            row_result = {target_column: column_value}
            
            poid_sequence = [row[columns[0]]] 
            
            bracket_content = str(column_value)[str(column_value).find('(')+1:str(column_value).find(')')]
            genes_in_operon = bracket_content.split(',')
            
            other_values = []
            for val in row:
                if pd.notna(val) and 'operon' not in str(val):
                    other_values.extend(str(val).split(','))

            homo_pg_list = list(set(genes_in_operon).union(set(other_values)))
            search_needed = True

            while search_needed:
                search_needed = False  
                
                for column in strains:
                    matched_entries = []
                    
                    for i, search_value in enumerate(POdata[column]):
                        search_str = str(search_value)
                        if pd.notna(search_value) and 'operon' in search_str:
                            matching_values = [val for val in homo_pg_list if val in search_str]
                            if matching_values:
                                b_start = search_str.find('(')+1
                                b_end = search_str.find(')')
                                bracket_values = search_str[b_start:b_end].split(',')

                                if set(matching_values) == set(bracket_values) or set(bracket_values).issubset(set(matching_values)):
                                    matching_first_col_value = POdata.iloc[i, 0]
                                    poid_sequence.append(matching_first_col_value)
                                    matched_entries.append(f"{search_str}")

                                elif not set(bracket_values).issubset(set(matching_values)):
                                    homo_pg_list = list(set(homo_pg_list) | set(bracket_values))  
                                    column_value = f"operon({','.join(homo_pg_list)})"
                                    row_result[target_column] = column_value  
                                    search_needed = True  

                    if matched_entries:
                        row_result[column] = '|'.join(matched_entries)
                    else:
                        row_result[column] = '-'  

            if len(poid_sequence) > 1:
                row_result['poid'] = '-'.join(str(x) for x in poid_sequence)
            else:
                row_result['poid'] = ''  

            output_data.append(row_result)

    output_df = pd.DataFrame(output_data)
    output_df = output_df.drop_duplicates(subset=output_df.columns[:len(strains)])  
    out_csv = os.path.join(output_dir, f'{target_column}.csv')
    output_df.to_csv(out_csv, index=False)


# Step 2: Merge and deduplicate

print("\n=== Step 2: Merging networks and deduplicating ===")
csv_files = [os.path.join(output_dir, file) for file in os.listdir(output_dir) if file.endswith('.csv') and file.replace('.csv', '') in strains]

df_list = [pd.read_csv(file_path) for file_path in csv_files]
merged_df = pd.concat(df_list, ignore_index=True)

merged_df = merged_df.drop_duplicates(subset=merged_df.columns[:len(strains)])
merged_csv_path = os.path.join(output_dir, 'merged_deduplicated.csv')
merged_df.to_csv(merged_csv_path, index=False)
print(f"Merge complete. Saved to: {merged_csv_path}")


# Step 3: Match PO_ID dictionary

print("\n=== Step 3: Updating PO_ID mapping ===")
POdict = {}
for index, row in POdata.iterrows():
    for column, cell_value in row.items():
        if isinstance(cell_value, str) and 'operon' in cell_value:
            POdict[cell_value] = str(row[columns[0]])
                
POdata1 = pd.read_csv(merged_csv_path)
POdata1['poid'] = ''

for index, row in POdata1.iterrows():
    poid_list = []
    for column_name1, cell_value1 in row.items():
        if isinstance(cell_value1, str):
            if 'operon' in cell_value1 and '|' not in cell_value1:             
                if cell_value1 in POdict:
                    poid_list.append(POdict[cell_value1])
            elif 'operon' in cell_value1 and '|' in cell_value1:
                multi_operon = cell_value1.split('|')
                multi_operon_list = []
                for operon_item in multi_operon:
                    if operon_item in POdict:
                        multi_operon_list.append(POdict[operon_item])
                poid_list.append('|'.join(multi_operon_list))
            
    if poid_list:
        POdata1.at[index, 'poid'] = '-'.join(poid_list)

updated_csv_path = os.path.join(output_dir, 'updated_merged_deduplicated.csv')
POdata1.to_csv(updated_csv_path, index=False)


# Step 4: Combine broken homologous operons

print("\n=== Step 4: Processing broken homologous operons ===")
podata = pd.read_csv(updated_csv_path)
new_data = []

podata_poid_list = podata['poid'].astype(str).tolist()

for index1, row1 in podata.iterrows():
    poid = str(row1['poid'])
    
    if '-' in poid or '1PO' in poid:
        new_data.append(row1.to_dict())
    else:
        for i in range(2, len(strains) + 1):
            if f'{i}PO' in poid:
                poid_count = podata_poid_list.count(poid)
                if poid_count == i:  
                    matching_rows = podata[podata['poid'] == poid]
                    combined_row = row1.to_dict()

                    for index2, row2 in matching_rows.iterrows():
                        if index2 != index1:  
                            for column, cell_value in row2.items():
                                if 'operon' in str(cell_value):
                                    combined_row[column] = cell_value

                    combined_row['poid'] += '*'  
                    new_data.append(combined_row)
                break  

new_df = pd.DataFrame(new_data)
new_df = new_df.drop_duplicates()

final_csv_path = os.path.join(output_dir, 'network_final.csv')
new_df.to_csv(final_csv_path, index=False)

print(f"\nStage 3: Execution complete! Final network file saved to: {final_csv_path}")
