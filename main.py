import pandas as pd
import scipy
from scipy import stats
import glob
import numpy as np

# Load the genes of interest
genes_of_interest = pd.read_csv('data/genes_of_interest.tsv', sep='\t')

# Load Mettl3 samples
mett_samples = []
for file in glob.glob('data/Mettl3_*.tsv'):
    sample_data = pd.read_csv(file, sep='\t')
    sample_data['Condition'] = 'Mettl3'
    mett_samples.append(sample_data)

# Load Mettl3 control samples
mett_ctrl_samples = []
for file in glob.glob('data/Mettl3_ctrl_*.tsv'):
    ctrl_data = pd.read_csv(file, sep='\t')
    ctrl_data['Condition'] = 'Mettl3_ctrl'
    mett_ctrl_samples.append(ctrl_data)

# Combine all samples into one DataFrame
all_samples = pd.concat(mett_samples + mett_ctrl_samples, ignore_index=True)

# Remove rows with zero Good_depth
all_samples = all_samples[all_samples['Good_depth'] > 0]

# Calculate total reads and editing levels
all_samples['Total'] = all_samples[['Count_A', 'Count_C', 'Count_G', 'Count_T']].sum(axis=1)
all_samples['Editing_Level'] = all_samples['Count_G'] / all_samples['Total']

# Prepare a DataFrame for the results
results = {}

# Loop through each gene of interest
for index, gene in genes_of_interest.iterrows():
    gene_id = gene['gene_id']
    # Filter samples for the genomic region of interest
    gene_samples = all_samples[
        (all_samples['#CHR'] == gene['chr']) &
        (all_samples['POS'] >= gene['start']) &
        (all_samples['POS'] <= gene['end'])
    ]
    
    if not gene_samples.empty:
        # Calculate mean editing levels for each condition
        mean_editing = gene_samples.groupby('Condition')['Editing_Level'].mean()
        
        # Store mean editing change
        mean_change = mean_editing.get('Mettl3', 0) - mean_editing.get('Mettl3_ctrl', 0)
        results[gene_id] = {
            'mean_editing_change': mean_change
        }
        
        # Perform statistical testing only if both conditions have data
        if len(gene_samples[gene_samples['Condition'] == 'Mettl3']) > 0 and len(gene_samples[gene_samples['Condition'] == 'Mettl3_ctrl']) > 0:
            
            # Prepare data for statistical testing
            editing_levels_mettl3 = gene_samples[gene_samples['Condition'] == 'Mettl3'][['POS', 'Editing_Level']]
            editing_levels_ctrl = gene_samples[gene_samples['Condition'] == 'Mettl3_ctrl'][['POS', 'Editing_Level']]

            # Merge the two DataFrames on the POS column to find common positions
            merged_data = pd.merge(editing_levels_mettl3, editing_levels_ctrl, on='POS', suffixes=('_Mettl3', '_Ctrl'))

            # Check if there are sufficient paired data points
            if len(merged_data) > 0:
                # Perform the paired t-test
                t_stat, p_value = stats.ttest_rel(merged_data['Editing_Level_Mettl3'], merged_data['Editing_Level_Ctrl'])
                results[gene_id]['p_value'] = p_value
            else:
                print(f"No sufficient paired data for {gene_id}.")

# Open a file to write the results
with open('change_results.txt', 'w') as file:
    file.write("Human-Readable Summary of Editing Change Results\n")
    file.write("=" * 50 + "\n")
    
    for gene_id, data in results.items():
        mean_editing_change = data['mean_editing_change']
        p_value = data['p_value']
        
        # Handle NaN values for p-value
        p_value_str = "NaN" if np.isnan(p_value) else f"{p_value:.4f}"
        
        file.write(f"Gene ID: {gene_id}\n")
        file.write(f"- Mean Editing Change: {mean_editing_change:.4f}\n")
        file.write(f"- p-value: {p_value_str}\n")
        
        # Interpretation
        if not np.isnan(p_value):
            if p_value < 0.05:
                file.write("- Statistically significant difference\n\n")
            else:
                file.write("- No significant difference\n\n")
        else:
            file.write("- Not enough data for statistical comparison\n\n")

print("Results written to change_results.txt.")