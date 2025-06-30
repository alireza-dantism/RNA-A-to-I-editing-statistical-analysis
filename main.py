import pandas as pd
import scipy
from scipy import stats
import glob
import numpy as np
import matplotlib.pyplot as plt

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
    strand = gene['strand']
    
    # Filter samples for the genomic region of interest and strand
    gene_samples = all_samples[
        (all_samples['#CHR'] == gene['chr']) &
        (all_samples['POS'] >= gene['start']) &
        (all_samples['POS'] <= gene['end'])
    ]
    
    if not gene_samples.empty:
        # Calculate mean and standard deviation of editing levels for each condition and strand
        mean_editing = gene_samples.groupby(['Condition'])['Editing_Level'].mean()
        std_editing = gene_samples.groupby(['Condition'])['Editing_Level'].std()

        # Store results for both strands
        results[gene_id] = {
            'strand': strand,
            'mean_editing': mean_editing.to_dict(),
            'std_dev': std_editing.to_dict()
        }

        # Perform statistical testing only if both conditions have data
        if (gene_samples['Condition'] == 'Mettl3').any() and (gene_samples['Condition'] == 'Mettl3_ctrl').any():
            # Prepare data for statistical testing
            editing_levels_mettl3 = gene_samples[gene_samples['Condition'] == 'Mettl3'][['POS', 'Editing_Level']]
            editing_levels_ctrl = gene_samples[gene_samples['Condition'] == 'Mettl3_ctrl'][['POS', 'Editing_Level']]

            # Merge the two DataFrames on the POS column to find common positions
            merged_data = pd.merge(editing_levels_mettl3, editing_levels_ctrl, on='POS', suffixes=('_Mettl3', '_Ctrl'))

            # Check if there are sufficient paired data points
            if len(merged_data) > 1:  # at least 2 points needed for a paired t-test
                # Perform the paired t-test
                t_stat, p_value = stats.ttest_rel(merged_data['Editing_Level_Mettl3'], merged_data['Editing_Level_Ctrl'])
                results[gene_id]['p_value'] = p_value
            else:
                results[gene_id]['p_value'] = np.nan  # Not enough data for comparison
        else:
            results[gene_id]['p_value'] = np.nan  # Not enough data for comparison

# Open a file to write the results
with open('change_results_strand_specific.txt', 'w') as file:
    file.write("Human-Readable Summary of Results (Strand-Specific)\n")
    file.write("=" * 50 + "\n")
    
    for gene_id, data in results.items():
        strand = data['strand']
        mean_editing = data['mean_editing']
        std_dev = data['std_dev']
        p_value = data.get('p_value', np.nan)

        # Handle NaN values for p-value
        p_value_str = "NaN" if np.isnan(p_value) else f"{p_value:.4f}"
        
        file.write(f"Gene ID: {gene_id} (Strand: {strand})\n")
        for condition, mean in mean_editing.items():
            std = std_dev.get(condition, 0)
            file.write(f"- Mean Editing Level ({condition}): {mean:.4f} (SD: {std:.4f})\n")
        file.write(f"- p-value: {p_value_str}\n")
        
        # Interpretation
        if not np.isnan(p_value):
            if p_value < 0.05:
                file.write("- Statistically significant difference\n\n")
            else:
                file.write("- No significant difference\n\n")
        else:
            file.write("- Not enough data for statistical comparison\n\n")

# Visualization of Editing Levels
for gene in results.keys():
    gene_data = all_samples[
        (all_samples['#CHR'] == genes_of_interest.loc[genes_of_interest['gene_id'] == gene, 'chr'].values[0]) &
        (all_samples['POS'] >= genes_of_interest.loc[genes_of_interest['gene_id'] == gene, 'start'].values[0]) &
        (all_samples['POS'] <= genes_of_interest.loc[genes_of_interest['gene_id'] == gene, 'end'].values[0])
    ]
    
    if not gene_data.empty:
        plt.figure(figsize=(8, 5))
        plt.title(f"Editing Levels for {gene} (Strand: {results[gene]['strand']})")
        for condition in gene_data['Condition'].unique():
            plt.boxplot(gene_data[gene_data['Condition'] == condition]['Editing_Level'], 
                        positions=[1 if condition == 'Mettl3' else 2], 
                        widths=0.5, 
                        labels=[condition])
        
        plt.ylabel("Editing Level")
        plt.xticks([1, 2], ['Mettl3', 'Mettl3_ctrl'])
        plt.savefig(f"{gene}_editing_levels_strand_specific.png")
        plt.close()

print("Strand-specific results written to change_results_strand_specific.txt and visualizations saved.")
