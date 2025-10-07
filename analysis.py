"""
RT-qPCR siRNA Validation Analysis
Processes qPCR data from siRNA knockdown experiments to calculate knockdown efficiency
Author: Joanna Renaut (adapted from an R script by Dr. Robert Zach)
Date: 7th October 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
import warnings
import os

warnings.filterwarnings('ignore')

# ============================================================================
# USER CONFIGURATION - EDIT THESE PARAMETERS
# ============================================================================

# Paths
path_data = Path('data')
path_results = Path('results')
path_figures = Path('figures')

# Generate output filename
output_name = "sample_output" # output name for .csv results and figures (.png/.pdf)

# create results output
os.makedirs(path_results, exist_ok=True)
os.makedirs(path_figures, exist_ok=True)

# Sample names for controls (will be excluded from plots but used for normalization)
control_sample_ntc = 'NTC'  # No template control
control_sample_nrt = 'NRT'  # No reverse transcriptase control
control_sample_scrambled = 'Scrm'  # Scrambled siRNA control

# Reference genes (housekeeping genes - won't be plotted as knockdown targets)
reference_gene = 'GAPDH'  # Primary reference gene for normalization
housekeeping_genes = ['ACTB']  # Additional housekeeping genes to exclude from analysis

# Column name mappings (if your CSV uses different column names)
# Leave as-is if using standard format
column_name_sample = 'siRNA Sample'  # Column containing sample names
column_name_primers = 'Primers'      # Column containing primer/gene names
column_name_cell_line = 'Cell Line'  # Column containing cell line names
column_name_cq = 'Cq'                # Column containing Cq values

# Genes to analyze (leave empty list [] to automatically detect all genes in data)
genes_to_plot = []  # e.g., ['BRCA']

# Cell lines to analyze (leave empty list [] to automatically detect all cell lines)
cell_lines_to_plot = []  # e.g., ['GENE_KO']

# Figure settings
y_axis_max = 1.5  # Maximum y-axis value (adjust to see your data range)
label_threshold = 0.2  # Values below this will be labeled with their exact value
fig_width_per_column = 0.8  # Width in inches per cell line column
fig_height_per_row = 0.4    # Height in inches per gene row

# Spacing between plots
horizontal_spacing = 0.10  # Space between cell line columns
vertical_spacing = 0.15    # Space between gene rows

# Colour palette for genes (one color per gene row)
colour_palette = ["#D4D3CF", "#DC6B83", "#75B1CE", "#D8C367", "#526C94",
                  "#000000", "#ccdba2", "#889466", "#E69F00", "#56B4E9",
                  "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442"]

# ============================================================================
# END USER CONFIGURATION
# ============================================================================

# Build control samples list
control_samples = [control_sample_ntc, control_sample_nrt]

# Set up matplotlib/seaborn styling
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 6
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.25
plt.rcParams['ytick.major.width'] = 0.25
plt.rcParams['xtick.major.size'] = 2
plt.rcParams['ytick.major.size'] = 2

# Create results directory if it doesn't exist
path_results.mkdir(exist_ok=True)

# ============================================================================
# DATA IMPORT AND PROCESSING
# ============================================================================

print("=" * 60)
print("RT-qPCR siRNA Validation Analysis")
print("=" * 60)

# Import all CSV files from data folder
data_raw = pd.DataFrame()

for file in path_data.glob("*.csv"):
    tmp_data = pd.read_csv(file)

    # Strip whitespace from column names (handles trailing spaces)
    tmp_data.columns = tmp_data.columns.str.strip()

    # Replace the column mapping section (around lines 96-121) with this:

    # Rename columns to standardized names using user-defined mappings
    column_mapping = {}
    for col in tmp_data.columns:
        if column_name_sample in col:
            column_mapping[col] = 'sample'
        elif column_name_primers in col:
            column_mapping[col] = 'primers'
        elif column_name_cell_line in col:
            column_mapping[col] = 'cell_line'

    # Handle Cq column separately to avoid mapping both 'Cq' and 'Average Cq'
    if column_name_cq in tmp_data.columns:
        column_mapping[column_name_cq] = 'Cq'

    tmp_data = tmp_data.rename(columns=column_mapping)

    # Check for required columns
    if 'primers' not in tmp_data.columns or 'sample' not in tmp_data.columns or 'Cq' not in tmp_data.columns:
        print(f"⚠ Warning: Required columns not found in {file.name}")
        print(f"   Need: 'sample', 'primers', 'Cq'")
        print(f"   Have: {tmp_data.columns.tolist()}")
        continue

    tmp_data['experiment'] = file.name
    tmp_data['Cq'] = pd.to_numeric(tmp_data['Cq'], errors='coerce')

    # Select only needed columns
    cols_to_keep = ['experiment', 'sample', 'primers', 'Cq']
    if 'cell_line' in tmp_data.columns:
        cols_to_keep.append('cell_line')

    tmp_data = tmp_data[cols_to_keep]

    # Remove rows with missing data
    tmp_data = tmp_data.dropna(subset=['sample', 'primers'])
    if 'cell_line' in tmp_data.columns:
        tmp_data = tmp_data.dropna(subset=['cell_line'])

    # Group by and average technical replicates
    group_cols = ['experiment', 'sample', 'primers']
    if 'cell_line' in tmp_data.columns:
        group_cols.append('cell_line')

    tmp_data = tmp_data.groupby(group_cols, as_index=False)['Cq'].mean()

    data_raw = pd.concat([data_raw, tmp_data], ignore_index=True)

print(f"\n✓ Data loaded: {len(data_raw)} measurements")

# Remove NaN Cq values
data_processed = data_raw.dropna(subset=['Cq']).copy()
print(f"✓ After filtering NaN: {len(data_processed)} measurements")

# Calculate RNA levels from Cq values (2^-Cq)
data_processed['RNA'] = 2 ** (-data_processed['Cq'])

# ============================================================================
# NORMALIZATION TO REFERENCE GENE
# ============================================================================

# Extract Scrm+reference gene values (universal reference from Controls section)
scrm_ref_ref = data_processed[
    (data_processed['sample'] == control_sample_scrambled) &
    (data_processed['primers'] == reference_gene)
].copy()

if 'cell_line' in scrm_ref_ref.columns:
    scrm_ref_ref = scrm_ref_ref.groupby(['experiment', 'cell_line'], as_index=False)['RNA'].mean()
    scrm_ref_ref = scrm_ref_ref.rename(columns={'RNA': f'{control_sample_scrambled}_{reference_gene}_RNA'})
else:
    scrm_ref_ref = scrm_ref_ref.groupby(['experiment'], as_index=False)['RNA'].mean()
    scrm_ref_ref = scrm_ref_ref.rename(columns={'RNA': f'{control_sample_scrambled}_{reference_gene}_RNA'})

print(f"\n✓ {control_sample_scrambled}+{reference_gene} reference values extracted")

# Get reference gene values for each sample
id_cols = ['experiment', 'sample']
if 'cell_line' in data_processed.columns:
    id_cols.append('cell_line')

ref_gene_values = data_processed[data_processed['primers'] == reference_gene][id_cols + ['RNA']].copy()
ref_gene_values = ref_gene_values.rename(columns={'RNA': f'{reference_gene}_RNA'})

# Get all gene measurements (exclude reference and housekeeping genes)
genes_to_exclude = [reference_gene] + housekeeping_genes
gene_measurements = data_processed[
    ~data_processed['primers'].isin(genes_to_exclude)
].copy()

# Merge gene measurements with their reference gene values
gene_measurements = gene_measurements.merge(ref_gene_values, on=id_cols, how='left')

# Calculate relative expression (gene/reference)
gene_measurements['RNA_rel'] = gene_measurements['RNA'] / gene_measurements[f'{reference_gene}_RNA']

# ============================================================================
# NORMALIZATION TO SCRAMBLED CONTROL
# ============================================================================

# Get scrambled control values for each gene (for normalization)
merge_cols = ['experiment']
if 'cell_line' in gene_measurements.columns:
    merge_cols.append('cell_line')

scrm_gene_values = gene_measurements[
    gene_measurements['sample'] == control_sample_scrambled
][merge_cols + ['primers', 'RNA_rel']].copy()
scrm_gene_values = scrm_gene_values.rename(columns={'RNA_rel': f'{control_sample_scrambled}_gene_RNA_rel'})

# Merge with main data
gene_measurements = gene_measurements.merge(
    scrm_gene_values,
    on=merge_cols + ['primers'],
    how='left'
)

# Normalize to scrambled control (Scrm becomes 1.0, knockdowns become fraction of Scrm)
gene_measurements['RNA_rel_norm'] = gene_measurements['RNA_rel'] / gene_measurements[f'{control_sample_scrambled}_gene_RNA_rel']

print(f"✓ Normalized to {control_sample_scrambled} control: {len(gene_measurements)} gene measurements")

# ============================================================================
# DETERMINE GENES AND CELL LINES TO PLOT
# ============================================================================

# Auto-detect or use specified genes
if len(genes_to_plot) == 0:
    genes = sorted([g for g in gene_measurements['primers'].unique()
                   if g not in control_samples])
    print(f"\n✓ Auto-detected genes: {genes}")
else:
    genes = genes_to_plot
    print(f"\n✓ Using specified genes: {genes}")

# Auto-detect or use specified cell lines
if 'cell_line' in gene_measurements.columns:
    if len(cell_lines_to_plot) == 0:
        cell_lines = sorted(gene_measurements['cell_line'].unique())
        print(f"✓ Auto-detected cell lines: {cell_lines}")
    else:
        cell_lines = cell_lines_to_plot
        print(f"✓ Using specified cell lines: {cell_lines}")
else:
    cell_lines = ['All']
    print("✓ No cell line column detected")

# ============================================================================
# EXPORT RESULTS TO CSV
# ============================================================================

# Create summary table with knockdown efficiencies
summary_data = []

for gene in genes:
    for cell_line in cell_lines:
        # Filter data for this gene and cell line
        if 'cell_line' in gene_measurements.columns and cell_line != 'All':
            plot_data = gene_measurements[
                (gene_measurements['primers'] == gene) &
                (gene_measurements['cell_line'] == cell_line)
                ].copy()
        else:
            plot_data = gene_measurements[
                gene_measurements['primers'] == gene
                ].copy()

        # Get Scrm and siRNA values
        scrm_val = plot_data[plot_data['sample'] == control_sample_scrambled]['RNA_rel_norm'].mean()
        sirna_val = plot_data[plot_data['sample'] == gene]['RNA_rel_norm'].mean()

        # Calculate knockdown percentage
        if not pd.isna(scrm_val) and not pd.isna(sirna_val) and scrm_val > 0:
            remaining_pct = (sirna_val / scrm_val) * 100
            knockdown_pct = 100 - remaining_pct
        else:
            remaining_pct = np.nan
            knockdown_pct = np.nan

        summary_data.append({
            'Gene': gene,
            'Cell_Line': cell_line,
            'Scrm_Normalized': scrm_val,
            'siRNA_Normalized': sirna_val,
            'Remaining_%': remaining_pct,
            'Knockdown_%': knockdown_pct
        })

# Create DataFrame and save
summary_df = pd.DataFrame(summary_data)
summary_filename = output_name.replace('graph_siRNA', 'knockdown_efficiency') + '.csv'
summary_df.to_csv(path_results / summary_filename, index=False)

print(f"✓ Summary saved: {path_results / summary_filename}")

# ============================================================================
# CREATE FIGURE
# ============================================================================

print(f"\n{'=' * 60}")
print("Generating figure...")
print(f"{'=' * 60}")

n_rows = len(genes)
n_cols = len(cell_lines)

# Calculate figure size based on number of panels and user settings
fig_width = 0.5 + (n_cols * fig_width_per_column)
fig_height = 0.5 + (n_rows * fig_height_per_row)
fig = plt.figure(figsize=(fig_width, fig_height))

# Create grid for subplots using user-defined spacing
gs = GridSpec(n_rows, n_cols, figure=fig,
              hspace=vertical_spacing,
              wspace=horizontal_spacing,
              left=0.18, right=0.90, top=0.92, bottom=0.08)

# Plot each gene × cell line combination
for i, gene in enumerate(genes):
    for j, cell_line in enumerate(cell_lines):
        ax = fig.add_subplot(gs[i, j])

        # Filter data for this specific gene and cell line
        if 'cell_line' in gene_measurements.columns and cell_line != 'All':
            plot_data = gene_measurements[
                (gene_measurements['primers'] == gene) &
                (gene_measurements['cell_line'] == cell_line)
            ].copy()
        else:
            plot_data = gene_measurements[
                gene_measurements['primers'] == gene
            ].copy()

        # Determine which samples are present (Scrm and/or gene knockdown)
        samples_present = [control_sample_scrambled, gene]
        samples_present = [s for s in samples_present if s in plot_data['sample'].values]

        x_pos = np.arange(len(samples_present))

        # Calculate mean values for each sample
        values = []
        for s in samples_present:
            sample_data = plot_data[plot_data['sample'] == s]['RNA_rel_norm']
            if len(sample_data) > 0:
                values.append(sample_data.mean())
            else:
                values.append(0)

        # Create bar plot
        if len(values) > 0:
            ax.bar(x_pos, values, width=0.8,
                  color=colour_palette[i % len(colour_palette)],
                  edgecolor='#000000',
                  linewidth=0.25)

        # Add value labels on small bars (for readability)
        for idx, (x, val) in enumerate(zip(x_pos, values)):
            if val < label_threshold:
                ax.text(x, val + 0.05, f'{val:.2f}',
                       ha='center', va='bottom', fontsize=4)

        # Set y-axis limits and ticks using user-defined max
        ax.set_ylim(0, y_axis_max)
        y_ticks = np.linspace(0, y_axis_max, 4)
        ax.set_yticks(y_ticks)
        ax.set_xticks(x_pos)

        # Style all spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)
            spine.set_color('#000000')

        # Column headers (cell line names) - only on top row
        if i == 0:
            ax.set_title(cell_line, fontsize=6, pad=2)

        # Y-axis labels - only on leftmost column
        if j == 0:
            ax.tick_params(axis='y', labelsize=5, length=2, width=0.25, pad=1)
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis='y', length=2, width=0.25)

        # Row labels (gene names) - only on rightmost column
        if j == n_cols - 1:
            ax.text(1.12, 0.5, gene, transform=ax.transAxes,
                   ha='left', va='center', fontsize=6, rotation=270)

        # X-axis labels - only on bottom row
        if i == n_rows - 1:
            x_labels = [control_sample_scrambled, 'siRNA'] if len(samples_present) == 2 else [control_sample_scrambled]
            ax.set_xticklabels(x_labels, rotation=90, ha='center', va='top', fontsize=5.5)
            ax.tick_params(axis='x', labelsize=5.5, length=2, width=0.25, pad=0.5)
        else:
            ax.set_xticklabels([])
            ax.tick_params(axis='x', length=0, width=0.25)

# Add overall figure labels
if 'cell_line' in gene_measurements.columns:
    fig.text(0.5, 0.97, 'Cell Line', ha='center', va='bottom', fontsize=6)
else:
    fig.text(0.5, 0.97, 'siRNA Knockdown', ha='center', va='bottom', fontsize=6)

fig.text(0.02, 0.5, f'RNA relative to {reference_gene} (normalised to {control_sample_scrambled})',
         ha='center', va='center', fontsize=6, rotation=90)

# ============================================================================
# SAVE FIGURE
# ============================================================================

if len(cell_lines_to_plot) > 0:
    output_name += f"_{'_'.join(cell_lines_to_plot).replace(' ', '_')}"

# Save in multiple formats
plt.savefig(path_figures / f"{output_name}.pdf", dpi=600, bbox_inches='tight')
plt.savefig(path_figures / f"{output_name}.png", dpi=600, bbox_inches='tight')
plt.show()

print(f"\n✓ Figure saved: {path_figures / output_name}.pdf")
print(f"✓ Figure saved: {path_figures / output_name}.png")
print(f"\nGenes plotted: {genes}")
print(f"Cell lines plotted: {cell_lines}")
print(f"\n{'=' * 60}")
print("Analysis complete!")
print(f"{'=' * 60}\n")