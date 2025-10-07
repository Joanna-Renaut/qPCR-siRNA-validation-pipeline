# RT-qPCR siRNA Validation Analysis

Python script for analyzing RT-qPCR data from siRNA knockdown experiments. Fully configurable for different experimental setups - simply edit the configuration section at the top of the script.

## Overview

This script automates the analysis of qPCR data to:
- Calculate relative gene expression normalized to a housekeeping gene (default: GAPDH)
- Normalize knockdown samples to scrambled siRNA controls (default: Scrm)
- Generate publication-quality multi-panel figures
- Export knockdown efficiency data to CSV
- Automatically average technical replicates
- Handle multiple cell lines and genes in a single run

## Requirements

Python Packages (install with pip):
- pandas
- numpy
- matplotlib

Python version: 3.12+

## Quick Start

1. Place your CSV file(s) in the data/ folder
2. Edit configuration at the top of analysis.py
3. Run: python analysis.py
4. Find figures in figures/ folder and results CSV in results/ folder

## Project Structure

    project/
    ├── analysis.py
    ├── data/               # Place your CSV files here
    │   └── experiment1.csv
    ├── results/            # CSV output files saved here
    │   └── knockdown_efficiency.csv
    └── figures/            # PNG and PDF figures saved here
        ├── graph_siRNA.pdf
        └── graph_siRNA.png

Folders are automatically created if they don't exist.

## Configuration Guide

All settings are in the USER CONFIGURATION section at the top of analysis.py.

### 1. Paths and Output Name

    path_data = Path('data')
    path_results = Path('results')
    path_figures = Path('figures')
    output_name = "graph_siRNA"

- path_data: Folder containing input CSV files
- path_results: Folder for CSV output
- path_figures: Folder for figure output
- output_name: Base name for all output files

Example outputs with output_name = "my_experiment":
- my_experiment.pdf
- my_experiment.png
- my_experiment_efficiency.csv

### 2. Control Sample Names

Edit these to match YOUR control names:

    control_sample_ntc = 'NTC'
    control_sample_nrt = 'NRT'
    control_sample_scrambled = 'Scrm'

Examples:
- control_sample_scrambled = 'siControl' if you use siControl
- control_sample_ntc = 'H2O' if you use water as negative control
- control_sample_nrt = 'noRT' if you use different naming

### 3. Reference Genes

Edit to match YOUR housekeeping genes:

    reference_gene = 'GAPDH'
    housekeeping_genes = ['ACTB']

Examples:
- reference_gene = '18S' for 18S ribosomal RNA
- reference_gene = 'HPRT1' for HPRT1
- housekeeping_genes = ['ACTB', 'B2M'] for multiple housekeeping genes

### 4. Column Names

Edit if YOUR CSV uses different column headers:

    column_name_sample = 'siRNA Sample'
    column_name_primers = 'Primers'
    column_name_cell_line = 'Cell Line'
    column_name_cq = 'Cq'

Examples:
- column_name_cq = 'Ct' if using Ct instead of Cq
- column_name_sample = 'Treatment' if samples called treatments
- column_name_primers = 'Target' if primers called targets

### 5. Which Genes/Cell Lines to Plot

Leave empty to auto-detect all:

    genes_to_plot = []
    cell_lines_to_plot = []

Or specify exactly which to plot:

    genes_to_plot = ['BRCA1', 'TP53', 'PTEN']
    cell_lines_to_plot = ['U2OS', 'HeLa']

### 6. Figure Appearance

Adjust visual settings:

    y_axis_max = 1.5              # Max y-axis value
    label_threshold = 0.2          # Label values below this
    fig_width_per_column = 0.8     # Width per column (inches)
    fig_height_per_row = 0.4       # Height per row (inches)
    horizontal_spacing = 0.10      # Space between columns
    vertical_spacing = 0.15        # Space between rows

Adjustments:
- To see values up to 2.0: y_axis_max = 2.0
- For wider plots: fig_width_per_column = 1.0
- For more compact: fig_width_per_column = 0.6

### 7. Color Palette

Customize colors for each gene (one color per row):

    colour_palette = ["#D4D3CF", "#DC6B83", "#75B1CE", ...]

The script cycles through these colors for each gene row.

## CSV Format Requirements

Your CSV needs these columns (names configurable):
- Sample identifier (e.g., NTC, Scrm, BRCA1)
- Primer/gene name (e.g., GAPDH, BRCA1)
- Cq value (numeric or "no cq")
- Cell line (optional)

Example CSV:

    Cell Line,siRNA Sample,Primers,Cq
    U2OS,NTC,GAPDH,no cq
    U2OS,Scrm,GAPDH,14.96
    U2OS,Scrm,GAPDH,14.98
    U2OS,Scrm,GAPDH,15.01
    U2OS,Scrm,BRCA1,22.85
    U2OS,Scrm,BRCA1,22.91
    U2OS,BRCA1,GAPDH,13.79
    U2OS,BRCA1,GAPDH,13.82
    U2OS,BRCA1,BRCA1,25.60
    U2OS,BRCA1,BRCA1,25.55

Key points:
- Technical replicates: Multiple rows for same condition are averaged
- Controls section: Include Scrm with reference gene for each cell line
- Per-gene measurements: Include Scrm + gene, gene siRNA + gene, both with reference gene

## Understanding Results

### Figure Output

The figure shows:
- Rows = genes
- Columns = cell lines
- Bars = Scrm (control at 1.0) and siRNA (knockdown efficiency)
- Small values (<0.2) are labeled with exact numbers

### CSV Output

The knockdown_efficiency.csv contains:

    Gene,Cell_Line,Scrm_Normalized,siRNA_Normalized,Remaining_%,Knockdown_%
    BRCA1,U2OS,1.0,0.07,7.0,93.0
    TP53,U2OS,1.0,0.12,12.0,88.0

Columns explained:
- Gene: Target gene
- Cell_Line: Cell line tested
- Scrm_Normalized: Scrambled control value (always 1.0)
- siRNA_Normalized: Knockdown sample value
- Remaining_%: Percentage of mRNA remaining
- Knockdown_%: Knockdown efficiency percentage

### Interpreting Values

    Bar Value | Remaining | Knockdown %
    1.0       | 100%      | 0%
    0.50      | 50%       | 50%
    0.20      | 20%       | 80%
    0.10      | 10%       | 90%
    0.05      | 5%        | 95%
    0.01      | 1%        | 99%

Example: Bar at 0.07 means 7% remaining = 93% knockdown efficiency

## Common Configuration Examples

### Example 1: Different Control Names

    control_sample_scrambled = 'siControl'
    control_sample_ntc = 'H2O'

### Example 2: Different Reference Gene

    reference_gene = '18S'
    housekeeping_genes = ['ACTB', 'GAPDH']

### Example 3: Custom Column Names

    column_name_sample = 'Treatment'
    column_name_primers = 'Target'
    column_name_cq = 'Ct'

### Example 4: Specific Genes Only

    genes_to_plot = ['BRCA1', 'BRCA2']
    cell_lines_to_plot = []

### Example 5: One Cell Line

    genes_to_plot = []
    cell_lines_to_plot = ['U2OS']

Output files:
- graph_siRNA_U2OS.pdf
- graph_siRNA_U2OS.png
- knockdown_efficiency_U2OS.csv

### Example 6: Custom Output Names

    output_name = "BRCA_screen_plate1"

Output files:
- BRCA_screen_plate1.pdf
- BRCA_screen_plate1.png
- BRCA_screen_plate1_efficiency.csv

## Troubleshooting

### "Required columns not found"

Problem: Script can't find necessary columns
Solution: Update column name settings to match your CSV:

    column_name_sample = 'YOUR_COLUMN_NAME'
    column_name_cq = 'YOUR_CQ_COLUMN'

The error message will show which columns are available.

### No data for a gene

Problem: Missing measurements
Solution: Ensure you have:
- Scrm with gene primers
- Gene siRNA with gene primers
- Both with reference gene primers

### Wrong control name in figure

Problem: Control sample name doesn't match
Solution: Update control name:

    control_sample_scrambled = 'YOUR_CONTROL_NAME'

### CSV output shows NaN values

Problem: Missing data for that gene/cell line combination
Solution: Check that you have complete measurements for that combination

### Figures saving to wrong location

Problem: Can't find output files
Solution: Check path settings and ensure folders are created:

    path_figures = Path('figures')

The script automatically creates folders if they don't exist.

## Tips for Best Results

1. Use consistent naming across all plates
2. Include 3 technical replicates for each condition
3. Test with one plate first to verify settings
4. Check console output for detected genes and cell lines
5. Ensure every sample has reference gene measurements
6. Use the CSV output to document knockdown efficiencies for methods sections
7. Keep a copy of your configuration settings with your data

## Output Files

### Figures (in figures/ folder)
- **PDF format**: High-quality vector graphics (recommended for publications)
- **PNG format**: Raster image at 600 DPI (good for presentations)

### CSV (in results/ folder)
- **knockdown_efficiency.csv**: Table with all knockdown percentages

All filenames include the output_name you specified, plus any cell line filters applied.

## Version History

- v1.2: Added CSV export of knockdown efficiencies
- v1.1: Fully configurable parameters for any experimental setup
- v1.0: Initial release

---

Last Updated: 7th October 2025