# RT-qPCR siRNA Validation Analysis

Python pipeline for analyzing RT-qPCR data from siRNA knockdown experiments. Includes tools for importing Agilent Aria data and generating publication-quality analysis figures.

## Overview

This pipeline provides two complementary workflows:

1. **Data Import** (`format_aria_data.py`): Import Agilent Aria qPCR data into a standardized platemap format
2. **Data Analysis** (`analysis.py`): Analyze knockdown efficiency and generate figures

### What it does:
- Import raw Agilent Aria qPCR data and map to plate layout
- Generate visual plate layout figures for quality control
- Calculate relative gene expression normalized to housekeeping genes
- Normalize knockdown samples to scrambled siRNA controls
- Generate publication-quality multi-panel figures
- Export knockdown efficiency data to CSV
- Automatically average technical replicates
- Handle multiple cell lines and genes in a single run

## Requirements

Python Packages (install with pip):
- pandas
- numpy
- matplotlib
- openpyxl (for Excel file support)

Python version: 3.7+

## Quick Start

### Option A: Full Pipeline (Aria Import + Analysis)

1. Fill in `qPCR_results_platemap.xlsx` with your plate layout
2. Export data from Agilent Aria (Excel format)
3. Run: `python format_aria_data.py`
4. Run: `python analysis.py`
5. Find results in `results/` and `figures/` folders

### Option B: Analysis Only (Already have CSV data)

1. Place your formatted CSV file(s) in the `data/` folder
2. Edit configuration in `analysis.py`
3. Run: `python analysis.py`
4. Find results in `results/` and `figures/` folders

## Project Structure

    project/
    ├── format_aria_data.py        # Import Aria data
    ├── platemap_figures.py        # Generate plate layout visualizations
    ├── analysis.py                # Main analysis script
    ├── qPCR_results_platemap.xlsx # Platemap template
    ├── agilent-data/              # Place Aria exports here
    │   └── your_aria_export.xlsx
    ├── data/                      # Formatted CSV files
    │   └── experiment1.csv
    ├── plate_layouts/             # QC plate layout figures
    │   ├── plate_layout_samples.png
    │   └── plate_layout_primers.png
    ├── results/                   # CSV output files
    │   └── knockdown_efficiency.csv
    └── figures/                   # Final analysis figures
        ├── graph_siRNA.pdf
        └── graph_siRNA.png

Folders are automatically created if they don't exist.

## Part 1: Data Import from Agilent Aria

### The Platemap Template

The `qPCR_results_platemap.xlsx` file is your experimental design template. Fill in these columns:

- **Gene of Interest**: Category/grouping (e.g., "Controls", "BRCA1", "TP53")
- **Cell Line**: Cell line name (e.g., "U2OS", "HeLa")
- **siRNA Sample**: Sample type (e.g., "NTC", "NRT", "Scrm", or gene name)
- **Primers**: Which primers used (e.g., "GAPDH", "ACTB", or gene name)
- **Well**: Well position (e.g., "A1", "B2")
- **Cq**: Leave blank - auto-filled from Aria
- **Melting Point**: Leave blank - auto-filled from Aria

**Flexible Layout**: You can organize your plate however you want:
- Controls (NTC, NRT, Scrm) with housekeeping genes (GAPDH, ACTB)
- Normalization samples (Scrm with each gene's primers)
- Gene knockdowns with their own primers
- Gene knockdowns with housekeeping gene primers

The import script will map Aria Cq values to the correct wells regardless of your layout.

### Configuring the Aria Import

Edit the top of `format_aria_data.py`:

    # Input files
    aria_export_file = 'agilent-data/your_export.xlsx'
    platemap_template = 'qPCR_results_platemap.xlsx'
    
    # Output
    output_dir = Path('data')
    output_file = 'experiment1.csv'
    
    # Aria sheet name (usually "Tabular Results")
    aria_sheet_name = 'Tabular Results'
    
    # Column names in Aria export
    aria_cq_column = 'Cq (∆R)'
    aria_tm_column = "Tm Product 1 (-R'(T))"

### Running the Import

    python format_aria_data.py

This will:
1. Read your Aria Excel export
2. Match wells between Aria data and your platemap
3. Fill in Cq and Melting Point values
4. Save formatted CSV to `data/` folder
5. Generate plate layout figures in `plate_layouts/` for QC

**Quality Control**: Check the generated plate layout figures to verify:
- Samples are distributed correctly
- Primers are mapped to the right wells
- No missing or misplaced data

## Part 2: Analysis Configuration

All settings are in the USER CONFIGURATION section at the top of `analysis.py`.

### 1. Paths and Output Name

    path_data = Path('data')
    path_results = Path('results')
    path_figures = Path('figures')
    output_name = "sample_output"

- `path_data`: Folder with formatted CSV files
- `path_results`: Where CSV results are saved
- `path_figures`: Where figures are saved
- `output_name`: Base name for all output files

### 2. Control Sample Names

**Edit these to match your platemap:**

    control_sample_ntc = 'NTC'
    control_sample_nrt = 'NRT'
    control_sample_scrambled = 'Scrm'

Examples:
- `control_sample_scrambled = 'siControl'` if you use siControl
- `control_sample_ntc = 'H2O'` if you use water
- `control_sample_nrt = 'noRT'` for different naming

### 3. Reference Genes

**Edit to match your housekeeping genes:**

    reference_gene = 'GAPDH'
    housekeeping_genes = ['ACTB']

Examples:
- `reference_gene = '18S'` for 18S ribosomal RNA
- `housekeeping_genes = ['ACTB', 'B2M']` for multiple housekeeping genes

### 4. Column Names

**Edit if your CSV uses different headers:**

    column_name_sample = 'siRNA Sample'
    column_name_primers = 'Primers'
    column_name_cell_line = 'Cell Line'
    column_name_cq = 'Cq'

Examples:
- `column_name_cq = 'Ct'` if using Ct instead of Cq
- `column_name_sample = 'Treatment'` for different naming

### 5. Which Genes/Cell Lines to Plot

Leave empty to auto-detect:

    genes_to_plot = []
    cell_lines_to_plot = []

Or specify exactly which to plot:

    genes_to_plot = ['BRCA1', 'TP53', 'PTEN']
    cell_lines_to_plot = ['U2OS', 'HeLa']

### 6. Figure Appearance

    y_axis_max = 1.5
    label_threshold = 0.2
    fig_width_per_column = 0.8
    fig_height_per_row = 0.4
    horizontal_spacing = 0.10
    vertical_spacing = 0.15

### 7. Color Palette

    colour_palette = ["#D4D3CF", "#DC6B83", "#75B1CE", ...]

One color per gene row.

## Data Format Requirements

### For Direct CSV Input (Skip Aria Import)

Your CSV needs these columns:
- **siRNA Sample**: Sample identifier
- **Primers**: Primer/gene name
- **Cq**: Cq value (numeric or "no cq")
- **Cell Line**: (Optional) Cell line name

Example:

    Cell Line,siRNA Sample,Primers,Cq
    U2OS,NTC,GAPDH,no cq
    U2OS,Scrm,GAPDH,14.96
    U2OS,Scrm,GAPDH,14.98
    U2OS,Scrm,BRCA1,22.85
    U2OS,BRCA1,GAPDH,13.79
    U2OS,BRCA1,BRCA1,25.60

**Key requirements:**
- 3 technical replicates (script averages them)
- Scrm with reference gene for each cell line
- For each gene: Scrm + gene primers, gene siRNA + gene primers, gene siRNA + reference gene primers

## Understanding Results

### Plate Layout Figures (QC)

Generated during import in `plate_layouts/`:
- **plate_layout_samples.png**: Color-coded by sample type
- **plate_layout_primers.png**: Color-coded by primer type

Use these to verify your plate setup before analysis.

### Analysis Figure

Shows knockdown efficiency:
- **Rows**: One per gene
- **Columns**: One per cell line
- **Bars**: Scrm (1.0) and siRNA (knockdown value)
- **Labels**: Values <0.2 show exact numbers

### CSV Output

`knockdown_efficiency.csv` contains:

    Gene,Cell_Line,Scrm_Normalized,siRNA_Normalized,Remaining_%,Knockdown_%
    BRCA1,U2OS,1.0,0.07,7.0,93.0
    TP53,U2OS,1.0,0.12,12.0,88.0

### Interpreting Values

    Bar Value | Remaining | Knockdown %
    1.0       | 100%      | 0%
    0.50      | 50%       | 50%
    0.20      | 20%       | 80%
    0.10      | 10%       | 90%
    0.05      | 5%        | 95%
    0.01      | 1%        | 99%

## Complete Workflow Example

### Step 1: Design Your Plate

Fill `qPCR_results_platemap.xlsx`:
- A1-A3: Controls (NTC) with GAPDH
- B1-B3: Controls (Scrm) with GAPDH
- C1-C3: Controls (Scrm) with ACTB
- D1-D3: BRCA1 (Scrm) with BRCA1 primers
- E1-E3: BRCA1 (knockdown) with GAPDH
- F1-F3: BRCA1 (knockdown) with BRCA1 primers
- etc.

### Step 2: Run qPCR

Run your plate on Agilent Aria, export as Excel.

### Step 3: Import Data

Edit `format_aria_data.py`:
    
    aria_export_file = 'agilent-data/myplate.xlsx'
    output_file = 'myplate.csv'

Run:

    python format_aria_data.py

Check `plate_layouts/` figures for QC.

### Step 4: Analyze

Edit `analysis.py`:

    output_name = "myplate_analysis"

Run:

    python analysis.py

Find results in `results/` and `figures/`.

## Common Use Cases

### Use Case 1: Full Pipeline (Aria Import)

    # 1. Fill platemap template
    # 2. Export from Aria
    python format_aria_data.py
    # 3. Check plate_layouts/ for QC
    python analysis.py

### Use Case 2: Analysis Only (Have CSV)

    # 1. Put CSV in data/
    python analysis.py

### Use Case 3: Batch Multiple Plates

    # Import all plates
    python format_aria_data.py  # plate1
    python format_aria_data.py  # plate2
    python format_aria_data.py  # plate3
    
    # Analyze all together
    python analysis.py

The analysis script automatically processes all CSVs in `data/`.

### Use Case 4: Different Experimental Setups

The platemap is flexible - organize however you want:
- Different cell lines across columns
- Multiple genes
- Different siRNA concentrations
- Time course experiments

Just fill the platemap template with your layout, and the scripts handle the rest.

## Troubleshooting

### Aria Import Issues

**"Column not found in Aria export"**
- Check the exact column names in your Aria Excel file
- Update `aria_cq_column` and `aria_tm_column` in config

**"Well not found"**
- Verify well names match exactly (e.g., "A1" not "a1")
- Check for extra spaces in well names

**Plate layout figures look wrong**
- Review your platemap Excel file
- Ensure Well, siRNA Sample, and Primers columns are filled correctly

### Analysis Issues

**"Required columns not found"**
- Update `column_name_*` settings to match your CSV headers

**No data for a gene**
- Ensure Scrm measured with gene primers
- Check gene siRNA measured with both gene and reference primers

**Wrong control name in figure**
- Update `control_sample_scrambled` to match your data

### General Tips

1. **Use example data**: Run with `example_data.csv` first to verify setup
2. **Check QC figures**: Always review plate layouts after import
3. **Consistent naming**: Use same control/gene names across all plates
4. **Technical replicates**: Include 3 replicates for reliable results
5. **Console output**: Read messages for detected genes/cell lines

## Output Files Summary

### From Aria Import:
- `data/*.csv`: Formatted data ready for analysis
- `plate_layouts/plate_layout_samples.png`: Sample distribution QC
- `plate_layouts/plate_layout_primers.png`: Primer distribution QC

### From Analysis:
- `results/*_efficiency.csv`: Knockdown percentages table
- `figures/*.pdf`: Vector graphics figure (publication)
- `figures/*.png`: Raster image (600 DPI, presentations)

## Version History

- v1.3: Added Aria data import pipeline and plate layout visualization
- v1.2: Added CSV export of knockdown efficiencies
- v1.1: Fully configurable parameters
- v1.0: Initial release

## Citation

If you use this pipeline in your research, please acknowledge the tool and cite relevant publications from your lab.

---

Last Updated: 7th October 2025
