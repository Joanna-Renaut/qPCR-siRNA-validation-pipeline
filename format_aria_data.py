"""
Agilent Aria Data Formatter
Takes raw Aria export and fills in a platemap template based on well positions
Author: Joanna Renaut
Date: 7th October 2025
"""

import pandas as pd
from pathlib import Path
import sys
from platemap_figures import generate_plate_figures

# ============================================================================
# USER CONFIGURATION
# ============================================================================

# Input files
aria_export_file = 'path-to-your-agilent-aria-data.xlsx' # .xlsx export from aria agilent software
platemap_template = 'qPCR_results_platemap.xlsx'  # filled in platemap .xlsx

# Output dir
output_dir = Path('data')
figures_dir = Path('plate_layouts')

# Output file
output_file = 'example.csv'

# Aria Excel sheet name
aria_sheet_name = 'Tabular Results'

# Column mappings from Aria export
aria_well_column = 'Well'
aria_cq_column = 'Cq (∆R)'
aria_tm_column = "Tm Product 1 (-R'(T))"

# Column names in your platemap
platemap_well_column = 'Well'
platemap_cq_column = 'Cq'
platemap_tm_column = 'Melting Point'
platemap_sample_column = 'siRNA Sample'
platemap_primer_column = 'Primers'

# Generate plate layout figures
generate_figures = True

# ============================================================================
# END USER CONFIGURATION
# ============================================================================

def format_aria_data(aria_file, template_file, output_file):
    """
    Fill platemap template with Aria export data based on well positions
    """

    print("=" * 60)
    print("Agilent Aria Data Formatter")
    print("=" * 60)

    # Read Aria export (Excel file - specific sheet)
    try:
        aria_data = pd.read_excel(aria_file, sheet_name=aria_sheet_name)
        print(f"\n✓ Loaded Aria export: {aria_file}")
        print(f"  Sheet: {aria_sheet_name}")
        print(f"  Found {len(aria_data)} wells")
    except Exception as e:
        print(f"\n✗ Error reading Aria file: {e}")
        return

    # Read platemap template (Excel file)
    try:
        platemap = pd.read_excel(template_file)  # Changed to read_excel
        print(f"✓ Loaded platemap template: {template_file}")
        print(f"  Found {len(platemap)} rows")
    except Exception as e:
        print(f"\n✗ Error reading platemap template: {e}")
        return

    # Check if required columns exist in Aria data
    if aria_well_column not in aria_data.columns:
        print(f"\n✗ Error: Column '{aria_well_column}' not found in Aria export")
        print(f"  Available columns: {aria_data.columns.tolist()}")
        return

    aria_cq_found = aria_cq_column in aria_data.columns
    aria_tm_found = aria_tm_column in aria_data.columns

    if not aria_cq_found:
        print(f"\n⚠ Warning: Column '{aria_cq_column}' not found in Aria export")
    if not aria_tm_found:
        print(f"⚠ Warning: Column '{aria_tm_column}' not found in Aria export")

    # Create lookup dictionary from Aria data
    aria_lookup = {}
    for _, row in aria_data.iterrows():
        well = row[aria_well_column]
        aria_lookup[well] = {
            'Cq': row.get(aria_cq_column, None) if aria_cq_found else None,
            'Tm': row.get(aria_tm_column, None) if aria_tm_found else None
        }

    print(f"\n✓ Created lookup for {len(aria_lookup)} wells")

    # Fill in the platemap
    filled_cq_count = 0
    filled_tm_count = 0

    for idx, row in platemap.iterrows():
        well = row[platemap_well_column]

        if pd.notna(well) and well in aria_lookup:
            # Fill Cq value
            if platemap_cq_column in platemap.columns and aria_cq_found:
                cq_value = aria_lookup[well]['Cq']
                if pd.notna(cq_value):
                    platemap.at[idx, platemap_cq_column] = cq_value
                    filled_cq_count += 1

            # Fill Melting Point value
            if platemap_tm_column in platemap.columns and aria_tm_found:
                tm_value = aria_lookup[well]['Tm']
                if pd.notna(tm_value):
                    platemap.at[idx, platemap_tm_column] = tm_value
                    filled_tm_count += 1

    # Save filled platemap as CSV
    output_path = output_dir / output_file
    output_dir.mkdir(exist_ok=True)
    platemap.to_csv(output_path, index=False)

    print(f"\n✓ Filled {filled_cq_count} Cq values")
    print(f"✓ Filled {filled_tm_count} Melting Point values")
    print(f"✓ Saved to: {output_path}")

    # Generate plate layout figures
    if generate_figures:
        print(f"\nGenerating plate layout figures...")
        try:
            sample_fig, primer_fig = generate_plate_figures(
                output_path,  # Use the CSV output for figure generation
                figures_dir,
                sample_column=platemap_sample_column,
                primer_column=platemap_primer_column,
                well_column=platemap_well_column
            )
            print(f"✓ Plate layout figures saved to: {figures_dir}")
        except Exception as e:
            print(f"⚠ Warning: Could not generate figures: {e}")

    print(f"\n{'=' * 60}")
    print("Formatting complete!")
    print(f"{'=' * 60}\n")


if __name__ == "__main__":
    # Check if files are provided as command line arguments
    if len(sys.argv) > 1:
        aria_export_file = sys.argv[1]
    if len(sys.argv) > 2:
        platemap_template = sys.argv[2]
    if len(sys.argv) > 3:
        output_file = sys.argv[3]

    format_aria_data(aria_export_file, platemap_template, output_file)