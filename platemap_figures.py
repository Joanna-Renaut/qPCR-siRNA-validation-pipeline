"""
Generate Plate Layout Figures
Creates visual representations of sample and primer layouts from platemap
Author: Joanna Renaut
Date: 7th October 2025
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from pathlib import Path

def parse_well(well):
    """Convert well (e.g., 'A1') to row and column indices"""
    if pd.isna(well) or not well:
        return None, None
    row = ord(well[0].upper()) - ord('A')
    col = int(well[1:]) - 1
    return row, col

def get_unique_colors(items):
    """Generate a color map for unique items"""
    unique_items = sorted([item for item in items if pd.notna(item)])

    # Use a nice color palette
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_items)))
    color_map = {item: colors[i] for i, item in enumerate(unique_items)}
    color_map[np.nan] = (1, 1, 1, 1)  # White for empty

    return color_map

def create_plate_layout(platemap, column_name, title, output_file,
                       well_column='Well'):
    """
    Create a 96-well plate layout figure

    Parameters:
    -----------
    platemap : pd.DataFrame
        DataFrame containing plate data
    column_name : str
        Column to visualize (e.g., 'siRNA Sample', 'Primers')
    title : str
        Figure title
    output_file : str or Path
        Output filename
    well_column : str
        Name of the well column (default: 'Well')
    """

    fig, ax = plt.subplots(figsize=(14, 10))

    # 96-well plate: 8 rows (A-H) x 12 columns (1-12)
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    cols = range(1, 13)

    # Get unique items and assign colors
    items = platemap[column_name].unique()
    color_map = get_unique_colors(items)

    # Create lookup dictionary from platemap
    well_data = {}
    for _, row in platemap.iterrows():
        well = row[well_column]
        if pd.notna(well):
            well_data[well] = row[column_name]

    # Draw wells
    for i, row_label in enumerate(rows):
        for j, col in enumerate(cols):
            well = f"{row_label}{col}"

            # Get data for this well
            item = well_data.get(well, np.nan)
            color = color_map.get(item, (1, 1, 1, 1))

            # Draw circle
            circle = patches.Circle((j, 7-i), 0.4,
                                    linewidth=1.5,
                                    edgecolor='black',
                                    facecolor=color)
            ax.add_patch(circle)

            # Add text label if there's data
            if pd.notna(item) and item != '':
                # Truncate long names
                label = str(item)[:8] if len(str(item)) > 8 else str(item)
                ax.text(j, 7-i, label,
                       ha='center', va='center',
                       fontsize=7, fontweight='bold')

    # Row labels (A-H)
    for i, row_label in enumerate(rows):
        ax.text(-0.8, 7-i, row_label,
               ha='center', va='center',
               fontsize=12, fontweight='bold')

    # Column labels (1-12)
    for j, col in enumerate(cols):
        ax.text(j, 7.8, str(col),
               ha='center', va='center',
               fontsize=12, fontweight='bold')

    ax.set_xlim(-1.2, 12)
    ax.set_ylim(-0.8, 8.5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)

    # Create legend
    legend_items = sorted([item for item in color_map.keys() if pd.notna(item)])
    legend_patches = [patches.Patch(color=color_map[item], label=item)
                     for item in legend_items]

    ax.legend(handles=legend_patches,
             loc='center left',
             bbox_to_anchor=(1.02, 0.5),
             fontsize=10,
             frameon=True)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file}")
    plt.close()

def generate_plate_figures(platemap_file, output_dir, experiment,
                          sample_column='siRNA Sample',
                          primer_column='Primers',
                          well_column='Well'):
    """
    Generate both sample and primer layout figures from a platemap

    Parameters:
    -----------
    platemap_file : str or Path
        Path to platemap CSV file
    output_dir : str or Path
        Directory to save output figures
    sample_column : str
        Column name for samples (default: 'siRNA Sample')
    primer_column : str
        Column name for primers (default: 'Primers')
    well_column : str
        Column name for wells (default: 'Well')

    Returns:
    --------
    tuple : (sample_layout_path, primer_layout_path)
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    # Read platemap
    platemap = pd.read_csv(platemap_file)

    # Generate sample layout
    sample_output = output_dir / f'{experiment}_plate_layout_samples.png'
    create_plate_layout(
        platemap,
        sample_column,
        '96-Well Plate Layout - Sample Distribution',
        sample_output,
        well_column=well_column
    )

    # Generate primer layout
    primer_output = output_dir / f'{experiment}_plate_layout_primers.png'
    create_plate_layout(
        platemap,
        primer_column,
        '96-Well Plate Layout - Primer Distribution',
        primer_output,
        well_column=well_column
    )

    return sample_output, primer_output