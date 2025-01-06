
# SEQ4EPI: ClusterFinder Pipeline for Viral Genomic Analysis

The **ClusterFinder** pipeline is a comprehensive solution for viral genomic data analysis, integrating sequence alignment, distance calculation, clustering, and result table generation. It supports the analysis of **SARS-CoV-2**, **RSV-A**, **RSV-B**, and **Influenza** viruses (H1N1, H3N2, Type B). This pipeline combines tools written in Python, R, and Shell scripts to provide a detailed understanding of viral evolutionary dynamics.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Details](#pipeline-details)
- [Script Descriptions](#script-descriptions)
- [Examples](#examples)
- [Contact](#contact)

## Prerequisites

Ensure the following tools and packages are installed on your system:

### Python 3.x
- Required packages: `pandas`, `argparse`, `logging`

### R
- Required packages: `tidyverse`, `ggplot2`

### Other Tools
- **IQ-TREE**: For phylogenetic tree construction.
- **MAFFT**: For sequence alignment.
- **Nextalign**: For sequence alignment (alternative to MAFFT).
- **Nextclade**: For sequence analysis (Nextclade datasets must be pre-installed).
- **TreeCluster**: For clustering sequences based on phylogenetic trees (clustering thresholds and methods are customizable).

## Installation

Clone the repository to your local machine:

```bash
git clone https://github.com/yourusername/ClusterFinder.git
cd ClusterFinder
```

## Usage

The pipeline can be executed depending on the virus type using the provided shell scripts:

- For **SARS-CoV-2**: `ClusterFinder_NCOV.sh`
- For **Influenza**: `ClusterFinder_Flu.sh`
- For **RSV**: `ClusterFinder_RSV.sh`

## Pipeline Details

The **ClusterFinder** pipeline follows these steps for viral genomic analysis:

1. **Sequence Alignment**: Sequences are aligned using **Nextalign** or **MAFFT** with a reference genome.
2. **Sequence Analysis**: Sequence analysis is conducted with **Nextclade** to identify mutations and variants.
3. **Phylogenetic Tree Construction**: The phylogenetic tree is constructed using **IQ-TREE** with the best-fit substitution model and bootstrap resampling.
4. **Sequence Clustering**: **TreeCluster** groups sequences based on the phylogenetic tree, applying clustering thresholds tailored to each virus’s evolution.
5. **Metadata Integration**: Clusters are combined with metadata (e.g., date, location) for a comprehensive analysis.
6. **Final Table Generation**: Results are combined in a final table using **R scripts**.
7. **Visualization**: The data is prepared for visualization with **iTOL**.

## Script Descriptions

### Python Scripts

- **ClusterFinder.py**: Main script for sequence alignment, analysis, and clustering.

  **Usage**:
  ```bash
  ClusterFinder.py [-h] --input_fasta INPUT_FASTA --output_prefix OUTPUT_PREFIX --genome_length GENOME_LENGTH --threshold THRESHOLD --reference REFERENCE --annotation ANNOTATION --metadata_file METADATA_FILE
  ```

  - `--input_fasta`: Path to the input FASTA file containing unaligned sequences.
  - `--output_prefix`: Prefix for the output files.
  - `--genome_length`: Length of the reference genome.
  - `--threshold`: Threshold for clustering (adjustable based on analysis goals).
  - `--reference`: Reference genome file for **Nextalign**.
  - `--annotation`: GFF3 annotation file for **Nextalign**.
  - `--metadata_file`: Metadata CSV file.

### R Scripts

- **final_table_of_clusters.R**: Generates the final table of clusters after merging metadata.

  **Usage**:
  ```bash
  Rscript final_table_of_clusters.R input_metadata.csv output_final_table.csv
  ```

  - `input_metadata.csv`: Path to the merged metadata CSV file.
  - `output_final_table.csv`: Path to the output final table CSV file.

- **Itol.R**: Prepares the visualization data for **iTOL**.

  **Usage**:
  ```bash
  Rscript Itol.R input_treefile.tree output_final_table.csv output_visualization_folder
  ```

  - `input_treefile.tree`: Path to the input tree file in Newick format.
  - `output_final_table.csv`: Path to the output final table CSV file.
  - `output_visualization_folder`: Path to the folder for storing visualization files.

### Shell Scripts

- **ClusterFinder_NCOV.sh**: For **SARS-CoV-2** clustering.

  **Usage**:
  ```bash
  ./ClusterFinder_NCOV.sh input_fasta prefix genome_length threshold reference gff3 metadata final_table treefile output_visualisation
  ```

  - `input_fasta`: Path to the input FASTA file.
  - `prefix`: Prefix for the output files.
  - `genome_length`: Length of the reference genome.
  - `threshold`: Threshold for clustering (modifiable).
  - `reference`: Path to the reference genome FASTA file.
  - `gff3`: Path to the GFF3 annotation file.
  - `metadata`: Path to the metadata CSV file.
  - `final_table`: Path to the final table CSV file.
  - `treefile`: Path to the output tree file in Newick format.
  - `output_visualisation`: Path to the folder for storing visualization files.

- **ClusterFinder_Flu.sh**: For **Influenza** analysis (HA and NA fragments).

  **Usage**:
  ```bash
  bash ClusterFinder_Flu.sh       /path/to/input_HA-fragment.fasta       /path/to/input_NA-fragment.fasta       /path/to/output_prefix       1700 1400       0.00035 0.00025       /path/to/nextclade_dataset_HA-fragment/       /path/to/metadata_HA-fragment.csv       /path/to/nextclade_dataset_NA-fragment/       /path/to/metadata_NA-fragment.csv       /path/to/output_HA-fragment.treefile       /path/to/output_NA-fragment.treefile       /path/to/visualisation_HA-fragment/       /path/to/visualisation_NA-fragment
  ```

  - `input_HA-fragment.fasta` and `input_NA-fragment.fasta`: Paths to HA and NA FASTA files.
  - `output_prefix`: Prefix for the output files.
  - `1700`, `1400`: Genome lengths for HA and NA.
  - `0.00035`, `0.00025`: Clustering thresholds for HA and NA (modifiable).
  - Paths to pre-installed **Nextclade** datasets and metadata files for HA and NA.
  - Output paths for tree files and visualization folders.

- **ClusterFinder_RSV.sh**: For **RSV** clustering.

  **Usage**:
  ```bash
  bash ClusterFinder_RSV.sh       /path/to/input.fasta       /path/to/output_prefix       15400       0.00007       /path/to/reference.fasta       /path/to/genome_annotation.gff3       /path/to/metadata.csv       /path/to/output.treefile       /path/to/output_visualisation_folder
  ```

  - `input.fasta`: Path to the input FASTA file.
  - `output_prefix`: Prefix for the output files.
  - `15400`: Genome length for RSV.
  - `0.00007`: Clustering threshold (modifiable).
  - Paths to reference genome, annotation file, and metadata.
  - Output paths for tree files and visualization folder.

## Examples

### Running SARS-CoV-2 Pipeline

```bash
./ClusterFinder_NCOV.sh   /path/to/input.fasta   /path/to/output_prefix   29903   0.000085   /path/to/reference.fasta   /path/to/genome_annotation.gff3   /path/to/metadata.csv   /path/to/final_table_merged_metadata.csv   /path/to/output.treefile   /path/to/output_visualisation_folder
```

### Running Influenza Pipeline

```bash
bash ClusterFinder_Flu.sh     /path/to/input_HA-fragment.fasta     /path/to/input_NA-fragment.fasta     /path/to/output_prefix     1700 1400     0.00035 0.00025     /path/to/nextclade_dataset_HA-fragment/     /path/to/metadata_HA-fragment.csv     /path/to/nextclade_dataset_NA-fragment/     /path/to/metadata_NA-fragment.csv     /path/to/output_HA-fragment.treefile     /path/to/output_NA-fragment.treefile     /path/to/visualisation_HA-fragment/     /path/to/visualisation_NA-fragment
```

### Running RSV Pipeline

```bash
bash ClusterFinder_RSV.sh     /path/to/input.fasta     /path/to/output_prefix     15400     0.00007     /path/to/reference.fasta     /path/to/genome_annotation.gff3     /path/to/metadata.csv     /path/to/output.treefile     /path/to/output_visualisation_folder
```

## Contact

For questions or support, please contact [stephanie.dan@chu-lyon.fr].
