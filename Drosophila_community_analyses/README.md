# Drosophila Community Analyses - User Guide

## Data Availability

Data are available on the European Nucleotide Archive (ENA) under the study accession number **PRJEB112398**. 

The dataset contains:
- **Pooled samples**: Sequenced for 16S and ITS markers
- **Single-fly samples**: Sequenced for COI (Cytochrome Oxidase I)

## Getting Started

### Step 1: Download the Data and Scripts

1. Download the folder `Drosophila_community_analyses` from GitHub
2. Download the raw sequencing data from ENA (PRJEB112398)
3. Add the sequencing data files to this folder

**Note:** Pooled samples can remain in the directory without interfering with the scripts.

## COI Analysis Pipeline

### Step 2: Run Bioinformatic Processing

Execute the bash scripts contained in the `COI_Bioinformatic` folder to preprocess the raw sequencing data.

```bash
# Run the bioinformatic pipeline scripts
bash COI_Bioinformatic/script_name.sh
```

### Step 3: Prepare Data for Taxonomic Affiliation

Run the R Markdown script to prepare and format the COI data:

```R
# In R/RStudio, run:
rmarkdown::render("Creating_excel_COI_reads.Rmd")
```

**Output:** Creates formatted Excel files and FASTA sequences ready for BLAST analysis.

### Step 4: Perform Taxonomic Affiliation

Run the BLAST analysis script to assign taxonomy to COI sequences:

```R
# In R/RStudio, run:
rmarkdown::render("COI_Blast.Rmd")
```

**Output:** Produces taxonomic assignments and identifies potential contaminations in pooled samples.

### Step 5: Statistical Analyses and Visualization

Run the final R script to generate statistical analyses, figures, and publication-ready results:

```R
# In R/RStudio, run:
source("Typing_Drosophila_results.R")
```

**Outputs include:**
- Graphical abstracts (bar plots and distribution maps)
- Species composition by locality
- Statistical tests (Chi-squared, Fisher's exact test)
- Distribution maps of Drosophila species in French Guiana

## Required R Packages

Before running the scripts, ensure the following packages are installed:

```R
# Bioinformatics packages
install.packages("BiocManager")
BiocManager::install("Biostrings")

# Data manipulation and visualization
install.packages(c(
  "readxl", "writexl", "dplyr", "tidyverse", "tidyr", "purrr",
  "ggplot2", "ggpubr", "rstatix", "gtsummary", "cowplot", "patchwork"
))

# Spatial visualization
install.packages(c("sf", "rnaturalearth", "rnaturalearthdata", "scatterpie"))

# Statistical analysis
install.packages(c("lme4", "lmerTest", "lmtest", "car", "emmeans", "entropy", "questionr"))

# Phylogenetic analysis
BiocManager::install("ape")
BiocManager::install("msa")
BiocManager::install("phangorn")
```

## Script Descriptions

### Creating_excel_COI_reads.Rmd
Converts raw FASTA consensus sequences into structured Excel tables for downstream analysis.

### COI_Blast.Rmd
Performs BLAST searches against reference sequences and assigns taxonomic identifications to Drosophila species.

### Typing_Drosophila_results.R
Generates comprehensive statistical analyses, publication-quality figures, and species distribution maps.

## Output Files

The pipeline generates the following main output files:

- `COI_consensus_sequence_table.xlsx` - Consensus sequences with metadata
- `Blast_results_filtered.xlsx` - BLAST results with taxonomic assignments
- `Miss_identify_sample.xlsx` - Potential contaminations detected
- `Identification_results_vf.csv` - Final species identifications
- `Graphical_abstract_category_barplots_Cayenne_vs_other_localities.pdf` - Summary figure
- `camemberts_par_localite.pdf` - Distribution maps with pie charts
- `barplots_with_legend_below_map.pdf` - Distribution maps with bar charts

## Troubleshooting

- **BLAST errors:** Ensure all input FASTA files are properly formatted
- **Missing packages:** Run the installation commands above before executing scripts
- **File path issues:** Set working directory to the folder containing the scripts using `setwd()`

## Citation

If you use these scripts and data, please cite the ENA study: **PRJEB112398**

## Contact

For questions about the analysis or data, please contact the repository maintainers.
