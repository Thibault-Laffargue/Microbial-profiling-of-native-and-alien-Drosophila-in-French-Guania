# Drosophila microbiome Analyses - User Guide

## Data Availability

Data are available on the European Nucleotide Archive (ENA) under the study accession number **PRJEB112398**. 

The dataset contains:
- **Pooled samples**: Sequenced for 16S and ITS markers
- **Single-fly samples**: Sequenced for COI (Cytochrome Oxidase I)

## Getting Started - Download the Data and Scripts

1. Download the repository `Microbial-profiling-of-native-and-alien-Drosophila-in-French-Guania` from GitHub and open the folder `Microbiome_analyses`
2. Download the raw sequencing data from ENA accession number PRJEB112398. (Optionnal if you want to skip Bioinformatic part)
3. Add the sequencing data files to this folder. (Optionnal if you want to skip Bioinformatic part)

**Note:** COI samples can remain in the directory without interfering with the scripts.

## Bioinformatics

In a terminal run Sequence_moving_script. 

### Bacterial 16s rRNA

![Workflow](Pipeline_Logigram.png)

1. Open the folder `bacterial_16s_rRNA_metabarcoding`
2. In a terminal run `Pipeline_part1`
3. In Rstudio run Stats_N_reads.R
4. In a terminal run `Pipeline_part2` 


### Fungal ITS

1. Open the folder `Fungi_ITS_metabarcoding`
2. In a terminal run `Pipeline_part1`
3. In Rstudio run Stats_N_reads.R
4. In a terminal run `Pipeline_part2` 

**Note:** For next step you will need `16s_OTU_table_without_chim.tsv`, `16s_metadata.csv`, `16s_tax_table_without_chim.tsv`,`centroids.tre`, `OL4_OTU_table_without_chim.tsv`, `OL4_metadata.csv` and `OL4_tax_table_without_chim.tsv`. If you skipped Bioinformatic part you will find this files in the compressed archive `Bio_informatics_results_file.zip`.

## Statistical analyses.

If you skipped Bioinformatic part pleas unzip `Bio_informatics_results_file.zip`.

Next steps will be done on Rstudio.

1. Run `1-Control+Phyloseq.Rmd`, this script allow to prepare the phyloseqs object for analyses. 
2. RUN `2-Alpha_diversity.Rmd`, this script allow to calculate alpha diversity indexes and to statistically compare the diversity accross Locality, HostSpecies and Invasion statut.
3. Run `3-Ordination.Rmd`, this script allow to represent and compare statically the structure and the composition of the microbiome accross Locality, HostSpecies and Invasion statut.
4. Run `4-Investigating_on_OTUs(unknow_treatment_version).Rmd`, this script allow the detection of taxa exclusive too Invasive alien Drosophila or Native species. 
5. Run `4.1-Class_relative_abundance`, this script allow the description of bacteria and fungi class distribution accross Locality, HostSpecies and Invasion statut.
6. Run `4.2-Pool_investigation.Rmd`, this script clarify pool composition for the pool with a suspect class composition.
7. Run `5-Core_microbiome.Rmd`, this script allow to describe core microbiome compoistion and to identify taxa unevenly distributed between Invasive alien Drosophila and Native Drosophila.
8. Run `6-comparison.Rmd`, this script allow comparison of the data first with excluded data in order to verify if mistake was made and then to reference dataset. The reference dataset is build with Staubach et al. (2013), Wang et al. (2020), and Brown et al. (2023) dataset.

**Note:** The file `Miss_identify_sample.xlsx`is a results of COI analyses. 
