# Microbes in Tibetan permafrost
[![DOI](https://zenodo.org/badge/642855012.svg)](https://zenodo.org/doi/10.5281/zenodo.11498006)

> This project repository associated with the following manuscript:

* Luyao Kang, Yutong Song, Rachel Mackelprang, Dianye Zhang, Shuqi Qin, Leiyi Chen, Linwei Wu, Yunfeng Peng and Yuanhe Yang*. Metagenomic insights into microbial community structure and metabolism in alpine permafrost on the Tibetan Plateau.


## Structure
```
.
├── data1
    ├── 16S
        ├── metadata.txt                                            <-- Metadata including the climatic, vegetable for samples
        ├── otus.nwk                                                <-- Phylogenetic tree for ASVs
        ├── otutab.txt                                              <-- ASV table
        ├── otutab_rare.txt                                         <-- Rarefy ASV table
        └── taxonomy.txt                                            <-- Taxonomic table
    ├── Functional genes
        ├── eggnog.KEGG_ko.raw.counts.txt                           <-- KEGG KO table with row counts
        └── eggnog.KEGG_ko.raw.tpm.txt                              <-- KEGG KO table with TPM tansformation
    └── Metabolic
        ├── Metabolic_Sankey_diagram                                <-- Resuts of metabolic weighted scores for the metagenomic assembly genomes 
            ├── Metabolic_Sankey_diagram_sur_input_rare.txt         <-- The contributions of microbial groups to individual biogeochemical processes in surface layer
            ├── Metabolic_Sankey_diagram_sub_input_rare.txt         <-- The contributions of microbial groups to individual biogeochemical processes in subsurface layer
            └── Metabolic_Sankey_diagram_pl_input_rare.txt          <-- The contributions of microbial groups to individual biogeochemical processes in permafrost layer 
        ├── bin_abundance_table.tab                                 <-- Abundance table for each bin
        ├── METABOLIC_result.xlsx                                   <-- Summary of the presence/absence of functional genes on each genome
        ├── MW_score_result.csv                                     <-- Resuts of metabolic weighted scores for the metagenomic assembly genomes
        └── tax.txt                                                 <-- Taxnomic identificcation for each bin
└── script
    ├── read_data.R                                                 <-- Codes for data input
    ├── amplicon_analysis.Rmd                                       <-- Microbial diversity and composition analysis based on the amplicon data
    ├── SPEC-OCCU.Rmd                                               <-- Specificity-occupancy analysis
    ├── env_effect_on_diver_comp.Rmd                                <-- Test the environmentaL effects on the microbial community
    ├── iCAMP_analysis.R                                            <-- null model analysis using iCAMP framework
    ├── selected_ko_analysis.Rmd                                    <-- Codes for the difference analysis in genes involved in C, N, S, and other elements among soil layers
    └── Metagenome_analysis.Rmd                                     <-- Codes for reproducing the results about the MAGs
```
All `.html` files in `script` folder are produced by Rmarkdown files (`.Rmd` files) using knitr package in R.
