# Microbes in Tibetan permafrost

> This project repository associated with the following manuscript:

* Luyao Kang, Yutong Song, Rachel Mackelprang, Dianye Zhang, Shuqi Qin, Leiyi Chen, Linwei Wu, Yunfeng Peng and Yuanhe Yang*. Metagenomic insights into microbial structure and metabolism in alpine permafrost on the Tibetan Plateau.


## Structure

```
.
├── data
    ├── uploading later...
    └── 
└── script
    ├── read_data.R                             <-- codes for data input
    ├── amplicon_analysis.Rmd                   <-- microbial diversity and composition analysis based on the amplicon data
    ├── SPEC-OCCU.Rmd                           <-- specificity-occupancy analysis
    ├── env_effect_on_diver_comp.Rmd            <-- test the environmentaL effects on the microbial community
    ├── iCAMP_analysis.R                        <-- null model analysis using iCAMP framework
    ├── selected_ko_analysis.Rmd                <-- codes for the difference analysis in genes involved in C, N, S, and other elements among soil layers
    └── Metagenome_analysis.Rmd                 <-- codes for reproducing the results about the MAGs

```
All `.html` files in `script` folder are produced by Rmarkdown files (`.Rmd` files) using knitr package in R.

[![DOI](https://zenodo.org/badge/642855012.svg)](https://zenodo.org/doi/10.5281/zenodo.11498006)
