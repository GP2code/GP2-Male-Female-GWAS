# Shared and Sex-Specific Genetic Risk for Parkinson's Disease Risk Across European Populations

`GP2 ❤️ Open Science 😍`

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Last Updated:** December 2025

## Summary

This is the online repository for the manuscript titled **"Shared and Sex-Specific Genetic Risk for Parkinson's Disease Risk Across European Populations"**. This study presents the largest sex-stratified genome-wide association study (GWAS) meta-analysis of Parkinson's disease (PD) risk to date, including 226,196 individuals (18,145 female PD cases, 95,558 female controls; 28,747 male PD cases, and 83,746 male controls) from European populations. We identified 57 genome-wide significant association signals, including five novel risk loci with sex-specific effects, and provide functional annotation and colocalization analyses.

## Citation

*(pending publication)*


## Data Statement

* All GP2 data are hosted in collaboration with the Accelerating Medicines Partnership in Parkinson's disease and are available via application at https://amp-pd.org/register-for-amp-pd. For up-to-date information on GP2 data acquisition, access, and policies, visit https://gp2.org/.
    * This study utilized GP2 releases 9 ([10.5281/zenodo.14510099](https://doi.org/10.5281/zenodo.14510099)) and 10 ([10.5281/zenodo.15748014](https://doi.org/10.5281/zenodo.15748014)). Genotyping imputation, quality control, ancestry prediction, and processing were performed using GenoTools (v1.0.0), publicly available on GitHub.
* UK Biobank data were accessed under Application Number 33601 (https://www.ukbiobank.ac.uk/use-our-data/apply-for-access/).
* Fox Insight data are available at https://foxinsight.michaeljfox.org/register.

### Helpful Links

- [GP2 Website](https://gp2.org/)
  - [GP2 Cohort Dashboard](https://gp2.org/cohort-dashboard-advanced/)
- [Introduction to GP2](https://movementdisorders.onlinelibrary.wiley.com/doi/10.1002/mds.28494)
  - [Other GP2 Manuscripts (PubMed)](https://pubmed.ncbi.nlm.nih.gov/?term=%22global+parkinson%27s+genetics+program%22)


## Repository Orientation

The `analyses/` directory includes analyses discussed in the manuscript, organized by analysis type:

```
THIS_REPO/
├── LICENSE
├── README.md
└── analyses/
    ├── GWAS/
    │   └── GWAS.ipynb
    ├── beta-beta/
    │   ├── betabeta_MF.ipynb
    │   └── gp2_vars_v2.txt
    ├── LDSC/
    │   └── male_female_LDSC.ipynb
    ├── forest_plots/
    │   ├── forest_plots.ipynb
    │   └── combined_forest_data.csv
    ├── locusCompare/
    │   └── locusCompare_setAxes_R.r
    └── coloc/
        ├── 01_define_loci.ipynb
        ├── 02_coloc_results.ipynb
        ├── coloc_GTEx_cisQTL_v10.v4.R
        ├── topSNPs.txt
        └── loci_500kb_topSNPs.tsv
```

## Analysis Notebooks

### Languages: Python, bash, and R

| Directory | Notebook/Script | Description |
|-----------|----------------|-------------|
| `GWAS/` | `GWAS.ipynb` | Sex-stratified data prep and GWAS meta-analyses across GP2 as an example |
| `beta-beta/` | `betabeta_MF.ipynb` | Beta-beta plots comparing effect sizes between male and female meta-analyses using 157 variants from GP2 et al. 2025; correlation analysis of sex-specific effects |
|  | `gp2_vars_v2.txt` | List of 157 PD risk variants from GP2 et al. 2025 used for beta-beta comparison |
| `LDSC/` | `male_female_LDSC.ipynb` | Linkage disequilibrium score regression (LDSC) analysis |
| `forest_plots/` | `forest_plots.ipynb` | Forest plot generation |
| | `combined_forest_data.csv` | Summary statistics for forest plot generation including beta estimates, standard errors, and p-values for both sexes |
| `locusCompare/` | `locusCompare_setAxes_R.r` | R script for creating locusCompare plots comparing GWAS signals between males and females at specific loci  |
| `coloc/` | `01_define_loci.ipynb` | Definition of genomic loci for colocalization analysis; extraction of lead SNPs and surrounding regions from sex-stratified GWAS |
|  | `02_coloc_results.ipynb` | Processing and filtering of colocalization results; identification of sex-specific and sex-shared eQTL colocalization signals |
|  | `coloc_GTEx_cisQTL_v10.v4.R` | Bayesian colocalization analysis using COLOC method with GTEx v10 brain eQTL data across 13 brain tissues; calculation of posterior probabilities for shared causal variants |
|  | `topSNPs.txt` | List of lead SNPs for colocalization analysis |
|  | `loci_500kb_topSNPs.tsv` | Genomic coordinates (±500kb) around lead SNPs for colocalization analysis |


## Software

| **Software** | **Version(s)** | **Resource URL** | **RRID** | **Notes** |
|--------------|----------------|------------------|----------|-----------|
| PLINK | 1.9 and 2.0 | http://www.cog-genomics.org/plink/ | RRID:SCR_001757 | Used for GWAS and meta-analyses |
| LDSC | 1.0.1 | https://github.com/bulik/ldsc | RRID:SCR_017248 | Used for genetic correlation and heritability estimation |
| FUMA | v1.5.3 | https://fuma.ctglab.nl/ | RRID:SCR_017346 | Used for functional annotation and LD clumping |
| COLOC | 5.2.3 | https://github.com/chr1swallace/coloc | - | Used for colocalization analysis with GTEx eQTLs |
| GenoTools | 1.0.0 | https://github.com/GP2code/GenoTools | - | Used for QC, ancestry prediction, and imputation |
| KING | 2.2.5 | https://www.kingrelatedness.com/ | RRID:SCR_018787 | Used for relatedness estimation |
| Python | 3.10.15 | http://www.python.org/ | RRID:SCR_008394 | pandas, numpy, seaborn, matplotlib, scipy; used for data wrangling and visualization |
| R | 4.3.3 | http://www.r-project.org/ | RRID:SCR_001905 | forestplot, coloc, ggplot2; used for statistical analysis and visualization |

## Contact

For questions about this analysis, please contact:
- Data and Working Group Committee (dawg@gp2.org)
