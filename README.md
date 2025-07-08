
# title: "Comorbidity Analysis of Systemic Sclerosis and Cancer using Bioinformatics"


# 📘 Overview

This repository contains the full bioinformatics workflow developed for the publication:

**"Network based systems biology approach to identify diseasome and comorbidity associations of Systemic Sclerosis with cancers"**  
📌 *Heliyon, Volume 8, Issue 4, 2022, e08892*  
🔗 [https://doi.org/10.1016/j.heliyon.2022.e08892](https://doi.org/10.1016/j.heliyon.2022.e08892)

Systemic Sclerosis (SSc) is a complex autoimmune disease increasingly associated with cancers such as **lung cancer**, **leukemia**, and **lymphoma**. To explore these comorbidities, we developed a modular R-based pipeline that performs gene expression analysis, enrichment testing, and network-based integration of multi-omics evidence.

This pipeline is generalizable and can be reused to explore comorbid relationships between any two diseases using transcriptomic datasets.

---

# 📁 Repository Structure

```bash
comorbidity/
├── data/                 # Input RNA-seq datasets from GEO
├── scripts/              # R scripts for each analysis module
├── results/              # Output DEG tables, GSEA results, pathway lists, similarity matrices
├── figures/              # Enrichment plots, DAGs, network diagrams
└── README.Rmd            # This file
