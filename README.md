# Isolated or Connected? Genetic Structure of the Arctic Fox ( *Vulpes lagopus* ) 🦊❄️

This repository contains code and figures for a population-genetic study of **Arctic fox** (*Vulpes lagopus*) populations from **Greenland, Canada and Siberia**.

The Arctic fox is a highly mobile, circumpolar species that relies on **sea ice** to travel long distances. Because sea ice is shrinking with climate change, patterns of **migration and gene flow** between populations are expected to change.

This project uses genome-wide SNP data from **47 individuals** to ask:

* Are Arctic fox populations **genetically isolated** or **well connected**?
* How much **gene flow** occurs between Greenland, Canada and Siberia?
* Are there signs of **bottlenecks** or reduced diversity in Greenlandic populations?

---

## 📂 Repository structure

```text
.
├─ data/
│  └─ README_data.md           # description of genotype files (raw data NOT included)
├─ scripts/
│  ├─ 01_pca.R                 # PCA of SNP genotypes
│  ├─ 02_admixture_evaladmix.R # ADMIXTURE + evalAdmix visualisation
│  ├─ 03_fst.R                 # pairwise FST between regions / clusters
│  └─ 04_nj_tree.R             # IBS matrix and neighbour-joining tree
├─ figures/
│  ├─ pca_pc1_pc2.png
│  ├─ pca_pc2_pc3.png
│  ├─ admixture_k5.png
│  ├─ admixture_k6.png
│  ├─ admixture_k7.png
│  ├─ fst_regions.png
│  ├─ fst_clusters.png
│  └─ nj_tree.png
└─ README.md
```

---

## 🧬 Data

* **47 Arctic fox individuals**
* **7 regions**:

  * Greenland: **Qanisartuut, Zackenberg, Scoresbysund, Kangerlussuaq**
  * Canada: **Bylot Island, Karrak Lake**
  * Siberia: **Taymyr**
* Genotypes: high-density **SNP data** in PLINK format

Raw genotype files are **not included** in this repository (size & privacy).
The analysis assumes PLINK files like:

* `AF.imputed.thin.bed`
* `AF.imputed.thin.bim`
* `AF.imputed.thin.fam`

---

## 🔧 Methods (overview)

### 1. Quality control

* Relationship checking using **KING robust kinship** (via PLINK 2).
* Checked for **first- and second-degree relatives**; none were found, so all 47 individuals were retained.
* Basic SNP filters (e.g. call rate, MAF) applied during pre-processing.

### 2. Population structure

**Principal Component Analysis (PCA)**

* PCA on SNP genotypes (PLINK v1.9).
* PC1–PC3 used to visualise separation between regions/countries.
* PCA plots show:

  * clear separation of **Greenland vs Canada/Siberia**
  * Greenlandic subpopulations forming distinct clusters.

**ADMIXTURE + evalAdmix**

* Model-based clustering with **ADMIXTURE** for multiple K values (e.g. K=5–7).
* **evalAdmix** used to assess model fit and residual correlations.
* Barplots and residual heatmaps reveal:

  * homogeneous ancestry blocks in some regions
  * individuals with **admixed ancestry** across Greenland, Canada and Siberia.

### 3. Genetic differentiation (FST)

* Pairwise **Weir & Cockerham FST** estimated in R.
* Analyses performed:

  * by **geographic region**
  * by **clusters** defined from ADMIXTURE.
* Heatmaps show:

  * **low FST between Canadian and Siberian populations** (high gene flow/shared ancestry)
  * **higher FST** for some Greenlandic regions, especially **Qanisartuut**.

### 5. Gene flow and clustering

**Neighbour-joining (NJ) tree**

* Genetic distance (e.g. IBS matrix) used to build a **NJ tree**.
* The tree topology supports:

  * tight clustering of Canada + Siberia
  * **Greenlandic Qanisartuut** forming a long branch (outgroup-like position)
  * internal structure among Greenland sites.

---

## 🧾 Key results

* The analysis reveals **clear population structure** in Arctic foxes despite their long-distance dispersal ability.
* **Canadian and Siberian populations are highly similar genetically**, forming a homogeneous cluster in PCA, ADMIXTURE and the NJ tree.
* **Greenlandic populations are more differentiated**:

  * **Qanisartuut** repeatedly appears as a **distinct, isolated subpopulation** with:

    * higher FST to other groups
    * unique ancestry components in ADMIXTURE
  * Other Greenlandic regions (e.g. Zackenberg, Scoresbysund, Kangerlussuaq) show varying degrees of similarity to each other and to the Canada–Siberia cluster.
* Evidence of **recent migration** is visible in ADMIXTURE, particularly in some individuals from **Zackenberg**.
* Overall, **Greenlandic foxes have lower genetic diversity** than their Canadian and Siberian counterparts, suggesting **historical bottlenecks or reduced connectivity**.

---

## 🌍 Interpretation & climate context

* Sea ice acts as a **bridge** connecting Arctic fox populations across oceans.
* The **close relationship** between Canadian and Siberian foxes is consistent with **frequent migration across sea ice**.
* As **climate change reduces sea-ice cover**, migration routes may close, increasing isolation and reducing genetic diversity in coastal or island populations.

These findings highlight that, even for a highly mobile species, **population structure and local isolation can emerge**, with important implications for **conservation and management** under ongoing climate change.

---

## 🙏 Acknowledgements

This repository is based on the report:

> *“Isolated or connected? Genetic structure dynamics of the Arctic fox (Vulpes lagopus)”*
> NBIA09043U Population Genetics 2025, University of Copenhagen.

Analyses use:

* **PLINK** (QC and PCA)
* **ADMIXTURE** + **evalAdmix**
* **R** + tidyverse / adegenet / ape (for FST, heterozygosity, DAPC and NJ).

---

