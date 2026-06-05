# Data Directory

This directory is intentionally kept empty in the GitHub repository due to file size limitations and data sharing policies.

## Required Datasets

Please download and organize the following datasets according to the structure described in [`config.R`](../config.R) before running the analysis pipeline.

### 1. Single-cell RNA-seq (scRNA-seq)

Six published LUAD scRNA-seq datasets were used in this study:

| Dataset | Source | Access |
|---------|--------|--------|
| Bischoff et al. | 3CA | [Curated Cancer Cell Atlas](https://www.weizmann.ac.il/sites/3CA) |
| Kim et al. | 3CA | 3CA |
| Laughney et al. | 3CA | 3CA |
| Qian et al. | 3CA | 3CA |
| Xing et al. | 3CA | 3CA |
| Zhu et al. | GEO | [GSE189357](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189357) |

**Expected format**: Seurat objects (`.rds`) saved after initial filtering.

### 2. Spatial Transcriptomics

| Dataset | Source | Access |
|---------|--------|--------|
| Zhu et al. | GEO | [GSE189487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189487) |

Six samples: TD5 (AIS), TD8 (AIS), TD3 (MIA), TD6 (MIA), TD1 (IAC), TD2 (IAC).

### 3. Bulk Transcriptomics

| Cohort | Platform | Access |
|--------|----------|--------|
| TCGA-LUAD | RNA-seq | [TCGA Data Portal](https://portal.gdc.cancer.gov/) |
| GSE13213 | Microarray | GEO |
| GSE31210 | Microarray | GEO |
| GSE41271 | Microarray | GEO |
| GSE26939 | Microarray | GEO |
| GSE30219 | Microarray | GEO |
| GSE72094 | Microarray | GEO |
| GSE11969 | Microarray | GEO |

### 4. Reference Gene Sets

- **MSigDB v2024.1**: Download from [GSEA MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/)
- Place under `msigdb_v2024.1.Hs_GMTs/`

### 5. SCENIC Databases

- **cisTarget databases**: Download from [SCENIC resources](https://resources.aertslab.org/cistarget/)
- Required files: `hs_hgnc_curated_tfs.txt`, motif rankings, motif annotations

---

**Note**: Due to patient privacy and data use agreements, we cannot redistribute raw or processed expression data. All datasets are publicly available from the sources listed above.
