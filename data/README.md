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

**Note**: This code package does not redistribute expression matrices. Obtain source datasets and observe their original terms of use.

## Additional manuscript inputs

The seven **diagnostic** cohorts differ from the survival cohorts: TCGA-LUAD, GSE10072, GSE32863, GSE63459, GSE19188, GSE68571 and GSE31210. `DIAGNOSTIC_DIR/all_datasets_merged_list.rds` is a named list of seven sample-by-gene tables. `TYPE` is `Y` for tumor and `N` for normal; `SAMPLE` and `STAGE` are optional metadata excluded from expression. Only GSE68571 receives the source workflow's `log2(pmax(x,0)+1)` transformation. Other inputs must already be on their study preprocessing scales.

Survival/modeling uses `bulk/TCGA_GEO7.rdata` containing `expr_list` (genes x samples) and `clin_list` (sample metadata with survival_time and vital_status). The eight-cohort adjusted Cox script additionally uses `CLINICAL_FILE` plus `bulk/TCGA/TCGA_sur_data.Rdata`. Inspect source clinical encodings and sample IDs before substituting data; the workflow preserves their historical handling.

FOXM1 dependency uses **DepMap 24Q4 Public, version 1**, article 27993248. `scripts/10_dependency/01_download_depmap.py` retrieves `Model.csv`, `CRISPRGeneEffect.csv` and the release README, validates checksums and sizes and records a download manifest. Outputs are placed under `Result_8/results_data/FOXM1_DepMap/`.

Compound screening requires a six-target FASTA (FOXM1, TFDP1, FOSL2, KLF4, CEBPG, ATF1), an ID/SMILES TSV, the three pretrained GraphBAN checkpoints, ESM and ChemBERTa resources. The manuscript library has 11,167 compounds; the preparation script reports if the supplied library count differs. Full DrugBank content and model weights are not included. Receptor structures, prepared PDBQTs, search boxes and ligand conversion tools are separate docking inputs. Historical sources are documented in `docs/provenance.md`.

Original processed input objects and all upstream download/preprocessing steps are not fully contained in this repository. Downloading public accessions alone does not recreate the exact curated RDS/RData objects without the recorded annotations and preprocessing.
