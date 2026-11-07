# LUAD-scMetaPrograms

Analysis code for **Identifying consensus tumor cell states and characterizing their roles in lung adenocarcinoma progression and clinical outcomes**.

Shaobo Kang, Yijing Zhang, Renjie Dou, Yansong Jin, and Yanyan Ping.

This repository follows the manuscript's Results 3.1–3.7: establish the single-cell atlas and nine consensus states; characterize transcriptional regulation and communication; examine progression and trajectories; evaluate clinical relevance; interpret the prognostic model and prioritize compounds. The manuscript itself has not been edited.

The final prognostic signature contains **eight TFs: FOXM1, TFDP1, FOSL2, KLF4, CEBPG, ATF1, POU6F1, XBP1**. The final model is **StepCox (forward) + elastic net (alpha = 0.7)**. Compound screening uses **six targets: FOXM1, TFDP1, FOSL2, KLF4, CEBPG, ATF1**, as confirmed by the author. VPS72 remains part of the state-regulation analyses in section 3.2; it is excluded from the compound-screening workflow.

## Manuscript order

| Results | Analysis | Code | Main figures |
|---|---|---|---|
| 3.1 | Atlas, epithelial malignancy, NMF, consensus MPs, state characterization | `scripts/01_preprocessing/`, `02_meta_programs/`, `03_state_characterization/` | Fig. 1; S1–S2 |
| 3.2 | SCENIC, state-associated TFs, TF–target networks, VPS72 characterization | `scripts/04_transcription_factors/` | Fig. 2; S3 |
| 3.3 | TME subtyping, CellChat, L–R expression and survival | `scripts/05_tme_communication/` | Fig. 3; S4–S5 |
| 3.4 | State dynamics across stages, bulk and spatial deconvolution | `scripts/02_meta_programs/06_mp_stage_association.R`, `scripts/06_deconvolution/` | Fig. 4; S6–S7 |
| 3.5 | Monocle3, Slingshot, CNV/stemness, TF dynamics, GeneSwitches | `scripts/07_trajectory/` | Fig. 5; S8 |
| 3.6 | MP survival, Scissor, diagnostic ROC, eight-TF model, adjusted survival | `scripts/08_clinical/`, `scripts/09_diagnostic/` | Fig. 6; S9 |
| 3.7 | Exact SHAP, FOXM1 DepMap, six-target GraphBAN, docking | `scripts/08_clinical/14_exact_shap.R`, `scripts/10_dependency/`, `scripts/11_compounds/` | Fig. 7; S10 |

See [the complete Methods crosswalk](docs/manuscript_crosswalk.md), [ordered workflow](workflow/manuscript_order.tsv), [reproduction guide](docs/reproduction.md), and [provenance and limitations](docs/provenance.md). Module numbers preserve the established repository paths; the workflow table gives the manuscript reading order and cross-module prerequisites.

## Data and configuration

```bash
git clone https://github.com/ShaoboKang/LUAD-scMetaPrograms.git
cd LUAD-scMetaPrograms
export LUAD_ROOT=/path/to/LUAD_project
export LUAD_CODE_ROOT="$PWD"
```

Run R scripts **from the repository root**. [config.R](config.R) separates the external data/results directory from the code repository. Inputs, large Seurat objects, pretrained weights, DrugBank structures, reference databases and manuscript files are not bundled. See [data/README.md](data/README.md) for input contracts and accession numbers.

The published study uses six single-cell cohorts (432,276 input cells from 88 samples), malignant-cell programs from 45 tumors, TCGA-LUAD plus seven GEO survival cohorts, seven diagnostic cohorts, six spatial samples, and DepMap 24Q4 Public v1. These study counts are reference values, not assertions that a fresh run has regenerated every result.

## Run selected analyses

```bash
# Atlas and consensus MPs: run preprocessing scripts first, then:
Rscript scripts/02_meta_programs/01_nmf_input_preparation.R
Rscript scripts/02_meta_programs/02_nmf_per_sample.R
Rscript scripts/02_meta_programs/03_mp_signature_identification.R

# Exact model inputs are saved alongside the Mime result for later SHAP.
Rscript scripts/08_clinical/08_ml_tf_selection.R
Rscript scripts/08_clinical/13_risk_multivariable_eight_cohorts.R
Rscript scripts/08_clinical/14_exact_shap.R
Rscript scripts/09_diagnostic/01_mp_diagnostic_auc.R

# FOXM1 dependency (pinned release; downloading the matrix is a separate step):
python3 scripts/10_dependency/01_download_depmap.py
Rscript scripts/10_dependency/02_foxm1_dependency.R
```

These examples are individual stages, not a one-command end-to-end pipeline. The [reproduction guide](docs/reproduction.md) identifies manually curated annotations, intermediate objects, spatial prerequisites and optional exploratory scripts. NMF, SCENIC, deconvolution and model training require substantial computation and method-specific environments.

## State identifiers

The original saved gene-list identifiers are preserved. Apply this mapping once when displaying manuscript states:

| Legacy MP | Manuscript MP | State |
|---|---|---|
| MP_1 | MP1 | G2M |
| MP_8 | MP2 | G1S |
| MP_2 | MP3 | NFkB (NFκB) |
| MP_7 | MP4 | MHC-II |
| MP_4 | MP5 | Translation |
| MP_6 | MP6 | OXPHOS |
| MP_5 | MP7 | Adhesion |
| MP_3 | MP8 | EMT |
| MP_9 | MP9 | Glycolysis |

The MP gene list must not be indexed using `Cluster_*` keys. The consensus identification script now checks the expected nine `MP_*` keys before saving. Existing serialized objects should retain their original labels; avoid remapping an already remapped object.

## Validation and evidence

[validation/](validation/) records checks performed for this release. [evidence/](evidence/) contains selected small numerical tables from the research project, including the FOXM1 dependency summary and eight retained docking interactions. Full datasets and a complete end-to-end rerun are outside the packaging checks.

Source files are catalogued with SHA-256 hashes, package references and input/output calls in [docs/source_inventory.json](docs/source_inventory.json). Installed virtual environments and generated outputs are excluded from that source inventory. The vendored GraphBAN inference source retains its upstream MIT license; see [third-party notes](docs/third_party.md). No new license is imposed on the remaining research code in this update.

Please cite the manuscript using its final publication details when available. For questions, open an issue in this repository.
