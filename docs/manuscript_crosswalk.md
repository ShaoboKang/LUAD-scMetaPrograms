# Manuscript-to-code crosswalk

Organized against the supplied manuscript's Methods 2.2.1–2.2.18 and Results 3.1–3.7. The manuscript file is unchanged. Directory numbers are stable implementation paths, not a replacement for manuscript section order.

| Methods | Results | Scripts relative to `scripts/` | Inputs / outputs and notes |
|---|---|---|---|
| 2.2.1 Preprocessing | 3.1 | `01_preprocessing/01_merge_and_qc.R` | `sc_data/*.rds` → merged QC object; decontX, doublets, gene filters |
| 2.2.2 Clustering and annotation | 3.1 | `01_preprocessing/02_cell_annotation.R`, `03_cell_proportions.R`, `04_epithelial_subtyping.R`; `05_tme_communication/01`–`05` | Sample metadata and manually reviewed cluster labels; eight major lineages and TME subtypes |
| 2.2.3 Malignant epithelial cells | 3.1 | `01_preprocessing/05_infercnv.R` | Epithelial object plus genomic gene positions → inferCNV and CNV scores |
| 2.2.4 Consensus MPs | 3.1 | `02_meta_programs/01_nmf_input_preparation.R`, `02_nmf_per_sample.R`, `03_mp_signature_identification.R`, `utils/` | >100 cells/sample; ranks 4–9, top 50 genes; robust overlaps 35/50 and 10/50; MP signatures and assignments |
| 2.2.5 Dominant MP states | 3.1 | `02_meta_programs/03_mp_signature_identification.R`, `05_mp_state_assignment.R` | Assignment is performed in script 03; script 05 visualizes assignments. Legacy identifiers are mapped for display |
| 2.2.6 Functional enrichment | 3.1 | `02_meta_programs/04_mp_functional_annotation.R`; `03_state_characterization/01`–`03` | MP gene lists / state DEGs → GO, Hallmark and pathway plots |
| 2.2.7 State-associated TFs | 3.2 | `04_transcription_factors/01_pyscenic_pipeline.R`, `02_vps72_characterization.R`, `03_vps72_spatial.R`, `04_tf_target_network.R` | Counts, motif databases, loom → RSS TF table, networks and VPS72 expression; spatial object prerequisite below |
| 2.2.8 Communication | 3.3 | `05_tme_communication/01`–`10` | Annotated TME populations + MP cells → CellChat object, L–R plots, spatial and survival support |
| 2.2.9 Bulk ssGSEA and deconvolution | 3.4 | `02_meta_programs/06_mp_stage_association.R`, `07_mp_prognosis_bulk.R`; `06_deconvolution/01_iobr_cibersort.R`, `02_bayesprism.R`, `03_tumor_ratio_clustering.R` | Bulk expression / single-cell reference → scores, fractions and stage comparisons |
| 2.2.10 Spatial RCTD | 3.4 | `06_deconvolution/04_rctd.R`, `05_rctd_visualization.R` | Six spatial samples + MP reference → RCTD fractions and cancer-region plots. Historical `CARD_list.rds` carries spatial annotations used by existing figures |
| 2.2.11 Trajectory and GeneSwitches | 3.5 | `07_trajectory/00_prepare_harmony.R`, `02_monocle3.R`, `03_cnv_stemness_scores.R`, `01_slingshot.R`, `04_tf_pseudotime_heatmap.R`, `05_clustergvis_all_mp_cells.R`, `06_gene_switch_analysis.R` | Integrated MP cells, CNV, stemness and loom → rooted CDS, lineages, TF pseudotime, switches (cutoff 0.2) |
| 2.2.12 Scissor | 3.6 | `08_clinical/02_scissor_analysis.R`, `utils/scissor_functions.R` | Bulk survival + integrated MP cells → Scissor labels; alpha 0.05 |
| 2.2.13 Diagnostic ROC | 3.6 | `09_diagnostic/01_mp_diagnostic_auc.R` | Seven diagnostic sample tables → ssGSEA, fixed-direction ROC/AUC with CIs, plots |
| 2.2.14 Prognostic model | 3.6 | `08_clinical/01_univariate_survival.R`, `08_ml_tf_selection.R`, `12_mp_survival_eight_cohorts.R`, `13_risk_multivariable_eight_cohorts.R` | State-associated TFs, TCGA+7 GEO → model comparison, final eight-TF model, KM and adjusted Cox. Script 01 includes preliminary LASSO exploration; its LASSO object is not the final model |
| 2.2.15 SHAP | 3.7 | `08_clinical/14_exact_shap.R` | Saved Mime fit + exact post-ComBat inputs → coefficient table, SHAP, prediction/additivity audits, summary figure |
| 2.2.16 CRISPR dependency | 3.7 | `10_dependency/01_download_depmap.py`, `02_foxm1_dependency.R`, `utils/ranked_landscape.R` | Pinned DepMap v1 with checksums → strict non-engineered LUAD subset, coverage and FOXM1 effects |
| 2.2.17 Compound screening | 3.7 | `11_compounds/01_prepare_screen.py`, `02_predict_graphban.sh`, `03_consensus_candidates.py`; `vendor/GraphBAN/inductive_mode/` | Six author-confirmed targets x compound library → three prediction tables → intersection requiring every score >0.9; VPS72 excluded |
| 2.2.18 Docking | 3.7 | `11_compounds/docking/01_smiles_to_3d.py`, `02_retry_ligands.py`, `03_run_vina.py`; `11_compounds/04_summarize_docking.py` | Ligand conformers, receptor/box manifest → preflight / Vina results → six-target pairs below −4 kcal/mol |

Exact panel reproduction also depends on the original plot settings, input versions and curated annotations. Section-level coverage does not claim a verified match for every panel. Historical differences that must not be silently rewritten are described in [provenance.md](provenance.md).
