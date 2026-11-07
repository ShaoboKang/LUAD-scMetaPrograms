# Reproduction guide

## Environments and paths

Use an R environment with Seurat v5, Harmony, decontX, scDblFinder, infercnv, NMF, clusterProfiler, GSVA, survival, glmnet and plotting packages. Individual steps additionally need CellChat, BayesPrism, IOBR, spacexr, Monocle3 (manuscript v1.3.7), Slingshot, CytoTRACE, GeneSwitches, Scissor, Mime1 and SCENIC loom support. These packages span CRAN, Bioconductor and upstream repositories; the old universal Bioconductor installation snippet was removed. An automatically inferred dependency list is in `docs/r_packages.tsv`; it is an inventory, not a tested environment lockfile.

Use separate Python environments for pySCENIC, GraphBAN and docking. GraphBAN needs its upstream PyTorch/DGL, ESM, transformers and feature dependencies. Docking needs RDKit (manuscript v2025.3.5), pandas, Meeko and AutoDock Vina. PyMOL visualization was external to the recorded batch scripts. There is no complete historical environment lockfile in the source project, so version-exact end-to-end reproduction is not asserted.

Set `LUAD_ROOT` to the data/results tree and `LUAD_CODE_ROOT` to this repository. Run R scripts from the repository root. English environment aliases `LUAD_DIAGNOSTIC_DIR` and `LUAD_CLINICAL_FILE` can override locally named diagnostic and clinical inputs without renaming them.

## Ordered analysis

Use [workflow/manuscript_order.tsv](../workflow/manuscript_order.tsv) as the executable-file index. It is deliberately not an automatic all-steps launcher: several stages require manually reviewed annotations and external tools.

1. Run preprocessing in order; supply the inferCNV gene-position file and inspect malignant/normal annotations. Original cluster label choices are preserved.
2. Prepare per-tumor NMF input, run NMF, identify consensus MPs, then run functional/state validation. Script `03_mp_signature_identification.R` both identifies MPs and assigns cells; `05_mp_state_assignment.R` is downstream visualization. Saved MP signatures retain legacy keys. Never apply the display remapping twice.
3. Run SCENIC input preparation, external pySCENIC commands and R RSS/network steps. The supplied motif databases, loom files and TF reference lists must match the study. TME subtyping can run in parallel with these analyses after atlas annotation.
4. Prepare annotated spatial objects before VPS72 and L–R spatial plots. The historical script `scripts/06_deconvolution/legacy/04_card.R` creates `CARD_list.rds`; it requires a compatible Seurat/CARD environment. Its manual cluster-to-region assignments are specific to the original six samples. Existing saved objects may be used directly. This preparation dependency occurs earlier than the manuscript's main RCTD presentation.
5. Run CellChat/L–R characterization, bulk scoring, CIBERSORT/BayesPrism, then RCTD and its visualization. Supply external reference matrices, genomic annotations and spatial image slots as expected by each method.
6. Run `07_trajectory/00_prepare_harmony.R` to produce the shared integrated MP object. Construct Monocle3 CDS before Slingshot, which reads the CDS. Root selection is fixed to the recorded cluster; use CNV/stemness analyses to review its biological justification. TF pseudotime and GeneSwitches need the SCENIC loom and recorded MP signatures.
7. Run MP survival, Scissor and diagnostic ROC. For the prognostic signature, run `08_clinical/08_ml_tf_selection.R`, then the eight-cohort adjusted Cox script and exact SHAP. The preliminary LASSO in `01_univariate_survival.R` is not the final eight-TF model.
8. Download the pinned DepMap release separately, then run dependency analysis. The download validates MD5 and file sizes from release metadata. The primary subset is non-engineered LUAD cell lines with observed FOXM1 (2305) scores.
9. Run six-target GraphBAN and docking as below.

Scripts `08_clinical/04`–`07` and `09`–`11`, tumor-ratio exploration and legacy CARD analyses retain exploratory/supplementary work. Their presence does not designate them as substitutes for the final model or main spatial results.

## Six-target compound workflow

The author-confirmed target panel is FOXM1, TFDP1, FOSL2, KLF4, CEBPG and ATF1. The screening input validator rejects VPS72 and any missing/extra target. The eight-TF prognostic model is unchanged.

```bash
python3 scripts/11_compounds/01_prepare_screen.py proteins_six_targets.fasta drugbank_id_smiles.tsv /path/to/df_test_in.csv
export GRAPHBAN_ROOT="$PWD/vendor/GraphBAN"
export GRAPHBAN_INPUT=/path/to/df_test_in.csv
export GRAPHBAN_CHEMBERTA=/path/to/ChemBERTa-77M-MTR
export GRAPHBAN_BIOSNAP_MODEL=/path/to/biosnap_checkpoint.pth
export GRAPHBAN_KIBA_MODEL=/path/to/kiba_checkpoint.pth
export GRAPHBAN_BINDINGDB_MODEL=/path/to/bindingdb_checkpoint.pth
bash scripts/11_compounds/02_predict_graphban.sh
python3 scripts/11_compounds/03_consensus_candidates.py biosnap.csv kiba.csv bindingdb.csv candidates.csv
```

The vendored predictor uses its YAML `RESULT.OUTPUT_DIR` and writes the output basename there. Inspect that directory for the three CSVs. Supply checkpoints and feature resources separately. Protein embeddings retain the original 1,022-residue truncation. `Y=1` is an input-format placeholder, not a measured binding label. The filter requires identical pair sets across models, unique pair identifiers, consistent SMILES and all three scores strictly greater than 0.9. Historical extra targets are retained only in a separate all-target audit table.

Create a separate docking work directory with `input/GraphBAN_candidates.csv`, `input/receptors.tsv`, `receptors_prepared/`, ligand structures and `results/`. Copy `scripts/11_compounds/docking/receptors_six_targets.tsv` as the receptor manifest **after reviewing the actual structural files and boxes**. Its entries preserve the historical structural sources. The optional `legacy_prepare_alphafold_receptors.sh` reproduces the historical AlphaFold preparation, with VPS72 removed from its target list.

```bash
export LUAD_DOCKING_DIR=/path/to/docking
python3 scripts/11_compounds/docking/01_smiles_to_3d.py "$LUAD_DOCKING_DIR/input/GraphBAN_candidates.csv"
python3 scripts/11_compounds/docking/02_retry_ligands.py
# Default is preflight; this does not perform docking:
python3 scripts/11_compounds/docking/03_run_vina.py
# Execute only once receptors, ligand PDBQTs and boxes pass preflight:
python3 scripts/11_compounds/docking/03_run_vina.py --run --workers 2 --cpu-per-job 2 --seed 2026
```

Ligand SDF-to-PDBQT conversion and receptor review use Meeko/external preparation. The batch script reports missing or invalid inputs rather than fabricating results. `04_summarize_docking.py` filters the original standardized Protein/Ligand/Component/Affinity table; it does not parse pose files into new experimental claims.

## Reusing the existing fitted model

The source project stores its historical Mime result and exact input list under the TF-selection work directory using a legacy filename. To run the new exact SHAP script without retraining, place or link these two files under `Result_8/Mime_TF_selection/` as:

- `Mime_res_batch_corrected.rds`
- `data_list_for_SHAP.rds`

They must be the paired originals. The SHAP function accepts a custom `model_dir` and `output_dir` when loaded interactively. It checks the selected eight genes, sample matching, saved-risk reconstruction and SHAP additivity. It neither recomputes ComBat nor refits the model.
