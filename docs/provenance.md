# Provenance and limits

## Source of truth

The attached manuscript supplies the section order and scientific narrative. Research scripts and existing results come from the author's `LUAD_project` directory. The author explicitly confirmed **six compound-screening targets, excluding VPS72**, and requested that the manuscript remain unchanged. The eight-gene prognostic signature and six-target compound screen are separate analyses.

The pre-existing repository was retained and extended with diagnostic ROC, eight-cohort MP survival and adjusted Cox, exact SHAP, pinned FOXM1 DepMap, GraphBAN inference preparation/filtering and docking scripts. This update does not retrain the prognostic model or replace the manuscript's results with newly generated scientific claims.

`source_inventory.json` records 122 source files found before packaging, including the restored third-party GraphBAN source. It records paths, line counts, SHA-256 hashes, package references and selected input/output calls. Installed `.venv` packages, binary model weights, the large GraphBAN zip and generated outputs are excluded. Inventory inclusion is not a claim that all scripts are used by the manuscript.

## Preserved distinctions

- The README presents only the final nine meta-programs and their states. Internal saved gene-list keys remain compatible with the existing analysis inputs. The broken `MP_list[Cluster_*]` indexing in the old packaged consensus script was removed; gene-list identities are retained for downstream explicit mapping. The saved RData now explicitly contains `MP_list` and `Cluster_list`; named arguments to `save()` did not rename the original lowercase variables.
- The original assignment implementation treats a maximum score **equal to 0.8** as active (`max_score < 0.8` is inactive), while the manuscript uses “exceeded 0.8.” The original boundary behavior is retained. No existing assignments were recomputed.
- The original QC implementation retains genes expressed in **more than 10** cells, whereas the manuscript excludes genes detected in fewer than 10 cells. Historical filtering is retained rather than silently changing the dataset.
- The MP clustering code uses an initial seed-neighborhood stopping rule, `sorted_intersection[1] > min_group_size`. This is not identical to a simple final cluster-size filter of at least five. The historical algorithm is retained; reproducing the nine MPs depends on the original inputs and tie handling.
- The packaged Monocle3 workflow roots the recorded analysis at cluster 2. CNV and CytoTRACE are supporting analyses; the script does not automatically discover a biologically optimal root in new datasets.
- Existing spatial expression and region plots read `CARD_list.rds`, which stores SCT-normalized objects and the historical region labels. RCTD is the main spatial deconvolution analysis. The historical CARD preparation helper is retained under `scripts/06_deconvolution/legacy/` to make this dependency explicit.
- The original model-selection workflow applies ComBat jointly to the cohorts. This preprocessing is retained for reproduction; it is not a frozen prospective deployment recipe. Exact post-ComBat inputs are now saved before modeling so SHAP can verify the same predictions.
- The eight-cohort adjusted Cox script preserves its source clinical endpoint handling, including the recorded TCGA endpoint choice and GSE30219 month-to-day conversion. It models stage/T/N/M as numeric terms when available, performs complete-case selection and records non-estimable terms. These are not newly inferred clinical annotations.

## Compound-screening history

Historical source files and receptor manifests include VPS72. Current screening preparation, consensus output and the Vina target list use only FOXM1, TFDP1, FOSL2, KLF4, CEBPG and ATF1. The consensus script retains a separate all-target audit, allowing historical prediction tables to be checked without presenting VPS72 as a current screening target.

The historical receptor manifest records FOXM1 PDB 3G73 and AlphaFold structures for other targets. Preparation scripts and evidence retain those actual structural sources. No receptor source was relabeled as PDB, no new docking calculation was substituted and the manuscript was not edited. `legacy_prepare_alphafold_receptors.sh` is explicitly historical preparation, not evidence of experimental structures for every TF.

`evidence/Table_2_2_molecular_docking.csv` and `evidence/vina_batch_summary.csv` preserve source tables, including historical VPS72 rows. `evidence/manuscript_docking_pairs.csv` is the current six-target view: eight interactions involving six TFs and five compounds, all below −4 kcal/mol. The energies match the supplied manuscript. Filtering this evidence table does not re-execute Vina or validate binding experimentally.

## Files and licensing

Raw expression data, full fitted models, pretrained weights, motif databases, licensed DrugBank libraries, receptor coordinates and manuscript files are not uploaded by this update. Selected small result tables are included as provenance. Patient-level clinical input tables remain external. No credential is copied into a release archive. GraphBAN retains its included upstream MIT license; an old unsupported repository-wide MIT badge was removed rather than inventing a license decision for the author's research code.
