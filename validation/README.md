# Release validation

Checks performed on 2026-11-07:

- All **72 R files** parsed successfully, including configuration and the serialization check.
- All **21 Python source files** parsed successfully; shell scripts passed `bash -n`.
- All **56 indexed workflow stages** resolve to existing scripts and cover Results 3.1–3.7.
- `check_release.py` verified the six-target FASTA requirement, rejection of VPS72, strict three-model `>0.9` filtering (a score of exactly 0.9 is excluded), rejection of duplicate pair identifiers and the eight docking-evidence pairs across six TFs/five compounds.
- `check_mp_serialization.R` exercised the actual MP reorder/serialization statements on sentinel data: no gene lists become NULL, legacy `MP_*` identities are retained, cluster display order is preserved and the saved RData exports `MP_list` / `Cluster_list` as required downstream.
- **Actual saved-model SHAP validation** used the original paired Mime model and post-ComBat input objects, without retraining. All eight cohorts passed, totaling **1,701 samples**. Maximum saved-risk reconstruction error was approximately **9.94e-15**; maximum SHAP additivity error was approximately **8.88e-16**. See `shap_prediction_audit.csv`.
- Existing FOXM1 evidence was checked: **48/51** scored LUAD cell lines had negative gene effects. Existing docking evidence reproduced the eight six-target pairs and their recorded energies; Vina was not rerun.
- Changed files passed `git diff --check`. Release file contents were scanned for common credential/private-key patterns; `.git` is excluded from the archive.

Re-run lightweight checks from the repository root:

```bash
python3 validation/check_release.py
Rscript validation/check_mp_serialization.R
Rscript -e 'fs <- list.files(".", pattern="[.][Rr]$", recursive=TRUE, full.names=TRUE); for(f in fs) parse(f)'
```

These checks do **not** establish a fresh end-to-end reproduction of QC, NMF, SCENIC, communication, deconvolution, trajectories, model training, GraphBAN inference or molecular docking. Those stages require their original large datasets, curated annotations, external references and compatible environments. The manuscript was not edited. Historical implementation differences and external prerequisites remain documented in `docs/provenance.md` and `docs/reproduction.md`.
