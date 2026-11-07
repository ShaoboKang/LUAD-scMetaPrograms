# Third-party code

The NMF/consensus approach was adapted in the research project from PanCan_scMetab and the Gavish et al. pan-cancer framework. The scoring helpers under `scripts/02_meta_programs/utils/` derive from scalop. SCENIC, Scissor and Mime helper files retain their existing source comments and scientific dependencies; packaging does not establish new ownership or relicense these components.

The restored GraphBAN source identifies the upstream checkout as `HamidHadipour-GraphBAN-d6da0a7`. Its `inductive_mode` Python/YAML files are included under `vendor/GraphBAN/` with the original README and LICENSE. In `predict.py`, the machine-specific ChemBERTa path is replaced by the required `GRAPHBAN_CHEMBERTA` environment variable. Model architecture and inference calculations are otherwise retained. Copied text files have normalized line endings and trailing whitespace. Pretrained weights, ESM/ChemBERTa downloads and training datasets remain external.

The author-specific GraphBAN input preparation, score intersection and docking helpers live in `scripts/11_compounds/`. Current target selection excludes VPS72 as requested by the author. Original project source hashes are recorded in `docs/source_inventory.json`.
