#!/usr/bin/env bash
# Run in the existing GraphBAN environment; pretrained weights remain external.
set -euo pipefail
: "${GRAPHBAN_ROOT:?Set GRAPHBAN_ROOT to the restored GraphBAN source directory}"
: "${GRAPHBAN_INPUT:?Set GRAPHBAN_INPUT to the absolute prepared CSV path}"
cd "$GRAPHBAN_ROOT"
python inductive_mode/predict.py --test_path "$GRAPHBAN_INPUT" --save_dir df_test_out_biosnap.csv --trained_model "${GRAPHBAN_BIOSNAP_MODEL:-$GRAPHBAN_ROOT/inductive_mode/trained_models/biosnap/inductive/seed12/result/best_model_epoch_45.pth}"
python inductive_mode/predict.py --test_path "$GRAPHBAN_INPUT" --save_dir df_test_out_kiba.csv --trained_model "${GRAPHBAN_KIBA_MODEL:-$GRAPHBAN_ROOT/inductive_mode/trained_models/kiba/inductive/seed12/result/best_model_epoch_7.pth}"
python inductive_mode/predict.py --test_path "$GRAPHBAN_INPUT" --save_dir df_test_out_bindingdb.csv --trained_model "${GRAPHBAN_BINDINGDB_MODEL:-$GRAPHBAN_ROOT/inductive_mode/trained_models/bindingdb/inductive/seed12/result/best_model_epoch_13.pth}"
