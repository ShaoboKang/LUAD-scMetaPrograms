# Coding Standards for LUAD-scMetaPrograms

## 1. Language Policy

- **ALL** comments, variable names, function names, documentation, and output messages must be in **English**.
- No Chinese characters (including in strings, file paths within code, comments, or variable names).
- Original Chinese path segments in `config.R` should use English aliases.

## 2. Function-Based Architecture

Every analysis script must follow this template:

```r
# ============================================================================
# Script: <script_name>.R
# Module: <module_number>_<module_name>
# Description: <one-line summary>
# ============================================================================

# ---- 1. Configuration & Dependencies ----
source("config.R")
library(Seurat)
# ... other packages

# ---- 2. Helper Functions ----

#' Brief description of what this function does.
#'
#' @param param1 Description of param1.
#' @param param2 Description of param2.
#' @return Description of return value.
#' @examples
#' result <- my_helper_function(x, y)
my_helper_function <- function(param1, param2) {
  # implementation
}

# ---- 3. Main Analysis Functions ----

#' Main analysis function for this script.
#'
#' @param input_path Path to input data file.
#' @param output_dir Directory for output files.
#' @return List of output file paths.
run_analysis <- function(input_path, output_dir) {
  # Step 1: Load data
  # Step 2: Process
  # Step 3: Save results
  # Step 4: Return metadata
}

# ---- 4. Entry Point ----

if (!interactive()) {
  message("[INFO] Running analysis pipeline...")
  run_analysis(
    input_path = file.path(SC_DATA_DIR, "input.rds"),
    output_dir = file.path(RESULT_1_DIR, "results_data")
  )
}
```

## 3. Naming Conventions

| Element | Convention | Example |
|---------|-----------|---------|
| Functions | `snake_case`, verb-first | `merge_and_qc`, `assign_mp_states` |
| Variables | `snake_case` | `cell_count`, `mp_scores` |
| Constants | `UPPER_SNAKE_CASE` | `MIN_SCORE`, `PROGRAM_SIZE` |
| File names | `snake_case.R` | `01_merge_and_qc.R` |
| Directories | `snake_case` | `01_preprocessing` |

## 4. Documentation Standards

- Every function must have a roxygen2-style header.
- Every script must have a module header with description.
- Complex logic blocks must have inline comments in English.
- Use `message()` / `warning()` for runtime feedback, never `cat()`.

## 5. Path Management

- All paths must use `config.R` variables: `SC_DATA_DIR`, `BULK_DIR`, `RESULT_X_DIR`, etc.
- Never hard-code absolute paths like `~/LUAD_project/...` or `F:/...`.
- Use `file.path()` for cross-platform compatibility.

## 6. Error Handling

- Use `tryCatch()` for operations that may fail (e.g., missing genes in Cox regression).
- Validate inputs with `stopifnot()` or explicit `if (...) stop("...")`.
- Error messages must be in English.

## 7. Output Organization

- Results saved to `RESULT_X_DIR/results_data/`
- Figures saved to `RESULT_X_DIR/results_figure/`
- Use descriptive filenames: `01_forestplot.pdf`, `mp_scores_heatmap.pdf`
