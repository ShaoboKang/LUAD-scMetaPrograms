# Validate the persisted MP interface without running NMF or loading dependencies.
code <- parse("scripts/02_meta_programs/03_mp_signature_identification.R")
functions <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
  length(x) == 3L && is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function")), code)
e <- new.env()
for (x in functions) eval(x, e)
statements <- as.list(body(e$run_analysis))[-1]
e$cluster_list <- setNames(lapply(1:9, function(i) paste0("program", i)), paste0("Cluster_", 1:9))
e$mp_list <- setNames(lapply(1:9, function(i) paste0("GENE", i)), paste0("MP_", 1:9))
e$nmf_programs_original <- matrix(0, 1, 1)
e$nmf_intersect_original <- matrix(0, 1, 1)
e$output_data_dir <- tempfile("luad_mp_interface_")
dir.create(e$output_data_dir)
selected <- Filter(function(x) {
  t <- paste(deparse(x), collapse = " ")
  grepl('^new_order <-|^cluster_list <- cluster_list\\[new_order\\]|^stopifnot\\(setequal\\(names\\(mp_list\\)|^Cluster_list <-|^MP_list <-|^save\\(Cluster_list', t)
}, statements)
stopifnot(length(selected) == 6L)
for (x in selected) eval(x, e)
saved <- new.env()
load(file.path(e$output_data_dir, "02_Malig_final_MP_top50.RData"), saved)
stopifnot(all(c("MP_list", "Cluster_list") %in% ls(saved)),
  identical(names(saved$MP_list), paste0("MP_", 1:9)),
  identical(saved$MP_list$MP_8, "GENE8"),
  identical(names(saved$Cluster_list), paste0("Cluster_", c(1,8,2,7,4,6,5,3,9))))
unlink(e$output_data_dir, recursive = TRUE)
message("PASS: MP gene identities, ordered clusters and RData variable names preserved")
