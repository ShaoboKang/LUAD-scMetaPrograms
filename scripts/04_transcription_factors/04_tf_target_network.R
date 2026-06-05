#' Build TF-Target Network from pySCENIC Output
# Script: 04_tf_target_network
#'
#' @description
#' Integrates pySCENIC regulons, MP gene lists, and cluster DEGs to build
#' a filtered TF-target network (edges and nodes) for visualization.
#'
#' @return Invisibly returns NULL. Writes CSV files to RESULT_4_DIR.
#'
#' @examples
#' \dontrun{
#' run_tf_target_network()
#' }
# Auto-locate config.R
config_path <- "config.R"
if (!file.exists(config_path)) {
  config_path <- file.path(dirname(dirname(dirname(getwd()))), "config.R")
}
if (!file.exists(config_path)) {
  config_path <- "~/LUAD-scMetaPrograms/config.R"
}
source(config_path)
suppressMessages({
  library(SCENIC)
  library(SCopeLoomR)
  library(dplyr)
  library(Seurat)
})

#' Build TF-Target Network
#'
#' @return Invisibly returns NULL.
run_tf_target_network <- function() {
  message("[04_tf_target_network] Starting network construction...")

  mp_rdata <- file.path(RESULT_2_DIR, "results_data", "02_Malig_final_MP_top50.RData")
  load(mp_rdata)

  tf_rds <- file.path(RESULT_4_DIR, "results_data", "01_top5_TF.rds")
  TF_names <- readRDS(tf_rds)
  TF_names <- rev(levels(TF_names$Topic))

  TFall <- data.frame(
    MP = c(
      rep("MP_1", 5), rep("MP_2", 6), rep("MP_3", 9), rep("MP_4", 7),
      rep("MP_5", 9), rep("MP_7", 6), rep("MP_8", 9), rep("MP_9", 5)
    ),
    TF = gsub("\\(.*\\)", "", TF_names)
  )

  loom_file <- file.path(RESULT_4_DIR, "scenic", "output", "out_SCENIC.loom")
  loom <- open_loom(loom_file)
  regulons <- regulonsToGeneLists(get_regulons(loom, column.attr.name = "Regulons"))
  regulons <- regulons[TF_names]
  regulons <- stack(regulons)
  colnames(regulons) <- c("Gene", "TF")
  close_loom(loom)

  adj_file <- file.path(RESULT_4_DIR, "scenic", "output", "adj.sample.tsv")
  TF_target <- read.delim(file = adj_file, sep = "\t")
  TF_target <- TF_target[TF_target$TF %in% TFall$TF, ]
  TF_target <- TF_target[TF_target$target %in% regulons$Gene, ]

  TF_target_with_MP <- merge(TF_target, TFall, by = "TF", all.x = TRUE)

  TF_target_with_MP <- TF_target_with_MP[TF_target_with_MP$target %in% unlist(MP_list), ]
  TF_target_with_MP <- TF_target_with_MP[TF_target_with_MP$importance > 0.3, ]

  mp_rds <- file.path(RESULT_2_DIR, "results_data", "04_LUAD_Epi_assigned_MP_final.rds")
  LUAD_Epi_assigned_MP <- readRDS(mp_rds)
  Idents(LUAD_Epi_assigned_MP) <- "MP"
  LUAD_Epi_assigned_MP$MP <- factor(
    LUAD_Epi_assigned_MP$MP,
    levels = paste0("MP_", seq_len(length(MP_list)))
  )

  DEGs <- FindAllMarkers(LUAD_Epi_assigned_MP, min.pct = 0.25, logfc.threshold = 1, only.pos = TRUE)
  DEGs <- DEGs[which(DEGs$p_val_adj < 0.05), ]

  TF_target_with_MP <- TF_target_with_MP[TF_target_with_MP$target %in% unique(DEGs$gene), ]

  edges <- TF_target_with_MP %>%
    dplyr::select(source = TF, target = target, weight = importance, MP) %>%
    mutate(interaction = "regulation")

  tf_nodes <- data.frame(
    name = unique(TF_target_with_MP$TF),
    type = "TF",
    MP = TF_target_with_MP$MP[match(unique(TF_target_with_MP$TF), TF_target_with_MP$TF)]
  )

  target_nodes <- data.frame(
    name = setdiff(unique(TF_target_with_MP$target), unique(TF_target_with_MP$TF)),
    type = "target",
    MP = TF_target_with_MP$MP[
      match(
        setdiff(unique(TF_target_with_MP$target), unique(TF_target_with_MP$TF)),
        TF_target_with_MP$target
      )
    ]
  )

  nodes <- bind_rows(tf_nodes, target_nodes)

  write.csv(edges, file = file.path(RESULT_4_DIR, "results_data", "edges.csv"), row.names = FALSE, quote = FALSE)
  write.csv(nodes, file = file.path(RESULT_4_DIR, "results_data", "nodes.csv"), row.names = FALSE, quote = FALSE)

  message("[04_tf_target_network] Completed successfully.")
  invisible(NULL)
}

if (!interactive()) {
  run_tf_target_network()
}
