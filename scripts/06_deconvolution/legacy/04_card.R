source("config.R")

# dir.create("~/R/SeuratV4")
# install.packages('Seurat', repos = c('https://satijalab.r-universe.dev'), lib = "~/R/SeuratV4")
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggsci)
library(CARD)
library(MuSiC)
library(ggpubr)
library(lemon)
library(patchwork)


LUAD_Epi_assigned_MP  <- readRDS(file.path(RESULT_2_DIR,"results_data","04_LUAD_Epi_assigned_MP_final.rds"))

dir <- c("seurat_ST1" ,"seurat_ST2","seurat_ST3","seurat_ST5","seurat_ST6","seurat_ST8")
stage <- c('IAC','IAC','MIA','AIS','MIA','AIS')

CARD_list <- list()

for (i in 1:length(dir)) {
  load(file.path(ST_DATA_DIR,paste0(dir[i],".rda")))
  print(paste0(dir[i],"_",stage[i]))

  ###PCA
  seurat <- SCTransform(seurat, assay = "Spatial", verbose = FALSE)
  seurat <- RunPCA(seurat, assay = "SCT", verbose = FALSE)
  #ElbowPlot(seurat, ndims = 50)
  seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:15)
  seurat <- FindClusters(seurat, verbose = FALSE,resolution = 0.5)
  table(seurat$seurat_clusters)
  #SpatialDimPlot(seurat, cells.highlight = CellsByIdentities(object = seurat, ), facet.highlight = TRUE, ncol = 4)
  #SpatialFeaturePlot(seurat, features = c("SFTPC","HLA-DRA","COL14A1","MDK"))

  seurat$subregion <- "subregion"

  if(dir[i] == "seurat_ST1"){
    #TD1 IAC
    seurat$subregion[seurat$seurat_clusters %in% c("0",'5',"7")]= "Cancer"
    seurat$subregion[seurat$seurat_clusters %in% c("6")]= "Stromal"
    seurat$subregion[seurat$seurat_clusters %in% c("1","2","8")]= "Normal epithelium"
    seurat$subregion[seurat$seurat_clusters %in% c("3","4")]= "Lymph"
  }

  if(dir[i] == "seurat_ST2"){
    #TD2 IAC
    seurat$subregion[seurat$seurat_clusters %in% c("0","1","2","3","6" )]= "Cancer"
    seurat$subregion[seurat$seurat_clusters %in% c("7")]= "Stromal"
    seurat$subregion[seurat$seurat_clusters %in% c( "8")]= "Normal epithelium"
    seurat$subregion[seurat$seurat_clusters %in% c("4","5")]= "Lymph"
  }

  if(dir[i] == "seurat_ST3"){
    #TD3 MIA
    seurat$subregion[seurat$seurat_clusters %in% c("0","1","3","4"  )]= "Cancer"
    seurat$subregion[seurat$seurat_clusters %in% c("2")]= "Stromal"
    seurat$subregion[seurat$seurat_clusters %in% c("5")]= "Lymph"

  }

  if(dir[i] == "seurat_ST5"){
    #TD5 AIS
    seurat$subregion[seurat$seurat_clusters %in% c("0" )]= "Cancer"
    seurat$subregion[seurat$seurat_clusters %in% c("2","4")]= "Stromal"
    seurat$subregion[seurat$seurat_clusters %in% c( "3","5","6","7")]= "Lymph"
    seurat$subregion[seurat$seurat_clusters %in% c( "1","8")]= "Normal epithelium"
  }

  if(dir[i] == "seurat_ST6"){
    #TD6 MIA
    seurat$subregion[seurat$seurat_clusters %in% c("0","3","4","5","6" )]= "Cancer"
    seurat$subregion[seurat$seurat_clusters %in% c("2")]= "Stromal"
    seurat$subregion[seurat$seurat_clusters %in% c("1" )]= "Lymph"
  }

  if(dir[i] == "seurat_ST8"){
    #TD8 AIS
    seurat$subregion[seurat$seurat_clusters %in% c("1","4" )]= "Cancer"
    seurat$subregion[seurat$seurat_clusters %in% c("6","7")]= "Stromal"
    seurat$subregion[seurat$seurat_clusters %in% c("2","3","5")]= "Lymph"
    seurat$subregion[seurat$seurat_clusters %in% c( "0" )]= "Normal epithelium"
  }
  spatial_location = seurat@images[["slice1"]]@coordinates[,c("row","col")]
  colnames(spatial_location) = c("x","y")

  CARD_obj = createCARDObject(
    sc_count = GetAssayData(LUAD_Epi_assigned_MP,layer = "counts"),
    sc_meta = LUAD_Epi_assigned_MP@meta.data,
    spatial_count = seurat@assays[["Spatial"]]@counts,
    spatial_location =  spatial_location,
    ct.varname = "MP",
    ct.select = unique(LUAD_Epi_assigned_MP$MP),
    sample.varname = "Sample" )

  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

  seurat <- AddMetaData(object = seurat,
                            metadata = CARD_obj@Proportion_CARD )

  CARD_list[[paste0(dir[i],"_",stage[i])]] <-  seurat

}

saveRDS(CARD_list, file =  file.path(RESULT_6_DIR,"results_data","CARD_list.rds")  )

source(file = file.path(CODE_DIR,"scripts/06_deconvolution/legacy/card_visualize.R"))
colors =c( "#E64B35FF" ,"#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
           "#91D1C2FF" ,"#DC0000FF", "#7E6148FF")

for (i in names(CARD_list)) {
  print(i)
  p1 = SpatialPlot(CARD_list[[i]] , group.by = "subregion",label = TRUE, label.size = 5)+
    scale_color_npg()+scale_fill_npg()+
    ggplot2::theme(legend.position = "none")
  spatial_location = CARD_list[[i]]@images[["slice1"]]@coordinates[,c("row","col")]
  colnames(spatial_location) = c("x","y")
  p2 = CARD.visualize.pie2(proportion = CARD_list[[i]]@meta.data[,paste0("MP_",1:9)],
                           spatial_location = spatial_location, colors = colors, round_size=0.7)

  pdf(file.path(RESULT_6_DIR,"results_figure",paste0(i,"_ratio.pdf")),width = 5,height = 8)
  print(p1+p2+ plot_layout(widths = c(1, 2)) )
  dev.off()
}

SpatialFeaturePlot(CARD_list$seurat_ST1_IAC, features = c("VPS72"))


CARD_list <- readRDS(file =  file.path(RESULT_6_DIR,"results_data","CARD_list.rds") )

a <- list()
for (i in names(CARD_list)) {

  a[[i]] = CARD_list[[i]]@meta.data[which(CARD_list[[i]]@meta.data$subregion == "Cancer"),paste0("MP_",1:9)]
  a[[i]]$group = i
}

a <- a[c( "seurat_ST5_AIS","seurat_ST6_MIA","seurat_ST2_IAC")]
a <- a[c( "seurat_ST8_AIS","seurat_ST3_MIA","seurat_ST1_IAC")]

combined_df <- bind_rows(a, .id = NULL)


combined_df$group <- factor(combined_df$group,levels = c(  "seurat_ST8_AIS","seurat_ST3_MIA","seurat_ST1_IAC"))
group=levels(factor( combined_df$group ))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}

boxplots <- list()
for (i in colnames(combined_df)[-10]) {

  boxplots[[i]] <- ggboxplot(combined_df,"group", i,  fill = "group",  outlier.shape = NA,title = i) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'none',
          axis.title.x = element_blank() )+ stat_compare_means(method = "t.test", hide.ns = FALSE,
                                                               comparisons = my_comparisons,  label = "p.adj.format",
                                                               p.adjust.method = "BH",vjust=0.02,bracket.size=0.6)  +
    scale_fill_manual(values= c( "#79B99D", "#4DBBD5FF", "#E64B35FF"  ))

}

pdf(file.path(RESULT_6_DIR,"results_figure","ST_MP_cellratio_CARD2.pdf"),width = 12,height = 12)
wrap_plots(boxplots, ncol = 3, nrow = 3) +
  plot_layout(guides = "collect") & theme(legend.position = "top")
dev.off()
