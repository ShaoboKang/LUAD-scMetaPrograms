# FOXM1 loss-of-function fitness in lung adenocarcinoma cell lines.
# First run 08_download_FOXM1_DepMap.py to download pinned source files + MD5 checks.
# This is a descriptive analysis of DepMap Chronos estimates, not percent viability.
# No per-cell-line P values / confidence intervals are invented from single estimates.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
source("config.R")
root <- BASE_DIR
out<-file.path(root,"Result_8/results_data/FOXM1_DepMap")
rawdir<-file.path(out,"raw");fig<-file.path(root,"Result_8/results_figure/FOXM1_DepMap")
for(d in c(out,fig))dir.create(d,recursive=TRUE,showWarnings=FALSE)
files<-file.path(rawdir,c("Model.csv","CRISPRGeneEffect.csv","README.txt"))
stopifnot(all(file.exists(files)))
meta<-as.data.frame(fread(files[1],check.names=FALSE))
stopifnot(!anyDuplicated(meta$ModelID),all(c("ModelID","CellLineName","StrippedCellLineName","OncotreeCode","ModelType")%in%names(meta)))
header<-names(fread(files[2],nrows=0,check.names=FALSE))
# Resolve the FOXM1 symbol AND Entrez identifier; avoid confusing FOXM1 paralogs.
fox_col<-which(header=="FOXM1 (2305)")
stopifnot(length(fox_col)==1L)
effect<-as.data.frame(fread(files[2],select=c(1L,fox_col),check.names=FALSE))
names(effect)<-c("ModelID","FOXM1_GeneEffect")
stopifnot(!anyDuplicated(effect$ModelID),is.numeric(effect$FOXM1_GeneEffect))
stopifnot(all(effect$ModelID%in%meta$ModelID))
# Exact histology and model-type filter; NSCLC-unspecified, LUSC, SCLC and LCLC excluded.
luad<-meta[!is.na(meta$OncotreeCode)&meta$OncotreeCode=="LUAD"&
             !is.na(meta$ModelType)&meta$ModelType=="Cell Line",,drop=FALSE]
stopifnot(nrow(luad)>0L)
n_annotated_luad<-nrow(luad)
engineering_tag<-trimws(ifelse(is.na(luad$EngineeredModel),"",as.character(luad$EngineeredModel)))
engineered<-luad[nzchar(engineering_tag),,drop=FALSE]
engineered$FOXM1_GeneEffect<-effect$FOXM1_GeneEffect[match(engineered$ModelID,effect$ModelID)]
write.csv(engineered,file.path(out,"FOXM1_LUAD_engineered_models_excluded.csv"),row.names=FALSE)
luad<-luad[!nzchar(engineering_tag),,drop=FALSE]
luad$FOXM1_GeneEffect<-effect$FOXM1_GeneEffect[match(luad$ModelID,effect$ModelID)]
luad$In_CRISPR_matrix<-luad$ModelID%in%effect$ModelID
luad$Has_FOXM1_score<-is.finite(luad$FOXM1_GeneEffect)
luad$Score_below_zero<-ifelse(luad$Has_FOXM1_score,luad$FOXM1_GeneEffect<0,NA)
# Thresholds below are descriptive cuts only, not validated dependency calls.
luad$Score_le_minus_0_5<-ifelse(luad$Has_FOXM1_score,luad$FOXM1_GeneEffect<=-.5,NA)
luad$Score_le_minus_1<-ifelse(luad$Has_FOXM1_score,luad$FOXM1_GeneEffect<=-1,NA)
luad<-luad[order(!luad$Has_FOXM1_score,luad$FOXM1_GeneEffect,luad$ModelID),,drop=FALSE]
luad$Effect_rank<-NA_integer_;luad$Effect_rank[luad$Has_FOXM1_score]<-seq_len(sum(luad$Has_FOXM1_score))
rownames(luad)<-NULL
write.csv(luad,file.path(out,"FOXM1_LUAD_all_models_with_coverage.csv"),row.names=FALSE)
z<-luad[luad$Has_FOXM1_score,,drop=FALSE]
stopifnot(nrow(z)>0L)
write.csv(z,file.path(out,"FOXM1_LUAD_GeneEffect_ranked.csv"),row.names=FALSE)
all_effect<-cbind(effect,meta[match(effect$ModelID,meta$ModelID),setdiff(names(meta),"ModelID"),drop=FALSE])
write.csv(all_effect,file.path(out,"FOXM1_all_screened_models.csv"),row.names=FALSE)
values<-z$FOXM1_GeneEffect
summary_table<-data.frame(Release="DepMap 24Q4 Public, version 1",Gene="FOXM1",Entrez=2305,
  Annotated_LUAD_models=n_annotated_luad,Engineered_models_excluded=nrow(engineered),
  LUAD_cell_lines_in_primary_filter=nrow(luad),LUAD_with_FOXM1_score=nrow(z),
  LUAD_without_FOXM1_score=nrow(luad)-nrow(z),Median=median(values),
  Q1=as.numeric(quantile(values,.25)),Q3=as.numeric(quantile(values,.75)),
  Minimum=min(values),Maximum=max(values),Below_zero=sum(values<0),
  At_or_below_minus_0_5=sum(values<=-.5),At_or_below_minus_1=sum(values<=-1))
write.csv(summary_table,file.path(out,"FOXM1_LUAD_summary.csv"),row.names=FALSE)
# Explicitly report availability for commonly used candidate LUAD models.
requested<-c("A549","HCC827","NCIH1975","PC9","NCIH1650","NCIH358")
featured<-do.call(rbind,lapply(requested,function(name){
  m<-meta[!is.na(meta$StrippedCellLineName)&meta$StrippedCellLineName==name,,drop=FALSE]
  if(!nrow(m))return(data.frame(Requested_name=name,ModelID=NA,CellLineName=NA,
     OncotreeCode=NA,FOXM1_GeneEffect=NA_real_,Status="Not in this release metadata"))
  v<-effect$FOXM1_GeneEffect[match(m$ModelID,effect$ModelID)]
  data.frame(Requested_name=name,ModelID=m$ModelID,CellLineName=m$CellLineName,
    OncotreeCode=m$OncotreeCode,FOXM1_GeneEffect=v,
    Status=ifelse(m$OncotreeCode!="LUAD","Not annotated LUAD: excluded from strict LUAD filter",
       ifelse(is.finite(v),"Available","No FOXM1 score in this release")))
}))
write.csv(featured,file.path(out,"FOXM1_selected_cell_lines.csv"),row.names=FALSE)

# Ranked plot: one model-level effect estimate per cell line, never technical replicates.
z$Cell_label<-paste0(z$CellLineName,"  [",z$ModelID,"]")
z$Cell_label<-factor(z$Cell_label,levels=rev(z$Cell_label))
z$Value_label<-sprintf("%.3f",z$FOXM1_GeneEffect)
span<-max(diff(range(c(values,-1,0))),.5)
x_low<-min(c(values,-1))-.12*span
x_high<-max(c(values,0))+.2*span
caption<-paste("DepMap 24Q4 Public v1 | CRISPRGeneEffect.csv | FOXM1 (2305)",
  "0: median nonessential reference; -1: median common-essential reference. More negative = stronger growth disadvantage.",
  "Scores are normalized fitness effects, not viability percentages. No score-level significance or cancer specificity is inferred.",sep="\n")
base_theme<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),
  plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0,size=8.5))
save_plot<-function(p,name,width,height){
  ggsave(file.path(fig,paste0(name,".pdf")),p,width=width,height=height,bg="white",limitsize=FALSE)
  ggsave(file.path(fig,paste0(name,".png")),p,width=width,height=height,dpi=220,bg="white",limitsize=FALSE)
}
source(file.path(CODE_DIR,"scripts/10_dependency/utils/ranked_landscape.R"),local=new.env(parent=environment()))

p_dist<-ggplot(z,aes(x="LUAD",y=FOXM1_GeneEffect))+
  geom_hline(yintercept=0,color="grey40",linewidth=.5)+
  geom_hline(yintercept=-1,color="#B2182B",linetype=2,linewidth=.5)+
  geom_boxplot(width=.3,outlier.shape=NA,fill="#EAF0F6",color="#335775")+
  geom_point(position=position_jitter(width=.10,height=0,seed=2026),alpha=.8,size=2,color="#335775")+
  labs(title="FOXM1 gene-effect distribution across LUAD cell lines",
    subtitle=sprintf("n = %d; median %.3f; IQR %.3f to %.3f",nrow(z),median(values),quantile(values,.25),quantile(values,.75)),
    x=NULL,y="Chronos gene effect",caption=caption)+base_theme
save_plot(p_dist,"FOXM1_LUAD_GeneEffect_distribution",8,6.5)

ff<-featured[featured$OncotreeCode=="LUAD"&is.finite(featured$FOXM1_GeneEffect),,drop=FALSE]
if(nrow(ff)){
  ff<-ff[order(ff$FOXM1_GeneEffect),]
  ff$CellLineName<-factor(ff$CellLineName,levels=rev(ff$CellLineName))
  p_selected<-ggplot(ff,aes(FOXM1_GeneEffect,CellLineName))+
    geom_vline(xintercept=0,color="grey40")+
    geom_vline(xintercept=-1,linetype=2,color="#B2182B")+
    geom_segment(aes(x=0,xend=FOXM1_GeneEffect,yend=CellLineName),color="grey70")+
    geom_point(size=3.3,color="#335775")+
    geom_text(aes(label=sprintf("%.3f",FOXM1_GeneEffect)),nudge_y=.25,size=3.6)+
    labs(title="FOXM1 knockout effect in selected LUAD cell lines",
      x="Chronos gene effect",y=NULL,caption=caption)+base_theme
  save_plot(p_selected,"FOXM1_selected_LUAD_cell_lines",9,5.5)
}
saveRDS(list(release="DepMap 24Q4 Public v1",doi="10.25452/figshare.plus.27993248.v1",
  luad=luad,engineered_excluded=engineered,summary=summary_table,selected=featured),file.path(out,"FOXM1_LUAD_DepMap_analysis.rds"))
# Human-readable report with the actual estimates; do not translate into percent viability.
lines<-c("# FOXM1 CRISPR knockout fitness in LUAD cell lines",
  "", "Source: DepMap 24Q4 Public, version 1.",
  "DOI: https://doi.org/10.25452/figshare.plus.27993248.v1", "",
  sprintf("Initial LUAD / Cell Line annotations: %d models. Excluded %d explicitly engineered models; retained %d models, %d with FOXM1 gene-effect estimates and %d without scores.",n_annotated_luad,nrow(engineered),nrow(luad),nrow(z),nrow(luad)-nrow(z)),
  sprintf("Median %.4f; IQR %.4f to %.4f; range %.4f to %.4f.",median(values),quantile(values,.25),quantile(values,.75),min(values),max(values)),
  sprintf("%d/%d scores below zero; %d/%d <= -0.5; %d/%d <= -1. These cuts are descriptive, not dependency probabilities or significance calls.",sum(values<0),length(values),sum(values<=-.5),length(values),sum(values<=-1),length(values)),
  "", "## Selected cell lines", "", "| Cell line | ModelID | FOXM1 gene effect | Status |", "|---|---|---:|---|")
for(i in seq_len(nrow(featured)))lines<-c(lines,sprintf("| %s | %s | %s | %s |",featured$Requested_name[i],featured$ModelID[i],
  ifelse(is.finite(featured$FOXM1_GeneEffect[i]),sprintf("%.4f",featured$FOXM1_GeneEffect[i]),"NA"),featured$Status[i]))
lines<-c(lines,"","## Interpretation and limits","",
  "The assay estimates fitness from pooled CRISPR knockout screens. Negative scores are consistent with growth inhibition and/or cell loss after gene perturbation; proliferation and death are not separated here.",
  "The scale is normalized to median nonessential controls (0) and median common-essential controls (-1). -1 is NOT 100% cell killing; score magnitude is not percent viability.",
  "No per-cell-line replicate uncertainty or P values can be inferred from this matrix alone. Small negative values may not represent a meaningful dependency.",
  "These are model-level integrated estimates. Expression of FOXM1, LUAD prognostic association and FOXM1 knockout fitness are different quantities.",
  "No normal lung controls were included, so tumor selectivity, therapeutic window and effects on normal cells cannot be concluded.",
  "The release is pinned for reproducibility and is not the latest quarterly release. Missing scores do not mean no biological effect.",
  "Explicitly engineered models (nonempty EngineeredModel annotation) are excluded from the primary analysis and saved separately. Absence of this tag is not proof of an entirely unmodified genetic background.",
  "Do not infer that every LUAD model depends on FOXM1 or select a cell line solely from a single score. Follow-up assays need independent guides and rescue/orthogonal validation.",
  "","## Reproduce","","1. python3 Result_8/code/08_download_FOXM1_DepMap.py","2. Rscript Result_8/code/09_FOXM1_DepMap_LUAD_CRISPR.R")
writeLines(lines,file.path(out,"FOXM1_LUAD_report.md"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(summary_table);print(featured,row.names=FALSE)
message("Analysis complete. Results: ",out,"; figures: ",fig)
