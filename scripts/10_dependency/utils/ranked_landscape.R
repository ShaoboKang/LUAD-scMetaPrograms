# Landscape, publication-style rendering of saved FOXM1 DepMap results.
# Only reorders and redraws existing estimates; no data filtering or recomputation.
suppressPackageStartupMessages(library(ggplot2))
root<-path.expand(Sys.getenv("LUAD_ROOT","~/LUAD_project"))
dat<-read.csv(file.path(root,"Result_8/results_data/FOXM1_DepMap/FOXM1_LUAD_GeneEffect_ranked.csv"),check.names=FALSE)
stopifnot(!anyDuplicated(dat$ModelID),all(is.finite(dat$FOXM1_GeneEffect)))
dat<-dat[order(dat$FOXM1_GeneEffect,dat$ModelID),]
stopifnot(!anyDuplicated(dat$CellLineName))
dat$CellLine<-factor(dat$CellLineName,levels=dat$CellLineName)
# Muted blue stems and darker blue points; no untested significance categories.
p<-ggplot(dat,aes(x=CellLine,y=FOXM1_GeneEffect))+
  geom_hline(yintercept=0,color="#424242",linewidth=.5)+
  geom_hline(yintercept=-1,color="#999999",linetype="dashed",linewidth=.45)+
  geom_segment(aes(xend=CellLine,y=0,yend=FOXM1_GeneEffect),
               color="#9DB6C9",linewidth=.75,lineend="round")+
  geom_point(shape=21,size=2.9,stroke=.35,fill="#32678C",color="#234D68")+
  geom_text(aes(label=sprintf("%.3f",FOXM1_GeneEffect)),
            nudge_y=-.025,angle=90,hjust=1,vjust=.5,size=4,color="#333333")+
  scale_y_continuous(breaks=c(-1.25,-1,-.75,-.5,-.25,0,.25),
    expand=expansion(mult=c(.10,.10)))+
  scale_x_discrete(expand=expansion(add=.7))+
  labs(title="FOXM1 knockout fitness effect in LUAD cell lines",
       x=NULL,y="Chronos gene effect")+
  theme_classic(base_size=17,base_family="sans")+
  theme(plot.title=element_text(face="bold",size=22,hjust=.5,margin=margin(b=15)),
        axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=13,color="#333333"),
        axis.text.y=element_text(size=16,color="#333333"),
        axis.title.y=element_text(size=18,margin=margin(r=10)),
        axis.line=element_line(color="#424242",linewidth=.45),
        axis.ticks=element_line(color="#424242",linewidth=.4),
        plot.margin=margin(14,18,12,14),legend.position="none")
figdir<-file.path(root,"Result_8/results_figure/FOXM1_DepMap")
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(figdir,"FOXM1_LUAD_GeneEffect_ranked.pdf"),p,width=18,height=6.5,bg="white")
ggsave(file.path(figdir,"FOXM1_LUAD_GeneEffect_ranked.png"),p,width=18,height=6.5,dpi=300,bg="white")
message("Landscape figure updated: ",figdir)
