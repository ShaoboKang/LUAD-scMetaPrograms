# RS + clinical features: eight multivariable Cox models, one combined forest.
# Only code and final PDF/PNG are saved; fits/tables stay in memory.
# Endpoint source is the clinical RDS (not the RS file). User previously confirmed
# TCGA OS TRUE=death/FALSE=censored, despite disagreement with older endpoints.
# Raw RS is continuous: HR per one original RS unit (no dichotomization/rescaling).
# Stage=1..4, T=1..4, N=0..3 and M=0/1 are single numeric terms.
# Their coefficients assume a constant log-hazard increment per stage/category step.
# TX/NX/MX remain missing; sex female->male and EGFR WT->MUT are binary factors.
suppressPackageStartupMessages({library(survival);library(ggplot2)})
source("config.R")
root <- BASE_DIR
figdir<-file.path(root,"Result_8/results_figure")
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
risk<-readRDS(file.path(root,"Result_8/results_data/Model_riskscore.rds"))
clinical<-readRDS(CLINICAL_FILE)
cohorts<-c("TCGA","GSE30219","GSE72094","GSE26939","GSE41271","GSE13213","GSE31210","GSE11969")
# Original AJCC labels are used only to verify T/N/M code meanings, not endpoints.
original<-new.env();load(file.path(root,"bulk/TCGA/TCGA_sur_data.Rdata"),original)
key<-function(x,ds){x<-toupper(trimws(as.character(x)));if(ds=="TCGA")substr(x,1,12) else x}
fits<-list();results<-list();audits<-list();ph_checks<-list()
for(ds in cohorts){
  d<-clinical[[paste0(if(ds=="TCGA")"tcga" else ds,"_cps")]]
  d<-d[,!grepl("^(CPS|CSS)($|_)",names(d),ignore.case=TRUE),drop=FALSE]
  rs<-risk[[ds]]
  stopifnot(all(c("ID","RS")%in%names(rs)))
  ids<-key(d$ID,ds);rids<-key(rs$ID,ds)
  stopifnot(!anyNA(ids),!anyNA(rids),!anyDuplicated(ids),!anyDuplicated(rids))
  ix<-match(ids,rids)
  d$RS<-suppressWarnings(as.numeric(as.character(rs$RS[ix])))
  d$time_days<-as.numeric(d$OS.time)*if(ds=="GSE30219")30 else 1
  d$event<-if(is.logical(d$OS))as.integer(d$OS) else as.integer(as.character(d$OS))
  stopifnot(all(na.omit(d$event)%in%c(0,1)))
  sex<-tolower(as.character(d$Sex));sex[sex=="f"]<-"female";sex[sex=="m"]<-"male"
  stopifnot(all(na.omit(sex)%in%c("female","male")))
  d$Sex<-factor(sex,levels=c("female","male"))
  covars<-character()
  if("Stage"%in%names(d)){
    stopifnot(all(na.omit(as.character(d$Stage))%in%as.character(1:4)))
    d$Stage<-as.numeric(as.character(d$Stage))
    covars<-c(covars,"Stage")
  }
  if(ds=="TCGA"){
    old<-original$clinical;oldids<-key(old$case_submitter_id,ds)
    stopifnot(!anyDuplicated(oldids));j<-match(ids,oldids);stopifnot(!anyNA(j))
    for(v in c("T","N","M")){
      labels<-sub("[abc]$","",as.character(old[[paste0("ajcc_pathologic_",tolower(v))]][j]))
      for(code in unique(na.omit(d[[v]])))stopifnot(length(unique(labels[!is.na(d[[v]])&d[[v]]==code]))==1L)
      labels[labels%in%c("TX","NX","MX")]<-NA_character_
      d[[v]]<-as.numeric(sub(paste0("^",v),"",labels))
    }
    covars<-c(covars,"T","N","M")
  }
  egfr<-intersect(c("EGFR","EGFR_status"),names(d))
  if(length(egfr)){
    vals<-toupper(as.character(d[[egfr[1]]]))
    stopifnot(all(na.omit(vals)%in%c("WT","MUT")))
    d$EGFR<-factor(vals,levels=c("WT","MUT"));covars<-c(covars,"EGFR")
  }
  if("Age"%in%names(d)){d$Age<-as.numeric(d$Age);covars<-c(covars,"Age")}
  covars<-c(covars,"Sex")
  has_rs<-is.finite(d$RS)
  valid_surv<-is.finite(d$time_days)&d$time_days>0&!is.na(d$event)
  complete<-complete.cases(d[,covars,drop=FALSE])
  dd<-d[has_rs&valid_surv&complete,,drop=FALSE]
  stopifnot(nrow(dd)>0,sum(dd$event==1)>0,sd(dd$RS)>0)
  for(v in covars)if(is.factor(dd[[v]]))dd[[v]]<-droplevels(dd[[v]])
  constant<-covars[vapply(dd[,covars,drop=FALSE],function(z)length(unique(z))<2L,logical(1))]
  adjusted<-setdiff(covars,constant)
  fmla<-reformulate(c("RS",adjusted),response="Surv(time_days,event)")
  warn<-character()
  fit<-withCallingHandlers(coxph(fmla,data=dd,ties="efron",x=TRUE,y=TRUE,
      control=coxph.control(iter.max=100)),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart("muffleWarning")})
  sm<-summary(fit)
  z<-data.frame(Dataset=ds,Term=rownames(sm$coefficients),
    HR=sm$conf.int[,"exp(coef)"],Lower=sm$conf.int[,"lower .95"],Upper=sm$conf.int[,"upper .95"],
    P=sm$coefficients[,"Pr(>|z|)"],row.names=NULL)
  z$Label<-z$Term
  # Show short variable names; preserve the modeled contrasts in the code.
  for(v in adjusted)if(is.factor(dd[[v]])){
    i<-startsWith(z$Term,v)
    z$Label[i]<-v
  }
  z$Estimable<-is.finite(z$HR)&is.finite(z$Lower)&is.finite(z$Upper)&z$Lower>0&is.finite(z$P)
  fits[[ds]]<-fit;results[[ds]]<-z
  ph_checks[[ds]]<-tryCatch(cox.zph(fit),error=function(e)conditionMessage(e))
  audits[[ds]]<-data.frame(Dataset=ds,Clinical_N=nrow(d),RS_matched=sum(has_rs),
    No_RS=sum(!has_rs),Invalid_survival=sum(has_rs&!valid_surv),
    Missing_clinical=sum(has_rs&valid_surv&!complete),Model_N=nrow(dd),
    Events=sum(dd$event==1),Censored=sum(dd$event==0),
    Covariates=paste(adjusted,collapse=" + "),Constant_omitted=paste(constant,collapse=", "),
    Warning=paste(warn,collapse="; "))
}
cox_results<-do.call(rbind,results);rownames(cox_results)<-NULL
sample_audit<-do.call(rbind,audits);rownames(sample_audit)<-NULL
print(sample_audit,row.names=FALSE)
cat("\nRS results (HR per 1 raw RS unit):\n")
print(cox_results[cox_results$Term=="RS",],row.names=FALSE)
if(any(!cox_results$Estimable)){
  cat("\nNot independently estimable terms (NE in figure):\n")
  print(cox_results[!cox_results$Estimable,c("Dataset","Label")],row.names=FALSE)
}
# Forest uses actual HR on a logarithmic axis (no negative HR values).
# Exact intervals are printed; arrows indicate values outside the display range.
axis_low<-.125;axis_high<-32
hx<-function(v)4.35+3.15*(log(v)-log(axis_low))/(log(axis_high)-log(axis_low))
plot_cohort<-function(ds){
  z<-results[[ds]];a<-audits[[ds]];row_count<-nrow(z)
  n<-max(vapply(results,nrow,integer(1)))
  bottom<-n-row_count
  z$Y<-n+1-seq_len(row_count)
  z$L<-pmax(axis_low,pmin(axis_high,z$Lower));z$U<-pmax(axis_low,pmin(axis_high,z$Upper))
  z$H<-pmax(axis_low,pmin(axis_high,z$HR))
  z$HRtext<-ifelse(z$Estimable,sprintf("%.2f (%.2f-%.2f)",z$HR,z$Lower,z$Upper),"Not estimable")
  z$Ptext<-ifelse(!z$Estimable,"NE",ifelse(z$P<.001,"<0.001",sprintf("%.3f",z$P)))
  z$Color<-ifelse(z$Term=="RS","#B2182B","black")
  ticks<-data.frame(HR=c(.25,1,2,4,8,16,32))
  ggplot(z)+
    annotate("segment",x=0,xend=8.75,y=n+.5,yend=n+.5,linewidth=.5)+
    annotate("text",x=0,y=n+1.05,label="Variable",hjust=0,fontface="bold",size=4.6)+
    annotate("text",x=3,y=n+1.05,label="HR (95% CI)",fontface="bold",size=4.6)+
    annotate("text",x=8.3,y=n+1.05,label="P value",fontface="bold",size=4.6)+
    annotate("segment",x=hx(1),xend=hx(1),y=bottom+.6,yend=n+.4,linetype=2,color="grey60")+
    geom_text(aes(x=0,y=Y,label=Label,color=Color),hjust=0,size=4.4)+
    geom_text(aes(x=3,y=Y,label=HRtext),size=4.4)+
    geom_text(aes(x=8.3,y=Y,label=Ptext),size=4.4)+
    geom_segment(data=z[z$Estimable,],aes(x=hx(L),xend=hx(U),y=Y,yend=Y),color="#497FA8",linewidth=.85)+
    geom_point(data=z[z$Estimable,],aes(x=hx(H),y=Y),shape=15,color="#CF3034",size=2.8)+
    geom_point(data=z[z$Estimable&z$Lower<axis_low,],aes(x=hx(axis_low),y=Y),shape=60,color="#497FA8",size=3)+
    geom_point(data=z[z$Estimable&z$Upper>axis_high,],aes(x=hx(axis_high),y=Y),shape=62,color="#497FA8",size=3)+
    annotate("segment",x=hx(axis_low),xend=hx(axis_high),y=bottom+.4,yend=bottom+.4,linewidth=.5)+
    geom_text(data=ticks,aes(x=hx(HR),y=bottom-.1,label=HR),inherit.aes=FALSE,size=3.8)+
    annotate("text",x=hx(1),y=bottom-.65,label="Hazard ratio",size=4.4)+
    scale_color_identity()+coord_cartesian(xlim=c(-.05,8.75),ylim=c(-1,n+1.65),expand=FALSE,clip="off")+
    labs(title=ds)+
    theme_void()+theme(plot.title=element_text(face="bold",hjust=0,size=18),plot.margin=margin(8,14,10,14))
}
plots<-lapply(cohorts,plot_cohort)
figure<-gtable::gtable_matrix("RS_clinical_forests",grobs=matrix(lapply(plots,ggplotGrob),nrow=2,ncol=4,byrow=TRUE),
  widths=grid::unit(rep(1,4),"null"),heights=grid::unit(rep(1,2),"null"))
figure<-gtable::gtable_add_rows(figure,grid::unit(.45,"in"),pos=0)
figure<-gtable::gtable_add_grob(figure,grid::textGrob("Multivariable Cox: risk score and clinical covariates",
  gp=grid::gpar(fontsize=24,fontface="bold")),t=1,l=1,r=4)
# No bottom caption, following the user's established figure preference.
pdf_file<-file.path(figdir,"All_8Cohorts_RS_clinical_multivariable_Cox.pdf")
ggsave(pdf_file,figure,width=28,height=8.5,bg="white")
ggsave(sub("pdf$","png",pdf_file),figure,width=28,height=8.5,dpi=180,bg="white")
stopifnot(length(fits)==8L,sum(cox_results$Term=="RS")==8L,
          all(cox_results$Estimable[cox_results$Term=="RS"]))
message("Final figure saved: ",pdf_file)
# Objects fits, cox_results, sample_audit and ph_checks remain in memory when sourced.
