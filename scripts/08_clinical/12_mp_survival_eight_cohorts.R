# Same MP state across eight cohorts: 9 panels in a 3 x 3 combined forest.
# ssGSEA and univariable Cox match 06-MP_prognosis.R, with explicit sample-ID checks.
suppressPackageStartupMessages({library(GSVA);library(BiocParallel);library(survival);library(ggplot2)})
source("config.R")
root <- BASE_DIR
work <- file.path(root,"Result_2/results_data/MP_univariable_Cox_by_state")
out <- file.path(root,"Result_8/results_figure")
for(d in c(work,out)) dir.create(d,recursive=TRUE,showWarnings=FALSE)
cohort_order <- c("TCGA","GSE11969","GSE13213","GSE26939","GSE30219","GSE31210","GSE41271","GSE72094")
mp_order <- c("G2M","G1S","NFkB","MHC-II","Translation","OXPHOS","Adhesion","EMT","Glycolysis")
mp_ids <- c("MP_1","MP_8","MP_2","MP_7","MP_4","MP_6","MP_5","MP_3","MP_9")
e <- new.env();load(file.path(root,"Result_2/results_data/02_Malig_final_MP_top50.RData"),e)
gene_sets <- e$MP_list[mp_ids];names(gene_sets)<-mp_order
files <- setNames(file.path(root,"bulk/sur_vaild_data",ifelse(cohort_order=="TCGA","TCGA_sur.Rdata",paste0(cohort_order,"_sur_data.Rdata"))),cohort_order)
manifest <- list(files=tools::md5sum(c(files,file.path(root,"bulk/TCGA/TCGA_sur_data.Rdata"))),
                 gene_sets=gene_sets,GSVA=as.character(packageVersion("GSVA")))
cache <- file.path(work,"scores_and_survival.rds")
if(file.exists(cache)) {
  saved <- readRDS(cache)
  stopifnot(identical(saved$manifest,manifest))
  inputs <- saved$inputs
} else {
  inputs <- list()
  for(ds in cohort_order) {
    message("Computing ssGSEA: ",ds)
    e <- new.env();load(files[[ds]],e)
    expr <- as.matrix(e$GSE_expr);cl <- as.data.frame(e$clinical)
    ids <- colnames(expr)
    stopifnot(!anyDuplicated(ids),!is.null(ids))
    match_note <- "Clinical rownames matched to expression sample IDs"
    if(ds=="TCGA") {
      orig <- new.env();load(file.path(root,"bulk/TCGA/TCGA_sur_data.Rdata"),orig)
      ot <- ifelse(as.character(orig$clinical$vital_status)=="Dead",
        as.numeric(orig$clinical$days_to_death),as.numeric(orig$clinical$days_to_last_follow_up))
      stopifnot(identical(ids,as.character(orig$clinical$samples_id)),
        isTRUE(all.equal(expr,as.matrix(orig$mRNA_exp),check.attributes=FALSE)),
        identical(as.character(cl$vital_status),as.character(orig$clinical$event)),
        isTRUE(all.equal(as.numeric(cl$survival_time),ot)))
      rownames(cl)<-orig$clinical$samples_id
      match_note <- "TCGA original samples_id; expression, event, survival verified before matching"
      rm(orig)
    }
    ix <- match(ids,rownames(cl));stopifnot(!anyNA(ix))
    cl <- cl[ix,,drop=FALSE];stopifnot(identical(rownames(cl),ids))
    raw_event <- as.character(cl$vital_status)
    if(all(na.omit(raw_event)%in%c("0","1"))) event<-as.integer(raw_event) else
      if(all(na.omit(raw_event)%in%c("1","2"))) event<-as.integer(raw_event)-1L else stop(ds,": unknown event coding")
    score <- GSVA::gsva(GSVA::ssgseaParam(exprData=expr,geneSets=gene_sets),
                        BPPARAM=BiocParallel::SerialParam(),verbose=FALSE)
    stopifnot(identical(colnames(score),ids),setequal(rownames(score),mp_order))
    inputs[[ds]] <- list(score=score[mp_order,,drop=FALSE],
      survival=data.frame(sample=ids,time=as.numeric(as.character(cl$survival_time)),event=event,row.names=ids),
      matching=match_note)
    rm(e,expr,cl,score);gc()
  }
  saveRDS(list(manifest=manifest,inputs=inputs),cache)
}



# Each score is fitted separately, matching the univariable analysis in script 06.
results <- list(); models <- list(); k <- 0L
for(ds in cohort_order) {
  z <- inputs[[ds]]$score
  for(mp in mp_order) {
    values <- as.numeric(z[mp,])
    stopifnot(all(is.finite(values)),sd(values)>0)
    dat <- inputs[[ds]]$survival
    dat$MP_score <- as.numeric(scale(values))
    dat <- dat[is.finite(dat$time)&dat$time>0&!is.na(dat$event),,drop=FALSE]
    stopifnot(sum(dat$event==1)>=5)
    fit <- coxph(Surv(time,event)~MP_score,data=dat,ties="efron")
    sm <- summary(fit)
    k <- k+1L
    results[[k]] <- data.frame(Dataset=ds,MP=mp,N=nrow(dat),Events=sum(dat$event==1),
      HR=as.numeric(sm$conf.int[1,"exp(coef)"]),CI_lower=sm$conf.int[1,"lower .95"],
      CI_upper=sm$conf.int[1,"upper .95"],P_value=sm$coefficients[1,"Pr(>|z|)"])
    models[[paste(ds,mp,sep="__")]] <- fit
  }
  message("Univariable Cox complete: ",ds)
}
res <- do.call(rbind,results);rownames(res)<-NULL
res$FDR <- p.adjust(res$P_value,method="BH")
write.csv(res,file.path(work,"All_9States_8Cohorts_univariable_Cox.csv"),row.names=FALSE)
saveRDS(list(results=res,models=models,matching=lapply(inputs,`[[`,"matching")),
        file.path(work,"All_9States_8Cohorts_univariable_Cox.rds"))
df <- as.data.frame(res)

# ============================================================
# 2. DATA VALIDATION
# ============================================================

required_cols <- c(
  "Dataset",
  "MP",
  "HR",
  "CI_lower",
  "CI_upper",
  "P_value"
)

missing_cols <- setdiff(
  required_cols,
  colnames(df)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

numeric_cols <- c(
  "HR",
  "CI_lower",
  "CI_upper",
  "P_value"
)

for (nm in numeric_cols) {
  df[[nm]] <- as.numeric(
    as.character(df[[nm]])
  )
}

if (
  any(!is.finite(as.matrix(df[, numeric_cols]))) ||
  any(df$HR <= 0) ||
  any(df$CI_lower <= 0) ||
  any(df$CI_upper < df$CI_lower) ||
  any(df$HR < df$CI_lower) ||
  any(df$HR > df$CI_upper) ||
  any(df$P_value < 0 | df$P_value > 1)
) {
  stop("Invalid Cox results.")
}

# ============================================================
# 3. DATASET ORDER
# ============================================================

cohort_order <- c(
  "TCGA",
  "GSE11969",
  "GSE13213",
  "GSE26939",
  "GSE30219",
  "GSE31210",
  "GSE41271",
  "GSE72094"
)

mp_order <- c(
  "G2M",
  "G1S",
  "NFkB",
  "MHC-II",
  "Translation",
  "OXPHOS",
  "Adhesion",
  "EMT",
  "Glycolysis"
)

mp_labels <- c(
  "G2/M state",
  "G1/S state",
  "NF-kB state",
  "MHC-II state",
  "Translation state",
  "OXPHOS state",
  "Adhesion state",
  "EMT state",
  "Glycolysis state"
)

if (!setequal(unique(as.character(df$Dataset)),
              cohort_order)) {
  stop("Unexpected cohort names.")
}

if (!setequal(unique(as.character(df$MP)),
              mp_order)) {
  stop("Unexpected MP names.")
}

if (nrow(df) != 72) {
  stop("Expected exactly 72 Cox results.")
}

if (anyDuplicated(df[, c("Dataset", "MP")])) {
  stop("Duplicated Dataset-MP combinations.")
}

df$Dataset <- factor(
  as.character(df$Dataset),
  levels = cohort_order
)

df$MP <- as.character(df$MP)

# Each panel is an MP state; rows are eight cohorts.
df$MP_label <- as.character(df$Dataset)
df$y <- 9 - match(as.character(df$Dataset), cohort_order)
df$MP <- factor(df$MP,levels=mp_order)

# ============================================================
# 4. FORMAT PVALUES
# ============================================================

format_pvalue <- function(p) {

  vapply(
    p,
    function(x) {

      if (x < 0.0001) {

        formatC(
          x,
          format = "e",
          digits = 2
        )

      } else if (x < 0.001) {

        sprintf("%.4f", x)

      } else {

        sprintf("%.3f", x)
      }

    },
    FUN.VALUE = character(1)
  )
}

df$P_label <- format_pvalue(
  df$P_value
)

# ============================================================
# 5. STATE COLORS
# ============================================================

# Colors corresponding to nine MP states

mp_colors <- c(
  "G2M"         = "#D7191C",
  "G1S"         = "#E34A33",
  "NFkB"        = "#279B8F",
  "MHC-II"      = "#40A6C0",
  "Translation" = "#435D8D",
  "OXPHOS"      = "#949AB9",
  "Adhesion"    = "#E69A85",
  "EMT"         = "#806B5B",
  "Glycolysis"  = "#8CC5B5"
)

# Reference figure:
# Blue CI + red square
CI_COLOR <- "#497FA8"
HR_COLOR <- "#CF3034"

# ============================================================
# 6. FIGURE GEOMETRY
# ============================================================

# The complete plot uses a numerical coordinate system.
# The left part is reserved for MP and Pvalue.
# The right part is reserved for HR and CI.

MP_X <- 0.08
P_X <- 1.45

SEPARATOR_X <- 1.95

# HR = 1 position
HR_CENTER <- 3.22

# Log2 transform allows equivalent spacing
# between HR = 0.5, 1 and 2.
HR_SCALE <- 0.95

hr_to_x <- function(hr) {
  HR_CENTER + HR_SCALE * log2(hr)
}

# HR axis breaks
HR_BREAKS <- c(
  0.5,
  1,
  1.5,
  2,
  3
)

# This range covers your current 72 results.
if (
  min(df$CI_lower) < 0.4 ||
  max(df$CI_upper) > 3
) {
  stop(
    "Some confidence intervals fall outside ",
    "the preset HR axis (0.4 to 3). ",
    "Please expand HR_BREAKS and plot limits."
  )
}

# ============================================================
# 7. FOREST PLOT FUNCTION
# ============================================================

make_forest <- function(
    plot_data,
    combined = FALSE,
    cohort_name = NULL
) {

  p <- ggplot(
    plot_data
  ) +

    # --------------------------------------------------------
  # Column separator
  # --------------------------------------------------------

  annotate(
    "segment",
    x = SEPARATOR_X,
    xend = SEPARATOR_X,
    y = 0.73,
    yend = 8.43,
    color = "#D0D0D0",
    linewidth = 0.55
  ) +

    # --------------------------------------------------------
  # HR = 1 reference line
  # --------------------------------------------------------

  annotate(
    "segment",
    x = hr_to_x(1),
    xend = hr_to_x(1),
    y = 0.70,
    yend = 8.43,
    color = "#B0B0B0",
    linewidth = 0.65,
    linetype = "dashed"
  ) +

    # --------------------------------------------------------
  # MP state names
  # --------------------------------------------------------

  geom_text(
    data = plot_data[as.character(plot_data$MP) %in% mp_order[c(1,4,7)], , drop = FALSE],
    aes(
      y = y,
      label = MP_label
    ),
    color = "black",
    x = MP_X,
    hjust = 0,
    size = 3.3,
    show.legend = FALSE
  ) +

    # --------------------------------------------------------
  # Pvalue column
  # --------------------------------------------------------

  geom_text(
    aes(
      y = y,
      label = P_label
    ),
    x = P_X,
    hjust = 0.5,
    size = 3.2,
    color = "black"
  ) +

    # --------------------------------------------------------
  # 95% confidence interval
  # --------------------------------------------------------

  geom_segment(
    aes(
      x = hr_to_x(CI_lower),
      xend = hr_to_x(CI_upper),
      y = y,
      yend = y
    ),
    color = CI_COLOR,
    linewidth = 1.25,
    lineend = "round"
  ) +

    # --------------------------------------------------------
  # CI left cap
  # --------------------------------------------------------

  geom_segment(
    aes(
      x = hr_to_x(CI_lower),
      xend = hr_to_x(CI_lower),
      y = y - 0.11,
      yend = y + 0.11
    ),
    color = CI_COLOR,
    linewidth = 0.65
  ) +

    # --------------------------------------------------------
  # CI right cap
  # --------------------------------------------------------

  geom_segment(
    aes(
      x = hr_to_x(CI_upper),
      xend = hr_to_x(CI_upper),
      y = y - 0.11,
      yend = y + 0.11
    ),
    color = CI_COLOR,
    linewidth = 0.65
  ) +

    # --------------------------------------------------------
  # HR point: red square
  # --------------------------------------------------------

  geom_point(
    aes(
      x = hr_to_x(HR),
      y = y
    ),
    shape = 15,
    size = 3.8,
    color = HR_COLOR
  ) +

    # --------------------------------------------------------
  # MP header
  # --------------------------------------------------------

  geom_text(
    data = unique(plot_data[as.character(plot_data$MP) %in% mp_order[c(1,4,7)], "MP", drop = FALSE]),
    x = MP_X,
    y = 8.70,
    label = "Dataset",
    hjust = 0,
    fontface = "bold",
    size = 3.6,
    color = "black"
  ) +

    # --------------------------------------------------------
  # Pvalue header
  # --------------------------------------------------------

  annotate(
    "text",
    x = P_X,
    y = 8.70,
    label = "Pvalue",
    fontface = "bold",
    size = 3.6,
    color = "black"
  ) +

    # --------------------------------------------------------
  # HR-axis baseline
  # --------------------------------------------------------

  annotate(
    "segment",
    x = hr_to_x(min(HR_BREAKS)),
    xend = hr_to_x(max(HR_BREAKS)),
    y = 0.62,
    yend = 0.62,
    linewidth = 0.65,
    color = "black"
  ) +

    # --------------------------------------------------------
  # HR-axis tick marks
  # --------------------------------------------------------

  geom_segment(
    data = data.frame(
      tick_x = hr_to_x(HR_BREAKS)
    ),
    aes(
      x = tick_x,
      xend = tick_x,
      y = 0.62,
      yend = 0.72
    ),
    inherit.aes = FALSE,
    linewidth = 0.65,
    color = "black"
  ) +

    # --------------------------------------------------------
  # HR-axis numerical labels
  # --------------------------------------------------------

  geom_text(
    data = data.frame(
      tick_x = hr_to_x(HR_BREAKS),
      tick_label = as.character(HR_BREAKS)
    ),
    aes(
      x = tick_x,
      y = 0.34,
      label = tick_label
    ),
    inherit.aes = FALSE,
    size = 3.0,
    color = "black"
  ) +

    # --------------------------------------------------------
  # HR axis title
  # --------------------------------------------------------

  annotate(
    "text",
    x = hr_to_x(1.4),
    y = 0.08,
    label = "HR",
    size = 3.4,
    color = "black"
  ) +

    # --------------------------------------------------------
  # State colors
  # --------------------------------------------------------

  scale_color_manual(
    values = mp_colors,
    guide = "none"
  ) +

    # --------------------------------------------------------
  # Plot coordinates
  # --------------------------------------------------------

  coord_cartesian(
    xlim = c(0, 5.00),
    ylim = c(0, 9.05),
    expand = FALSE,
    clip = "off"
  ) +

    theme_void(
      base_size = 11
    ) +

    theme(
      plot.background = element_rect(
        fill = "white",
        color = NA
      ),

      panel.background = element_rect(
        fill = "white",
        color = NA
      ),

      panel.spacing = grid::unit(
        1.4,
        "lines"
      ),

      plot.margin = margin(
        10, 12, 10, 12
      ),

      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14,
        margin = margin(b = 6)
      ),

      strip.text = element_text(
        face = "bold",
        size = 13,
        color = "black"
      ),

      strip.background = element_rect(
        fill = "#F2F2F2",
        color = NA
      )
    )

  # ----------------------------------------------------------
  # Combined plot
  # ----------------------------------------------------------

  if (combined) {

    p <- p +

      facet_wrap(
        ~ MP,
        ncol = 3,
        labeller = as_labeller(setNames(paste0("MP",seq_along(mp_order)," ",mp_order),mp_order)),
        drop = FALSE
      ) +

      labs(
        title = paste(
          "Univariate Cox analysis: each MP state across eight cohorts"
        ),

        caption = paste(
          "Red squares: HR;",
          "blue lines: 95% CI;",
          "Pvalue: unadjusted Cox Wald P value;",
          "HR per 1-SD increase in MP score."
        )
      ) +

      theme(
        plot.title = element_text(
          face = "bold",
          size = 17,
          hjust = 0.5
        ),

        plot.caption = element_text(
          size = 10,
          hjust = 0.5,
          margin = margin(t = 12)
        )
      )

  } else {

    p <- p +

      labs(
        title = cohort_name
      )

  }

  return(p)
}

# ============================================================

p_combined <- make_forest(df,combined=TRUE)
# Color facet titles individually using native grid text grobs.
# All dataset labels are black; only first-column panels show them.
title_colors <- setNames(unname(mp_colors[mp_order]),paste0("MP",seq_along(mp_order)," ",mp_order))
color_title_grob <- function(g) {
  if (inherits(g,"text") && length(g$label)==1L &&
      as.character(g$label) %in% names(title_colors)) {
    g$gp$col <- unname(title_colors[as.character(g$label)])
  }
  if (!is.null(g$children)) {
    for (i in seq_along(g$children)) g$children[[i]] <- color_title_grob(g$children[[i]])
  }
  if (!is.null(g$grobs)) g$grobs <- lapply(g$grobs,color_title_grob)
  g
}
p_combined <- color_title_grob(ggplotGrob(p_combined))
pdf_file <- file.path(out,"All_9States_MP_Cox_forest_across_8Cohorts.pdf")
ggsave(pdf_file,p_combined,width=21,height=15,bg="white")
ggsave(sub("pdf$","png",pdf_file),p_combined,width=21,height=15,dpi=180,bg="white")
capture.output(sessionInfo(),file=file.path(work,"sessionInfo.txt"))
message("Saved: ",pdf_file)
