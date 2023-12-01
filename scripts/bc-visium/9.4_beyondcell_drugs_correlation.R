rm(list = ls())

library("beyondcell")
library("Seurat")
library("tidyverse")
library("tidygraph")
library("patchwork")
library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")

bc.ranked <- read_rds("./results/analysis/bc_ranked_all.rds")

top.diff <- as.data.frame(bc.ranked@ranks) %>%
  select(starts_with(match = "TCs_res.0.3.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("TCs_res.0.3.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  pull("signature") %>%
  unique()

distrib.enrich <- bc.ranked@normalized[top.diff,] 
dim(distrib.enrich)
distrib.enrich <- t(distrib.enrich)
dim(distrib.enrich)
distrib.enrich <- as.data.frame(distrib.enrich) %>%
  rownames_to_column("spot")
distrib.metadata <- bc.ranked@meta.data %>%
  select(Cancer.Epithelial, TCs_res.0.3) %>%
  rownames_to_column("spot")
distrib.total <- right_join(distrib.metadata,distrib.enrich, by = "spot")
rownames(distrib.total) <- distrib.total$spot

distrib.total.pivot <- distrib.total %>%
  pivot_longer(cols = -c(spot,Cancer.Epithelial,TCs_res.0.3), names_to = "signatures", values_to = "enrichment")

library(ggstatsplot)

#To cite package 'ggstatsplot' in publications use:
#  
#  Patil, I. (2021). Visualizations with statistical details: The
#'ggstatsplot' approach. Journal of Open Source Software, 6(61), 3167,
#doi:10.21105/joss.03167

plot.correlation.list <- function(df,drugs,idents,cluster){
  plots <- lapply(drugs, function(drug){
    data <- df %>%
      filter(!!sym(idents) == cluster) %>%
      select(Cancer.Epithelial, drug) %>%
      rename(enrichment = drug) %>%
      ggscatterstats(
        data = ., ## data frame from which variables are taken
        x = Cancer.Epithelial, ## predictor/independent variable
        y = enrichment, ## dependent variable
        marginal = F,
        bf.message = F,
        ggplot.component = list(ylim(-75,75),
                                xlim(0,1),
                                ggtitle(drug)))
    return(data)
  })
  return(plots)
}

plots.TC3 <- plot.correlation.list(df = distrib.total, drugs = top.diff, idents = "TCs_res.0.3", cluster = "TC-3")

TCs <- unique(levels(bc.ranked@meta.data$TCs_res.0.3))

plots.TCs <- lapply(TCs, function(TC){
  plot.TC <- plot.correlation.list(df = distrib.total,
                                   drugs = top.diff,
                                   idents = "TCs_res.0.3",
                                   cluster = TC)
})
#stats.TC3 <- lapply(seq_along(drugs), function(index){
#  drug <- drugs[index]
#  df <- as.data.frame(x = drug)
#  stats <- extract_stats(plots.TC3[[index]])
#  stats <- stats$subtitle_data
#  df <- merge(df,stats)
#  return(df)
#})
#stats.TC3.df <- do.call(rbind,stats.TC3)

stats.function <- function(drugs,plots.list,cluster){
  stats.cluster <- lapply(seq_along(drugs), function(index){
    drug <- drugs[index]
    df <- as.data.frame(x = drug)
    stats <- extract_stats(plots.list[[index]])
    stats <- stats$subtitle_data
    df <- merge(df,stats)
    return(df)
  })
  stats.cluster.df <- do.call(rbind,stats.cluster)
  stats.cluster.df$TC <- cluster
  return(stats.cluster.df)
}

stats.TC3 <- stats.function(drugs = top.diff, plots.list = plots.TC3, cluster = "TC-3")

stats.TCs <- lapply(seq_along(TCs), function(TC){
  plots <- plots.TCs[[TC]]
  stats <- stats.function(drugs = top.diff,
                          plots.list = plots,
                          cluster = TCs[TC])
  return(stats)
})

stats.TCs.global <- do.call(rbind,stats.TCs)

stats.TC1.df.split <- stats.TCs.global %>%
  filter(TC == "TC-1") %>%
  select(drug,estimate,p.value) %>%
  rename(estimate.TC1 = estimate,
         p.value.TC1 = p.value)

stats.split.list <- lapply(seq_along(TCs), function(index){
  split <- stats.TCs.global %>%
    filter(TC == TCs[index]) %>%
    select(drug, estimate, p.value) %>%
    rename_with(.fn = ~paste0("estimate.",TCs[index]), .cols = estimate) %>%
    rename_with(.fn = ~paste0("p.value.",TCs[index]), .cols = p.value) %>%
    mutate(drug = NULL)
})


stats.split.merged <- do.call(cbind,stats.split.list)
rownames(stats.split.merged) <- top.diff


matrix.cor <- stats.split.merged %>%
  select(starts_with("estimate")) %>%
  rename_with(.fn = ~gsub(pattern = "estimate.",replacement = "",.x)) %>%
  as.matrix()


pvalue.TC1 = stats.TC1.df.split$p.value.TC1
is_sig = pvalue.TC1 < 0.01
pch = rep("*", length(pvalue.TC1))
pch[!is_sig] = NA
pvalue_col_fun = colorRamp2(c(0, 2, 3), c("gree", "white", "red")) 

TC.colors <- c("TC-1" = "#00b2d7",
               "TC-2" = "#e5c22f",
               "TC-3" = "#903ca2",
               "TC-4" = "#3f8741",
               "TC-5" = "#ff7b00",
               "TC-6" = "#cb5c42")

ht.TC1 <- Heatmap(
  matrix = matrix.cor,
  name = "Corr Cancer-TC1",
  #row_labels = stats.TC1.df.split$drug,
  cluster_columns = F,
  show_column_names = F,
  cluster_rows = F,
  #row_labels = toupper(stats.global$preferred.drug.names),
  top_annotation = HeatmapAnnotation(TC = TCs,
                                     col = list("TC" = TC.colors)),
  #right_annotation = rowAnnotation(
  #  #pvalue = anno_simple(-log10(pvalue), col = pvalue_col_fun, pch = pch),
  #  corr = anno_numeric(round(stats.TC1.df.split$estimate.TC1,2), 
  #                         align_to = "right",
  #                         bg_gp = gpar(fill = "white", col = "white"),
  #                         labels_gp = NULL,
  #                         labels_offset = unit(0, "cm"),
  #                         width = unit(1.25, "cm")), 
  #  pvalue = anno_numeric(round((stats.TC1.df.split$p.value.TC1),3), 
  #                         align_to = "right",
  #                         bg_gp = gpar(fill = "white", col = "white"),
  #                         labels_gp = NULL,
  #                         labels_offset = unit(0, "cm"),
  #                         width = unit(1.25, "cm")),
  #  pvalue.sig = anno_simple(-log10(pvalue.TC1), 
  #                       col = pvalue_col_fun, 
  #                       pch = pch)),
  #col = colorRamp2(c(min, 0, max), c("blue", "white", "red")),
  #heatmap_legend_param = list(at = c(min, 0, max))
)
ht.TC1
