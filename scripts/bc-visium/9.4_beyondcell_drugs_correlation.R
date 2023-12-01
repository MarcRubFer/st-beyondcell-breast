rm(list = ls())

library("beyondcell")
library("Seurat")
library("tidyverse")
library("tidygraph")
library("patchwork")
library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")

bc.ranked <- readRDS("./results/analysis/bc_ranked_all.rds")

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

#plots.TC3 <- plot.correlation.list(df = distrib.total, drugs = top.diff, idents = "TCs_res.0.3", cluster = "TC-3")

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

#stats.TC3 <- stats.function(drugs = top.diff, plots.list = plots.TC3, cluster = "TC-3")

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

# Create annotations for p-values

pvalue_col_fun = colorRamp2(c(0, 2, 3), c("white", "white", "yellow")) 
p.value.df <- stats.split.merged %>%
  select(all_of(starts_with("p.value")))

# p.values
pvalue.TC1 = p.value.df$`p.value.TC-1`
is_sig.TC1 = pvalue.TC1 < 0.01
pch.TC1 = rep("*", length(pvalue.TC1))
pch.TC1[!is_sig] = NA

list <- lapply(c(1:6), function(x){
  p.value <- p.value.df[,x]
  is_sig = p.value < 0.01
  pch = rep("*", length(pvalue.TC1))
  pch[!is_sig] = NA
  l <- list(p.value = p.value,
            pch = pch)
  return(l)
})
names(list) <- TCs
names(list)
TC.colors <- c("TC-1" = "#00b2d7",
               "TC-2" = "#e5c22f",
               "TC-3" = "#903ca2",
               "TC-4" = "#3f8741",
               "TC-5" = "#ff7b00",
               "TC-6" = "#cb5c42")

# Collapsed moas for breast-SSc signature
## Import manual collapsed moas drugs
collapsed.moas <- read_tsv(file = "./data/tsv/collapsed.moas.top.differential.drugs - top.differential.drugs.tsv")
collapsed.moas <- as.data.frame(collapsed.moas)
rownames(collapsed.moas) <- collapsed.moas$top.diff

ht <- Heatmap(
  matrix = matrix.cor[,1],
  name = "Correlation",
  width = unit(5, "cm"),
  row_labels = toupper(collapsed.moas$preferred.drug.names),
  row_names_gp = gpar(fontsize = 7),
  cluster_columns = F,
  show_column_names = F,
  cluster_rows = F,
  top_annotation = HeatmapAnnotation(TC = TCs[1],
                                     col = list("TC" = TC.colors)),
  right_annotation = rowAnnotation(
    pvalue.TC1 = anno_simple(-log10(list[["TC-1"]][["p.value"]]), 
                             col = pvalue_col_fun, 
                             pch = list[["TC-1"]][["pch"]])),
  layer_fun = function(j, i, x, y, width, height, fill) {
    # since grid.text can also be vectorized
    grid.text(sprintf("%.1f", pindex(matrix.cor, i, j)), x, y, 
              gp = gpar(fontsize = 10))
    }
)
ht
ht1 <- ht[, 1] + rowAnnotation(
  pvalue.TC1 = anno_simple(-log10(list[["TC-1"]][["p.value"]]), 
                           col = pvalue_col_fun, 
                           pch = list[["TC-1"]][["pch"]])
)
ht2 <- ht[, 2]+ rowAnnotation(
  pvalue.TC2 = anno_simple(-log10(list[["TC-2"]][["p.value"]]), 
                           col = pvalue_col_fun, 
                           pch = list[["TC-2"]][["pch"]])
  )
ht3 <- ht[, 3]+ rowAnnotation(
  pvalue.TC3 = anno_simple(-log10(list[["TC-3"]][["p.value"]]), 
                           col = pvalue_col_fun, 
                           pch = list[["TC-3"]][["pch"]])
)
ht4 <- ht[, 4]+ rowAnnotation(
  pvalue.TC4 = anno_simple(-log10(list[["TC-4"]][["p.value"]]), 
                           col = pvalue_col_fun, 
                           pch = list[["TC-4"]][["pch"]])
)
ht5 <- ht[, 5]+ rowAnnotation(
  pvalue.TC5 = anno_simple(-log10(list[["TC-5"]][["p.value"]]), 
                           col = pvalue_col_fun, 
                           pch = list[["TC-5"]][["pch"]])
)
ht6 <- ht[, 6]+ rowAnnotation(
  pvalue.TC6 = anno_simple(-log10(list[["TC-6"]][["p.value"]]), 
                           col = pvalue_col_fun, 
                           pch = list[["TC-6"]][["pch"]])
)


figure.correlation <- ht1 + ht2 + ht3 + ht4 + ht5 + ht6
figure.correlation

