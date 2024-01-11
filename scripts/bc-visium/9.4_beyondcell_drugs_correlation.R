rm(list = ls())

library("beyondcell")
library("Seurat")
library("tidyverse")
library("tidygraph")
library("patchwork")
library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")

out.dir <- "./results"
dir.create(path = out.dir, recursive = TRUE)

# establish seed
set.seed(1)

# Load beyondcell object
bc.ranked <- readRDS("./results/analysis/bc_ranked_all.rds")

# Extract top High/Low-Differential Drugs 
top.diff <- as.data.frame(bc.ranked@ranks) %>%
  select(starts_with(match = "TCs_res.0.3.group.")) %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = starts_with("TCs_res.0.3.group."), names_to = "cluster", values_to = "group") %>%
  filter(group != is.na(group),
         grepl("Differential", group)) %>%
  pull("signature") %>%
  unique()

# Extract data about of enrichment (bcScore) of top.diff
distrib.enrich <- bc.ranked@normalized[top.diff,] 
dim(distrib.enrich)
distrib.enrich <- t(distrib.enrich)
dim(distrib.enrich)
distrib.enrich <- as.data.frame(distrib.enrich) %>%
  rownames_to_column("spot")

#Extract metadata for Cancer epithelial proportion and TCs
distrib.metadata <- bc.ranked@meta.data %>%
  select(Cancer.Epithelial, TCs_res.0.3) %>%
  rownames_to_column("spot")

# Merge by spot-ID both data
distrib.total <- right_join(distrib.metadata,distrib.enrich, by = "spot")
rownames(distrib.total) <- distrib.total$spot

#distrib.total.pivot <- distrib.total %>%
#  pivot_longer(cols = -c(spot,Cancer.Epithelial,TCs_res.0.3), names_to = "signatures", values_to = "enrichment")

library(ggstatsplot)

#To cite package 'ggstatsplot' in publications use:
#  
#  Patil, I. (2021). Visualizations with statistical details: The
#'ggstatsplot' approach. Journal of Open Source Software, 6(61), 3167,
#doi:10.21105/joss.03167

# This function create correlation plots (ggstatplot library). Correlation between
# top.diff drug and its enrichment score
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

# Extract levels of TCs
TCs <- unique(levels(bc.ranked@meta.data$TCs_res.0.3))

# Create a list of plots with the information of correlation drug-enrichment in 
# each TC
plots.TCs <- lapply(TCs, function(TC){
  plot.TC <- plot.correlation.list(df = distrib.total,
                                   drugs = top.diff,
                                   idents = "TCs_res.0.3",
                                   cluster = TC)
})

# This function extract the information of each correlation drug-enrichment plot
# and create a data frame
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

# Create a list with the extracted information from correlation plots
stats.TCs <- lapply(seq_along(TCs), function(TC){
  plots <- plots.TCs[[TC]]
  stats <- stats.function(drugs = top.diff,
                          plots.list = plots,
                          cluster = TCs[TC])
  return(stats)
})

# Merge the information of all TCs in a unique data frame
stats.TCs.global <- do.call(rbind,stats.TCs)
stats.TCs.global <- stats.TCs.global %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))


write.table(stats.TCs.global,
            file = "./results/tables/table_correlation_drugs_enrichment_global.tsv",
            sep = "\t")



stats.split.list <- lapply(seq_along(TCs), function(index){
  split <- stats.TCs.global %>%
    filter(TC == TCs[index]) %>%
    select(drug, estimate, p.value, p.adj) %>%
    rename_with(.fn = ~paste0("estimate.",TCs[index]), .cols = estimate) %>%
    rename_with(.fn = ~paste0("p.value.",TCs[index]), .cols = p.value) %>%
    rename_with(.fn = ~paste0("p.adj.",TCs[index]), .cols = p.adj) %>%
    mutate(drug = NULL)
})
stats.split.merged <- do.call(cbind,stats.split.list)
rownames(stats.split.merged) <- top.diff

write.table(stats.split.merged,
            file = "./results/tables/table_correlation_drugs_enrichment_global.tsv",
            sep = "\t")

matrix.cor <- stats.split.merged %>%
  select(starts_with("estimate")) %>%
  rename_with(.fn = ~gsub(pattern = "estimate.",replacement = "",.x)) %>%
  as.matrix()

# Create annotations for p-values

pvalue_col_fun = colorRamp2(c(0, -log10(0.05), -log10(0.01)), c("white", "white", "white")) 
#pvalue_col_fun = colorRamp2(breaks = c(0, -log10(0.05), -log10(0.01)), 
#                            hcl_palette = "YlOrBr")
#                            ,
#                            c("white", "yellow", "orange")) 
pvalue_col_fun = colorRamp2(c(0, -log10(0.05), -log10(0.01)), c("green", "white", "red"))
pvalue_col_fun = colorRamp2(c(0, -log10(0.05), -log10(0.01), -log10(0.001)), hcl_palette = "Purp", reverse = F)

p.adj.df <- stats.split.merged %>%
  select(all_of(starts_with("p.adj")))

# Create a list with contain significance of p.values and named
list.pvalue <- lapply(c(1:6), function(x){
  p.adj <- p.adj.df[,x]
  is_sig = p.adj < 0.05
  pch = rep("*", length(p.adj))
  pch[!is_sig] = NA
  l <- list(p.adj = p.adj,
            pch = pch)
  return(l)
})
names(list.pvalue) <- TCs
names(list.pvalue)

# Colors for TCs legend
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


# Create Heatmap
ht <- Heatmap(
  matrix = matrix.cor,
  name = "Correlation",
  width = unit(20, "cm"),
  row_labels = toupper(collapsed.moas$preferred.drug.names),
  row_names_gp = gpar(fontsize = 7),
  cluster_columns = F,
  show_column_names = F,
  cluster_rows = F,
  top_annotation = HeatmapAnnotation(TC = TCs,
                                     col = list("TC" = TC.colors)),
  right_annotation = rowAnnotation(
    pval.adj.TC1 = anno_simple(-log10(list.pvalue[["TC-1"]][["p.adj"]]), 
                           col = pvalue_col_fun, 
                           pch = list.pvalue[["TC-1"]][["pch"]]),
    pval.adj.TC2 = anno_simple(-log10(list.pvalue[["TC-2"]][["p.adj"]]), 
                           col = pvalue_col_fun, 
                           pch = list.pvalue[["TC-2"]][["pch"]]),
    pval.adj.TC3 = anno_simple(-log10(list.pvalue[["TC-3"]][["p.adj"]]), 
                           col = pvalue_col_fun, 
                           pch = list.pvalue[["TC-3"]][["pch"]]),
    pval.adj.TC4 = anno_simple(-log10(list.pvalue[["TC-4"]][["p.adj"]]), 
                           col = pvalue_col_fun, 
                           pch = list.pvalue[["TC-4"]][["pch"]]),
    pval.adj.TC5 = anno_simple(-log10(list.pvalue[["TC-5"]][["p.adj"]]), 
                           col = pvalue_col_fun, 
                           pch = list.pvalue[["TC-5"]][["pch"]]),
    pval.adj.TC6 = anno_simple(-log10(list.pvalue[["TC-6"]][["p.adj"]]), 
                             col = pvalue_col_fun, 
                             pch = list.pvalue[["TC-6"]][["pch"]])
    ),
  layer_fun = function(j, i, x, y, width, height, fill) {
    # since grid.text can also be vectorized
    grid.text(sprintf("%.1f", pindex(matrix.cor, i, j)), x, y, 
              gp = gpar(fontsize = 10))
    }
)
ht <- draw(ht, merge_legend = TRUE)
ht
# Save heatmap as png
png(filename = "./results/plots/Beyondcell_oct23_DrugRank/Correlation_heatmap.png",
    width = 48,
    height = 24,
    units = "cm",
    res = 320)
draw(ht)
dev.off()

# now we generate two legends, one for the p-value
# see how we define the legend for pvalue
lgd_pvalue = Legend(title = "p-value.adj", col_fun = pvalue_col_fun, at = c(0,-log10(0.1), -log10(0.01), -log10(0.001)), 
                    labels = c("1","0.1", "0.01", "0.001"))
# and one for the significant p-values
lgd_sig = Legend(title = "p-value.adj", pch =  "*", type = "points", labels = "< 0.05")
# these two self-defined legends are added to the plot by `annotation_legend_list`
ht2 <- draw(ht, annotation_legend_list = list(lgd_sig))
ht2

# Save heatmap as png
png(filename = "./results/plots/Beyondcell_oct23_DrugRank/Correlation_heatmap_new.png",
    width = 48,
    height = 24,
    units = "cm",
    res = 320)
draw(ht2)
dev.off()
