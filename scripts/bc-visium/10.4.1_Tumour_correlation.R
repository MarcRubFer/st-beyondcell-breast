distrib.enrich <- bc.TCs.Tumour@normalized[top.diff.Tumour,] 
dim(distrib.enrich)
distrib.enrich <- t(distrib.enrich)
dim(distrib.enrich)
distrib.enrich <- as.data.frame(distrib.enrich) %>%
  rownames_to_column("spot")
distrib.metadata <- bc.TCs.Tumour@meta.data %>%
  select(Cancer.Epithelial, TCs_res.0.3) %>%
  rownames_to_column("spot")
distrib.total <- left_join(distrib.metadata,distrib.enrich, by = "spot")
rownames(distrib.total) <- distrib.total$spot
distrib.total$spot <- NULL
distrib.total <- distrib.total %>%
  pivot_longer(cols = -c(Cancer.Epithelial,TCs_res.0.3), names_to = "signatures", values_to = "enrichment")

corr.cancer <- distrib.total %>%
  group_by(TCs_res.0.3) %>%
  summarise(cor = cor(Cancer.Epithelial,enrichment))

ggplot(distrib.total, aes(x=Cancer.Epithelial, y=enrichment)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(vars(TCs_res.0.3))

library(corrplot)
library(ggstatsplot)

#To cite package 'ggstatsplot' in publications use:
#  
#  Patil, I. (2021). Visualizations with statistical details: The
#'ggstatsplot' approach. Journal of Open Source Software, 6(61), 3167,
#doi:10.21105/joss.03167

corr.drugs <- distrib.total %>%
  group_by(TCs_res.0.3,signatures) %>%
  summarise(cor = cor(Cancer.Epithelial,enrichment)) 

ggplot(corr.drugs, aes(x = signatures, y=cor)) +
  geom_bar(stat = "identity") +
  facet_grid(vars(TCs_res.0.3))

ggplot(corr.drugs, aes(x=cor,y=signatures)) +
  geom_bar(stat = "identity") +
  facet_grid(vars(TCs_res.0.3))


###############################################################################

drugs <- unique(as.factor(distrib.total$signatures))

df <- distrib.total %>%
  filter(TCs_res.0.3 == "TC-3")
df1 <- df %>%
  select(Cancer.Epithelial, "AM-580_PRISM_K06854232")
ggscatterstats(
  data = df1, ## data frame from which variables are taken
  x = Cancer.Epithelial, ## predictor/independent variable
  y = "AM-580_PRISM_K06854232", ## dependent variable
  marginal = F,
  bf.message = F,
  ggplot.component = list(ylim(-75,75),
                          xlim(0,1))
)
###############################################################################
plots.TC3 <- lapply(drugs, function(drug){
  df <- distrib.total %>%
    filter(TCs_res.0.3 == "TC-3") %>%
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
    
})
wrap_plots(plots.TC3, ncol = 5)
stats.TC3 <- lapply(seq_along(drugs), function(index){
  drug <- drugs[index]
  df <- as.data.frame(x = drug)
  stats <- extract_stats(plots.TC3[[index]])
  stats <- stats$subtitle_data
  df <- merge(df,stats)
  return(df)
})
stats.TC3.df <- do.call(rbind,stats.TC3)

stats.TC3.df.split <- stats.TC3.df %>%
  select(drug,estimate,p.value) %>%
  rename(estimate.TC3 = estimate,
         p.value.TC3 = p.value)
matrix.cor.TC3 <- stats.TC3.df.split %>%
  select(estimate.TC3)
matrix.cor.TC3 <- as.matrix(matrix.cor.TC3)

pvalue.TC3 = stats.TC3.df.split$p.value.TC3
is_sig = pvalue < 0.01
pch = rep("*", length(pvalue.TC3))
pch[!is_sig] = NA
pvalue_col_fun = colorRamp2(c(0, 2, 3), c("white", "white", "green")) 

ht.TC3 <- Heatmap(
  matrix = matrix.cor.TC3,
  name = "Corr Cancer-TC1",
  row_labels = stats.TC3.df.split$drug,
  cluster_columns = F,
  show_column_names = F,
  cluster_rows = F,
  #row_labels = toupper(stats.global$preferred.drug.names),
  top_annotation = HeatmapAnnotation(TC = c("TC-3"),
                                     col = list("TC" = TC.colors)),
  right_annotation = rowAnnotation(
    #pvalue = anno_simple(-log10(pvalue), col = pvalue_col_fun, pch = pch),
    numeric = anno_numeric(round(stats.TC3.df.split$estimate.TC3,2), 
                           align_to = "right",
                           bg_gp = gpar(fill = "white", col = "white"),
                           labels_gp = NULL,
                           labels_offset = unit(0, "cm"),
                           width = unit(1.25, "cm")), 
    pvalue = anno_simple(-log10(pvalue), 
                         col = pvalue_col_fun, 
                         pch = pch)),
  #col = colorRamp2(c(min, 0, max), c("blue", "white", "red")),
  #heatmap_legend_param = list(at = c(min, 0, max))
)
ht.TC3
###############################################################################
plots.TC4 <- lapply(drugs, function(drug){
  df <- distrib.total %>%
    filter(TCs_res.0.3 == "TC-4") %>%
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
  
})
wrap_plots(plots.TC4, ncol = 5)
stats.TC4 <- lapply(seq_along(drugs), function(index){
  drug <- drugs[index]
  df <- as.data.frame(x = drug)
  stats <- extract_stats(plots.TC4[[index]])
  stats <- stats$subtitle_data
  df <- merge(df,stats)
  return(df)
})
stats.TC4.df <- do.call(rbind,stats.TC4)

stats.TC4.df.split <- stats.TC4.df %>%
  select(drug,estimate,p.value) %>%
  rename(estimate.TC4 = estimate,
         p.value.TC4 = p.value)
matrix.cor.TC4 <- stats.TC4.df.split %>%
  select(estimate.TC4)
matrix.cor.TC4 <- as.matrix(matrix.cor.TC4)

pvalue.TC4 = stats.TC4.df.split$p.value.TC4
is_sig = pvalue < 0.01
pch = rep("*", length(pvalue.TC4))
pch[!is_sig] = NA
pvalue_col_fun = colorRamp2(c(0, 2, 3), c("white", "white", "green")) 

ht.TC4 <- Heatmap(
  matrix = matrix.cor.TC4,
  name = "Corr Cancer-TC1",
  row_labels = stats.TC4.df.split$drug,
  cluster_columns = F,
  show_column_names = F,
  cluster_rows = F,
  #row_labels = toupper(stats.global$preferred.drug.names),
  top_annotation = HeatmapAnnotation(TC = c("TC-4"),
                                     col = list("TC" = TC.colors)),
  right_annotation = rowAnnotation(
    #pvalue = anno_simple(-log10(pvalue), col = pvalue_col_fun, pch = pch),
    numeric = anno_numeric(round(stats.TC4.df.split$estimate.TC4,2), 
                           align_to = "right",
                           bg_gp = gpar(fill = "white", col = "white"),
                           labels_gp = NULL,
                           labels_offset = unit(0, "cm"),
                           width = unit(1.25, "cm")), 
    pvalue = anno_simple(-log10(pvalue), 
                         col = pvalue_col_fun, 
                         pch = pch)),
  #col = colorRamp2(c(min, 0, max), c("blue", "white", "red")),
  #heatmap_legend_param = list(at = c(min, 0, max))
)
ht.TC4
###############################################################################
plots.TC5 <- lapply(drugs, function(drug){
  df <- distrib.total %>%
    filter(TCs_res.0.3 == "TC-5") %>%
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
  
})
wrap_plots(plots.TC5, ncol = 5)
stats.TC5 <- lapply(seq_along(drugs), function(index){
  drug <- drugs[index]
  df <- as.data.frame(x = drug)
  stats <- extract_stats(plots.TC5[[index]])
  stats <- stats$subtitle_data
  df <- merge(df,stats)
  return(df)
})
stats.TC5.df <- do.call(rbind,stats.TC5)

stats.TC5.df.split <- stats.TC5.df %>%
  select(drug,estimate,p.value) %>%
  rename(estimate.TC5 = estimate,
         p.value.TC5 = p.value)
matrix.cor.TC5 <- stats.TC5.df.split %>%
  select(estimate.TC5)
matrix.cor.TC5 <- as.matrix(matrix.cor.TC5)

pvalue.TC5 = stats.TC5.df.split$p.value.TC5
is_sig = pvalue < 0.01
pch = rep("*", length(pvalue.TC5))
pch[!is_sig] = NA
pvalue_col_fun = colorRamp2(c(0, 2, 3), c("white", "white", "green")) 

ht.TC5 <- Heatmap(
  matrix = matrix.cor.TC5,
  name = "Corr Cancer-TC1",
  row_labels = stats.TC5.df.split$drug,
  cluster_columns = F,
  show_column_names = F,
  cluster_rows = F,
  #row_labels = toupper(stats.global$preferred.drug.names),
  top_annotation = HeatmapAnnotation(TC = c("TC-5"),
                                     col = list("TC" = TC.colors)),
  right_annotation = rowAnnotation(
    #pvalue = anno_simple(-log10(pvalue), col = pvalue_col_fun, pch = pch),
    numeric = anno_numeric(round(stats.TC5.df.split$estimate.TC5,2), 
                           align_to = "right",
                           bg_gp = gpar(fill = "white", col = "white"),
                           labels_gp = NULL,
                           labels_offset = unit(0, "cm"),
                           width = unit(1.25, "cm")), 
    pvalue = anno_simple(-log10(pvalue), 
                         col = pvalue_col_fun, 
                         pch = pch)),
  #col = colorRamp2(c(min, 0, max), c("blue", "white", "red")),
  #heatmap_legend_param = list(at = c(min, 0, max))
)
ht.TC5
###############################################################################
plots.TC6 <- lapply(drugs, function(drug){
  df <- distrib.total %>%
    filter(TCs_res.0.3 == "TC-6") %>%
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
  
})
wrap_plots(plots.TC6, ncol = 5)
stats.TC6 <- lapply(seq_along(drugs), function(index){
  drug <- drugs[index]
  df <- as.data.frame(x = drug)
  stats <- extract_stats(plots.TC6[[index]])
  stats <- stats$subtitle_data
  df <- merge(df,stats)
  return(df)
})
stats.TC6.df <- do.call(rbind,stats.TC6)

stats.TC6.df.split <- stats.TC6.df %>%
  select(drug,estimate,p.value) %>%
  rename(estimate.TC6 = estimate,
         p.value.TC6 = p.value)
matrix.cor.TC6 <- stats.TC6.df.split %>%
  select(estimate.TC6)
matrix.cor.TC6 <- as.matrix(matrix.cor.TC6)

pvalue.TC6 = stats.TC6.df.split$p.value.TC6
is_sig = pvalue < 0.01
pch = rep("*", length(pvalue.TC6))
pch[!is_sig] = NA
pvalue_col_fun = colorRamp2(c(0, 2, 3), c("white", "white", "green")) 

ht.TC6 <- Heatmap(
  matrix = matrix.cor.TC6,
  name = "Corr Cancer-TC1",
  row_labels = stats.TC6.df.split$drug,
  cluster_columns = F,
  show_column_names = F,
  cluster_rows = F,
  #row_labels = toupper(stats.global$preferred.drug.names),
  top_annotation = HeatmapAnnotation(TC = c("TC-6"),
                                     col = list("TC" = TC.colors)),
  right_annotation = rowAnnotation(
    #pvalue = anno_simple(-log10(pvalue), col = pvalue_col_fun, pch = pch),
    numeric = anno_numeric(round(stats.TC6.df.split$estimate.TC6,2), 
                           align_to = "right",
                           bg_gp = gpar(fill = "white", col = "white"),
                           labels_gp = NULL,
                           labels_offset = unit(0, "cm"),
                           width = unit(1.25, "cm")), 
    pvalue = anno_simple(-log10(pvalue), 
                         col = pvalue_col_fun, 
                         pch = pch)),
  #col = colorRamp2(c(min, 0, max), c("blue", "white", "red")),
  #heatmap_legend_param = list(at = c(min, 0, max))
)
ht.TC6
###############################################################################

ht <- ht.TC3 + ht.TC4 + ht.TC5 + ht.TC6

valnemulin <- distrib.total %>%
  select(2:3,"valnemulin_PRISM_K33813875") %>%
  rename(enrichment = "valnemulin_PRISM_K33813875") %>%
  grouped_ggscatterstats(data = .,
                         x = Cancer.Epithelial,
                         y = enrichment,
                         grouping.var = TCs_res.0.3,
                         plotgrid.args = list(nrow = 4,
                                              ncol = 1),
                         marginal = F,
                         bf.message = F,
                         ggplot.component = list(ylim(-50,50),
                                                 xlim(0.3,1)
                         ),
                         annotation.args = list(
                           title = "VALNEMULIN",
                           subtitle = "Correlation Cancer Epithelial proportion vs BCS enrichment"
                         ))
valnemulin

simvastatin <- distrib.total %>%
  select(2:3,"simvastatin_CTRP_K22134346") %>%
  rename(enrichment = "simvastatin_CTRP_K22134346") %>%
  grouped_ggscatterstats(data = .,
                         x = Cancer.Epithelial,
                         y = enrichment,
                         grouping.var = TCs_res.0.3,
                         plotgrid.args = list(nrow = 4,
                                              ncol = 1),
                         marginal = F,
                         bf.message = F,
                         ggplot.component = list(ylim(-50,50),
                                                 xlim(0.3,1)
                         ),
                         annotation.args = list(
                           title = "SINVASTATIN",
                           subtitle = "Correlation Cancer Epithelial proportion vs BCS enrichment"
                         ))
simvastatin
