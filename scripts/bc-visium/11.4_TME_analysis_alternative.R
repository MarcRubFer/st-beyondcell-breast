TC.colors <- c("TC-1" = "#00b2d7",
               "TC-2" = "#e5c22f",
               "TC-3" = "#903ca2",
               "TC-4" = "#3f8741",
               "TC-5" = "#ff7b00",
               "TC-6" = "#cb5c42")

# Analysis of TME TCs (TCs 1 and 2)
TCs.TME <- levels(bc.recomputed@meta.data$TCs_res.0.3)[1:2]
cells.TCs.TME <- bc.recomputed@meta.data %>%
  filter(TCs_res.0.3 %in% TCs.TME) %>%
  rownames_to_column(var = "spots") %>%
  pull(spots)

bc.TCs.TME <- bcSubset(bc.recomputed, cells = cells.TCs.TME)

head(bc.TCs.TME@meta.data)

df <- bc.TCs.TME@meta.data %>%
  select(B.cells:T.cells, spot.collapse, TCs_res.0.3) %>%
  pivot_longer(cols = B.cells:T.cells, names_to = "decon.cell.type", values_to = "proportion")

head(df)

ggplot(df, aes(x=spot.collapse, y=proportion)) +
  geom_boxplot(aes(fill = decon.cell.type)) +
  facet_wrap(vars(TCs_res.0.3))

df2 <- df %>%
  filter(TCs_res.0.3 == "TC-2")

ggplot(df2, aes(x=spot.collapse, y=proportion)) +
  geom_boxplot(aes(fill = decon.cell.type)) +
  facet_wrap(vars(TCs_res.0.3))

df3 <- df %>%
  filter(decon.cell.type == "Cancer.Epithelial")

ggplot(df3, aes(x=TCs_res.0.3, y=proportion, fill = TCs_res.0.3)) +
  geom_boxplot() +
  ggtitle(label = "Proportion of Cancer Epithelial by TC", subtitle = "TC-1 and TC-2")

distrib.enrich <- bc.TCs.TME@normalized 
dim(distrib.enrich)
distrib.enrich <- t(distrib.enrich)
distrib.enrich <- as.data.frame(distrib.enrich) %>%
  rownames_to_column("spot")
distrib.metadata <- bc.TCs.TME@meta.data %>%
  select(TCs_res.0.3) %>%
  rownames_to_column("spot")

distrib.total <- left_join(distrib.metadata,distrib.enrich, by = "spot")
head(distrib.total)
distrib.total.TC2 <- distrib.total %>%
  filter(TCs_res.0.3 == "TC-2") %>%
  mutate(spot = NULL) %>%
  pivot_longer(cols = -TCs_res.0.3, names_to = "signatures", values_to = "enrichment")

ggplot(distrib.total.TC2, aes(x=enrichment)) +
  geom_histogram()


distrib.enrich.pathway <- bc.TME.reactome@normalized 
distrib.enrich.pathway <- t(distrib.enrich.pathway)
distrib.enrich.pathway <- as.data.frame(distrib.enrich.pathway) %>%
  rownames_to_column("spot")
distrib.total.pathway <- left_join(distrib.metadata,distrib.enrich.pathway, by = "spot")
distrib.total.TC2.pathway <- distrib.total.pathway %>%
  filter(TCs_res.0.3 == "TC-2") %>%
  mutate(spot = NULL) %>%
  pivot_longer(cols = -TCs_res.0.3, names_to = "signatures", values_to = "enrichment")

ggplot(distrib.total.TC2.pathway, aes(x=enrichment)) +
  geom_histogram()


distrib.enrich <- bc.TCs.TME@normalized[top.diff.TME,] 
dim(distrib.enrich)
distrib.enrich <- t(distrib.enrich)
dim(distrib.enrich)
distrib.enrich <- as.data.frame(distrib.enrich) %>%
  rownames_to_column("spot")
distrib.metadata <- bc.TCs.TME@meta.data %>%
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

corr.drugs.TC1 <- corr.drugs %>%
  filter(TCs_res.0.3 == "TC-1")

ggplot(corr.drugs.TC1, aes(x = signatures, y = 1, size = abs(cor), fill = cor)) +
  geom_point(shape = 21) +  # Burbujas
  scale_size(range = c(5, 15)) +  # Rango de tamaños de burbujas
  geom_text(aes(label = round(cor, 2)), size = 3, color = "black", vjust = -5) +  # Etiquetas de correlación dentro de las burbujas
  labs(x = "Variables", y = "Cancer Epithelial", title = "Correlation Drugs-Cancer Epithelial") +
  theme_minimal() +
  guides(fill = guide_legend(title = "Correlation")) +  # Leyenda de colores
  scale_fill_gradient(low = "blue", high = "red") + # Colores de la escala
  theme(panel.grid =  element_blank(),
        axis.text.x = element_text(angle = 90),axis.title.y = element_blank(),  # Eliminar título del eje Y
        axis.text.y = element_blank(),  # Eliminar etiquetas del eje Y
        axis.line.y = element_blank(),  # Eliminar línea del eje Y
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 0, r = 0.5, b = 0.5, l = 0.5, unit = "cm"))


##############################################################################

ggplot(distrib.total, aes(x=Cancer.Epithelial, y=enrichment)) +
  geom_point(aes(colour=TCs_res.0.3)) +
  geom_point(colour = "grey90", size = 0.2) +
  geom_smooth(aes(colour=TCs_res.0.3),method = "lm") +
  facet_grid(vars(TCs_res.0.3)) +
  facet_wrap(vars(signatures))

importazole <- distrib.total %>%
  filter(signatures == "importazole_CTRP_A02481876") %>%
  grouped_ggscatterstats(data = .,
                       x = Cancer.Epithelial,
                       y = enrichment,
                       grouping.var = TCs_res.0.3,
                       plotgrid.args = list(nrow = 2,
                                            ncol = 1),
                       marginal = F,
                       bf.message = F,
                       ggplot.component = list(ylim(-50,50),
                                               xlim(0,1)
                                               ),
                       annotation.args = list(
                         title = "IMPORTAZOLE",
                         subtitle = "Correlation Cancer Epithelial proportion vs BCS enrichment"
                       ))
stats <- extract_stats(importazole)
stats$subtitle_data

staurosporine <- distrib.total %>%
  filter(signatures == "Staurosporine_GDSC_1034") %>%
  grouped_ggscatterstats(data = .,
                         x = Cancer.Epithelial,
                         y = enrichment,
                         grouping.var = TCs_res.0.3,
                         plotgrid.args = list(nrow = 2,
                                              ncol = 1),
                         marginal = F,
                         bf.message = F,
                         ggplot.component = list(ylim(-75,75),
                                                 xlim(0,1)),
                         annotation.args = list(
                           title = "STAUROSPORINE",
                           subtitle = "Correlation Cancer Epithelial proportion vs BCS enrichment"
                         ))
importazole | staurosporine

TC1.drugs <- distrib.total %>%
  filter(TCs_res.0.3 == "TC-1")

drugs <- levels(as.factor(TC1.drugs$signatures))

plots.TC1 <- lapply(drugs, function(drug){
  df <- TC1.drugs %>%
    filter(signatures == drug) %>%
    grouped_ggscatterstats(data = .,
                           x = Cancer.Epithelial,
                           y = enrichment,
                           grouping.var = TCs_res.0.3,
                           plotgrid.args = list(nrow = 2,
                                                ncol = 1),
                           marginal = F,
                           bf.message = F,
                           ggplot.component = list(ylim(-75,75),
                                                   xlim(0,1),
                                                   ggtitle(drug)),
                           annotation.args = list(
                             title = drug,
                             subtitle = "Correlation Cancer Epithelial proportion vs BCS enrichment"
                           ))
})
wrap_plots(plots.TC1, ncol = 5)
extract_stats(plots.TC1[[1]])


stats.TC1 <- lapply(seq_along(drugs), function(index){
  drug <- drugs[index]
  df <- as.data.frame(x = drug)
  stats <- extract_stats(plots.TC1[[index]])
  stats <- stats$subtitle_data
  df <- merge(df,stats)
  return(df)
})

stats.TC1.df <- do.call(rbind,stats.TC1)


TC2.drugs <- distrib.total %>%
  filter(TCs_res.0.3 == "TC-2")

drugs.TC2 <- levels(as.factor(TC2.drugs$signatures))

plots.TC2 <- lapply(drugs.TC2, function(drug){
  df <- TC2.drugs %>%
    filter(signatures == drug) %>%
    grouped_ggscatterstats(data = .,
                           x = Cancer.Epithelial,
                           y = enrichment,
                           grouping.var = TCs_res.0.3,
                           plotgrid.args = list(nrow = 2,
                                                ncol = 1),
                           marginal = F,
                           bf.message = F,
                           ggplot.component = list(ylim(-75,75),
                                                   xlim(0,1),
                                                   ggtitle(drug)),
                           annotation.args = list(
                             title = drug,
                             subtitle = "Correlation Cancer Epithelial proportion vs BCS enrichment"
                           ))
})
wrap_plots(plots.TC2, ncol = 5)

stats.TC2 <- lapply(seq_along(drugs), function(index){
  drug <- drugs[index]
  df <- as.data.frame(x = drug)
  stats <- extract_stats(plots.TC2[[index]])
  stats <- stats$subtitle_data
  df <- merge(df,stats)
  return(df)
})

stats.TC2.df <- do.call(rbind,stats.TC2)

stats.TC1.df.split <- stats.TC1.df %>%
  select(drug,estimate,p.value) %>%
  rename(estimate.TC1 = estimate,
         p.value.TC1 = p.value)

stats.TC2.df.split <- stats.TC2.df %>%
  select(drug,estimate,p.value) %>%
  rename(estimate.TC2 = estimate,
         p.value.TC2 = p.value)

stats.global <- left_join(stats.TC1.df.split, stats.TC2.df.split, by = "drug")

matrix.cor <- stats.global %>%
  select(estimate.TC1,estimate.TC2)
matrix.cor <- as.matrix(matrix.cor)

pvalue = stats.global$p.value.TC1
is_sig = pvalue < 0.01
pch = rep("*", length(pvalue))
pch[!is_sig] = NA

pvalue.TC2 = stats.global$p.value.TC2
is_sig2 = pvalue.TC2 < 0.01
pch2 = rep("*", length(pvalue.TC2))
pch2[!is_sig2] = NA
# color mapping for -log10(pvalue)
pvalue_col_fun = colorRamp2(c(0, 2, 3), c("white", "white", "green")) 

min <- round(min(apply(matrix.cor, 1, function(row) min(row))),2)
max <- round(max(apply(matrix.cor, 1, function(row) max(row))),2)

row.namex <-collapsed.moas.TME %>%
  select(preferred.drug.names) %>%
  rownames_to_column("drug")

stats.global <- left_join(stats.global,row.namex, by = "drug")

ht.TC1 <- Heatmap(
  matrix = matrix.cor[,"estimate.TC1"],
  name = "Corr Cancer-TC1",
  row_labels = stats.global$drug,
  cluster_columns = F,
  show_column_names = F,
  cluster_rows = F,
  #row_labels = toupper(stats.global$preferred.drug.names),
  top_annotation = HeatmapAnnotation(TC = c("TC-1"),
                                     col = list("TC" = TC.colors)),
  right_annotation = rowAnnotation(
    #pvalue = anno_simple(-log10(pvalue), col = pvalue_col_fun, pch = pch),
    numeric = anno_numeric(round(stats.global$estimate.TC1,2), 
                           align_to = "right",
                           bg_gp = gpar(fill = "white", col = "white"),
                           labels_gp = NULL,
                           labels_offset = unit(0, "cm"),
                           width = unit(1.25, "cm")), 
    pvalue = anno_simple(-log10(pvalue), 
                         col = pvalue_col_fun, 
                         pch = pch)),
  col = colorRamp2(c(min, 0, max), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(min, 0, max))
)
ht.TC1
ht.TC2 <- Heatmap(
  matrix = matrix.cor[,"estimate.TC2"],
  name = "Corr Cancer-TC2",
  row_labels = stats.global$drug,
  cluster_columns = F,
  show_column_names = F,
  cluster_rows = F,
  top_annotation = HeatmapAnnotation(TC = c("TC-2"),
                                     col = list("TC" = TC.colors)),
  right_annotation = rowAnnotation(
    #pvalue = anno_simple(-log10(pvalue), col = pvalue_col_fun, pch = pch),
    numeric = anno_numeric(round(stats.global$estimate.TC2,2), 
                           align_to = "right",
                           bg_gp = gpar(fill = "white", col = "white"),
                           labels_gp = NULL,
                           labels_offset = unit(0, "cm"),
                           width = unit(1.25, "cm")), 
    pvalue = anno_simple(-log10(pvalue.TC2), 
                         col = pvalue_col_fun, 
                         pch = pch2)),
  col = colorRamp2(c(min, 0, max), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(min, 0, max))
)
ht.TC2


ht.TC1 + ht.TC2




drugs.total <- levels(as.factor(distrib.total$signatures))
plots.total <- lapply(drugs.total, function(drug){
  df <- distrib.total %>%
    filter(signatures == drug) %>%
    grouped_ggscatterstats(data = .,
                           x = Cancer.Epithelial,
                           y = enrichment,
                           grouping.var = TCs_res.0.3,
                           plotgrid.args = list(nrow = 2,
                                                ncol = 1),
                           marginal = F,
                           bf.message = F,
                           ggplot.component = list(ylim(-75,75),
                                                   xlim(0,1),
                                                   ggtitle(drug)),
                           annotation.args = list(
                             title = drug,
                             subtitle = "Correlation Cancer Epithelial proportion vs BCS enrichment"
                           ))
})
wrap_plots(plots.total, ncol = 5)
extract_stats(plots.total[[1]])
