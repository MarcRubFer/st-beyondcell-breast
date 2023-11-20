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

corr.drugs <- distrib.total %>%
  group_by(TCs_res.0.3,signatures) %>%
  summarise(cor = cor(Cancer.Epithelial,enrichment)) 

ggplot(corr.drugs, aes(x = signatures, y=cor)) +
  geom_bar(stat = "identity") +
  facet_grid(vars(TCs_res.0.3))

ggplot(corr.drugs, aes(x=cor,y=signatures)) +
  geom_bar(stat = "identity") +
  facet_grid(vars(TCs_res.0.3))

ggplot(distrib.total, aes(x=Cancer.Epithelial, y=enrichment, colour = TCs_res.0.3)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(vars(signatures))

importazole <- distrib.total %>%
  filter(signatures == "importazole_CTRP_A02481876") 
importazole
grouped_ggscatterstats(data = importazole,
                       x = Cancer.Epithelial,
                       y = enrichment,
                       grouping.var = TCs_res.0.3,
                       plotgrid.args = list(nrow = 2,
                                            ncol = 1),
                       marginal = F,
                       bf.message = F,
                       ggplot.component = list(ylim(-50,50),
                                               xlim(0,1)),
                       annotation.args = list(
                         title = "IMPORTAZOLE",
                         subtitle = "Correlation Cancer Epithelial proportion vs BCS enrichment"
                       ))


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
staurosporine

TC1.drugs <- distrib.total %>%
  filter(TCs_res.0.3 == "TC-1")

drugs <- levels(as.factor(TC1.drugs$signatures))

plots <- lapply(drugs, function(drug){
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
wrap_plots(plots, ncol = 5)
