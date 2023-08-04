head(bc.recomputed@meta.data)

tc.composition <- bc.recomputed@meta.data %>%
  select(spot.composition.filter, bc_clusters_res.0.3) 

tc.composition2 <- tc.composition %>%
  count(bc_clusters_res.0.3, spot.composition.filter) %>%
  group_by(bc_clusters_res.0.3) %>%
  mutate(percentage = n / sum(n) * 100)
head(tc.composition2)

ggplot(data = tc.composition2, mapping = aes(x=bc_clusters_res.0.3, y = percentage)) +
  geom_col(aes(fill = spot.composition.filter))


new_df <- data.frame()
for (i in 1:nrow(tc.composition)) {
  spot.composition <- tc.composition$spot.composition.filter[i]
  
  if (grepl(pattern = "\\+", spot.composition)) {
    terms <- strsplit(spot.composition, "\\+")[[1]]
    
    for (term in terms) {
      new_row <- tc.composition[i, ]
      new_row$cell.composition <- term
      new_df <- rbind(new_df, new_row)
    }
  } else {
    tc.composition[i, "cell.composition"] <- spot.composition
    new_df <- rbind(new_df, tc.composition[i, ])
  }
}
new_df
rownames(new_df) <- NULL
head(new_df)

tc.composition3 <- new_df %>%
  mutate(spot.composition.filter = NULL) %>%
  count(bc_clusters_res.0.3, cell.composition) %>%
  group_by(bc_clusters_res.0.3) %>%
  mutate(percentage = n / sum(n) * 100)
head(tc.composition3)

ggplot(data = tc.composition3, mapping = aes(x=bc_clusters_res.0.3, y = percentage)) +
  geom_col(aes(fill = cell.composition), col = "black", linewidth = 0.2) +
  ggtitle("Cell type composition within each Beyondcell cluster")

VlnPlot(seurat.spatial, features = "percent.mt", group.by = "bc_clusters_res.0.3")
seurat.spatial@meta.data
