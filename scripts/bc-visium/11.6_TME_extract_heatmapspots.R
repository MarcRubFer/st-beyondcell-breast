# Draw Heatmap
heatmap.drugs.TME.cancerepith <- Heatmap(
  drugs.matrix.TME2,
  name = "bcScore",
  cluster_columns = T,
  top_annotation = HeatmapAnnotation("TCs" = col.order.TME2$TCs_res.0.3,
                                     "Cancer Epith" = col.order.TME2$Cancer.Epithelial,
                                     "Lymphoid" = col.order.TME2$Lymphoid,
                                     "B cells" = col.order.TME2$B.cells,
                                     "T cells" = col.order.TME2$T.cells,
                                     "Cell type" = col.order.TME2$spot.collapse,
                                     col = list("TCs" = TC.colors[1:2],
                                                "Cell type" = colors.categories,
                                                "Cancer Epith" = scale.cancer,
                                                "Lymphoid" = scale.lymphoid,
                                                "B cells" = scale.bcells,
                                                "T cells" = scale.tcells)),
  right_annotation = rowAnnotation(TC1.sens = collapsed.moas.TME$TC1.sensitivity,
                                   TC2.sens = collapsed.moas.TME$TC2.sensitivity,
                                   MoA = collapsed.moas.TME$collapsed.MoAs,
                                   col = list(MoA = cols.drugs.TME,
                                              TC1.sens = col.sensitivity,
                                              TC2.sens = col.sensitivity)),
  show_column_names = FALSE,
  #column_split = col.order.TME2$TCs_res.0.3,
  column_km = 6,
  row_names_gp = gpar(fontsize = 6),
  row_labels = toupper(collapsed.moas.TME$preferred.drug.names),
  cluster_rows = T,
  row_split =  5,
  col = colorRamp2(c(drugs.matrix.TME.min, 0, drugs.matrix.TME.max), c("blue", "white", "red")),
  heatmap_legend_param = list(at = c(drugs.matrix.TME.min, 0, drugs.matrix.TME.max))
)      
heatmap.drugs.TME.cancerepith
c.order.list <- column_order(heatmap.drugs.TME.cancerepith)
c.dend <- column_dend(heatmap.drugs.TME.cancerepith)

lapply(c.order.list, function(x) length(x))
