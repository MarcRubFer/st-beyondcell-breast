# Signatures in bc.
sigs <- rownames(bc.recomputed@normalized)
# Cells in bc.
cells <- colnames(bc.recomputed@normalized)

idents <- "TCs_res.0.3"
resm.cutoff = c(0.05,0.95)
sp.cutoff = c(0.1, 0.4, 0.6, 0.9)
# Keep column to group by.
meta <- bc.recomputed@meta.data %>%
  tibble::rownames_to_column("cells") %>%
  dplyr::select(cells, all_of(idents)) %>%
  dplyr::rename(group.var := !!idents) %>%
  dplyr::mutate(group.var = factor(group.var)) %>%
  unique()
lvls <- levels(meta$group.var)
# Column to order by.
order.col <- paste0("rank.", levels(meta$group.var)[1])
cols.additional <- c("median", "sd", "variance", "min", "max", "prop.na")
cols.stats <- c("rank", "rank.median", "switch.point", "mean", cols.additional, 
                "residuals.mean", "residuals.median", "group.mean", "group.median")

cols.stats.level <- tidyr::expand_grid(lvls, cols.stats) %>%
  dplyr::mutate(col.name = paste(cols.stats, lvls, sep = ".")) %>%
  dplyr::pull(col.name)

sp <- data.frame(switch.point = bc.recomputed@switch.point) %>%
  tibble::rownames_to_column("IDs")

# Compute long normalized BCS.
normalized.long <- bc.recomputed@normalized %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cells") %>%
  tidyr::pivot_longer(cols = all_of(sigs), names_to = "IDs", 
                      values_to = "enrichment", values_drop_na = FALSE)

# Add grouping information and switch point.
normalized.long <- normalized.long %>%
  dplyr::inner_join(sp, by = "IDs") %>%
  dplyr::inner_join(meta, by = "cells")

# Compute mean BCS and residual's mean per signature.
stats.long <- normalized.long %>%
  dplyr::group_by(IDs) %>%
  dplyr::mutate(mean = round(mean(enrichment, na.omit = TRUE), digits = 2)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(resid = enrichment - mean) %>%
  dplyr::group_by(IDs, group.var) %>%
  dplyr::mutate(residuals.mean = round(mean(resid, na.rm = TRUE), digits = 2),
                residuals.median = round(median(resid, na.rm = TRUE), digits = 2)) %>%
  dplyr::ungroup()

stats.long <- stats.long %>%
  dplyr::group_by(IDs) %>%
  dplyr::mutate(median = round(median(enrichment, na.rm = TRUE), digits = 2),
                sd = round(sd(enrichment, na.rm = TRUE), digits = 2),
                variance = round(var(enrichment, na.rm = TRUE), digits = 2),
                min = round(min(enrichment, na.rm = TRUE), digits = 2),
                max = round(max(enrichment, na.rm = TRUE), digits = 2),
                prop.na = round(sum(is.na(enrichment))/length(cells), 
                                digits = 2)) %>%
  dplyr::ungroup()

# Residual's deciles. (mean)
res.decil.mean <- stats.long %>%
  dplyr::group_by(group.var) %>%
  dplyr::group_modify(~as.data.frame(t(quantile(.$residuals.mean, resm.cutoff, 
                                                na.rm = TRUE)))) %>%
  dplyr::ungroup()

res.decil.median <- stats.long %>%
  dplyr::group_by(group.var) %>%
  dplyr::group_modify(~as.data.frame(t(quantile(.$residuals.median, resm.cutoff, 
                                                na.rm = TRUE)))) %>%
  dplyr::ungroup()
colnames(res.decil.mean)[2:3] <- c("Pmin", "Pmax")
colnames(res.decil.median)[2:3] <- c("Pmin", "Pmax")

stats.long <- stats.long %>%
  dplyr::select(-cells, -enrichment) %>%
  unique() %>%
  dplyr::inner_join(res.decil.mean, by = "group.var") %>%
  dplyr::inner_join(res.decil.median, by = "group.var")

stats.long.annotated <- stats.long %>%
  dplyr::mutate(
    group.mean = dplyr::case_when(switch.point < sp.cutoff[1] & 
                               residuals.mean > Pmax.x ~ 
                               "TOP-HighSensitivity",
                             switch.point > sp.cutoff[4] & 
                               residuals.mean < Pmin.x ~ 
                               "TOP-LowSensitivity",
                             switch.point > sp.cutoff[2] & 
                               switch.point < sp.cutoff[3] & 
                               residuals.mean < Pmin.x ~ 
                               "TOP-Differential-LowSensitivity",
                             switch.point > sp.cutoff[2] &
                               switch.point < sp.cutoff[3] & 
                               residuals.mean > Pmax.x ~ 
                               "TOP-Differential-HighSensitivity",
                             TRUE ~ NA_character_),
    group.median = dplyr::case_when(switch.point < sp.cutoff[1] & 
                                    residuals.median > Pmax.y ~ 
                                    "TOP-HighSensitivity",
                                  switch.point > sp.cutoff[4] & 
                                    residuals.median < Pmin.y ~ 
                                    "TOP-LowSensitivity",
                                  switch.point > sp.cutoff[2] & 
                                    switch.point < sp.cutoff[3] & 
                                    residuals.median < Pmin.y ~ 
                                    "TOP-Differential-LowSensitivity",
                                  switch.point > sp.cutoff[2] &
                                    switch.point < sp.cutoff[3] & 
                                    residuals.median > Pmax.y ~ 
                                    "TOP-Differential-HighSensitivity",
                                  TRUE ~ NA_character_))

# Order.
rank <- stats.long.annotated %>%
  dplyr::mutate(in.range = switch.point > sp.cutoff[2] & 
                  switch.point < sp.cutoff[3],
                sp.rank = switch.point * as.numeric(in.range)) %>%
  dplyr::select(IDs, group.var, sp.rank, residuals.mean, in.range) %>%
  unique() %>%
  dplyr::group_split(group.var)
rank <- lapply(rank, FUN = function(x) {
  dt <- data.table::as.data.table(x)
  dt[, rank := data.table::frank(dt, -sp.rank, -residuals.mean, 
                                 ties.method = "dense")]
  return(dt)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(rank = dplyr::if_else(in.range, rank, NA_integer_)) %>%
  dplyr::select(IDs, group.var, rank) %>%
  unique()

rank.median <- stats.long.annotated %>%
  dplyr::mutate(in.range = switch.point > sp.cutoff[2] & 
                  switch.point < sp.cutoff[3],
                sp.rank = switch.point * as.numeric(in.range)) %>%
  dplyr::select(IDs, group.var, sp.rank, residuals.median, in.range) %>%
  unique() %>%
  dplyr::group_split(group.var)
rank.median <- lapply(rank.median, FUN = function(x) {
  dt <- data.table::as.data.table(x)
  dt[, rank.median := data.table::frank(dt, -sp.rank, -residuals.median, 
                                 ties.method = "dense")]
  return(dt)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(rank.median = dplyr::if_else(in.range, rank.median, NA_integer_)) %>%
  dplyr::select(IDs, group.var, rank.median) %>%
  unique()

stats.long.ranked <- stats.long.annotated %>%
  dplyr::inner_join(rank, by = c("IDs", "group.var")) %>%
  dplyr::inner_join(rank.median, by = c("IDs", "group.var"))

final.stats <- stats.long.ranked %>%
  dplyr::select(IDs, group.var, all_of(cols.stats)) %>%
  unique() %>%
  tidyr::pivot_wider(names_from = group.var, values_from = all_of(cols.stats),
                     names_sep = ".")


# Add Drug name and MoA to final.stats.
info <- drugInfo$IDs %>%
  dplyr::filter(IDs %in% final.stats$IDs) %>%
  dplyr::select(IDs, preferred.drug.names, studies) %>%
  dplyr::left_join(y = drugInfo$MoAs[, c("IDs", "MoAs")], by = "IDs",
                   relationship = "many-to-many") %>%
  dplyr::left_join(y = drugInfo$Targets, by = "IDs",
                   relationship = "many-to-many") %>%
  dplyr::left_join(y = drugInfo$Synonyms, by = "IDs",
                   relationship = "many-to-many")
if (dim(info)[1] > 0) {
  info <- aggregate(.~ IDs, data = info, na.action = NULL, FUN = function(x) {
    paste(na.omit(unique(x)), collapse = "; ")
  })
  cols.druginfo <- c("drugs", "preferred.drug.names", "MoAs", "targets", 
                     "studies")
} else {
  info <- data.frame(IDs = rownames(bc.recomputed@normalized))
  cols.druginfo <- NULL
}

final.stats <- final.stats %>%
  #dplyr::left_join(info, by = "IDs") %>%
  tibble::column_to_rownames("IDs") %>%
  unique()

# Order by rank and reorder columns.
final.stats <- final.stats[order(final.stats[, order.col], decreasing = FALSE),
                           c(cols.stats.level)]


# Add to beyondcell object.
bc.prueba <- bc.recomputed
bc.prueba@ranks[[idents]] <- final.stats

bc.prueba@ranks$TCs_res.0.3


