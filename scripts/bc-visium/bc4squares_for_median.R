bc4Squares.median <- function(bc, idents, lvl = NULL, top = 3, topnames = NULL, 
                       force = 1, alpha = 0.7, pt.size = 3, 
                       x.cutoff = NULL, y.cutoff = NULL, ...) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check idents.
  if (length(idents) != 1) stop('Idents must be a single metadata column.')
  if (!(idents %in% names(bc@ranks))) {
    stop(paste0('$', idents, ' not found in bc@ranks.'))
  } else if (idents == "general") {
    stop('General rank can\'t be used in bc4Squares(). All residuals are 0.')
  }
  # Check lvl.
  if (is.null(lvl)) {
    lvl <- levels(as.factor(bc@meta.data[, idents]))
    in.lvl <- rep(TRUE, times = length(lvl))
  } else {
    in.lvl <- lvl %in% unique(bc@meta.data[, idents])
    if (all(!in.lvl)) {
      stop(paste0('None of the specified levels were found in ', idents, '.'))
    } else if (any(!in.lvl)) {
      warning(paste0('The following levels were not found in ', idents, ': ',
                     paste0(lvl[!in.lvl], collapse = ", "), '.'))
    }
  }
  # Check top.
  if (length(top) != 1 | top[1]%%1 != 0 | top[1] < 0) {
    stop('top must be a single integer >= 0.')
  }
  # Check topnames.
  if (!is.null(topnames)) {
    info <- drugInfo$IDs %>%
      dplyr::select(IDs, preferred.drug.names, studies) %>%
      dplyr::left_join(y = drugInfo$MoAs[, c("IDs", "MoAs")], by = "IDs") %>%
      dplyr::left_join(y = drugInfo$Targets, by = "IDs") %>%
      dplyr::left_join(y = drugInfo$Synonyms, by = "IDs") %>%
      dplyr::mutate(sources = studies) %>%
      as.data.frame()
    in.topnames <- toupper(topnames) %in% info$drugs |
      tolower(topnames) %in% info$IDs |
      toupper(topnames) %in% toupper(rownames(bc@normalized))
    if (all(!in.topnames)) {
      warning('None of the specified topname drugs were found in bc.')
    } else if (any(!in.topnames)) {
      warning(paste0('The following topname drugs were not found in bc: ' ,
                     paste0(topnames[!in.topnames], collapse = ", "), '.'))
    }
  } else {
    in.topnames <- NULL
  }
  # Check force.
  if (length(force) != 1 | force[1] < 0) {
    stop('force must be a single number >= 0.')
  }
  # Check alpha.
  if (length(alpha) != 1 | alpha[1] < 0 | alpha[1] > 1) {
    stop('alpha must be a single number between 0 and 1.')
  }
  # Check pt.size.
  if (length(pt.size) != 1 | !is.numeric(pt.size)) {
    stop('pt.size must be a single number.')
  }
  # Check x.cutoff.
  if (!is.null(x.cutoff)) {
    if (length(x.cutoff) != 2 | !is.numeric(x.cutoff)) {
      stop('x.cutoff must be a numeric vector of length 2.')
    }
    if (x.cutoff[2] < x.cutoff[1]) {
      warning(paste('Upper x cut-off is smaller than lower x cut-off.',
                    'Sorting x cut-offs in increasing order.'))
      x.cutoff <- sort(x.cutoff, decreasing = FALSE)
    }
  }
  # Check y.cutoff.
  if (!is.null(y.cutoff)) {
    if (length(y.cutoff) != 4 | !is.numeric(y.cutoff)) {
      stop('y.cutoff must be a numeric vector of length 4.')
    }
    if (any(y.cutoff < 0 | y.cutoff > 1)) {
      stop('y.cutoff must contain 4 switch point values between 0 and 1.')
    }
    sorted.y.cutoff <- sort(y.cutoff, decreasing = FALSE)
    if (!identical(y.cutoff, sorted.y.cutoff)) {
      warning(paste('Sorting y cut-offs in increasing order.'))
      y.cutoff <- sorted.y.cutoff
    }
  } else {
    y.cutoff <- c(0.1, 0.4, 0.6, 0.9)
  }
  # --- Code ---
  # Get info about drugs (their corresponding name in bc, the preferred name
  # used by beyondcell and the MoA).
  info <- FindDrugs(bc, x = rownames(bc@scaled), na.rm = FALSE)
  # Switch points.
  sp <- data.frame(switch.point = bc@switch.point[info$bc.names],
                   row.names = info$bc.names)
  # One plot per level.
  p4s <- lapply(lvl[in.lvl], function(l) {
    ### Subset residuals' medians and switch points.
    res <- bc@ranks[[idents]][info$bc.names,
                              paste0("residuals.median.", l), drop = FALSE]
    colnames(res) <- "residuals.median"
    df <- transform(merge(res, sp, by = 0), row.names = Row.names, Row.names = NULL)
    if (is.null(x.cutoff)) {
      ### Residual's deciles.
      res.decil <- quantile(as.numeric(res$residuals.median), na.rm = TRUE,
                            prob = seq(from = 0, to = 1, length = 11))
      x.cutoff <- c(res.decil[["10%"]], res.decil[["90%"]])
      x.cutoff.caption <- "first and last deciles"
    } else {
      residuals.limits <- c(min(res$residuals.median), max(res$residuals.median))
      if (x.cutoff[1] < residuals.limits[1]) {
        x.cutoff[1] <- residuals.limits[1]
        warning(paste0('For lvl = ', l, ': Lower x cut-off is outside of range.',
                       'Setting it equal to the minimum residual\'s median.'))
      }
      if (x.cutoff[2] > max(res$residuals.median)) {
        x.cutoff[2] <- residuals.limits[2]
        warning(paste0('For lvl = ', l, ': Upper x cut-off is outside of range.',
                       'Setting it equal to the maximum residual\'s median.'))
      }
      x.cutoff.caption <- paste0(x.cutoff, collapse = " and ")
    }
    ### y cut-off caption.
    y.cutoff.caption <- paste(paste0(y.cutoff[1:3], collapse = ", "), "and",
                              y.cutoff[4])
    ### Drug annotation.
    sp_cutoff1 <- as.numeric(df$switch.point) < y.cutoff[1]
    sp_cutoff2 <- as.numeric(df$switch.point) > y.cutoff[2]
    sp_cutoff3 <- as.numeric(df$switch.point) < y.cutoff[3]
    sp_cutoff4 <- as.numeric(df$switch.point) > y.cutoff[4]
    res_cutoff1 <- as.numeric(df$residuals.median) < x.cutoff[1]
    res_cutoff2 <- as.numeric(df$residuals.median) > x.cutoff[2]
    df$annotation <- rep("no", times = nrow(df))
    df$annotation[sp_cutoff1 & res_cutoff2] <- "TOP-HighSensitivityDrugs"
    df$annotation[sp_cutoff4 & res_cutoff1] <- "TOP-LowSensitivityDrugs"
    df$annotation[sp_cutoff2 & sp_cutoff3 &
                    res_cutoff1] <- "TOP-Differential-LowSensitivityDrugs"
    df$annotation[sp_cutoff2 & sp_cutoff3 &
                    res_cutoff2] <- "TOP-Differential-HighSensitivityDrugs"
    ### Drug labels.
    df$labels <- rep(NA, times = nrow(df))
    decreasing_order <- c("TOP-Differential-HighSensitivityDrugs",
                          "TOP-HighSensitivityDrugs")
    unique.annotations <- unique(df$annotation[df$annotation != "no"])
    sel.labels <- unlist(sapply(unique.annotations, function(x) {
      sub.df <- subset(df, subset = df$annotation == x)
      if (x %in% decreasing_order) {
        sub.df <- sub.df[order(sub.df$residuals.median, sub.df$switch.point,
                               decreasing = TRUE), ]
      } else {
        sub.df <- sub.df[order(sub.df$residuals.median, sub.df$switch.point,
                               decreasing = FALSE), ]
      }
      return(rownames(sub.df)[1:min(top, nrow(sub.df))])
    }))
    df[sel.labels, "labels"] <- info$preferred.and.sigs[match(sel.labels,
                                                              info$bc.names)]
    ### Topnames.
    if(length(topnames[in.topnames]) > 0) {
      topnames <- FindDrugs(bc, x = topnames[in.topnames], na.rm = FALSE)
      df[match(topnames$bc.names,
               table = rownames(df)), "labels"] <- topnames$preferred.and.sigs
    }
    ### Colours and names.
    colors <- c("#1D61F2", "#DA0078", "orange", "#C7A2F5", "grey80", "black")
    names <- c("TOP-LowSensitivityDrugs", "TOP-HighSensitivityDrugs",
               "TOP-Differential-HighSensitivityDrugs",
               "TOP-Differential-LowSensitivityDrugs", "no", "black")
    ### Circle's borders colour.
    df$borders <- df$annotation
    df$borders[df$labels != ""] <- "black"
    ### Reorder df so labeled points are plotted on top.
    df <- rbind(subset(df, subset = df$borders != "black"),
                subset(df, subset = df$borders == "black"))
    ### Plot.
    p <- ggplot(df, aes(x = as.numeric(residuals.median),
                        y = as.numeric(switch.point), color = borders,
                        fill = annotation)) +
      geom_point(shape = 21, alpha = alpha, size = pt.size) +
      scale_color_manual(values = setNames(colors, names)) +
      scale_fill_manual(values = setNames(colors, names), breaks = names[1:4],
                        drop = FALSE) + theme_classic() +
      geom_vline(xintercept = x.cutoff[1], linetype = "dotted") +
      geom_vline(xintercept = x.cutoff[2], linetype = "dotted") +
      geom_hline(yintercept = y.cutoff[4], linetype = "dotted") +
      geom_hline(yintercept = y.cutoff[1], linetype = "dotted") +
      geom_hline(yintercept = y.cutoff[2], linetype = "dotted") +
      geom_hline(yintercept = y.cutoff[3], linetype = "dotted") + ylim(0, 1) +
      labs(title = paste(idents, "=", l),
           caption = paste0("x cut-offs: ", x.cutoff.caption, "; y cut-offs: ",
                            y.cutoff.caption)) + 
      xlab("Residuals' median") + ylab("Switch Point") +
      ggrepel::geom_text_repel(label = df$labels, force = force, na.rm = TRUE,
                               ...) +
      guides(fill = guide_legend(title = "Drug Annotation"), color = "none") +
      cowplot::theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
    return(p)
  })
  # If p4s has only one plot, don't return a list.
  if (length(p4s) == 1) p4s <- p4s[[1]]
  return(p4s)
}

prueba4.median <- bc4Squares.median(bc.prueba, idents = "TCs_res.0.3")
wrap_plots(prueba4.median)
