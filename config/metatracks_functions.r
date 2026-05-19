
# metatracks_functions.r

overlaps_f <- function(list_peaks,regions){
  regions_gr <- GRanges(regions)
  overlaps_list <- list()
  for (i in names(list_peaks)){
    genes <- copy(list_peaks[[i]])
    genes_gr <- GRanges(genes)
    #---------------------------------------------------------------------------
    hits <- findOverlaps(genes_gr, regions_gr)
    #---------------------------------------------------------------------------
    hits_dt <- data.table(
      chrDomains = regions$chr[hits@to],
      startDomains = start(regions_gr)[hits@to],
      endDomains = end(regions_gr)[hits@to],
      domain = regions$id[hits@to],
      startGenes = start(genes_gr)[hits@from],
      endGenes = end(genes_gr)[hits@from],
      gene = genes$id[hits@from],
      score = genes$score[hits@from],
      upBoundary = regions$upBoundary[hits@to],
      downBoundary = regions$downBoundary[hits@to]
    )
    #-------------------------------------------------------------------------------
    overlaps_list[[i]] <- hits_dt
  }
  return(overlaps_list)
}

density_f <- function(
  dt,
  line_color,
  color,
  axis_ls,
  scale_label,
  ymin = NA,
  ymax = NA) {
  
  ymin <- if (is.na(ymin)) min(dt$signal, na.rm = TRUE) else ymin
  ymax <- if (is.na(ymax)) max(dt$signal, na.rm = TRUE) else ymax
  
  gradient <- scale_fill_gradient2(
    low = "white",
    high = color,
    limits = c(ymin, ymax),
    oob = scales::squish,
    name = scale_label
    )
  # gradient + 
  densityPlot <- ggplot(dt, aes(x = bin)) +
    geom_ribbon(aes(ymin = ymin, ymax = signal),
    fill = color, alpha = alpha_level) +
    geom_line(aes(y = signal),
    linewidth = 0.5, color = line_color) +
    plot_theme + axis_ls +
    scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0))
  
  return(densityPlot)
}

graficar_density <- function(plot, mark) {
  nplots <- length(plot)
  if (formato == "png") {
    path.apaplot <- paste0(dir_path_cluster,
        "/",
        nomen,
        ".",
        mark,
        ".",
        nbins,
        "_bins.png")
    device.plot <- png(
        path.apaplot,
        width = density_plot_width,
        height = 5 * nplots,
        units = "cm",
        res = 500
    )
  } else if (formato == "pdf") {
    path.apaplot <- paste0(
        dir_path_cluster,
        "/",
        nomen,
        ".",
        mark,
        ".",
        nbins,
        "_bins.pdf")
    device.plot <- pdf(
        path.apaplot,
        width = density_plot_width / 2.54,
        height = (5 * nplots) / 2.54)
  } else {
    stop("Image format not available. Use instead 'pdf' or 'png'")
  }
  path.apaplot
  device.plot
  print(plot)
  dev.off()
}
