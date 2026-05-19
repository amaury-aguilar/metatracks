#boundariesNpeaks_meta.v14.r
#=============================================================
#------------------------------Libraries----------------------
#=============================================================
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(parallel)
library(scales)
library(patchwork)
library(yaml)
library(rtracklayer)
library(matrixStats)
#=============================================================
#--------------------------------Env Configuration------------
#=============================================================
start.time <- Sys.time()
options(scipen=999)
dir_path <- getwd()
setwd(dir_path)

config <- yaml::read_yaml("config/metatracks_config.yaml")
source("config/metatracks_functions.r")
#=============================================================
# GLOBAL
#=============================================================
(density_plot_width <- config$global$density_plot_width)
(ncores <- config$global$ncores)
(alpha_level <- config$global$alpha_level)
(output_dir <- config$global$output_dir)
(chr_prefix <- config$global$chromosome_prefix_in_bw)
#=============================================================
# DIMENSIONS
#=============================================================
(nbins <- config$global$dimensions$nbins)

(in_boundaries_pb_up <- config$global$dimensions$in_boundaries_pb_up)
(in_boundaries_pb_down <- config$global$dimensions$in_boundaries_pb_down)

(near_boundaries_pb_up <- config$global$dimensions$near_boundaries_pb_up)
(near_boundaries_pb_down <- config$global$dimensions$near_boundaries_pb_down)

(out_boundaries_pb_up <- config$global$dimensions$out_boundaries_pb_up)
(out_boundaries_pb_down <- config$global$dimensions$out_boundaries_pb_down)

(anchor_bin_size <- config$global$dimensions$anchor_bin_size)

(anchor_pb <- out_boundaries_pb_up + near_boundaries_pb_up + in_boundaries_pb_up + in_boundaries_pb_down + near_boundaries_pb_down + out_boundaries_pb_down)

(resol <- config$global$dimensions$resol)
(TAD_extend <- config$global$dimensions$TAD_extend)
(distance_To_Nearest_TAD <- config$global$dimensions$distance_To_Nearest_TAD)

#=============================================================
# OUTPUT
#=============================================================
(formato <- config$global$output$formato)
(row_height <- config$global$output$row_height)
(base_height <- config$global$output$base_height)
#=============================================================
# AGGREGATION
#=============================================================
(aggreagation_stat_label <- config$global$aggreagation_stat$label)
(aggreagation_stat_method <- config$global$aggreagation_stat$method)
#=============================================================
#------------------------------Make_directories---------------
#=============================================================
dir.create(output_dir)
#output_sample_dir <- paste0(output_dir, "/", sample_group)
#dir.create(output_sample_dir)
#=============================================================
#----------------- Set pb based bins limits ------------------
#=============================================================
(in_boundaries_up_bins <-
  in_boundaries_pb_up / anchor_bin_size)
(in_boundaries_down_bins <-
  in_boundaries_pb_down / anchor_bin_size)
(near_boundaries_up_bins <-
  near_boundaries_pb_up / anchor_bin_size)
(near_boundaries_down_bins <-
  near_boundaries_pb_down / anchor_bin_size)
(out_boundaries_up_bins <-
  out_boundaries_pb_up / anchor_bin_size)
(out_boundaries_down_bins <-
  out_boundaries_pb_down / anchor_bin_size)

(near_boundaries_up_start_d <-
  out_boundaries_up_bins + near_boundaries_up_bins)
(in_boundaries_up_start_d <-
  near_boundaries_up_start_d + in_boundaries_up_bins)
(in_boundaries_down_start_d <-
  in_boundaries_up_start_d + in_boundaries_down_bins)
(near_boundaries_down_start_d <-
  in_boundaries_down_start_d + near_boundaries_down_bins)
(out_boundaries_down_start_d <-
  near_boundaries_down_start_d + out_boundaries_down_bins)
(out_boundaries_up_end_d <-
  out_boundaries_down_start_d + out_boundaries_up_bins)
(near_boundaries_up_end_d <-
  out_boundaries_up_end_d + near_boundaries_up_bins)
(in_boundaries_up_end_d <-
  near_boundaries_up_end_d + in_boundaries_up_bins)
(in_boundaries_down_end_d <-
  in_boundaries_up_end_d + in_boundaries_down_bins)
(near_boundaries_down_end_d <-
  in_boundaries_down_end_d + near_boundaries_down_bins)
(out_boundaries_down_end_d <-
  near_boundaries_down_end_d + out_boundaries_down_bins)
#=============================================================
(tad_center <- nbins / 2)
(ext_bins <- ceiling(nbins * TAD_extend / (1 + 2 * TAD_extend)))
(mid_bins <- nbins - 2 * ext_bins)
#=============================================================
#---------------------------- Plot Axis and theme ------------
#=============================================================
axis_tads_ls <- list(
  scale_x_continuous(
    limits = c(1, nbins),
    breaks = c(1, ext_bins, ext_bins + mid_bins, nbins),
    labels = c(paste0("-", TAD_extend, "X"), paste0("Up"),
    paste0("Down"), paste0("+", TAD_extend, "X")),
    expand = c(0, 0)),
  geom_vline(
    xintercept = ext_bins,
    linewidth = 0.7,
    linetype = "dashed",
    alpha = 0.6),
  geom_vline(
    xintercept = ext_bins + mid_bins,
    linewidth = 0.7,
    linetype = "dashed",
    alpha = 0.6))

axis_anchors_ls <- list(
  scale_x_continuous(
    limits = c(1, out_boundaries_down_end_d),
    breaks = c(
      1,
      out_boundaries_up_bins,
      near_boundaries_up_start_d,
      in_boundaries_up_start_d,
      in_boundaries_down_start_d,
      near_boundaries_down_start_d,
      out_boundaries_down_start_d,
      out_boundaries_up_end_d,
      near_boundaries_up_end_d,
      in_boundaries_up_end_d,
      in_boundaries_down_end_d,
      near_boundaries_down_end_d,
      out_boundaries_down_end_d),
    labels = c(
      paste0(in_boundaries_up_bins +
      near_boundaries_up_bins +
      out_boundaries_up_bins),
      paste0(in_boundaries_up_bins +
      near_boundaries_up_bins),
      paste0(in_boundaries_up_bins),
      c(paste0("0", "\n", "Up")),
      paste0(in_boundaries_down_bins),
      paste0(in_boundaries_down_bins +
      near_boundaries_down_bins),
      paste0(in_boundaries_down_bins +
      near_boundaries_down_bins +
      out_boundaries_down_bins),
      paste0(in_boundaries_up_bins +
      near_boundaries_up_bins),
      paste0(in_boundaries_up_bins),
      c(paste0("0", "\n", "Down")),
      paste0(in_boundaries_down_bins),
      paste0(in_boundaries_down_bins +
      near_boundaries_down_bins),
      paste0(in_boundaries_down_bins +
      near_boundaries_down_bins +
      out_boundaries_down_bins)
      ),
    expand = c(0, 0)),
  geom_vline(
    xintercept = c(
      near_boundaries_up_start_d,
      in_boundaries_down_start_d,
      near_boundaries_up_end_d,
      in_boundaries_down_end_d),
    linewidth = 0.5,
    linetype = "dashed",
    alpha = 0.6),
  geom_vline(
    xintercept = out_boundaries_down_start_d,
    linewidth = 0.5,
    linetype = "dotted",
    alpha = 0.25))
#
plot_theme <- list(
  theme_minimal(),
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14, angle = 0, hjust = 0.75),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14, angle = 0, hjust = 0.5),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.line = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length.x = unit(3, "pt"),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length.y = unit(3, "pt"),
    panel.border = element_rect(linetype = "solid", fill = NA),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.title = element_text(hjust = 0.5)))

remove_x <- theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank())

#=============================================================
#-------------- Iteration over cluster --------------------
#=============================================================
anchor_plots_strips_ls <- list()
for (c in names(config$samples)) {
  sample_group <- c
#=============================================================
#---------------------------- Genes---------------------------
#=============================================================
(genes <- config$regions[[2]][1])
(genes_cols <- strsplit(config$regions[[2]][2], ",")[[1]])
(genes_names <- names(config$regions)[[2]])
genes_dt <- fread(genes, col.names = genes_cols)
genes_dt$chr <- gsub("chr", "", genes_dt$chr)
grouping_col <- config$regions[[2]][3]
genes_cluster <- genes_dt[get(grouping_col) %in% sample_group]

(genes_ctrl <- config$regions[[3]][1])
(genes_ctrl_cols <- strsplit(config$regions[[3]][2], ",")[[1]])
(genes_ctrl_names <- names(config$regions)[[3]])
genes_ctrl_dt <- fread(genes_ctrl,  col.names = genes_ctrl_cols)
genes_ctrl_dt$chr <- gsub("chr", "", genes_ctrl_dt$chr)
genes_rndm <- genes_ctrl_dt[!id %in% genes_dt$id]

genes_ls <- list()
genes_ls[[genes_names]] <- genes_cluster
genes_ls[[genes_ctrl_names]] <- genes_rndm

#=============================================================
#--------------------------- Domains -------------------------
#=============================================================
(domains <- config$regions[[1]][1])
(domains_cols <- strsplit(config$regions[[1]][2], ",")[[1]])
(domains_names <- names(config$regions)[[1]])
domains_bedpe_dt <- fread(domains, col.names = domains_cols)

domains_dt <- data.table(
  chr = gsub("chr", "", domains_bedpe_dt$chr1),
  start = floor(domains_bedpe_dt$end1 / as.numeric(resol)) * resol,
  end = floor(domains_bedpe_dt$start2 / as.numeric(resol)) * resol,
  id = domains_bedpe_dt$id
)
domains_dt$length <- domains_dt$end - domains_dt$start
domains_dt$extension <- domains_dt$length * TAD_extend

domains_ext_dt <- data.table(
  chr = domains_dt$chr,
  start = domains_dt$start - domains_dt$extension,
  end = domains_dt$end + domains_dt$extension,
  id = domains_dt$id,
  upBoundary = domains_dt$start,
  downBoundary = domains_dt$end
)
#=================================================================
#--------------- Filtering domains overlapping genes --------------
#=================================================================
domains_overlapping_genes_ls <- overlaps_f(genes_ls, domains_ext_dt)

domains_filtered_gr <- GRanges(
  seqnames = paste0(chr_prefix, domains_overlapping_genes_ls$DEGs$chrDomains),
  ranges = IRanges(
    start = domains_overlapping_genes_ls$DEGs$startDomains,
    end = domains_overlapping_genes_ls$DEGs$endDomains
  ),
  domain = domains_overlapping_genes_ls$DEGs$domain
)

domains_filtered_ctrl_gr <- GRanges(
  seqnames = paste0(chr_prefix, domains_overlapping_genes_ls$Ctrl$chrDomains),
  ranges = IRanges(
    start = domains_overlapping_genes_ls$Ctrl$startDomains,
    end = domains_overlapping_genes_ls$Ctrl$endDomains
  ),
  domain = domains_overlapping_genes_ls$Ctrl$domain
)

domains_filtered_gr <- unique(domains_filtered_gr)
domains_filtered_ctrl_gr <- unique(domains_filtered_ctrl_gr)

domains_overlapping_genes_gr_ls <- list()
domains_overlapping_genes_gr_ls[[genes_names]] <-domains_filtered_gr
domains_overlapping_genes_gr_ls[[genes_ctrl_names]] <-domains_filtered_ctrl_gr

#=================================================================
#--------------------------- Anchors -----------------------------
#=================================================================
anchors_bed_dt <- data.table(
  chr = domains_bedpe_dt$chr1,
  start = (domains_bedpe_dt$start1 + domains_bedpe_dt$end1) / 2,
  end = (domains_bedpe_dt$start2 + domains_bedpe_dt$end2) / 2,
  id = domains_dt$id
)
anchors_up_dt <- data.table(
  chr = anchors_bed_dt$chr,
  start = anchors_bed_dt$start - in_boundaries_pb_up,
  end = anchors_bed_dt$start + in_boundaries_pb_down,
  id = anchors_bed_dt$id
)
anchors_down_dt <- data.table(
  chr = anchors_bed_dt$chr,
  start = anchors_bed_dt$end - in_boundaries_pb_up,
  end = anchors_bed_dt$end + in_boundaries_pb_down,
  id = anchors_bed_dt$id
)

up_anchors_overlapping_genes_ls <-
  overlaps_f(genes_ls, anchors_up_dt)
down_anchors_overlapping_genes_ls <-
  overlaps_f(genes_ls, anchors_down_dt)

up_gene_anchors <- up_anchors_overlapping_genes_ls$DEGs$domain
down_gene_anchors <- down_anchors_overlapping_genes_ls$DEGs$domain
anchors_overlapping_genes <-
  unique(c(up_gene_anchors, down_gene_anchors))

up_gene_ctrl_anchors <- up_anchors_overlapping_genes_ls$Ctrl$domain
down_gene_ctrl_anchors <- down_anchors_overlapping_genes_ls$Ctrl$domain
anchors_overlapping_ctrl_genes <- unique(c(up_gene_ctrl_anchors, down_gene_ctrl_anchors))

anchors_overlapping_ctrl_genes <-setdiff(anchors_overlapping_ctrl_genes, anchors_overlapping_genes)

anchors_overlapping_ctrl_genes_same_size <-sample(anchors_overlapping_ctrl_genes,
length(anchors_overlapping_genes),replace = FALSE)

anchors_overlapping_genes_dt <-
  anchors_bed_dt[id %in% anchors_overlapping_genes]
anchors_overlapping_ctrl_genes_dt <-
  anchors_bed_dt[id %in% anchors_overlapping_ctrl_genes_same_size]

anchors_up_dt <- data.table(
  chr = paste0(chr_prefix, anchors_overlapping_genes_dt$chr),
  start = anchors_overlapping_genes_dt$start -
    in_boundaries_pb_up -
    near_boundaries_pb_up -
    out_boundaries_pb_up,
  end = anchors_overlapping_genes_dt$start +
    in_boundaries_pb_down +
    near_boundaries_pb_down +
    out_boundaries_pb_down,
  domain = anchors_overlapping_genes_dt$id
)
anchors_down_dt <- data.table(
  chr = paste0(chr_prefix, anchors_overlapping_genes_dt$chr),
  start = anchors_overlapping_genes_dt$end -
    in_boundaries_pb_up -
    near_boundaries_pb_up -
    out_boundaries_pb_up,
  end = anchors_overlapping_genes_dt$end +
    in_boundaries_pb_down +
    near_boundaries_pb_down +
    out_boundaries_pb_down,
  domain = anchors_overlapping_genes_dt$id
)

anchors_up_ctrl_dt <- data.table(
  chr = paste0(chr_prefix, anchors_overlapping_ctrl_genes_dt$chr),
  start = anchors_overlapping_ctrl_genes_dt$start -
    in_boundaries_pb_up -
    near_boundaries_pb_up -
    out_boundaries_pb_up,
  end = anchors_overlapping_ctrl_genes_dt$start +
    in_boundaries_pb_down +
    near_boundaries_pb_down +
    out_boundaries_pb_down,
  domain = anchors_overlapping_ctrl_genes_dt$id
)
anchors_down_ctrl_dt <- data.table(
  chr = paste0(chr_prefix, anchors_overlapping_ctrl_genes_dt$chr),
  start = anchors_overlapping_ctrl_genes_dt$end -
    in_boundaries_pb_up -
    near_boundaries_pb_up -
    out_boundaries_pb_up,
  end = anchors_overlapping_ctrl_genes_dt$end +
    in_boundaries_pb_down +
    near_boundaries_pb_down +
    out_boundaries_pb_down,
  domain = anchors_overlapping_ctrl_genes_dt$id
)

anchors_up_gr <- GRanges(anchors_up_dt)
anchors_down_gr <- GRanges(anchors_down_dt)
anchors_up_ctrl_gr <- GRanges(anchors_up_ctrl_dt)
anchors_down_ctrl_gr <- GRanges(anchors_down_ctrl_dt)
#=================================================================
#-------------------------- Bigwigs ------------------------------
#=================================================================
bw_files <- config$marks$files
marks_colors <- config$marks$colors
#=================================================================
#----------------- Import Bigwig Signal for Domains---------------
#=================================================================
all_domains <- domains_overlapping_genes_gr_ls[[genes_names]]
# Split
domain_split <- split(all_domains, all_domains$domain)
# tile once
bins_ls <- lapply(domain_split, function(gr){
  tile(gr, n = config$global$dimensions$nbins)[[1]]
})
# combine all bins
all_bins <- unlist(GRangesList(bins_ls), use.names = FALSE)
# store domain labels
domain_id <- rep(names(bins_ls), lengths(bins_ls))
bw_signal_ls <- list()
for(mark in names(bw_files)){
  cat("Processing", mark, "\n")
  # import ONCE
  bw <- import(bw_files[[mark]], which = reduce(all_domains))
  # overlaps all bins at once
  hits <- findOverlaps(all_bins, bw, ignore.strand = TRUE)
  # initialize
  signal_vec <- numeric(length(all_bins))
  # aggregate means
  dt <- data.table(
    bin = queryHits(hits),
    score = mcols(bw)$score[subjectHits(hits)]
  )

  agg_dt <- dt[, .(signal = get(aggreagation_stat_method)(score)), by = bin]

  signal_vec[agg_dt$bin] <- agg_dt$signal
  # convert to matrix by domain
  signal_mat <- matrix(signal_vec,
  nrow = config$global$dimensions$nbins,
  byrow = FALSE
  )
  # agg across domains
  signal_v <- switch(
    aggreagation_stat_method,
    mean = rowMeans(signal_mat, na.rm = TRUE),
    median = rowMedians(signal_mat, na.rm = TRUE),
    sum = rowSums(signal_mat, na.rm = TRUE),
    stop("Unsupported aggregation method")
  )

  signal_dt <- data.table(
    bin = 1:config$global$dimensions$nbins,
    signal = signal_v
  )
  bw_signal_ls[[mark]] <- signal_dt
}


#=============================================================
#---------------- Import Bigwig Signal for anchors -----------
#=============================================================
anchor_bins <- anchor_pb / anchor_bin_size
up_anchor_bins <- seq(1, anchor_bins, 1)
down_anchor_bins <- seq(anchor_bins + 1, anchor_bins * 2, 1)
# Split
up_anchors_split <- split(anchors_up_gr, anchors_up_gr$domain)
down_anchors_split <-
split(anchors_down_gr, anchors_down_gr$domain)
# tile once
bins_ls <- lapply(names(up_anchors_split), function(id){

  up_gr <- up_anchors_split[[id]]
  down_gr <- down_anchors_split[[id]]

  up_bins <- tile(up_gr, n = anchor_bins)[[1]]
  down_bins <- tile(down_gr, n = anchor_bins)[[1]]

  up_bins$domain <- id
  down_bins$domain <- id

  up_bins$anchor <- "up"
  down_bins$anchor <- "down"

  up_bins$bin <- up_anchor_bins
  down_bins$bin <- down_anchor_bins

  c(up_bins, down_bins)
})

names(bins_ls) <- names(up_anchors_split)
# combine all bins
all_bins <- unlist(GRangesList(bins_ls), use.names = FALSE)

template_dt <- unique(data.table(
  domain = all_bins$domain,
  bin = all_bins$bin
))

anchor_signal_by_mark_ls <- list()
for (mark in names(bw_files)) {
  cat("Processing", mark, "\n")
  # import ONCE
  bw <- import(bw_files[[mark]], which = reduce(all_bins))
  # overlaps all bins at once
  hits <- findOverlaps(all_bins, bw, ignore.strand = TRUE)
  #
  dt <- data.table(
    domain = all_bins$domain[queryHits(hits)],
    bin = all_bins$bin[queryHits(hits)],
    score = mcols(bw)$score[subjectHits(hits)]
    )

  signal_dt <- dt[,
    .(signal = get(aggreagation_stat_method)(score, na.rm = TRUE)), by = .(domain, bin)]
  signal_dt <- merge(template_dt, signal_dt,
    by = c("domain", "bin"), all.x = TRUE)

  signal_dt[is.na(signal), signal := 0]

  anchor_signal_by_mark_ls[[mark]] <- signal_dt
}

anchors_m_signal_by_mark_ls <-
lapply(anchor_signal_by_mark_ls, function(dt) {
  dt[, .(signal = get(aggreagation_stat_method)(signal, na.rm = TRUE)),
    by = bin][order(bin)]
})
#=================================================================
#---------------- Flip anchors ---------------
#=================================================================
up_domains <- unique(up_anchors_overlapping_genes_ls$DEGs$domain)
down_domains <-unique(down_anchors_overlapping_genes_ls$DEGs$domain)

# Domains to flip, present in down, but not in up
domains_to_flip <- setdiff(down_domains, up_domains)

anchor_signal_by_mark_flipped_ls <- lapply(
  anchor_signal_by_mark_ls,
  function(dt){
    dt <- copy(dt)
    # Flip only selected domains
    dt[domain %in% domains_to_flip,
     bin := max(bin) - bin + 1, by = domain]
    # Reorder after flipping
    setorder(dt, domain, bin)
    dt
  }
)

  anchors_agg_signal_by_mark_flipped_ls <-
  lapply(anchor_signal_by_mark_flipped_ls, function(dt) {
    dt[, .(signal = get(aggreagation_stat_method)(signal, na.rm = TRUE)),
      by = bin][order(bin)]
  })

  plot_anchors_ls <- Map(
    function(signal, color) {
      density_f(
        signal,
        "black",
        color,
        axis_anchors_ls,
        aggreagation_stat_label
      )
  }, anchors_agg_signal_by_mark_flipped_ls,
    config$marks$colors[names(anchors_agg_signal_by_mark_flipped_ls)]
  )

#=============================================================
#----------------  ---------------
#=============================================================

# saveRDS()
#=============================================================
#----------------  Plots by cluster ---------------
#=============================================================

plot_anchors_ls[-length(plot_anchors_ls)] <-
  lapply(plot_anchors_ls[-length(plot_anchors_ls)],
   function(p) p + remove_x)

anchor_plots_strips_ls[[sample_group]] <- wrap_plots(plot_anchors_ls, ncol = 1)

#wrap_plots(plot_anchors_strip, plot_anchors_strip, ncol = 2)
#graficar_density(plots_density_sel_wraped, name)
}

wrap_plots(anchor_plots_strips_ls, ncol = length(anchor_plots_strips_ls))
