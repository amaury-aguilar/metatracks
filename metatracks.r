#boundariesNpeaks_meta.v14.r
#------------------------------Libraries----------------------------------------
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(parallel)
library(scales)
library(patchwork)
library(yaml)
library(rtracklayer)
#--------------------------------Env Configuration------------------------------
start.time <- Sys.time()
options(scipen=999)
dir_path <- getwd()
setwd(dir_path)

config <- yaml::read_yaml("config/metatracks_config.yaml")
source("config/metatracks_functions.r")

#=================================================================
# GLOBAL
#=================================================================
(density_plot_width <- config$global$density_plot_width)
(ncores <- config$global$ncores)
(alpha_level <- config$global$alpha_level)
(normalize_by <- config$global$normalize_by)
(reference_bed <- config$global$reference_bed)
(segment_overlaps <- config$global$segment_overlaps)
(output_dir <- config$global$output_dir)
chr_prefix <- config$global$chromosome_prefix_in_bw

#=================================================================
# DIMENSIONS
#=================================================================
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

#=================================================================
# OUTPUT
#=================================================================
(formato <- config$global$output$formato)
(row_height <- config$global$output$row_height)
(base_height <- config$global$output$base_height)

#=================================================================
# AGGREGATION
#=================================================================
(aggreagation_stat_label <- config$global$aggreagation_stat$label)
(aggreagation_stat_method <- config$global$aggreagation_stat$method)
#=================================================================
#------------------------------Make_directories-------------------
#=================================================================
sample_group <- names(config$samples)[1]

dir.create(output_dir)
output_sample_dir <- paste0(output_dir, "/", sample_group)
dir.create(output_sample_dir)
#=================================================================
#----------------- Set pb based bins limits ----------
#=================================================================
(boundary_width_pb <- in_boundaries_pb_up + in_boundaries_pb_down)
(cross_boundaries_bins <- near_boundaries_pb_up / resol)
(outside_boundaries_bins <- out_boundaries_pb_up / resol)
(boundaries_bins_in <- in_boundaries_pb_up / resol)
(boundaries_bins_out <- in_boundaries_pb_down / resol)
(boundary_width_bin <- boundaries_bins_in + boundaries_bins_out)
(upBoundary_Start_bin <- outside_boundaries_bins)
(upBoundary_End_bin <- upBoundary_Start_bin + boundary_width_bin)
(middle_region_boundaries_bin <- upBoundary_Start_bin + boundary_width_bin + cross_boundaries_bins)
(downBoundary_Start_bin <- middle_region_boundaries_bin + cross_boundaries_bins)
(downBoundary_End_bin <- downBoundary_Start_bin + boundary_width_bin)
(end_region_boundaries_bin <- downBoundary_End_bin + outside_boundaries_bins)
tad_center <- nbins / 2
ext_bins <- ceiling(nbins * TAD_extend / (1 + 2 * TAD_extend))
mid_bins <- nbins - 2 * ext_bins
#=================================================================
#---------------------------- Plot Axis and theme ----------------
#=================================================================
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

in_boundaries_up_bins <-
  in_boundaries_pb_up / anchor_bin_size
in_boundaries_down_bins <-
  in_boundaries_pb_down / anchor_bin_size
near_boundaries_up_bins <-
  near_boundaries_pb_up / anchor_bin_size
near_boundaries_down_bins <-
  near_boundaries_pb_down / anchor_bin_size
out_boundaries_up_bins <-
  out_boundaries_pb_up / anchor_bin_size
out_boundaries_down_bins <-
  out_boundaries_pb_down / anchor_bin_size

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
#=================================================================
#---------------------------- Genes-------------------------------
#=================================================================
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

#=================================================================
#--------------------------- Domains -----------------------------
#=================================================================
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
  start = anchors_overlapping_genes_dt$start - in_boundaries_pb_up - near_boundaries_pb_up - out_boundaries_pb_up,
  end = anchors_overlapping_genes_dt$start + in_boundaries_pb_down + near_boundaries_pb_down + out_boundaries_pb_down,
  domain = anchors_overlapping_genes_dt$id
)
anchors_down_dt <- data.table(
  chr = paste0(chr_prefix, anchors_overlapping_genes_dt$chr),
  start = anchors_overlapping_genes_dt$end - in_boundaries_pb_up -near_boundaries_pb_up - out_boundaries_pb_up,
  end = anchors_overlapping_genes_dt$end + in_boundaries_pb_down + near_boundaries_pb_down + out_boundaries_pb_down,
  domain = anchors_overlapping_genes_dt$id
)

anchors_up_ctrl_dt <- data.table(
  chr = paste0(chr_prefix, anchors_overlapping_ctrl_genes_dt$chr),
  start = anchors_overlapping_ctrl_genes_dt$start - in_boundaries_pb_up - near_boundaries_pb_up - out_boundaries_pb_up,
  end = anchors_overlapping_ctrl_genes_dt$start + in_boundaries_pb_down + near_boundaries_pb_down + out_boundaries_pb_down,
  domain = anchors_overlapping_ctrl_genes_dt$id
)
anchors_down_ctrl_dt <- data.table(
  chr = paste0(chr_prefix, anchors_overlapping_ctrl_genes_dt$chr),
  start = anchors_overlapping_ctrl_genes_dt$end - in_boundaries_pb_up -near_boundaries_pb_up - out_boundaries_pb_up,
  end = anchors_overlapping_ctrl_genes_dt$end + in_boundaries_pb_down + near_boundaries_pb_down + out_boundaries_pb_down,
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
bw_signal_ls <- list()
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

  mean_dt <- dt[, .(signal = mean(score)), by = bin]

  signal_vec[mean_dt$bin] <- mean_dt$signal
  # convert to matrix by domain
  signal_mat <- matrix(
    signal_vec,
    nrow = config$global$dimensions$nbins,
    byrow = FALSE
  )
  # mean across domains
  signal_v <- rowMeans(signal_mat)
  signal_dt <- data.table(
    bin = 1:config$global$dimensions$nbins,
    signal = signal_v
  )
  bw_signal_ls[[mark]] <- signal_dt
}

plot_test <- density_f(
  bw_signal_ls$H3K4me1,
  "black",
  config$marks$colors$H3K4me1,
  axis_tads_ls)
plot_test

#=================================================================
#---------------- Import Bigwig Signal for anchors ---------------
#=================================================================

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
    score = mcols(bw)$score[subjectHits(hits)])

  signal_dt <- dt[,
    .(signal = mean(score, na.rm = TRUE)), by = .(domain, bin)]
  signal_dt <- merge(template_dt, signal_dt,
    by = c("domain", "bin"), all.x = TRUE)

  signal_dt[is.na(signal), signal := 0]

  anchor_signal_by_mark_ls[[mark]] <- signal_dt
}

anchor_metaprofile_ls <- lapply(anchor_signal_by_mark_ls, function(dt) {
  dt[, .(signal = mean(signal, na.rm = TRUE)),
    by = bin][order(bin)]
})
#----------------------------------------------------

plot_test <- density_f(
  anchor_metaprofile_ls$H3K4me1,
  "black",
  config$marks$colors$H3K4me1,
  axis_anchors_ls)
plot_test

ggplot(anchor_metaprofile_ls$H3K27me3, aes(x = bin)) +
    geom_ribbon(aes(ymin = 0, ymax = signal, fill = signal),
    alpha = alpha_level) +
    geom_line(aes(y = signal, linewidth = 0.5), color = "black") +
    scale_y_continuous(expand = c(0, 0)) + plot_theme + axis_anchors_ls

#
#------------------------------ Compute Fisher Statistics --------------------
#
DEGs_inTADs_InsideOut_dt <-
  identify_InsideOut_f(domains_overlapping_genes_ls[[genes_names]])
genesRndm_inTADs_InsideOut_dt <-
  identify_InsideOut_f(domains_overlapping_genes_ls[[genes_ctrl_names]])

DEGs_inTADs_InsideOut_flip_vectors <- 
  flip_regions_f(DEGs_inTADs_InsideOut_dt,mode="TAD")
setkey(DEGs_inTADs_InsideOut_dt, region)
setkey(DEGs_inTADs_InsideOut_flip_vectors, region)
DEGs_inTADs_InsideOut_dt[
  DEGs_inTADs_InsideOut_flip_vectors[,.(region,flip)],
  flip := i.flip,
  on = "region"
]

genesRndm_inTADs_InsideOut_vectors <- 
  flip_regions_f(genesRndm_inTADs_InsideOut_dt,mode="TAD")
setkey(genesRndm_inTADs_InsideOut_dt, region)
setkey(genesRndm_inTADs_InsideOut_vectors, region)
genesRndm_inTADs_InsideOut_dt[
  genesRndm_inTADs_InsideOut_vectors[,.(region,flip)],
  flip := i.flip,
  on = "region"
]

DEGs_fisher <- run_fisher_fast(200,DEGs_inTADs_InsideOut_dt,
                                genesRndm_inTADs_InsideOut_dt,
                                "inside_domain",
                                "inTAD","outTAD","peak",
                                ncores)

DEGs_fisher_anchors <- run_fisher_fast(200,DEGs_inTADs_InsideOut_dt,
                                        genesRndm_inTADs_InsideOut_dt,
                                        "inside_boundaries",
                                        "inBoundaries","outBoundaries",
                                        "peak",ncores)

DEGs_fisher_dt <- DEGs_fisher$results
DEGs_fisher_anchors_dt <- DEGs_fisher_anchors$results
odds_ratio_DEGs <- 
  plot_fisher_f(
    DEGs_fisher$results[,.(iter,p_value,odds_ratio)],
    "odds_ratio","Odds Ratio","DEGs")
p_value_DEGs <- 
  plot_fisher_f(
    DEGs_fisher$results[,.(iter,p_value,odds_ratio)],
    "p_value","p value","DEGs")
odds_ratio_DEGs_anchors <- 
  plot_fisher_f(
    DEGs_fisher_anchors$results[,.(iter,p_value,odds_ratio)],
    "odds_ratio","Odds Ratio","DEGs")
p_value_DEGs_anchors <- 
  plot_fisher_f(
    DEGs_fisher_anchors$results[,.(iter,p_value,odds_ratio)],
    "p_value","p value","DEGs")

saveRDS(odds_ratio_DEGs,
        file = paste0(dir_path_cluster,"/odds_ratio_DEGs"))
saveRDS(p_value_DEGs,
        file = paste0(dir_path_cluster,"/p_value_DEGs"))
    
#
#----------------- Re-sample Rndm genes for density plots------------------
#
p_values <- DEGs_fisher$results$p_value
p_values_sorted <- sort(p_values)
middle_p_value <- p_values_sorted[ceiling(length(p_values_sorted)/2)]
rndm_DEGs_dt <- DEGs_fisher$results[p_value==middle_p_value][sample(.N,1)]

setkey(genesRndm_inTADs_InsideOut_dt,peak)
genesRndm_inTADs_InsideOut_sample_dt <- 
  genesRndm_inTADs_InsideOut_dt[rndm_DEGs_dt$sampled_ids]

(n_degs_peaks <- uniqueN(DEGs_inTADs_InsideOut_dt$peak))
(n_tads_peaks <- uniqueN(DEGs_inTADs_InsideOut_dt$region))
(n_Rndm_degs <- uniqueN(genesRndm_inTADs_InsideOut_sample_dt$peak))
(n_Rndm_tads <- uniqueN(genesRndm_inTADs_InsideOut_sample_dt$region))
#
#------------------------- Plot density plots --------------------
#  
DEGs_inTADs_bins <- binmaker_redim_f(DEGs_inTADs_InsideOut_dt,"score",
                                      aggreagation_stat[2],nbins,FALSE)

genesRndm_inTADs_InsideOut_sample_dt[,log_TPMs:= log10(score+1)]
rndmGenes_inTADs_bins <- binmaker_redim_f(genesRndm_inTADs_InsideOut_sample_dt,"log_TPMs",
                                          aggreagation_stat[2],nbins,FALSE)
#
#----------------------- Plot ANCHOR density plots -------------
#
segment_overlaps_ls <- c("farOutside_Up_left",
                          "outside_Up_left",
                          "inside_UpBoundary",
                          "crossing_Up_right",
                          "farCrossing_Up_right",
                          "farCrossing_Down_left",
                          "crossing_Down_left",
                          "inside_DownBoundary",
                          "outside_Down_right",
                          "farOutside_Down_right")
segment_overlaps_ls <- c("inside_UpBoundary",
                          "inside_DownBoundary")

DEGs_inTADs_InsideOut_anchors_ls <- 
  binmaker_pb_f(DEGs_inTADs_InsideOut_dt,
                "score",aggreagation_stat[2],
                end_region_boundaries_bin,
                segment_overlaps_ls,FALSE)
DEGs_inTADs_InsideOut_anchors_bins <- DEGs_inTADs_InsideOut_anchors_ls$bins

genesRndm_inTADs_InsideOut_sample_anchors_ls <- 
  binmaker_pb_f(genesRndm_inTADs_InsideOut_sample_dt,
                "log_TPMs",aggreagation_stat[2],
                end_region_boundaries_bin,
                segment_overlaps_ls,FALSE)
genesRndm_inTADs_InsideOut_sample_anchors_bins <-
  genesRndm_inTADs_InsideOut_sample_anchors_ls$bins

(n_degs_peaks_anchors <- DEGs_inTADs_InsideOut_anchors_ls$peaks)
(n_tads_peaks_anchors <- DEGs_inTADs_InsideOut_anchors_ls$regions)
(n_Rndm_degs_anchors <- genesRndm_inTADs_InsideOut_sample_anchors_ls$peaks)
(n_Rndm_tads_anchors <- genesRndm_inTADs_InsideOut_sample_anchors_ls$regions)
#
#---------------------- Flip DEGs density plots ---------------
#
DEGs_inTADs_InsideOut_FLIPPED_bins <- binmaker_redim_f(
  DEGs_inTADs_InsideOut_dt,"score",
  aggreagation_stat[2],nbins,TRUE)

genesRndm_inTADs_InsideOut_sample_FLIPPED_bins <- binmaker_redim_f(
  genesRndm_inTADs_InsideOut_sample_dt,"log_TPMs",
  aggreagation_stat[2],nbins,TRUE)
#
#---------------- Flip DEGs for ANCHOR density plots --------
#
DEGs_inTADs_InsideOut_anchors_FLIPPED_dt <- copy(DEGs_inTADs_InsideOut_dt)
DEGs_inTADs_InsideOut_anchors_FLIPPED_dt$flip <- NULL
DEGs_inTADs_InsideOut_anchors_flip_vectors <- 
  flip_regions_f(DEGs_inTADs_InsideOut_anchors_FLIPPED_dt,mode="boundary")
setkey(DEGs_inTADs_InsideOut_anchors_FLIPPED_dt, region)
setkey(DEGs_inTADs_InsideOut_anchors_flip_vectors, region)
DEGs_inTADs_InsideOut_anchors_FLIPPED_dt[
  DEGs_inTADs_InsideOut_anchors_flip_vectors[,.(region,flip)],
  flip := i.flip,
  on = "region"
]

DEGs_inTADs_InsideOut_FLIPPED_bins_anchors_ls <- binmaker_pb_f(
  DEGs_inTADs_InsideOut_anchors_FLIPPED_dt,"score",
  aggreagation_stat[2],end_region_boundaries_bin,
  segment_overlaps_ls,TRUE)
DEGs_inTADs_InsideOut_FLIPPED_bins_anchors_dt <- 
  DEGs_inTADs_InsideOut_FLIPPED_bins_anchors_ls$bins

genesRndm_inTADs_InsideOut_sample_anchors_FLIPPED_dt <- 
  copy(genesRndm_inTADs_InsideOut_sample_dt)
genesRndm_inTADs_InsideOut_sample_anchors_FLIPPED_dt$flip <- NULL
genesRndm_inTADs_InsideOut_sample_anchors_vectors <- 
  flip_regions_f(genesRndm_inTADs_InsideOut_sample_anchors_FLIPPED_dt,mode="boundary")
setkey(genesRndm_inTADs_InsideOut_sample_anchors_FLIPPED_dt, region)
setkey(genesRndm_inTADs_InsideOut_sample_anchors_vectors, region)
genesRndm_inTADs_InsideOut_sample_anchors_FLIPPED_dt[
  genesRndm_inTADs_InsideOut_sample_anchors_vectors[,.(region,flip)],
  flip := i.flip,
  on = "region"
]

genesRndm_inTADs_InsideOut_sample_anchors_ls <- 
  binmaker_pb_f(genesRndm_inTADs_InsideOut_sample_dt,
                "log_TPMs",aggreagation_stat[2],
                end_region_boundaries_bin,
                segment_overlaps_ls,FALSE)
genesRndm_inTADs_InsideOut_sample_anchors_bins <-
  genesRndm_inTADs_InsideOut_sample_anchors_ls$bins

genesRndm_inTADs_InsideOut_sample_FLIPPED_bins_anchors_ls <- 
  binmaker_pb_f(
    genesRndm_inTADs_InsideOut_sample_anchors_FLIPPED_dt,
    "score",aggreagation_stat[2],
    end_region_boundaries_bin,
    segment_overlaps_ls,TRUE)
genesRndm_inTADs_InsideOut_sample_FLIPPED_bins_anchors_dt <- 
  genesRndm_inTADs_InsideOut_sample_FLIPPED_bins_anchors_ls$bins
#
#-------------------- List of Bins and statistics -------------
#
statistics_ls <- list()
statistics_ls[["Genes"]] <- 
  data.table(Plot=c("Domains","Anchors"),
              Color=c("lfc","lfc"),
              p_value=c(format(DEGs_fisher$summary$median_p,
                              scientific = TRUE,digits=2),
                        format(DEGs_fisher_anchors$summary$median_p,
                              scientific = TRUE,digits=2)),
              OR=c(round(DEGs_fisher$summary$median_OR,3),
                  round(DEGs_fisher_anchors$summary$median_OR,3)),
              rGenes=c(n_Rndm_degs, n_Rndm_degs_anchors),
              rTADs=c(n_Rndm_tads, n_Rndm_tads_anchors),
              DEGs=c(n_degs_peaks, n_degs_peaks_anchors),
              dTADs=c(n_tads_peaks, n_tads_peaks_anchors))

marks_inTADs_with_DEGs_InsideOut_bins_ls <- list()
marks_in_TADs_w_rndGenes_sample_InsideOut_bins_ls <- list()
marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED_ls <- list()
marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED_ls <- list()
marks_inTADs_with_DEGs_InsideOut_anchors_ls_ls <- list()
marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls_ls <- list()
marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls_ls <- list()
marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_FLIPPED_ls_ls <- list()

marks_inTADs_with_DEGs_InsideOut_bins_ls[["Genes"]] <- 
  DEGs_inTADs_bins
marks_in_TADs_w_rndGenes_sample_InsideOut_bins_ls[["Genes"]] <- 
  rndmGenes_inTADs_bins
marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED_ls[["Genes"]] <- 
  DEGs_inTADs_InsideOut_FLIPPED_bins
marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED_ls[["Genes"]] <- 
  genesRndm_inTADs_InsideOut_sample_FLIPPED_bins
marks_inTADs_with_DEGs_InsideOut_anchors_ls_ls[["Genes"]] <- 
  DEGs_inTADs_InsideOut_anchors_bins
marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls_ls[["Genes"]] <- 
  genesRndm_inTADs_InsideOut_sample_anchors_bins
marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls_ls[["Genes"]] <- 
  DEGs_inTADs_InsideOut_FLIPPED_bins_anchors_dt
marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_FLIPPED_ls_ls[["Genes"]] <- 
  genesRndm_inTADs_InsideOut_sample_FLIPPED_bins_anchors_dt
#
#----------------------- Plot density for DEGs ----------------
#
# Domain plots
print(paste0("Ploting ", j, " gene distributions over domains and anchors"))

plots_density_DEGs_overTADs <- 
  density_f(marks_inTADs_with_DEGs_InsideOut_bins_ls[["Genes"]],
            density_plot_y, aggreagation_stat[1],
            peaks_density_limits_ls$DEGs[1],peaks_density_limits_ls$DEGs[2],
            "black",peaks_density_limits_ls$DEGs[3],
            peaks_density_limits_ls$DEGs[4],"lfc",
            paste0("p=",statistics_ls[["Genes"]][Plot=="Domains"]$p_value,
                    "\n OR=",statistics_ls[["Genes"]][Plot=="Domains"]$OR),
            4,axis_tads_ls)
plots_density_sameRndm_Genes_overTADs <- 
  density_f(marks_in_TADs_w_rndGenes_sample_InsideOut_bins_ls[["Genes"]],
            density_plot_y, aggreagation_stat[1],
            peaks_density_limits_Rnd_ls$DEGs[1],peaks_density_limits_Rnd_ls$DEGs[2],
            "black",peaks_density_limits_Rnd_ls$DEGs[3],
            peaks_density_limits_Rnd_ls$DEGs[4],"lfc",
            paste0("Rnd_Genes:",statistics_ls[["Genes"]][Plot=="Domains"]$rGenes,
                    "\n Rnd_TADs:",statistics_ls[["Genes"]][Plot=="Domains"]$rTADs,
                    "\n Ref_DEGs:",statistics_ls[["Genes"]][Plot=="Domains"]$DEGs,
                    "\n Ref_TADs:",statistics_ls[["Genes"]][Plot=="Domains"]$dTADs),
            4,axis_tads_ls)

# Domain flipped plots
DEGs_overlaps_TADs_flipped_dt <- 
  density_f(marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED_ls[["Genes"]],
            density_plot_y, aggreagation_stat[1],
            peaks_density_limits_ls$DEGs[1],peaks_density_limits_ls$DEGs[2],
            "black",peaks_density_limits_ls$DEGs[3],
            peaks_density_limits_ls$DEGs[4],"lfc",
            paste0("p=",statistics_ls[["Genes"]][Plot=="Domains"]$p_value,
                    "\n OR=",statistics_ls[["Genes"]][Plot=="Domains"]$OR),
            4,axis_tads_ls)
genes_all_Rndm_overlaps_TADs_sample_flipped_dt <- 
  density_f(marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED_ls[["Genes"]],
            density_plot_y, aggreagation_stat[1],
            peaks_density_limits_Rnd_ls$DEGs[1],
            peaks_density_limits_Rnd_ls$DEGs[2],
            "black",peaks_density_limits_Rnd_ls$DEGs[3],
            peaks_density_limits_Rnd_ls$DEGs[4],"lfc",
            paste0("Rnd_Genes:",statistics_ls[["Genes"]][Plot=="Domains"]$rGenes,
                    "\n Rnd_TADs:",statistics_ls[["Genes"]][Plot=="Domains"]$rTADs,
                    "\n Ref_DEGs:",statistics_ls[["Genes"]][Plot=="Domains"]$DEGs,
                    "\n Ref_TADs:",statistics_ls[["Genes"]][Plot=="Domains"]$dTADs),
            4,axis_tads_ls)
# Anchors plots
plots_density_DEGs_overTADs_inAnchors <- 
  density_f(marks_inTADs_with_DEGs_InsideOut_anchors_ls_ls[["Genes"]],
            density_plot_y, aggreagation_stat[1],
            peaks_density_limits_ls$DEGs[1],
            peaks_density_limits_ls$DEGs[2],
            "black", peaks_density_limits_ls$DEGs[3],
            peaks_density_limits_ls$DEGs[4], "lfc",
            paste0("p=",statistics_ls[["Genes"]][Plot=="Anchors"]$p_value,
                    "\n OR=",statistics_ls[["Genes"]][Plot=="Anchors"]$OR),
            4,axis_anchors_ls)
plots_density_sameRndm_Genes_overTADs_inAnchors <- 
  density_f(marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls_ls[["Genes"]],
            density_plot_y,aggreagation_stat[1],
            peaks_density_limits_ls$DEGs[1],
            peaks_density_limits_ls$DEGs[2],
            "black",peaks_density_limits_ls$DEGs[3],
            peaks_density_limits_ls$DEGs[4],"lfc",
            paste0("Rnd_Genes:",statistics_ls[["Genes"]][Plot=="Anchors"]$rGenes,
                    "\n Rnd_TADs:",statistics_ls[["Genes"]][Plot=="Anchors"]$rTADs,
                    "\n Ref_DEGs:",statistics_ls[["Genes"]][Plot=="Anchors"]$DEGs,
                    "\n Ref_TADs:",statistics_ls[["Genes"]][Plot=="Anchors"]$dTADs),4,
            axis_anchors_ls)
# Anchors flipped plots
plots_density_DEGs_overTADs_inAnchors_flipped <- 
  density_f(marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls_ls[["Genes"]],
            density_plot_y, aggreagation_stat[1],
            peaks_density_limits_ls$DEGs[1],
            peaks_density_limits_ls$DEGs[2],
            "black", peaks_density_limits_ls$DEGs[3],
            peaks_density_limits_ls$DEGs[4], "lfc",
            paste0("p=",statistics_ls[["Genes"]][Plot=="Anchors"]$p_value,
                    "\n OR=",statistics_ls[["Genes"]][Plot=="Anchors"]$OR),
            4,axis_anchors_ls)
plots_density_sameRndm_Genes_overTADs_inAnchors_flipped <- 
  density_f(marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_FLIPPED_ls_ls[["Genes"]],
            density_plot_y,aggreagation_stat[1],
            peaks_density_limits_ls$DEGs[1],
            peaks_density_limits_ls$DEGs[2],
            "black",peaks_density_limits_ls$DEGs[3],
            peaks_density_limits_ls$DEGs[4],"lfc",
            paste0("Rnd_Genes:",statistics_ls[["Genes"]][Plot=="Anchors"]$rGenes,
                    "\n Rnd_TADs:",statistics_ls[["Genes"]][Plot=="Anchors"]$rTADs,
                    "\n Ref_DEGs:",statistics_ls[["Genes"]][Plot=="Anchors"]$DEGs,
                    "\n Ref_TADs:",statistics_ls[["Genes"]][Plot=="Anchors"]$dTADs),
            4, axis_anchors_ls)
#
#------------------------ ITERATION OVER MARKS ------------------
#
plots_density_peaks_overTADs_ls <- list()
plots_density_Rndm_peaks_overTADs_ls <- list()
plots_density_peaks_overTADs_flipped_ls <- list()
plots_density_Rndm_peaks_overTADs_flipped_ls <- list()
plots_density_Marks_inTADs_with_DEGs_atAnchors_ls <- list()
plots_density_Marks_inTADs_with_rndGenes_atAnchors_ls <- list()
plots_density_Marks_inTADs_with_DEGs_atAnchors_flipped_ls <- list()
plots_density_Marks_inTADs_with_rndGenes_atAnchors_flipped_ls <- list()

odds_ratio_Marks_ls <- list()
odds_ratio_Marks_anchors_ls <- list()
p_value_Marks_ls <- list()
p_value_Marks_anchors_ls <- list()

# i <- "H3K27ac_Enh"
# i <- "H3K27me3"
# i <- "H3K4me3"
side_peaks_ls <- names(peaks_list)[!names(peaks_list)%in%c("DEGs","genes_rndm_All")]
# side_peaks_ls <- "H3K27me3_Super"
for (i in side_peaks_ls){
  #
  #----------------------- Copy Marks to Plot ----------------
  #
  print(paste0("Ploting ",i," peaks over domains and anchors of ", j, " genes"))
  
  TADs_with_DEGs <- unique(peaks_overlaps_TADs_ls[["DEGs"]]$region)
  side_peaks_inTADs_with_DEGs_dt <- 
    peaks_overlaps_TADs_ls[[i]][region%in%TADs_with_DEGs]
  marks_inTADs_with_DEGs_InsideOut_dt <- 
    identify_InsideOut_f(side_peaks_inTADs_with_DEGs_dt)
  
  setkey(marks_inTADs_with_DEGs_InsideOut_dt, region)
  setkey(DEGs_inTADs_InsideOut_flip_vectors, region)
  marks_inTADs_with_DEGs_InsideOut_dt[
    DEGs_inTADs_InsideOut_flip_vectors[
      ,.(region,flip)],
    flip := i.flip,
    on = "region"]
  #
  #------------------------- Copy Marks to Export ---------------
  #                         AND STATISTICAL TEST
  #Filtering peaks in TADs, to TADs having active genes and marks
  
  setkey(peaks_overlaps_TADs_ls[[i]],region)
  setkey(genesRndm_inTADs_InsideOut_sample_dt,peak)
  rnd_genes <- unique(unlist(DEGs_fisher$results$sampled_ids))
  TADs_w_rnd_Genes <- 
    genesRndm_inTADs_InsideOut_sample_dt[rnd_genes,unique(region)]
  marks_in_TADs_w_rndGenes_sample_dt <- 
    peaks_overlaps_TADs_ls[[i]][
      ,.(region,peak,score,
          chrMark,startMark,endMark,
          upBoundary,downBoundary,
          regStart,regEnd)][TADs_w_rnd_Genes,nomatch = 0]
  
  #---------------------------------------------------------------------------
  rndm_regions <- unique(genesRndm_inTADs_InsideOut_dt$region)
  rndmGenes_inTADs_withMarks_dt <- 
    peaks_overlaps_TADs_ls[[i]][region%in%rndm_regions]
  
  rndmGenes_inTADs_withMarks_dt[,uniqueN(peak)]
  rndmGenes_inTADs_withMarks_dt[,uniqueN(region)]
  marks_in_TADs_w_rndGenes_sample_dt[,uniqueN(peak)]
  marks_in_TADs_w_rndGenes_sample_dt[,uniqueN(region)]
  side_peaks_inTADs_with_DEGs_dt[,uniqueN(peak)]
  side_peaks_inTADs_with_DEGs_dt[,uniqueN(region)]
  
  rndmGenes_inTADs_withMarks_InsideOut_dt <- 
    identify_InsideOut_f(rndmGenes_inTADs_withMarks_dt)
  
  marks_in_TADs_w_rndGenes_sample_InsideOut_dt <- 
    identify_InsideOut_f(marks_in_TADs_w_rndGenes_sample_dt)
  
  setkey(marks_in_TADs_w_rndGenes_sample_InsideOut_dt, region)
  setkey(genesRndm_inTADs_InsideOut_vectors, region)
  marks_in_TADs_w_rndGenes_sample_InsideOut_dt[
    genesRndm_inTADs_InsideOut_vectors[,.(region,flip)],
    flip := i.flip,
    on = "region"
  ]
  
  Marks_fisher <- run_fisher_fast(200,marks_inTADs_with_DEGs_InsideOut_dt,
                                  rndmGenes_inTADs_withMarks_InsideOut_dt,
                                  "inside_domain",
                                  "inTAD","outTAD","peak",ncores)
  Marks_fisher_anchors <- run_fisher_fast(200,marks_inTADs_with_DEGs_InsideOut_dt,
                                          rndmGenes_inTADs_withMarks_InsideOut_dt,
                                          "inside_boundaries",
                                          "inBoundaries","outBoundaries",
                                          "peak",ncores)
  
  odds_ratio_Marks_ls[[i]] <- 
    plot_fisher_f(
      Marks_fisher$results[,.(iter,p_value,odds_ratio)],
      "odds_ratio","Odds Ratio",i)
  p_value_Marks_ls[[i]] <- 
    plot_fisher_f(
      Marks_fisher$results[,.(iter,p_value,odds_ratio)],
      "p_value","p value",i)
  odds_ratio_Marks_anchors_ls[[i]] <- 
    plot_fisher_f(
      Marks_fisher_anchors$results[,.(iter,p_value,odds_ratio)],
      "odds_ratio","Odds Ratio",i)
  p_value_Marks_anchors_ls[[i]] <- 
    plot_fisher_f(
      Marks_fisher_anchors$results[,.(iter,p_value,odds_ratio)],
      "p_value","p value",i)
  #
  #----------------------- Lineal density of marks in TADs -------------------
  #
  (nTADs_with_Marks <- 
      uniqueN(side_peaks_inTADs_with_DEGs_dt$region))
  (nRndm_TADs_with_Marks <- 
      uniqueN(marks_in_TADs_w_rndGenes_sample_dt$region))
  (nPeaks_in_TADs <- 
      uniqueN(side_peaks_inTADs_with_DEGs_dt$peak))
  (nRndm_Peaks_in_TADs <- 
      uniqueN(marks_in_TADs_w_rndGenes_sample_dt$peak))
  
  marks_inTADs_with_DEGs_InsideOut_dt$log_Score <-
    log10(marks_inTADs_with_DEGs_InsideOut_dt$score+1)
  marks_inTADs_with_DEGs_InsideOut_bins <- 
    binmaker_redim_f(
      marks_inTADs_with_DEGs_InsideOut_dt,"log_Score",
      aggreagation_stat[2],nbins,FALSE)
  
  marks_in_TADs_w_rndGenes_sample_InsideOut_dt$log_Score <-
    log10(marks_in_TADs_w_rndGenes_sample_InsideOut_dt$score+1)
  marks_in_TADs_w_rndGenes_sample_InsideOut_bins <- 
    binmaker_redim_f(
      marks_in_TADs_w_rndGenes_sample_InsideOut_dt,"log_Score",
      aggreagation_stat[2],nbins,FALSE)
  #
  #---------------------- Flip TADs for Marks density plots ----------------
  #
  marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED <- 
    binmaker_redim_f(
      marks_inTADs_with_DEGs_InsideOut_dt,"log_Score",
      aggreagation_stat[2],nbins,TRUE)
  
  marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED <- 
    binmaker_redim_f(
      marks_in_TADs_w_rndGenes_sample_InsideOut_dt,"log_Score",
      aggreagation_stat[2],nbins,TRUE)
  #
  #---------------- Density of marks at anchors -----------
  #
  
  marks_inTADs_with_DEGs_InsideOut_anchors_ls <- 
    binmaker_pb_f(marks_inTADs_with_DEGs_InsideOut_dt,
                  "log_Score",aggreagation_stat[2],
                  end_region_boundaries_bin,
                  segment_overlaps_ls,FALSE)
  
  marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls <- 
    binmaker_pb_f(marks_in_TADs_w_rndGenes_sample_InsideOut_dt,
                  "log_Score",aggreagation_stat[2],
                  end_region_boundaries_bin,
                  segment_overlaps_ls,FALSE)
  #
  #---------------------- flip marks regions at anchors --------------------
  #
  marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls <- 
    binmaker_pb_f(marks_inTADs_with_DEGs_InsideOut_dt,
                  "log_Score",aggreagation_stat[2],
                  end_region_boundaries_bin,
                  segment_overlaps_ls,TRUE)
  
  marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_FLIPPED_ls <- 
    binmaker_pb_f(marks_in_TADs_w_rndGenes_sample_InsideOut_dt,
                  "log_Score",aggreagation_stat[2],
                  end_region_boundaries_bin,
                  segment_overlaps_ls,TRUE)
  #
  #---------------------- Bins files --------------------
  #
  marks_inTADs_with_DEGs_InsideOut_bins_ls[[i]] <-
    marks_inTADs_with_DEGs_InsideOut_bins
  marks_in_TADs_w_rndGenes_sample_InsideOut_bins_ls[[i]] <-
    marks_in_TADs_w_rndGenes_sample_InsideOut_bins
  marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED_ls[[i]] <-
    marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED
  marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED_ls[[i]] <-
    marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED
  marks_inTADs_with_DEGs_InsideOut_anchors_ls_ls[[i]] <-
    marks_inTADs_with_DEGs_InsideOut_anchors_ls$bins
  marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls_ls[[i]] <-
    marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls$bins
  marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls_ls[[i]] <-
    marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls$bins
  marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_FLIPPED_ls_ls[[i]] <-
    marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_FLIPPED_ls$bins
  #
  #---------------------- Density Plots --------------------
  #
  statistics_ls[[i]] <- 
    data.table(Plot=c("Domains","Anchors"),
                Color=c(colors_ls$hex[[i]],colors_ls$hex[[i]]),
                p_value=c(format(Marks_fisher$summary$median_p,
                                scientific = TRUE,digits=2),
                          format(Marks_fisher_anchors$summary$median_p,
                                scientific = TRUE,digits=2)),
                OR=c(round(Marks_fisher$summary$median_OR,3),
                    round(Marks_fisher_anchors$summary$median_OR,3)),
                rGenes=c(nRndm_Peaks_in_TADs,
                        marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls$peaks),
                rTADs=c(nRndm_TADs_with_Marks,
                        marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls$regions),
                DEGs=c(nPeaks_in_TADs,
                      marks_inTADs_with_DEGs_InsideOut_anchors_ls$peaks),
                dTADs=c(nTADs_with_Marks,
                        marks_inTADs_with_DEGs_InsideOut_anchors_ls$regions))
  # 1 Domains
  plots_density_peaks_overTADs_ls[[i]] <- 
    density_f(marks_inTADs_with_DEGs_InsideOut_bins_ls[[i]],
              density_plot_y, aggreagation_stat[1],
              peaks_density_limits_ls[[i]][1],
              peaks_density_limits_ls[[i]][2],
              "black",peaks_density_limits_ls[[i]][3],
              peaks_density_limits_ls[[i]][4],
              statistics_ls[[i]][Plot=="Domains"]$Color,
              paste0("p=",statistics_ls[[i]][Plot=="Domains"]$p_value,
                      "\n OR=",statistics_ls[[i]][Plot=="Domains"]$OR),
              4,axis_tads_ls)
  plots_density_Rndm_peaks_overTADs_ls[[i]] <- 
    density_f(marks_in_TADs_w_rndGenes_sample_InsideOut_bins_ls[[i]],
              density_plot_y,
              aggreagation_stat[1],
              peaks_density_limits_ls[[i]][1],
              peaks_density_limits_ls[[i]][2],
              "black",peaks_density_limits_ls[[i]][3],
              peaks_density_limits_ls[[i]][4],
              statistics_ls[[i]][Plot=="Domains"]$Color,
              paste0("rGenes=",statistics_ls[[i]][Plot=="Domains"]$rGenes,
                      "\n rTADs=",statistics_ls[[i]][Plot=="Domains"]$rTADs,
                      "\n DEGs=",statistics_ls[[i]][Plot=="Domains"]$DEGs,
                      "\n dTADs=",statistics_ls[[i]][Plot=="Domains"]$dTADs),
              4,axis_tads_ls)
  
  # 2 Domains_Flipped
  plots_density_peaks_overTADs_flipped_ls[[i]] <- 
    density_f(marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED_ls[[i]],
              density_plot_y, aggreagation_stat[1],
              peaks_density_limits_ls[[i]][1],
              peaks_density_limits_ls[[i]][2],
              "black",peaks_density_limits_ls[[i]][3],
              peaks_density_limits_ls[[i]][4],
              statistics_ls[[i]][Plot=="Domains"]$Color,
              paste0("p=",statistics_ls[[i]][Plot=="Domains"]$p_value,
                      "\n OR=",statistics_ls[[i]][Plot=="Domains"]$OR),
              4,axis_tads_ls)
  plots_density_Rndm_peaks_overTADs_flipped_ls[[i]] <-
    density_f(marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED_ls[[i]],
              density_plot_y,
              aggreagation_stat[1],
              peaks_density_limits_ls[[i]][1],
              peaks_density_limits_ls[[i]][2],
              "black",peaks_density_limits_ls[[i]][3],
              peaks_density_limits_ls[[i]][4],
              statistics_ls[[i]][Plot=="Domains"]$Color,
              paste0("rGenes=",statistics_ls[[i]][Plot=="Domains"]$rGenes,
                      "\n rTADs=",statistics_ls[[i]][Plot=="Domains"]$rTADs,
                      "\n DEGs=",statistics_ls[[i]][Plot=="Domains"]$DEGs,
                      "\n dTADs=",statistics_ls[[i]][Plot=="Domains"]$dTADs),
              4,axis_tads_ls)
  
  # Anchors
  plots_density_Marks_inTADs_with_DEGs_atAnchors_ls[[i]] <- 
    density_f(marks_inTADs_with_DEGs_InsideOut_anchors_ls_ls[[i]],
              density_plot_y, aggreagation_stat[1],
              peaks_density_limits_ls[[i]][1],
              peaks_density_limits_ls[[i]][2],
              "black",peaks_density_limits_ls[[i]][3],
              peaks_density_limits_ls[[i]][4],
              statistics_ls[[i]][Plot=="Anchors"]$Color,
              paste0("p=",statistics_ls[[i]][Plot=="Anchors"]$p_value,
                      "\n OR=",statistics_ls[[i]][Plot=="Anchors"]$OR),
              4,axis_anchors_ls)
  plots_density_Marks_inTADs_with_rndGenes_atAnchors_ls[[i]] <- 
    density_f(marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls_ls[[i]],
              density_plot_y, aggreagation_stat[1],
              peaks_density_limits_ls[[i]][1],
              peaks_density_limits_ls[[i]][2],
              "black",peaks_density_limits_ls[[i]][3],
              peaks_density_limits_ls[[i]][4],
              statistics_ls[[i]][Plot=="Anchors"]$Color,
              paste0("rGenes=",statistics_ls[[i]][Plot=="Anchors"]$rGenes,
                      "\n rTADs=",statistics_ls[[i]][Plot=="Anchors"]$rTADs,
                      "\n DEGs=",statistics_ls[[i]][Plot=="Anchors"]$DEGs,
                      "\n dTADs=",statistics_ls[[i]][Plot=="Anchors"]$dTADs),
              4,axis_anchors_ls)
  # Anchors Flipped
  plots_density_Marks_inTADs_with_DEGs_atAnchors_flipped_ls[[i]] <- 
    density_f(marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls_ls[[i]],
              density_plot_y, aggreagation_stat[1],
              peaks_density_limits_ls[[i]][1],
              peaks_density_limits_ls[[i]][2],
              "black",peaks_density_limits_ls[[i]][3],
              peaks_density_limits_ls[[i]][4],
              statistics_ls[[i]][Plot=="Anchors"]$Color,
              paste0("p=",statistics_ls[[i]][Plot=="Anchors"]$p_value,
                      "\n OR=",statistics_ls[[i]][Plot=="Anchors"]$OR),
              4,axis_anchors_ls)
  plots_density_Marks_inTADs_with_rndGenes_atAnchors_flipped_ls[[i]] <- 
    density_f(marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_FLIPPED_ls_ls[[i]],
              density_plot_y, aggreagation_stat[1],
              peaks_density_limits_ls[[i]][1],
              peaks_density_limits_ls[[i]][2],
              "black",peaks_density_limits_ls[[i]][3],
              peaks_density_limits_ls[[i]][4],
              statistics_ls[[i]][Plot=="Anchors"]$Color,
              paste0("rGenes=",statistics_ls[[i]][Plot=="Anchors"]$rGenes,
                      "\n rTADs=",statistics_ls[[i]][Plot=="Anchors"]$rTADs,
                      "\n DEGs=",statistics_ls[[i]][Plot=="Anchors"]$DEGs,
                      "\n dTADs=",statistics_ls[[i]][Plot=="Anchors"]$dTADs),
              4,axis_anchors_ls)
  
}

saveRDS(marks_inTADs_with_DEGs_InsideOut_bins_ls,
        file = paste0(dir_path_cluster,"/marks_inTADs_with_DEGs_InsideOut_bins_ls"))
saveRDS(marks_in_TADs_w_rndGenes_sample_InsideOut_bins_ls,
        file = paste0(dir_path_cluster,"/marks_in_TADs_w_rndGenes_sample_InsideOut_bins_ls"))
saveRDS(marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED_ls,
        file = paste0(dir_path_cluster,"/marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED_ls"))
saveRDS(marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED_ls,
        file = paste0(dir_path_cluster,"/marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED_ls"))
saveRDS(marks_inTADs_with_DEGs_InsideOut_anchors_ls_ls,
        file = paste0(dir_path_cluster,"/marks_inTADs_with_DEGs_InsideOut_anchors_ls_ls"))
saveRDS(marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls_ls,
        file = paste0(dir_path_cluster,"/marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls_ls"))
saveRDS(marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls_ls,
        file = paste0(dir_path_cluster,"/marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls_ls"))
saveRDS(marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_FLIPPED_ls_ls,
        file = paste0(dir_path_cluster,"/marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_FLIPPED_ls_ls"))
saveRDS(statistics_ls,
        file = paste0(dir_path_cluster,"/stats_ls"))
saveRDS(odds_ratio_Marks_ls,
        file = paste0(dir_path_cluster,"/odds_ratio_Marks_ls"))
saveRDS(p_value_Marks_ls,
        file = paste0(dir_path_cluster,"/p_value_Marks_ls"))
saveRDS(odds_ratio_Marks_anchors_ls,
        file = paste0(dir_path_cluster,"/odds_ratio_Marks_anchors_ls"))
saveRDS(p_value_Marks_anchors_ls,
        file = paste0(dir_path_cluster,"/p_value_Marks_anchors_ls"))
}
#
#---------------------------------STRIP_PLOTS-----------------------------
#
print(paste0("Ploting strips plots for ", j))
plots_density_TADs_and_peaks_ls <- 
append(plots_density_peaks_overTADs_ls,
        list(DEGs=plots_density_DEGs_overTADs),
        after = 0)
plots_density_Rndm_TADs_and_peaks_ls <- 
append(plots_density_Rndm_peaks_overTADs_ls,
        list(DEGs=plots_density_sameRndm_Genes_overTADs),
        after = 0)

plots_density_TADs_and_peaks_flipped_ls <-
append(plots_density_peaks_overTADs_flipped_ls,
        list(DEGs=DEGs_overlaps_TADs_flipped_dt),
        after = 0)
plots_density_Rndm_TADs_and_peaks_flipped_ls <- 
append(plots_density_Rndm_peaks_overTADs_flipped_ls,
        list(DEGs=genes_all_Rndm_overlaps_TADs_sample_flipped_dt),
        after = 0)

plots_density_TADs_and_peaks_atAnchors_ls <- 
append(plots_density_Marks_inTADs_with_DEGs_atAnchors_ls,
        list(DEGs=plots_density_DEGs_overTADs_inAnchors),
        after = 0)
plots_density_Rndm_TADs_and_peaks_atAnchors_ls <- 
append(plots_density_Marks_inTADs_with_rndGenes_atAnchors_ls,
        list(DEGs=plots_density_sameRndm_Genes_overTADs_inAnchors),
        after = 0)

plots_density_TADs_and_peaks_atAnchors_flipped_ls <-
append(plots_density_Marks_inTADs_with_DEGs_atAnchors_flipped_ls,
        list(DEGs=plots_density_DEGs_overTADs_inAnchors_flipped),
        after = 0)
plots_density_Rndm_TADs_and_peaks_atAnchors_flipped_ls <- 
append(plots_density_Marks_inTADs_with_rndGenes_atAnchors_flipped_ls,
        list(DEGs=plots_density_sameRndm_Genes_overTADs_inAnchors_flipped),
        after = 0)

plots_Odds_ratio_ls <- append(odds_ratio_Marks_ls,
                            list(DEGs=odds_ratio_DEGs),
                            after = 0)
plots_p_value_ls <- append(p_value_Marks_ls,
                          list(DEGs=p_value_DEGs),
                          after = 0)
print(paste0("Making list of strips plots for ", j))
marks_to_plot <- 
list(
  "density_strip"=list(side_peaks_ls,plots_density_TADs_and_peaks_ls,"no_X"),
  "density_strip_SELECTION"=list(strip_plot_selection,plots_density_TADs_and_peaks_ls,"no_X"),
  "density_strip_RNDM_SELECTION"=list(strip_plot_selection,plots_density_Rndm_TADs_and_peaks_ls,"no_X"),
  "density_strip_flipped"=list(side_peaks_ls,plots_density_TADs_and_peaks_flipped_ls,"no_X"),
  "density_strip_flipped_SELECTION"=list(strip_plot_selection,plots_density_TADs_and_peaks_flipped_ls,"no_X"),
  "density_strip_RNDM_flipped_SELECTION"=list(strip_plot_selection,plots_density_Rndm_TADs_and_peaks_flipped_ls,"no_X"),
  "density_strip_AtAnchors"=list(side_peaks_ls,plots_density_TADs_and_peaks_atAnchors_ls,"no_X"),
  "density_strip_AtAnchors_SELECTION"=list(strip_plot_selection,plots_density_TADs_and_peaks_atAnchors_ls,"no_X"),
  "density_strip_RNDM_AtAnchors_SELECTION"=list(strip_plot_selection,plots_density_Rndm_TADs_and_peaks_atAnchors_ls,"no_X"),
  "density_strip_AtAnchors_Flipped"=list(side_peaks_ls,plots_density_TADs_and_peaks_atAnchors_flipped_ls,"no_X"),
  "density_strip_AtAnchors_Flipped_SELECTION"=list(strip_plot_selection,plots_density_TADs_and_peaks_atAnchors_flipped_ls,"no_X"),
  "density_strip_RNDM_AtAnchors_Flipped_SELECTION"=list(strip_plot_selection,plots_density_Rndm_TADs_and_peaks_atAnchors_flipped_ls,"no_X"),
  "odds_ratio_strip"=list(strip_plot_selection,plots_Odds_ratio_ls,"X"),
  "p_value_strip"=list(strip_plot_selection,plots_p_value_ls,"X")
)
print(paste0("Ploting list of strip plots for ", j))
strip_f(marks_to_plot)
print(paste0(j, " genes successfully plotted"))