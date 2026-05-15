#boundariesNpeaks_meta.v14.r
#------------------------------Libraries----------------------------------------
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(parallel)
library(scales)
library(patchwork)
#--------------------------------Env Configuration------------------------------
start.time <- Sys.time()
options(scipen=999)
dir_path <- getwd()
setwd(dir_path)

#---------------------------------Configuration file ---------------------------
variables <- commandArgs(trailingOnly = TRUE)
configuracion <- read.table(as.character(variables[1]),header = T,row.names = 1, sep=" ")
#configuracion <- read.table("config/marks_loops_overlaps.density.byDomain.config", header = T,row.names = 1)
#configuracion <- read.table("config/marks_TADs_Down_p05_lfc0.density_heatmaps", header = T,row.names = 1)
#--------------------------------- Functions------------------------------------
extract_colors_f <- function(list_peaks){
  colors_rgb_ls <- list()
  colors_hex_ls <- list()
  for (i in names(list_peaks)){
    
    marks <- copy(list_peaks[[i]])
    rgb_vals <- tstrsplit(marks$color[1], ",", fixed = TRUE)
    color_hex <- rgb(as.integer(rgb_vals[[1]]),
                     as.integer(rgb_vals[[2]]),
                     as.integer(rgb_vals[[3]]),
                     maxColorValue = 255)
    colors_rgb_ls[[i]] <- marks$color[1]
    colors_hex_ls[[i]] <- color_hex
  }
  return(list(rgb=colors_rgb_ls,hex=colors_hex_ls))
}

overlaps_f <- function(list_peaks,regions){
  regions$chr <- gsub("chr","",regions$chr)
  regions_gr <- GRanges(regions)
  overlaps_list <- list()
  for (i in names(list_peaks)){
    marks <- copy(list_peaks[[i]])
    marks$chr <- gsub("chr","",marks$chr)
    marks_gr <- GRanges(marks)
    #---------------------------------------------------------------------------
    hits <- findOverlaps(marks_gr, regions_gr)
    #---------------------------------------------------------------------------
    hits_dt <- data.table(
      chrMark=marks$chr[hits@from],
      startMark = start(marks_gr)[hits@from],
      endMark = end(marks_gr)[hits@from],
      peak = marks$peak[hits@from],
      score = marks$score[hits@from],
      region = regions$id[hits@to],
      regStart = start(regions_gr)[hits@to],
      regEnd = end(regions_gr)[hits@to],
      upBoundary = regions$upBoundary[hits@to],
      downBoundary = regions$downBoundary[hits@to]
    )
    #-------------------------------------------------------------------------------
    overlaps_list[[i]] <- hits_dt
  }
  return(overlaps_list)
}

identify_InsideOut_f <- function(pb_dt){
  pb_dt <- pb_dt[, .(region,peak,score, 
                     chrMark,startMark,endMark,
                     regStart,regEnd,
                     upBoundary,downBoundary)]
  pb_dt[, `:=`(
    upBoundary_start   = upBoundary - boundaries_pb_out,
    upBoundary_end     = upBoundary + boundaries_pb_in,
    downBoundary_start = downBoundary - boundaries_pb_in,
    downBoundary_end   = downBoundary + boundaries_pb_out
  )]
  
  pb_dt[,length_domain := pmax(0, downBoundary_start - upBoundary_end)]
  
  pb_dt[, `:=`(
    upTAD_out   = upBoundary_start - (length_domain/2),
    downTAD_out = downBoundary_end + (length_domain/2)
  )]
  
  #Filter_out_peaks
  pb_dt <- pb_dt[endMark >= upTAD_out & startMark <= downTAD_out]
  
  #Inside TADs
  pb_dt[, `:=`(
    inTAD = startMark > upBoundary_end & endMark < downBoundary_start,
    outTAD =
      (startMark > upTAD_out & endMark < upBoundary_start) |
      (startMark > downBoundary_end & endMark < downTAD_out),
    inside_UpBoundary =
      startMark <= upBoundary_end & endMark >= upBoundary_start,
    inside_DownBoundary =
      startMark <= downBoundary_end & endMark >= downBoundary_start,
    outside_Up_left =
      startMark <= upBoundary_start &
      endMark >= (upBoundary_start - boundary_width_pb),
    outside_Down_right =
      startMark <= (downBoundary_end + boundary_width_pb) &
      endMark >= downBoundary_end,
    crossing_Up_right =
      startMark <= (upBoundary_end + boundary_width_pb) &
      endMark >= upBoundary_end,
    crossing_Down_left =
      startMark <= downBoundary_start &
      endMark >= (downBoundary_start - boundary_width_pb),
    farOutside_Up_left =
      startMark <= (upBoundary_start - boundary_width_pb) &
      endMark >= (upBoundary_start - (2*boundary_width_pb)),
    farCrossing_Up_right =
      startMark <= (upBoundary_end + (2*boundary_width_pb)) &
      endMark >= (upBoundary_end + boundary_width_pb),
    farCrossing_Down_left =
      startMark <= (downBoundary_start - boundary_width_pb) &
      endMark >= (downBoundary_start - (2*boundary_width_pb)),
    farOutside_Down_right =
      startMark <= (downBoundary_end + (2*boundary_width_pb)) &
      endMark >= (downBoundary_end + boundary_width_pb)
  )]
  #Inside domains
  pb_dt[, inside_domain := fcase(
    inTAD, "inTAD",
    outTAD, "outTAD",
    default = "Faraway"
  )]
  #Inside Boundaries
  pb_dt[, inside_boundaries := fcase(
    inside_UpBoundary | inside_DownBoundary, "inBoundaries",
    outside_Up_left | outside_Down_right, "outBoundaries",
    default = "Faraway"
  )]
  #Near Boundaries
  pb_dt[, near_boundaries_outTADs :=
          (outside_Up_left | outside_Down_right) &
          !(inside_UpBoundary | inside_DownBoundary)]
  pb_dt[, near_boundaries_inTADs :=
          (crossing_Up_right | crossing_Down_left) &
          !(inside_UpBoundary | inside_DownBoundary)]
  pb_dt[, near_boundaries := fcase(
    (near_boundaries_outTADs | near_boundaries_inTADs), "nearBoundaries",
    (farOutside_Up_left | farOutside_Down_right |
       farCrossing_Up_right | farCrossing_Down_left), "farBoundaries",
    default = "Faraway"
  )]
  #Location
  pb_dt[, location := fcase(
    inside_boundaries == "inBoundaries", "inBoundaries",
    inTAD, "inTAD",
    default = "outside"
  )]
  pb_dt[,`:=`(inTAD=NULL,
              outTAD=NULL)]
  return(pb_dt)
}

run_fisher_N_resamplings_genes <- 
  function(n_iter = 100,DEGs_dt,random_source_dt,boolean_column,var,nvar,sampling_id) {
    results <- vector("list", n_iter)
    rndG <- data.table()
    for (i in seq_len(n_iter)) {
      setkeyv(random_source_dt,sampling_id)
      setkeyv(DEGs_dt,sampling_id)
      # 1. Sample random peaks (same number as observed)
      sampled_rndm <- random_source_dt[
        sample(unique(get(as.character(sampling_id))), 
               size = uniqueN(DEGs_dt[[as.character(sampling_id)]]),
               replace = FALSE)]
      # 2. Compute counts for observed and random data
      obsPeaks_in  <- DEGs_dt[get(boolean_column) == var,  uniqueN(peak)]
      obsPeaks_out <- DEGs_dt[get(boolean_column) == nvar, uniqueN(peak)]
      rndPeaks_in  <- sampled_rndm[get(boolean_column) == var,  uniqueN(peak)]
      rndPeaks_out <- sampled_rndm[get(boolean_column) == nvar, uniqueN(peak)]
      
      obs_inTADs  <- DEGs_dt[get(boolean_column) == var,  uniqueN(region)]
      obs_outTADs <- DEGs_dt[get(boolean_column) == nvar, uniqueN(region)]
      rnd_inTADs  <- sampled_rndm[get(boolean_column) == var,  uniqueN(region)]
      rnd_outTADs <- sampled_rndm[get(boolean_column) == nvar, uniqueN(region)]
      
      # 3. Contingency table
      mat <- matrix(c(obsPeaks_in, obsPeaks_out, rndPeaks_in, rndPeaks_out),
                    nrow = 2, byrow = TRUE)
      # 4. Fisher test
      ft <- fisher.test(mat)
      # 5. Store results
      results[[i]] <- data.table(iter = i,p_value = ft$p.value,
                                 odds_ratio = ft$estimate,
                                 obsPeaks_inTADs=obsPeaks_in,obsPeaks_outTADs=obsPeaks_out,
                                 rndPeaks_inTADs=rndPeaks_in,rndPeaks_outTADs=rndPeaks_out,
                                 obs_inTADs=obs_inTADs,obs_outTADs=obs_outTADs,
                                 rnd_inTADs=rnd_inTADs,rnd_outTADs=rnd_outTADs)
      
      rndG <- rbind(rndG, data.table(p_value = ft$p.value,
                                     peaks = list(sampled_rndm[,unique(peak)])))
    }
    # Combine results
    results_dt <- rbindlist(results)
    # Summary
    summary_dt <- results_dt[, .(
      median_p = median(p_value),
      mean_p = mean(p_value),
      median_OR = median(odds_ratio),
      mean_OR = mean(odds_ratio)
    )]
    p_values <- rndG$p_value
    p_values_sorted <- sort(p_values)
    middle_p_value <- p_values_sorted[ceiling(length(p_values_sorted)/2)]
    rndm_DEGs_dt <- rndG[p_value==middle_p_value][sample(.N,1)]
    
    return(list(results = results_dt,summary = summary_dt,peaks_dt=rndG,rndm_DEGs_median_p=rndm_DEGs_dt))
  }
run_fisher_fast <- 
  function(n_iter = 100,DEGs_dt,
           random_source_dt,boolean_column,var,nvar, sampling_id,
           n_cores = parallel::detectCores() - 1) {
    
    # ---- Observed ----
    obsPeaks_in  <- DEGs_dt[get(boolean_column) == var,  uniqueN(peak)]
    obsPeaks_out <- DEGs_dt[get(boolean_column) == nvar, uniqueN(peak)]
    
    # ---- Precompute flags ----
    random_source_dt[, in_flag  := get(boolean_column) == var]
    random_source_dt[, out_flag := get(boolean_column) == nvar]
    
    setkeyv(random_source_dt, sampling_id)
    
    sampling_ids <- unique(random_source_dt[[sampling_id]])
    n_sample <- uniqueN(DEGs_dt[[sampling_id]])
    
    results_list <- mclapply(seq_len(n_iter), function(i) {
      
      sampled_ids <- sample(sampling_ids, n_sample)
      
      sampled_rndm <- random_source_dt[.(sampled_ids)]
      
      rndPeaks_in  <- uniqueN(sampled_rndm[in_flag == TRUE, peak])
      rndPeaks_out <- uniqueN(sampled_rndm[out_flag == TRUE, peak])
      
      mat <- matrix(c(obsPeaks_in, obsPeaks_out,
                      rndPeaks_in, rndPeaks_out),
                    nrow = 2, byrow = TRUE)
      
      ft <- fisher.test(mat)
      
      list(
        iter = i,
        p_value = ft$p.value,
        odds_ratio = ft$estimate,
        sampled_ids = list(sampled_ids) 
      )
      
    }, mc.cores = n_cores)
    
    results_dt <- rbindlist(results_list)
    
    summary_dt <- results_dt[, .(
      median_p = median(p_value),
      mean_p   = mean(p_value),
      median_OR = median(odds_ratio),
      mean_OR   = mean(odds_ratio)
    )]
    
    # ---- Recover peaks correctly ----
    middle_row <- results_dt[order(p_value)][ceiling(.N/2)]
    selected_ids <- middle_row$sampled_ids[[1]]
    
    selected_peaks <- random_source_dt[
      .(selected_ids), unique(peak)
    ]
    
    return(list(
      results = results_dt,
      summary = summary_dt,
      selected_peaks = selected_peaks
    ))
  }
plot_fisher_f <- function(fisher_results,variable,title,x_lable){
  plot_fisher <- 
    ggplot(fisher_results, aes(x = get(variable))) +
    geom_histogram(bins = 40) +
    geom_vline(xintercept = median(fisher_results[[variable]]),
               color = "red", linetype = "dashed") + plot_theme_statistics +
    geom_text(x = Inf, y = Inf, 
              label = paste0(title,"\n Median: ",
                             format(median(fisher_results[[variable]]),
                                    scientific = TRUE,digits=2)),
              color = "red", hjust = 1.1, vjust = 1.5,size = 4)+
    labs(x = title, y = paste0(x_lable,"\n Frequency"))
  return(plot_fisher)
}

binmaker_redim_f <- function(hits_dt,measure_fill,agg,maxbins,do_flip = FALSE){
  hits_dt <- hits_dt[,.(region,peak,var=get(measure_fill),
                        chrMark,startMark,endMark,
                        regStart,regEnd,
                        upBoundary,downBoundary,flip)]
  hits_dt[, peakStart_rel := (startMark - regStart) / (regEnd - regStart)]
  hits_dt[, peakEnd_rel := (endMark - regStart) / (regEnd - regStart)]
  #-------------------------------------------------------------------------------
  hits_dt[, start_bin := round(peakStart_rel * nbins,0) + 1]
  hits_dt[, end_bin := round(peakEnd_rel * nbins,0) + 1]
  hits_dt[start_bin > nbins, start_bin := nbins]
  hits_dt[end_bin > nbins, end_bin := nbins]
  hits_dt[start_bin < 0, start_bin := 1]
  hits_dt[end_bin < 0, end_bin := 1]
  hits_dt[start_bin == 0, start_bin := 1]
  hits_dt[end_bin == 0, end_bin := 1]
  #-------------------------------------------------------------------------------
  hits_dt <- hits_dt[,.(bin = seq(start_bin[1], end_bin[1])),
                     by = .(region,peak,var,flip)]
  if(do_flip){
    hits_dt <- hits_dt[flip == TRUE, bin := maxbins - bin + 1]
  }
  all_bins <- data.table(bin = 1:maxbins) 
  bin_stats <-
    hits_dt[,.(N = .N,agg_score = get(agg)(var, na.rm = TRUE)),by = bin]
  setkey(bin_stats, bin)
  bin_stats <- bin_stats[all_bins]
  bin_stats[is.na(N), N := 0]
  bin_stats[is.na(agg_score), agg_score := 0]
  bin_stats[, density := N / uniqueN(hits_dt$region)]
  
  return(bin_stats)
}

density_f <- function(bin_stats, measure, agg_label,
                      min_score, max_score, line_color,
                      ymin,ymax,scala,label_text,sizeL,axis_ls) {
  
  ymin <- if (is.na(ymin)) min(bin_stats[[measure]], na.rm = TRUE) else ymin
  ymax <- if (is.na(ymax)) max(bin_stats[[measure]], na.rm = TRUE) else ymax
  
  range_vals <- range(bin_stats[,agg_score], na.rm = TRUE)
  min_limit <- if (is.na(min_score)) range_vals[1] else min_score
  max_limit <- if (is.na(max_score)) range_vals[2] else max_score
  
  if (scala == "lfc"){
    gradient <- scale_fill_gradient2(low = "blue4",mid = "white",high = "red4",
                                     midpoint = 0,limits = c(min_limit,max_limit),
                                     oob = scales::squish,name=agg_label)
  } else {
    gradient <- scale_fill_gradient2(low = "white",
                                     high = scala,
                                     limits = c(min_limit,max_limit),
                                     oob = scales::squish,name=agg_label)
  }
  
  densityPlot <- ggplot(bin_stats, aes(x = bin)) +
    geom_ribbon(aes(ymin = ymin, ymax = get(measure),
                    fill = agg_score),alpha = alpha_level) +
    geom_line(aes(y = get(measure)),linewidth = 0.5, color=line_color) +
    gradient + plot_theme + axis_ls + 
    scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0))+
    annotate("text", x = Inf, y = Inf, label = label_text,
             hjust = 1, vjust = 1, size = sizeL)
  
  return(densityPlot)
}

binmaker_pb_f <- function(datatable_bins,
                          measure_fill,agg,maxbins,
                          filter_list,do_flip = FALSE){
  datatable_bins <- datatable_bins[
    rowSums(datatable_bins[, ..segment_overlaps_ls]) > 0]
  
  datatable_bins <- datatable_bins[,.(
    region, peak, chrMark, startMark, endMark,
    var=get(measure_fill), upBoundary_start, upBoundary_end,
    downBoundary_start, downBoundary_end,flip
  )]
  datatable_bins[, `:=`(
    up_win_start   = upBoundary_start - (upBoundary_Start_bin * resol),
    up_win_end     = upBoundary_end   + (upBoundary_Start_bin * resol),
    
    down_win_start = downBoundary_start - (upBoundary_Start_bin * resol),
    down_win_end   = downBoundary_end   + (upBoundary_Start_bin * resol)
  )]
  
  datatable_bins[, overlaps_up :=
                   endMark >= up_win_start & startMark <= up_win_end]
  datatable_bins[, overlaps_down :=
                   endMark >= down_win_start & startMark <= down_win_end]
  # ----------- Trim -----------
  # Up window
  datatable_bins[overlaps_up == TRUE, `:=`(
    up_start_trim = pmax(startMark, up_win_start),
    up_end_trim   = pmin(endMark,   up_win_end)
  )]
  # Down window
  datatable_bins[overlaps_down == TRUE, `:=`(
    down_start_trim = pmax(startMark, down_win_start),
    down_end_trim   = pmin(endMark,   down_win_end)
  )]
  
  # Compute bins
  datatable_bins[overlaps_up == TRUE, `:=`(
    bin_up_start = floor((up_start_trim - up_win_start) / resol) + 1,
    bin_up_end   = floor((up_end_trim   - up_win_start) / resol) + 1
  )]
  datatable_bins[overlaps_down == TRUE, `:=`(
    bin_down_start = floor((down_start_trim - down_win_start) / resol) + 1,
    bin_down_end   = floor((down_end_trim   - down_win_start) / resol) + 1
  )]
  # Shift down
  datatable_bins[overlaps_down == TRUE, `:=`(
    bin_down_start = bin_down_start + middle_region_boundaries_bin,
    bin_down_end   = bin_down_end   + middle_region_boundaries_bin
  )]
  
  # ---------------- FIX: clamp bins ----------------
  # Up bins 
  datatable_bins[overlaps_up == TRUE, `:=`(
    bin_up_start = pmax(1,  pmin(middle_region_boundaries_bin, bin_up_start)),
    bin_up_end   = pmax(1,  pmin(middle_region_boundaries_bin, bin_up_end))
  )]
  # Down bins
  datatable_bins[overlaps_down == TRUE, `:=`(
    bin_down_start = pmax(middle_region_boundaries_bin+1, pmin(end_region_boundaries_bin, bin_down_start)),
    bin_down_end   = pmax(middle_region_boundaries_bin+1, pmin(end_region_boundaries_bin, bin_down_end))
  )]
  
  dt_up <- if (datatable_bins[overlaps_up == TRUE, .N] > 0) {
    datatable_bins[
      overlaps_up == TRUE,
      .(bin = seq(bin_up_start, bin_up_end)),
      by = .(region, peak, var, flip)
    ]
  } else {
    data.table(region=character(), peak=character(), var=numeric(), flip=logical(), bin=integer())
  }
  dt_down <- if (datatable_bins[overlaps_down == TRUE, .N] > 0) {
    datatable_bins[
      overlaps_down == TRUE,
      .(bin = seq(bin_down_start, bin_down_end)),
      by = .(region, peak, var, flip)
    ]
  } else {
    data.table(region=character(), peak=character(), var=numeric(), flip=logical(), bin=integer())
  }
  datatable_bins_expanded <- rbind(dt_up, dt_down)
  #-----------------------------------------------------------------------------
  if(do_flip){
    datatable_bins_expanded[flip == TRUE, bin := maxbins - bin + 1]
  }
  #-----------------------------------------------------------------------------
  all_bins <- data.table(bin = 1:maxbins) 
  bin_stats <-
    datatable_bins_expanded[
      ,.(N = .N,agg_score = get(agg)(var, na.rm = TRUE)),by = bin]
  setkey(bin_stats, bin)
  bin_stats <- bin_stats[all_bins]
  bin_stats[is.na(N), N := 0]
  bin_stats[is.na(agg_score), agg_score := 0]
  n_regions <- uniqueN(datatable_bins_expanded$region)
  bin_stats[, density := if (n_regions > 0) N / n_regions else 0]
  #-----------------------------------------------------------------------------
  return(list(bins=bin_stats,
              peaks=uniqueN(datatable_bins_expanded$peak),
              regions=uniqueN(datatable_bins_expanded$region)))
}

flip_regions_f <- function(dt, mode=c("TAD","boundary")){
  
  mode <- match.arg(mode)
  
  if(mode=="TAD"){
    
    dt_counts <- dt[, {
      center <- (regStart[1] + regEnd[1]) / 2
      peak_center <- (startMark + endMark) / 2
      
      .(
        left  = sum(peak_center < center, na.rm=TRUE),
        right = sum(peak_center >= center, na.rm=TRUE)
      )
      
    }, by=region]
    
  } else if(mode=="boundary"){
    
    dt_counts <- dt[, .(
      left  = sum(inside_UpBoundary  | crossing_Up_right  |
                    outside_Up_left    | farOutside_Up_left |
                    farCrossing_Up_right, na.rm=TRUE),
      
      right = sum(inside_DownBoundary | crossing_Down_left |
                    outside_Down_right  | farOutside_Down_right |
                    farCrossing_Down_left, na.rm=TRUE)
      
    ), by=region]
  }
  
  dt_counts[, flip := right > left]
  
  return(dt_counts)
}

graficar_density <- function(plot,mark){
  nplots <- length(plot)
  if(formato=="png"){
    path.apaplot <- paste0(dir_path_cluster,"/",nomen,".",mark,".",nbins,"_bins.png")
    device.plot <- png(path.apaplot,width=density_plot_width,
                       height=5*nplots,units="cm",res=500)
  } else if (formato=="pdf"){
    path.apaplot <- paste0(dir_path_cluster,"/",nomen,".",mark,".",nbins,"_bins.pdf")
    device.plot <- pdf(path.apaplot, width = density_plot_width / 2.54, 
                       height = (5*nplots) / 2.54) 
  } else{
    stop("Image format not available. Use instead 'pdf' or 'png'")
  } 
  path.apaplot
  device.plot
  print(plot)
  dev.off()
}

strip_f <- function(lista){
  for (name in names(lista)){
    plots_density_sel <- lista[[name]][[2]][lista[[name]][[1]]]
    
    if(lista[[name]][[3]]=="no_X"){
      plots_no_x <- head(names(plots_density_sel),
                         length(names(plots_density_sel))-1)
      for (i in plots_no_x){
        plots_density_sel[[i]] <- plots_density_sel[[i]]+remove_x
      }
    }
    plots_density_sel_wraped <- wrap_plots(plots_density_sel, ncol = 1)
    graficar_density(plots_density_sel_wraped,name)
  }
}


#--------------------------------- Iterations ----------------------------------
# j <- "Monocytes"
just_plot <- strsplit(configuracion["just_plot:",1]," ")[[1]]
for (j in just_plot){
  nomen <- j
  #------------------------------Variables--------------------------------------
  cluster_sense <- strsplit(configuracion["cluster_sense:", 1]," ")[[1]]
  plot_select <- as.character(configuracion["plot_select:", 1])
  strip_plot_selection <- strsplit(configuracion["strip_plot_selection:", 1]," ")[[1]]
  density_plot_y <- as.character(configuracion["density_plot_y:", 1])
  density_plot_width <- as.numeric(configuracion["density_plot_width:", 1])
  sample <- as.character(configuracion["sample:", 1])
  tads_dir <- as.character(configuracion["tads_dir:", 1])
  tads_file <- as.character(configuracion["tads:", 1])
  genes_dir <- as.character(configuracion["genes_dir:", 1])
  genes_file <- as.character(configuracion["genes_file:", 1])
  genes_ctrl <- as.character(configuracion["genes_ctrl:", 1])
  #Defining how TAD regions are ordered
  aggreagation_stat <- strsplit(configuracion["aggreagation_stat:", 1]," ")[[1]]
  alpha_level <- as.numeric(configuracion["alpha_level:", 1])
  ncores <- as.numeric(configuracion["ncores:", 1])
  #Define dimensions
  nbins <- as.numeric(configuracion["nbins:", 1])
  boundaries_pb_in <- as.numeric(configuracion["boundaries_pb_in:", 1])
  boundaries_pb_out <- as.numeric(configuracion["boundaries_pb_out:", 1])
  cross_boundaries_pb <- as.numeric(configuracion["cross_boundaries_pb:", 1])
  outside_boundaries_pb <- as.numeric(configuracion["outside_boundaries_pb:", 1])
  resol <- as.numeric(configuracion["resol:", 1])
  TAD_extend <- as.numeric(configuracion["TAD_extend:", 1])
  distance_To_Nearest_TAD <- as.numeric(configuracion["distance_To_Nearest_TAD:", 1])
  formato <- as.character(configuracion["formato:", 1])
  row_height <- as.numeric(configuracion["row_height:", 1])
  base_height <- as.numeric(configuracion["base_height:", 1])
  peaks_density_limits <- 
    as.list(strsplit(configuracion["peaks_density_limits:", 1]," ")[[1]])
  peaks_density_limits_Rnd <- 
    as.list(strsplit(configuracion["peaks_density_limits_Rnd:", 1]," ")[[1]])
  peaks_cols <- strsplit(configuracion["peaks_cols:", 1]," ")[[1]]
  marks_names <- strsplit(configuracion["marks_names:", 1]," ")[[1]]
  marks_files <- as.list(strsplit(configuracion["marks_files_list:", 1]," ")[[1]])
  normalize_by <- as.character(configuracion["normalize_by:", 1])
  (runMode <- as.character(configuracion["runMode:",j]))
  (bins_paths <- as.character(configuracion["bins_paths:",j]))
  #
  #------------------------------Make_directories-------------------------------
  #
  dir.create(file.path(dir_path,"metaplots"))
  dir_path_sample <- file.path(dir_path,"metaplots",sample)
  dir.create(dir_path_sample)
  dir_path_cluster <- file.path(dir_path_sample,j)
  dir.create(dir_path_cluster)
  dir_path_cluster_beds <- file.path(dir_path_cluster,"beds")
  dir.create(dir_path_cluster_beds) 
  #
  #------------------------------Set pb based bins limits ----------------------
  #
  (boundary_width_pb <- boundaries_pb_in + boundaries_pb_out)
  (cross_boundaries_bins <- cross_boundaries_pb / resol)
  (outside_boundaries_bins <- outside_boundaries_pb / resol)
  (boundaries_bins_in <- boundaries_pb_in / resol)
  (boundaries_bins_out <- boundaries_pb_out / resol)
  (boundary_width_bin <- boundaries_bins_in+boundaries_bins_out)
  (upBoundary_Start_bin <- outside_boundaries_bins)
  (upBoundary_End_bin <- upBoundary_Start_bin+boundary_width_bin)
  (middle_region_boundaries_bin <- upBoundary_Start_bin+boundary_width_bin+cross_boundaries_bins)
  (downBoundary_Start_bin <- middle_region_boundaries_bin+cross_boundaries_bins)
  (downBoundary_End_bin <- downBoundary_Start_bin+boundary_width_bin)
  (end_region_boundaries_bin <- downBoundary_End_bin+outside_boundaries_bins)
  tad_center <- nbins / 2 
  ext_bins <- ceiling(nbins * TAD_extend / (1 + 2 * TAD_extend))
  mid_bins <- nbins - 2 * ext_bins
  #
  #---------------------------- Plot Axis and theme ------------------------
  #
  axis_tads_ls <- list(scale_x_continuous(limits = c(1, nbins),
                                          breaks=c(1,ext_bins,ext_bins+mid_bins,nbins),
                                          labels = c(paste0("-",TAD_extend,"X"),
                                                     paste0("Up"),
                                                     paste0("Down"),
                                                     paste0("+",TAD_extend,"X")),
                                          expand = c(0, 0)),
                       geom_vline(xintercept = ext_bins,linewidth = 0.7,linetype = "dashed",alpha=0.6),
                       geom_vline(xintercept = ext_bins+mid_bins,linewidth = 0.7,linetype = "dashed",alpha=0.6))
  
  axis_bined_boundaries_ls <- list(scale_x_continuous(limits = c(0, end_region_boundaries_bin-1),
                                                      breaks=c(0,
                                                               upBoundary_Start_bin,
                                                               upBoundary_End_bin,
                                                               downBoundary_Start_bin,
                                                               downBoundary_End_bin,
                                                               end_region_boundaries_bin-1),
                                                      labels = c(paste0("-",outside_boundaries_pb/1000,"kb"),
                                                                 paste0("Up"),
                                                                 paste0("Boundary"),
                                                                 paste0("Down"),
                                                                 paste0("Boundary"),
                                                                 paste0("+",outside_boundaries_pb/1000,"kb")),
                                                      expand = c(0, 0)),
                                   geom_vline(xintercept = c(upBoundary_Start_bin,
                                                             upBoundary_End_bin,
                                                             downBoundary_Start_bin,
                                                             downBoundary_End_bin),
                                              linewidth = 0.5, linetype = "dashed",alpha=0.6),
                                   geom_vline(xintercept = middle_region_boundaries_bin,linewidth = 0.5,
                                              linetype = "dotted",alpha=0.25))
  #
  plot_theme <- list(theme_minimal(),theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                           axis.text.x = element_text(size = 14,angle = 0, hjust = 0.75),
                                           axis.title.x = element_blank(),
                                           axis.text.y = element_text(size = 14,angle = 0, hjust = 0.5),
                                           axis.title.y = element_blank(),
                                           plot.title = element_text(hjust = 0.5,size=14),
                                           axis.line = element_blank(),
                                           axis.ticks.x = element_line(color = "black", linewidth = 0.5),
                                           axis.ticks.length.x = unit(3, "pt"),
                                           axis.ticks.y = element_line(color = "black", linewidth = 0.5),
                                           axis.ticks.length.y = unit(3, "pt"),
                                           panel.border = element_rect(linetype = "solid", fill = NA),
                                           plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
                                           legend.title = element_text(hjust = 0.5)))
  plot_theme_statistics <- list(theme_minimal(),theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                                      axis.text.x = element_text(size = 14,angle = 0, hjust = 0.75),
                                                      axis.title.x = element_text(size = 14,angle = 0, hjust = 0.5),
                                                      axis.text.y = element_text(size = 14,angle = 0, hjust = 0.5),
                                                      axis.title.y = element_text(size = 14,angle = 90, hjust = 0.5),
                                                      plot.title = element_text(hjust = 0.5,size=14),
                                                      axis.line = element_blank(),
                                                      axis.ticks.x = element_line(color = "black", linewidth = 0.5),
                                                      axis.ticks.length.x = unit(3, "pt"),
                                                      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
                                                      axis.ticks.length.y = unit(3, "pt"),
                                                      panel.border = element_rect(linetype = "solid", fill = NA),
                                                      plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
                                                      legend.title = element_text(hjust = 0.5)))
  remove_x <- theme(axis.title.x = element_blank(),axis.text.x  = element_blank(),
                    axis.ticks.x = element_blank())
  #
  #---------------------------Set plot Limits ----------------------------------
  #
  degs_marks_names <- c("DEGs",marks_names)
  names(peaks_density_limits) <- degs_marks_names
  names(peaks_density_limits_Rnd) <- degs_marks_names
  
  peaks_density_limits_ls <- list()
  for (n in c("DEGs",marks_names)){
    peaks_density_limits_ls[[n]] <- as.numeric(
      strsplit(peaks_density_limits[[n]],",")[[1]]
    )
  }
  
  peaks_density_limits_Rnd_ls <- list()
  for (n in c("DEGs",marks_names)){
    peaks_density_limits_Rnd_ls[[n]] <- as.numeric(
      strsplit(peaks_density_limits_Rnd[[n]],",")[[1]]
    )
  }
  #--------------------------------- Choose Mode -------------------------------
  if(runMode=="just_plots"){
    #---------------------------------Import bins-------------------------------
    marks_inTADs_with_DEGs_InsideOut_bins_ls <-
      readRDS(paste0(bins_paths,"/marks_inTADs_with_DEGs_InsideOut_bins_ls"))
    marks_in_TADs_w_rndGenes_sample_InsideOut_bins_ls <-
      readRDS(
        paste0(bins_paths,"/marks_in_TADs_w_rndGenes_sample_InsideOut_bins_ls"))
    marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED_ls <-
      readRDS(
        paste0(bins_paths,"/marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED_ls"))
    marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED_ls <-
      readRDS(
        paste0(bins_paths,"/marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED_ls"))
    marks_inTADs_with_DEGs_InsideOut_anchors_ls_ls <-
      readRDS(
        paste0(bins_paths,"/marks_inTADs_with_DEGs_InsideOut_anchors_ls_ls"))
    marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls_ls <-
      readRDS(
        paste0(bins_paths,"/marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls_ls"))
    marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls_ls <-
      readRDS(
        paste0(bins_paths,"/marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls_ls"))
    marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_FLIPPED_ls_ls <-
      readRDS(
        paste0(bins_paths,"/marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_FLIPPED_ls_ls"))
    statistics_ls <-
      readRDS(paste0(bins_paths,"/stats_ls"))
    odds_ratio_Marks_ls <- 
      readRDS(paste0(bins_paths,"/odds_ratio_Marks_ls"))
    odds_ratio_DEGs <-
      readRDS(paste0(bins_paths,"/odds_ratio_DEGs"))
    p_value_Marks_ls <-
      readRDS(paste0(bins_paths,"/p_value_Marks_ls"))
    p_value_DEGs <-
      readRDS(paste0(bins_paths,"/p_value_DEGs"))
    #---------------------------------Plot Genes--------------------------------
    if (normalize_by=="byGenes"){
      normalizador_domains <- statistics_ls[["Genes"]][[Plot=="Domains"]]$DEGs
      normalizador_domains_Rnd <- statistics_ls[["Genes"]][[Plot=="Domains"]]$rGenes
      normalizador_anchors <- statistics_ls[["Genes"]][[Plot=="Anchors"]]$DEGs
      normalizador_anchors_Rnd <- statistics_ls[["Genes"]][[Plot=="Anchors"]]$rGenes
    } else if (normalize_by=="byDomain"){
      normalizador_domains <- statistics_ls[["Genes"]][[Plot=="Domains"]]$dTADs
      normalizador_domains_Rnd <- statistics_ls[["Genes"]][[Plot=="Domains"]]$rTADs
      normalizador_anchors <- statistics_ls[["Genes"]][[Plot=="Anchors"]]$dTADs
      normalizador_anchors_Rnd <- statistics_ls[["Genes"]][[Plot=="Anchors"]]$rTADs
    } else {
      normalizador_domains <- 1
      normalizador_domains_Rnd <- 1
      normalizador_anchors <- 1
      normalizador_anchors_Rnd <- 1
    }
    
    marks_inTADs_with_DEGs_InsideOut_bins_ls[["Genes"]][,N:=N/normalizador_domains]
    marks_inTADs_with_DEGs_InsideOut_bins_ls[["Genes"]][,agg_score:=agg_score/normalizador_domains]
    marks_in_TADs_w_rndGenes_sample_InsideOut_bins_ls[["Genes"]][,N:=N/normalizador_domains_Rnd]
    marks_in_TADs_w_rndGenes_sample_InsideOut_bins_ls[["Genes"]]$agg_score/normalizador_domains_Rnd
    marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED_ls[["Genes"]][,N:=N/normalizador_domains]
    marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED_ls[["Genes"]][,agg_score:=agg_score/normalizador_domains]
    marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED_ls[["Genes"]][,N:=N/normalizador_domains_Rnd]
    marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED_ls[["Genes"]][,agg_score:=agg_score/normalizador_domains_Rnd]
    marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls_ls[["Genes"]][,N:=N/normalizador_anchors_Rnd]
    marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls_ls[["Genes"]][,agg_score:=agg_score/normalizador_anchors_Rnd]
    marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls_ls[["Genes"]][,N:=N/normalizador_anchors]
    marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls_ls[["Genes"]][,agg_score:=agg_score/normalizador_anchors]
    
    # Domain plots
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
    marks_inTADs_with_DEGs_InsideOut_anchors_ls_ls[["Genes"]][,N:=N/normalizador_anchors]
    marks_inTADs_with_DEGs_InsideOut_anchors_ls_ls[["Genes"]][,agg_score:=agg_score/normalizador_anchors]
    plots_density_DEGs_overTADs_inAnchors <- 
      density_f(marks_inTADs_with_DEGs_InsideOut_anchors_ls_ls[["Genes"]],
                density_plot_y, aggreagation_stat[1],
                peaks_density_limits_ls$DEGs[1],
                peaks_density_limits_ls$DEGs[2],
                "black", peaks_density_limits_ls$DEGs[3],
                peaks_density_limits_ls$DEGs[4], "lfc",
                paste0("p=",statistics_ls[["Genes"]][Plot=="Anchors"]$p_value,
                       "\n OR=",statistics_ls[["Genes"]][Plot=="Anchors"]$OR),
                4,axis_bined_boundaries_ls)
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
                axis_bined_boundaries_ls)
    
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
                4,axis_bined_boundaries_ls)
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
                4, axis_bined_boundaries_ls)
    #---------------------- Plot Marks -----------------------------------------
    for (i in side_peaks_ls){
      marks_inTADs_with_DEGs_InsideOut_bins_ls[[i]][,N:=N/normalizador_domains]
      marks_inTADs_with_DEGs_InsideOut_bins_ls[[i]][,agg_score:=agg_score/normalizador_domains]
      marks_in_TADs_w_rndGenes_sample_InsideOut_bins_ls[[i]][,N:=N/normalizador_domains_Rnd]
      marks_in_TADs_w_rndGenes_sample_InsideOut_bins_ls[[i]]$agg_score/normalizador_domains_Rnd
      marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED_ls[[i]][,N:=N/normalizador_domains]
      marks_inTADs_with_DEGs_InsideOut_bins_FLIPPED_ls[[i]][,agg_score:=agg_score/normalizador_domains]
      marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED_ls[[i]][,N:=N/normalizador_domains_Rnd]
      marks_in_TADs_w_rndGenes_sample_InsideOut_bins_FLIPPED_ls[[i]][,agg_score:=agg_score/normalizador_domains_Rnd]
      marks_inTADs_with_DEGs_InsideOut_anchors_ls_ls[[i]][,N:=N/normalizador_anchors]
      marks_inTADs_with_DEGs_InsideOut_anchors_ls_ls[[i]][,agg_score:=agg_score/normalizador_anchors]
      marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls_ls[[i]][,N:=N/normalizador_anchors_Rnd]
      marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_ls_ls[[i]][,agg_score:=agg_score/normalizador_anchors_Rnd]
      marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls_ls[[i]][,N:=N/normalizador_anchors]
      marks_inTADs_with_DEGs_InsideOut_anchors_FLIPPED_ls_ls[[i]][,agg_score:=agg_score/normalizador_anchors]
      marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_FLIPPED_ls_ls[[i]][,N:=N/normalizador_anchors_Rnd]
      marks_in_TADs_w_rndGenes_sample_InsideOut_anchors_FLIPPED_ls_ls[[i]][,agg_score:=agg_score/normalizador_anchors_Rnd]
      
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
                  4,axis_bined_boundaries_ls)
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
                  4,axis_bined_boundaries_ls)
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
                  4,axis_bined_boundaries_ls)
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
                  4,axis_bined_boundaries_ls)
    }
  } else {
    #
    #----------------------------Import TADs and DEGs-----------------------------
    #
    tads <- fread(paste0(tads_dir,tads_file),select = c(1:7))
    degs <- fread(paste0(genes_dir,genes_file),
                  select = c(1:9,11), col.names = c(peaks_cols,"cluster"))
    genes_ctrl_bed <- fread(paste0(genes_dir,genes_ctrl),
                            select = c(1:9,13), col.names = c(peaks_cols,"cluster"))
    genes_rndm <- genes_ctrl_bed[!peak%in%degs$peak]
    #
    #------------------------------Import_Marks-----------------------------------
    #
    names(marks_files) <- marks_names
    #
    marks_ls <- list()
    for (mark in marks_names){
      marks_ls[[mark]] <- fread(marks_files[[mark]],col.names = peaks_cols)
    }
    #
    #------------------------------Make Marks List--------------------------------
    #
    degs_cluster <- degs[cluster%in%c(nomen)]
    degs_cluster_list <- list(DEGs=degs_cluster)
    genes_rndm_list <- list(genes_rndm_All=genes_rndm)
    peaks_list <- c(degs_cluster_list,genes_rndm_list,marks_ls)
    #
    #------------------------------ Make boundaries ------------------------------
    #
    tads_dt <- data.table(chr=tads[[1]],
                          start=floor(tads[[3]]/as.numeric(resol))*resol,
                          end=floor(tads[[5]]/as.numeric(resol))*resol,
                          id=tads[[7]])
    tads_dt$length <- tads_dt$end-tads_dt$start
    #
    #------------------------------ Make TAD extension -----------------------
    #
    tads_dt$extension <- tads_dt$length*TAD_extend
    tads_ext_dt <- data.table(chr=tads_dt$chr,
                              start=tads_dt$start-tads_dt$extension,
                              end=tads_dt$end+tads_dt$extension,
                              id=tads_dt$id,
                              sense=".",
                              upBoundary=tads_dt$start,
                              downBoundary=tads_dt$end)
    #
    #------------------------------ Overlaps Marks ---------------------------
    #
    colors_ls <- extract_colors_f(marks_ls)
    peaks_overlaps_TADs_ls <- overlaps_f(peaks_list,tads_ext_dt)
    #
    #------------------------------ Compute Fisher Statistics --------------------
    #
    DEGs_inTADs_InsideOut_dt <-
      identify_InsideOut_f(peaks_overlaps_TADs_ls[["DEGs"]])
    genesRndm_inTADs_InsideOut_dt <-
      identify_InsideOut_f(peaks_overlaps_TADs_ls[["genes_rndm_All"]])
    
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
    #------------------------- Plot ANCHOR density plots --------------------
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
    #------------------------- Flip DEGs density plots ---------------------------
    #
    DEGs_inTADs_InsideOut_FLIPPED_bins <- binmaker_redim_f(
      DEGs_inTADs_InsideOut_dt,"score",
      aggreagation_stat[2],nbins,TRUE)
    
    genesRndm_inTADs_InsideOut_sample_FLIPPED_bins <- binmaker_redim_f(
      genesRndm_inTADs_InsideOut_sample_dt,"log_TPMs",
      aggreagation_stat[2],nbins,TRUE)
    #
    #--------------------- Flip DEGs for ANCHOR density plots --------------------
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
    #---------------------- List of Bins and statistics --------------------------
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
    #-------------------------- Plot density for DEGs ----------------------------
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
                4,axis_bined_boundaries_ls)
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
                axis_bined_boundaries_ls)
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
                4,axis_bined_boundaries_ls)
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
                4, axis_bined_boundaries_ls)
    #
    #
    #---------------------------ITERATION OVER MARKS--------------------------
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
      #--------------------------- Copy Marks to Plot ----------------------
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
      #--------------------------- Copy Marks to Export ----------------------
      #                            AND STATISTICAL TEST
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
                  4,axis_bined_boundaries_ls)
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
                  4,axis_bined_boundaries_ls)
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
                  4,axis_bined_boundaries_ls)
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
                  4,axis_bined_boundaries_ls)
      
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
}