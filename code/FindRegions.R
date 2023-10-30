###############################################################################
###############################################################################
## Example code for regions detection 
## Input files and code are limited to chr19 
###############################################################################
###############################################################################

library(pacman)
p_load(data.table, dplyr, GenomicRanges, gUtils, viridis, rtracklayer, BSgenome.Mmusculus.UCSC.mm10, foreach, doParallel)
options(scipen=999)
genome <- BSgenome.Mmusculus.UCSC.mm10


get_dens_reg <- function(regionai, densitis) {
  # Simple function calculates median density per region
  ov <- findOverlaps(regionai, densitis)  
  
  ret <- cbind(
    densitis[subjectHits(ov)] %>% 
      as.data.frame() %>% 
      as.data.table() %>% 
      .[, .(score)],
    regionai[queryHits(ov)] %>% 
      as.data.frame() %>% 
      as.data.table() %>% 
      .[, ID := paste0(seqnames, "_", start, "_", end)] %>% 
      .[, .(seqnames, start, end, peakID, testi_up, testi_down)] 
  ) %>% as.data.table %>% 
    .[, lapply(.SD, median), by=c("seqnames", "start", "end", "peakID")] %>% 
    makeGRangesFromDataFrame(., keep.extra.column=TRUE)
  return(ret)  
}

get_regs <- function(pikai, d2test, th) {
  # Recursive function evaluates and extends potential regions
  if(sum(mcols(pikai)$testi_up > -2) > 0 | sum(mcols(pikai)$testi_down > -2) > 0) {
    upai <- resize(pikai, 10, fix="start")
    downai <- resize(pikai, 10, fix="end")

    upai_stat <- get_dens_reg(upai, d2test) 
    prideti <- upai[(!mcols(upai)$peakID %in% mcols(upai_stat)$peakID)]
    if(length(prideti) > 0) {mcols(prideti)$score <- -1}
    upai_stat <- c(upai_stat, prideti)

    downai_stat <- get_dens_reg(downai, d2test)
    prideti <- downai[(!mcols(downai)$peakID %in% mcols(downai_stat)$peakID)]
    if(length(prideti) > 0) {mcols(prideti)$score <- -1}
    downai_stat <- c(downai_stat, prideti)

    new_Peaks <- merge(
      upai_stat  %>% as.data.frame() %>% as.data.table() %>% .[, .(peakID, score, testi_up)] %>% setnames("score", "upscore"),
      downai_stat  %>% as.data.frame() %>% as.data.table() %>% .[, .(peakID, score, testi_down)] %>% setnames("score", "downscore"),
      by="peakID") %>% 
      merge(., pikai %>% as.data.frame() %>% as.data.table() %>% .[, .(seqnames, start, end, peakID)]) %>% 
      as.data.table() %>% 
      .[testi_up > -2, start := start - 10] %>% 
      .[testi_down > -2, end := end + 10 ] %>% 
      .[upscore <= th, testi_up := testi_up - 1] %>% 
      .[upscore > th, testi_up := 0] %>% 
      .[downscore <= th, testi_down := testi_down - 1] %>% 
      .[downscore > th, testi_down := 0] %>% 
      .[, .(peakID, testi_up, testi_down, seqnames, start, end)] %>% 
      makeGRangesFromDataFrame(., keep.extra.columns=TRUE)  
    return(get_regs(new_Peaks, d2test, th))} else {return(pikai)}   
}

regions_perSample <- function(x, samplasName, targetName, thas, thas1, cov2test, chr2do, q1, q2) {
  xx <- foreach(chr=chr2do) %dopar% {
      set.seed(1987)
      dens2test <- x[seqnames(x) == chr]
      # Identify summits and reduce/collapse nearby
      denSummits <- dens2test[mcols(dens2test)$score>thas, ]
      denSummits <- reduce(denSummits, min.gapwidth=10)
      mcols(denSummits)$testi_up <- 0
      mcols(denSummits)$testi_down <- 0
      mcols(denSummits)$peakID <- paste0(seqnames(denSummits), "_", start(denSummits), "_", end(denSummits));
      rezu <- get_regs(denSummits, dens2test, thas1)
      # resize regions by 20 bp because we extened this before stopping.
      regionai <- reduce(rezu-20)
      # Identify first and last targets from coverage. Regions should be limited by cov>0 target. 
      i <- countOverlaps(regionai, cov2test[elementMetadata(cov2test)[, samplasName] > 0] )
      regionai <- regionai[i>=2]
      regionai <- reduce(regionai)
      ov_f <- findOverlaps(regionai, cov2test[elementMetadata(cov2test)[, samplasName] > 0], select="first")
      ov_l <- findOverlaps(regionai, cov2test[elementMetadata(cov2test)[, samplasName] > 0], select="last")
    
      modR <- cbind(regionai %>% as.data.table, 
        cov2test[elementMetadata(cov2test)[, samplasName] > 0][ov_f] %>% as.data.table %>% .[, .(start)] %>% setnames(., "pradzia")) %>% 
      cbind(.,    
        cov2test[elementMetadata(cov2test)[, samplasName] > 0][ov_l] %>% as.data.table %>% .[, .(start)] %>% setnames(., "pabaiga")) %>% 
        .[, start := pradzia] %>% 
        .[, end := pabaiga] %>% 
        .[start < end, ] %>% 
        makeGRangesFromDataFrame()

      modR <- reduce(modR, min.gapwidth=75)
      modR <- modR[width(modR) >= 25]
      mcols(modR)$peakID <- paste0(seqnames(modR), "_", start(modR), "_", end(modR))
      mcols(modR)$score <- 1
      #rtracklayer::export.bw(modR, con =  paste0(dir2save))
  }}

findRengions <- function(targetai, megas, chr2do, quant1, quant2) {
  # This function calculates initial thresholds and calls other functions to identify regions per sample
  # targetai <- "open"
  # megas <- "TT_S0"; sample ID
  # 75 and 75 - first one is threshold for seeds and second one is extension threshold. These numbers do not have to be the same. 

  coveragai <- readRDS("input/cov2test.RDS")
  dens2do <-  import.bw("input/dens2test.bw")
  
  # Thresholds should be set using  data from the whole genome and not just a  single chromosome
  set.seed(1987)
  thas <- dens2do[sample(1:length(dens2do), 1000000)] %>% 
      as.data.frame() %>% 
      as.data.table() %>% 
      .[, score] %>% 
      quantile(quant1/100) %>% 
      round(., 2)

  set.seed(1987)
  thas1 <- dens2do[sample(1:length(dens2do), 1000000)] %>% 
      as.data.frame() %>% 
      as.data.table() %>% 
      .[, score] %>% 
      quantile(quant2/100) %>% 
      round(., 2)      
  # When analysing other stages than S0, density of specific stage should be read here.     
  dens2do <- import.bw("input/dens2test.bw")
  gc()
  regions_perSample(dens2do, megas, targetai, thas, thas1, coveragai, chr2do, quant1, quant2);
}

findRengions("open", "TT_S0", "19", 75, 75)