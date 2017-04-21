################################################################################
####Kyle Halliwill 04 18 17
####Converting translocation and segmentation data into a single, unified frame 
####for plotting
##This script is designed to take two data frames, one consisting of segmented
##data from a segmentation algorithm (in this case BICseq), and formatted 
##translocation data (in this case from lumpy), and return a single merged 
##sheet engineered for plotting. 


################################################################################
####Libraries
require(reshape)


################################################################################
####Function
combBICLPY <- function(trans_data,
                       seg_data, 
                       chroms = NULL,
                       lbound = NULL,
                       rbound = NULL,
                       increase_fraction_base = 0.1, 
                       increase_fraction_interval = 0.05){
  
  ##getting chromosomes, if none are provided, grab all available
  if(is.null(chroms)){
    chroms <- unique(c(seg_data$chrom, trans_data$p1_chr, trans_data$p2_chr))
  }
  
  ##The maximum value in the segmentation data, again to identify where 
  ##rearrangements should appear. This is a global parameter, not influenced by 
  ##which chromosome we are currently working on
  dMax <- max(as.numeric(seg_data$log2.copyRatio))
  
  ##Building a results container
  resFRAME <- data.frame("sample" = NA,
                         "type" = NA,
                         "chrom1" = NA,
                         "chrom2" = NA,
                         "x_pos1" = NA,
                         "x_pos2" = NA,
                         "y_value1" = NA,
                         "y_value2" = NA)[ -1, ]
  
  ##Looping over by chromosome, attaching formatted data to the results container
  for(chr in chroms){
    ##This is the fraction above the maximum segmentation value that 
    ##rearrangements will appear
    incFac <- increase_fraction_base
    
    ##Subsetting to each chromosome
    ctDat <- trans_data[ which(trans_data$p1_chr == chr), ]
    csDat <- seg_data[ which(seg_data$chrom == chr), ]
    
    ##Getting the median value of the chromosome. This helps draw prettier 
    ##rearrangements
    tMin <- median(csDat$log2.copyRatio, na.rm = T)
    
    ##Setting boundaries. This is designed to play nice with both multiple 
    ##chromosomes and single boundaries, but won't really make sense to approach
    ##things that way. 
    if(is.null(lbound)){
      lbound <- min(c(ctDat$p1_pos, csDat$start), na.rm = T)
    } 
    if(is.null(rbound)){
      rbound <- max(c(ctDat$p2_pos, csDat$end), na.rm = T)
    } 
    
    ##Initially formatting the segmentation data
    seg_ss <- csDat[ which(csDat$end >= lbound & csDat$start <= rbound), ]
    seg_melt <- melt(seg_ss[ , c("chrom", "start", "end", "log2.copyRatio")], 
                     id.vars = c("chrom", "log2.copyRatio"),
                     measure.vars = c("start", "end"))
    seg_melt <- seg_melt[ order(seg_melt$value), ]
    seg_melt$value[ which(seg_melt$value < lbound)] <- lbound
    seg_melt$value[ which(seg_melt$value > rbound)] <- rbound
    seg_melt$chrom2 <- seg_melt$chrom
    seg_melt$pos2 <- seg_melt$value
    seg_melt$type <- "COPY_NUMBER"
    seg_melt$sample <- seg_data$tmr[1]
    seg_melt$value2 <- seg_melt$log2.copyRatio
    seg_form <- seg_melt[ , c("sample", "type", "chrom", "chrom2", "value", "pos2", "log2.copyRatio", "value2")]
    colnames(seg_form) <- colnames(resFRAME)
    
    ##Now attaching this to the eventual output and formatting
    resFRAME <- rbind(resFRAME, seg_form)
    
    ##Alright, great, now working on the rearrangment data. First only intra
    ##chromosomal events
    trans_cis <- ctDat[ which(ctDat$p1_pos >= lbound &
                                ctDat$p2_pos < rbound & 
                                ctDat$p1_chr == chr &
                                ctDat$p2_chr == chr), ]
    for(tType in c("DEL","INV","ITX")){
      tDatRaw <- trans_cis[ which(trans_cis$var_type == tType), ]
      tDat <- tDatRaw[ duplicated(tDatRaw$var_id_base) == FALSE, ]
      yhighVals <- dMax + dMax*incFac
      if(nrow(tDat) != 0){
        for(i in 1:nrow(tDat)){
          itDat <- tDat[ i, ]
          resFRAME[ nrow(resFRAME) + 1, ] <- c(itDat$samp, paste(itDat$var_type, "_LLEG", sep = ""), 
                                               as.character(itDat$p1_chr), as.character(itDat$p2_chr), 
                                               itDat$p1_pos, itDat$p1_pos, 
                                               tMin, yhighVals)
          resFRAME[ nrow(resFRAME) + 1, ] <- c(itDat$samp, paste(itDat$var_type, "_ARC", sep = ""),
                                               as.character(itDat$p1_chr), as.character(itDat$p2_chr), 
                                               itDat$p1_pos, itDat$p2_pos, 
                                               yhighVals, yhighVals)
          resFRAME[ nrow(resFRAME) + 1, ] <- c(itDat$samp, paste(itDat$var_type, "_RLEG", sep = ""), 
                                               as.character(itDat$p1_chr), as.character(itDat$p2_chr), 
                                               itDat$p2_pos, itDat$p2_pos, 
                                               yhighVals, tMin)
        }
      }
      incFac <- incFac + increase_fraction_interval
    }
    
    ##Ok, and finally, the inter chromosomal events. Only going to be taking the
    ##portion of the event that corresponds to this chromosome
    trans_trans <- ctDat[ which(ctDat$var_type == "CTX" & ctDat$p1_chr == chr), ]
    yhighVals <- dMax + dMax*incFac
    if(nrow(trans_trans) != 0){
      for(i in 1:nrow(trans_trans)){
        itDat <- trans_trans[ i, ]
        resFRAME[ nrow(resFRAME) + 1, ] <- c(itDat$samp, paste(itDat$var_type, "_LLEG", sep = ""), 
                                             as.character(itDat$p1_chr), as.character(itDat$p1_chr), 
                                             itDat$p1_pos, itDat$p1_pos, 
                                             tMin, yhighVals)
        resFRAME[ nrow(resFRAME) + 1, ] <- c(itDat$samp, paste(itDat$var_type, "_ULEG", sep = ""),
                                             as.character(itDat$p1_chr), as.character(itDat$p1_chr),
                                             as.numeric(itDat$p1_pos), as.numeric(itDat$p1_pos) + ifelse(as.numeric(itDat$p2_chr) > as.numeric(itDat$p1_chr),
                                                                                                         (rbound - lbound)/20,
                                                                                                         -(rbound - lbound)/20),
                                             yhighVals, yhighVals + yhighVals/10)
      }
    }
  }
  resFRAME$x_pos1 <- as.numeric(resFRAME$x_pos1)
  resFRAME$x_pos2 <- as.numeric(resFRAME$x_pos2)
  resFRAME$y_value1 <- as.numeric(resFRAME$y_value1)
  resFRAME$y_value2 <- as.numeric(resFRAME$y_value2)
  resFRAME$chrom1 <- factor(resFRAME$chrom1, 
                            levels= c(1:19, "X"),
                            labels = c(1:19, "X"))
  resFRAME$chrom2 <- factor(resFRAME$chrom2, 
                            levels= c(1:19, "X"),
                            labels = c(1:19, "X"))
  
  return(resFRAME)
}





