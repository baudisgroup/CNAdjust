# fit state-specific distribution 
fit.total <- function(seg,label){
  idx <- which(seg$label == label)
  if (length(idx) == 0){return()}
  sel_seg <- seg[idx,]
  sel_seg_len <- ceiling((sel_seg[,4]-sel_seg[,3]+1)/1000000)
  value <- rep(sel_seg[,6],sel_seg_len)
  fit <- MASS::fitdistr(value, "normal")
  return(fit)
}

compute.fit <- function(seg){
  low_dup_fit <- fit.total(seg,'+1')
  high_dup_fit <- fit.total(seg,'+2')
  low_del_fit <- fit.total(seg,'-1')
  high_del_fit <- fit.total(seg,'-2')
  normal_fit <- fit.total(seg,'0')
  
  fit_lst <- list()
  fit_lst[['+1']] <- low_dup_fit
  fit_lst[['+2']] <- high_dup_fit
  fit_lst[['-1']] <- low_del_fit
  fit_lst[['-2']] <- high_del_fit
  fit_lst[['0']] <- normal_fit
  return(fit_lst)
}

compute.pointprob <- function(fit, point, reffit=NULL,label=NULL){
  # if only one called CNV for specific level
  if (is.na(fit$estimate[2])) fit <- reffit[[label]]
  ## P(X <= di)
  if (label %in% c("+1","+2")){
    prob <- pnorm(point, mean=fit$estimate[1], sd=fit$estimate[2], lower.tail = T) 
  ## P(X >= di)
  } else if (label %in% c("-1","-2")){
    prob <- pnorm(point, mean=fit$estimate[1], sd=fit$estimate[2], lower.tail = F)
  ## 1 - (two-sided tail probability of di) 
  } else{
    prob <- 2 * pnorm(fit$estimate[1]-abs(point-fit$estimate[1]), mean=fit$estimate[1], sd=fit$estimate[2], lower.tail = T) 
  }  
  return(prob)
}

# calculate CNA frequency

extract.bin.feature <- function(data,bins,genome='hg38',exclude.sex.chrom = NULL){
  # transform sex chrom
  data <- unify.sexchro(data) 
  bins <- unify.sexchro(bins)

  if (!is.null(exclude.sex.chrom)){
    bins <- bins[!bins[,2] %in% exclude.sex.chrom,]
    data <- data[!data[,2] %in% exclude.sex.chrom,]
  }

  
  feature_dup_low <- list()
  feature_del_low <- list()
  feature_dup_high <- list()
  feature_del_high <- list()
  # for each sample
  for (idx in c(1:length(unique(data[,1])))){
    ind.data <- data[data[,1] %in% unique(data[,1])[idx],]
    ind_dup_low <- rep(0,dim(bins)[1])
    ind_del_low <- rep(0,dim(bins)[1])
    ind_dup_high <- rep(0,dim(bins)[1])
    ind_del_high <- rep(0,dim(bins)[1])
    # for each segment
    for (j in c(1:dim(ind.data)[1])){
      ind.seg.start <- ind.data[j,3]
      ind.seg.end <- ind.data[j,4]
      # find bins where the segment spans
      sel.bin <- which(bins[,2] == ind.data[j,2] & bins[,4] > ind.seg.start & bins[,3] < ind.seg.end)
      if (length(sel.bin) == 0){next}
      ind_dup_high[sel.bin] <-   ind_dup_high[sel.bin] + as.numeric(ind.data[j,'label'] == '+2')
      ind_del_high[sel.bin] <- ind_del_high[sel.bin] + as.numeric(ind.data[j,'label'] == '-2')
      ind_dup_low[sel.bin] <- ind_dup_low[sel.bin] + as.numeric(ind.data[j,'label'] == '+1')
      ind_del_low[sel.bin] <-  ind_del_low[sel.bin] + as.numeric(ind.data[j,'label'] == '-1')
    } # binary representation
    ind_dup_high[ind_dup_high > 1] <- 1
    ind_del_high[ind_del_high > 1] <- 1
    ind_dup_low[ind_dup_low > 1] <- 1
    ind_del_low[ind_del_low > 1] <- 1
    feature_dup_high[[idx]] <- ind_dup_high
    feature_del_high[[idx]] <- ind_del_high
    feature_dup_low[[idx]] <- ind_dup_low
    feature_del_low[[idx]] <- ind_del_low
  }
  feature_dup_high <- do.call(rbind,feature_dup_high)
  feature_del_high <- do.call(rbind,feature_del_high)
  feature_dup_low <- do.call(rbind,feature_dup_low)
  feature_del_low <- do.call(rbind,feature_del_low)
  
  if (length(unique(data[,1])) > 1){
    rownames(feature_dup_low) <- unique(data[,1])
    rownames(feature_dup_high) <- unique(data[,1])
    rownames(feature_del_low) <- unique(data[,1])
    rownames(feature_del_high) <- unique(data[,1])
  }
  
  
  feature.list <- list()
  feature.list[['low.dup']] <- feature_dup_low
  feature.list[['high.dup']] <- feature_dup_high
  feature.list[['low.del']] <- feature_del_low
  feature.list[['high.del']] <- feature_del_high
  feature.list[['samples']] <- unique(data[,1])
  return(feature.list)
}

get.freq <- function(data,bins,output_path=NULL,save_freq=F,genome='hg38',exclude_sex_chrom = NULL,mergelevel=F){
  # special: return flat prior
  if (is.null(data)){
    freq.lst <- list()
    freq.lst[['+1']] <- rep(1,nrow(bins))
    freq.lst[['+2']] <- rep(1,nrow(bins))
    freq.lst[['-1']] <- rep(1,nrow(bins))
    freq.lst[['-2']] <- rep(1,nrow(bins))
    freq.lst[['0']] <- rep(1,nrow(bins))
    return(freq.lst)
  }
  
  features <- extract.bin.feature(data,bins,genome,exclude_sex_chrom)

  freq.lst <- list()
  # merge low-level and high-level calling
  if (mergelevel){
    freq_dup <- features$low.dup+features$high.dup
    freq_dup[freq_dup > 1] <- 1
    freq_dup <- colSums(freq_dup)/nrow(freq_dup)
      
    freq_del <- features$low.del+features$high.del
    freq_del[freq_del > 1] <- 1
    freq_del <- colSums(freq_del)/nrow(freq_del)

    freq.lst[['+1']] <- freq_dup
    freq.lst[['+2']] <- freq_dup
    freq.lst[['-1']] <- freq_del
    freq.lst[['-2']] <- freq_del
    freq.lst[['0']] <- 1-colSums(rbind(freq_dup,freq_del))

  } else{
    freq_low_del <-  features$low.del
    freq_low_del <- colSums(freq_low_del)/nrow(freq_low_del)
  
    freq_low_dup <-  features$low.dup
    freq_low_dup <- colSums(freq_low_dup)/nrow(freq_low_dup)
  
    freq_high_del <-  features$high.del
    freq_high_del <- colSums(freq_high_del)/nrow(freq_high_del)
  
    freq_high_dup <- features$high.dup
    freq_high_dup <- colSums(freq_high_dup)/nrow(freq_high_dup)
  

    freq.lst[['+1']] <- freq_low_dup
    freq.lst[['+2']] <- freq_high_dup
    freq.lst[['-1']] <- freq_low_del
    freq.lst[['-2']] <- freq_high_del
    freq.lst[['0']] <- 1-colSums(rbind(freq_low_dup,freq_low_del, freq_high_dup,freq_high_del))
  }

  if (any(freq.lst[['0']] < 0)) freq.lst[['0']][which(freq.lst[['0']] < 0)] <- 0
  
  if(save_freq){   
    saveRDS(freq.lst,output_path)
  }
  
  return(freq.lst)
}

# posterior calculation
compute.posterior <- function(output_dir,reference_freq_df,prior_code,series.seg,seg,bins_info,fit,reffit,return_details=F,genome='hg38',calling = F,priorshift='n',shift_samples=NULL,shiftnum=NULL){
  if (length(prior_code) > 0 & !is.na(prior_code) & !is.null(reference_freq_df)){
    reference_freq <- reference_freq_df[[prior_code]]
    prior <- list()
    prior[['+1']] <- reference_freq[c(1:(length(reference_freq)/2))]
    prior[['+2']] <- prior[['+1']] 
    prior[['-1']] <- reference_freq[c(((length(reference_freq)/2)+1):length(reference_freq))]
    prior[['-2']] <- prior[['-1']] 
    prior[['0']] <- 1- prior[['+1']] - prior[['-1']]
    # in case diverse cnv regions in single bin
    if (any(prior[['0']] < 0)) prior[['0']][which(prior[['0']] < 0)] <- 0
  } else{
    # if metadata is not available, use avg freq of all samples, which means all samples from the input segment data should be similar 
    freq_obj_path <- file.path(output_dir,'freq',paste0('CNV',"-",priorshift,'-',shiftnum,'-freq.rds'))
    if (file.exists(freq_obj_path)){
      prior <- readRDS(freq_obj_path)
    } else{
      dir.create(file.path(output_dir,'freq'),showWarnings = F)
      data <- series.seg
      ## consider situations where prior needs to be adjusted
      ### if the series only contains 1 sample
      if (length(unique(data[,1])) == 1){
        data <- NULL
        ### if whole shift is considered
      } else if (priorshift != 'n'){
        data[data[,1] %in% shift_samples,] <- labelSeg::labelseg(data[data[,1] %in% shift_samples,],baseshift = priorshift,genome=genome,shiftnum=shiftnum,labeled=!calling)
      }
      ## compute 
      prior <- get.freq(data,bins_info,output_path = freq_obj_path,save_freq = T,genome=genome,mergelevel=T)
    }
  }
  # missing levels for adjusted callings
  if (any(!unique(seg$label) %in% names(fit))){
    misslevels <- unique(seg$label)[!unique(seg$label) %in% names(fit)]
    for (level in misslevels){
      fit[[level]] <- reffit[[level]]
    }
  } 
  
  seg <- unify.sexchro(seg) 
  bins_info <- unify.sexchro(bins_info)
  
  posterior <- vapply(seq_len(dim(seg)[1]), function(i){
    x <- seg[i,]
    ind_fit <- fit[[x$label]]
    lld <- compute.pointprob(ind_fit, x[,6], reffit, x$label) 
    ind_prior <- prior[[x$label]]
    ind_bins <- bins_info[bins_info[,2] == x[,2],]
    ind_intersect_bins <- ind_bins[ind_bins[,3] < x[,4] & ind_bins[,4] > x[,3],'index']
    ind_avg_prior <- mean(ind_prior[ind_intersect_bins])
    return((ceiling((x[,4]-x[,3])/1000000)) * log10(ind_avg_prior * lld+1e-25))
  },numeric(1))
  
  if (return_details){
    return(list(prior=prior,fit=fit,posterior=posterior))
  } else{
    return(sum(posterior))
  }
}
