# determine if whole shift is applied (shift all samples with bias in one direction and leave all other samples unchanged, including samples with bias in the opposite direction)
compare.whole.shift <- function(output_dir,report,reference_freq_df,prior_code,ori_label_seg,shift_samples,total_shift_sample_prop,baseshift,ori_fit,ref_fit,bins_lst,genome,calling){
  if (length(shift_samples) == 0) return()
  # original data 
  shift_sample_len <- length(shift_samples)
  shift_sampleids <- rep(NA, shift_sample_len)

  total_ori_posterior <- rep(0,length(unique(ori_label_seg[,1])))
  for (i in seq_len(shift_sample_len)){
    shift_sampleids[i] <- report$sample_id[report$sample_name == shift_samples[i]]
    total_ori_posterior[i] <- compute.posterior(output_dir,reference_freq_df,prior_code,ori_label_seg,ori_label_seg[ori_label_seg[,1] == shift_sampleids[i],],bins_lst,ori_fit,ref_fit,genome=genome,priorshift = 'n')
  }  

  opposite_shift_sampleids <- unique(ori_label_seg[,1])[!unique(ori_label_seg[,1]) %in% shift_sampleids]
  
  for (j in seq_len(length(opposite_shift_sampleids))){
    total_ori_posterior[i+j]  <- compute.posterior(output_dir,reference_freq_df,prior_code,ori_label_seg,ori_label_seg[ori_label_seg[,1] == opposite_shift_sampleids[j],],bins_lst,ori_fit,ref_fit,genome=genome,priorshift = 'n')
  }

  total_ori_posterior <- sum(total_ori_posterior)

  # relabeled data 
  current_total_shift_posterior <- -1e100
  max_total_shift_posterior <- -Inf
  current_seg <- ori_label_seg
  # determine shift strength which is indicated by k
  k <- 1
  while (current_total_shift_posterior > max_total_shift_posterior){
      max_total_shift_posterior <- current_total_shift_posterior
      finalseg <- current_seg
      new_label_seg <- ori_label_seg
      new_label_seg[new_label_seg[,1] %in% shift_sampleids,] <- labelSeg::labelseg(new_label_seg[new_label_seg[,1] %in% shift_sampleids,],baseshift = baseshift,genome=genome,shiftnum=k,labeled=!calling)
      # if most samples have abnormal baselines, fit itself doesn't make sense for comparing different callings
      if (total_shift_sample_prop > 0.5){
        shift_fit <- ref_fit
      } else{
        shift_fit <- compute.fit(new_label_seg)
      }
      total_shift_posterior <- rep(0,length(unique(new_label_seg[,1])))

      for (i in seq_len(length(unique(new_label_seg[,1])))){
        total_shift_posterior[i] <- compute.posterior(output_dir,reference_freq_df,prior_code,ori_label_seg,new_label_seg[new_label_seg[,1] == unique(new_label_seg[,1])[i],],bins_lst,shift_fit,ref_fit,genome=genome,calling=calling,priorshift = baseshift,shift_samples=shift_sampleids,shiftnum=k)
      }

      current_total_shift_posterior <- sum(total_shift_posterior)
      current_seg <- new_label_seg
      k <- k+1
  }


  if (total_ori_posterior  >= max_total_shift_posterior){
    shift <- "n"
    shiftnum <- NULL
    shiftsampleids <- NULL
    posterior <- total_ori_posterior 
    fit <- ori_fit
  } else{
    shift <- baseshift
    shiftnum <- k-2
    shiftsampleids <- shift_sampleids
    posterior <- max_total_shift_posterior
    fit <- compute.fit(finalseg)
  }
  
  compare_list <- list()
  compare_list[['shift']] <- shift
  compare_list[['shiftnum']] <- shiftnum
  compare_list[['shiftids']] <- shiftsampleids
  compare_list[['posterior']] <- posterior
  compare_list[['fit']] <- fit
  return(compare_list)
}

# shift the baseline of specific sample, and update plot and report accordingly
shift.baseline <- function(output_dir,cohort_seg,oriseg,prior_code,reference_freq_df,sample,report,bins_lst,fit,ref_fit,genome,calling,shift,whole_shift,shift_ids,shiftnum,outputfilename,outputplotname){
  labelplot_path <- file.path(output_dir ,sample, outputplotname)
  renameplot_path <- file.path(output_dir ,sample,gsub("\\.","_before_shift.",outputplotname))
  
  ori_posterior <- compute.posterior(output_dir,reference_freq_df,prior_code,cohort_seg,oriseg,bins_lst,fit,ref_fit,genome=genome,calling=calling,priorshift=whole_shift,shift_samples=shift_ids,shiftnum=shiftnum)

  relabelSeg <- oriseg
  current_relabel_posterior <- -1e100
  max_relabel_posterior <- -Inf
  k <- 1
  while (current_relabel_posterior > max_relabel_posterior){
      max_relabel_posterior <-  current_relabel_posterior
      finalrelabelSeg <- relabelSeg
      relabelSeg <- labelSeg::labelseg(data = oriseg,baseshift = substr(shift,1,1),genome=genome,shiftnum=k,labeled=!calling)
      current_relabel_posterior <- compute.posterior(output_dir,reference_freq_df,prior_code,cohort_seg,relabelSeg,bins_lst,fit,ref_fit,genome=genome,calling=calling,priorshift=whole_shift,shift_samples=shift_ids,shiftnum=shiftnum)
      k <- k+1
  }

  if (max_relabel_posterior > ori_posterior){
      if (file.exists(labelplot_path)) system(sprintf("mv %s %s",labelplot_path,renameplot_path))
      # if high-level calling doesn't exist in original calling
      if (!any(c("+2","-2") %in% unique(cohort_seg$label))){
        finalrelabelSeg$label[finalrelabelSeg$label == "+2"] <- "+1"
        finalrelabelSeg$label[finalrelabelSeg$label == "-2"] <- "-1"
      }
      cov <- seg.cov(finalrelabelSeg)
      report <- update.seg(output_dir=output_dir,sample=sample,oriseg=oriseg,relabelSeg=finalrelabelSeg,report=report,score=cov,problem=paste0("shift-baseline-",shift),genome=genome,outputfilename=outputfilename,outputplotname=outputplotname)
  } else{
      idx <- which(report$sample_name == sample)
      report$note[idx] <- paste0(report$note[idx],";shift-baseline-",shift,"_reverse")
  }
  
  return(report)
}
