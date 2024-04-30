cutoff.calling <- function(seg,low_thre,high_thre){
  lrr <- seg[,6]
  recall <- rep('0',length(lrr))
  recall[lrr >= low_thre] <- '+1'
  recall[lrr <= -low_thre] <- '-1'
  recall[lrr >= high_thre] <- '+2'
  recall[lrr <= -high_thre] <- '-2'
  seg$label <- recall
  return(seg)
}

process.noise <- function(output_dir,cohort_seg,oriseg,prior_code,reference_freq_df,sample,report,bins_lst,ori_fit,ref_fit,genome,calling,whole_shift,shift_id,shiftnum,outputfilename,outputplotname, lowthres, highthres){
  assess <- compute.posterior(output_dir,reference_freq_df,prior_code,cohort_seg,oriseg,bins_lst,ori_fit,ref_fit,genome=genome,calling=calling,priorshift=whole_shift,shift_samples=shift_id,shiftnum=shiftnum)
  target <- oriseg
  method <- 'nochange'
  for ( l_thre in lowthres){
    for (h_thre in highthres){
      relabelSeg <- cutoff.calling(oriseg,l_thre,h_thre)
      reassess <-  compute.posterior(output_dir,reference_freq_df,prior_code,cohort_seg,relabelSeg,bins_lst,ori_fit,ref_fit,genome=genome,calling=calling,priorshift=whole_shift,shift_samples=shift_id,shiftnum=shiftnum)
      if (reassess > assess){
        assess <- reassess
        target <- relabelSeg
        method <- paste0('lowthre_',l_thre,'_highthre_',h_thre)
      }
    }
  }
  if (method != 'nochange'){
    # if high-level calling doesn't exist in original calling
    if (!any(c("+2","-2") %in% unique(cohort_seg$label))){
        target$label[target$label == "+2"] <- "+1"
        target$label[target$label == "-2"] <- "-1"
    }
    cov <- seg.cov(target)
    report <- update.seg(output_dir=output_dir,sample=sample,oriseg=oriseg,relabelSeg=target,report=report,score=cov,problem=paste0("noisy_",method),genome=genome,outputfilename=outputfilename,outputplotname=outputplotname)
  } else{
    idx <- which(report$sample_name == sample)
    report$note[idx] <- paste0(report$note[idx],";noisy_",method)
  }
  
  return(report)
}

