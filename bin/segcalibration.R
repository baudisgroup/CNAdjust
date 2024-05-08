#! /usr/bin/env Rscript

args = R.utils::commandArgs(asValue = T)

work_dir <- args$workdir
data_dir <- args$inputdir
output_dir <- args$outputdir
usetumorref <- as.logical(args$usetumorref)
tumorcodefilename <- args$cohortfile
priorfilename <- args$priorfile
use_custom_region <- as.logical(args$useregion)
regionfilename <- args$regionfile
genome <- args$genome
outputfilename <- args$outputfile
outputplotname <- args$outputplot

logrsd_thre <- as.numeric(args$logrsd)
lowfrac_thre <- as.numeric(args$cnafrac)
callratio_thre <- as.numeric(args$dupdelratio)
lowfrac_thre2 <- as.numeric(args$cnafrac2)
callratio_thre2 <- as.numeric(args$dupdelratio2)
lowfrac_thre3 <- as.numeric(args$cnafrac3)
segnum_thre <- as.numeric(args$segnum)
lowcutoff <- as.numeric(unlist(strsplit(args$lowthre,split=',')))
highcutoff <- as.numeric(unlist(strsplit(args$highthre,split=',')))

# source helper functions
source(file.path(work_dir,'bin','utils','general-utils.R'))
source(file.path(work_dir,'bin','utils','computePosterior-utils.R'))
source(file.path(work_dir,'bin','utils','shiftBaseline-utils.R'))
source(file.path(work_dir,'bin','utils','processNoisy-utils.R'))
plotfiles_sources <- list.files(path = file.path(work_dir,'bin','utils','plotseg-utils'),recursive = T,full.names = T)
invisible(sapply(plotfiles_sources,source))
load(file.path(work_dir,'data','plot_seg_genome.rda'))
source(file.path(work_dir,'bin','plotseglabel.R'))


# read qc report
report_path <- checkfile(filepath=output_dir,filetype="QC report",pattern="data_quality_report.txt")
report <- read.table(report_path,sep = '\t',header=T)
# read tumor code annotation file and tumor referecne frequency file 
if (usetumorref){
    series <- basename(output_dir)
    prior_code_file <- checkfile(filepath=file.path(data_dir,series),filetype="cohort assign file",pattern=tumorcodefilename)
    reference_prior_file <- checkfile(filepath=file.path(data_dir,series), filetype="cohort CNA prior file", pattern=priorfilename)  
    prior_code_df <- read.table(prior_code_file,sep = '\t',header=T,check.names = F)
    prior_prob_df <-  read.table(reference_prior_file,sep = '\t',header=T,check.names = F)
    # check if these two files are matched
    unique_prior_codes <- unique(prior_code_df[,2])
    if (!all(unique_prior_codes %in% colnames(prior_prob_df))) stop("prior is not found for ",unique_prior_codes[!unique_prior_codes %in% colnames(prior_prob_df)])
    if (!all(report$sample_name %in% prior_code_df[,1])) stop("sample mapping is failed between segment data and cohort assign")
} else{
    prior_code_df <- NULL
    prior_prob_df <- NULL
}

if (use_custom_region){
    bins_info_path <- checkfile(filepath=file.path(data_dir,series),filetype="genome region file",pattern=regionfilename)    
} else{
    bins_info_path <- checkfile(filepath=file.path(work_dir,'data'),filetype="genome region file",pattern=paste0(genome,"_bin.txt"))
}

bins_info <- read.table(bins_info_path ,sep="\t",header=T)

# identify problematic samples
fail_label_idx <- which(report$lowCNA_cov == 0 & report$highCNA_cov == 0 & report$normal_cov == 0)

noisy_idx <- which(report$segment_num > segnum_thre)
noisy_idx <-  setdiff(noisy_idx ,fail_label_idx)

h_shift_idx1 <- which(report$LLR_segsd > logrsd_thre & report$lowCNA_cov > lowfrac_thre & report$lowCNA_ratio >= callratio_thre)
h_shift_idx2 <- which(report$lowCNA_cov > lowfrac_thre2 & report$lowCNA_ratio >= callratio_thre2)
h_shift_idx3 <- which(report$lowCNA_cov > lowfrac_thre3 & report$lowCNA_ratio >= 1)
h_shift_idx <- unique(c(h_shift_idx1, h_shift_idx2, h_shift_idx3))
h_shift_idx <- setdiff(setdiff(h_shift_idx,fail_label_idx),noisy_idx)

l_shift_idx1 <- which(report$LLR_segsd > logrsd_thre  & report$lowCNA_cov > lowfrac_thre & report$lowCNA_ratio <= 1/callratio_thre)
l_shift_idx2 <- which(report$lowCNA_cov > lowfrac_thre2 & report$lowCNA_ratio <= 1/callratio_thre2)
l_shift_idx3 <- which(report$lowCNA_cov > lowfrac_thre3 & report$lowCNA_ratio <= 1)
l_shift_idx <- unique(c(l_shift_idx1, l_shift_idx2, l_shift_idx3))
l_shift_idx <- setdiff(setdiff(l_shift_idx,fail_label_idx),noisy_idx)

noisy_samples <- report$sample_name[noisy_idx]
h_shift_samples <- report$sample_name[h_shift_idx]
l_shift_samples <- report$sample_name[l_shift_idx]


# read labeled segment data
ori_label_seg <- get.labelseg(output_dir,outputfilename)
if (is.null(ori_label_seg)) stop("labeled data is not found in output dir")

num_cohorts <- ifelse(is.null(prior_code_df),1,length(unique(prior_code_df[,2])))
ref_fit <- readRDS(file.path(work_dir,'data',"reference_fit.rds"))

for (i in seq_len(num_cohorts)){
    if (is.null(prior_code_df)){
        cohort_code <- NA
        cohort_samples <- report$sample_name
      } else{
        cohort_code <- unique(prior_code_df[,2])[i]
        cohort_samples <- prior_code_df[,1][prior_code_df[,2] == cohort_code]
    }
    cohort_noisy_samples <- noisy_samples[noisy_samples %in% cohort_samples] 
    cohort_h_shift_samples <- h_shift_samples[h_shift_samples %in% cohort_samples]
    cohort_l_shift_samples <- l_shift_samples[l_shift_samples %in% cohort_samples]
    cohort_oriseg <- ori_label_seg[ori_label_seg[,1] %in% report$sample_id[report$sample_name %in% cohort_samples],]
    
    if (dim(cohort_oriseg)[1] == 0) next
    
    # fit level-specific logR distribution for the same cohort 
    total_adjust_sample_prop <- sum(length(cohort_h_shift_samples),length(cohort_l_shift_samples),length(cohort_noisy_samples))/length(unique(cohort_oriseg[,1]))
    ## if majority of samples have potential problems, fit itself doesn't make sense for adjustment
    if (total_adjust_sample_prop > 0.5){
        fit <- ref_fit
    } else{
        fit <- compute.fit(cohort_oriseg)
    } 
       
    whole_shift <- "n"
    whole_shift_num <- NULL
    whole_shift_ids <- NULL
    total_shift_sample_prop <- sum(length(cohort_h_shift_samples),length(cohort_l_shift_samples))/length(unique(cohort_oriseg[,1]))
    # if some samples have abnormal baselines, data likelihood and prior (when use frequency) are changed. Then we need to compare posterior between new calling by shifting all of the sample and original callings 
    if (total_shift_sample_prop >= 0.25){
        compare_1 <- compare.whole.shift(output_dir,report,prior_prob_df,cohort_code,cohort_oriseg,cohort_h_shift_samples,total_shift_sample_prop,baseshift = 'h',fit,ref_fit,bins_info,genome,calling=FALSE)
        compare_2 <- compare.whole.shift(output_dir,report,prior_prob_df,cohort_code,cohort_oriseg,cohort_l_shift_samples,total_shift_sample_prop,baseshift = 'l',fit,ref_fit,bins_info,genome,calling=FALSE)
        if (is.null(compare_1)) comapre <- compare_2
        if (is.null(compare_2)) comapre <- compare_1
        if (!is.null(compare_1) & !is.null(compare_2)){
            if (compare_2$posterior > compare_1$posterior){
                compare <- compare_2
            } else{
                compare <- compare_1
            }
        }

        whole_shift <- compare$shift
        whole_shift_num <- compare$shiftnum
        fit <- compare$fit
        whole_shift_ids <- compare$shiftids
    } 
    
    # shift the baseline higher
    if (length(cohort_h_shift_samples) > 0){
        for (sample in cohort_h_shift_samples){
            ind_seg <- cohort_oriseg[cohort_oriseg[,1] == report$sample_id[report$sample_name == sample],]
            report <- shift.baseline(output_dir,cohort_oriseg,ind_seg,cohort_code,prior_prob_df,sample,report,bins_info,fit,ref_fit,genome,calling=FALSE,shift="higher",whole_shift = whole_shift,shift_id=whole_shift_ids,shiftnum=whole_shift_num,outputfilename=outputfilename,outputplotname=outputplotname)
        }
    }
    
    # shift the baseline lower
    if (length(cohort_l_shift_samples) > 0){
        for (sample in cohort_l_shift_samples){
            ind_seg <- cohort_oriseg[cohort_oriseg[,1] == report$sample_id[report$sample_name == sample],]
            report <- shift.baseline(output_dir,cohort_oriseg,ind_seg,cohort_code,prior_prob_df,sample,report,bins_info,fit,ref_fit,genome,calling=FALSE,shift="lower",whole_shift = whole_shift,shift_id=whole_shift_ids,shiftnum=whole_shift_num,outputfilename=outputfilename,outputplotname=outputplotname)
        }
    }

    # cut-off using different thresholds
    if (length(cohort_noisy_samples) > 0){
        for (sample in cohort_noisy_samples){
            ind_seg <- cohort_oriseg[cohort_oriseg[,1] == report$sample_id[report$sample_name == sample],]
            report <- process.noise(output_dir,cohort_oriseg,ind_seg,cohort_code,prior_prob_df,sample,report,bins_info,fit,ref_fit,genome,calling=FALSE,whole_shift = whole_shift,shift_id=whole_shift_ids,shiftnum=whole_shift_num,outputfilename=outputfilename,outputplotname=outputplotname,lowthres=lowcutoff,highthres=highcutoff)
        }
    }
}


write.table(report, file=file.path(output_dir,'data_quality_report.txt'),sep = '\t',quote = F,row.names = F)

