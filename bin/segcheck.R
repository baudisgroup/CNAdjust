#! /usr/bin/env Rscript

args = R.utils::commandArgs(asValue = T)

series <- args$series
work_dir <- args$workdir
data_dir <- args$inputdir
idmapping <- as.logical(args$useidmapping)
genome <- args$genome
outputfilename <- args$outputfile
outputplotname <- args$outputplot
mappingfile <- args$idmappingfile

# source helper functions
source(file.path(work_dir,'bin','utils','general-utils.R'))
plotfiles_sources <- list.files(path = file.path(work_dir,'bin','utils','plotseg-utils'),recursive = T,full.names = T)
invisible(sapply(plotfiles_sources,source))
load(file.path(work_dir,'data','plot_seg_genome.rda'))
source(file.path(work_dir,'bin','plotseglabel.R'))

# check if input exist

if (!dir.exists(file.path(data_dir,series))) stop(series, " is not found in the input folder")

# create output dir
dir.create(path= file.path(series),showWarnings = F)
output_dir <- file.path(series)

# read segment data
data_file <- checkfile(filepath=file.path(data_dir,series),filetype="input segment data file", pattern="\\.seg")

data <- read.csv(data_file,header = T, sep='\t',colClasses = c("label"="character"))
## check data content
if (any(data[,5] <= 0)) stop("aberrant probe number in the input")
if (! "label" %in% colnames(data)) stop("label column is missing")
if (!all(data$label %in% c("+1","+2","0","-1","-2"))) stop ("CNA calling label is invalid, it should be the following values: '-1','+1','0','-2','+2'") 


# read sample id name mapping file
sample_ids <-  unique(data[,1])
if (idmapping){
  output_sample_names_file <- checkfile(filepath=file.path(data_dir,series),filetype="sample id mapping file",pattern=mappingfile)
  output_sample_names <- read.csv(output_sample_names_file,sep = '\t',header=T)
  output_sample_names <- output_sample_names[,2][match(sample_ids, output_sample_names[,1])]
} else{
  output_sample_names <- sample_ids
}

report <- list()
for (i in seq_len(length(sample_ids))){
  ind_sample_id <-  sample_ids[i]
  ind_sample_name <- output_sample_names[i]
  ind_sample_seg <- data[data[,1] == ind_sample_id,]
  dir.create(path= file.path(output_dir,ind_sample_name),showWarnings = F)
  # compute QC metrics
  ## LRR standard deviation of segments weighted by probe number
  lrrsegsd <- round(sd(rep(ind_sample_seg[,6],ind_sample_seg[,5])),3)
  ## segment number
  segNum <- dim(ind_sample_seg)[1]
  if (dim(ind_sample_seg)[1] > 0){
    ## called coverage
    cov <- seg.cov(ind_sample_seg)
    lowCov <- cov$lowCNA_cov
    highCov <-  cov$highCNA_cov
    normalCov <- cov$normal_cov
    lowCov_ratio <- cov$lowCNA_ratio
    
    # plot and save image of labelled segments
    plotseglabel(filepath = file.path(output_dir,ind_sample_name),filename=outputplotname,data = ind_sample_seg, assembly = genome)
    # write labelled segment file
    write.table(ind_sample_seg, file=file.path(output_dir,ind_sample_name,outputfilename), sep = '\t',row.names = F,quote=F)
    note <- "initial"
    
  } else{
    # if the segment data is failed to label
    lowCov <- 0
    highCov <- 0
    normalCov <- 0
    lowCov_ratio <- 0
    plotseglabel(filepath = file.path(output_dir,ind_sample_name),filename=outputplotname,data = ind_sample_seg,assembly = genome,no_label = T)
    write.table(ind_sample_seg, file=file.path(output_dir,ind_sample_name,outputfilename), sep = '\t',row.names = F,quote=F)
    note <- "failed-to-label"
  }
  
  report[[i]] <- data.frame(sample_id=ind_sample_id,sample_name=ind_sample_name, segment_num=segNum, LLR_segsd = lrrsegsd, 
                            normal_cov= normalCov, lowCNA_cov=lowCov, lowCNA_ratio=lowCov_ratio , highCNA_cov= highCov, note=note)
}
report <- do.call(rbind, report)
# write report 
write.table(report, file=file.path(output_dir,'analyse_report.txt'),sep = '\t',quote = F,row.names = F)

