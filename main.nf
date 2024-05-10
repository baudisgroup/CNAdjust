#!/usr/bin/env nextflow

nextflow.enable.dsl=2


def helpMessage() {
    log.info"""
    ============================================================
    CNAdjust 
    ============================================================
    Usage:
    
    The typical command for running the pipeline is as follows:
    nextflow run hangjiaz/CNAdjust --inputdir <path to input directory> --series <name of series to be analyzed> --outputdir <path to output directory>
    
    Mandatory parameters:
      --inputdir                  Path to the input data directory.
      --series                    Name of the series to be analyzed. It should be folder names of the subfolders under `inputdir` and concatenated with commas.

    Other important parameters (optional):
      --outputdir                 Path to the output data directory. By default, it's in the current folder with the folder name 'output'.
      --use_external_ref          Logical value to determine if external information is used to adjust CNA. Default is true.
      --use_idmapping             Logical value to determine if sample identidier mapping is applied given the sampleid mapping file. Default is false.
      --use_custom_region         Logical value to determine if a custom prior is applied given the genomic bin location file. Default is false.      
      --genome                    Reference genome assembly version used in segment data. Available options are "hg38", "hg19", and "hg18". Default is "hg38".


   Filename parameters (optional):
      --cohort_assign_file        Filename of the input cohort assignment file. Default is "cohort-assignment.txt".
      --prior_file                Filename of the input cohort CNA occurrence file. Default is "cohort-cna-pattern.txt".
      --idmapping_file            Filename of the input sample id mapping file. Default is "sampleid-mapping.txt".
      --region_file               Filename of the input genomic region file. Default is "cohort-cna-region.txt".
      --output_seg                Filename of the output segment data. Default is "result.seg.txt".
      --output_plot               Filename of the output plot of segment data. Default is "segments.pdf".
    
   Potential problematic profile identification parameters (optional):
      --logrsd                    Upper limit of weighted standard deviation of logR by the number of involved markers to identify potential problematic samples in the strict criteria. Default is 0.35.
      --cnafrac                   Upper limit of CNA fraction to identify potential problematic samples in the strict criteria. Default is 0.5.
      --dupdelratio               Upper limit of the fraction ratio and its reciprocal as lower limit to identify potential problematic samples in the strict criteria. Default is 3.
      --cnafrac2                  Upper limit of CNA fraction to identify potential problematic samples in the more relaxed criteria 2. Default is 0.2.
      --dupdelratio2              Upper limit of the fraction ratio and its reciprocal as lower limit to identify potential problematic samples in the more relaxed criteria 2. Default is 5.
      --cnafrac3                  Upper limit of CNA fraction to identify potential problematic samples in the more relaxed criteria 3. Default is 0.7.
      --segnum                    Upper limit of segment number to identify noisy profiles. Defalut is 1000.
      --lowthre                   Thresholds to call low-level CNA in noisy profiles. It should be a positive value used to call duplication, and its opposite is used to call deletion. It can be multiple values concatenated with commas. Default is 0.1,0.15,0.3.
      --highthre                  Thresholds to call high-level CNA in noisy profiles. It should be a positive value used to call duplication, and its opposite is used to call deletion. It can be multiple values concatenated with commas. Default is 1,1.5,2.

    """.stripIndent()
}


// show helper message
if (params.help){
    helpMessage()
    exit 0
}

// Print log info
log.info """\

         CNAdjust  
         ==================================================
         series:              $params.series
         inputdir:            $params.inputdir
         outputdir:           $params.outputdir
         use_external_ref:    $params.use_external_ref
         cohort_assign_file:  $params.cohort_assign_file
         prior_file:          $params.prior_file
         use_custom_region:   $params.use_custom_region
         region_file:         $params.region_file
         use_idmapping:       $params.use_idmapping
         idmapping_file:      $params.idmapping_file
         genome:              $params.genome 
         output_seg:          $params.output_seg
         output_plot:         $params.output_plot
         criteria1_logrSD:    $params.logrsd
         criteria1_cnaFrac:   $params.cnafrac
         criteria1_fracRatio: $params.dupdelratio
         criteria2_cnaFrac:   $params.cnafrac2
         criteria2_fracRatio: $params.dupdelratio2
         criteria3_cnaFrac:   $params.cnafrac3
         noisy_segNum:        $params.segnum
         noisy_lowCutoff:     $params.lowthre 
         noisy_highCutoff:    $params.highthre
         """
         .stripIndent()  


 include {segcheck} from './modules/segcheck.nf'
 include {segcalibration} from './modules/segcalibration.nf'


workflow {
   series = Channel.of(params.series.split(',')) 

   segcheck(series,params.inputdir,params.use_idmapping,params.genome,params.output_seg,params.output_plot,params.idmapping_file)
      
   segcalibration(segcheck.out,params.inputdir,params.use_external_ref,params.cohort_assign_file,params.prior_file,params.use_custom_region,params.region_file,params.genome,
   params.output_seg,params.output_plot,params.logrsd,params.cnafrac,params.dupdelratio,params.cnafrac2,params.dupdelratio2,params.cnafrac3,params.segnum,params.lowthre,params.highthre)          
}
