#!/usr/bin/env nextflow

nextflow.enable.dsl=2


log.info """\

         CNAdjust  
         ===================================
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
