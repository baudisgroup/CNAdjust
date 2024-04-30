#!/usr/bin/env nextflow

nextflow.enable.dsl=2


log.info """\

         CNV CALIBRATION PIPELINE    
         ===================================
         series:           $params.series
         inputdir:         $params.inputdir
         outputdir:        $params.outputdir
         idmapping:        $params.idmapping
         usetumorref:      $params.usetumorref
         genome:           $params.genome 
         outputfile:       $params.outputfilename
         outputplot:       $params.outputplotname
         idmappingfile:    $params.idmappingfilename
         cohortassignfile: $params.tumorcodefilename
         cohortpriorfile:  $params.priorfilename
         """
         .stripIndent()  


 include {segcheck} from './modules/segcheck.nf'
 include {segcalibration} from './modules/segcalibration.nf'


workflow {
   series = Channel.of(params.series.split(',')) 

   segcheck(series,params.inputdir,params.idmapping,params.genome,params.outputfilename,params.outputplotname,params.idmappingfilename)
      
   segcalibration(segcheck.out,params.inputdir,params.usetumorref,params.genome,params.outputfilename,params.outputplotname,params.tumorcodefilename,params.priorfilename,
   params.logrsd,params.lowfrac,params.dupdelratio,params.lowfrac2,params.dupdelratio2,params.lowfrac3,params.segnum,params.lowthre,params.highthre)          
}
