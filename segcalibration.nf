/* 
 * calibrate CNV calling 
 * 
 */
process segcalibration{
    tag "$outputdir"

    publishDir params.outputdir, mode: 'copy'

    input:
    path  outputdir
    path  inputdir
    val   usetumorref
    val   genome
    val   outputfilename
    val   outputplotname
    val   tumorcodefilename
    val   priorfilename
    val   logrsd
    val   lowfrac
    val   dupdelratio
    val   lowfrac2
    val   dupdelratio2
    val   lowfrac3
    val   segnum
    val   lowthre
    val   highthre
    output:
    path outputdir

    script:       
    """
    segcalibration.R -workdir $projectDir -outputdir $outputdir -inputdir $inputdir -usetumorref $usetumorref -genome $genome -outputfilename $outputfilename -outputplotname $outputplotname -tumorcodefilename $tumorcodefilename -priorfilename $priorfilename -logrsd $logrsd -lowfrac $lowfrac -dupdelratio $dupdelratio -lowfrac2 $lowfrac2 -dupdelratio2 $dupdelratio2 -lowfrac3 $lowfrac3  -segnum $segnum -lowthre $lowthre -highthre $highthre
    """
}