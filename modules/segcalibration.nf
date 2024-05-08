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
    val   cohortfile
    val   priorfile
    val   useregion
    val   regionfile
    val   genome
    val   outputfile
    val   outputplot
    val   logrsd
    val   cnafrac
    val   dupdelratio
    val   cnafrac2
    val   dupdelratio2
    val   cnafrac3
    val   segnum
    val   lowthre
    val   highthre
    output:
    path outputdir

    script:       
    """
    segcalibration.R -workdir $projectDir -outputdir $outputdir -inputdir $inputdir -usetumorref $usetumorref -cohortfile $cohortfile -priorfile $priorfile -useregion $useregion -regionfile $regionfile -genome $genome -outputfile $outputfile -outputplot $outputplot  -logrsd $logrsd -cnafrac $cnafrac -dupdelratio $dupdelratio -cnafrac2 $cnafrac2 -dupdelratio2 $dupdelratio2 -cnafrac3 $cnafrac3  -segnum $segnum -lowthre $lowthre -highthre $highthre
    """
}