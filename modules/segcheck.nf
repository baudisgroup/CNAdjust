/* 
 * label and evaluate segments
 * 
 */
process segcheck{    
    tag {"$series"}
    
    input:
    val  series
    path inputdir 
    val  useidmapping
    val  genome
    val  outputfile
    val  outputplot
    val  idmappingfile
    output:
    path "$series"

    script:       
    """
    segcheck.R -series $series -workdir $projectDir -inputdir $inputdir -useidmapping $useidmapping -genome $genome -outputfile $outputfile -outputplot $outputplot -idmappingfile $idmappingfile
    """
}

