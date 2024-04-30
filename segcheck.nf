/* 
 * label and evaluate segments
 * 
 */
process segcheck{    
    tag {"$series"}
    
    input:
    val  series
    path inputdir 
    val  idmapping
    val  genome
    val  outputfilename
    val  outputplotname
    val  idmappingfilename
    output:
    path "$series"

    script:       
    """
    segcheck.R -series $series -workdir $projectDir -inputdir $inputdir -idmapping $idmapping -genome $genome -outputfilename $outputfilename -outputplotname $outputplotname -idmappingfilename $idmappingfilename
    """
}

