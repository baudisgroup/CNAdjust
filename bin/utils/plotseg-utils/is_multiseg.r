#Function to check if segments come from pcf-routine or multipcf-routine


##Input:
### segments: a data frame with segmentation results from pcf, multipcf or aspcf

## Output:
### multi: logical value indicating whether segments comes from multipcf or not

##Required by:
## checkSegments
## getUnisegFormat
## selectSegments

##Requires:
## none


is.multiseg <- function(segments){
	#If not multisegment, the first column name should be "sampleID"
	multi <- ifelse(colnames(segments)[1]=="sampleID",FALSE,TRUE)
	
	return(multi)
}#end is.multiseg

