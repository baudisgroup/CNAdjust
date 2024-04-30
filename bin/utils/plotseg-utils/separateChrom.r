#Function that returns the index where each chromosome starts (and the last chromosome ends)

##Input:
### v: a vector of chromosome numbers

## Output:
### cp: indeces for start of each chromosome and end of last chromosome

##Required by:
### addChromlines
### adjustSeg

##Requires:
### none


separateChrom <- function(v){
	d <- diff(v)   #get difference between value (i+1) and value i in vector v
	cp <- which(d!=0)+1  #get changepoints
	
	#Add start of vector and (stop+1) of the whole vector
	cp <- c(1,cp,(length(v)+1))

	return(cp)
}#end separateChrom
