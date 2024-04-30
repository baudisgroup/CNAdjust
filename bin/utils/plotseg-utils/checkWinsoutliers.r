# Function that check that input in winsoutliers has same dimension as data (in plots)

##Input:
### winsoutliers: data frame with outliers statueses
### data: data frame with logR data


##Required by:
### checkAndRetrievePlotInput


##Requires:
### pullOutContent

checkWinsoutliers <- function(winsoutliers,data){
	
	  winsoutliers <- pullOutContent(winsoutliers,what="wins.outliers")
		#Check that winsoutliers has same dimension as data:
		if(!all(dim(winsoutliers)==dim(data))){
			stop("winsoutliers must have the same number of rows and columns as data",call.=FALSE)
		}
		return(winsoutliers)
}#end checkWinsoutliers
