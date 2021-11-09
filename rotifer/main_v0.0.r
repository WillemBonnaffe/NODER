############
## main.r ##
############

## goal: load time series and apply the NODERMT analysis 

##############
## INITIATE ##
##############

## goal: initiate the RMT algorithm

## replicates
for(m in 1:3)
{
	## load data
	TS = read.table(paste("TS_",m,".csv",sep=""),sep=";",header=T)
	
	## order data
	# TS = TS[,c("t","G","B","R")]
	
	## transform real variables
	for(i in 2:ncol(TS)){TS[,i][which(TS[,i]<0.001)] = 0.005}
	TS[,c("R","G","B")] = log(TS[,c("R","G","B")])
	
	## apply NODERMT anlaysis
	pdf(paste("TS",m,".pdf",sep=""))
	plot(TS)
	source("NODERMT_v0.18.r")
	dev.off()
}

#
###

