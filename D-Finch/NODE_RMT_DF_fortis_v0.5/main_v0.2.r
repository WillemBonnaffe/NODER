############
## main.r ##
############

## goal: load time series and apply the NODERMT analysis 

## update
# 10-03-2021 - created version 0.0
# 11-03-2021 - created version 0.1
#            - implemented figures for the paper 
# 29-03-2021 - created version 0.2
#            - combined figures for the paper
# 24-05-2021 - created version 0.3
#            - adapted the script to analyse the Darwin's Finch time series

##############
## INITIATE ##
##############

## goal: initiate the RMT algorithm

## visualise the time series
pdf("TS.pdf")
m = 1
par(mfrow=c(3,1))#,mar=c(4,4,2,2)+0.1)
	for(m in 1:1)
	{
		## load data
		TS = read.table(paste("TS_",m,".csv",sep=""),sep=";",header=T)
		
		## order data
		# TS = TS[,c("t","G","B","R")]
		
		## transform real variables
		# for(i in 2:ncol(TS)){TS[,i][which(TS[,i]<0.001)] = 0.005}
		# TS[,c("R","G","B")] = log(TS[,c("R","G","B")])
		
		## plot 
		plot(TS[,1],TS[,2],xlab="Time (days)",ylab="Density (normalised)")
		for(j in 1:6)
		{
			points(TS[,1],TS[,j+1])
			lines(TS[,1],TS[,j+1],col=rainbow(3)[j])
		}
	}
par(mfrow=c(1,1))#,mar=c(5,4,4,2)+0.1)
dev.off()

## run the NODERMT analysis for each time series
m = 1
for(m in 1:1)
{
	## load data
	TS = read.table(paste("TS_",m,".csv",sep=""),sep=";",header=T)
	
	## order data
	# TS = TS[,c("t","G","B","R")]
	
	## transform real variables
	for(i in 2:ncol(TS)){TS[,i][which(TS[,i]<0.001)] = 0.005}
	TS[,c("R","G","B")] = log(TS[,c("R","G","B")])
	
	## apply NODERMT anlaysis
	# pdf(paste("TS",m,".pdf",sep=""),height=18,width=6)
	# plot(TS)
	source("NODERMT_v0.20.r")
	save(outputList,file=paste("TS",m,".RData",sep=""))
	# dev.off()
}

#
###
