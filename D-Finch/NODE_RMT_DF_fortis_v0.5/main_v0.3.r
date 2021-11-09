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
# 25-05-2021 - created version 0.3
#            - modularised code to work with any single time series file

##############
## INITIATE ##
##############

## goal: initiate the RMT algorithm

## load NODE functions
source("NODERMT_v0.25.r")

## visualise the time series
## load data
TS = read.table("TS.csv",sep=";",header=T)

## transform real variables
TS[,c("rain","small.seeds","large.seeds","fortis.N")] = log(TS[,c("rain","small.seeds","large.seeds","fortis.N")])

## subset and order data
TS = TS[,c("year","rain","small.seeds","large.seeds","fortis.N","fortis.PC1.beak","fortis.PC2.beak")] 

# ## plot time series
# plot(TS[,1],TS[,2],xlab="Time (days)",ylab="Density (normalised)")
# for(j in 1:6)
# {
# 	points(TS[,1],TS[,j+1])
# 	lines(TS[,1],TS[,j+1],col=rainbow(6)[j])
# }

## apply NODERMT anlaysis
for(p in seq(0.02,0.1,0.02))
{
	pdf(paste("result_",p,".pdf",sep=""))
	#
	print(p)
	outputList = fit.NODE(sd2_o=0.1,sd2_p=p,N_e=500) 
	save(outputList,file=paste("TS_",p,".RData",sep=""))
	#
	dev.off()
}

## plot interaction networks
pdf("DIN.pdf")
for(p in seq(0.02,0.1,0.02))
{
	load(paste("TS_",p,".RData",sep=""))
    EMat = outputList[[2]]
	CMat = outputList[[3]]
	# CMat = t(apply(CMat,1,function(x)x/sum(x)))
	CMat[1,] = CMat[1,]*0
	# CMat = CMat/max(CMat)
	# CMat = CMat*(CMat>0.1)
	.plot.DIN(EMat,CMat,colnames(TS)[-1])
}
dev.off()

#
###
