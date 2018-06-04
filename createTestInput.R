#Create a synthetic example

#need a list for each core, that has an entry for each layer
#that includes,
#thickness
#OC%
#UC%
#marker layer ID

#example
source("~/Dropbox/VarveModel/syntheticVarves.R")
source("~/Dropbox/VarveModel/simulateOverAndUndercounting.R")
library(ggplot2)
nYr=2500
nML=10;
trueVarves=syntheticVarves(nYr,0.5)

OCP=.05#overcounting probability
UCP=.05#undercounting probability

while(TRUE){
markerLayers = sort(sample.int(nYr,size=nML))
if (min(diff(markerLayers))<100){break}
}

core1=trueVarves+syntheticVarves(2500,0.1)#add some varvey noise
core2=trueVarves+syntheticVarves(2500,0.1)#
core3=trueVarves+syntheticVarves(2500,0.1)

mCore1=simulateOverAndUndercounting(core1,OCP,UCP)

mCore2=simulateOverAndUndercounting(core2,OCP,UCP)

mCore3=simulateOverAndUndercounting(core3,OCP,UCP)

ML1=matrix(NA,nrow=length(markerLayers))
ML2=ML1
ML3=ML1

for(i in 1:length(markerLayers)){
  ML1[i]=which(abs(mCore1$old2new-markerLayers[i])==min(abs(mCore1$old2new-markerLayers[i])))
  ML2[i]=which(abs(mCore2$old2new-markerLayers[i])==min(abs(mCore2$old2new-markerLayers[i])))
  ML3[i]=which(abs(mCore3$old2new-markerLayers[i])==min(abs(mCore3$old2new-markerLayers[i])))
}
tp1=data.frame(x=1:length(mCore1$newThicks),c1=mCore1$newThicks)
tp2=data.frame(x=1:length(mCore2$newThicks),c2=mCore2$newThicks)
tp3=data.frame(x=1:length(mCore3$newThicks),c3=mCore3$newThicks)
class(tp1)
ggplot(data = tp1)+geom_line(aes(x=x,y=c1),colour="red") +geom_line(data=tp2,aes(x=x,y=c2),colour="blue")+geom_line(data=tp3,aes(x=x,y=c3))
#create marker layer column
ml1Col=matrix(NA,nrow=length(mCore1$newThicks))
ml1Prior =ml1Col
ucPrior1 =matrix(UCP,nrow=length(mCore1$newThicks))
ocPrior1 =matrix(OCP,nrow=length(mCore1$newThicks))
ml1Col[ML1]=paste0("ML",seq(1:length(ML1)))
ml1Prior[ML1]=0.9


ml2Col=matrix(NA,nrow=length(mCore2$newThicks))
ml2Prior =ml2Col
ucPrior2 =matrix(UCP,nrow=length(mCore2$newThicks))
ocPrior2 =matrix(OCP,nrow=length(mCore2$newThicks))
ml2Col[ML2]=paste0("ML",seq(1:length(ML2)))
ml2Prior[ML2]=0.9

ml3Col=matrix(NA,nrow=length(mCore3$newThicks))
ml3Prior =ml3Col
ucPrior3 =matrix(UCP,nrow=length(mCore3$newThicks))
ocPrior3 =matrix(OCP,nrow=length(mCore3$newThicks))
ml3Col[ML3]=paste0("ML",seq(1:length(ML3)))
ml3Prior[ML3]=0.9

#create "Truth columns"
OCTruth1=matrix(0,nrow=nYr)
OCTruth1[mCore1$OCi]=1;

UCTruth1=matrix(0,nrow=nYr)
UCTruth1[mCore1$UCi]=1;

OCTruth2=matrix(0,nrow=nYr)
OCTruth2[mCore2$OCi]=1;

UCTruth2=matrix(0,nrow=nYr)
UCTruth2[mCore2$UCi]=1;

OCTruth3=matrix(0,nrow=nYr)
OCTruth3[mCore3$OCi]=1;

UCTruth3=matrix(0,nrow=nYr)
UCTruth3[mCore3$UCi]=1;

OCTruth=data.frame(OCTruth1,OCTruth2,OCTruth3)
UCTruth=data.frame(UCTruth1,UCTruth2,UCTruth3)
thickTruth=data.frame(orig=trueVarves,core1,core2,core3)

#recover lost wasUCI data
test = simulateOverAndUndercounting()




truth=list(thickness=thickTruth,OCTruth=OCTruth,UCTruth=UCTruth)





#create input for model
core1 = data.frame(year = 1:length(mCore1$newThicks), thickness = mCore1$newThicks,markerLayers = ml1Col,mlPrior = ml1Prior,ocPrior = ocPrior1 , ucPrior=ucPrior1,OverCountedTruth=OCTruth1,UnderCountedTruth=UCTruth1)

core2 = data.frame(year = 1:length(mCore2$newThicks), thickness = mCore2$newThicks,markerLayers = ml2Col,mlPrior = ml2Prior,ocPrior = ocPrior2 , ucPrior=ucPrior2,OverCountedTruth=OCTruth2,UnderCountedTruth=UCTruth2)

core3 = data.frame(year = 1:length(mCore3$newThicks), thickness = mCore3$newThicks,markerLayers = ml3Col,mlPrior = ml3Prior,ocPrior = ocPrior3 , ucPrior=ucPrior3,OverCountedTruth=OCTruth3,UnderCountedTruth=UCTruth3)

cores = list(core1,core2,core3)
save(cores,file="~/Dropbox/VarveModel/testCores.Rdata")
save(truth,file="~/Dropbox/VarveModel/testCoresTruth.Rdata")

save