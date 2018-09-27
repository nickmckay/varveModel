load("output0106.Rdata")
library(here)
library(tidyverse)
#trim the rows of all NAs
good= which(!apply(is.na(compMat),1,all))
compMat=compMat[good,]
cMs=t(apply(compMat,1,sort,na.last=TRUE))

corenames = c("2","3","4","20","21")

#which iterations have been run:
goodCol= which(!apply(is.na(compMat),2,all))

nSim = max(goodCol)

medianComp = cMs[,round(nSim*.5)]
c975= cMs[,floor(nSim*.975)]
c025= cMs[,ceiling(nSim*.025)]

plot(medianComp,type="l")
lines(c975)
lines(c025)

#plot(compMat[,250],type="l")
ageModels = vector(mode="list",length=5)
origCores = cores
Depth.22BP = c(601.90	,621.14	,511.00,	874.71	,874.84)


#extrapolate age models?
corePost[[1]]$thicks




for(c in 1:5){
  #plot age vs depth:
  
  #look at cumsums
  good= which(!apply(is.na(corePost[[c]]$thicks),1,all))
  trimmed = corePost[[c]]$thicks[good,]
  trimmed.filled = trimmed
  ###IDEA? Populate missing areas in trimmed with other iterations so that we have more iteations for cumsumming #NO BAD IDEA
  for(nc in 1:ncol(trimmed)){
    thisIt =trimmed[,nc]
    newSec = TRUE
    for(ri in 1:nrow(trimmed)){
      if(is.na(thisIt[ri])){
        if(newSec){#if it's a new section
          #sample a column that has data for that row
          sampCol = sample(which(!is.na(trimmed[ri,])),1)
        }
        sampled = trimmed[ri,sampCol]
        if(is.na(sampled)){
          sampCol = sample(which(!is.na(trimmed[ri,])),1)
        }
        trimmed.filled[ri,nc]=sampled
      }else{
        next
      }
    }
  }
  
  
  csMat = (apply(trimmed,2,cumsum) + Depth.22BP[c])/1000 #account for missing top
  
  #remove columns of all NAs
  goodCol = which(!apply(is.na(csMat),2,all))
  csMat = csMat[,goodCol]
  
  mCs=apply(csMat,1,quantile,.5,na.rm=TRUE)
  cs975=apply(csMat,1,quantile,.975,na.rm=TRUE)
  cs025=apply(csMat,1,quantile,.025,na.rm=TRUE)
  year = seq_along(mCs) + 21
  
  
  
  orig.thick = (cumsum(cores[[c]]$thickness) + Depth.22BP[c])/1000
  orig.year = seq_along(orig.thick)
  
  if(c < 5){
    #load in and prepare original counts and uncertainty
    countUncertainty = na.omit(read.csv(here("Eklutna","Input0106",paste0("replicateCounting",corenames[c],".csv")),header = FALSE,strip.white = T)[,1:3])
    

    unc.hi = matrix(nrow = length(orig.thick))
    unc.lo = unc.hi
    for(d in 1:nrow(countUncertainty)){
      q = which(orig.thick>(countUncertainty[d,1]/100) & orig.thick<=(countUncertainty[d,2]/100))
      # unc.hi[q] = countUncertainty[d,1]+seq_along(q)*(1+abs(countUncertainty[d,3])/100)
      # unc.lo[q] = countUncertainty[d,1]+seq_along(q)*(1-abs(countUncertainty[d,3])/100)
      unc.hi[q] = orig.year[q]*(1+abs(countUncertainty[d,3])/100)
      unc.lo[q] = orig.year[q]*(1-abs(countUncertainty[d,3])/100)
    }
    
    orig.year = orig.year+21
    unc.hi = unc.hi+21
    unc.lo = unc.lo +21
    
    
    #remove reversals in unc
    na.lo = which(is.na(unc.lo))
    unc.lo = unc.lo[-na.lo]
    year.lo = orig.thick[-na.lo]
    # while(any(diff(unc.lo)<0)){
    #   unc.lo[which(diff(unc.lo)<0)] = NA
    #   na.lo = which(is.na(unc.lo))
    #   unc.lo = unc.lo[-na.lo]
    #   year.lo = year.lo[-na.lo]
    # }
    # 
    # k = 101
    # unc.lo.min = caTools::runmin(unc.lo,k = k,endrule = "NA")
    
    #remove reversals in unc
    na.hi = which(is.na(unc.hi))
    unc.hi = unc.hi[-na.hi]
    year.hi = orig.thick[-na.hi]
    
    #make better
    unc.hi = unc.hi[order(unc.hi)]
    unc.lo = unc.lo[order(unc.lo)]
    
    meandiff = rowMeans(cbind(abs(unc.hi - orig.year[-na.hi]),abs(unc.lo - orig.year[-na.lo])))
    unc.hi = orig.year[-na.hi]+meandiff
    unc.lo = orig.year[-na.lo]-meandiff
    
    # 
    unc.hi = unc.hi[order(unc.hi)]
    unc.lo = unc.lo[order(unc.lo)]              
    
    
    # unc.hi.max = caTools::runmax(unc.hi,k = k,endrule = "NA")
    
    # while(any(diff(unc.hi)<0)){
    #   unc.hi[which(diff(unc.hi)<0)] = NA
    #   na.hi = which(is.na(unc.hi))
    #   unc.hi = unc.hi[-na.hi]
    #   year.hi = year.hi[-na.hi]
    # }
    # 
    
  }else{
    orig.year = orig.year+21
  }
  
  na.with.year= apply(is.na(csMat),1,sum)/ncol(csMat)
  #year.to.stop = max(which(na.with.year<.5))
  year.to.stop = max(orig.year)
  library(ggplot2)
  age.depth.plot = ggplot()+
    geom_ribbon(aes(x = year,ymin = cs025,ymax = cs975,colour = "95% HDR",fill = "95% HDR"))+
    geom_line(aes(x = year,y=mCs,colour = "Median modeled age"))+
    
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  #add in y-scale
  if(c <= 3){
    age.depth.plot = age.depth.plot +    scale_y_reverse(name = "Depth (mblf)",limits = c(16,0),expand = c(0,0))+scale_x_continuous(name = "Age (yr BP)",expand = c(0,0),limits = c(0,2400))
  }else{
    age.depth.plot = age.depth.plot +    scale_y_reverse(name = "Depth (mblf)",limits = c(10,0),expand = c(0,0),breaks = seq(10,0,by=-2))+scale_x_continuous(name = "Age (yr BP)",expand = c(0,0),limits = c(0,950))
  }
  
  if(c < 5){
    #add in original counts
    age.depth.plot = age.depth.plot+
      geom_line(aes(x = orig.year,y=orig.thick,colour = "Original Counts"))+
      geom_line(aes(x = unc.hi,y=year.hi,colour = "Repeat Counter Uncertainty"),linetype = "dashed")+
      geom_line(aes(x = unc.lo,y=year.lo,colour = "Repeat Counter Uncertainty"),linetype = "dashed")+
      # geom_line(aes(x = unc.hi[order(unc.hi)],y=year.hi,colour = "Repeat Counter Uncertainty"),linetype = "dashed")+
      # geom_line(aes(x = unc.lo[order(unc.lo)],y=year.lo,colour = "Repeat Counter Uncertainty"),linetype = "dashed")+
      scale_colour_manual(name="Legend",
                          values=c("gray70","black","red","orange"))+
      scale_fill_manual(name = "Legend",values = c("gray70"))+
      ggtitle(paste("Core",corenames[c]))+
      theme(legend.position="none")
    
  }else{
    age.depth.plot = age.depth.plot+
      geom_line(aes(x = orig.year,y=orig.thick,colour = "Original Counts"))+
      scale_colour_manual(name="Legend",
                          values=c("gray70","black","red"))+
      scale_fill_manual(name = "Legend",values = c("gray70"))+
      ggtitle(paste("Core",corenames[c]))+
      theme(legend.position="none")
  }
  
  #for core 3, add in 14C age distribution
  #scaler (what percent of range for distribution?)
  scaler = 2
  if(c %in% c(1,2,4)){
    setwd("~/Dropbox/vplR/clam/")
    source("clam.R")
    
    if(c==2){
      depths = c(242.3,466.7,723.2,709.1,918.2)/100
      ages = c(98,  346,  1157,  1154,  1559)
      age.unc = c(29,
                  20,
                  31,
                  44,
                  28)
    }
    
    
    if(c==1){
      depths = c(52.8)/100
      ages = c(236)
      age.unc = c(30)
    }
    if(c==4){
      depths = c(825.2)/100
      ages = c(932)
      age.unc = c(26)
    }
    
    
    plot.range = c(0,max(cs975,na.rm = TRUE))
    
    for(a in 1:length(ages)){
      print(a)
      calibrate(cage = ages[a],error = age.unc[a],storedat = T)
      
      #if 
      if(max(diff(dat$calib[,1]))>1){
        new1 = min(dat$calib[,1]):max(dat$calib[,1])
        new2 = approx(x = dat$calib[,1],y = dat$calib[,2],xout = new1,method = "constant")
        dat$calib = (cbind(new1,new2$y))
      }
      
      #turn probabilities into depth units
      plot.vals = (dat$calib[,2]/sum(dat$calib[,2]))*abs(diff(plot.range))*scaler+0.001
      
      ageData = data.frame(year = dat$calib[,1],pmax = depths[a]+plot.vals,pmin = depths[a]-plot.vals)
      ageData$pmin[ageData$pmin < 0] = 0
      
      age.depth.plot = age.depth.plot + 
        geom_line(data = ageData,aes(x = year,y = pmax))+
        geom_ribbon(data = ageData,aes(x = year,ymax = pmax, ymin =pmin) ,
                    colour = "purple",fill = "purple",alpha = 0.2)
      
      if(c == 1){#save for later
        age.depth.plot1 <- age.depth.plot
      }
      #load("~/Dropbox/VarveModel/output4.Rdata")
      cores = origCores
    }
  }
  
  print(age.depth.plot)
  ggsave(filename = paste0("~/Dropbox/VarveModel/EklutnaFigures/AgeDepth.Core",corenames[c],".pdf"),width = 5,height = 4,units = "in")  
  
  #age.depth.plot+geom_line(aes(x=ageModels[[c]]$ageMedian,y=depthSequence/100),colour="green")
  csMat.mm = csMat*1000
  #create depth -age ensembles
  depthSequence = seq_len(median(apply(csMat.mm,2,max,na.rm=T))) #vector of depth in mm
  func = function(x) approx(x,y=seq_len(length(x)), xout = depthSequence)$y
  ageModels[[c]]$ageEnsemble = apply(csMat.mm,2,func)+21
  ageModels[[c]]$depth = depthSequence
  
  if(c==1){#
    #extrapolate
    newDepths <- readr::read_csv(here("Eklutna","Input0106","depths_deep_turbidites_loc2.csv")) %>% 
      mutate(depth = `depth (in cm)`*10)
    nd <- seq(max(ageModels[[c]]$depth)+1,to=max(newDepths$depth)+1)
    mf <- function(x){Hmisc::approxExtrap(ageModels[[c]]$depth,x,nd)$y}
    
    mf2ind <- round(runif(1000)*(nrow(ageModels[[c]]$ageEnsemble)-length(nd)))
mf2 <- function(x){
  mf2ind <- round(runif(1)*(nrow(ageModels[[c]]$ageEnsemble)-(length(nd)+604)))+604
  fake <- x[mf2ind:(mf2ind+length(nd)-1)]
  fake <-  fake - x[mf2ind-1] + max(x,na.rm = T)
  return(fake)
}   

mf3 <- function(x){
  gx <- diff(na.omit(x))
  
  return(cumsum(sample(gx,size = length(nd),replace = T))+max(x,na.rm = T))
}

    newAges <- apply(ageModels[[c]]$ageEnsemble,MARGIN = 2,mf3)
    ageModels[[c]]$ageEnsemble <- rbind(ageModels[[c]]$ageEnsemble,newAges)
    ageModels[[c]]$depth <- c(ageModels[[c]]$depth,nd)
    }
  
  
  ageModels[[c]]$ageMedian = apply(ageModels[[c]]$ageEnsemble,1,median,na.rm=T)
  ageModels[[c]]$IQR = apply(ageModels[[c]]$ageEnsemble,1,IQR,na.rm=T)
  ageModels[[c]]$c975 = apply(ageModels[[c]]$ageEnsemble,1,quantile,na.rm=T,probs = 0.975)
  ageModels[[c]]$c025 = apply(ageModels[[c]]$ageEnsemble,1,quantile,na.rm=T,probs = 0.025)
  
  ageModels[[c]]$st.dev = apply(ageModels[[c]]$ageEnsemble,1,sd,na.rm=T)
  
  ageModelToExport = as.data.frame(cbind(ageModels[[c]]$depth,ageModels[[c]]$ageMedian,ageModels[[c]]$st.dev,ageModels[[c]]$IQR))
  
  names(ageModelToExport) = c("Depth (cm)","Median Age (yr BP {1950})","Age Standard Deviation","Age Inter-quartile Range")
  write.csv(ageModelToExport,file = paste0("~/Dropbox/VarveModel/AgeModelData-core",corenames[c],".csv"))

    #add extrapolated lines in
  if(c == 1){
    newDepth <- nd
    ndi <- which(ageModels[[c]]$depth >= min(nd))
    newMed <- ageModels[[c]]$ageMedian[ndi]
    new975 <- ageModels[[c]]$c975[ndi]
    new025 <- ageModels[[c]]$c025[ndi]
    
    age.depth.plot.extrap <- age.depth.plot1 + xlim(c(0,max(new975))) +
      geom_line(aes(x = newMed,y = newDepth/1000),colour = "blue",linetype = 2) +
      geom_line(aes(x = new975,y = newDepth/1000),colour = "blue",linetype = 3) +
      geom_line(aes(x = new025,y = newDepth/1000),colour = "blue",linetype = 3)
    ggsave(age.depth.plot.extrap,filename = here("ExtrapolatedAgeModel.pdf"))
  }
}




# 
# plot(mCs/1000,type="l",xlab="varve year", ylab="depth (m)")
# lines(cs975/1000)
# lines(cs025/1000)
# #lines(csMs[,1]/1000)
# #lines(csMs[,1000]/1000)
# 
# lines(cumsum(cores[[c]]$thickness)/1000,col="blue")
# 
# depth.offset = -89.2 #in cm
# depths = (c(242.3,466.7,723.2,709.1,918.2) + depth.offset)/100
# 
# ages = c(111,387,1074,1072,1467)+age.offset
# points(ages,depths)
# title(main="Eklutna - Core 3")###
# #legend("bottomright",inset=.05,c("model results","original counts","14C ages"),col=c("black","red","blue"))
# c=5
# #indiviual core ensembles
# plot((corePost[[c]]$thicks[,1]),type="l")
# for(i in 1:10){
#   
#   lines((corePost[[c]]$thicks[,i]))
# }



#create histograms by depth here....

#histogram by depth in mm
#hist(ageMat[1400,]) 




library(ggplot2)
for(ML in 1:length(allMarkerLayers)){
  #marker layer age distributions
  test = sapply(sapply(corePost,"[[","mlCum"),function(x){which(x==allMarkerLayers[ML],arr.ind = TRUE)[,1]})
  if(is.list(test)){
    mlDist = rowMeans(as.data.frame(test[sapply(test,length)==max(sapply(test,length))]),na.rm=TRUE)
  }else{
    mlDist = rowMeans(test,na.rm=TRUE)
  }
  
  
  if(length(mlDist)>1){
    mlDist = mlDist+19 #add twenty one for varve year
    ml.hist = ggplot()+
      geom_histogram(aes(x = mlDist,y=..density..),bins = 15,fill = "grey90",colour = "grey30")+
      scale_x_continuous(name = "Age (yr BP)")+
      scale_y_continuous(name = "Probability Density",limits = c(0,NA),expand=c(0,0))+
      theme_bw()+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())+
      ggtitle(paste("Marker Layer",as.character(allMarkerLayers[ML])))
    
    
    #add in individual cores estimates
    countYears = sapply(sapply(cores,"[[","markerLayers"),function(x){which(x==allMarkerLayers[ML])})
    
    if(is.list(countYears)){
      countYears = unlist(countYears)
      countYears[(length(countYears)+1):5] = NA
    }
    
    #jitter countYears slightly
    t = 0
    for(j in 1:5){
      while(any(na.omit(countYears[j]==countYears[-j]))){
        t = t+1
        if (t%%2==0){
          countYears[j] = countYears[j]-.1
        }else if(t%%3==0){
          countYears[j] = countYears[j]-.15
        }else{
          countYears[j] = countYears[j]+.1
        }
      }
    }
    
    # countYears = countYears + rnorm(n = 5,sd = .5)
    
    #and account for varve offset
    countYears = countYears + 19
    
    names(countYears) = c("2","3","4","20","21")
    
    ml.hist = ml.hist + geom_vline(aes(xintercept = countYears,colour = names(countYears)),size = 1,alpha = 0.8)+
      theme(legend.position="none")
    
    print(ml.hist)
    ggsave(ml.hist,units = "in",width = 2.5,height = 2.5,filename = paste("~/Dropbox/VarveModel/EklutnaFigures/MarkerLayerHistograms/",paste0("MarkerLayer",as.character(allMarkerLayers[ML]),".pdf")))
  }
}
# 
# c=2
# #Individual core cumulatives
# #trim the rows of all NAs
# good= which(!apply(is.na(corePost[[c]]$thicks),1,all))
# cs1= apply(t(corePost[[c]]$thicks[good,]),MARGIN=1,cumsum)
# cs1sd=apply(cs1,1,sd,na.rm=TRUE)
# plot(cs1[,c],type="l")
# #indiviual core ensembles
# plot(cs1[,1]+19,type="l")
# for(i in 2:10){
#   
#   lines(cs1[,i]+19)
# }



#Tephra plots
tephfile = readr::read_csv(here("Eklutna","Input0106","TephraDepths.csv"))
tephdata <- tephfile[,-1]
ageData =vector(mode = "list",length = nrow(tephdata))

for(c in 1:nCores){
  ageData[[c]] =vector(mode = "list",length = nCores)
  
  for(t in 1:length(ageData)){
    ageData[[c]][[t]] = NA
    if(!is.na(tephdata[t,c])){
      adjDepth = tephdata[t,c]*10
      if(adjDepth < nrow(ageModels[[c]]$ageEnsemble)){
    ageData[[c]][[t]] = as.vector(ageModels[[c]]$ageEnsemble[as.numeric(round(adjDepth)),])
      }
    }
  }
  
}
#names(ageData) = corenames
core.names = colnames(tephdata)
tephra.names = tephfile$`Tephra ID/ Core ID`

# meltedAges = melt(ageData,na.rm = TRUE)
# 
# for(n in 1:5){
# meltedAges[meltedAges[,2]==n,2] = core.names[n]
# }

tephPlotHollow =list()
tephPlot= list()
tephPlotStack =list()

for(t in 1:nrow(tephfile)){
  rm("tdf")
  for(c in 1:5){
    
    if(!all(is.na(ageData[[c]][[t]]))){
      if(exists("tdf")){
        tdf = rbind(tdf,data.frame(age = ageData[[c]][[t]],Core = core.names[c]))
      }else{
        tdf = data.frame(age = ageData[[c]][[t]],Core = core.names[c])
        
      }
    }
    
  }
  tephPlotHollow[[t]]=ggplot(tdf)+
    geom_density(aes(x = age,colour = Core),alpha = 0.5)+
    #geom_freqpoly(aes(x = age,fill = Core))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    scale_colour_brewer(palette = "Set1")+
    scale_y_continuous(expand = c(0,0))+
    ggtitle(paste("Tephra",tephra.names[t]))+
    xlab("Age (yr BP)")+
    ylab("Probability Density")
  print( tephPlotHollow[[t]])
  
  tephPlot[[t]]=ggplot(tdf)+
    geom_density(aes(x = age,fill = Core),colour = "black",alpha = 0.3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    scale_fill_brewer(palette = "Set1")+
    scale_y_continuous(expand = c(0,0))+
    ggtitle(paste("Tephra",tephra.names[t]))+
    xlab("Age (yr BP)")+
    ylab("Probability Density")
  print( tephPlot[[t]])
  
  tephPlotStack[[t]]=ggplot(tdf)+
    geom_density(aes(x = age,fill = Core),colour = "black",alpha = 1,position = "stack")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    scale_fill_brewer(palette = "Accent")+
    scale_y_continuous(expand = c(0,0))+
    ggtitle(paste("Tephra",tephra.names[t]))+
    xlab("Age (yr BP)")+
    ylab("Probability Density")
  print( tephPlotStack[[t]])
  
    
}

library(gridExtra)

tephraOut = grid.arrange(grobs = tephPlot,ncol = 1)
ggsave(plot = tephraOut,filename = here("EklutnaFigures","TephraPlot.pdf"),width = 5,height = 11, units = "in")


tephraOutHollow = grid.arrange(grobs = tephPlotHollow,ncol = 1)
ggsave(plot = tephraOutHollow,filename = here("EklutnaFigures","TephraPlot-hollow.pdf"),width = 5,height = 11, units = "in")


tephraOutStack = grid.arrange(grobs = tephPlotStack,ncol = 1)
ggsave(plot = tephraOutStack,filename = here("EklutnaFigures","TephraPlot-Stack.pdf"),width = 5,height = 11, units = "in")


#Nore's events

#Tephra plots
tephfile = readr::read_csv(here("Eklutna","Input0106","eventsNore.csv"))
tephdata <- tephfile[,-1]/10
ageData =vector(mode = "list",length = nrow(tephdata))

for(c in 1:nCores){
  ageData[[c]] =vector(mode = "list",length = nCores)
  
  for(t in 1:length(ageData)){
    ageData[[c]][[t]] = NA
    if(!is.na(tephdata[t,c])){
      adjDepth = tephdata[t,c]*10
      if(adjDepth < nrow(ageModels[[c]]$ageEnsemble)){
        ageData[[c]][[t]] = as.vector(ageModels[[c]]$ageEnsemble[as.numeric(round(adjDepth)),])
      }
    }
  }
  
}
#names(ageData) = corenames
core.names = colnames(tephdata)
tephra.names = tephfile$`EVENT NAME`


tephPlotHollow =list()
tephPlot= list()
tephPlotStack =list()

for(t in 1:nrow(tephfile)){
  rm("tdf")
  for(c in 1:5){
    
    if(!all(is.na(ageData[[c]][[t]]))){
      if(exists("tdf")){
        tdf = rbind(tdf,data.frame(age = ageData[[c]][[t]],Core = core.names[c]))
      }else{
        tdf = data.frame(age = ageData[[c]][[t]],Core = core.names[c])
        
      }
    }
    
  }
  if(!exists("tdf")){
   #no data
    next
  }
  tephPlotHollow[[t]]=ggplot(tdf)+
    geom_density(aes(x = age,colour = Core),alpha = 0.5)+
    #geom_freqpoly(aes(x = age,fill = Core))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    scale_colour_brewer(palette = "Set1")+
    scale_y_continuous(expand = c(0,0))+
    ggtitle(paste("Tephra",tephra.names[t]))+
    xlab("Age (yr BP)")+
    ylab("Probability Density")
  ggsave( tephPlotHollow[[t]],filename = here("EklutnaFigures","EventHistograms",paste0(tephra.names[t],"-histHollow.pdf")))
  
  tephPlot[[t]]=ggplot(tdf)+
    geom_density(aes(x = age,fill = Core),colour = "black",alpha = 0.3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    scale_fill_brewer(palette = "Set1")+
    scale_y_continuous(expand = c(0,0))+
    ggtitle(paste("Event",tephra.names[t]))+
    xlab("Age (yr BP)")+
    ylab("Probability Density")
  ggsave( tephPlot[[t]],filename = here("EklutnaFigures","EventHistograms",paste0(tephra.names[t],"-hist.pdf")))
  
  tephPlotStack[[t]]=ggplot(tdf)+
    geom_density(aes(x = age,fill = Core),colour = "black",alpha = 1,position = "stack")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    scale_fill_brewer(palette = "Accent")+
    scale_y_continuous(expand = c(0,0))+
    ggtitle(paste("Tephra",tephra.names[t]))+
    xlab("Age (yr BP)")+
    ylab("Probability Density")
  ggsave( tephPlotStack[[t]],filename = here("EklutnaFigures","EventHistograms",paste0(tephra.names[t],"-histStack.pdf")))
  
  
}

