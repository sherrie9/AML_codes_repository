remove(list=ls())

suppressPackageStartupMessages({
  # library(flowWorkspaceData)
  library(flowWorkspace)
  library(ggdendro)
  library(scamp)
  library(ggplot2)
  library(cowplot)
  library(knitr)
  library(dplyr)
  library(tidyr)
  library(faust)
  library(grDevices)
})

setwd("/home/rstudio/")

source("~/code/helperFunc_S.R")
source("~/code/helperFunc.R")
source("~/code/IMPC_Functions.R")

{library('flowCore')
library('flowCut')
library('stringr')
library('Cairo')
library('flowBin')
library('flowFP')
library('flowDensity')
library('pracma')
library('fields')
library('e1071')
library('flowStats')
library('MASS')
library('KernSmooth')
library('gridExtra')
library('plyr')
library('doMC')
library('alphahull')
library(igraph)}


no_cores <- 10
registerDoMC(no_cores)

store.allFCS.AML.normal <- readRDS("~/results/Preprocessing/store.allFCS.AML.normal.rds") # preprocessed metadata
days <- c("Day22","last before 2nd ind","last before cons") 
tubes <- c("001","002","003","004","005")


st <- 3 #st is starting populations, 1=cd45-ssc-+, 2=cd45+ssc+, 3=cd45+ssc-, 4=cd45++ssc-
        #5 is singlets

#----------------other days--------------------------------------------------------
# for(t in 1:length(tubes)){
results <- foreach (t = 3:4, .combine = rbind, .maxcombine = length(tubes), .multicombine = T) %:%
# for(d in 1:length(days)){
 foreach (d = 3, .combine = rbind, .maxcombine = length(days), .multicombine = T) %dopar% {   
   #===========================load files==============================================
   load(file=paste0("~/results/GatingSet_tube",tubes[t],"_Normal.RData"))
   fs.cd45negsscnegpos.normal <- llply(as(cd45negsscnegpos.flowD,Class = "list"),function(fD){
     f <- fD@flow.frame
     f@exprs <- na.omit(f@exprs)
     return(f)})
   fs.cd45possscpos.normal <- llply(as(cd45possscpos.flowD,Class = "list"),function(fD){
     f <- fD@flow.frame
     f@exprs <- na.omit(f@exprs)
     return(f)})
   fs.cd45possscneg.normal <- llply(as(cd45possscneg.flowD,Class = "list"),function(fD){
     f <- fD@flow.frame
     f@exprs <- na.omit(f@exprs)
     return(f)})
   fs.cd45pospossscneg.normal <- llply(as(cd45pospossscneg.flowD,Class = "list"),function(fD){
     f <- fD@flow.frame
     f@exprs <- na.omit(f@exprs)
     return(f)})
   singlet.fs.normal <- singlets.fs
   lv2.normal <- lv2
   
   
   load(file = paste0("~/results/GatingSet_tube",tubes[t],"_",days[d], ".RData"))

    if(st==1){
      fs <- llply(as(cd45negsscnegpos.flowD,Class = "list"),function(fD){
        f <- fD@flow.frame
        f@exprs <- na.omit(f@exprs)
        return(f)})
      fs.normal <- fs.cd45negsscnegpos.normal
      SN <- "CD45-SSC-+"
    }
    if(st==2){
      fs <- llply(as(cd45possscpos.flowD,Class = "list"),function(fD){
        f <- fD@flow.frame
        f@exprs <- na.omit(f@exprs)
        return(f)})
      fs.normal <- fs.cd45possscpos.normal
      SN <- "CD45+SSC+"
    }
    if(st==3){
      fs <- llply(as(cd45possscneg.flowD,Class = "list"),function(fD){
        f <- fD@flow.frame
        f@exprs <- na.omit(f@exprs)
        return(f)})
      fs.normal <- fs.cd45possscneg.normal
      SN <- "CD45+SSC-"
    }
    if(st==4){
      fs <- llply(as(cd45pospossscneg.flowD,Class = "list"),function(fD){
        f <- fD@flow.frame
        f@exprs <- na.omit(f@exprs)
        return(f)})
      fs.normal <- fs.cd45pospossscneg.normal
      SN <- "CD45++SSC-"
    }
    if(st==5|d==1){
      fs <- singlets.fs
      SN <- "Singlets"
      fs.normal <- singlet.fs.normal
    }
  if(class(fs)=="flowSet"){
    fs <- llply(as(fs,Class="list"),function(f)return(f))
    fs.normal <- llply(as(fs.normal,Class="list"),function(f)return(f))
    
  }
  
    # size.relapse <- length(which(lv2$st=="R"))
    # size.nonrelaspe <- length(which(lv2$st == "0"))
    # set.seed(1234)
    # f.ind <- c(sample(which(lv2$st=="R"), min(size.relapse,size.nonrelaspe)),
    #           sample(which(lv2$st=="0"), min(size.relapse,size.nonrelaspe)))
    # 
    # fs <- fs[f.ind]
    # lv2 <- lv2[f.ind,]
    # if(days[d] != "Day22"){
    #   set.seed(1234)
    #   normal.ind <- sample(1:nrow(lv2.normal), min(size.relapse,size.nonrelaspe))
    #   fs.normal <- fs.normal[normal.ind]
    #   lv2 <- rbind(lv2,lv2.normal[normal.ind,])
    # }
    
    # lv2 <- rbind(lv2,lv2.normal[normal.ind])  
  
    #---------------------------------Faust----------------------------------------------------
   
    names(fs.normal) <- paste0("V",c(41:60))
    # if(t==1){
      fs <- flowSet(c(fs,fs.normal))
      gs <- GatingSet(fs)
    # }else{
      # fs <- flowSet(fs)
      # gs <- GatingSet(fs)
    # }
    
    # fs <- flowSet(fs)
  
    
    
    #########################  FAUST ANALYSIS   ########################################
    # flowWorkspaceDataInfo()
    # gsPath <- system.file("extdata","gs_bcell_auto", package = "flowWorkspaceData")
    # gs_test <- load_gs(gsPath)
    
    startingNode <- "root"
    scat.chans <- c(grep(colnames(fs),pattern = "FSC*"), grep(colnames(fs),pattern = "SSC*"))
    names(scat.chans) <- colnames(fs)[scat.chans]
    
    activeChannelsIn <- markernames(gs)
    
    if(class(activeChannelsIn)=="list" && length(activeChannelsIn) >= 2 &&
       t==4|t==5){
      # fsApply(fs,function(f){return(na.omit(f@parameters@data$desc))})
      if(t==4){
        fs <- fsApply(fs,function(f){
          f@parameters@data$desc[which(f@parameters@data$desc == "CD38H7")] <- "CD38"
          # f@parameters@data$desc[which(f@parameters@data$desc == "CD11a")] <- "CD11A"
          return(f)
        }) 
      }
        
      if(t==5){
        fs <- fsApply(fs,function(f){
          f@parameters@data$desc[which(f@parameters@data$desc == "CD11a")] <- "CD11A"
          return(f)
        }) 
      }
      
      gs <- GatingSet(fs)
      activeChannelsIn <- markernames(gs)
    }
    # activeChannelsIn <- c(colnames(fs)[scat.chans],markernames(gs))
    
    # channelBoundsIn <- matrix(0,nrow=2,ncol=length(activeChannelsIn))
    # colnames(channelBoundsIn) <- activeChannelsIn
    # rownames(channelBoundsIn) <- c("Low","High")
    # channelBoundsIn["High",] <- 3.3
    # channelBoundsInpaste0(projPath,"/Sample_",timepoints[t],startingNode,".pdf")
    
    projPath <- paste0("~/results/Faust2/Faust-",SN,"_",days[d],"_tube",tubes[t])
    suppressWarnings(dir.create(projPath, recursive = TRUE))
    
    faust(
      gatingSet = gs,
      experimentalUnit = "name",
      activeChannels = activeChannelsIn,
      # channelBounds = channelBoundsIn,
      startingCellPop = startingNode,
      projectPath = projPath,
      depthScoreThreshold = 0.05,
      selectionQuantile = 1.0,
      debugFlag = FALSE,
      #set this to the number of threads you want to use on your system
      threadNum = parallel::detectCores() - 3 ,
      nameOccuranceNum=1,
      seedValue = 271828,
      annotationsApproved = FALSE # set to false before we inspect the scores plots.
    )
  
    path <- file.path(projPath,"faustData","plotData","scoreLines.pdf")
    cat(paste0("<img src='",path,"'>"))
    
    root_selectedChannels <- readRDS( file.path(projPath,"faustData","gateData","root_selectedChannels.rds"))
    start.faust <- Sys.time()
    faust(
      gatingSet = gs,
      experimentalUnit = "name",
      activeChannels = root_selectedChannels,
      # channelBounds = channelBoundsIn,
      startingCellPop = startingNode,
      projectPath = projPath,
      depthScoreThreshold = 0.05,
      selectionQuantile = 1.0,
      debugFlag = FALSE,
      #set this to the number of threads you want to use on your system
      threadNum = parallel::detectCores()- 3,
      nameOccuranceNum=1,
      seedValue = 271828,
      annotationsApproved = TRUE
    )
    time.faust <- TimeOutput(start.faust)
  
    
   
    write.table(time.faust, file = "~/results/Faust/Faust-runtime.txt",append = TRUE)
    
    rm(gs)
     

  
}
  


















