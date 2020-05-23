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
  library(xtable)
  library(flowType)
  library(RchyOptimyx)
})
require(emmeans)
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
  library(igraph)
  # library("VennDiagram")
}
no_cores <- 10
registerDoMC(no_cores)


store.allFCS.AML.normal <- readRDS("~/results/Preprocessing/store.allFCS.AML.normal.rds") # preprocessed metadata
days <- c("Day22","last before 2nd ind","last before cons") 
tubes <- c("001","002","003","004","005")
st <- 3 #3is CD45+SSC- 2 IS CD45+SSC+

suppressWarnings(dir.create("~/results/flowType/", recursive = TRUE))
suppressWarnings(dir.create("~/results/flowType/CD45+SSC-/", recursive = TRUE))
suppressWarnings(dir.create("~/results/flowType/CD45+SSC+/", recursive = TRUE))

results <- foreach (t = c(2:5), .combine = rbind, .maxcombine = length(tubes), .multicombine = T) %dopar% { 
  # for(d in 1:length(days)){
  # foreach (d = c(2), .combine = rbind, .maxcombine = length(days), .multicombine = T) %dopar% { 
    # if(d==1){
    #   SN <- "Singlets"
    # }else{
    #   if(st == 3){
    #     SN <- "CD45+SSC-"
    #   }
    #   if(st == 2){
    #     SN <- "CD45+SSC+"
    #   }
    # }
  
    
  load(file = paste0("~/results/GatingSet_tube",tubes[t],"_",days[2], ".RData"))
  fs <- singlets.fs
  markers_list <- as.character(fs[[1]]@parameters@data[1:ncol(fs[[1]]), 2])
  markers_list <- Find.markers(fs[[1]],na.omit(markers_list))
  scat.chans <- c(grep(colnames(fs),pattern = "FSC*"), grep(colnames(fs),pattern = "SSC*"))
  names(scat.chans) <- colnames(fs)[scat.chans]
  PropMarkers <- c(scat.chans[1],markers_list[-which(names(markers_list) == "CD45")]) 
  MarkerNames <- names(PropMarkers)
  
  #load day22 flowType results
  SN <- "CD45+SSC-"
  projPath_fT <- paste0("~/results/flowType/",SN,"/tube",tubes[t])
  load(file = paste0(projPath_fT,"/Tresholds_",days[1],".RData"))
  load(file = paste0(projPath_fT,"/flowTypeResults_",days[1],".RData"))
  
  tresholds.list.day22 <- tresholds.list
  ResList.day22 <- ResList
  all.proportions.day22 <- matrix(0,length(ResList[[1]]@CellFreqs),60)
  for (i in 1:length(ResList)){
    all.proportions.day22[,i] = ResList[[i]]@CellFreqs / ResList[[i]]@CellFreqs[1]
  }
  rownames(all.proportions.day22) <- unlist(lapply(ResList[[1]]@PhenoCodes,
                                             function(x){return(
                                               flowType::decodePhenotype(x,
                                                                         MarkerNames,
                                                                         ResList[[1]]@PartitionsPerMarker))}))
  #load last before 2nd ind flowType results
  SN <- "CD45+SSC-"
  projPath_fT <- paste0("~/results/flowType/",SN,"/tube",tubes[t])
  load(file = paste0(projPath_fT,"/Tresholds_",days[2],".RData"))
  load(file = paste0(projPath_fT,"/flowTypeResults_",days[2],".RData"))
  
  tresholds.list.bef2ndind <- tresholds.list
  ResList.bef2ndind <- ResList
  all.proportions.bef2ndind <- matrix(0,length(ResList[[1]]@CellFreqs),60)
  for (i in 1:length(ResList)){
    all.proportions.bef2ndind[,i] = ResList[[i]]@CellFreqs / ResList[[i]]@CellFreqs[1]
  }
  rownames(all.proportions.bef2ndind) <- unlist(lapply(ResList[[1]]@PhenoCodes,
                                                   function(x){return(
                                                     flowType::decodePhenotype(x,
                                                                               MarkerNames,
                                                                               ResList[[1]]@PartitionsPerMarker))}))
  
  #load last before cons results
  projPath_fT <- paste0("~/results/flowType/",SN,"/tube",tubes[t])
  suppressWarnings(dir.create(paste0(projPath_fT,"/TrendPlots"), recursive = TRUE))
  load(file = paste0(projPath_fT,"/Tresholds_",days[3],".RData"))
  load(file = paste0(projPath_fT,"/flowTypeResults_",days[3],".RData"))
  
  tresholds.list.befcons <- tresholds.list
  ResList.befcons<- ResList
  all.proportions.befcons <- matrix(0,length(ResList[[1]]@CellFreqs),60)
  
  for (i in 1:length(ResList)){
    all.proportions.befcons[,i] = ResList[[i]]@CellFreqs / ResList[[i]]@CellFreqs[1]
  }
  rownames(all.proportions.befcons) <- unlist(lapply(ResList[[1]]@PhenoCodes,
                                                       function(x){return(
                                                         flowType::decodePhenotype(x,
                                                                                   MarkerNames,
                                                                                   ResList[[1]]@PartitionsPerMarker))}))
  
  #combine all proportions for all time points. only with relapsed vs non relapsed groups
  all.proportions <- cbind(all.proportions.day22[,c(1:40)],
                           all.proportions.bef2ndind[,c(1:40)],
                           all.proportions.befcons[,c(1:40)])
  all.proportions.relapsed <- cbind(all.proportions.day22[,c(1:20)],
                                    all.proportions.bef2ndind[,c(1:20)],
                                    all.proportions.befcons[,c(1:20)])
  all.proportions.nonrelapsed <- cbind(all.proportions.day22[,c(21:40)],
                                    all.proportions.bef2ndind[,c(21:40)],
                                    all.proportions.befcons[,c(21:40)])
   all.proportions.relapsed <- all.proportions.relapsed[-1,]
   all.proportions.nonrelapsed <- all.proportions.relapsed[-1,]
   
   all.proportions.relapsed <- data.frame(all.proportions.relapsed)
   all.proportions.nonrelapsed <- data.frame(all.proportions.nonrelapsed)
   
   p.vals <- NULL
   for(n in 1:nrow(all.proportions)){
     
   
   data.row <- data.frame(
     Phenotype=rep(rownames(all.proportions)[n], ncol(all.proportions)),
     Value = all.proportions[n,],
     Days = c(rep(1,40),rep(2,40),rep(3,40)),
     AML_Status =  c(rep("R",20),rep("0",20),rep("R",20),rep("0",20),
                       rep("R",20),rep("0",20)))

 

  
   fit  <- data.row %>% group_by(AML_Status) %>% do(model =lm(Value ~ Days, data=.))
   
   fit <- lm(Value ~ Days * AML_Status, data=data.row)
   # anova(fit)
 
   slopes <- emtrends(fit, 'AML_Status', var = 'Days')
   p.value <- summary(pairs(slopes))[6]
   
   p.vals <- c(p.vals,p.value)
   
   }
   
   #select pvals < 0.05
   
   selected <- unname(which(p.vals < 0.01))
   print(selected)
   p <- list()
   
   MyTable=cbind(rownames(all.proportions)[selected],
                 format(p.vals[selected], digits=2),
                 format(p.adjust(p.vals)[selected],digits=3),
                 format(rowMeans(all.proportions[selected,c(1:20)]),
                        digits=3),
                 format(rowMeans(all.proportions[selected,c(21:40)]),
                        digits=3),
                 format(rowMeans(all.proportions[selected,c(41:60)]),
                        digits=3))
   # MyTable <- MyTable[-which(rowMeans(all.proportions[selected,]) < 0.001),]
   colnames(MyTable)=c('Phenotype', 'p-value','adjusted p-value','mean relapse cell frequency',
                       'mean nonrelapse cell frequency','mean normal cell frequency')
   
   
   write.table(MyTable,file=paste0(projPath_fT,"/Trends_",tubes[t],
                                   ".csv"),row.names = FALSE,
               sep=",")
   
   #plotting regeneration linear models
   for(k in 1:length(selected)){
     plot.data <- data.frame(
       Phenotype=rep(rownames(all.proportions)[selected[k]], ncol(all.proportions)),
       Value = all.proportions[selected[k],],
       Days = c(rep(1,40),rep(2,40),rep(3,40)),
       AML_Status =  c(rep("R",20),rep("0",20),rep("R",20),rep("0",20),
                       rep("R",20),rep("0",20)))
     p[[k]] <- ggplot(plot.data, aes(x=Days,y=Value, color=AML_Status, fill=AML_Status)) + 
       geom_smooth(method="lm")+
       geom_point() + theme_bw() + ylab(paste0(unique(plot.data$Phenotype)," Proportions - % of singlets"))
     ggsave(filename = paste0(projPath_fT,"/TrendPlots/trends_",unique(plot.data$Phenotype),".png"),p[[k]])
    
   }
   
   
   
 
  
  
 
  

}
