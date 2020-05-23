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


store.allFCS.AML.normal <- readRDS("~/results/Preprocessing/store.allFCS.AML.normal.rds") # read preprocessed metadata
days <- c("Day22","last before 2nd ind","last before cons") 
tubes <- c("001","002","003","004","005")
st <- 3 #3is CD45+SSC- 2 IS CD45+SSC+

suppressWarnings(dir.create("~/results/flowType/", recursive = TRUE))
suppressWarnings(dir.create("~/results/flowType/CD45+SSC-/", recursive = TRUE))
suppressWarnings(dir.create("~/results/flowType/CD45+SSC+/", recursive = TRUE))

results <- foreach (t = c(2:5), .combine = rbind, .maxcombine = length(tubes), .multicombine = T) %:%
    
  foreach (d = c(2), .combine = rbind, .maxcombine = length(days), .multicombine = T) %dopar% { 

      if(st == 3){
        SN <- "CD45+SSC-"
      }
      if(st == 2){
        SN <- "CD45+SSC+"
      }
   
    load(file = paste0("~/results/GatingSet_tube",tubes[t],"_",days[d], ".RData")) #load preprocessed data
 
      fs <- llply(as(cd45possscneg.flowD,Class = "list"),function(fD){
        f <- fD@flow.frame
        f@exprs <- na.omit(f@exprs)
        return(f)})
      fs.singlets.p <- llply(as(singlets.fs,Class = "list"),function(f){return(f)})
   

    
    singlets.counts <- fsApply(singlets.fs,function(f)return(nrow(f@exprs)))

    load(file=paste0("~/results/GatingSet_tube",tubes[t],"_Normal.RData")) #load preprocessed data
  
      if(st==3){
        fs.normal <- llply(as(cd45possscneg.flowD,Class = "list"),function(fD){
          f <- fD@flow.frame
          f@exprs <- na.omit(f@exprs)
          return(f)})
        fs.singlets.normal.p <- llply(as(singlets.fs,Class = "list"),function(f){return(f)})
       
      }
      if(st == 2){
        fs.normal <- llply(as(cd45possscpos.flowD,Class = "list"),function(fD){
          f <- fD@flow.frame
          f@exprs <- na.omit(f@exprs)
          return(f)})
        fs.singlets.normal.p <- llply(as(singlets.fs,Class = "list"),function(f){return(f)})
      }
     
    
    
    
      names(fs.normal) <- paste0("V",41:60)
      normals.counts <- fsApply(singlets.fs,function(f)return(nrow(f@exprs)))
      names(fs.singlets.normal.p) <- paste0("v",41:60)
      fs.singlets.p <- flowSet(c(fs.singlets.p,fs.singlets.normal.p))
      singlets.counts <- rbind(singlets.counts,normals.counts)
      lv2 <- data.frame(
        sample=c(paste0("V",1:60)),
        AML_status=c(rep("R",20),rep("0",20),rep("Normal",20))
                        )
    
  
      names(fs) <- paste0("V",1:40)
      # load(file = paste0("~/results/GatingSet_tube",tubes[t],"_",days[d], ".RData"))
      fs <- flowSet(c(fs,fs.normal))
    
      if(class(fs)=="list"){
        fs <- flowSet(fs)
      }
      markers_list <- as.character(fs[[1]]@parameters@data[1:ncol(fs[[1]]), 2])
      markers_list <- Find.markers(fs[[1]],na.omit(markers_list))
      scat.chans <- c(grep(colnames(fs),pattern = "FSC*"), grep(colnames(fs),pattern = "SSC*"))
      names(scat.chans) <- colnames(fs)[scat.chans]
      PropMarkers <- c(scat.chans[1],markers_list[-which(names(markers_list) == "CD45")]) 
      MarkerNames <- names(PropMarkers) 
    
      
      #------------------flowType analysis---------------------------------------------
      projPath_fT <- paste0("~/results/flowType/",SN,"/tube",tubes[t])
      suppressWarnings(dir.create(paste0(projPath_fT,"/gatingStrat"), recursive = TRUE))
     
      #finding automated gating thresholds for each marker
      tresholds.list <- lapply(1:length(fs),function(i){
        tresholds <- NULL
        #   ###############################
        #   #tube code insert
        #   ###############################
          for(n in colnames(fs)[PropMarkers]){
            f <- fs[[i]]
            tresholds[n] <- deGate(fs[[i]],n, percentile = 0.95)
            if(n == "FSC-A"){
              peaks <- getPeaks(fs[[i]],n)
                if(length(peaks$Peaks) == 4){
                  peaks$Peaks <- peaks$Peaks[-which(peaks$Peaks < 50000)]
                }
      
      
                if(length(peaks$Peaks)>=3){
                  allcuts <- deGate(f, channel = c("FSC-A"),all.cuts = T) #find gates at head
                  gt <- min(allcuts)
                  if(gt > 100000){
                    gt <- deGate(f, channel = c("FSC-A"), use.upper = T, upper=F)
                  }
          
      
                }else if(length(peaks$Peaks)==2){
      
                  allcuts <- deGate(f, channel = c("FSC-A"),all.cuts = T,use.percentile = FALSE,
                                    percentile = 0.99)
                  gt <- min(allcuts)
                  if(length(allcuts) == 2){
                    if(allcuts[1] < 50000 && allcuts[2] < 50200){
                      gt <- allcuts[2]
                    }
                  }
      
                  gt_perc <- deGate(f, channel = c("FSC-A"),all.cuts = T,use.percentile = TRUE,
                                    percentile = 0.99)
                  if(gt == gt_perc ){
                    gt <- deGate(f, channel = c("FSC-A"), use.upper = T, upper=F)
                  }
                  if(gt > 100000){
                    gt <- deGate(f, channel = c("FSC-A"), use.upper = T, upper=F) #find gates at head
                  }
      
                }else{
                  gt <- deGate(f, channel = c("FSC-A"), use.upper = T, upper=F) #find gates at head
      
                }
      
              tresholds[n] <- gt
            }
            if(n == "HLADR" && d==1 ){
              tresholds[n] <- 2
            }
          }
        return(tresholds)
      
      })
      ### running flowType analysis
      ResList <-
        lapply(1:length(fs),function(i){
          flowType::flowType(fs[[i]],PropMarkers,Methods ='thresholds',
                             Thresholds= as.list(tresholds.list[[i]]),
                             MarkerNames,
                             MaxMarkersPerPop = 4)
          })
      save(ResList,file = paste0(projPath_fT,"/flowTypeResults_",days[d],".RData"))
      load(file = paste0(projPath_fT,"/flowTypeResults_",days[d],".RData"))
      Proportions <-ResList[[1]]@CellFreqs
      Proportions <- Proportions / max(Proportions)
      # PropMarkers = c(3,6,7,8,9)
      names(Proportions) <- unlist(lapply(ResList[[1]]@PhenoCodes,
                                          function(x){return(flowType::decodePhenotype(
                                            x,MarkerNames, partitions.per.marker = rep(2,length(MarkerNames))
                                            # ,
                                            # ,
                                            # c(F,F,T,F,F,T,T,T,T,F)
                                            # ResList[[1]]@PartitionsPerMarker
                                          ))}))
      
      #----------------Pvalues from flowType-----------------------------------------------
      
      suppressWarnings(dir.create(projPath_fT, recursive = TRUE))
      all.proportions <- matrix(0,length(ResList[[1]]@CellFreqs),60)
      for (i in 1:length(ResList)){
        all.proportions[,i] = ResList[[i]]@CellFreqs / ResList[[i]]@CellFreqs[1]
      }
      rownames(all.proportions) <- unlist(lapply(ResList[[1]]@PhenoCodes,
                                                 function(x){return(
                                                   flowType::decodePhenotype(x,
                                                                   MarkerNames,
                                                                   ResList[[1]]@PartitionsPerMarker))}))
      Pvals <- vector()
      EffectSize <- vector()
      Labels <- lv2['AML_status']
      
      for (i in 1:dim(all.proportions)[1]){
      
        if (length(which(all.proportions[i,]!=1))==0){
          Pvals[i]=1
          EffectSize[i]=0
          next
        }
        temp = t.test(all.proportions[i,Labels=="R"],
                           all.proportions[i,Labels=="Normal"],
                           paired = FALSE, alternative = "two.sided")
      
        Pvals[i] <- temp$p.value
      
      
      }
      
      Pvals[is.nan(Pvals)]=1
      
      names(Pvals)=rownames(all.proportions)
      
      selected <- which(p.adjust(Pvals)<0.05)
      print(selected)
      
      Pvals2 <- vector()
      EffectSize2 <- vector()
      
      
      for (i in 1:dim(all.proportions)[1]){
      
        if (length(which(all.proportions[i,]!=1))==0){
          Pvals2[i]=1
          EffectSize2[i]=0
          next
        }
        temp = t.test(all.proportions[i,Labels=="0"],
                      all.proportions[i,Labels=="Normal"],
                      paired = FALSE, alternative = "two.sided")
      
        Pvals2[i] <- temp$p.value
      
      
      }
      
      Pvals2[is.nan(Pvals2)]=1
      
      names(Pvals2)= rownames(all.proportions)
      
      selected2 <- which(p.adjust(Pvals2)<0.05)
      print(selected2)
      
      #----------plotting venn diagram -----------------------------
      png(filename = paste0(projPath_fT,"/",days[d],"_Venn_Diagram.png"),
          width = 700, height = 500, units = "px", pointsize = 14)
      draw.pairwise.venn(area1=length(selected),area2 = length(selected2),
                         cross.area = length(intersect(selected,selected2)),
                         category = c("Relapsed Vs Normal","Nonrelapsed Vs Normal"))
      
      dev.off()
      
      selected_rel <-  setdiff(selected, intersect(selected,selected2))
      selected_non <-  setdiff(selected2, intersect(selected,selected2))
      lo <- 1
      for(lo in 1:2){ #iterate twice over both relapsed and non-relapsed
        if(lo==1){
          selected <- selected_rel
          Pvals<- Pvals
        }
        if(lo == 2){
          selected <- selected_non
          Pvals <- Pvals2
        }
      
      names(selected) <- rownames(all.proportions[selected,])
      print(selected)    
      
      if(length(selected) > 0){
        # write and save pvalue tables --------------------------------------------------------
        MyTable=cbind(rownames(all.proportions)[selected],
                      format(Pvals[selected], digits=2),
                      format(p.adjust(Pvals)[selected],digits=3),
                      format(rowMeans(all.proportions[selected,c(1:20)]),
                             digits=3),
                      format(rowMeans(all.proportions[selected,c(21:40)]),
                             digits=3),
                      format(rowMeans(all.proportions[selected,c(41:60)]),
                             digits=3))
        MyTable <- MyTable[-which(rowMeans(all.proportions[selected,]) < 0.001),]
        colnames(MyTable)=c('Phenotype', 'p-value','adjusted p-value','mean relapse cell frequency',
                            'mean nonrelapse cell frequency','mean normal cell frequency')
        
        
        write.table(MyTable,file=paste0(projPath_fT,"/CellPopulations2_",days[d],
                                        ifelse(lo==1,"withinRelapsed","withinNonrelapsed"),
                                        ".csv"),row.names = FALSE,
                    sep=",")
        
        
        MyTable <- MyTable[-which(rowMeans(all.proportions[selected,]) < 0.001),]
        plot.data <- data.frame(MyTable)
        
        plot.data <-  plot.data %>% gather(meanfrequency,value,c(4:6), factor_key = TRUE)
        plot.data$value<- round(as.numeric(plot.data$value), 3)
          p <- plot.data %>% ggplot(aes(fill=meanfrequency, x= Phenotype,y=value)) +
          geom_bar(position="dodge", stat="identity") +
          theme_bw() +
          theme(axis.text.x= element_text(angle=25,hjust=0.9, size=5),
                            plot.title = element_text(size=7),
                            plot.margin = unit(c(0.2,0.2,-0.5,0.1),"cm")) +
          xlab('Phenotype') + ylab('Cell proportion - % of singlets')
        ggsave(filename =paste0(projPath_fT,"/cellproportions_",days[d],".png"))

        phenotypes <- unname(MyTable[,"Phenotype"])
        suppressWarnings(dir.create(paste0(projPath_fT,"/RchyOp/"), recursive = TRUE))
        suppressWarnings(dir.create(paste0(projPath_fT,"/RchyOp/",days[d],"/"), recursive = TRUE))

        
       #----Plotting phenotypes --------------------------------------------------------------------------
        for(po in 1:length(phenotypes)){
          pheno <- unlist(str_split(phenotypes[po],"\\W+"))
          pheno <- pheno[pheno!=""]
          if(pheno[1]=="FSC"){
            pheno <- c("FSC-A",pheno[3:length(pheno)])
          }
          if(length(pheno) == 4){
        
        
        
            res.Rchy<-RchyOptimyx(pheno.codes=ResList[[1]]@PhenoCodes, phenotypeScores=-log10(Pvals),
                                  startPhenotype=ResList[[1]]@PhenoCodes[selected[po]], factorial(6), FALSE)
        
        
        
            old.par <- par(mar = c(0, 0, 0, 0))
            par(old.par)
            png(filename = paste0(projPath_fT,"/RchyOp/",days[d],"/",
                                  names(selected[po]),
                                  "_Rchy.png"),
                width = 700, height = 500, units = "px", pointsize = 12)
            plot(res.Rchy, phenotypeScores=-log10(Pvals), phenotypeCodes=ResList[[1]]@PhenoCodes,
                 marker.names=MarkerNames, ylab='-log10(Pvalue)')
            dev.off()
          }
        }
        
        
        if(length(phenotypes) > 10){
          phenotypes <- phenotypes[1:10]
        }

        for(p in 1:length(phenotypes)){
          pop_to_plot_name <- phenotypes[p]
          pheno <- unlist(str_split(phenotypes[p],"\\W+"))
          pheno <- pheno[pheno!=""]
          quat <- unlist(str_split(phenotypes[p],"\\w+"))
          quat <- quat[quat!=""]

          if(pheno[1]=="FSC"){
            pheno <- c("FSC-A",pheno[3:length(pheno)])
            quat <- quat[-1]
          }
          if(length(quat) > 2){
            maintitle <- paste0(pheno[1],quat[1],pheno[2],quat[2])
            if(paste0(quat[1:2],collapse = "") == "++"){
              fDposition <- c(T,T)
            }
            if(paste0(quat[1:2],collapse = "") == "+-"){
              fDposition <- c(T,F)
            }
            if(paste0(quat[1:2],collapse = "") == "-+"){
              fDposition <- c(F,T)
            }
            if(paste0(quat[1:2],collapse = "") == "--"){
              fDposition <- c(F,F)
            }
          }
          pop_to_plot <- PropMarkers[pheno]
          if(length(pop_to_plot)==1){
            pop_to_plot <- c(pop_to_plot,scat.chans["SSC-A"])
          }

          suppressWarnings(dir.create(paste0(projPath_fT,"/gatingStrat/",days[d],"/",pop_to_plot_name,"/"), recursive = TRUE))

          a <- lapply(1:length(fs.singlets.p),function(i){

            tresholds <- tresholds.list[[i]]
            chans <- c(markers_list["CD45"],scat.chans["SSC-A"])
            png(filename = paste0(projPath_fT,"/gatingStrat/",days[d],"/",pop_to_plot_name,"/Sample_V",i,".png"),
                width = 700, height = 700, units = "px", pointsize = 12)
            par(mfrow = c(2,2), mar = (c(5, 5, 4, 2) + 0.1))
            plotDens(fs.singlets.p[[i]],channels = chans, main="singlets")
            points(fs[[i]]@exprs[,chans],pch=".",col="red")
            plotDens(fs[[i]], channels = c(pop_to_plot[1],pop_to_plot[2]),main=SN,density.overlay = c(T,T))
            abline(v=tresholds[colnames(fs)[pop_to_plot][1]], lwd=2)
            abline(h=tresholds[colnames(fs)[pop_to_plot][2]], lwd=2)
            if(length(pop_to_plot) ==3){
              temp <- flowDensity(fs[[i]],channels =pop_to_plot[1:2],position = fDposition,
                                  gate=c(tresholds[colnames(fs)[pop_to_plot][1]],
                                         tresholds[colnames(fs)[pop_to_plot][2]]))
              if(temp@cell.count > 1 ){
                plotDens(temp,channels =c(pop_to_plot[3],scat.chans[3]),main=maintitle,
                         density.overlay = c(T,F))
                abline(v=tresholds[colnames(fs)[pop_to_plot[3]]], lwd=2)
              }

            }
            if(length(pop_to_plot) ==4){
              temp <- flowDensity(fs[[i]],channels =pop_to_plot[1:2],position = fDposition,
                                  gate=c(tresholds[colnames(fs)[pop_to_plot][1]],
                                         tresholds[colnames(fs)[pop_to_plot][2]]))
              if(temp@cell.count > 1){
                plotDens(temp,channels =c(pop_to_plot[3],pop_to_plot[4]),main=maintitle,
                         density.overlay = c(T,T))
                abline(v=tresholds[colnames(fs)[pop_to_plot[3]]],type="l", lwd=2)
                abline(h=tresholds[colnames(fs)[pop_to_plot[4]]],type="l", lwd=2)
              }
            }

            dev.off()
          })
        }
        
     
        
      }
        
  
      }
      
      




