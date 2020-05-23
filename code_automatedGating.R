remove(list=ls())
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
library('Rtsne')
library("gplots")
library('fields')
library('e1071')
library('flowStats')
library('ggplot2')
library('MASS')
library('KernSmooth')
library('gridExtra')
library('plyr')
library('doMC')
library('alphahull')
library(igraph)}

no_cores <- 10
registerDoMC(no_cores)

suppressWarnings(dir.create("~/results/Plots/AutomatedGating"))
suppressWarnings(dir.create("~/results/Summary/AutomatedGating"))

skip_plot <- T
store.allFCS.AML.normal <- readRDS("~/results/Preprocessing/store.allFCS.AML.normal.rds") # preprocessed metadata

all.gthres <- list()
#---------------------tube 1 gating----------------------------------------------

tube1fcs <- store.allFCS.AML.normal[which(store.allFCS.AML.normal$Tube =="001"),]
# patients <- unique(tube1fcs$PatientID) #don't do patients[14] ==13201700010, transformation is wrong
scatter_chans <- c("FSC-A","FSC-H","SSC-A")
# tubes <- c("001","002","003","004","005","006","007","008")
days <- c("Day22","last before 2nd ind","last before cons", "Normal")

tubes <- c("001")
suppressWarnings(dir.create("~/results/Summary/AutomatedGating/", recursive = T))
for(d in 2:length(days)){
  
  
  
  print(as.character(days[d]))
  # if(days[x] == "Day22"){
  #   tube1fcs_1 <- tube1fcs[which(tube1fcs$Status == days[x]),]
  # }else{
  #   tube1fcs_1 <- tube1fcs[c(which(tube1fcs$Status == days[x]),which(tube1fcs$Status=="Normal")),]
  # }
  # 
  # 
  # paths <- apply(tube1fcs_1, 1, function(r){
  #   if(is.na(r['Path3'])){
  #     path <- paste0(r['Path1'],"/",r['Path2'],"/",r['FCS.files'])
  #   }else{
  #     path <- paste0(r['Path1'],"/",r['Path2'],"/",r['Path3'],"/",r['FCS.files'])
  #   } 
  #   return(path)
  # })
  # 
  # f_list <- preprocessFCS_files(paths,scatter_chans)
  # 
  # # lv2 <- data.frame(id=tube1fcs_1[intersect(names(f_list),rownames(tube1fcs_1)),"PatientID"],
  # #                   st=tube1fcs_1[intersect(names(f_list),rownames(tube1fcs_1)),"AML_Status"])
  # lv2 <- tube1fcs_1[intersect(names(f_list),rownames(tube1fcs_1)),]
  # 
  # fs.raw <- flowSet(f_list)
  # fs <- fs.raw
  # scat.chans <- c(grep(colnames(fs),pattern = "FSC*"), grep(colnames(fs),pattern = "SSC*"))
  # names(scat.chans) <- colnames(fs)[scat.chans]
  # time.loc <- which(tolower(colnames(fs)) == "time"); names(time.loc) <- NULL
  # 
  # fs <- fsApply(fs,function(f) removeMargins(f,chans = scat.chans, verbose = T))
  # 
  # fs <- fsApply(fs, function(f) removeMargins(f,chans = scat.chans,debris = T, neg = T, verbose = T))
  # 
  # fs <- fsApply(fs, function(f) compensate(f,f@description$SPILL))
  # 
  # gf <- globalFrame(fs)
  # lgl <- gf$lgl
  # 
  # fs <- fsApply(fs, function(f) transform(f, lgl))
  # 
  # 
  # fs.preflowCut <- fs
  # 
  # fcut <- fsApply(fs, function(f){
  #   channels.to.clean <- which(complete.cases(f@parameters@data$desc) == TRUE)
  #   time.loc <- unique(c(grep("time", tolower(colnames(f)))))
  #   
  #   f.Clean <- flowCut(f, Channels = channels.to.clean,
  #                      Plot = 'None', PrintToConsole = FALSE)
  #   
  #   return(list(ind = f.Clean$ind, passed.flowCut = f.Clean$data['Has the file passed',], f=f.Clean$frame))
  # })
  # 
  # 
  # fs <- flowSet(lapply(fcut, function(f){return(f$f)}))
  # 
  # 
  # 
  # markers_list <- as.character(fs[[1]]@parameters@data[1:ncol(fs[[1]]), 2])
  # markers_list <- Find.markers(fs[[1]],na.omit(markers_list))
  # channel_names <-  as.character(fs[[1]]@parameters@data[1:ncol(fs[[1]]), 1])
  # 
  # 
  # 
  # #===============================gating live===================================================
  # live.flowD <- llply(as(fs, Class="list"), function(f){
  #   
  #   peaks <- getPeaks(f,channel =  c("FSC-A"))
  #   if(length(peaks$Peaks) == 4){
  #     peaks$Peaks <- peaks$Peaks[-which(peaks$Peaks < 50000)]
  #   }
  # 
  #   
  #   if(length(peaks$Peaks)>=3){
  #     allcuts <- deGate(f, channel = c("FSC-A"),all.cuts = T) #find gates at head
  #     gt <- min(allcuts)
  #     if(gt > 100000){
  #       gt <- deGate(f, channel = c("FSC-A"), use.upper = T, upper=F)
  #     }
  #     # if(gt > peaks$Peaks[2]){
  #     #   if(50000 > peaks$Peaks[1]){
  #     #     gt <- 50000
  #     #   }else{
  #     #     gt <- deGate(f, channel = c("FSC-A"), use.upper = T, upper=F, tinypeak.removal = 0.1) 
  #     #   }
  #     # }
  #     
  #   }else if(length(peaks$Peaks)==2){
  #     
  #     allcuts <- deGate(f, channel = c("FSC-A"),all.cuts = T,use.percentile = FALSE,
  #                       percentile = 0.99)
  #     gt <- min(allcuts)
  #     if(length(allcuts) == 2){
  #       if(allcuts[1] < 50000 && allcuts[2] < 50200){
  #         gt <- allcuts[2] 
  #       }
  #     }
  #     
  #     gt_perc <- deGate(f, channel = c("FSC-A"),all.cuts = T,use.percentile = TRUE,
  #                       percentile = 0.99) 
  #     if(gt == gt_perc ){
  #       gt <- deGate(f, channel = c("FSC-A"), use.upper = T, upper=F)
  #     }
  #     if(gt > 100000){
  #       gt <- deGate(f, channel = c("FSC-A"), use.upper = T, upper=F) #find gates at head
  #     }
  #     
  #    
  #     # if(peaks$Peaks[2] < 150000){
  #     #   allcuts <- deGate(f, channel = c("FSC-A"),all.cuts = T)
  #     #   gt <- min(allcuts)
  #     # }else{
  #     #   gt <- deGate(f, channel = c("FSC-A"), use.upper = T, upper=F) #find gates at head
  #     # }
  #   
  #   }else{
  #     gt <- deGate(f, channel = c("FSC-A"), use.upper = T, upper=F) #find gates at head
  #     
  #   }
  #   
  #   # if(gt < 28000){
  #   #   gt <- deGate(f, channel = c("FSC-A"))
  #   #   if(gt > 100000){
  #   #     gt <- 50000
  #   #   }
  #   # }
  #   fD <- flowDensity(f, channels = c("FSC-A","SSC-A"), position=c(T,NA), gates = c(gt,NA))
  #   return(fD)
  # }, .parallel=T)
  # 
  # par(mfrow=c(2,3))
  # a <- lapply(1:length(live.flowD), function(x){
  #   plotDens(fs[[x]],channels=c("FSC-A","SSC-A"))
  #   lines(live.flowD[[x]]@filter)
  # })
  # 
  # live.fs <- flowSet(lapply(live.flowD, function(fD){getflowFrame(fD)}))
  # 
  # 
  # #============================gating singlets=================================================
  # 
  # if(length(grep('FSC-H',colnames(live.fs))) == 0){
  #   print("no FSC-H channels, skipping singlets gating")  
  #   singlets.fs <- live.fs
  # }else{
  #   singlets.flowD <- llply(as(live.fs, Class = "list"),function(f){
  #     if(length(grep('FSC-H', colnames(f))) == 0){
  #       print("no FSC-H channels, skipping singlets gating")  
  #       return(f)
  #     }else{
  #       flowD <- GateSinglets(f,scat.chans,channels.ind)
  #       if((length(flowD@index)/nrow(f)) < 0.7){
  #         theta <- pi/5
  #         chans <- c(grep("FSC-A", colnames(f)), grep("FSC-H", colnames(f)))
  #         rot <- rotate.data(f, chans = chans, theta = theta)$data
  #         gate <- deGate(rot, channel = chans[2], use.upper = T, upper = F, alpha = .01)
  #         flowD <- flowDensity(rot, channels = chans, gates = c(NA, gate), position = c(NA,T) )
  #         flowD@flow.frame <- rotate.data(flowD@flow.frame, chans = chans, theta = -theta)$data
  #         flowD@filter <- rotate.data(flowD@filter, chans = chans, theta = -theta)$data
  #       }
  #       return(flowD)
  #     }
  #   }, .parallel=T)
  #   
  #   singlets.fs <- flowSet(lapply(singlets.flowD, function(fD){getflowFrame(fD)}))
  #   
  #  
  # }
  
  #========================== gating CD33+CD117+B mast cells===============================================
  # cd33poscd117posB.flowD <- fsApply(singlets.fs,function(f){
  #   # gt.cd117 <- mean(f@exprs[,markers_list["CD117"]]) + 5*sd(f@exprs[,markers_list["CD117"]])
  #   all.cuts <- deGate(f, c(markers_list["CD117"]), upper = T, all.cuts = T, tinypeak.removal = 0.0001)
  #   
  #   denf <- density(f@exprs[,markers_list["CD117"]])
  #   x.max <- denf$x[which.max(denf$y)]
  #   gt.cd117 <- min(intersect(denf$x[denf$y < 0.01],
  #                             denf$x[denf$x > x.max]))
  #   gt.cd117 <- all.cuts[min(which(all.cuts > gt.cd117))]
  #   flowD <- flowDensity(f, channels = c(markers_list["CD117"],markers_list["CD33"]), 
  #                        position = c(T,F),gates=c(gt.cd117,3.55))
  #   return(flowD)
  # })
  # 
  # par(mfrow=c(2,3))
  # a <- lapply(1:length(cd33poscd117posB.flowD), function(x){
  #   plotDens(singlets.fs[[x]],channels=c(markers_list["CD117"],markers_list["CD33"]))
  #   lines(cd33poscd117posB.flowD[[x]]@filter)
  # })
  # 
  # cd33poscd117posB.fs <- flowSet(lapply(cd33poscd117posB.flowD, function(fD){getflowFrame(fD)}))
  
 
 
  
  load(file = paste0("~/results/GatingSet_tube",tubes,"_",days[d], ".RData")) 
  fs <- singlets.fs
  markers_list <- as.character(fs[[1]]@parameters@data[1:ncol(fs[[1]]), 2])
  markers_list <- Find.markers(fs[[1]],na.omit(markers_list))
  channel_names <-  as.character(fs[[1]]@parameters@data[1:ncol(fs[[1]]), 1])
  scat.chans <- c(grep(colnames(fs),pattern = "FSC*"), grep(colnames(fs),pattern = "SSC*"))
  names(scat.chans) <- colnames(fs)[scat.chans]
  time.loc <- which(tolower(colnames(fs)) == "time"); names(time.loc) <- NULL
  #===========================gate HLA-DR mast cells======================================
  par(mfrow=c(2,3))
  hladrmastcell.flowD <- llply(as(singlets.fs, Class = "list"), function(f){
    chans <- c(markers_list[grep(paste("HLA-DR|HLADR", collapse = "|"),names(markers_list))],
               markers_list[grep(paste("CD117|CD117BC",collapse = "|"),names(markers_list))])
    gate1 <- deGate(f,chans[1],upper= TRUE, use.upper = T, tinypeak.removal = 0.1,alpha = 0.01)
    
    # if(gate1 > 2.5+0.5){
    #   if(center=="karolinska"){
    #     gate1 <- 3.2
    #   }else{
    #     gate1 <- 2.5
    #   }
    # }
  
    temp <- flowDensity(f,chans, position = c(F,NA), gates= c(gate1,NA))

    rot <- rotate.data(getflowFrame(temp),chans,theta=0.18)$data
    # rot <- rotate.data(f,chans,theta=0.18)$data
    
    gate <- deGate(rot,chans[2],upper= TRUE, use.upper = T, tinypeak.removal = 0.01,alpha = 0.01)
    gate2 <- deGate(rot,chans[2],upper= TRUE, use.upper = T, tinypeak.removal = 0.01,alpha = 0.9)
    
    gate <- c(gate,gate2)
    
  
    gate <- gate[which.min(abs(gate-3.3))] 
      if(gate < 3.1){
        gate <- 3.2
      }

      if(abs(gate-3.3) > 0.5){
        gate <- 3.3
      }
    
    if(f@description$`EXPERIMENT NAME` == "a 280503-6071_5570_17_Dag29_eInd1"){
      gate = deGate(rot,chans[2],upper= TRUE, use.upper = T, tinypeak.removal = 0.01,alpha = 0.05)
      gate <- gate
      }
    
    flowD <- flowDensity(rot,chans,position = c(NA,T), gates= c(NA,gate))
    flowD@flow.frame <- rotate.data(flowD@flow.frame,chans,theta=-0.18)$data
    
    # gate3 <-  deGate(flowD@flow.frame,chans[1],upper= TRUE, use.upper = T, tinypeak.removal = 0.1,alpha = 0.9)
    # 
    # if(gate3 > 2.5+0.5){
    #   if(center=="karolinska"){
    #     gate3 <- 3.2
    #   }else{
    #     gate3 <- 2.5
    #   }
    # }
    # 
    # temp <- flowDensity(flowD,chans,position = c(F,NA), gates=c(gate3,NA))
    # 
    # flowD@filter <- temp@filter
    # 
    flowD@filter <- rotate.data(flowD@filter,chans, theta=-0.18)$data
 
    if(flowD@cell.count <=1){
      return(flowD)
    }else{
      if(any(is.infinite(flowD@filter))|any(is.nan(flowD@filter))){
        return(flowD)
      }
      
    }
    
    
    #check if the filter has missed any cells right outside its borders
    test.filter <- flowD@filter
    test.flowD <- flowDensity(f,chans,position = c(F,T),gates=c(max(test.filter[,1])+0.2, min(test.filter[,2])))
    
    if(test.flowD@cell.count > flowD@cell.count * 2 &&
       test.flowD@cell.count < flowD@cell.count *25
       ){#if with 0.2 in x-expansion,the filter captures more than 
                                                    #twice the cell populations, then it can be inferred there
                                                    #are cells outside the filter that was missed previously
      
      test.temp <- flowDensity(f,chans,position = c(NA,T),gates=c(NA,min(test.filter[,2])))
      test.gate <- deGate(test.temp,chans[1],upper=T)
      if(test.gate > 4){test.gate <- 2.5}
        
      flowD <- flowDensity(test.temp,chans,position=c(F,NA),gates=c(test.gate,NA))
    }

    
    
    
    return(flowD)
    
  }, .parallel=T)
  
  hladrmastcell.fs <- flowSet(lapply(hladrmastcell.flowD, function(fD){getflowFrame(fD)}))
  chans <- c(markers_list[grep(paste("HLA-DR|HLADR", collapse = "|"),names(markers_list))],
             markers_list[grep(paste("CD117|CD117BC",collapse = "|"),names(markers_list))])
  par(mfrow=c(2,3))
  a <- lapply(1:length(hladrmastcell.fs), function(x){
    plotDens(singlets.fs[[x]],channels=chans)
    lines(hladrmastcell.flowD[[x]]@filter)
  })
  
  all.gthres[[1]] <- llply(as(hladrmastcell.flowD,Class="list"), function(fD){
    chans <- c(markers_list[grep(paste("HLA-DR|HLADR", collapse = "|"),names(markers_list))],
               markers_list[grep(paste("CD117|CD117BC",collapse = "|"),names(markers_list))])
    
    ft <- fD@filter
    colnames(ft) <- c(channel_names[chans[1]], channel_names[chans[2]])
  
    return(ft)
  })
  
  # names(all.gthres[[1]]) <- rep(paste0("HLADR.mast.filter"),length(hladrmastcell.flowD))
 
  # llply(as(singlets.fs, Class = "list"), function(f){
  #   chans <- c(markers_list["CD13"],markers_list["CD33"])
  #   flowPeaks.Res <-  flowPeaks::flowPeaks(f@exprs[,chans],tol=0.1,h0=1,h=1.5)
  #   plot(flowPeaks.Res)
  #   })
  # 
   
  #===============================gating CD13+B CD33- =========================================
  cd13posbcd33neg.flowD <- llply(as(singlets.fs, Class = "list"), function(f){
    
    chans <- c(markers_list["CD13"],markers_list["CD33"])
    f.org <- f
    temp <- f
    temp <- flowDensity(temp,chans,position = c(NA,T), gate=c(NA,0))
    f <-temp@flow.frame
    
      rot <- rotate.data(f,chans, theta = 1.2)$data
    
    
    gate <- deGate(rot,chans[2],upper = F, tinypeak.removal = 0.0001)
    # gate <- deGate(rot,chans[2],upper = F, tinypeak.removal = 0.001,use.upper = T)
    
    if(abs(2.9-abs(gate)) > 0.5){
      if(identifier(f)=="10"){
        gate <- -3.3
      }else{
        gate <- -2.9
      }
     
    }
    
    temp <- flowDensity(rot,chans,position = c(NA,F),gates=c(NA,gate))

    if(temp@proportion > 10){
      gate1 <- deGate(temp,chans[2],use.upper = F,tinypeak.removal = 0.0001)
      temp1 <- flowDensity(temp,chans,position = c(NA,F),gates=c(NA,gate1))
      temp <- temp1
    }
  
      temp@filter <- rotate.data(temp@filter, chans,theta = -1.2)$data
      # temp@gates <- rotate.data(temp@gates,chans,theta=-1.2)$data
      temp@flow.frame <- rotate.data(getflowFrame(temp), chans,theta = -1.2)$data
      
    temp2 <- flowDensity(temp,chans,position = c(NA,F), gate=c(NA,2.5))
    
    # plotDens(f,chans)
    # lines(temp@filter)
    
    
    # gate1 <- deGate(f,chans[1], tinypeak.removal = 0.01, upper = T,use.upper = T)
    return(temp2)
  }, .parallel=T)
  
  cd13posbcd33neg.fs <- flowSet(lapply(cd13posbcd33neg.flowD, function(fD){getflowFrame(fD)}))
  # cd13posbcd33neg.flowD1 <- fsApply(cd13posbcd33neg.fs,function(f){
  #   
  # })
  par(mfrow=c(2,3))
  chans <- c(markers_list["CD13"],markers_list["CD33"])
  a <- lapply(1:length(cd13posbcd33neg.fs), function(x){
    plotDens(singlets.fs[[x]],channels=chans)
    lines(cd13posbcd33neg.flowD[[x]]@filter)
  })
  
  all.gthres[[2]] <- llply(as(cd13posbcd33neg.flowD,Class="list"), function(fD){
   
    ft <- fD@filter
    colnames(ft) <- c(channel_names[chans[1]], channel_names[chans[2]])
    
    return(ft)
  })
  
  # names(all.gthres[[2]]) <- rep(paste0("CD33CD13.filter"),length(hladrmastcell.flowD))
  
  #=============================gate CD45-/CD45+D cells=======================================
  
  cd45negcd45posD.flowD <- llply(as(cd13posbcd33neg.fs, Class = "list"), function(f){
    chans <- c(markers_list[grep("CD45",names(markers_list))], scat.chans["SSC-A"])
    gate <- deGate(f,chans[1],upper = T,all.cuts = T,tinypeak.removal = 0.01)
    if(length(gate) > 1){#if there are several gates, pick the one that is closest to 3.5
      dens <- density(f@exprs[,chans[1]])
      xmax <- dens$x[which.max(dens$y)]
      # gate <- min(gate[which(gate>xmax)])
      gate <- gate[which.min(abs(gate-3.5))]
      if(gate < xmax){#if all gates are smaller than highest density point, then neeed to find gate at the tail
        gate <- deGate(f,chans[1], upper=T, use.upper = T)
      }
    }else{
      if(abs(3.5-gate) < 0.5){
        gate <- gate
      }else{
        dens <- density(f@exprs[,chans[1]])
        xmax <- dens$x[which.max(dens$y)]
        if(gate < xmax){
          gate <- deGate(f,chans[1],upper = T,use.upper=T)
          # if(gate==-Inf){
          #   gate <- 1000000
          # }
        }
      }
      
      
    }
    
    
    flowD <- flowDensity(f,chans,position = c(F,NA),gates=c(gate,NA))

    return(flowD)
  }, .parallel=T)
  par(mfrow=c(2,3))
  chans <- c(markers_list[grep("CD45",names(markers_list))], scat.chans["SSC-A"])
  a <- lapply(1:length(cd45negcd45posD.flowD), function(x){
    plotDens(cd13posbcd33neg.fs[[x]],channels=chans,pch=20)
    abline(v=cd45negcd45posD.flowD[[x]]@gates)
  })
  
  all.gthres[[3]] <- llply(as(cd45negcd45posD.flowD,Class="list"), function(fD){
    
    ft <- fD@gates[1]
    # colnames(ft) <- c(channel_names[chans[1]], channel_names[chans[2]])
    
    return(ft)
  })
  
  
  
  cd45negcd45posD.fs <- flowSet(lapply(cd45negcd45posD.flowD, function(fD){getflowFrame(fD)}))
  # cd45negcd45posD.gates <- unlist(lapply(cd45negcd45posD.flowD,function(fD){fD@gates[1]}))
  
  #=========================gate CD34-mesenchymal stromacells/fibroblasts===================
  par(mfrow=c(2,3))
  cd34negstroma.flowD <- llply(as(cd45negcd45posD.fs, Class = "list"), function(f){
    chans <- c(markers_list[grep(paste("CD117|CD117BC|CD117",collapse = "|"),names(markers_list))],
               markers_list[grep("CD34",names(markers_list))])
    
    # if(nrow(f@exprs)<=3){
    #   gate <- 100000
    # }else{
      gate <- deGate(f,chans[2],upper = T,use.upper = T, tinypeak.removal = 0.01,alpha = 0.9)
      if(gate > 4.2){
        gate <-3
      }
    # }
   
    if(gate==-Inf){
      gate <- 100000
    }
   
    # plot(density(na.omit(f@exprs[,chans[2]])))
    # plot(density(na.omit(f@exprs[,chans[1]])))
    gate1 <- deGate(f,chans[1],upper = T,use.upper = T, tinypeak.removal = 0.01,alpha = 0.9)
    if(gate1==-Inf){
      gate1 <- 100000
    }
    flowD <- flowDensity(f,chans,position = c(F,F),gates=c(gate1,gate))

    return(flowD)
  }, .parallel=T)
  
  
  cd34negstroma.fs <- flowSet(lapply(cd34negstroma.flowD, function(fD){getflowFrame(fD)}))
  
  par(mfrow=c(2,3))
  chans <-  chans <- c(markers_list[grep(paste("CD117|CD117BC|CD117",collapse = "|"),names(markers_list))],
                       markers_list[grep("CD34",names(markers_list))])
  a <- lapply(1:length(cd34negstroma.flowD), function(x){
    plotDens(cd45negcd45posD.fs[[x]],channels=chans,pch=20)
    # abline(h=cd34negstroma.flowD[[x]]@gates[2])
    # abline(v=cd34negstroma.flowD[[x]]@gates[1])
    lines(cd34negstroma.flowD[[x]]@filter)
  })
  
  all.gthres[[4]] <- llply(as(cd34negstroma.flowD,Class="list"), function(fD){
    
    ft <- fD@filter
    colnames(ft) <- c(channel_names[chans[1]], channel_names[chans[2]])
    
    return(ft)
  })
  
  # all.gthres[[5]] <- llply(as(cd34negstroma.flowD,Class="list"), function(fD){
  #   
  #   ft <- fD@filter
  #   colnames(ft) <- c(channel_names[chans[1]], channel_names[chans[2]])
  #   
  #   return(ft)
  # })
  
  
  #============================gating CD33+HLA-DR+==============================================
  # if(days[x]=="Day22"){
  #   cd33poshladrpos.flowD<- llply(as(singlets.fs, Class="list"),function(f){
  #     chans <- c(markers_list[grep(paste("HLA-DR|HLADR", collapse = "|"),names(markers_list))],markers_list["CD33"])
  #     set.seed(1234)
  #     flowPeaks.Res <- flowPeaks::flowPeaks(f@exprs[,chans],tol=0.1,h0=1,h=1.5)
  #     # plot(flowPeaks.Res)
  #     pop.cid <- flowPeaks.Res$peaks$cid[intersect(which(flowPeaks.Res$peaks$mu[,2] > 1.9),
  #                                                  which(flowPeaks.Res$peaks$mu[,1] > 2.2))]
  #     pop.f <- f
  #     pop.f@exprs <- f@exprs[which(flowPeaks.Res$peaks.cluster %in% pop.cid), ]
  #     flowD <- flowDensity(pop.f,chans,position = c(T,T),gate=c(0,0))
  #     if( min(flowD@filter[,2]) < 1){
  #       flowD <- flowDensity(pop.f,chans,position = c(NA,T),gate=c(NA,2))
  #     }
  #     
  #     return(flowD)
  #   })
  # }else{
    
    cd33poshladrpos.flowD<- llply(as(singlets.fs, Class="list"),function(f){
      
      chans <- c(markers_list[grep(paste("HLA-DR|HLADR", collapse = "|"),names(markers_list))],markers_list["CD33"])
      #check for marginal events on CD33
      # remove cells with CD33 > 4
      if(length(which(f@exprs[,chans[2]] > 4)) > 0){
        f@exprs <- f@exprs[-which(f@exprs[,chans[2]] > 4),]
      }
      
      
      gate <- deGate(f,chans[2],tinypeak.removal = 0.001,alpha=0.01)
      if(gate < 1.3){
        gate <- deGate(f,chans[2],tinypeak.removal = 0.001,alpha=0.1)
        if(gate < 1.2){
          gate <- 2.3
        }
      }
      if(gate > 3.6){
        gate <- deGate(f,chans[2],tinypeak.removal = 0.01,alpha=0.01)
      }
      
      temp <- flowDensity(f,chans,position = c(NA,T),gates=c(NA,gate))
      
      
      if(temp@proportion < 5){#if unable to find good gate on x direction, switch to another
        gate <- deGate(f,chans[1],tinypeak.removal = 0.001,alpha=0.01)
        if(gate > 3.4){
          gate <- deGate(f,chans[1],tinypeak.removal = 0.01,alpha=0.01)
        } 
        if(gate < 1.4){
          gate <- 2.3
        }
        temp <- flowDensity(f,chans,position = c(T,NA),gates=c(gate,NA))
        gate1 <- deGate(temp,chans[2],tinypeak.removal = 0.01,alpha=0.1)
        if(gate1 < 1.4){
          gate1 <- 2.3
        }
        if(gate1 > 3.5){
          gate1 <- 2.3
        }
        temp1 <- flowDensity(temp,chans,position = c(NA,T),gates=c(NA,gate1))
        gate2 <- deGate(temp1,chans[2],tinypeak.removal = 0.01,use.percentile = TRUE,percentile = 0.99)
        
        temp2<- flowDensity(temp1,chans,position=c(NA,F),gates=c(NA,gate2))
        if(temp2@proportion==0|temp2@proportion==100){#if there is no cells being capured, create a standard filter
          gate3 <- deGate(f,chans[2],tinypeak.removal = 0.01,use.percentile = TRUE,percentile = 0.9999)
          tempx <- flowDensity(f,chans,position = c(NA,F),gates=c(NA,gate3))
          temp2 <- flowDensity(tempx, chans, position = c(T,T),gates=c(2,1.8))
        }
        return(temp2)
      }else{
        if(days[d]=="Day22"){
          gate1 <- deGate(temp,chans[1],tinypeak.removal = 0.001, alpha=0.9, upper = F, use.upper = T)
          if(gate1 < 1.5){
            gate1 <- 2
          }
        }else{
          gate1 <- deGate(temp,chans[1],tinypeak.removal = 0.001,alpha=0.01,all.cuts = T)
          
          if(length(gate1) > 1){
            denx <- density(getflowFrame(temp)@exprs[,chans[1]])
            x.max <- denx$x[which.max(denx$y)]
            gate1 <- min(gate1[which(gate1>x.max)])
          }
        }
       
        gate1_1 <- deGate(f,chans[1],upper = F,tinypeak.removal = 0.01)
        
        gates <- c(gate1,gate1_1)
        
        gate1 <- gates[which.min(abs(gates-2))]
        
        if(abs(gate1-2) > 1.2){
          gate1 <- 2
        }
        if(days[d]=="Day22"){
          temp3 <- flowDensity(temp, chans, position = c(T, NA), gates = c(gate1,NA), ellip.gate = T)
          
        }else{
          gate2 <- deGate(temp,chans[2],tinypeak.removal = 0.001,alpha=0.01,
                          use.percentile = T, percentile = 0.999)
          temp1 <- flowDensity(getflowFrame(temp),
                               channels =chans,
                               position = c(T,F),
                               gates=c(gate1,gate2))
          gate3<-deGate(temp1,chans[2], upper=F,use.upper = T,tinypeak.removal = 0.8)
          temp3 <- flowDensity(getflowFrame(temp1),
                               channels =chans,
                               position = c(NA,T),
                               gates=c(NA,gate3))
        }
       
    
        
        return(temp3)
      }
    }, .parallel=T)
  # }

  
  
  par(mfrow=c(2,3))
  a <- lapply(1:length(cd33poshladrpos.flowD), function(x){
    plotDens(singlets.fs[[x]],channels=c(markers_list[grep(paste("HLA-DR|HLADR", collapse = "|"),names(markers_list))],markers_list["CD33"]))
    lines(cd33poshladrpos.flowD[[x]]@filter)
  })
  
  chans <- c(markers_list[grep(paste("HLA-DR|HLADR", collapse = "|"),names(markers_list))],markers_list["CD33"])
  
  cd33poshladrpos.fs <- flowSet(lapply(cd33poshladrpos.flowD, function(fD){getflowFrame(fD)}))
  all.gthres[[5]] <- llply(as(cd33poshladrpos.flowD,Class="list"), function(fD){
    
    ft <- fD@filter
    colnames(ft) <- c(channel_names[chans[1]], channel_names[chans[2]])
    
    return(ft)
  })
  
  #=============================gate CD11b+============================================
  
  cd11bpos.flowD <- llply(as(cd33poshladrpos.fs, Class = "list"),function(f){
    chans <- c(markers_list["CD56"],
               markers_list[grep(paste("CD11b|CD11B",collapse = "|"),names(markers_list))])
    gate <- deGate(f,chans[2],tinypeak.removal = 0.1,upper = F, alpha=0.9)
    if(days[d]!="Day22"){
      if(gate > 2.75){
        gate <- 2.25
      }
    }else{
      if(gate > 3){
        gate <- 2.5
      }
      if(gate < 0.5){
        gate <- mean(f@exprs[,chans[2]])
      }
    }
    # if(gate < 0.5){
    #   peak.cd56 <- getPeaks(f,chans[2])
    #   sd.cd56 <- sd(f@exprs[,chans[2]])
    #   gate <- peak.cd56$Peaks + sd.cd56
    # }
    flowD <- flowDensity(f,chans,position=c(NA,T),gates= c(NA,gate))
    
  
    return(flowD)
  }, .parallel=T)
  
  par(mfrow=c(2,3))
  a <- lapply(1:length(cd33poshladrpos.fs), function(x){
    plotDens(cd33poshladrpos.fs[[x]],channels=c(markers_list["CD56"],
                                                markers_list[grep(paste("CD11b|CD11B",collapse = "|"),names(markers_list))]))
    lines(cd11bpos.flowD[[x]]@filter)
  })
  chans <- c(markers_list["CD56"],
             markers_list[grep(paste("CD11b|CD11B",collapse = "|"),names(markers_list))])
  cd11bpos.fs <- flowSet(lapply(cd11bpos.flowD, function(fD){getflowFrame(fD)}))
  
  all.gthres[[6]] <- llply(as(cd11bpos.flowD,Class="list"), function(fD){
    
    ft <- fD@gates[2]
    # colnames(ft) <- c(channel_names[chans[1]], channel_names[chans[2]])
    
    return(ft)
  })
  
  
  #=========================== gate SSC+=============================================

  sscpos.flowD <- llply(as(cd11bpos.fs, Class = "list"),function(f){
    chans <- c(markers_list[grep("CD45",names(markers_list))],scat.chans["SSC-A"])
    gate <- deGate(f,chans[2],tinypeak.removal = 0.01, alpha=0.1,upper = TRUE)
  
    temp <- flowDensity(f,chans,position = c(NA,F),gates=c(NA,gate))
    
    
    
    gate1 <- deGate(temp,chans[2],tinypeak.removal = 0.01,alpha=0.1,upper = TRUE, twin.factor = 0.7)

    if(gate1 > 150000){
      gate1 <- 100000
      # if(center == "karolinska"){
      #   gate1 <- 150000
      # }
    }
    if(gate1 <30000){
      gate1 <- gate
      if(gate1 < 30000){
        gate1 <-deGate(f,chans[2], alpha=0.1,upper = TRUE) 
      }
    }
    if(gate1 == -Inf){
      gate1 <- 50000
    }
    
    flowD <- flowDensity(f,chans,position=c(NA,T),gates=c(NA,gate1))
    notsscpos <- getflowFrame(notSubFrame(f,chans,position = c(NA,T),gates = c(NA,gate1)))
    
    return(list(sscpos.flowD=flowD, notsscpos=notsscpos))
  }, .parallel=T)
  
  par(mfrow=c(2,3))
  a <- lapply(1:length(cd11bpos.flowD), function(x){
    chans <- c(markers_list[grep("CD45",names(markers_list))],scat.chans["SSC-A"])
    plotDens(cd11bpos.fs[[x]],chans)
    lines(sscpos.flowD[[x]]$sscpos.flowD@filter)
  })
  
  
  
  sscpos.fs <- flowSet(lapply(1:length(sscpos.flowD), function(idx){
    fD<- sscpos.flowD[[idx]]
    return(fD$sscpos.flowD@flow.frame)}))
  
  notsscpos.fs <- flowSet(lapply(1:length(sscpos.flowD), function(idx){
    fD<- sscpos.flowD[[idx]]
    return(fD$notsscpos)}))
  chans <- c(markers_list[grep("CD45",names(markers_list))],scat.chans["SSC-A"])
  all.gthres[[7]] <- llply(as(sscpos.flowD,Class="list"), function(fD){
    
    ft <- fD$sscpos.flowD@gates[2]
    # colnames(ft) <- c(channel_names[chans[1]], channel_names[chans[2]])
    
    return(ft)
  })
  #==================================gate CD34-CD117-=======================================
  
  cd34negcd117neg.flowD <- llply(as(sscpos.fs, Class = "list"),function(f){
    chans <- c(markers_list[grep(paste("CD117|CD117BC",collapse = "|"),names(markers_list))],
               markers_list[grep("CD34",names(markers_list))])
    gate <- deGate(f,chans[2],twin.factor = 0.95,upper=T, alpha = 0.9)
    if(gate < 1.6){
      gate <- deGate(f,chans[2],upper=T,use.upper = T)
      
    }
    if(gate==-Inf){
      gate <- 3
    }
    gate1 <- deGate(f,chans[1],use.upper = TRUE,upper = TRUE)
    if(gate1 > 3){
      gate1 <- deGate(f,chans[1])
      if(gate1 > 2.5){
        gate1 <- 2.5
      }
      if(gate1 <1.5){
        gate1 <- 2.5
      }
    }
    # if(gate1 > 2.5){
    #   gate1 <- 2.5
    # }
    if(gate1==-Inf){
      gate1 <- 2.5
    }
    flowD <- flowDensity(f,chans,position = c(F,F),gates=c(gate1,gate))
    
   
    return(flowD)
  }, .parallel=T)
  
  par(mfrow=c(2,3))
  a <- lapply(1:length(sscpos.flowD), function(x){
    plotDens(sscpos.fs[[x]],channels=c(markers_list[grep(paste("CD117|CD117BC",collapse = "|"),names(markers_list))],
                                       markers_list[grep("CD34",names(markers_list))]))
    lines(cd34negcd117neg.flowD[[x]]@filter)
  })
  
  
  cd34negcd117neg.fs <- flowSet(lapply(cd34negcd117neg.flowD,function(fD){getflowFrame(fD)}))
  chans <- c(markers_list[grep(paste("CD117|CD117BC",collapse = "|"),names(markers_list))],
             markers_list[grep("CD34",names(markers_list))])
  
  all.gthres[[8]] <- llply(as( cd34negcd117neg.flowD,Class="list"), function(fD){

    ft <- fD@gates[1]
    # colnames(ft) <- c(channel_names[chans[1]], channel_names[chans[2]])
    
    return(ft)
  })
  all.gthres[[9]] <- llply(as( cd34negcd117neg.flowD,Class="list"), function(fD){
    
    ft <- fD@gates[2]
    # colnames(ft) <- c(channel_names[chans[1]], channel_names[chans[2]])
    
    return(ft)
  })
  
  #====================================gate singlet monocyte derived cells high ssc================
  if(length(grep("FSC-H",colnames(singlets.fs))) == 0){
    
    print("no FSC-H channeL, skipping gating singlet monocyte derived cells high ssc..")
    singmonohighssc.fs <- sscpos.fs
    
  }else{
  
    singmonohighssc.flowD <- llply(as(sscpos.fs, Class = "list"),function(f){
      chans <- c(scat.chans["FSC-A"],scat.chans["FSC-H"])
      if(nrow(na.omit(f@exprs)) <3 ){#creating a dummy gate if no cells from parent population
        
        res <- flowDensity(f,chans,position = c(T,T),gates=c(-10,-10))
        res@flow.frame <- f
      }else{
        res <- GateSinglets(f,scat.chans,channels.ind)
        plotDens(f,chans,pch=20)
        lines(res@filter)
      }
      
      
      # flowSOM.res <- ReadInput(f)
      # flowSOM.res <- BuildSOM(flowSOM.res,7:14)
      # flowSOM.res <- BuildMST(flowSOM.res)
      # 
      # labels_pre <- flowSOM.res$map$mapping[, 1]
      # 
      # k <- 3
      # 
      # seed <- 1234
      # out <- FlowSOM::metaClustering_consensus(flowSOM.res$map$codes, k = k, seed = seed)
      # 
      # labels <- out[labels_pre]
      
     
      return(res)
      
    }, .parallel=T)
    
    singmonohighssc.fs <- flowSet(lapply(singmonohighssc.flowD,function(fD){getflowFrame(fD)}))
    par(mfrow=c(2,3))
    a <- lapply(1:length(sscpos.flowD), function(x){
      plotDens(sscpos.fs[[x]],channels=c(scat.chans["FSC-A"],scat.chans["FSC-H"]))
      lines(singmonohighssc.flowD[[x]]@filter)
    })
    
  }
  chans <- c(scat.chans["FSC-A"],scat.chans["FSC-H"])
  all.gthres[[10]] <- llply(as(singmonohighssc.flowD,Class="list"), function(fD){
    
    ft <- fD@filter
    colnames(ft) <- c(channel_names[chans[1]], channel_names[chans[2]])
    
    return(ft)
  })
  
  Gates.list <- NULL

    # for(l in 1:length(all.gthres)){
    #   Gates.list[[l]] <-all.gthres[[l]][[k]]
    # }
    for(k in 1:length(singlets.fs)){
      Gates.list[[k]] <- lapply(all.gthres, function(gl){gl[[k]]})
    }
    
    for(l in 1:length(singlets.fs)){
      names(Gates.list[[l]]) <- c("HLADR.mastcells.filter","CD13BCD33.filter","CD45.gate",
                                  "CD34stroma.filter","CD33HLADR.filter",
                                  "CD11B.gate","SSC.gate",
                                  "CD117.gate","CD34.gate",
                                  "monocyteshighssc.filter")
    }
 
  
  
  save(singlets.fs,Gates.list,lv2,file=paste0("~/results/gatethreslist_",days[d],".RData"))
  
  ##------------------------------------------------------------------
  ##-----------     PLOTTING                   -----------------------
  ##------------------------------------------------------------------
  if(!skip_plot){
    
  
  dir.create("~/results/Plots/AutomatedGating/", recursive = T)
  l_ply(1:length(fs.preflowCut), function(idx){
    f <- fs[[idx]]
    fileName <- gsub('.fcs',".png",
      unlist(strsplit(f@description$FILENAME,"/"))[length(unlist(strsplit(f@description$FILENAME,"/")))]
    )
    
    png(file=paste0("~/results/Plots/AutomatedGating/",fileName), width = 2000*4/4, height = 2000*4/4)
    par(mfrow = c(3,4), mar = (c(5, 5, 4, 2) + 0.1))
    
    
    plotDens(fs[[idx]],channels=c("FSC-A","SSC-A"), main=paste0("CLOG-: "),
             cex.lab = 2, cex.axis = 2, cex.main = 2)
    lines(live.flowD[[idx]]@filter,type="l",lwd=2)
    
    if(length(grep("FSC-H",colnames(singlets.fs))) == 0){}else{
      plotDens(live.fs[[idx]],channels = c("FSC-A","FSC-H"),main=paste0("Live: "),
               cex.lab = 2, cex.axis = 2, cex.main = 2)
      lines(singlets.flowD[[idx]]@filter,type="l",lwd=2)
    }
    
    
    # plotDens(singlets.fs[[idx]],channels=c(markers_list["CD117"],markers_list["CD33"]),main=paste0("Singlets: ",basename(description(f)$FIL)),
    #          cex.lab = 2, cex.axis = 2, cex.main = 2)
    # lines(cd33poscd117posB.flowD[[x]]@filter,type="l",lwd=2)
    
    
    plotDens(singlets.fs[[idx]],c(markers_list[grep(paste("HLA-DR|HLADR", collapse = "|"),names(markers_list))],markers_list[grep(paste("CD117|CD117BC",collapse = "|"),names(markers_list))]),main=paste0("Singlets: "),
             cex.lab = 2, cex.axis = 2, cex.main = 2)
    lines(hladrmastcell.flowD[[idx]]@filter,type="l",lwd=2)
    
    
    plotDens(singlets.fs[[idx]], c(markers_list["CD13"],markers_list["CD33"]),main=paste0("Singlets: "),
             cex.lab = 2, cex.axis = 2, cex.main = 2)
    lines(cd13posbcd33neg.flowD[[idx]]@filter,type="l",lwd=2)
    
    plotDens(cd13posbcd33neg.fs[[idx]],c(markers_list[grep("CD45",names(markers_list))],scat.chans["SSC-A"]),main=paste0("CD13+B CD33-: "),
             cex.lab = 2, cex.axis = 2, cex.main = 2, pch=20)
    abline(v=cd45negcd45posD.flowD[[idx]]@gates, lwd=2)
    
    plotDens(cd45negcd45posD.fs[[idx]], c(markers_list[grep(paste("CD117|CD117BC",collapse = "|"),names(markers_list))],markers_list[grep("CD34",names(markers_list))]), main=paste0("CD45-/CD45+ D: "),
             cex.lab = 2, cex.axis = 2, cex.main = 2, pch=20)
    abline(h=cd34negstroma.flowD[[idx]]@gates[2], lwd=2)
    abline(v=cd34negstroma.flowD[[idx]]@gates[1], lwd=2)
    
    plotDens(singlets.fs[[idx]], c(markers_list[grep(paste("HLA-DR|HLADR", collapse = "|"),names(markers_list))],markers_list["CD33"]),main=paste0("Singlets: "),
             cex.lab = 2, cex.axis = 2, cex.main = 2)
    lines(cd33poshladrpos.flowD[[idx]]@filter,type="l", lwd=2)
    
    plotDens(cd33poshladrpos.fs[[idx]], c(markers_list["CD56"],markers_list[grep(paste("CD11b|CD11B",collapse = "|"),names(markers_list))]), main=paste0("CD33+HLA-DR+: "),
             cex.lab = 2, cex.axis = 2, cex.main = 2, pch=20)
    lines(cd11bpos.flowD[[idx]]@filter, type="l",lwd=2)
    
    plotDens(cd11bpos.fs[[idx]], c(markers_list[grep("CD45",names(markers_list))],scat.chans["SSC-A"]), main=paste0("CD11b+: "),
             cex.lab = 2, cex.axis = 2, cex.main = 2, pch=20)
    lines(sscpos.flowD[[idx]]$sscpos.flowD@filter, type="l",lwd=2)
    
    plotDens(sscpos.fs[[idx]], c(markers_list[grep(paste("CD117|CD117BC",collapse = "|"),names(markers_list))],markers_list[grep("CD34",names(markers_list))]), main=paste0("SSC+: "),
             cex.lab = 2, cex.axis = 2, cex.main = 2, pch=20)
    lines(cd34negcd117neg.flowD[[idx]]@filter, type="l",lwd=2)
    
    if(length(grep("FSC-H",colnames(singlets.fs))) == 0){
      
    }else{
      if(length(na.omit(sscpos.fs[[idx]]@exprs)) == 0){#no plot of plot blank if no cells in this population
     
      }else{
        plotDens(sscpos.fs[[idx]],c("FSC-A","FSC-H"),main=paste0("SSC+: "),
                 cex.lab = 2, cex.axis = 2, cex.main = 2, pch=20)
        lines(singmonohighssc.flowD[[idx]]@filter, type = "l", lwd=2)
      }
    
    
    }
    dev.off()
    
  }, .parallel=T)
  }
  #-----------------------------------
  # Cell counts--------
  #-----------------------------------
  
  # names.cell.counts <- c('FCS.files',"Passed flowCut","ALl events","Live","Singlets",
  #                        "HLA-DR- mast cells","CD13+B CD33-","CD45-/CD45+D","CD34- mesenchymal stromacells/fibroblasts",
  #                        "CD33+HLA-DR+","CD11b+","SSC+","CD34-CD117-","Singlets monocyte derived cells high SSC")
  # 
  names.cell.counts <- c('FCS.files',"Singlets",
                         "HLA-DR- mast cells","CD13+B CD33-","CD45-/CD45+D","CD34- mesenchymal stromacells/fibroblasts",
                         "CD33+HLA-DR+","CD11b+","SSC+","CD34-CD117-","Singlets monocyte derived cells high SSC")

  cell.counts <- ldply(1:length(fs), function(idx){
    f <- fs[[idx]]
    # cc.vector <- c(basename(description(f)$FIL),fcut[[idx]]$passed.flowCut,
    #                nrow(fs[[idx]]),nrow(na.omit(live.fs[[idx]]@exprs)),
    #                nrow(na.omit(singlets.fs[[idx]]@exprs)),
    #                hladrmastcell.flowD[[idx]]@cell.count,
    #                cd13posbcd33neg.flowD[[idx]]@cell.count,
    #                cd45negcd45posD.flowD[[idx]]@cell.count,
    #                cd34negstroma.flowD[[idx]]@cell.count,
    #                cd33poshladrpos.flowD[[idx]]@cell.count,
    #                cd11bpos.flowD[[idx]]@cell.count,
    #                sscpos.flowD[[idx]]$sscpos.flowD@cell.count,
    #                cd34negcd117neg.flowD[[idx]]@cell.count,
    #                nrow(na.omit(singmonohighssc.fs[[idx]]@exprs))
    # 
    #                )
    cc.vector <- c(basename(description(f)$FIL),
                   nrow(na.omit(singlets.fs[[idx]]@exprs)),
                   hladrmastcell.flowD[[idx]]@cell.count,
                   cd13posbcd33neg.flowD[[idx]]@cell.count,
                   cd45negcd45posD.flowD[[idx]]@cell.count,
                   cd34negstroma.flowD[[idx]]@cell.count,
                   cd33poshladrpos.flowD[[idx]]@cell.count,
                   cd11bpos.flowD[[idx]]@cell.count,
                   sscpos.flowD[[idx]]$sscpos.flowD@cell.count,
                   cd34negcd117neg.flowD[[idx]]@cell.count,
                   nrow(na.omit(singmonohighssc.fs[[idx]]@exprs))
                   
    )
    names(cc.vector) <- names.cell.counts
    return(cc.vector)
  }, .parallel=T)
  # cell.counts <- merge(tube1fcs_1,cell.counts,by=c('FCS.files'))
  if(days[d]!="Normal"){
    lv2 <- data.frame(sample=c(paste0("V",1:40)),
      AML_status=c(rep("R",20),rep("0",20)))
  }else{
    lv2 <- data.frame(sample=c(paste0("V",41:60)),
      AML_status=rep("Normal",20))
    
  }
  
  cell.counts <- cbind(cell.counts,lv2)
  # save(cell.counts,file=paste0("~/results/Summary/AutomatedGating/Patient_",unique(cell.counts$PatientID),"_Tube",unique(cell.counts$Tube)[1],"_cellcounts.RData"))

  save(cell.counts,file=paste0("~/results/Summary/AutomatedGating/Tube",tubes,days[d],"cellcounts.RData"))
  
  
  
  }
