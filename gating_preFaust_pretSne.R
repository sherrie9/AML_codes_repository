
source("~/code/helperFunc_S.R")
source("~/code/helperFunc.R")
source("~/code/IMPC_Functions.R")
# if(interactive()){
#   diva <- readDIVAFunc()
# }

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
  library('plotrix')
  library('e1071')
  library('flowStats')
  library('ggplot2')
  library('gridExtra')
  library('plyr')
  library('flowPeaks')
  library('dplyr')}

store.allFCS.AML.normal <- readRDS("~/results/Preprocessing/store.allFCS.AML.normal.rds") # preprocessed metadata

scatter_chans <- c("FSC-A","FSC-H","SSC-A")
tubes <- c("001","002","003","004","005")
days <- c("Day22","last before 2nd ind","last before cons","Normal") #only extract these dates
suppressWarnings(dir.create("~/results/Plots/gating_preFaust_pretSNE"))
seed <- 1234

for(t in 2:length(tubes)){
  
  FCS_1 <- store.allFCS.AML.normal[intersect(which(store.allFCS.AML.normal[,'Tube'] == tubes[t]),
                                             grep(paste(days, collapse = "|"), store.allFCS.AML.normal$Status)),]
  
 
  for(d in 1){
  # for(d in 3){
   
  FCS_11 <- FCS_1[which(FCS_1$Status == days[d]),]

  paths <- apply(FCS_11, 1, function(r){
    if(is.na(r['Path3'])){
      path <- paste0(r['Path1'],"/",r['Path2'],"/",r['FCS.files'])
    }else{
      path <- paste0(r['Path1'],"/",r['Path2'],"/",r['Path3'],"/",r['FCS.files'])
    } 
    return(path)
  })
  
  f_list <- preprocessFCS_files(paths,scatter_chans)
  
  lv2 <- data.frame(id=FCS_11[intersect(names(f_list),rownames(FCS_11)),"PatientID"],
                    st=FCS_11[intersect(names(f_list),rownames(FCS_11)),"AML_Status"])
  
  fs.raw <- flowSet(f_list)
  fs <- fs.raw
  fsApply(fs,function(f){return(na.omit(f@parameters@data$desc))})
  if(t==4&&d==2 | t==2&&d==2){
    fs <- fs[-10]
    lv2 <- lv2[-10,]
  }
  if(t==3&&d==2){
    fs <- fs[-9]
    lv2 <- lv2[-9,]
  }
  if(t==5&&d==2){
    fs <- fs[-11]
    lv2 <- lv2[-11,]
  }
  fsApply(fs,function(f){return(na.omit(f@parameters@data$desc))})
  
  scat.chans <- c(grep(colnames(fs),pattern = "FSC*"), grep(colnames(fs),pattern = "SSC*"))
  names(scat.chans) <- colnames(fs)[scat.chans]
  time.loc <- which(tolower(colnames(fs)) == "time"); names(time.loc) <- NULL
  
  fs <- fsApply(fs,function(f) removeMargins(f,chans = scat.chans, verbose = T))
  
  fs <- fsApply(fs, function(f) removeMargins(f,chans = scat.chans,debris = T, neg = T, verbose = T))
  
  fs <- fsApply(fs, function(f) compensate(f,f@description$SPILL))
  
  gf <- globalFrame(fs)
  lgl <- gf$lgl
  
  fs <- fsApply(fs, function(f) transform(f, lgl))
  
  
  fs.preflowCut <- fs
  
  fcut <- fsApply(fs, function(f){
    channels.to.clean <- which(complete.cases(f@parameters@data$desc) == TRUE)
    time.loc <- unique(c(grep("time", tolower(colnames(f)))))
    
    f.Clean <- flowCut(f, Channels = channels.to.clean,
                       Plot = 'None', PrintToConsole = FALSE)
    
    return(list(ind = f.Clean$ind, passed.flowCut = f.Clean$data['Has the file passed',], f=f.Clean$frame))
  })
  
  
  fs <- flowSet(lapply(fcut, function(f){return(f$f)}))
  
  
  # gs <- GatingSet(fs) 
  

  markers_list <- as.character(fs[[1]]@parameters@data[1:ncol(fs[[1]]), 2])
  markers_list <- Find.markers(fs[[1]],na.omit(markers_list))
  channel_names <-  as.character(fs[[1]]@parameters@data[1:ncol(fs[[1]]), 1])
  
  #===============================gating live===================================================
  live.flowD <- llply(as(fs, Class="list"), function(f){
    
    peaks <- getPeaks(f,channel =  c("FSC-A","SSC-A"))
    if(length(peaks$Peaks) == 4){
      peaks$Peaks <- peaks$Peaks[-which(peaks$Peaks < 50000)]
    }
    
    
    if(length(peaks$Peaks)>=3){
      allcuts <- deGate(f, channel = c("FSC-A"),all.cuts = T) #find gates at head
      gt <- min(allcuts)
      if(gt > 90000){
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
      if(gt > 99000){
        gt <- deGate(f, channel = c("FSC-A"), use.upper = T, upper=F) #find gates at head
      }
      
      
    }else{
      gt <- deGate(f, channel = c("FSC-A"), use.upper = T, upper=F) #find gates at head
      
    }
    
    fD <- flowDensity(f, channels = c("FSC-A","SSC-A"), position=c(T,NA), gates = c(gt,NA))
    return(fD)
  }, .parallel=T)
  
  
  live.fs <- flowSet(lapply(live.flowD, function(fD){getflowFrame(fD)}))
  

  #============================gating singlets=================================================
  
  if(length(grep('FSC-H',colnames(live.fs))) == 0){
    print("no FSC-H channels, skipping singlets gating")  
    singlets.fs <- live.fs
  }else{
    singlets.flowD <- llply(as(live.fs, Class = "list"),function(f){
      if(length(grep('FSC-H', colnames(f))) == 0){
        print("no FSC-H channels, skipping singlets gating")  
        return(f)
      }else{
        flowD <- GateSinglets(f,scat.chans,channels.ind)
        if((length(flowD@index)/nrow(f)) < 0.7){
          theta <- pi/5
          chans <- c(grep("FSC-A", colnames(f)), grep("FSC-H", colnames(f)))
          rot <- rotate.data(f, chans = chans, theta = theta)$data
          gate <- deGate(rot, channel = chans[2], use.upper = T, upper = F, alpha = .01)
          flowD <- flowDensity(rot, channels = chans, gates = c(NA, gate), position = c(NA,T) )
          flowD@flow.frame <- rotate.data(flowD@flow.frame, chans = chans, theta = -theta)$data
          flowD@filter <- rotate.data(flowD@filter, chans = chans, theta = -theta)$data
        }
        return(flowD)
      }
    }, .parallel=T)
    
    singlets.fs <- flowSet(lapply(singlets.flowD, function(fD){getflowFrame(fD)}))
  }
  

  
 

  
  chans <- c(markers_list["CD45"],scat.chans["SSC-A"])

  if(days[d] == "Normal"){
    singlets.resample <- randomResampling_Normal(singlets.fs,20)
    singlets.fs <- singlets.resample
    lv2.temp <- data.frame(cbind(sample(c(1:4000),(20-nrow(lv2))),"Normal"))
    colnames(lv2.temp) <- colnames(lv2)
    lv2 <- rbind(lv2,lv2.temp)
  }else{
    if(t==1){
      singlets.resample <- randomResampling(singlets.fs,20,lv2)
    }else{
      singlets.resample <- randomResampling(singlets.fs,20,lv2)
    }
    lv2 <- singlets.resample$lv2
    singlets.fs <- singlets.resample$singlets.fs
  }
  
  # if(days[d] == "Day22"){
  #   save(singlets.fs,
  #        lv2,file=paste0("~/results/GatingSet_tube",tubes[t],"_",days[d],".RData"))
  # }
  
  if(d == 1 ){
    fs.day22 <- singlets.fs
    cd45negsscnegpos.flowD <- llply(as(singlets.fs, Class = "list"), function(f){
      flowD <- flowDensity(f,channels = chans,position = c(F,NA), gate=c(1.5,NA))
      return(flowD)
    })
    # res.fs <- llply(as(1:length(fs), Class="list"), function(i){
    #   fs.day22[[i]]@exprs <- fs[[i]]@exprs[-cd45negsscnegpos.flowD[[i]]@index,] 
    #   return(fs[[i]])
    # })
    res.flowD <- llply(as(singlets.fs, Class = "list"), function(f){
      flowD <- notSubFrame(f,channels = chans,position = c(F,NA), gate=c(1.5,NA))
      return(flowD)
    })
    
    cd45pospossscneg.flowD <- llply(as(res.flowD, Class = "list"), function(f){
      temp <- f@flow.frame
      temp@exprs[,chans[2]] <- temp@exprs[,chans[2]] / 10000
      # flowD <- flowDensity(f,channels = chans,position = c(T,F),gate=c(2.6,50000),ellip.gate=TRUE)
      temp2 <- rotate.data(temp,chans,theta=-0.25)$data
      
      fD <- flowDensity(temp2,chans,position = c(T,NA),gate=c(2.4,NA))
      
      fD@flow.frame <- rotate.data(fD@flow.frame,chans,theta=0.25)$data
      fD@flow.frame@exprs[,chans[2]] <- fD@flow.frame@exprs[,chans[2]] * 10000
      fD@filter <- rotate.data(fD@filter,chans, theta=0.25)$data
      fD@filter[,2] <- fD@filter[,2] * 10000
      return(fD)
      })
    
    res.fs <- llply(as(res.flowD,Class="list"),function(f){return(f@flow.frame)})
 
    res.fs <- llply(as(1:length(res.fs), Class="list"), function(i){
      res.fs[[i]]@exprs <- res.fs[[i]]@exprs[-cd45pospossscneg.flowD[[i]]@index,]
      return(res.fs[[i]])
    })
   

    
    cd45possscpos.flowD <- llply(as(res.fs,Class="list"), function(f){
      temp <- f
      temp@exprs[,chans[2]] <- temp@exprs[,chans[2]] / 10000
      # flowD <- flowDensity(f,channels = chans,position = c(T,F),gate=c(2.6,50000),ellip.gate=TRUE)
      temp2 <- rotate.data(temp,chans,theta=-0.18)$data
      fD <- flowDensity(temp2,chans,position = c(F,NA),gate=c(1.5,NA))
      
      fD@flow.frame <- rotate.data(fD@flow.frame,chans,theta=0.18)$data
      fD@flow.frame@exprs[,chans[2]] <- fD@flow.frame@exprs[,chans[2]] * 10000
      fD@filter <- rotate.data(fD@filter,chans, theta=0.18)$data
      fD@filter[,2] <- fD@filter[,2] * 10000
      
      fD2 <- flowDensity(fD@flow.frame,chans,position = c(NA,T),gate=c(NA,50000))
      # fD2@filter[,2] <- fD2@filter[,2] * 10000
      return(fD2)
    })
    
    cd45possscneg.fs <- llply(as(1:length(res.fs), Class="list"), function(i){
      res.fs[[i]]@exprs <- res.fs[[i]]@exprs[-cd45possscpos.flowD[[i]]@index,]
      return(res.fs[[i]])
    })
    
    cd45possscneg.flowD <- llply(as(cd45possscneg.fs, Class="list"), function(f){
      fD <- flowDensity(f,chans,position = c(T,T), gate=c(0,0))
      return(fD)
    })
    
    # par(mfrow=c(2,3))
    # a <- lapply(1:length(singlets.fs), function(x){
    #   plotDens(singlets.fs[[1]],channels=c(markers_list["CD45"],scat.chans["SSC-A"]))
    #   # lines(cd45possscpos.flowD[[x]]@filter)
    #   points(cd45negsscnegpos.flowD[[1]]@flow.frame@exprs[,chans], pch=".",col="brown")
    #   points(cd45possscpos.flowD[[1]]@flow.frame@exprs[,chans], pch=".",col="green")
    #   points(cd45pospossscneg.flowD[[1]]@flow.frame@exprs[,chans], pch=".",col="red")
    #   points(cd45possscneg.flowD[[1]]@flow.frame@exprs[,chans],pch=".",col="grey")
    #   
    # })
    png(file=paste0("~/results/Plots/gating_preFaust_pretSNE/gating_tubes",tubes[t],"_",days[d]),
        width = 2000*4/4, height = 2000*4/4)
    par(mfrow = c(7,9), mar = (c(5, 5, 4, 2) + 0.1))
    
    
    for(p in 1:length(singlets.fs)){
      plotDens(singlets.fs[[p]],chans, main=paste0(lv2$id[p],"_",lv2$st[p]))
      # plotDens(singlets.fs[[p]],chans, main=paste0(lv2$id[p],"_",lv2$st[p]))
      # points(cd45negsscnegpos.flowD[[p]]$fD@flow.frame@exprs[,chans], pch=".",col="blue")
      points(cd45negsscnegpos.flowD[[p]]@flow.frame@exprs[,chans], pch=".",col="blue")
      points(cd45possscpos.flowD[[p]]@flow.frame@exprs[,chans], pch=".",col="green")
      points(cd45pospossscneg.flowD[[p]]@flow.frame@exprs[,chans], pch=".",col="red")
      points(cd45possscneg.flowD[[p]]@flow.frame@exprs[,chans],pch=".",col="grey")
    }
    
    dev.off()
    save(singlets.fs,
         cd45negsscnegpos.flowD,
         cd45possscpos.flowD,
         cd45possscneg.flowD,
         cd45pospossscneg.flowD,
         lv2,file=paste0("~/results/GatingSet_tube",tubes[t],"_",days[d],".RData"))
  }
  
  
  
    if(days[d]=="Normal" | days[d] == "last before cons" | days[d]=="last before 2nd ind"){
      cd45subpop.flowD <- llply(as(singlets.fs, Class="list"), function(f){
        if(days[d]=="last before 2nd ind"){#TODO MONDAY: remove cells with negative CD45 values
          fD.rmneg <- flowDensity(f,chans,position = c(T,NA),gate=c(0,NA))
          f <- fD.rmneg@flow.frame
          f@exprs <- na.omit(f@exprs)
        }
        
        
        ff <- f
        
        ff@exprs[,chans[2]] <- ff@exprs[,chans[2]] * 0.00001
        
        set.seed(seed)
        flowPeaks.Res <- flowPeaks(ff@exprs[,chans], tol = 0.1, h0=1,h=1.5)
        if(days[d]=="Normal" && length(flowPeaks.Res$peaks$cid) == 3){
          flowPeaks.Res <- flowPeaks(ff@exprs[,chans], tol = 0.08, h0=0.8,h=1.1)
        } 
        if(days[d]=="Normal" && length(flowPeaks.Res$peaks$cid) == 2){
          flowPeaks.Res <- flowPeaks(ff@exprs[,chans], tol = 0.06, h0=0.5,h=1.1)
        } 
        if(t==5 && days[d]=="last before cons" &&
           (basename(f@description$FILENAME)=="13201300004_25130928_d22post2ndind_13201300004_25130928_d22post2ndind_Tube_005.fcs" |
           basename(f@description$FILENAME) == "resample_7.fcs" |
           basename(f@description$FILENAME) =="resample_15.fcs")){
          min.cid <- flowPeaks.Res$peaks$cid[intersect(which(flowPeaks.Res$peaks$mu[,1] < 1),
                                                       which(flowPeaks.Res$peaks$mu[,2] < 1))]
          first.cid <- flowPeaks.Res$peaks$cid[intersect(which(flowPeaks.Res$peaks$mu[,1] < 2),
                                                         which(flowPeaks.Res$peaks$mu[,2] < 0.5))]
        }else{
          min.cid <- flowPeaks.Res$peaks$cid[intersect(which(flowPeaks.Res$peaks$mu[,1] < 2),
                                                       which(flowPeaks.Res$peaks$mu[,2] < 1.9))]
          if(d==1){
            min.cid <- flowPeaks.Res$peaks$cid[intersect(which(flowPeaks.Res$peaks$mu[,1] < 1.5),
                                                         which(flowPeaks.Res$peaks$mu[,2] < 1.5))]
          }
          first.cid <- flowPeaks.Res$peaks$cid[intersect(which(flowPeaks.Res$peaks$mu[,1] < 2),
                                                         which(flowPeaks.Res$peaks$mu[,2] < 0.5))]
        }
        
        # min.cid <- min.cid[which.min(flowPeaks.Res$peaks$mu[min.cid,2])]
        
        
        if(length(first.cid)== 0){
          peaks.cd45 <- getPeaks(f,chans[1])  
          if(length(peaks.cd45$Peaks) ==1){
            gate.cd45 <- deGate(f,chans[1],use.upper = T, upper = F, alpha= 0.0001)
            if(gate.cd45 < 1.5){
              gate.cd45 <- deGate(f,chans[1],use.upper = T, upper = F, alpha= 0.01)
            }
          }
          if(length(peaks.cd45$Peaks) >= 2){
            if(min(peaks.cd45$Peaks) > 2){
              gate.cd45 <- deGate(f,chans[1],use.upper = T, upper = F)
            }else{
              gate.cd45 <- deGate(f,chans[1])
              temp1 <- flowDensity(f,chans,position = c(F,NA),gate=c(gate.cd45,NA))
              gate.cd45 <- deGate(temp1,chans[1],use.upper = T, upper = T)
            }
            
          }
          
          if(gate.cd45 < 1){
            gate.cd45  <- 2
          }
          if(gate.cd45 > 2.5){#or at least for day22
            gate.cd45 <- 2.5
          }
          
          fD <- flowDensity(f,chans,position = c(F,NA), gates=c(gate.cd45,NA))
          first.pop.flowD <- fD
          notfD <- notSubFrame(f,chans,position = c(F,NA), gates=c(gate.cd45,NA))
          f <- notfD@flow.frame
        }else{
          first.pop <- f
          first.pop@exprs <- f@exprs[which(flowPeaks.Res$peaks.cluster %in% min.cid), ]
          
          first.pop.flowD <- flowDensity(first.pop, chans, position = c(T, T), gates = c(0,0))
          
          
        }
        
        
        res.cid <- setdiff(flowPeaks.Res$peaks$cid,min.cid)
        big.cid <- res.cid[intersect(which(flowPeaks.Res$peaks$mu[res.cid,1] > 2.8),
                           which(flowPeaks.Res$peaks$mu[res.cid,1] < 4.3))]
        lymph.cid <- big.cid[which.min(flowPeaks.Res$peaks$mu[big.cid,2])]
        
        lymph3 <- f
        lymph3@exprs <- f@exprs[which(flowPeaks.Res$peaks.cluster %in% lymph.cid), ]
        
        lymph3.flowD <- flowDensity(lymph3, chans, position = c(T, T), gates = c(0,0))
        
        res.cid <- setdiff(flowPeaks.Res$peaks$cid,c(min.cid,lymph.cid))
        
        gran.cid <- res.cid[intersect(which(flowPeaks.Res$peaks$mu[res.cid,2] > 0.75),
                                                      which(flowPeaks.Res$peaks$mu[res.cid,1] < 3.89))]
        gran <- f
        gran@exprs <- f@exprs[which(flowPeaks.Res$peaks.cluster %in% gran.cid), ]
        
        gran.flowD <- flowDensity(gran, chans, position = c(T, T), gates = c(0,0))
        
       
        cd45possscneg.cid  <- setdiff(flowPeaks.Res$peaks$cid, c(lymph.cid,gran.cid,min.cid))
        cd45possscneg <- f
        cd45possscneg@exprs <- matrix(f@exprs[which(flowPeaks.Res$peaks.cluster %in% cd45possscneg.cid), ],
                                      ncol=ncol(f@exprs))
        colnames(cd45possscneg@exprs) <- colnames(f@exprs)
        
        cd45possscneg.flowD <- flowDensity(cd45possscneg, chans, position = c(T, T), gates = c(0,0))
        
        return(list(cd45negsscnegpos.flowD = first.pop.flowD, cd45possscpos.flowD = gran.flowD,
                    cd45pospossscneg.flowD = lymph3.flowD,
                    cd45possscneg.flowD = cd45possscneg.flowD, flowPeaks.Res=flowPeaks.Res))
      })
      
      
    
    
    

    # par(mfrow=c(2,3))
    # a <- lapply(1:length(cd45subpop.flowD), function(x){
    #   plotDens(singlets.fs[[x]],channels=c(markers_list["CD45"],scat.chans["SSC-A"]))
    #   points(cd45subpop.flowD[[x]]$cd45negsscnegpos.flowD@flow.frame@exprs[,chans], pch=".",col="blue")
    #   points(cd45subpop.flowD[[x]]$cd45possscpos.flowD@flow.frame@exprs[,chans], pch=".",col="green")
    #   points(cd45subpop.flowD[[x]]$cd45pospossscneg.flowD@flow.frame@exprs[,chans], pch=".",col="red")
    #   points(cd45subpop.flowD[[x]]$cd45possscneg.flowD@flow.frame@exprs[,chans],pch=".",col="grey")
    # 
    # })

    cd45negsscnegpos.flowD <- llply(as(cd45subpop.flowD, Class="list"),function(fD){
      return(fD$cd45negsscnegpos.flowD)
    })

    cd45possscpos.flowD <- llply(as(cd45subpop.flowD, Class = "list"), function(fD){
      return(fD$cd45possscpos.flowD)
    })

    cd45pospossscneg.flowD <- llply(as(cd45subpop.flowD, Class = "list"), function(fD){
      return(fD$cd45pospossscneg.flowD)
    })

    cd45possscneg.flowD <- llply(as(cd45subpop.flowD, Class = "list"), function(fD){
      return(fD$cd45possscneg.flowD)
    })
    
  
  
  #------------------------plotting----------------------------------------------------
  png(file=paste0("~/results/Plots/gating_preFaust_pretSNE/gating_tubes",tubes[t],"_",days[d]),
      width = 2000*4/4, height = 2000*4/4)
  par(mfrow = c(7,9), mar = (c(5, 5, 4, 2) + 0.1))


  for(p in 1:length(singlets.fs)){
    plotDens(singlets.fs[[p]],chans, main=paste0(lv2$id[p],"_",lv2$st[p]))
    # plotDens(singlets.fs[[p]],chans, main=paste0(lv2$id[p],"_",lv2$st[p]))
    # points(cd45negsscnegpos.flowD[[p]]$fD@flow.frame@exprs[,chans], pch=".",col="blue")
    points(cd45subpop.flowD[[p]]$cd45negsscnegpos.flowD@flow.frame@exprs[,chans], pch=".",col="blue")
    points(cd45subpop.flowD[[p]]$cd45possscpos.flowD@flow.frame@exprs[,chans], pch=".",col="green")
    points(cd45subpop.flowD[[p]]$cd45pospossscneg.flowD@flow.frame@exprs[,chans], pch=".",col="red")
    points(cd45subpop.flowD[[p]]$cd45possscneg.flowD@flow.frame@exprs[,chans],pch=".",col="grey")
  }

  dev.off()

  # png(file=paste0("~/results/Plots/gating_preFaust_pretSNE/tubes",tubes[t],"_",days[d]),
  #     width = 2000*4/4, height = 2000*4/4)
  # par(mfrow = c(4,4), mar = (c(5, 5, 4, 2) + 0.1))
  # 
  # 
  # for(p in 1:length(singlets.fs)){
  #   plotDens(singlets.fs[[p]],chans, main=paste0(lv2$id[p],"_",lv2$st[p]))
  # 
  # }
  # 
  # dev.off()
  
  save(singlets.fs,
       cd45negsscnegpos.flowD,
       cd45possscpos.flowD,
       cd45possscneg.flowD,
       cd45pospossscneg.flowD,
       lv2,file=paste0("~/results/GatingSet_tube",tubes[t],"_",days[d],".RData"))
  
    }
  
  } 
 
}
