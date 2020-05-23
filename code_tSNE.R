remove(list=ls())
# setwd("/data/projects/Codes/MRD")



source("~/code/helperFunc_S.R")
source("~/code/helperFunc.R")
source("~/code/IMPC_Functions.R")

if(interactive()){
  center <- readCenterFunc()
}

addDay0 <- FALSE

library('flowCore')
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

inputPath <- paste0("~/data/from ",center)
inputPathNormal <- list(c("~/data/Normal BM DIVA 6"),c("~/data/Normal BM DIVA 8"))

suppressWarnings(dir.create("~/results/Plots/"))
suppressWarnings(dir.create("~/results/Plots/flowCut/"))
suppressWarnings(dir.create("~/results/Plots/tSNE/"))
suppressWarnings(dir.create("~/results/Summary/"))
suppressWarnings(dir.create("~/results/tSNE_Rdata"))
suppressWarnings(dir.create("~/results/Plots/normalization"))
suppressWarnings(dir.create("~/results/Summary/Cleaned_Files"))
suppressWarnings(dir.create("~/results/Plots/tSNE_normal"))
suppressWarnings(dir.create(paste0("~/results/Plots/tSNE_",center)))

flowCutdir <- "~/results/Plots/flowCut/"

store.allFCS <- preProcessFCS_center(inputPath,center)
store.allFCS.normal <- preProcessFCS_normal(inputPathNormal)


Patients <- unique(store.allFCS$PatientID)
Days <- c("Diagnos","Diagnosis","Day0","Day22",
          "day0","day22","day29","1HAM","last before 2nd induction","last before 2nd ind","last before consolidation","d22post2ndind",
          "d22","Last bef 2nd ind","Last bef consolidation","before Ind2","before cons1")
markersv <- c("CD34", "CD117","HLADR","HLA-DR","CD13","CD33","CD56","NG2","CD7","CD123","CD38","CD99","CD133")

for(x in 1:length(Patients)){

  
  if(addDay0){
    FCS_1 <- store.allFCS[which(store.allFCS[,'PatientID'] == Patients[x]), ]
  }else{
    FCS_1 <- store.allFCS[which(store.allFCS[,'PatientID'] == Patients[x]), ]
    if(length(grep(paste("Diagnos|Day0|day0",collapse = "|"), FCS_1[,"Status"]))>0){
      FCS_1 <- FCS_1[-grep(paste("Diagnos|Day0|day0",collapse = "|"), FCS_1[,"Status"]),]
    }else{
      print("There is no day 0 for this patient")
    }
   
  }
  
  
  
  
  
  tube <- unique(FCS_1[,'Tube'])
  
  
  for(i in 1:length(tube)){
    
    

    FCS_11 <- FCS_1[which(FCS_1[,'Tube'] == tube[i]),]
    
    if(nrow(FCS_11) == 1){
      next
    }
    
    days <- unique(FCS_11$Status)
    
    if(addDay0){
     
      days <- intersect(days, Days)
      FCS_11 <-FCS_11[unlist(sapply(days,function(h){grep(h,FCS_11[,"Status"])})),]
      if(center=="Israel"){
        if(length(grep(paste("1w|2w", collapse = "|"), FCS_11[,"Status"])) > 0){
          FCS_11 <-FCS_11[-grep(paste("1w|2w", collapse = "|"), FCS_11[,"Status"]),]
        }
        
        
      }
      
    }
    
    # fs.raw <- makefs(FCS_11)
    # fs.raw <- read.flowSet(paste0(FCS_11$Path1,"/",FCS_11$Path2,"/",FCS_11$Path3,"/",FCS_11$FCS.files))
    
    f_list <- list()
 
    for(m in 1:nrow(FCS_11)){
      path <- FCS_11[m,]
      f <- read.FCS(filename = paste0(path$Path1,"/",path$Path2,"/",path$Path3,"/",path$FCS.files))
  
      scat.chans <- c(grep(colnames(f),pattern = "FSC*"), grep(colnames(f),pattern = "SSC*"))
      names(scat.chans) <- colnames(f)[scat.chans]
      time.loc <- which(tolower(colnames(f)) == "time"); names(time.loc) <- NULL
      
      f <-  removeMargins(f,chans = scat.chans, verbose = T)
      
      f <- removeMargins(f,chans = scat.chans,debris = T, neg = T, verbose = T)
      
      # if(dim(f@description$SPILL)[1] != dim(f@description$SPILL)[2]){
      #   print(paste0("file ",path$FCS.files,"spill matrix is not square, skipping compensation.."))
      # }else{
        f <-compensate(f,f@description$SPILL)
      # }
      

      f_list[[m]] <- f
    }
  
    
    f_list <- makefs(f_list)
    
    fs <- flowSet(f_list)
    
    # scat.chans <- c(grep(colnames(fs),pattern = "FSC*"), grep(colnames(fs),pattern = "SSC*"))
    # names(scat.chans) <- colnames(fs)[scat.chans]
    # time.loc <- which(tolower(colnames(fs)) == "time"); names(time.loc) <- NULL
    # 
    # fs <- fsApply(fs,function(f) removeMargins(f,chans = scat.chans, verbose = T))
    # 
    # fs <- fsApply(fs, function(f) removeMargins(f,chans = scat.chans,debris = T, neg = T, verbose = T))
    # 
    # fs <- fsApply(fs, function(f) compensate(f,f@description$SPILL))
    
    gf <- globalFrame(fs)
    lgl <- gf$lgl
    
    fs <- fsApply(fs, function(f) transform(f, lgl))
    
    
    fs.preflowCut <- fs
    
    
    
    fcut <- fsApply(fs, function(f){
      channels.to.clean <- which(complete.cases(f@parameters@data$desc) == TRUE)
      time.loc <- unique(c(grep("time", tolower(colnames(f)))))
      # png(file=paste0("~/results/Plots/flowCut/",identifier(f),".png"),
      #     width = 5*300, height=2*300)
      
      f.Clean <- flowCut(f, Channels = channels.to.clean,Directory = flowCutdir,
                         FileID = identifier(f),Plot = 'None', PrintToConsole = FALSE)
      
      # dev.off()

      # if(!file.exists("~/results/Summary/flowCut_Summary.csv")){
      #   header <- c("File_Name", "Flag","Events_Removed(%)")
      #   write.csv((t(header)), row.names = FALSE, file="~/results/Summary/flowCut_Summary.csv")
      # }
      # 
      # write.table(t(c(identifier(f),f.Clean$data[17,1], f.Clean$data[13,1])),
      #             row.names = FALSE, file="~/results/Summary/flowCut_Summary.csv", sep=",",
      #             col.names = FALSE, append = TRUE)
      
      return(list(ind = f.Clean$ind, passed.flowCut = f.Clean$data['Has the file passed',]))
    })
    
    
    #--------write cleaned fs.raw to file-------------------------------------
    # fs.raw <- fsApply(fs.raw, function(f){
    #   if(length(fcut[[which(FCS_11$FCS.files == basename(description(f)$FIL))]]$ind) > 0){
    #     f@exprs <- f@exprs[-fcut[[which(FCS_11$FCS.files == basename(description(f)$FIL))]]$ind,] 
    #     write.FCS(f,filename = paste0 ("~/results/Summary/Cleaned_Files/",identifier(f)))
    #   }
    #   return(f)
    # })
    #-----------------------------------------
    
    
    fs <- fsApply(fs.preflowCut, function(f){
      if(length(fcut[[which(FCS_11$FCS.files == basename(description(f)$FIL))]]$ind) > 0){
        f@exprs <- f@exprs[-fcut[[which(FCS_11$FCS.files == basename(description(f)$FIL))]]$ind, ]
      }
      return(f)
    })
    markers_list <- as.character(fs[[1]]@parameters@data[1:ncol(fs[[1]]), 2])
    
    markers_list <- Find.markers(fs[[1]],na.omit(markers_list))
    #------------------------------------Gating live-----------------------------------------------
    
    live.flowD <- fsApply(fs, function(f){
      fD <- flowDensity(f, channels = c("FSC-A","SSC-A"), position=c(T,NA), gates = c(60000,NA))
    })
    
    par(mfrow=c(2,3))
    a <- lapply(1:length(live.flowD), function(x){
      plotDens(fs[[x]],channels=c("FSC-A","SSC-A"))
      lines(live.flowD[[x]]@filter)
    })
    
    live.fs <- flowSet(lapply(live.flowD, function(fD){fD@flow.frame}))
    
    # live.flowD <- fsApply(live.fs, function(f){
    #   fD <- flowDensity(f, channels = c("FSC-H","FSC-A"), position=c(NA,T), gates = c(NA,50000))
    # })
    # 
    # par(mfrow=c(2,3))
    # a <- lapply(1:length(live.flowD), function(x){
    #   plotDens(fs[[x]],channels=c("FSC-H","FSC-A"))
    #   lines(live.flowD[[x]]@filter)
    # })
    # 
    # live.fs <- flowSet(lapply(live.flowD, function(fD){fD@flow.frame}))
    
    #----------------------------Gating Singlets--------------------------------------------------------------
    
    if(any(grepl("FSC-H",names(scat.chans)))){
      singlets.flowD <- fsApply(live.fs,function(f){
        flowD <- GateSinglets(f,scat.chans,channels.ind)
        return(flowD)
      })
      
      par(mfrow=c(2,3))
      a <- lapply(1:length(singlets.flowD), function(x){
        plotDens(live.fs[[x]],channels=c("FSC-A","FSC-H"))
        lines(singlets.flowD[[x]]@filter)
      })
      
      singlets.fs <- flowSet(lapply(singlets.flowD, function(fD){fD@flow.frame}))
    }else{
      
      # singlets.fs <- live.fs
      singlets.fs <- fsApply(live.fs,function(f){
        f@exprs <- f@exprs[-which(is.na(f@exprs[,1])),]
        return(f)
      })
    }
    
    
    
    # #----------------------------------------------------------------------------------------------------------
    # #-----------------------------normalizing-------------------------------------------------
    # 
    # phData <- data.frame(FCS_111[,c("Status","PatientID","Tube")])
    # 
    # rownames(phData) <- FCS_111[,"FCS.files"]
    # 
    # pData(singlets.fs) <- phData
    # 
    # channels.to.clean <- which(complete.cases(singlets.fs[[1]]@parameters@data$desc) == TRUE)
    # pars <- colnames(fs)[channels.to.clean]
    # 
    # norm <- normalization(normFunction  = function(x, parameters, ...)
    #   warpSet(x, stains=parameters,...),
    #   parameters =  pars,
    #   normalizationId = "Warping"
    # )
    # 
    # wf <- workFlow(singlets.fs)
    # 
    # add(wf, norm)
    # 
    # # dev.off()
    # png(file=paste0("~/results/Plots/normalization/Patients",x,"_Tube00",i,"_Diva",
    #                 diva,"_BeforeNormalization.png"), width = 4*300, height = 8*300)
    # par(mfrow=c(8,4), mar=c(15,15,15,16))
    # 
    # densityplot(Status~., singlets.fs, channels=pars, main=("before normalization"),
    #             filter=lapply(par,curv1Filter)) #before normalizing
    # 
    # dev.off()
    # 
    # png(file=paste0("~/results/Plots/normalization/Patients",x,"_Tube00",i,"_Diva",
    #                 diva,"_AfterNormalization.png"), width = 4*300, height = 8*300)
    # par(mfrow=c(8,4), mar=c(15,15,15,16))
    # densityplot(Status~.,Data(wf[["Warping"]]), channels=pars, main=("after normalization"),
    #             filter=lapply(par,curv1Filter))
    # 
    # 
    # dev.off()
    # 
    # singlets.fs <- Data(wf[["Warping"]])
    
    #--------------------------------------run tSNE------------------------------------------------------------
    set.seed(1234)
    channels.ind <- setdiff(1:length(colnames(fs[[length(singlets.fs)]])), c(scat.chans,time.loc))
    
    fb_result <- fb_result_function(singlets.fs,8000)
    
    
    #-------------------flowSOM------------------------------------------------------
    data_flowSOM <- flowFrame(fb_result)
    out <- FlowSOM::ReadInput(data_flowSOM, transform = FALSE, scale=FALSE)
    out <- FlowSOM::BuildSOM(out,colsToUse = markers_list)
    out <- FlowSOM::BuildMST(out)
    
    FlowSOM::PlotStars(out)
    
    labels_pre <- out$map$mapping[, 1]
    
    k <- 40
    
    seed <- 1234
    out <- FlowSOM::metaClustering_consensus(out$map$codes, k = k, seed = seed)
    
    labels <- out[labels_pre]
    
    table(labels)
    
    res <- data.frame(cluster = labels)
    # fb_result_fs <- flowSet(flowFrame(fb_result[c(1:6000),]),flowFrame(fb_result[c(6001:12000),]),
    #                         flowFrame(fb_result[c(12001:18000),]),flowFrame(fb_result[c(18001:24000),]))
    
    lv <- NULL
    for(j in 1:length(days)){
      lv <- c(lv,rep(j,8000))
    }
    
    fb_result_f <- data.frame(cbind(fb_result,lv))
    
    
    # if(file.exists(file=paste0("~/results/tSNE_Rdata/Patient",Patients[x],"_tube",tube[i],
    #                            ifelse(addDay0, yes="_withDay0", no="_noDay0"), "_Rtsne.RData"))){
    #   
    #   load(file=paste0("~/results/tSNE_Rdata/Patient",Patients[x],"_tube",tube[i],
    #                    ifelse(addDay0, yes="_withDay0", no="_noDay0"), "_Rtsne.RData"))
    # }else{
      set.seed(1234)
      tsne_out <- Rtsne(fb_result[,channels.ind])
      print("end of tSNE")
      save(tsne_out, file=paste0("~/results/tSNE_Rdata/Patient",Patients[x],"_tube",tube[i],ifelse(addDay0, yes="_withDay0", no="_noDay0"), "_Rtsne.RData"))

    # }
    
    
    # load(file=paste0("~/results/tSNE_Rdata/Patient",Patients[x],"_tube",tube[i],ifelse(addDay0, yes="_withDay0", no="_noDay0"), "_Rtsne.RData"))
    
    data_plot <- tsne_out$Y
    
    #-----------------------------------flowSOM tSNE overlay--------------------------------
    
    
    
    #-------------------------------------normalizing marker expression for each day-----------
    
    # phData <- data.frame(FCS_111[,c("Status","PatientID","Tube")])
    # 
    # rownames(phData) <- rownames(pData(fb_result_fs))
    # 
    # pData(fb_result_fs) <- phData
    # 
    # # channels.to.clean <- which(complete.cases(singlets.fs[[1]]@parameters@data$desc) == TRUE)
    # pars <- colnames(fb_result_fs)
    # 
    # norm <- normalization(normFunction  = function(x, parameters, ...)
    #   warpSet(x, stains=parameters,...),
    #   parameters =  pars,
    #   normalizationId = "Warping"
    # )
    # 
    # wf <- workFlow(fb_result_fs)
    # 
    # add(wf, norm)
    # 
    # fb_result <- rbind(Data(wf[["Warping"]])[[1]]@exprs,
    #                    Data(wf[["Warping"]])[[2]]@exprs,
    #                    Data(wf[["Warping"]])[[3]]@exprs,
    #                    Data(wf[["Warping"]])[[4]]@exprs
    # )
    # 
    # fb_result_f <- data.frame(cbind(fb_result,lv))
    #------------------------------------- plot by marker--------------------------------------------------------------
    
    # z1 <- ceiling(sqrt(length(markers_list) + 2))
    
    
    if(addDay0){
      
      # if(length(intersect(names(markers_list), markersv)) == 0){
      #   next
      # }else{
      #   markers_list <-Find.markers(fs[[1]],intersect(names(markers_list),markersv))
      # }
      if(center == "copenhagen"){
        markers_list <- markers_list
      }else{
        markers_list <-Find.markers(fs[[1]],intersect(names(markers_list),markersv))
      }
     
    }
    markers_to_see <- c(scat.chans,markers_list)
    png(file=paste0("~/results/Plots/tSNE_",center,"/Patient_",Patients[x],"_allMarkers_Tube",tube[i],ifelse(addDay0, yes="_withDay0", no="_noDay0"),".png"),
        width = 2*300, height= ceiling(length(markers_to_see)/2) *300)
    par(mfrow=c(ceiling(length(markers_to_see)/2),2), mar=c(5,5,5,6))
    
   
    
    for(z in 1:length(markers_to_see)){
      flowresults_marker <- fb_result[,z]
      tsnecol <- c(NA)
      flowresults_scale <- c(NA)
      
      if(length(which(flowresults_marker < 0)) > 0){
        flowresults_marker[which(flowresults_marker < 0)] <- 0
      }
      
      for (j in 1:length(flowresults_marker)) {
        flowresults_scale[j] <- flowresults_marker[j] - min(flowresults_marker)
      }
      max <- max(flowresults_scale)
      
      for (j in 1:length(flowresults_scale)) {
        flowresults_scale[j] <- flowresults_scale[j]/max
      }
      
      # colorScale <- NULL
      # colorScale[which(flowresults_marker_temp < 0.2)] <- 0
      # colorScale[intersect(which(flowresults_marker_temp>=0.2),which(flowresults_marker_temp <0.4))] <- 1
      # colorScale[intersect(which(flowresults_marker_temp>=0.4),which(flowresults_marker_temp <0.6))] <- 2
      # colorScale[intersect(which(flowresults_marker_temp>=0.6),which(flowresults_marker_temp <0.8))] <- 3
      # colorScale[which(flowresults_marker_temp >= 0.8)] <- 4
      # 
      NumEvents <- nrow(fb_result)
      ColourPal <- rainbow(round(NumEvents*3/2))
      
      
      markernames <- colnames(fs[[1]]@exprs)
      tsnecol <- color.scale(flowresults_scale,
                             extremes=c("blue1", "green", "yellow", "orange", "red"), color.spec="rgb")
      tsnecol2 <- color.scale(1:length(flowresults_scale), extremes=c("blue1", "green", "yellow", "orange", "red"), color.spec="rgb")
      
      plot(data_plot, col=tsnecol, main=paste0(markernames[markers_to_see[z]], " ", names(markers_to_see)[z]),
           xlab="tsne1",ylab="tsne2", pch=19, cex=0.5, cex.lab=2, cex.main=2,
           cex.axis=2)
      
      image.plot( legend.only=TRUE, zlim= range(fb_result[,z]), col=tsnecol2[seq(1,NumEvents,by=1)],
                  legend.width=2, legend.shrink=1.0, legend.mar=8, nlevel=400, axis.args = list(cex.axis = 2))
      
      # d_tsne<- as.data.frame(cbind(data_plot, colorScale))
      # 
      # p <- ggplot(d_tsne, aes(x=V1,y=V2)) +
      #   geom_point(size=0.5, aes(colour= factor(colorScale))) +
      #   scale_color_manual(values=c("blue1","green","yellow","orange","red"))+
      #   guides(colour=guide_legend(override.aes=list(size=6), title = "Expression")) +
      #   xlab("") + ylab("") +
      #   ggtitle(paste0(markernames[markers_list[k]], " ",
      #                  names(markers_list)[k],"_",day[unique(fb_result_f$lv)[j]]))+
      #   # theme_light(base_size =20)+ 
      #   theme_bw()+
      #   theme(axis.text.x=element_blank(),
      #         axis.text.y=element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()
      #   )
      
    }
    dev.off()
    
    
    #------------------------
    
    if(length(fs)==1){
      next
    }
    
    #---------------------------------------plotting marker progression----------------------------
    
    
     png(file=paste0("~/results/Plots/tSNE_",center,"/Patient_",Patients[x],"_MarkersProgression_Tube",tube[i],ifelse(addDay0, yes="_withDay0", no="_noDay0"),".png"),
        width = length(days)*300, height=length(markers_to_see)*300)
    # par(mfrow=c(8,4), mar=c(5,5,5,6))
    
    plots <- list()
    u <- 0
    for(k in 1:length(markers_to_see)){
      
      u <- u +1
      flowresults_marker <- fb_result[,k]
      
      flowresults_scale <- c(NA)
      
      for (j in 1:length(flowresults_marker)) {
        flowresults_scale[j] <- flowresults_marker[j] - min(flowresults_marker)
      }
      max <- max(flowresults_scale)
      
      for (j in 1:length(flowresults_scale)) {
        flowresults_scale[j] <- flowresults_scale[j]/max
      }
      
      
      
      for(j in 1:length(unique(fb_result_f$lv))){
        
        
        flowresults_marker_temp <- flowresults_scale
        
        flowresults_marker_temp[which(fb_result_f$lv != unique(fb_result_f$lv)[j])] <- 0
        
        colorScale <- NULL
        colorScale[which(flowresults_marker_temp < 0.2)] <- 0
        colorScale[intersect(which(flowresults_marker_temp>=0.2),which(flowresults_marker_temp <0.4))] <- 1
        colorScale[intersect(which(flowresults_marker_temp>=0.4),which(flowresults_marker_temp <0.6))] <- 2
        colorScale[intersect(which(flowresults_marker_temp>=0.6),which(flowresults_marker_temp <0.8))] <- 3
        colorScale[which(flowresults_marker_temp >= 0.8)] <- 4
        
        # tsnecol <- color.scale(colorScale, extremes = c("grey","yellow","orange","red","dark red"))
        # tsnecol2 <- color.scale(1:length(flowresults_scale), extremes=c("grey","yellow","orange","red","dark red"), color.spec="rgb")
        # 
        # plot(data_plot, col=tsnecol, main=paste0(markernames[markers_list[k]], " ",
        #                                          names(markers_list)[k],"_",day[unique(fb_result_f$lv)[j]]),
        #      xlab="tsne1",ylab="tsne2", pch=19, cex=0.5, cex.lab=2, cex.main=2, cex.axis=2)
        # 
        # 
        # image.plot( legend.only=TRUE, zlim= range(flowresults_marker), col=tsnecol2[seq(1,NumEvents,by=1)],
        #             legend.width=2, legend.shrink=1.0, legend.mar=8, nlevel=400, axis.args = list(cex.axis = 2))
        
        
        
        
        
        
        d_tsne<- as.data.frame(cbind(data_plot, colorScale))
        
        p <- ggplot(d_tsne, aes(x=V1,y=V2)) +
          geom_point(size=0.5, aes(colour= factor(colorScale))) +
          scale_color_manual(values=c("grey","yellow","orange","red","dark red"))+
          guides(colour=guide_legend(override.aes=list(size=6), title = "Expression")) +
          xlab("") + ylab("") +
          ggtitle(paste0(markernames[markers_to_see[k]], " ",
                         names(markers_to_see)[k],"_",days[unique(fb_result_f$lv)[j]]))+
          # theme_light(base_size =20)+ 
          theme_bw()+
          theme(axis.text.x=element_blank(),
                axis.text.y=element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()
          )
        
        plots[[u]] <- p
        if(j == length(unique(fb_result_f$lv))){
          
        }else{
          u <- u + 1
        }
        
      }
      
      
      
      
    }
    
    
    do.call("grid.arrange", c(plots, ncol= length(days)))
    
    
    dev.off()
    
    
    
    
    }

  
}


#===========================================================================================
# the code below runs normal files alone and with patients files from last before 2nd induction
#===========================================================================================


# normalsamples <- unique(store.allFCS.normal$SampleID)
patientslb2ind <- store.allFCS[grep(paste("last before consolidation|d22post2ndind", collapse = "|"),store.allFCS$Status),]
tube_normal <- intersect(unique(store.allFCS.normal$Tube),unique(patientslb2ind$Tube))

# for(j in 1:length(normalsamples)){
#   
#   FCS_1n <- store.allFCS.normal[which(store.allFCS.normal[,"SampleID"] == normalsamples[j]),]
#   
#   tube_normal <- unique(FCS_1n$Tube)
  
for(t in 3:length(tube_normal)){
    patientslb2ind_1 <- as.data.frame(patientslb2ind[which(patientslb2ind[,"Tube"]== tube_normal[t]),])
    
    if(length(which(duplicated(patientslb2ind_1[,"PatientID"]) == "TRUE")) > 0 ){
      
      patientslb2ind_1 <- patientslb2ind_1[-c(which(duplicated(patientslb2ind_1[,"PatientID"]) == "TRUE")-1),]
      
    } 
    FCS_1n <- store.allFCS.normal[which(store.allFCS.normal[,"Tube"] ==tube_normal[t]),]
    
    f_normal <- apply(FCS_1n,1 , function(x){
      f <- read.FCS(filename = paste0(x["Path1"],"/",x["Path2"],"/",x["FCS.files"]))
      # f <- read.FCS(filename = paste0(x["Path1"],"/",x["Path2"],"/",x["Path3"], "/", x["FCS.files"]))
      colnames(f)[which(colnames(f) == "Horizon V450-A")] <- "V450-A"
      colnames(f)[which(colnames(f) == "Horizon V500-A")] <- "V500-A"
      colnames(f@description$SPILL)[which(colnames(f@description$SPILL) =="Horizon V450-A" )] <- "V450-A"
      colnames(f@description$SPILL)[which(colnames(f@description$SPILL) =="Horizon V500-A" )] <- "V500-A"
      return(f)
    })
    
for(x in 1:nrow(patientslb2ind_1)){
    # for(x in 1){  
    
    f_pt <- read.FCS(filename = paste0(patientslb2ind_1[x,]$Path1,"/",patientslb2ind_1[x,]$Path2,"/",
                                       patientslb2ind_1[x,]$Path3,"/",patientslb2ind_1[x,]$FCS.files))
    colnames(f_pt)[which(colnames(f_pt) == "Horizon V450-A")] <- "V450-A"
    colnames(f_pt)[which(colnames(f_pt) == "Horizon V500-A")] <- "V500-A"
    colnames(f_pt@description$SPILL)[which(colnames(f_pt@description$SPILL) =="Horizon V450-A" )] <- "V450-A"
    colnames(f_pt@description$SPILL)[which(colnames(f_pt@description$SPILL) =="Horizon V500-A" )] <- "V500-A"
    
    
  #-------------------------------------  
    fs <- f_normal
    fs[[length(f_normal) + 1]] <- f_pt
    names(fs)[length(fs)] <- as.character(patientslb2ind_1[x,"PatientID"])
    fs <- flowSet(fs)
    
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
      # png(file=paste0("~/results/Plots/flowCut/",identifier(f),".png"),
      #     width = 5*300, height=2*300)
      
      f.Clean <- flowCut(f, Channels = channels.to.clean,Directory = flowCutdir,
                         FileID = identifier(f),Plot = 'None', PrintToConsole = FALSE)
      
      # dev.off()
      
      # if(!file.exists("~/results/Summary/flowCut_Summary.csv")){
      #   header <- c("File_Name", "Flag","Events_Removed(%)")
      #   write.csv((t(header)), row.names = FALSE, file="~/results/Summary/flowCut_Summary.csv")
      # }
      # 
      # write.table(t(c(identifier(f),f.Clean$data[17,1], f.Clean$data[13,1])),
      #             row.names = FALSE, file="~/results/Summary/flowCut_Summary.csv", sep=",",
      #             col.names = FALSE, append = TRUE)
      
      return(list(ind = f.Clean$ind, passed.flowCut = f.Clean$data['Has the file passed',],
                  frame = f.Clean$frame))
    })
    
    
    fs <- lapply(fcut, function(f){f <- f$frame 
      return(f)
      })
    
    fs<- flowSet(fs)
    
    markers_list <- as.character(fs[[1]]@parameters@data[1:ncol(fs[[1]]), 2])
    
    markers_list <- Find.markers(fs[[1]],na.omit(markers_list))
  #-----------------------------------------GATING LIVE--------------------------------------------  
    live.flowD <- fsApply(fs, function(f){
      fD <- flowDensity(f, channels = c("FSC-A","SSC-A"), position=c(T,NA), gates = c(60000,NA))
    })
    
    par(mfrow=c(2,3))
    a <- lapply(1:length(live.flowD), function(x){
      plotDens(fs[[x]],channels=c("FSC-A","SSC-A"))
      lines(live.flowD[[x]]@filter)
    })
    
    live.fs <- flowSet(lapply(live.flowD, function(fD){fD@flow.frame}))
    
  #--------------------------------Gating Singlets-------------------------------------------------
    
    singlets.flowD <- fsApply(live.fs,function(f){
      flowD <- GateSinglets(f,scat.chans,channels.ind)
      return(flowD)
    })
    
    par(mfrow=c(2,3))
    a <- lapply(1:length(singlets.flowD), function(x){
      plotDens(live.fs[[x]],channels=c("FSC-A","FSC-H"))
      lines(singlets.flowD[[x]]@filter)
    })
    
    singlets.fs <- flowSet(lapply(singlets.flowD, function(fD){fD@flow.frame}))
    
  #-------------------------------tsne normal only------------------------------------------------
    if(x == 1){
      
    
    
    set.seed(1234)
    channels.ind <- setdiff(1:length(colnames(fs[[length(singlets.fs)]])), c(scat.chans,time.loc))
    # fb_result <- singlets.fs[[length(singlets.fs)]]@exprs[sample(1:nrow(singlets.fs[[19]]@exprs), 25000),channels.ind]
    fs_sub <- singlets.fs
    fs_sub <- sapply(1:length(singlets.fs), function(y){
      if(y==length(singlets.fs)){
        return(NULL)
        }
      return(singlets.fs[[y]])
    })
    fs_sub[sapply(fs_sub, is.null)] <- NULL
    
    fs_sub <- flowSet(fs_sub)
    
    # invLogicle <- inverseLogicleTransform(lgl)
    # fs_sub_inv <- fsApply(fs_sub,function(f) transform(f, invLogicle))
    
    fb_result <- fb_result_function(fs_sub, n=2000)
    
    
    tsne_out <- Rtsne(unique(fb_result[,channels.ind]))
    print("end of tSNE")
    save(tsne_out, file=paste0("~/results/tSNE_Rdata/normal_tube",tube_normal[t], "_Rtsne.RData"))
    
    # load(file=paste0("~/results/tSNE_Rdata/Patient",x,"_tube00",i,"_Rtsne.RData"))
    
    data_plot <- tsne_out$Y
    
    
    #--------------------------plotting---------------------------------------------------------
    markers_to_see <- c(scat.chans,markers_list)
    # markers_to_see <- markers_list
    png(file=paste0("~/results/Plots/tSNE_normal/normals_tube",tube_normal[t],".png"),
        width = 2*300, height= ceiling(length(markers_to_see)/2) *300)
    par(mfrow=c(ceiling(length(markers_to_see)/2),2), mar=c(5,5,5,6))
    
    
    
    for(z in 1:length(markers_to_see)){
      flowresults_marker <- fb_result[,z]
      tsnecol <- c(NA)
      flowresults_scale <- c(NA)
      
      if(length(which(flowresults_marker < 0)) > 0){
        flowresults_marker[which(flowresults_marker < 0)] <- 0
      }
      
      for (j in 1:length(flowresults_marker)) {
        flowresults_scale[j] <- flowresults_marker[j] - min(flowresults_marker)
      }
      max <- max(flowresults_scale)
      
      for (j in 1:length(flowresults_scale)) {
        flowresults_scale[j] <- flowresults_scale[j]/max
      }
      
      # colorScale <- NULL
      # colorScale[which(flowresults_marker_temp < 0.2)] <- 0
      # colorScale[intersect(which(flowresults_marker_temp>=0.2),which(flowresults_marker_temp <0.4))] <- 1
      # colorScale[intersect(which(flowresults_marker_temp>=0.4),which(flowresults_marker_temp <0.6))] <- 2
      # colorScale[intersect(which(flowresults_marker_temp>=0.6),which(flowresults_marker_temp <0.8))] <- 3
      # colorScale[which(flowresults_marker_temp >= 0.8)] <- 4
      # 
      NumEvents <- nrow(fb_result)
      ColourPal <- rainbow(round(NumEvents*3/2))
      
      
      markernames <- colnames(fs[[1]]@exprs)
      tsnecol <- color.scale(flowresults_scale,
                             extremes=c("blue1", "green", "yellow", "orange", "red"), color.spec="rgb")
      tsnecol2 <- color.scale(1:length(flowresults_scale), extremes=c("blue1", "green", "yellow", "orange", "red"), color.spec="rgb")
      
      plot(data_plot, col=tsnecol, main=paste0(markernames[markers_to_see[z]], " ", names(markers_to_see)[z]),
           xlab="tsne1",ylab="tsne2", pch=19, cex=0.5, cex.lab=2, cex.main=2,
           cex.axis=2)
      
      image.plot( legend.only=TRUE, zlim= range(fb_result[,z]), col=tsnecol2[seq(1,NumEvents,by=1)],
                  legend.width=2, legend.shrink=1.0, legend.mar=8, nlevel=400, axis.args = list(cex.axis = 2))
      
      
      
    }
    dev.off()
    }   
    #--------------------------------tsne normals overlay with patient-----------------------------------
    
    # set.seed(1234)
    
    # fts <- flowSet(singlets.fs[[1]],singlets.fs[[3]],singlets.fs[[5]],singlets.fs[[7]],
    #                singlets.fs[[9]],singlets.fs[[11]],singlets.fs[[13]],singlets.fs[[15]],
    #                singlets.fs[[17]],singlets.fs[[19]])
    
    # invLogicle <- inverseLogicleTransform(lgl)
    # singlets.fs_inv <- fsApply(singlets.fs,function(f) transform(f, invLogicle))
    # 
    fb_result <- fb_result_function(singlets.fs, n=4000)
    
    # fb_result_fs <- flowSet(flowFrame(fb_result[c(1:6000),]),flowFrame(fb_result[c(6001:12000),]),
    #                         flowFrame(fb_result[c(12001:18000),]),flowFrame(fb_result[c(18001:24000),]))
    
    lv <- NULL
    for(b in 1:length(singlets.fs) ){
      lv <- c(lv,rep(b,4000))
    }
    

    
    fb_result_f <- data.frame(cbind(fb_result,lv))
    
    set.seed(1234)
    tsne_out <- Rtsne(fb_result[,channels.ind],check_duplicates = FALSE)
    print("end of tSNE")
    save(tsne_out, file=paste0("~/results/tSNE_Rdata/normals+patient_",patientslb2ind_1[x,"PatientID"],"_tube",tube_normal[t],"_Rtsne.RData"))
    
    # load(file=paste0("~/results/tSNE_Rdata/Patient",x,"_tube00",i,"_Rtsne.RData"))
    
    data_plot <- tsne_out$Y
    
    markers_to_see <- c(scat.chans,markers_list)
    png(file=paste0("~/results/Plots/tSNE_normal/AllMarkers_normals+patient_",patientslb2ind_1[x,"PatientID"],"_tube",tube_normal[t],".png"),
        width = 2*300, height= ceiling(length(markers_to_see)/2) *300)
    par(mfrow=c(ceiling(length(markers_to_see)/2),2), mar=c(5,5,5,6))
    
    
    
    for(z in 1:length(markers_to_see)){
      flowresults_marker <- fb_result[,z]
      tsnecol <- c(NA)
      flowresults_scale <- c(NA)
      
      if(length(which(flowresults_marker < 0)) > 0){
        flowresults_marker[which(flowresults_marker < 0)] <- 0
      }
      
      for (j in 1:length(flowresults_marker)) {
        flowresults_scale[j] <- flowresults_marker[j] - min(flowresults_marker)
      }
      max <- max(flowresults_scale)
      
      for (j in 1:length(flowresults_scale)) {
        flowresults_scale[j] <- flowresults_scale[j]/max
      }
      
      # colorScale <- NULL
      # colorScale[which(flowresults_marker_temp < 0.2)] <- 0
      # colorScale[intersect(which(flowresults_marker_temp>=0.2),which(flowresults_marker_temp <0.4))] <- 1
      # colorScale[intersect(which(flowresults_marker_temp>=0.4),which(flowresults_marker_temp <0.6))] <- 2
      # colorScale[intersect(which(flowresults_marker_temp>=0.6),which(flowresults_marker_temp <0.8))] <- 3
      # colorScale[which(flowresults_marker_temp >= 0.8)] <- 4
      # 
      NumEvents <- nrow(fb_result)
      ColourPal <- rainbow(round(NumEvents*3/2))
      
      
      markernames <- colnames(fs[[1]]@exprs)
      tsnecol <- color.scale(flowresults_scale,
                             extremes=c("blue1", "green", "yellow", "orange", "red"), color.spec="rgb")
      tsnecol2 <- color.scale(1:length(flowresults_scale), extremes=c("blue1", "green", "yellow", "orange", "red"), color.spec="rgb")
      
      plot(data_plot, col=tsnecol, main=paste0(markernames[markers_to_see[z]], " ", names(markers_to_see)[z]),
           xlab="tsne1",ylab="tsne2", pch=19, cex=0.5, cex.lab=2, cex.main=2,
           cex.axis=2)
      
      image.plot( legend.only=TRUE, zlim= range(fb_result[,z]), col=tsnecol2[seq(1,NumEvents,by=1)],
                  legend.width=2, legend.shrink=1.0, legend.mar=8, nlevel=400, axis.args = list(cex.axis = 2))
      
    
      
     
    }
  
    for(p in 1:length(lv)){
      if(lv[p] == unique(lv)[length(unique(lv))]){
        # lv[p] <- as.character(patientslb2ind_1[t,"PatientID"])
        lv[p] <- "#CC0000" #patient 
      }else{
        lv[p] <- "#A0A0A0" #normal
      }
    }

    plot(data_plot, col=lv, main=paste0("normals+patient_",patientslb2ind_1[x,"PatientID"]),
         xlab="tsne1", ylab="tsne2",pch=19, cex=0.5,cex.lab=2, cex.main=2, cex.axis=2)
    
    legend(x="bottomright",legend = c("Normals", "Patient"), xpd=TRUE,horiz = TRUE, inset = c(0,0),
          col=c("#A0A0A0","#CC0000"), cex=0.8, lwd=1)
    
    dev.off()
    
    
}  
    
  }
  
# }


