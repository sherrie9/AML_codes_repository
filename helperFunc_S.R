readDIVAFunc <- function(){
  diva <- readline(prompt = "Enter DIVA number (6/8)")
  diva <- as.integer(diva)
  
  if(diva == 6 | diva == 8){
    return(diva)
  }else{
    print("Incorrect Diva number")
    diva <- readDIVAFunc()
    return(diva)
  }
  
}


readCenterFunc <- function(){
  center <- readline(prompt="Enter lab name: ")
  center <- tolower(center)
  if(center == "israel"){
    center <- "Israel"
  }
  
  return(center)
}


preProcessFCS_normal <- function(inputPathNormal){
  library("flowCore")
  library("stringr")
  
  store.allFCS <- NULL
  
  
  for(i in 1:length(inputPathNormal)){
    
    pathFCS <- unlist(inputPathNormal[i])
    
    allFCS <- dir(pathFCS, full.names = T, recursive = T, pattern = "*.fcs")
    
    if(length(allFCS)==0){
      next
    }
    store.allFCS.temp <- sapply(1:length(allFCS), function(x){pathFCS})
    
    store.allFCS.temp <- cbind(store.allFCS.temp, sapply(1:length(allFCS), function(x){unlist(strsplit(allFCS[x], split = "/"))[length(unlist(strsplit(allFCS[x], split = "/")))-1]}))
    
    # store.allFCS.temp <- cbind(store.allFCS.temp, sapply(1:length(allFCS), function(x){unlist(strsplit(allFCS[x], split = "/"))[length(unlist(strsplit(allFCS[x], split = "/")))-1]}))
    
    
    
    # store.allFCS.temp <- apply(store.allFCS.temp,1, function(x){paste0(x[1],x[2],x[3])})
    # if(all(grepl("Normal", store.allFCS.temp))){
    #   store.allFCS.temp <- cbind(store.allFCS.temp, "Normal")
    # }

    store.allFCS.temp <- cbind(store.allFCS.temp, sapply(1:length(allFCS), function(x){unlist(strsplit(allFCS[x], split = "/"))[length(unlist(strsplit(allFCS[x], split = "/")))]}))
    store.allFCS.temp <- cbind(store.allFCS.temp,gsub(".fcs","",substring(unname(sapply(store.allFCS.temp[,3], function(x){
      unlist(strsplit(x,split="Tube"))[length(unlist(strsplit(x,split="Tube")))]})),2)))
   
    if(all(grepl("Normal", store.allFCS.temp[,2]))){
      store.allFCS.temp <- cbind(store.allFCS.temp, "Normal")
      store.allFCS.temp <- cbind(store.allFCS.temp,unname(sapply(store.allFCS.temp[,2],
                                                                 function(x){unlist(strsplit(x,split=" "))[length(unlist(strsplit(x,split = " ")))]})))
    }
  
    
    # if(all(grepl("DIVA 6", store.allFCS.temp[,1], ignore.case = TRUE ))){
    #   store.allFCS.temp <- cbind(store.allFCS.temp,"DIVA 6")
    # }
    # if(all(grepl("DIVA 8", store.allFCS.temp[,1],  ignore.case = TRUE))){
    #   store.allFCS.temp <- cbind(store.allFCS.temp,"DIVA 8")
    # }
    
    colnames(store.allFCS.temp) <- c("Path1","Path2","FCS.files","Tube","Status","PatientID")
    store.allFCS.temp <- as.data.frame(store.allFCS.temp)
    
    store.allFCS <- rbind(store.allFCS, store.allFCS.temp)
  }#end of for loop
  
  # a lot of meta data are messed up, need to tidy them.
  # store.allFCS$Status[which(tolower(store.allFCS$Status) == "dag22")] <- "Day22"
  # store.allFCS$Status[which(store.allFCS$Status=="day22")] <- "Day22"
  # store.allFCS$Status[which(store.allFCS$Status=="before cons")] <- "before consolidation"
  # store.allFCS$Status[which(store.allFCS$Status=="last before cons")] <- "last before consolidation"
  # store.allFCS$Status[which(store.allFCS$Status=="last before 2nd ind")] <- "last before 2nd induction"
  # store.allFCS$Status[which(store.allFCS$Status=="diagnos")] <- "Diagnos"
  # store.allFCS$Status[which(store.allFCS$Status=="Day29")] <- "Day30"
  # store.allFCS$Status[which(store.allFCS$Status=="Day39")] <- "Day37"
  store.allFCS$Tube <- unlist(lapply(strsplit(as.character(store.allFCS$Tube),split="_"), function(l){
    if(length(l) > 1){
      return(l[1])
    }else{
      return(l)
    }
  })) 
  
  
  return(store.allFCS)
}

preProcessFCS_center <- function(inputPath, center){
  library("flowCore")
  library("stringr")
  
  store.allFCS <- NULL

    # pathFCS <- unlist(inputPath[i])
    
    allFCS <- dir(inputPath, full.names = T, recursive = T, pattern = "*.fcs")
    
    store.allFCS.temp <- sapply(1:length(allFCS), function(x){inputPath})
    
    store.allFCS.temp <- cbind(store.allFCS.temp, sapply(1:length(allFCS), function(x){unlist(strsplit(allFCS[x], split = "/"))[length(unlist(strsplit(allFCS[x], split = "/")))-2]}))
    
    store.allFCS.temp <- cbind(store.allFCS.temp, sapply(1:length(allFCS), function(x){unlist(strsplit(allFCS[x], split = "/"))[length(unlist(strsplit(allFCS[x], split = "/")))-1]}))
    
    
    
  
    
    store.allFCS.temp <- cbind(store.allFCS.temp, sapply(1:length(allFCS), function(x){unlist(strsplit(allFCS[x], split = "/"))[length(unlist(strsplit(allFCS[x], split = "/")))]}))
    
    if(center == "***"|center =="***"|center=="***"|center=="***"
       |center=="*****"){
      store.allFCS.temp <- cbind(store.allFCS.temp,gsub(".fcs","",substring(unname(sapply(store.allFCS.temp[,4], function(x){
        unlist(strsplit(tolower(x),split="tube"))[length(unlist(strsplit(tolower(x),split="tube")))]})),2)))
      
      
      
      store.allFCS.temp <- cbind(store.allFCS.temp, unname(sapply(store.allFCS.temp[,3],
                                                                  function(x){unlist(strsplit(x,split="_"))[length(unlist(strsplit(x,split="_")))]})))
      
      store.allFCS.temp <- cbind(store.allFCS.temp,store.allFCS.temp[,2])
      
      
      colnames(store.allFCS.temp) <- c("Path1","Path2","Path3","FCS.files","Tube","Status","PatientID")
      store.allFCS.temp <- as.data.frame(store.allFCS.temp)
      
      store.allFCS <- rbind(store.allFCS, store.allFCS.temp)
      
      
      if(center == "gothenburg"){
        
      
      # a lot of meta data are messed up, need to tidy them.
      store.allFCS$Status[which(tolower(store.allFCS$Status) == "dag22")] <- "Day22"
      store.allFCS$Status[which(store.allFCS$Status=="day22")] <- "Day22"
      store.allFCS$Status[which(store.allFCS$Status=="before cons")] <- "before consolidation"
      store.allFCS$Status[which(store.allFCS$Status=="LastBefore2ndInd")] <- "last before 2nd induction"
      store.allFCS$Status[which(store.allFCS$Status=="last before cons")] <- "last before consolidation"
      store.allFCS$Status[which(store.allFCS$Status=="Last before cons")] <- "last before consolidation"
      store.allFCS$Status[which(store.allFCS$Status=="before consolidation")] <- "last before consolidation"
      store.allFCS$Status[which(store.allFCS$Status=="last before 2nd ind")] <- "last before 2nd induction"
      store.allFCS$Status[which(store.allFCS$Status=="last before 2nd course")] <- "last before 2nd induction"
      store.allFCS$Status[which(store.allFCS$Status=="Lats")] <- "last before 2nd induction"
      store.allFCS$Status[which(store.allFCS$Status=="diagnos")] <- "Diagnos"
      store.allFCS$Status[which(store.allFCS$Status==" Day0")] <- "Day0"
      store.allFCS$Status[which(store.allFCS$Status=="d22 post 2nd ind")] <- "d22post2ndind"
      store.allFCS$Status[which(store.allFCS$Status=="d22 after 2nd ind")] <- "d22post2ndind"
      # store.allFCS$Status[which(store.allFCS$Status=="Day29")] <- "Day30"
      # store.allFCS$Status[which(store.allFCS$Status=="Day39")] <- "Day37"
      store.allFCS$Tube <- unlist(lapply(strsplit(as.character(store.allFCS$Tube),split="_"), function(l){
        if(length(l) > 1){
          return(l[1])
        }else{
          return(l)
        }
      }))
      }
      if(center == "umea"){
        store.allFCS$Status[which(store.allFCS$Status=="diagnosis")] <- "Diagnosis"
        store.allFCS$Status[which(store.allFCS$Status=="day22")] <- "Day22"
        store.allFCS$Status[which(store.allFCS$Status=="last before ind2")] <- "last before 2nd ind"
        store.allFCS$Status[which(store.allFCS$Status=="Last before 2nd induction")] <- "last before 2nd ind"
        store.allFCS$Status[which(store.allFCS$Status=="last before 2 ind")] <- "last before 2nd ind"
        # store.allFCS <- as.character(store.allFCS)
        store.allFCS$Status <- as.character(store.allFCS$Status)
        store.allFCS$Status[which(store.allFCS$Status=="day29")] <- "d22+1w"
        store.allFCS$Status[which(store.allFCS$Status=="day22+1w")] <- "d22+1w"
        store.allFCS$Status[which(store.allFCS$Status=="last before cons")] <- "last before consolidation"
        store.allFCS <- as.data.frame(store.allFCS)
        }
      
      if(center =="Israel"){
        store.allFCS$Status[which(store.allFCS$Status=="Last bef 2nd induction")] <- "Last bef 2nd ind"
      }
        
      if(center =="karolinska"){
        store.allFCS$Status <- as.character(store.allFCS$Status)
        store.allFCS$Status[which(store.allFCS$Status=="diagnosis")] <- "Diagnosis"
        store.allFCS$Status[which(store.allFCS$Status=="day29")] <- "d22+1w"
        store.allFCS$Status[which(store.allFCS$Status=="day22")] <- "Day22"
        store.allFCS$Status[which(store.allFCS$Status=="before2ind dx")] <- "last before 2nd ind"
        store.allFCS$Status[which(store.allFCS$Status=="last before ind2")] <- "last before 2nd ind"
        store.allFCS$Status[which(store.allFCS$Status=="before 2 ind")] <- "last before 2nd ind"
        
        store.allFCS$Status[which(store.allFCS$Status=="before1cons")] <- "last before consolidation"
        store.allFCS$Status[which(store.allFCS$Status=="before1 con")] <- "last before consolidation"
        
        store.allFCS$Status[which(store.allFCS$Status=="last before cons")] <- "last before consolidation"
        store.allFCS$Status[which(store.allFCS$Status=="first after 2nd ind dx")] <- "last before consolidation"
 
        }
      
    }
    
    if(center == "copenhagen"){
      store.allFCS.temp <- cbind(store.allFCS.temp,
                                 unname(sapply(store.allFCS.temp[,3], function(x){
                                   unlist(strsplit(x,split="_"))}[length(unlist(strsplit(x,split = "_")))])))

      store.allFCS.temp <- cbind(store.allFCS.temp, store.allFCS.temp[,2])
      colnames(store.allFCS.temp) <- c("Path1","Path2","Path3","FCS.files","Status","Tube","PatientID")
      store.allFCS.temp <- as.data.frame(store.allFCS.temp)
      
      
      # store.allFCS.temp$Status[which(store.allFCS.temp$Status == "1HAM")] <- "Kons 1HAM"
      # store.allFCS.temp$Tube <- sapply(store.allFCS.temp$Tube,function(x)sub("^\\s+", "", x))
      
      tubenames <- sapply(1:nrow(store.allFCS.temp), function(x){
        f <- read.FCS(filename = paste0(store.allFCS.temp[x,]$Path1,"/",store.allFCS.temp[x,]$Path2,"/",
                                        store.allFCS.temp[x,]$Path3,"/",store.allFCS.temp[x,]$FCS.files))
        # tubename <- paste0(na.omit(unname(f@parameters@data$desc)),collapse = "-")
        tubename <- na.omit(unname(f@parameters@data$desc))
        a <- ifelse(length(which(tubename == "HLA-DR V450")), yes=tubename[which(tubename =="HLA-DR V450")] <- "HLA-DR", no="")
        a <- ifelse(length(which(tubename == "CD117 PC7")), yes=tubename[which(tubename =="CD117 PC7")] <- "CD117", no="")
        a <- ifelse(length(which(tubename == "CD45 V500")), yes=tubename[which(tubename =="CD45 V500")] <- "CD45", no="")
        a <- ifelse(length(which(tolower(tubename) == "cd14 percp")), yes=tubename[which(tolower(tubename) =="cd14 percp")] <- "CD14", no="")
        a <- ifelse(length(which(tubename == "Propidiumiodid")), yes=tubename[which(tubename =="Propidiumiodid")] <- "propidiiumiodid", no="")
        a <- ifelse(length(which(tubename == "IgG1 V500")), yes=tubename[which(tubename =="IgG1 V500")] <- "mouse IgG1 V500", no="")
        a <- ifelse(length(which(tubename == "cyTdT")), yes=tubename[which(tubename =="cyTdT")] <- "CyTdT", no="")
        a <- ifelse(length(which(tubename == "cyCD3")), yes=tubename[which(tubename =="cyCD3")] <- "CyCD3", no="")
        
        
        # a <- ifelse(length(which(tubename == "CD45 H7")), yes=tubename[which(tubename =="CD45 H7")] <- "CD45", no="")
        tubename <- paste0(tubename,collapse = "-")
        return(tubename)
      })
      
      store.allFCS.temp$Tube <- tubenames
      
      store.allFCS <- store.allFCS.temp
      
      store.allFCS$Status[which(store.allFCS$Status == "before 2nd ind")] <- "before Ind2"
      store.allFCS$Status[which(store.allFCS$Status =="1HAM")] <-"before cons1"
    }
  
  
  return(store.allFCS)
}




preprocessFCS_files <- function(paths,scatter_chans){
  
 
  
  f_list <- lapply(paths,function(p){
    f <- read.FCS(p)
    return(f)})
  
  markers <- c("CD56","CD13","CD34|CD34*","CD117|CD117 BIO|CD117BC","CD33",
               "CD11b|CD11B","HLADR|HLA-DR","CD45|CD45V")
  colnames <- c("FSC-A*","SSC-A","FITC-A","PE-A","PerCP-Cy5-5-A|PerCP-A","PE-Cy7-A","APC-A","APC-H7-A|APC-Cy7-A",
                "V450-A|Horizon V450-A|HZ V450-A|Pacific Blue-A","Horizon V500-A|HZ V500-A|AmCyan-A|V500-A")
  
  f_list <- lapply(f_list, function(f){
    
    f@parameters@data$desc[grep(f@parameters@data$desc, pattern = "CD34|CD34 PC5")] <- "CD34"
    f@parameters@data$desc[grep(f@parameters@data$desc, pattern = "CD117|CD117 BIO|CD117BC")] <- "CD117"
    f@parameters@data$desc[grep(f@parameters@data$desc, pattern = "CD11b|CD11B")] <- "CD11B"
    f@parameters@data$desc[grep(f@parameters@data$desc, pattern = "HLADR|HLA-DR")] <- "HLADR"
    f@parameters@data$desc[grep(f@parameters@data$desc, pattern = "CD45|CD45V")] <- "CD45"
    
    colnames(f@exprs)[grep(colnames(f@exprs), pattern = "PerCP-Cy5-5-A|PerCP-A")] <- "PerCP-A"
    colnames(f@exprs)[grep(colnames(f@exprs), pattern = "APC-H7-A|APC-Cy7-A")] <- "APC-H7/Cy7-A"
    colnames(f@exprs)[grep(colnames(f@exprs), pattern = "V450-A|Horizon V450-A|HZ V450-A|Pacific Blue-A")] <- "V450-A"
    colnames(f@exprs)[grep(colnames(f@exprs), pattern = "Horizon V500-A|HZ V500-A|AmCyan-A|V500-A")] <- "V500-A"
    
    colnames(f)[grep(colnames(f), pattern = "PerCP-Cy5-5-A|PerCP-A")] <- "PerCP-A"
    colnames(f)[grep(colnames(f), pattern = "APC-H7-A|APC-Cy7-A")] <- "APC-H7/Cy7-A"
    colnames(f)[grep(colnames(f), pattern = "V450-A|Horizon V450-A|HZ V450-A|Pacific Blue-A")] <- "V450-A"
    colnames(f)[grep(colnames(f), pattern = "Horizon V500-A|HZ V500-A|AmCyan-A|V500-A")] <- "V500-A"
    
    colnames(f@description$SPILL)[grep(colnames(f@description$SPILL), pattern = "PerCP-Cy5-5-A|PerCP-A")] <- "PerCP-A"
    colnames(f@description$SPILL)[grep(colnames(f@description$SPILL), pattern = "APC-H7-A|APC-Cy7-A")] <- "APC-H7/Cy7-A"
    colnames(f@description$SPILL)[grep(colnames(f@description$SPILL), pattern = "V450-A|Horizon V450-A|HZ V450-A|Pacific Blue-A")] <- "V450-A"
    colnames(f@description$SPILL)[grep(colnames(f@description$SPILL), pattern = "Horizon V500-A|HZ V500-A|AmCyan-A|V500-A")] <- "V500-A"
    # return(f@parameters@data$desc[-which(is.na(f@parameters@data$desc))])
    return(f)
    })
  

  
  if(length(scatter_chans) == 2){
    f_list <- lapply(f_list,function(f){
      if(length(grep(colnames(f),pattern = "FSC-H")) * length(grep(colnames(f),pattern = "SSC-H")) > 0 ){
        f@exprs <- f@exprs[,-c(which(colnames(f@exprs)=="FSC-H"),which(colnames(f@exprs)=="SSC-H"))]
        f@parameters@data$name <- f@parameters@data$name[-c(grep(f@parameters@data$name,pattern="FSC-H"),
                                                            grep(f@parameters@data$name,pattern="SSC-H"))]
      }else if(length(grep(colnames(f),pattern = "FSC-H")) >0 ){
        f@exprs <- f@exprs[,-which(colnames(f@exprs)=="FSC-H")]
        f@parameters@data$name <- f@parameters@data$name[-grep(f@parameters@data$name,pattern="FSC-H")]
      }
      
      return(f)
    })
  }
  if(length(scatter_chans) == 3){#if user wants to use FSC-A, FSC-H SSC-A channels, then we need to remove files that don't have FSC-H
                                 #and remove channel SSC-H
    f_list <- lapply(f_list,function(f){
      if(length(grep(colnames(f),pattern = "FSC-H")) == 0){
        f <- NULL
      }else{
        if(length(grep(colnames(f),pattern = "SSC-H")) > 0 ){
          f@exprs <- f@exprs[,-c(which(colnames(f@exprs)=="SSC-H"))]
          f@parameters@data <- f@parameters@data[-grep(f@parameters@data$name,pattern="SSC-H"),]
          if(f@parameters@data$name[1] == "Time"){
            f@parameters@data <- rbind(f@parameters@data[-1,],f@parameters@data[1,])
            f@exprs <- cbind(f@exprs[,-1],f@exprs[,1])
            colnames(f@exprs) <- colnames(f)
          }
          # f@parameters@data$name <- f@parameters@data$name[-grep(f@parameters@data$name,pattern="SSC-H")]
        }
        
      }
      return(f)
    })
  }
  
  f_list[sapply(f_list, is.null)] <- NULL
  
  f_list <- lapply(f_list,function(f){
   if(colnames(f)[1] == "Time"){

     f@parameters@data <- rbind(f@parameters@data[-1,],f@parameters@data[1,])
     f@exprs <- cbind(f@exprs[,-1],f@exprs[,1])
     colnames(f@exprs) <- colnames(f)
     # colnames(f) <- f@parameters@data$name
   }

   return(f)
  })
  
  
  
  return(f_list)
  # f_colnames <- lapply(f_list,function(f){return(unname(colnames(f@exprs)))})
  
  
}


















makefs <- function(f_list){
  
   
   f_list<- sapply(f_list,function(f){
    if(length(grep("horizon v450-a",tolower(colnames(f_list[[1]])))) > 0){
      colnames(f)[grep("horizon v450-a",tolower(colnames(f_list[[1]])))] <- "V450-A"
    }
     
     if(length(grep("horizon v500-a",tolower(colnames(f_list[[1]])))) > 0){
       colnames(f)[grep("horizon v500-a",tolower(colnames(f_list[[1]])))] <- "V500-A"
     }
     
     return(f)
   })
   
   # fs<-flowSet(f_list)
   return(f_list)
}



is.equal <- function(mylist) {
  
  check.eq <- sapply(mylist[-1], function(x) {x == mylist[[1]]})
  
  as.logical(apply(check.eq, 1, prod))                   
  
}


globalFrame <- function(fs){
      library('plyr')
      print("Start printing Global Frame")
  

      
       gFrame <- fsApply(fs, function(f){
        # markers <- c("CD56","CD36","CD15","CD7","CD99","CD13","CD64","NG2","CD96","CD11a","CD34","CD117",
        #              "CD33","CD2","CD123","CD133","CD11b","CD14","CD19","CD38","CD4","HLADR","CD45")
        # channels.ind <- Find.markers(f,markers)
        
        # markers <- as.character(f@parameters@data[1:ncol(f), 2])
        scatter.channels <- grep("FS|SS",f@parameters@data[,'name'])
        time.loc <- which(tolower(f@parameters@data$name) == "time"); names(time.loc) <- NULL
        # channels.ind <- Find.markers(f,na.omit(markers))
        channels.ind <- setdiff(1:length(colnames(f)), c(scatter.channels,time.loc))
        
        if(nrow(f@exprs) < 8000){
          temp <- f@exprs[,channels.ind]
        }else{
          temp <- f@exprs[sample(1:length(f@exprs[,1]), 8000), channels.ind]
        }
        
        
        return(data.frame(temp,check.names = F))
      })
       

       # tempFrame <- gFrame[[1]]
       # 
       # for(q in 1:length(gFrame)){
       #   if(!identitical(names(gFrame[[q]]), names(tempFrame))){
       #     temp <- gFrame[[q]]
       #     ids <- unlist(sapply(names(tempFrame), function(n){grep(n,names(temp))}))
       #     gFrame[[q]] <- temp[,ids]
       #     names(gFrame[[q]]) <- names(tempFrame)
       #   }
       # }
       
       gFrame <- do.call(rbind, gFrame)
       
       fileswNAs <- NULL
       colswNAs <- which(is.na(gFrame), arr.ind = TRUE)
       if(length(colswNAs) > 0){
         rowswNAs <- unique(colswNAs[,'row'])
         colswNAs <- unique(colswNAs[,'col'])
         fileswNAs <- unique(gFrame[rowswNAs,1])
         gFrame <- gFrame[-rowswNAs,]
       }
       
       g <- fs[[1]]
       
       Transform.idx <- unlist(sapply(colnames(gFrame), function(x) {grep(x, colnames(g))})) 
       gexprs.temp <- matrix(0, nrow = nrow(gFrame), ncol = ncol(g@exprs))
       gexprs.temp[, Transform.idx] <- as.matrix(gFrame[, 1:ncol(gFrame)])
       g@exprs <- gexprs.temp
       colnames(g@exprs) <- colnames(g) 
       print("End of creating the Global Frame")
       
       
       print("Start computing the transform using the estimateLogicle()")
       lgl <- estimateLogicle(g, channels = colnames(g)[Transform.idx])
       print("End of computing the transform using the estimateLogicle()")
       
       return(list(globalFrame=g, globalFrame.Matrix=gFrame, lgl=lgl))
}




GateSinglets <- function(f, scat.chans,channels.ind){
  
  if(length(which(is.na(f@exprs[,"FSC-A"]))) > 0){
    f@exprs <- f@exprs[-which(is.na(f@exprs[,"FSC-A"])),]
  }
  
  temp <- f@exprs[, c(scat.chans['FSC-A'], scat.chans['FSC-H'])]
 
  temp.flowD <- flowDensity(f, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA,F), gates = c(NA, 0.98*max(f@exprs[,c(scat.chans['FSC-H'])])))
  #x0 <- 125000
  x0 <- quantile(f@exprs[, c(scat.chans['FSC-H'])], c(0.05, 0.99))
  x0 <- x0[1] + 0.75*(x0[2]-x0[1])
 
  if(length(which(abs(temp[ , scat.chans['FSC-A']] - x0) < 5000)) < 2){
    if(length(which(abs(temp[ , scat.chans['FSC-A']] - x0) < 20000)) < 2){
      
      temp <- temp
    }else{
      temp <- temp[which(abs(temp[ , scat.chans['FSC-A']] - x0) < 20000),]
    }
    
    
    
  }else{
    temp <- temp[which(abs(temp[ , scat.chans['FSC-A']] - x0) < 5000),]
  }
  
  
  maxDens <- density(temp[, 2])
  peak.lcn <- findpeaks(maxDens$y) # get all peaks
  peak.lcn <- peak.lcn[which(peak.lcn[,1] > 0.1*max(peak.lcn[,1])),] # choose only peaks with significant y values
  if(length(peak.lcn) > 4){
    peak.lcn <- maxDens$x[peak.lcn[,2]] # get x-values
    peak.lcn <- peak.lcn[which.min(abs(peak.lcn - x0))] # assume we are interested in peak closed to diagonal
  }else{
    peak.lcn <- maxDens$x[peak.lcn[2]]
  }
  theta0 <- -atan(peak.lcn/x0)
  
  # rotate data
  temp.flowD <- flowDensity(f, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA,F), gates = c(NA, 0.975*max(f@exprs[,c(scat.chans['FSC-H'])])))
  rot <- rotate.data(f, c(scat.chans['FSC-A'], scat.chans['FSC-H']), theta = -theta0)$data
  
  maxDens <- density(rot@exprs[,c(scat.chans['FSC-H'])])
  maxDens <- smooth.spline(maxDens$x, maxDens$y, spar = 0.25)
  
  gate <- deGate(rot, channel = c(scat.chans['FSC-H']), all.cuts = T, tinypeak.removal = 0.001, percentile = 0)
  gate0 <- gate
  singlets.peak <- maxDens$x[which.max(maxDens$y)]
  
  if(gate!=-Inf){
    if(is.null(names(gate))){ 
      idx <- which(gate < (singlets.peak - 10000))
      if(length(idx > 0)){
        gate <- gate[idx]
        p05percentile <- min(deGate(rot, channel = c(scat.chans['FSC-H']), use.percentile = T, percentile = 0.0005), 50000) 
        gate <- gate[which.min(abs(singlets.peak - 0.45*(singlets.peak - p05percentile) - gate))]
      }else{
        gate <- deGate(rot, channel = c(scat.chans['FSC-H']), use.upper = T, upper = F, tinypeak.removal = 0.9, alpha = 0.01)
      }
    }
    
  }
  
  

  if((singlets.peak - gate) > 50000){
    gate <- deGate(rot, channel = c(scat.chans['FSC-H']), use.upper = T, upper = F, tinypeak.removal = 0.9, alpha = 0.5) - 15000
  }
  
  temp <- flowDensity(rot, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA, T), gates = c(NA, gate))
  #   , ellip.gate = T, scale = .9999999);
  temp@filter <- rotate.data(temp@filter, c(scat.chans['FSC-A'], scat.chans['FSC-H']), theta = theta0)$data
  temp@flow.frame <- rotate.data(getflowFrame(temp), c(scat.chans['FSC-A'], scat.chans['FSC-H']),theta = theta0)$data
  
  # Gate really closesly - need this to tell if SSC-H/SSC-W channels are interchanged 
  temp2 <- flowDensity(rot, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA, T), gates = c(NA, singlets.peak - 2000))
  temp2@flow.frame <- rotate.data(getflowFrame(temp2), c(scat.chans['FSC-A'], scat.chans['FSC-H']),theta = theta0)$data
  
  FSCsinglets.flowD <- temp
  
  return(FSCsinglets.flowD)
}


fb_result_function <- function(fs,n){
  
  fb <- c()
  
  # markers_list <- as.character(fs[[1]]@parameters@data[1:ncol(fs[[1]]), 2])
  # 
  # markers_list <- Find.markers(fs[[1]],na.omit(markers_list))
  
  scatter.channels <- grep("FS|SS",fs[[1]]@parameters@data[,'name'])
  time.loc <- which(tolower(fs[[1]]@parameters@data$name) == "time"); names(time.loc) <- NULL
  # channels.ind <- Find.markers(f,na.omit(markers))
  channels.ind <- setdiff(1:length(colnames(fs[[1]])), c(scatter.channels,time.loc))
  channels.indfcsssc<- c(scatter.channels,channels.ind)
  for (k in 1: length(fs)){
   
    if(nrow(fs[[k]]@exprs) < n){
      expr <- fs[[k]]@exprs
      expr <- apply(expr, 2, function(l){
        rn <- runif(n-length(l), min = mean(l)-2*sd(l), max = mean(l)+ 2 *sd(l))
        l <- c(l, rn)
      })
      
    }else{
      expr <- fs[[k]]@exprs[sample(c(1:nrow(fs[[k]]@exprs)), n), ]
    }
   
    fb <- rbind(fb,expr)
  }
  
  
  
  
  
  fb <- fb[,channels.indfcsssc]
  
 
  
  
  
  return(fb)
}


###################### RANDOM SAMPLING ##############################################################
randomResampling <- function(fs,n.samples, lv2){
  ns <- n.samples
  sub.idx <- which(lv2$st == "R") 
  sub.fs_R <- fs[sub.idx] 
  
  n.samples <- n.samples - length(sub.fs_R)
  sub.fs_R_S <- llply(1:n.samples, function(i){ 
    
    idx <- sample(c(1:length(fs)),1)
    f <- fs[[idx]]
    pop.size <- round(0.9 * nrow(f@exprs),0)
    f@exprs <- as.matrix(sample_n(as.data.frame(f@exprs),pop.size))
    fileName <- paste0("resample_",i,".fcs")
    f@description$FILENAME <- fileName
    return(f)
  }) 
  
  ptsID_S <- paste0("13",sample(c(3000:4000),n.samples),"000",sample(c(10:99),n.samples))
  lv2_R <- data.frame(cbind(ptsID_S,"R"))
  colnames(lv2_R) <- colnames(lv2) 
  lv2_R <- rbind(lv2[sub.idx,],lv2_R)
  singlets.fs_R <- c(lapply(1:length(sub.fs_R),function(k)return(sub.fs_R[[k]])),
                           sub.fs_R_S)
  
  
  ##-------------------------------------------------------------------------
  sub.idx <- which(lv2$st == "0") 
  sub.fs_0 <- fs[sub.idx] 
  n.samples <- ns
  n.samples <- n.samples - length(sub.fs_0)
  sub.fs_0_S <- llply(1:n.samples, function(i){ 
    
    idx <- sample(c(1:length(fs)),1)
    f <- fs[[idx]]
    pop.size <- round(0.9 * nrow(f@exprs),0)
    f@exprs <- as.matrix(sample_n(as.data.frame(f@exprs),pop.size))
    fileName <- paste0("resample_",i,".fcs")
    f@description$FILENAME <- fileName
    return(f)
  }) 
  
  ptsID_S <- paste0("13",sample(c(4000:5000),n.samples),"000",sample(c(10:99),n.samples))
  lv2_0 <- data.frame(cbind(ptsID_S,"0"))
  colnames(lv2_0) <- colnames(lv2) 
  lv2_0 <- rbind(lv2[sub.idx,],lv2_0)
  singlets.fs_0 <- c(lapply(1:length(sub.fs_0),function(k)return(sub.fs_0[[k]])),
                             sub.fs_0_S)
  singlets.fs <- flowSet(c(singlets.fs_R, singlets.fs_0))
  
  lv2 <- rbind(lv2_R,lv2_0)
  
  return(list(singlets.fs = singlets.fs, lv2 = lv2))


  }

### Random resample normal #################################################
randomResampling_Normal <- function(fs, n.samples){
  
  ##-------------------------------------------------------------------------
  # sub.idx <- which(lv2$st == "Normal") 
  # sub.fs_n <- fs[sub.idx] 
  # n.samples <- ns
  n.samples <- n.samples - length(fs)
  sub.fs_n_S <- llply(1:n.samples, function(i){ 
    
    idx <- sample(c(1:length(fs)),1)
    f <- fs[[idx]]
    pop.size <- round(0.9 * nrow(f@exprs),0)
    f@exprs <- as.matrix(sample_n(as.data.frame(f@exprs),pop.size))
    fileName <- paste0("resample_",i,".fcs")
    f@description$FILENAME <- fileName
    return(f)
  }) 
  
  # ptsID_S <- paste0("13",sample(c(5000:6000),n.samples),"000",sample(c(10:99),n.samples))
  # lv2_n <- data.frame(cbind(ptsID_S,"Normal"))
  # colnames(lv2_n) <- colnames(lv2) 
  # lv2_n <- rbind(lv2[sub.idx,],lv2_n)
  singlets.fs_n <- flowSet(c(lapply(1:length(fs),function(k)return(fs[[k]])),
                     sub.fs_n_S))
  
  return(singlets.fs_n)
}

