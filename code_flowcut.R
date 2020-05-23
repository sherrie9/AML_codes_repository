flowCutdir <- "~/results/Plots/flowCut/"

source("~/code/helperFunc_S.R")
source("~/code/helperFunc.R")
source("~/code/IMPC_Functions.R")

suppressWarnings(dir.create("~/results/Summary/cleaned_not_trans"))

if(interactive()){
  center <- readCenterFunc()
}

store.allFCS <- preProcessFCS_center(inputPath,center)
store.allFCS.normal <- preProcessFCS_normal(inputPathNormal)
Patients <- unique(store.allFCS$PatientID)

for(x in 1:length(Patients)){
  

  FCS_1 <- store.allFCS[which(store.allFCS[,'PatientID'] == Patients[x]), ]
  tube <- unique(FCS_1[,'Tube'])
  

  for(i in 1:length(tube)){
    
    
    
    FCS_11 <- FCS_1[which(FCS_1[,'Tube'] == tube[i]),]
    
    # fs.raw <- makefs(FCS_11)
    # fs.raw <- read.flowSet(paste0(FCS_11$Path1,"/",FCS_11$Path2,"/",FCS_11$Path3,"/",FCS_11$FCS.files))
    
    f_list <- list()
    path <- list()
    f_org <- list()
    ind.marg.neg <- list()
    for(m in 1:nrow(FCS_11)){
      path[[m]] <- FCS_11[m,]
      f <- read.FCS(filename = paste0(path[[m]]$Path1,"/",
                                      path[[m]]$Path2,"/",
                                      path[[m]]$Path3,"/",path[[m]]$FCS.files))
      f_org[[m]] <- f
      scat.chans <- c(grep(colnames(f),pattern = "FSC*"), grep(colnames(f),pattern = "SSC*"))
      names(scat.chans) <- colnames(f)[scat.chans]
      time.loc <- which(tolower(colnames(f)) == "time"); names(time.loc) <- NULL
      
      res.remove <- removeMargNegs(f, chans=scat.chans, neg=0, verbose= F)
      f <- res.remove$f; ind.marg.neg[[m]] <- res.remove$ind.marg.neg; remove(res.remove)
      
      
  
      f <-compensate(f,f@description$SPILL)
      f_list[[m]] <- f
    }
     f_list <- makefs(f_list)
    
     fs <- flowSet(f_list)
      
     gf <- globalFrame(fs)
     lgl <- gf$lgl
    
     fs <- fsApply(fs, function(f) transform(f, lgl))
      
  
    
      for(m in 1:length(fs)){
        
      f <- fs[[m]]
      fileID <- gsub(".fcs","",path[[m]]$FCS.files)
      channels.to.clean <- which(complete.cases(f@parameters@data$desc) == TRUE)
    
      png(file=paste0("~/results/Plots/flowCut/",fileID,".png"),
          width = 5*300, height=2*300)
      
      f.Clean <- flowCut(f, Channels = channels.to.clean,Directory = flowCutdir,
                         FileID = fileID ,Plot = 'All', PrintToConsole = TRUE)
      
      dev.off()
      
      ind.clean <- f.Clean$ind
      ind.data <- f.Clean$data
      if(length(ind.clean) > 0 ){
        #f_cleaned_not_trans is the file obtained after flowcut and margin removal but not transformation
        f_cleaned_not_trans <- f_org[[m]]
        f_cleaned_not_trans <- f_cleaned_not_trans[-ind.marg.neg[[m]]][-ind.clean]
       
        suppressWarnings(write.FCS(f_cleaned_not_trans, file=paste0("~/results/Summary/cleaned_not_trans/", 
                                                                    path[[m]]$FCS.files)))
      }
      
      
      
      if(!file.exists("~/results/Summary/flowCut_Summary.csv")){
        header <- c("File_Name", "Flag","Events_Removed(%)")
        write.csv((t(header)), row.names = FALSE, file="~/results/Summary/flowCut_Summary.csv")
      }

      write.table(t(c(fileID,f.Clean$data[17,1], f.Clean$data[13,1])),
                  row.names = FALSE, file="~/results/Summary/flowCut_Summary.csv", sep=",",
                  col.names = FALSE, append = TRUE)
      
    
      }
    }
   
  }

