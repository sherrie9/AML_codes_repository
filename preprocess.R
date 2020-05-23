library(ggdendro)
library(scamp)
library(ggplot2)
library(cowplot)
library(knitr)
library(dplyr)
library(tidyr)


{library('flowCore')
  library('flowCut')
  library('stringr')
  library('Cairo')
  library('flowBin')
  library('flowDensity')
  library('pracma')
  library('fields')
  library('e1071')
  library('flowStats')
  library('gridExtra')
  library('plyr')
  library('doMC')}

source("~/code/helperFunc_S.R")
source("~/code/helperFunc.R")
source("~/code/IMPC_Functions.R")


#---------------Read in normal files----------------------------------
#---------------------------------------------------------------------
InPath1 <-"~/data/Normal BM DIVA 6/"
InPath2 <- "~/data/Normal BM DIVA 8/"
inputPath <- list(InPath1, InPath2)

store.FCS.normal <- preProcessFCS_normal(inputPath)
store.FCS.normal <- as.data.frame(store.FCS.normal)

store.FCS.normal <- cbind(store.FCS.normal,rep("Normal",nrow(store.FCS.normal)))
colnames(store.FCS.normal)[length(colnames(store.FCS.normal))] <- "AML_Status"

#-------------- Read in AML patient files----------------------------
#--------------------------------------------------------------------
centers <- c("copenhagen","gothenburg","Israel","umea","karolinska")

# inputPath <- list(InPath1, InPath2, InPath3, InPath4)
inpaths <- unname( sapply(centers, function(c)paste0("~/data/from ",c))) 

# store.allFCS <- preProcessFCS(inputPath)
store.allFCS <- NULL

for(i in 1:length(centers)){
  store.allFCS <- rbind(preProcessFCS_center(inpaths[i],centers[i]),store.allFCS)
}


# load(file="~/code/Rtsne_Results.RData")#patient 1 tube 1


flowCutdir <- "~/results/Plots/flowCut/"

Patients <- unique(store.allFCS$PatientID)

{store.allFCS$Tube <- as.character(store.allFCS$Tube)
  store.allFCS$Tube[c(which(store.allFCS$Tube=="1"),which(store.allFCS$Tube=="001"),
                      which(store.allFCS$Tube=="CD56-CD13-CD34-CD117-CD33-CD11b-HLA-DR-CD45"),
                      which(store.allFCS$Tube=="001_dx"))] <- "001"
  store.allFCS$Tube[c(which(store.allFCS$Tube=="2"),
                      which(store.allFCS$Tube=="CD36-CD64-CD34-CD117-CD33-CD14 H7-HLA-DR-CD45"),
                      which(store.allFCS$Tube=="002_dx"),which(store.allFCS$Tube=="002dx"))] <- "002"
  store.allFCS$Tube[c(which(store.allFCS$Tube=="3"),
                      which(store.allFCS$Tube=="CD15-NG2-CD34-CD117-CD2-CD19 H7-HLA-DR-CD45"),
                      which(store.allFCS$Tube=="003_dx"), which(store.allFCS$Tube=="003dx"),
                      which(store.allFCS$Tube=="003sin"))] <- "003"
  store.allFCS$Tube[c(which(store.allFCS$Tube=="4"),
                      which(store.allFCS$Tube=="CD7-CD96-CD34-CD117-CD123-CD38 H7-HLA-DR-CD45"),
                      which(store.allFCS$Tube=="004_dx"))] <- "004"
  store.allFCS$Tube[c(which(store.allFCS$Tube=="5"),
                      which(store.allFCS$Tube=="CD99-CD11a-CD34-CD117-CD133-CD38 H7-HLA-DR-CD45"),
                      which(store.allFCS$Tube=="005_dx"))] <- "005"
  store.allFCS$Tube <- as.factor(store.allFCS$Tube)
}
{store.allFCS$Status <- as.character(store.allFCS$Status)
  
  store.allFCS$Status[which(store.allFCS$Status == "Day0")] <- "Diagnosis"
  store.allFCS$Status[which(store.allFCS$Status == "day0")] <- "Diagnosis"
  store.allFCS$Status[which(store.allFCS$Status == "Diagnos")] <- "Diagnosis"
  store.allFCS$Status[which(store.allFCS$Status == "day22")] <- "Day22"
  store.allFCS$Status[which(store.allFCS$Status == "d22")] <- "Day22"
  store.allFCS$Status[which(store.allFCS$Status == "Day24")] <- "Day22"
  
  store.allFCS$Status[which(store.allFCS$Status == "Day30")] <- "d22+1w"
  store.allFCS$Status[which(store.allFCS$Status == "Day29")] <- "d22+1w"
  store.allFCS$Status[which(store.allFCS$Status == "day29")] <- "d22+1w"
  store.allFCS$Status[which(store.allFCS$Status == "day29 sin and dx")] <- "d22+1w"
  
  
  store.allFCS$Status[which(store.allFCS$Status == "day51")] <- "d22+3w"
  store.allFCS$Status[intersect(which(store.allFCS$PatientID=="13201400006"),
                                which(store.allFCS$Status == "day40"))] <- "d22+2w"
  
  store.allFCS$Status[which(store.allFCS$Status == "Day34")] <- "last before 2nd ind"
  store.allFCS$Status[which(store.allFCS$Status == "Day39")] <- "last before 2nd ind"
  store.allFCS$Status[which(store.allFCS$Status == "Day37")] <- "last before 2nd ind"
  store.allFCS$Status[which(store.allFCS$Status == "Day40")] <- "last before 2nd ind"
  
  
  store.allFCS$Status[which(store.allFCS$Status == "before Ind2")] <- "last before 2nd ind"
  store.allFCS$Status[which(store.allFCS$Status == "last before 2nd induction")] <- "last before 2nd ind"
  store.allFCS$Status[which(store.allFCS$Status == "Last bef 2nd ind")] <- "last before 2nd ind"
  store.allFCS$Status[which(store.allFCS$Status == "ind")] <- "last before 2nd ind"
  store.allFCS$Status[which(store.allFCS$Status == "before 2nd ind")] <- "last before 2nd ind"
  
  store.allFCS$Status[which(store.allFCS$Status == "last before consolidation")] <- "last before cons"
  store.allFCS$Status[which(store.allFCS$Status == "before cons1")] <- "last before cons"
  store.allFCS$Status[which(store.allFCS$Status == "cons")] <- "last before cons"
  store.allFCS$Status[which(store.allFCS$Status == "Last bef consolidation")] <- "last before cons"
  store.allFCS$Status[which(store.allFCS$Status == "day22post2ndind")] <- "d22post2ndind"
  store.allFCS$Status[which(store.allFCS$Status == "Last before consolidation")] <- "last before cons"
  store.allFCS$Status[which(store.allFCS$Status == "1")] <- "last before cons"
  
  #change 4 patients d22post2ndind status to last bef cons
  store.allFCS$Status[intersect(which(store.allFCS$PatientID == "13201300002"),
                                which(store.allFCS$Status=="d22post2ndind"))] <- "last before cons"
  store.allFCS$Status[intersect(which(store.allFCS$PatientID == "13201300004"),
                                which(store.allFCS$Status=="d22post2ndind"))] <- "last before cons"
  
  store.allFCS$Status[intersect(which(store.allFCS$PatientID == "13201400003"),
                                which(store.allFCS$Status=="d22post2ndind"))] <- "last before cons"
  store.allFCS$Status[intersect(which(store.allFCS$PatientID == "13201400011"),
                                which(store.allFCS$Status=="d22post2ndind"))] <- "last before cons"
  store.allFCS$Status[which(store.allFCS$Status == "day22after2ind")] <- "d22post2ndind"
  
  store.allFCS$Status <- as.factor(store.allFCS$Status)
  
}

tubes <- unique(store.allFCS$Tube)

AML_Status <- read.csv("~/code/AML_Stat.csv")
colnames(AML_Status) <- c("PatientID","AML_Status")
AML_Status[,2] <- as.character(AML_Status[,2])
AML_Status[,2][is.na(AML_Status[,2])] <- "NoData"
AML_Status <- as.data.frame(AML_Status)


store.allFCS.AML <- merge(store.allFCS, AML_Status, by="PatientID")
store.allFCS.AML <- store.allFCS.AML[-which(store.allFCS.AML$AML_Status == "NoData"),]

store.allFCS.AML.normal <- rbind.fill(store.allFCS.AML,store.FCS.normal)

suppressWarnings(dir.create("~/results/Preprocessing/"))
saveRDS(store.allFCS.AML.normal,"~/results/Preprocessing/store.allFCS.AML.normal.rds") 



