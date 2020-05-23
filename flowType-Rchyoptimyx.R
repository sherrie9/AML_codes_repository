library(flowType)
library(flowCore)
library(flowTypeFilter)
library(RchyOptimyx)

source("~/code/helperFunc_S.R")
source("~/code/helperFunc.R")
source("~/code/IMPC_Functions.R")


days <- c("Day22","last before 2nd ind","last before cons") 

for(d in 1:length(days)){
  load(file = paste0("~/results/GatingSet_tube",tubes[t],"_",days[1], ".RData"))
  
  markers_list <- as.character(singlets.fs[[1]]@parameters@data[1:ncol(singlets.fs[[1]]), 2])
  markers_list <- Find.markers(singlets.fs[[1]],na.omit(markers_list))
  
  MarkerNames <- c('HLADR-mast','CD13posBCD33neg','CD45',"CD34negstroma","CD33posHLADRneg",
                   "CD11B","SSC","CD117","CD34","singletsmonocytehighSSC")
  
  PropMarkers <- c('V500-A','APC-H7/Cy7-A','SSC-A','PE-Cy7-A','PerCP-A')
  
  
  ResList <- lapply(1:length(singlets.fs),function(i){
    Res <- flowTypeFilter::flowTypeFilter(singlets.fs[[i]],PropMarkers = PropMarkers,Methods="Filters",
                                          MaxMarkersPerPop = 5,PartitionsPerMarker = 2, 
                                          Thresholds=Gates.list[[i]], MarkerNames = MarkerNames)
    return(Res)
  })
  
  Proportions <-ResList[[1]]@CellFreqs
  Proportions <- Proportions / max(Proportions)
  # PropMarkers = c(3,6,7,8,9)
  names(Proportions) <- unlist(lapply(ResList[[1]]@PhenoCodes, 
                                      function(x){return(flowTypeFilter::decodePhenotype(
                                        x,MarkerNames, partitions.per.marker = rep(2,length(MarkerNames)),
                                        # , 
                                        c(F,F,T,F,F,T,T,T,T,F)
                                        # ResList[[1]]@PartitionsPerMarker
                                      ))}))
  
}




index=order(Proportions,decreasing=TRUE)[1:30]
bp=barplot(Proportions[index], axes=FALSE, names.arg=FALSE)
text(bp+0.2, par("usr")[3]+0.02, srt = 90, adj = 0, labels = names(Proportions[index]), xpd = TRUE, cex=0.8)
axis(2);
axis(1, at=bp, labels=FALSE);
title(xlab='Phenotype Names', ylab='Cell Proportion')

#----------------Pvalues from flowType-----------------------------------------------
all.proportions <- matrix(0,length(ResList[[1]]@CellFreqs),length(lv2[,1]))
for (i in 1:length(ResList)){
  all.proportions[,i] = ResList[[i]]@CellFreqs / ResList[[i]]@CellFreqs[1]
} 

Pvals <- vector()
EffectSize <- vector()
Labels <- lv2$AML_Status

for (i in 1:dim(all.proportions)[1]){
  
  if (length(which(all.proportions[i,]!=1))==0){
         Pvals[i]=1
         EffectSize[i]=0
         next
  }
  temp = wilcox.test(all.proportions[i,Labels=="R"], 
                     all.proportions[i,Labels=="0"], 
                     paired = FALSE, alternative = "two.sided")
  # temp=t.test(all.proportions[i, Labels=="R"],
  #             all.proportions[i, Labels=="Normal"])
  Pvals[i] <- temp$p.value 
  # EffectSize[i] <- abs(temp$statistic)
  
}

Pvals[is.nan(Pvals)]=1

names(Pvals)=names(Proportions)
Pvals1 <- vector()
for (i in 1:dim(all.proportions)[1]){
  
  if (length(which(all.proportions[i,]!=1))==0){
    Pvals1[i]=1
    EffectSize[i]=0
    next
  }
  temp = wilcox.test(all.proportions[i,Labels=="R"], 
                     all.proportions[i,Labels=="Normal"], 
                     paired = FALSE, alternative = "two.sided")
  # temp=t.test(all.proportions[i, Labels=="R"],
  #             all.proportions[i, Labels=="Normal"])
  Pvals1[i] <- temp$p.value 
  # EffectSize[i] <- abs(temp$statistic)
  
}
selected <- which(Pvals<0.006)
# names(Pvals1)=names(Proportions)
# selected1 <- which(Pvals1<0.05)
print(names(selected))
# selected <- intersect(selected,selected1)

rownames(all.proportions) <- names(Proportions)

select.mat <- all.proportions[selected,]

samplemedians <- matrix(0,nrow = dim(select.mat)[1], ncol = 3)
for (i in 1:dim(select.mat)[1]){
  samplemedians[i,] <- c(median(select.mat[i,Labels.faust =="R"]),
                         median(select.mat[i,Labels.faust =="0"]),
                         median(select.mat[i,Labels.faust =="Normal"]))
}
colnames(samplemedians) <- c("Relapsed","Non-relapsed","Normal")
rownames(samplemedians) <- rownames(select.mat)

data.plot <- data.frame(samplemedians)
# library(tidyverse)
data.plot <-  data.plot %>% rownames_to_column(var="Cell_Type") %>% gather(key="Status",value = "Median",Relapsed,Non.relapsed,Normal)


data.plot$Status <- as.factor(data.plot$Status)

ggplot(data=data.plot, aes(x=Cell_Type,y=Median, col=Status)) + geom_point() +
  scale_color_manual(values=c("green","blue","red"))+
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust=1, size = 8))

ggsave(filename = paste0("~/results/samplemedianbycelltype_flowtype.pdf"), height = 10,width= 15)


#----------------Pvalues from faust-----------------------------------------------
countMatrix <- readRDS(file=paste0("~/results/Faust-Singlets_",days[2],"/faustData/faustCountMatrix.rds"))
countMatrix <- countMatrix[,-ncol(countMatrix)]
lv2.faust <-  lv2[rownames(countMatrix),]
Labels.faust <- lv2.faust$AML_Status
countMatrix <- t(countMatrix)

#-------------------normalizing---------------------------------------------------
inds <- sapply(colnames(countMatrix),function(rn)grep(rn,rownames(lv2)))

singlets.faust.fs <- list(singlets.fs[[1]],singlets.fs[[5]],singlets.fs[[6]],singlets.fs[[13]],
                          singlets.fs[[11]],singlets.fs[[9]],singlets.fs[[2]],singlets.fs[[3]],
                          singlets.fs[[15]],singlets.fs[[20]],singlets.fs[[19]],singlets.fs[[18]])


for (i in 1:length(singlets.faust.fs)){
  countMatrix[,i] = countMatrix[,i] / nrow(singlets.faust.fs[[i]]@exprs)
} 


Pvals.faust <- vector()
for (i in 1:dim(countMatrix)[1]){
  

  temp = wilcox.test(countMatrix[i,Labels=="R"], 
                     countMatrix[i,Labels=="0"], 
                     paired = FALSE, alternative = "two.sided")

 
  Pvals.faust[i] <- temp$p.value 

  
}
Pvals.faust1 <- vector() 
for (i in 1:dim(countMatrix)[1]){
  

  temp = wilcox.test(countMatrix[i,Labels=="R"], 
                     countMatrix[i,Labels=="Normal"], 
                     paired = FALSE, alternative = "two.sided")
  
  
  Pvals.faust1[i] <- temp$p.value 
  
  
}
names(Pvals.faust) <- rownames(countMatrix)
selected <- which(Pvals.faust<0.05)
names(Pvals.faust1) <- rownames(countMatrix)
selected1 <- which(Pvals.faust1 < 0.05)

selected <- intersect(selected,selected1)

countMatrix <- countMatrix[selected,]

samplemedians <- matrix(0,nrow = dim(countMatrix)[1], ncol = 3)
for (i in 1:dim(countMatrix)[1]){
  samplemedians[i,] <- c(median(countMatrix[i,Labels.faust =="R"]),
                         median(countMatrix[i,Labels.faust =="0"]),
                         median(countMatrix[i,Labels.faust =="Normal"]))
}

colnames(samplemedians) <- c("Relapsed","Non-relapsed","Normal")
rownames(samplemedians) <- rownames(countMatrix)

data.plot <- data.frame(samplemedians)
library(tidyverse)
data.plot <-  data.plot %>% rownames_to_column(var="Cell_Type") %>% gather(key="Status",value = "Median",Relapsed,Non.relapsed,Normal)


data.plot$Status <- as.factor(data.plot$Status)

ggplot(data=data.plot, aes(x=Cell_Type,y=Median, col=Status)) + geom_point() +
  scale_color_manual(values=c("green","blue","red"))+
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust=1, size = 8))

ggsave(filename = paste0("~/results/samplemedianbycelltype_faust.pdf"), height = 10,width= 15)

#-----------Rchyoptimyx------
res<-RchyOptimyx(ResList[[1]]@PhenoCodes, -log10(Pvals),             
                   ResList[[1]]@PhenoCodes[selected[38]], factorial(6),FALSE)


plot(res, phenotypeScores=-log10(Pvals), 
     phenotypeCodes=ResList[[1]]@PhenoCodes, 
     marker.names=MarkerNames, ylab='-log10(Pvalue)')


#--------------------------------------------------------------------------
res<-RchyOptimyx(pheno.codes=ResList[[1]]@PhenoCodes, phenotypeScores=-log10(Pvals), 
                 startPhenotype=ResList[[1]]@PhenoCodes[4287], 2, FALSE)
plot(res, phenotypeScores=-log10(Pvals), phenotypeCodes=ResList[[1]]@PhenoCodes,
     marker.names=MarkerNames, ylab='-log10(Pvalue)')

