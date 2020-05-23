## Prompts user for centre name

readCentreFunc <- function()
{
  centre <- readline(prompt = "Enter the centre name (Sanger/TCP/CIPHE/BCM/Jax):")
  centre <- tolower(centre)
  if(centre == "sanger" | centre == "tcp" | centre == "ciphe" | centre == "bcm" | centre == "jax")
  {
    return(centre)
  }else{
    centre <- readline(prompt = "Incorrect centre name. Please enter again:")
    centre <- tolower(centre)
    return(centre)
  }
}

###############################################################################################################

## Prompts user for panel number

readPanelFunc <- function()
{
  panel <- readline(prompt = "Enter the Panel number (1/2):")
  panel <- as.integer(panel)
  if(panel == 1 | panel == 2)
  {
    return(panel)
  }else{
    panel <- readline(prompt = "Incorrect Panel number. Please enter again (1/2):")
    panel <- as.integer(panel)
    return(panel)
  }
}

##################################################################################################################

# Compensation function ---------------------------------------------------------------------------

compensateIMPC <- function(ff, fname, fpath, centre, panel.no){
  centre <- tolower(centre)
  
  if (centre == "bcm"){
    compMatrix <- dir(fpath, full.names = T, recursive = F, pattern = "Comp.*csv")
    compMatrix <- read.csv(compMatrix, stringsAsFactors = F)
    
    chanNames <- str_match(compMatrix[,1], pattern = '(?<=:: ).*|DAPI-A')
    idx <- which(is.na(chanNames))
    if(length(idx) < 7){
      chanNames[idx] <- compMatrix[idx, 1]
      chanNames <- gsub('-A','', chanNames)
      chanNames <- gsub('.*-','', chanNames)
      chanNames <- gsub('_', '/', chanNames)
      chanNames <- gsub('BV786', 'CD23', chanNames)  # There are actually 2 different BV786 labels in pData(parameters(f))$desc. One of which is CD23, I hope this is correct.
      channel.ind <- sapply(chanNames, function(x){
        if(x == 'CD4'){ 
          ind<- grep('CD4(?!4)', pData(parameters(ff))$desc, perl = T)
        }else{
          ind <- grep(x, pData(parameters(ff))$desc)
          if(length(ind) < 1){
            if(x == 'IA/E'){
              ind <- grep("MHCII", pData(parameters(ff))$desc)
            }
          }
        }
        return(ind)
      })
      channel.ind <- unlist(channel.ind)
    }else{ 
      # For P1, there are four days (150123, 150505, 150512, 160920)
      # where the compensation matrix row/column names have fluorophore names instead of Ab names
      # The following channel.ind values work for those days
      channel.ind <- 7:14
    }
    compMatrix <- compMatrix[,-c(1)]
    colnames(compMatrix) <- colnames(ff)[channel.ind]
    rownames(compMatrix) <- colnames(ff)[channel.ind]
    ff <- compensate(ff, compMatrix)
    
  }else if(centre == "tcp"){
    
    file.names.0909A <- c("PANEL_A_ABIX_223_C7C07009.fcs", "PANEL_A_ABOV_126_C6C06008.fcs",
                          "PANEL_A_ABQX_108_C4C04006.fcs", "PANEL_A_ABRB_82_C2C02004.fcs",
                          "PANEL_A_ABRB_85_C3C03005.fcs",  "PANEL_A_ABRL_111_C5C05007.fcs",
                          "PANEL_A_ABUN_78_C1C01003.fcs",  "PANEL_A_B6NC_C8C08010.fcs", 
                          "FMO_PANEL_A_CD5_D01011.fcs", "FMO_PANEL_A_CD8A_D05015.fcs", 
                          "FMO_PANEL_A_CD25_D03013.fcs", "FMO_PANEL_A_CD44_D02012.fcs", 
                          "FMO_PANEL_A_CD161_D04014.fcs")
    
    file.names.0909B <- c("PANEL_B_ABIX_223_F7F07022.fcs", "PANEL_B_ABOV_126_F6F06021.fcs",
                          "PANEL_B_ABQX_108_F4F04019.fcs", "PANEL_B_ABRB_82_F2F02017.fcs",
                          "PANEL_B_ABRB_85_F3F03018.fcs",  "PANEL_B_ABRL_111_F5F05020.fcs",
                          "PANEL_B_ABUN_78_F1F01016.fcs",  "PANEL_B_B6NC_F8F08023.fcs", 
                          "FMO_PANEL_B_CD5_G07030.fcs", "FMO_PANEL_B_CD11B_G03026.fcs",
                          "FMO_PANEL_B_CD11C_G04027.fcs", "FMO_PANEL_B_CD21-35_G02025.fcs",
                          "FMO_PANEL_B_CD23_G08031.fcs", "FMO_PANEL_B_CD161_G06029.fcs", 
                          "FMO_PANEL_B_LY6C_G01024.fcs", "FMO_PANEL_B_MHCII_G05028.fcs")
    
    file.names.0915 <- c("PANEL_A_ABMX_153_C7C07009.fcs", "PANEL_A_ABNL_105_C3C03005.fcs",
                         "PANEL_A_ABNL_106_C4C04006.fcs", "PANEL_A_ABRB_94_C1C01003.fcs",
                         "PANEL_A_ABSU_98_C2C02004.fcs",  "PANEL_A_ABTV_105_C5C05007.fcs",
                         "PANEL_A_ABTV_108_C6C06008.fcs", "PANEL_A_B6NC_736_C8C08010.fcs",
                         "PANEL_B_ABMX_153_F7F07021.fcs", "PANEL_B_ABNL_105_F3F03017.fcs",
                         "PANEL_B_ABNL_106_F4F04018.fcs", "PANEL_B_ABRB_94_F1F01015.fcs",
                         "PANEL_B_ABSU_98_F2F02016.fcs",  "PANEL_B_ABTV_105_F5F05019.fcs",
                         "PANEL_B_ABTV_108_F6F06020.fcs", "PANEL_B_B6NC_736_F8F08022.labelled.fcs", 
                         "FMO_PANEL_A_CD5_D01031.fcs", "FMO_PANEL_A_CD8A_D05014.fcs", 
                         "FMO_PANEL_A_CD25_D03012.fcs", "FMO_PANEL_A_CD44_D02032.fcs",
                         "FMO_PANEL_A_CD161_D04013.fcs", "FMO_PANEL_B_CD5_G07029.fcs",
                         "FMO_PANEL_B_CD11B_G03025.fcs", "FMO_PANEL_B_CD11C_G04026.fcs",
                         "FMO_PANEL_B_CD21-35_G02024.fcs", "FMO_PANEL_B_CD23_G08030.fcs",
                         "FMO_PANEL_B_CD161_G06028.fcs", "FMO_PANEL_B_LY6C_G01023.fcs",
                         "FMO_PANEL_B_MHCII_G05027.fcs")
    
    file.names.20160526A <- c("PANEL_A_ACWH_243_C1C01003.fcs", "PANEL_A_ACWH_244_C2C02004.fcs",
                              "PANEL_A_ACWH_245_C3C03011.fcs", "PANEL_A_ACWH_247_C4C04007.fcs",
                              "PANEL_A_ACWH_255_C5C05008.fcs", "PANEL_A_B6NC_950_C6C06009.fcs",
                              "PANEL_A_B6NC_955_C7C07010.fcs", "FMO_PANEL_A_CD5_D01005.fcs",
                              "FMO_PANEL_A_CD8A_D05014.fcs", "FMO_PANEL_A_CD25_D03012.fcs",
                              "FMO_PANEL_A_CD44_D02006.fcs", "FMO_PANEL_A_CD161_D04013.fcs")
    
    file.names.20160526B <- c("PANEL_B_ACWH_243_F1F01002.fcs", "PANEL_B_ACWH_244_F2F02003.fcs",
                              "PANEL_B_ACWH_245_F3F03004.fcs", "PANEL_B_ACWH_247_F4F04005.fcs",
                              "PANEL_B_ACWH_255_F5F05006.fcs", "PANEL_B_B6NC_950_F6F06007.fcs",
                              "PANEL_B_B6NC_955_F7F07008.fcs", "FMO_PANEL_B_CD5_G07015.fcs",
                              "FMO_PANEL_B_CD11B_G03011.fcs", "FMO_PANEL_B_CD11C_G04012.fcs",
                              "FMO_PANEL_B_CD21-35_G02010.fcs", "FMO_PANEL_B_CD23_G08017.fcs",
                              "FMO_PANEL_B_CD161_G06014.fcs", "FMO_PANEL_B_LY6C_G01009.fcs",
                              "FMO_PANEL_B_MHCII_G05013.fcs")
    
    file.names.20160623A <- c("PANEL_A_ACRC_43_C1C01003.fcs", "PANEL_A_ACRK_62_C2C02004.fcs", 
                              "PANEL_A_ACRK_65_C3C03005.fcs", "FMO_PANEL_A_CD5_D01007.fcs",
                              "FMO_PANEL_A_CD8A_D05010.fcs", "FMO_PANEL_A_CD25_D03008.fcs",
                              "FMO_PANEL_A_CD44_D02024.fcs", "FMO_PANEL_A_CD161_D04009.fcs")
    
    file.names.20160623B <- c("PANEL_B_ACRC_43_F1F01012.fcs", "PANEL_B_ACRK_62_F2F02013.fcs", 
                              "PANEL_B_ACRK_65_F3F03014.fcs", "FMO_PANEL_B_CD5_G07022.fcs",
                              "FMO_PANEL_B_CD11B_G03018.fcs", "FMO_PANEL_B_CD11C_G04019.fcs",
                              "FMO_PANEL_B_CD21-35_G02017.fcs", "FMO_PANEL_B_CD23_G08023.fcs",
                              "FMO_PANEL_B_CD161_G06021.fcs", "FMO_PANEL_B_LY6C_G01016.fcs",
                              "FMO_PANEL_B_MHCII_G05020.fcs")
    
    
    if(fname %in% file.names.0909A){
      SPILL.matrix <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/Tcp_extra/Compensations/Originals/COMPENSATION_CONTROL_2015-09-09_PanelA.csv", check.names = F)[,-1]
    }else if(fname %in% file.names.0909B){
      SPILL.matrix <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/Tcp_extra/Compensations/Originals/Compensation_Controls_2015-09-09_PanelB.csv", check.names = F)[,-1]
    }else if(fname %in% file.names.0915){
      SPILL.matrix <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/Tcp_extra/Compensations/Originals/COMPENSATION_2015-09-15_PANELA&B.csv", check.names = F)[,-1]
      if(panel.no == 1) # one channel name is different from Panel 1 to 2.
        colnames(SPILL.matrix)[8] <- "BV786-A"
    }else if(fname %in% file.names.20160526A){
      SPILL.matrix <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/Tcp_extra/Compensations/Originals/MATRIX_PANELA_2016-05-26.csv", check.names = F)[,-1]
    }else if(fname %in% file.names.20160526B){
      SPILL.matrix <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/Tcp_extra/Compensations/Originals/MATRIX_PANELB_2016-05-26.csv", check.names = F)[,-1]
      colnames(SPILL.matrix) <- gsub("(.*) ::.*", '\\1', colnames(SPILL.matrix))
      rownames(SPILL.matrix) <- colnames(SPILL.matrix)
    }else if(fname %in% file.names.20160623A){
      SPILL.matrix <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/Tcp_extra/Compensations/Originals/Compensation_Matrix_2016-06-23_PANEL_A.csv", check.names = F)[,-1]
    }else if(fname %in% file.names.20160623B){
      SPILL.matrix <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/Tcp_extra/Compensations/Originals/Compensation_Matrix_2016-06-23_PANEL_B.csv", check.names = F)[,-1]
      colnames(SPILL.matrix) <- gsub("(.*) ::.*", '\\1', colnames(SPILL.matrix))
      rownames(SPILL.matrix) <- colnames(SPILL.matrix)
    }else{
      if(det(ff@description$SPILL) == 1){
        cat(paste0("Check the spillover matrix for ", fname, "it's probably an identity matrix!", "\n"))
        SPILL.matrix <- ff@description$SPILL
      }else{
        SPILL.matrix <- ff@description$SPILL
      }
    }
    ff <- compensate(ff, SPILL.matrix)
    
  }else if(centre == "jax"){
    if(panel.no == 1){
      compMatrix <- dir(fpath, full.names = T, recursive = F, pattern = "CompMatrix_Tmem")
    }else{
      compMatrix <- dir(fpath, full.names = T, recursive = F, pattern = "CompMatrix_APC")
    }
    
    # Reading the file depends on whether the compensation matrix data is saved as a csv file or not
    if(length(grep(".csv", compMatrix)) > 0){ 
      compMatrix <- read.csv(compMatrix, stringsAsFactors = F)
      chanNames <- compMatrix[2, 2:ncol(compMatrix)]
      if(fpath == '/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/Jax/Jax_2017-03-29/PKG - JAXKOMP data files/16-0602  KOMP'){
        compMatrix <- compMatrix[3:(2 + length(chanNames)), 2:ncol(compMatrix)]
      }else{
        compMatrix <- compMatrix[5:(4 + length(chanNames)), 2:ncol(compMatrix)]
      }
      compMatrix <- sapply(compMatrix, as.numeric)
      rownames(compMatrix) <- chanNames
      colnames(compMatrix) <- chanNames
    }else{ 
      chanNames <- read.table(compMatrix, skip = 2, nrows = 1, sep ='\t', stringsAsFactors =F)
      compMatrix <- read.table(compMatrix, sep = '\t', skip = 3, nrows = length(chanNames), stringsAsFactors = F)
      rownames(compMatrix) <- chanNames
      colnames(compMatrix) <- chanNames
    }
    
    index.Remove.row <- which(is.na(compMatrix[,1]))
    index.Remove.col <- which(is.na(compMatrix[1,]))
    if(length(index.Remove.row) != 0 & length(index.Remove.col) != 0){
      compMatrix <- compMatrix[-index.Remove.row,]
      compMatrix <- compMatrix[,-index.Remove.col]
    }
    
    ff <- compensate(ff, compMatrix)
    
  }else if(centre == "ciphe" | centre == "sanger"){  # all other centers (Jax, CIPHE)
    
    if(det(ff@description$SPILL) == 1){
      cat(paste0("Check the spillover matrix, it's probably an identity matrix!", "\n"))
      ff <- compensate(ff, ff@description$SPILL)
    }else{
      ff <-compensate(ff, ff@description$SPILL)
    }
  }
  
  return(ff)
  
}


##################################################################################################

# ############################################################################################
# # removes all cells that exist on the edges of FSC-A, FSC-H, SSC-A, and SSC-H. 

removeMargins<- function(f,chans,sens=1, debris=FALSE,return.ind=F,neg=500, verbose = T)
{
  neg <-cbind(chans,neg)[,2]
  #Size is a vector of size 2, to be passed to mfrow in case of plotting
  data <- exprs(f)
  margins <- c()
  marg.list <-list()
  if(!debris)
  {
    for(chan in chans)
    {
      stain.max <-max(data[,chan])
      margins <- which ( data[, chan] >= stain.max*sens)
      marg.list <- c(marg.list, list(margins))
      data <- data[ -margins, ]
      if(verbose == T){print(paste(length(margins), "margin events in",colnames(f)[chan], "will be removed.",sep =" "))}
    }
    
  }else
  {
    for(i in 1:length(chans))
    {
      stain.min <-min(data[,chans[i]])
      margins <- which ( data[, chans[i]] <= stain.min*sens)
      if (neg[i]<500)
      {
        negs <- which ( data[, chans[i]] < neg[i])
        margins <- negs
      }
      marg.list <- c(marg.list, list(margins))
      if (length(margins)!=0){
        data <- data[ -margins, ]
      }
      if(verbose == T){print(paste(length(margins), "debris events in",colnames(f)[chans[i]], "will be removed.",sep =" "))}
    }
  }
  exprs(f) <- data
  if (!return.ind)
    return(f)
  else
    return(list(frame=f,ind=marg.list))
  
}

############################################################################################
# removes all cells that exist on the edges of FSC-A, FSC-H, SSC-A, and SSC-H. 
removeMarginsAndNegatives <-function(f.temp){
  
  
  for ( q1 in 1:length(colnames(f.temp))){ # 6 because there may be 3 FSC and 3 SSC columns 
    if (length(which(f.temp@exprs[,q1]==262143))>0){f.temp <- f.temp[-which(f.temp@exprs[,q1]==262143)]}
    if (length(which(f.temp@exprs[,q1]<=0))>0)     {f.temp <- f.temp[-which(f.temp@exprs[,q1]<=0)]}
  }
  return(f.temp)
  
}


############################################################################################
#Remove all cell types on the outside of an elipse from the FSC-A and FSC-H plot.
removeDoublets <-function(f.temp, temp.nameU, q){
  if (length(which(f.temp@parameters@data$name=="FSC-H"))>=1) {   
    png ( file = paste("Plots/Figures/",temp.nameU, "_", q, "_B4DoubMarg.png", sep=""))
    plotDens(f.temp, c(1,2))  
    dev.off()
    
    singlet <- flowDensity(f.temp,c(which(f.temp@parameters@data$name=="FSC-A"),which(f.temp@parameters@data$name=="FSC-H")),
                           position=c(T,F),percentile=c(.01,.99),use.percentile=c(T,T),ellip.gate=T)
    f.temp <- getflowFrame(singlet)
    
    png ( file = paste("Plots/Figures/",temp.nameU, "_", q, "_AftDoubMarg.png", sep=""))
    plotDens(f.temp, c(1,2))  
    dev.off()
  }
  return(f.temp)
}



############################################################################################
# Remove all cell types on the outside of an elipse from the FSC-A and FSC-H plot.
removeDoubletsSSC <-function(f.temp, name1, name2, directory, Plt){
  
  if (length(which(f.temp@parameters@data$name=="SSC-H"))>=1) {   
    singlet <- flowDensity(f.temp,c(which(f.temp@parameters@data$name=="SSC-A"),which(f.temp@parameters@data$name=="SSC-H")),
                           position=c(T,F),percentile=c(0.01,.99),use.percentile=c(T,T),ellip.gate=T,scale = 0.999999)
    
    f.temp <- getflowFrame(singlet)
  }
  return(f.temp)
}

#################################################################################
# rearranges the order of the flowFrame
reorderfSet <- function(f,commonMarker,columnPlacement){
  
  NumOfChannels <- length(f@parameters@data$desc)    
  temp.number <- which(f@parameters@data$desc==commonMarker)
  if ( length(temp.number)>1){ temp.number <- temp.number[length(temp.number)]}
  row.order <- NULL
  for(q in 1:NumOfChannels) {row.order <- c(row.order,q)}
  row.order[temp.number]    <- columnPlacement
  row.order[columnPlacement]<- temp.number
  f <- f[,row.order]
  return(f)
}


#################################################################################
# Prints out the time since start_time. Used for optimizing code.
TimeOutput <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
TimeOutput(Sys.Date())


#################################################################################
#SVD reduction for flowType to be used for grouping patients based on reduced matrix
#Kmeans used for grouping patients, any other method can be used

#Author: Mehrnoush Malek
#Date: April 2012

svd.reduction <- function(cell.prop,kmean.centres=3,kmeans.start=1000)
  #cell.prop is a matrix of size 3^k by m, where k is number of markers used in flowType and m is number of samples
{
  
  svdr.Result<-svd(t(cell.prop))
  #Find a threshold where samples get far from others
  x11();plot(sort(svdr.Result$d))
  
  inds<-which(svdr.Result$d > locator()$y)
  acf<-t(cell.prop) %*% (svdr.Result$v[,inds])
  kmeans.cluster<-kmeans(acf, centers=kmean.centres, nstart=kmeans.start)
  x11();plot(data.frame(acf), pch=16,col=kmeans.cluster$cluster)
  return(kmeans.cluster)
}


#######################################################################################

##Finds markers in the FCS file
Find.markers <- function(frame,marker.list)
{
  #Parameters:
  #*frame: a flowFrame in the flowSet
  #**marker.list: A vector of characters
  #Output:
  #*channels.ind: a vector of channels numbers in the frame  corresponding to marker.list
  channels.ind <- unlist(lapply(marker.list, function(x) {
    ind <- grep(x, frame@parameters@data[,2], ignore.case=T)
    ind_store <- ind
    if(length(ind)==0){
      warning(paste (x, "not found, check markers!"))
      return(NA)
    } else {
      if(length(ind)>1) {
        cnt <- 0
        repeat{
          cnt <- cnt + 1
          fs.markers<-unlist(lapply(frame@parameters@data[,2], function(x) unlist(strsplit(x," "))[cnt]))
          ind<-match(x,fs.markers)
          if (is.na(ind))
          {
            fs.markers<-unlist(lapply(frame@parameters@data[,2], function(x) unlist(strsplit(x,"-"))[cnt]))
            ind<-match(x,fs.markers)
            if(!is.na(ind))
              break;
          } else {
            break;
          }
          if(cnt >= 10) {
            
            if (length(ind_store) >= 2){
              ind <- ind_store[1]
              warning(paste (x, "found more than one, choosing first. Check markers!"))
            } else {
              warning(paste (x, "not found, check markers!"))
            }
            break;
          }
        }
      }
    }
    return(ind)
  }))
  names(channels.ind)<-marker.list
  #Removing NAs in the channel vector, as not all centres have all of these markers
  #Note that most centres should have Live/CD4/CD8/CD44/CD62/CD5/CD161/CD25 in their channels
  ind <- which (is.na(channels.ind))
  if (length(ind)!=0)
    channels.ind <- channels.ind[-ind]
  return(channels.ind)
}


###################################################################################################

rotate.data <- function(data, chans=NULL, theta=NULL)
{
  if (class(data)== "flowFrame" & !is.null(chans))
  {
    data.new <- exprs(data)[,chans]
    if (is.null(theta))
    {
      reg.slope <- atan(lm(data.new[,1] ~ data.new[,2])$coefficients[2])
      theta <- pi/2 - reg.slope
    }
    data.new <- data.new %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
    exprs(data)[,chans] <- data.new
  }else{
    data <- data %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
  }
  return(list(data=data,theta=theta))
}



# Written by Quentin

logiclTransformCiphe <- function(flow.frame, markers.transform){
  
  #######################################################################################################
  #
  #
  #
  #
  #
  #######################################################################################################
  ## These two lines were commented out and modfied by Sibyl
  #no.transform <- c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","Time","Flag")
  # markers.transform <- colnames(flow.frame)[colnames(flow.frame)%in%no.transform == FALSE]
  
  list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
  list.index <- gsub("N","", list.index)
  list.index <- gsub("\\$P","", list.index)
  
  if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]]))
  {
    r.values <- unlist(lapply(list.index, function(x) 
      as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
    )	
  } else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]]))
  {
    r.values <- unlist(lapply(list.index, function(x) 
      as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
    )	
  } else 
  {
    r.values <- rep(90, length(list.index))
  }
  
  w.values <- (4.5-log10(262143/abs(r.values)))/2
  w.values[which(w.values<0)] <- 0.5
  w.values[which(is.infinite(w.values))] <- 0.5
  
  for(t in 1:length(markers.transform)){
    lgcl <- logicleTransform(w=w.values[t])
    flow.frame <- transform(flow.frame, transformList(markers.transform[t],lgcl))
  }
  
  return(flow.frame)
}

making.graph <- function(ft,start.path,nodes.lbl, edges.from, edges.to,frame,gates,nodes.id,parent="CD45+")
{
  nxt<- tail(start.path,1)
  pop.path <- start.path
  while(nxt!=as.vector(edges.from[1]))
  {
    prnt <- which(edges.to==nxt)[1]
    nxt<- as.vector(edges.from[prnt])
    pop.path <-c( pop.path ,as.vector(nodes.id[which(nodes.id==nxt)]))
  }
  parent.names <- rev(nodes.lbl[match(pop.path,nodes.id)])
  parent.names <- parent.names[-length(parent.names)]
  pop.path <- rev(pop.path)[-1]
  
  pop.path <- lapply(pop.path,function(p) unlist(strsplit(p,"")))
  #par(mfrow=c(ceiling(sqrt(length(pop.path))),ceiling(sqrt(length(pop.path)))))
  ###You need to find a way to remove the 1/2 from old p when you go to next loop, as ind should only return 1 in each step
  mat <-matrix(c(1:(ceiling(sqrt(length(pop.path)))*2)),nrow =ceiling(sqrt(length(pop.path))) ,byrow = T)
  layout(mat)
  #layout.show()
  for (i1 in 1:length(pop.path))
  {
    print(i1)
    pop.vec <- pop.path[[i1]]
    ind <-which(pop.vec!="0")
    if (length(names(ft@Thresholds)) != 0){
      if (names(ft@Thresholds)[ind]!="gate")
      {
        plotDens(frame, c(names(ft@Thresholds)[ind],"SSC-A"),main=parent.names[i1])
        abline(v=gates[[ind]],lwd=3)
        f.temp <- frame
        if (pop.vec[ind]=="1")
          exprs(f.temp) <- exprs(f.temp)[which(exprs(f.temp)[,names(ft@Thresholds)[ind]]<=gates[[ind]]),]
        if (pop.vec[ind]=="2")
          exprs(f.temp) <- exprs(f.temp)[which(exprs(f.temp)[,names(ft@Thresholds)[ind]]>gates[[ind]]),]
        frame <- f.temp
        if (i1<length(pop.path))
          for (k1 in (i1+1):length(pop.path))
            pop.path[[k1]][ind]<-"0"
        ##You need to apply the gate and update frame for the next loop.
      }
    }else{
      plotDens(frame, names(gates)[ind],main=parent.names[i1],
               ylim=c(min(c(exprs(frame)[,names(gates)[ind]],
                            gates[[ind]]),na.rm = T),
                      max(c(exprs(frame)[,names(gates)[ind]],
                            gates[[ind]]),na.rm = T)),
               xlim=c(min(c(exprs(frame)[,names(gates)[ind]],
                            gates[[ind]] ),na.rm = T),
                      max(c(exprs(frame)[,names(gates)[ind]],
                            gates[[ind]]),na.rm = T)))
      lines(gates[[ind]], lwd=3,col=1)
      if (pop.vec[ind]=="1")
        frame <- getflowFrame(notSubFrame(frame, colnames(gates[[ind]]),position = c(T,T),filter=gates[[ind]]))
      if (pop.vec[ind]=="2")
        frame <- getflowFrame(flowDensity(frame, colnames(gates[[ind]]),position = c(T,T),filter=gates[[ind]]))
      if (i1<length(pop.path))
        for (k1 in (i1+1):length(pop.path))
          pop.path[[k1]][ind]<-"0"
    }
  }
}