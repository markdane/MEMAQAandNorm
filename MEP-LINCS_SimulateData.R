# Simulate MEMA data to evaluate normalizations and hit selection
#1/2016 Mark Dane

#Start with a simple approach with ECMp and Ligand effect rates
#and random noise
#Maintain the 8 well, 8 plate 700 spt MEMA structure of the LINCS experiments

library("limma")#read GAL file and strsplit2
library("MEMA")#merge, annotate and normalize functions
library("data.table")#fast file reads, data merges and subsetting
library(parallel)
library(RUVnormalize)
library(ruv)
library(ggplot2)
library(plotly)


source("~/Documents/MEP-LINCS/MEPLINCSFunctions.R")


addMarginValues <- function(dt, mValue, MValue){
  #browser()
  dt.l <- dt[,list(m.l = median(get(mValue), na.rm = TRUE),
                   M.l = median(get(MValue), na.rm= TRUE)),
             by="Ligand"]
  
  dte. <- dt[,list(me. = median(get(mValue), na.rm = TRUE),
                   Me. = median(get(MValue), na.rm= TRUE)),
             by="ECMp"]
  
  setkey(dt,Ligand)
  setkey(dt.l,Ligand)
  dtm <- merge(dt,dt.l)
  
  setkey(dtm,ECMp)
  setkey(dte.,ECMp)
  dtmm <- merge(dtm, dte.) 
  return(dtmm)
}

plotCenteredBoxes <- function(dt, value, groupBy, colourBy=NULL, titleExpression, yLimits=NULL){
  
  if(is.null(colourBy))  p <- ggplot(dt,aes_string(x=groupBy, y=value))
  else p <- ggplot(dt,aes_string(x=groupBy, y=value, colour=colourBy))
  
  p <- p + geom_boxplot()+ 
    ggtitle(titleExpression)+
    xlab("")+ylab("")+
    guides(fill=FALSE)+
    coord_cartesian(ylim=yLimits)+
    theme_grey(base_size = 12, base_family = "")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(1)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(1)), plot.title = element_text(size = rel(1)),legend.text=element_text(size = rel(.3)),legend.title=element_text(size = rel(1)))
  suppressWarnings(print(p))
}


plotLEHmap <- function(dt, fill, titleExpression, limits){
  p <- ggplot(dt, aes_string(x="Ligand", y="ECMp", fill = fill))+
    geom_tile()+
    scale_fill_gradient(low="white",high="red",oob = scales::squish,
                        limits=limits)+
    ggtitle(titleExpression)+
    xlab("")+ylab("")+
    guides(fill=FALSE)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(.6)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(.5)), plot.title = element_text(size = rel(1)),legend.text=element_text(size = rel(.3)),legend.title=element_text(size = rel(.7)),panel.grid.major = element_line(linetype = 'blank'),panel.background = element_rect(fill = "dimgray"))
  suppressWarnings(print(p))
}


#Use the phase 2 GAL file to set the layout
# Read and clean spotmetadata

#Read in the spot metadata from the gal file
spotMetadata <- readSpotMetadata("~/Documents/MEP-LINCS/GALFiles/20160120_LI8X001_2.gal")
#Relabel the column Name to ECMpAnnotID
setnames(spotMetadata, "Name", "ECMpAnnotID")

#Create a data table with the Actual-Simulated replacements
ECMReplace <- data.table(Actual=unique(spotMetadata$ECMpAnnotID),Simulated=unique(spotMetadata$ECMpAnnotID))
ECMReplace$Simulated[!grepl("Fiducial|PBS",ECMReplace$Simulated)]<-paste0("ECMp",1:48)
ECMReplace$ECMpER <- c(1,ECMpER=seq(.95,1.05,length.out = 48),0)

#replace the actual ECMp names
setkey(spotMetadata,ECMpAnnotID)
setkey(ECMReplace,Actual)
spotMetadata <- spotMetadata[ECMReplace]
spotMetadata <- spotMetadata[,ECMpAnnotID :=NULL]
setnames(spotMetadata,"Simulated","ECMpAnnotID")

#Make a rotated version of the spot metadata to match the print orientation
setkey(spotMetadata,Spot)
spotMetadata180 <- rotateMetadata(spotMetadata)

#Create the well metadata
wmd <-data.table(Barcode=rep(paste0("LI8S","0000",1:8), each=8*700),
                 Well=rep(wellAN(2,4), each=700),
                 CellLine="SimCell",
                 LigandAnnotID=rep(paste0("Ligand",1:64), each=700),
                 LigandER=rep(seq(.8,1.2,length.out = 64), each=700))
setkey(wmd,Barcode,Well)
#Create A row and B row array metadata 
AWells <- cbind(spotMetadata,wmd)
BWells <- cbind(spotMetadata180,wmd)
#Merge the correct metadata based on the row
l3 <-rbind(AWells[grepl("A",AWells$Well)],
            BWells[grepl("B",BWells$Well)])

#Create short display names
l3$ECMp <- gsub("_.*","",l3$ECMpAnnotID)
l3$Ligand <- gsub("_.*","",l3$LigandAnnotID)
l3$MEP <- paste(l3$ECMp,l3$Ligand,sep = "_")

#Delete the blank spots
l3 <- l3[!grep("PBS|Fiducial",l3$ECMp)]
#Add a convenience label for wells and ligands
l3$WL <- paste(l3$Well,l3$Ligand,sep = "_")

#Calculate the spot cell count based on the ECMp and Ligand effect rates
SCCSeed <- 80
set.seed(1234)
SCCNoise <-sample(-5:5,size=nrow(l3),replace=TRUE)
l3$Spot_PA_SpotCellCount <- as.integer(SCCSeed*l3$ECMpER*l3$LigandER+SCCNoise)

metadataNames <- "ObjectNumber|^Row$|^Column$|Block|^ID$|PrintOrder|Depositions|CellLine|Endpoint|WellIndex|Center|Array|ECMp|Ligand|MEP|Well_Ligand|ImageID|Sparse|Wedge|OuterCell|Spot_PA_Perimeter|Nuclei_PA_Cycle_State|_SE|ReplicateCount|LoessSCC|QAScore|^WL$"

#Save the un-normalized parameters to merge in later
mdDT <- l3[,grep(paste0(metadataNames,"|Barcode|^Well$|^Spot$"),colnames(l3),value=TRUE), with = FALSE]
#Identify parameters to be normalized
normParameters <- grep(metadataNames,colnames(l3),value=TRUE,invert=TRUE)

#Apply RUV-3 normalization to each feature
k=2
nDT <- normRUV3Dataset(l3[,normParameters, with = FALSE], k)

nDT$NormMethod <- "RUV3"
#Merge the normalized data with its metadata
setkey(nDT,Barcode,Well,Spot)
setkey(mdDT,Barcode,Well,Spot)
nmdDT <- merge(nDT,mdDT)

#merge spot level normalized and raw data
setkey(l3, Barcode, Well, Spot)
slDT <- merge(l3[,normParameters, with = FALSE], nmdDT)

p <- create8WellPseudoImage(slDT,pr = "Spot_PA_SpotCellCount",prDisplay = "SCC")
print(p)

DT <- l3[l3$Ligand=="Ligand1"]
p <- ggplot(DT, aes(x=ECMpAnnotID))+
  geom_bar(width=.8)+geom_hline(yintercept = mean(table(DT$ECMpAnnotID)), colour="blue")+
  ggtitle(" \n\nCount of Replicate ECM Proteins In Each MEMA")+
  xlab("Printed ECM Protein")+ylab("Number of spots")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(.8)),axis.title.x = element_text(size=rel(.8)),axis.title.y = element_text(size=rel(.8)),plot.title = element_text(size = rel(1)),legend.text=element_text(size = rel(.3)),legend.title=element_text(size = rel(.5)))

print(p)

#TEMP: Use format in QA Norm
x<-data.frame(Signal="SCC",Method="RUV3")
l3$SCCmel <- as.numeric(median(log2(l3$Spot_PA_SpotCellCount), na.rm=TRUE))

l1MEPs <- l3[,list(mel = median(eval(parse(text=paste0(unique(unique(x[["Signal"]])),"mel")))),
                   Mel = mad(eval(parse(text=paste0(unique(x[["Signal"]]),"mel"))))),
             by="Barcode,Well,Ligand,ECMp,MEP"]
#l1MEPS is the data.table of the MEP level, raw, transformed responses
l1MEPs <- addMarginValues(l1MEPs,"mel","Mel")

#calculate the overall median and MAD for the staining set
m.. <- median(l1MEPs$mel, na.rm = TRUE)
M.. <- median(l1MEPs$Mel, na.rm = TRUE)


#For an array  Unit of study there are 64 'samples'of 594 spots with 8 replicates
#Create a spot assignment that recognizes rotated B row arrays
l3$RSpot <- l3$Spot
l3$RSpot[grepl("B",l3$Well)] <- 701-l3$RSpot[grepl("B",l3$Well)]
#Add in ligand and ECMp names so they will carry through the normalization
l3$BWL <- paste(l3$Barcode,l3$Well,l3$Ligand,sep="_")
l3$SE <- paste(l3$RSpot,l3$ECMp,sep="_")
l3$BW <- paste(l3$Barcode,l3$Well,sep="_")
l3$WSE <- paste(l3$Well, l3$RSpot,l3$ECMp,sep="_")
l3$BWLSE <- paste(l3$Barcode, l3$Well, l3$Ligand, l3$RSpot,l3$ECMp,sep="_")

#Cast to get mel values by barcode_well_ligand rows and spot columns
#logit transform the values that are porportions in the [0,1] range
#Coerce missing values to have near 0 values before transformation
if(unique(unique(x[["Signal"]])) %in% c("EdU", "DNA2N", "Ecc")){
  fill <- log2(.01/(1-.01))
} else if(unique(x[["Signal"]]) %in% c("SCC")){
  fill <- 1
} else if(unique(x[["Signal"]]) %in% c("LineageRatioLog2")){
  fill <- log2(.001)
} else(stop(paste("Need fill value for",unique(x[["Signal"]])," signal"))) 

l3c <- dcast(l3, BWL~SE, value.var=paste0(unique(x[["Signal"]]),"mel"), fill = fill)
#Remove the BWL column and use it as rownames in the matrix
l3m <- l3c[,grep("BWL",colnames(l3c), value=TRUE, invert=TRUE), with=FALSE]
#Y is a numeric matrix of the raw transformed responses 
#of each spot extracted from the l3 dataset
YArray <- matrix(unlist(l3m), nrow=nrow(l3m), dimnames=list(l3c$BWL, colnames(l3m)))


#Setup data with plate as the unit
#There are 8 'samples' with 694 controls
#Order the spot level data by plate, well, spot number

#Cast to get mel values with Barcode rows and well+Spot+ECMp columns
#Use the names to hold the spot contents
#Coerce missing values to have near 0 proliferation signals
l3PlateC <- dcast(l3, Barcode~WSE, value.var=paste0(unique(x[["Signal"]]),"mel"),fill = fill)
#Remove the Barcode column and use it as rownames in the matrix
l3m <- l3PlateC[,grep("Barcode",colnames(l3PlateC), value=TRUE, invert=TRUE), with=FALSE]
YPlate <- matrix(unlist(l3m), nrow=nrow(l3m), dimnames=list(l3PlateC$Barcode, colnames(l3m)))


limits=quantile(l1MEPs[["mel"]], probs=c(.005, .995), na.rm=TRUE)
titleExpression <- bquote(.(unique(x[["Signal"]]))~Medians~(m["e,l"]))
plotLEHmap(l1MEPs, fill="mel", titleExpression, limits)

limits=quantile(l1MEPs[["Mel"]], probs=c(.005, .995), na.rm=TRUE)
titleExpression <- bquote(.(unique(x[["Signal"]]))~MADs~(M["e,l"]))
plotLEHmap(l1MEPs, fill="Mel", titleExpression, limits)
