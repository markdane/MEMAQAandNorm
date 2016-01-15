library("rmarkdown")
library("ruv")

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


medianNaiveRUV2 <- function(nu, Y, cIdx, k){
  nY <- naiveRandRUV(Y, cIdx, nu, k)
  #Median summarize the normalized values
  #melt matrix to have ECMp and Ligand columns
  nYm <- melt(nY, varnames=c("BWL","SE"))
  #Recreate a MEP column
  nYm$ECMp <- sub("([[:alnum:]]*_){1}","",nYm$SE)
  nYm$Ligand <- sub("([[:alnum:]]*_){2}","",nYm$BWL)
  nYm$MEP <- paste(nYm$ECMp,nYm$Ligand, sep="_")
  nYm$Barcode <- strsplit2(nYm$BWL,split="_")[,1]
  #Median summarize by MEP
  nYm <- data.table(nYm)
  nYm <- nYm[, BWL := NULL]
  nYm <- nYm[, SE := NULL]
  nYm <- nYm[,list(Value = median(value)), by="Barcode,MEP,ECMp,Ligand"]
  return(nYm)
}

madNaiveRUV2 <- function(nu, Y, cIdx, k){
  nY <- naiveRandRUV(Y, cIdx, nu, k)
  #MAD summarize the normalized values
  #melt matrix to have ECMp and Ligand columns
  nYm <- melt(nY, varnames=c("BWL","SE"))
  #Recreate a MEP column
  nYm$ECMp <- sub("([[:alnum:]]*_){1}","",nYm$SE)
  nYm$Ligand <- sub("([[:alnum:]]*_){2}","",nYm$BWL)
  nYm$MEP <- paste(nYm$ECMp,nYm$Ligand, sep="_")
  nYm$Barcode <- strsplit2(nYm$BWL,split="_")[,1]
  #Median summarize by MEP
  nYm <- data.table(nYm)
  nYm <- nYm[, BWL := NULL]
  nYm <- nYm[, SE := NULL]
  nYM <- nYm[,list(Value = mad(value)), by="Barcode,MEP,ECMp,Ligand"]
  return(nYM)
}

medianNaiveReplicateRUV2 <- function(k, Y, cIdx, scIdx){
  nY <- naiveReplicateRUV(Y, cIdx, scIdx, k)$cY
  #Median summarize the normalized values
  #melt matrix to have ECMp and Ligand columns
  nYm <- melt(nY, varnames=c("BWL","SE"))
  #Recreate a MEP column
  nYm$ECMp <- sub("([[:alnum:]]*_){1}","",nYm$SE)
  nYm$Ligand <- sub("([[:alnum:]]*_){2}","",nYm$BWL)
  nYm$MEP <- paste(nYm$ECMp,nYm$Ligand, sep="_")
  nYm$Barcode <- strsplit2(nYm$BWL,split="_")[,1]
  #Median summarize by MEP
  nYm <- data.table(nYm)
  nYm <- nYm[, BWL := NULL]
  nYm <- nYm[, SE := NULL]
  nYm <- nYm[,list(Value = median(value)), by="Barcode,MEP,ECMp,Ligand"]
  return(nYm)
}

madNaiveReplicateRUV2 <- function(k, Y, cIdx, scIdx){
  nY <- naiveReplicateRUV(Y, cIdx, scIdx, k)$cY
  #MAD summarize the normalized values
  #melt matrix to have ECMp and Ligand columns
  nYm <- melt(nY, varnames=c("BWL","SE"))
  #Recreate a MEP column
  nYm$ECMp <- sub("([[:alnum:]]*_){1}","",nYm$SE)
  nYm$Ligand <- sub("([[:alnum:]]*_){2}","",nYm$BWL)
  nYm$MEP <- paste(nYm$ECMp,nYm$Ligand, sep="_")
  nYm$Barcode <- strsplit2(nYm$BWL,split="_")[,1]
  #Median summarize by MEP
  nYm <- data.table(nYm)
  nYm <- nYm[, BWL := NULL]
  nYm <- nYm[, SE := NULL]
  nYM <- nYm[,list(Value = mad(value)), by="Barcode,MEP,ECMp,Ligand"]
  return(nYM)
}

medianRUVIIIPlate <- function(k, Y, M, cIdx){
  #browser()
  nY <- RUVIII(Y, M, cIdx, k)
  #melt matrix to have ECMp and Ligand columns
  nYm <- melt(nY, varnames=c("Barcode","WSE"))
  #Extract ECMp and Well values into separate columns
  nYm$ECMp <- sub("([[:alnum:]]*_){2}","",nYm$WSE)
  nYm$Well <- strsplit2(nYm$WSE,split="_")[,1]
  
  #Median summarize by MEP
  nYm <- data.table(nYm)
  nYm <- nYm[, WSE := NULL]
  nYm <- nYm[,list(Value = median(value)), by="Barcode,Well,ECMp"]
  return(nYm)
}

madRUVIIIPlate <- function(k, Y, M, cIdx){
  #browser()
  nY <- RUVIII(Y, M, cIdx, k)
  #melt matrix to have ECMp and Ligand columns
  nYm <- melt(nY, varnames=c("Barcode","WSE"))
  #Extract ECMp and Well values into separate columns
  nYm$ECMp <- sub("([[:alnum:]]*_){2}","",nYm$WSE)
  nYm$Well <- strsplit2(nYm$WSE,split="_")[,1]
  
  #MAD summarize the normalized values by MEP
  nYm <- data.table(nYm)
  nYm <- nYm[, WSE := NULL]
  nYM <- nYm[,list(Value = mad(value)), by="Barcode,Well,ECMp"]
  return(nYM)
}

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

RUVIII = function(Y, M, ctl, k=NULL, eta=NULL, average=FALSE, fullalpha=NULL)
{
  Y = RUV1(Y,eta,ctl)
  if (is.null(k))
  {
    ycyctinv = solve(Y[,ctl]%*%t(Y[,ctl]))
    newY = (M%*%solve(t(M)%*%ycyctinv%*%M)%*%(t(M)%*%ycyctinv)) %*% Y
    fullalpha=NULL
  }
  else if (k == 0) newY = Y
  else
  {
    m = nrow(Y)
    Y0 = residop(Y,M)
    fullalpha = t(svd(Y0%*%t(Y0))$u[,1:(m-ncol(M)),drop=FALSE])%*%Y
    alpha = fullalpha[1:k,,drop=FALSE]
    ac = alpha[,ctl,drop=FALSE]
    W = Y[,ctl]%*%t(ac)%*%solve(ac%*%t(ac))
    newY = Y - W%*%alpha
  }
  if (average) newY = ((1/apply(M,2,sum))*t(M)) %*% newY
  return(list(newY = newY, fullalpha=fullalpha))
}

medianNaiveRUV2Plate <- function(nu, Y, cIdx, k){
  #browser()
  nY <- naiveRandRUV(Y, cIdx, nu, k)
  #melt matrix to have ECMp and Ligand columns
  nYm <- melt(nY, varnames=c("Barcode","WSE"))
  #Extract ECMp and Well values into separate columns
  nYm$ECMp <- sub("([[:alnum:]]*_){2}","",nYm$WSE)
  nYm$Well <- strsplit2(nYm$WSE,split="_")[,1]
  
  #Median summarize by MEP
  nYm <- data.table(nYm)
  nYm <- nYm[, WSE := NULL]
  nYm <- nYm[,list(Value = median(value)), by="Barcode,Well,ECMp"]
  return(nYm)
}

madNaiveRUV2Plate <- function(nu, Y, cIdx, k){
  #browser()
  nY <- naiveRandRUV(Y, cIdx, nu, k)
  #melt matrix to have ECMp and Ligand columns
  nYm <- melt(nY, varnames=c("Barcode","WSE"))
  #Extract ECMp and Well values into separate columns
  nYm$ECMp <- sub("([[:alnum:]]*_){2}","",nYm$WSE)
  nYm$Well <- strsplit2(nYm$WSE,split="_")[,1]
  
  #MAD summarize the normalized values by MEP
  nYm <- data.table(nYm)
  nYm <- nYm[, WSE := NULL]
  nYM <- nYm[,list(Value = mad(value)), by="Barcode,Well,ECMp"]
  return(nYM)
}
#cellLine
#ss is staining set
#s is signal
#m is method

calcResidual <- function(x){
  mel <- mean(x, na.rm=TRUE)
  return(x-mel)
}

####################################

dataFiles <- data.frame(CellLine=rep(c("PC3", "YAPC","MCF7"),each=1),
                        StainingSet = rep(c("SS2","SS3","SS3","SS3","SS3"), each=3),
                        Signal=rep(c("EdU","LineageRatioLog2","DNA2N","SCC","Ecc"), each=3),
                        inputFileName=rep(c("../MEP-LINCS/PC3/SS2/AnnotatedData/PC3_SS2_Level1.txt", "../MEP-LINCS/YAPC/SS2/AnnotatedData/YAPC_SS2_Level1.txt", "../MEP-LINCS/MCF7/SS2/AnnotatedData/MCF7_SS2_Level1.txt","../MEP-LINCS/PC3/SS3/AnnotatedData/PC3_SS3_Level1.txt", "../MEP-LINCS/YAPC/SS3/AnnotatedData/YAPC_SS3_Level1.txt", "../MEP-LINCS/MCF7/SS3/AnnotatedData/MCF7_SS3_Level1.txt","../MEP-LINCS/PC3/SS3/AnnotatedData/PC3_SS3_Level1.txt", "../MEP-LINCS/YAPC/SS3/AnnotatedData/YAPC_SS3_Level1.txt", "../MEP-LINCS/MCF7/SS3/AnnotatedData/MCF7_SS3_Level1.txt","../MEP-LINCS/PC3/SS3/AnnotatedData/PC3_SS3_Level1.txt", "../MEP-LINCS/YAPC/SS3/AnnotatedData/YAPC_SS3_Level1.txt", "../MEP-LINCS/MCF7/SS3/AnnotatedData/MCF7_SS3_Level1.txt","../MEP-LINCS/PC3/SS3/AnnotatedData/PC3_SS3_Level1.txt", "../MEP-LINCS/YAPC/SS3/AnnotatedData/YAPC_SS3_Level1.txt", "../MEP-LINCS/MCF7/SS3/AnnotatedData/MCF7_SS3_Level1.txt"), each=1),stringsAsFactors = FALSE)

x <- dataFiles

callQALM <- function(x){
  render(paste0("Mep-LINCS_QALM.Rmd"),
         output_file = paste0("Mep-LINCS_QALM_",x[["CellLine"]],"_",x[["Signal"]],".html"),
         output_format = "html_document")
}


apply(dataFiles, 1, callQALM)
