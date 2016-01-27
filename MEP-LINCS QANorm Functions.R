#MEP-LINCS QA Norm functions

naiveRandRUV_MD <- function (Y, cIdx, nuCoeff = 0.001, k = nrow(Y)) 
{
  if (k == 0) {
    nY <- Y
  } else if (k==1){
    svdYc <- svd(Y[, cIdx])
    W <- svdYc$u[, 1:k] * svdYc$d[1]
    nu <- nuCoeff * svdYc$d[1]^2
    nY <- Y - W %*% solve(t(W) %*% W + nu * diag(k), t(W) %*% Y)
  } else {
    svdYc <- svd(Y[, cIdx])
    W <- svdYc$u[, 1:k] %*% diag(svdYc$d[1:k])
    nu <- nuCoeff * svdYc$d[1]^2
    nY <- Y - W %*% solve(t(W) %*% W + nu * diag(k), t(W) %*% Y)
  }
  return(nY)
}

plotLEHmap <- function(dt, fill, titleExpression, limits, xAxisSize=.6, yAxisSize=.5){
  p <- ggplot(dt, aes_string(x="Ligand", y="ECMp", fill = fill))+
    geom_tile()+
    scale_fill_gradient(low="white",high="red",oob = scales::squish,
                        limits=limits)+
    ggtitle(titleExpression)+
    xlab("")+ylab("")+
    guides(fill=FALSE)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(xAxisSize)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(yAxisSize)), plot.title = element_text(size = rel(1)),legend.text=element_text(size = rel(.3)),legend.title=element_text(size = rel(.7)),panel.grid.major = element_line(linetype = 'blank'),panel.background = element_rect(fill = "dimgray"))
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
  nY <- RUVIII(Y, M, cIdx, k)["newY"]
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
  nY <- RUVIII(Y, M, cIdx, k)["newY"]
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

medianRZS <- function(Y){
  #RZS Normalize each row of Y
  nY <- apply(Y,1, function(x){
    xCtrl <- x[grepl("A03", names(x))]
    xMedian <- median(xCtrl,na.rm=TRUE)
    #impute answers if xMAD is 0
    xMad <- mad(xCtrl,na.rm = TRUE)+.01
    xRZS <- (x-xMedian)/xMad
  })
  #melt matrix to have ECMp and Ligand columns
  nYm <- melt(nY, varnames=c("WSE","Barcode"))
  #Extract ECMp and Well values into separate columns
  nYm$ECMp <- sub("([[:alnum:]]*_){2}","",nYm$WSE)
  nYm$Well <- strsplit2(nYm$WSE,split="_")[,1]
  
  #Median summarize by MEP
  nYm <- data.table(nYm)
  nYm <- nYm[, WSE := NULL]
  nYm <- nYm[,list(Value = median(value, na.rm=TRUE)), by="Barcode,Well,ECMp"]
  return(nYm)
}

madRZS <- function(Y){
  #RZS Normalize each row of Y
  nY <- apply(Y,1, function(x){
    xCtrl <- x[grepl("A03", names(x))]
    xMedian <- median(xCtrl,na.rm=TRUE)
    #impute answers if xMAD is 0
    xMad <- mad(xCtrl,na.rm = TRUE)+.01
    xRZS <- (x-xMedian)/xMad
  })
  #melt matrix to have ECMp and Ligand columns
  nYm <- melt(nY, varnames=c("WSE","Barcode"))
  #Extract ECMp and Well values into separate columns
  nYm$ECMp <- sub("([[:alnum:]]*_){2}","",nYm$WSE)
  nYm$Well <- strsplit2(nYm$WSE,split="_")[,1]
  
  #Median summarize by MEP
  nYm <- data.table(nYm)
  nYm <- nYm[, WSE := NULL]
  nYm <- nYm[,list(Value = mad(value)), by="Barcode,Well,ECMp"]
  return(nYm)
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
  #browser()
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


calcResidual <- function(x){
  mel <- median(x, na.rm=TRUE)
  return(x-mel)
}

naiveRUV2Plate <- function(k, nu, Y, cIdx){
  #browser()
  nY <- naiveRandRUV_MD(Y, cIdx, nu, k)
  #melt matrix to have ECMp and Ligand columns
  nYm <- melt(nY, varnames=c("Barcode","WSE"))
  #Extract ECMp and Well values into separate columns
  nYm$ECMp <- sub("([[:alnum:]]*_){2}","",nYm$WSE)
  nYm$Spot <- strsplit2(nYm$WSE,split="_")[,2]
  nYm$Well <- strsplit2(nYm$WSE,split="_")[,1]
  nYm$k <- k
  
  #Median summarize by MEP
  nYm <- data.table(nYm)
  nYm <- nYm[, WSE := NULL]
  
  return(nYm)
}

naiveReplicateRUV2Plate <- function(k, Y, cIdx, scIdx){
  nY <- naiveReplicateRUV(Y, cIdx, scIdx, k)$cY
  #Median summarize the normalized values
  #melt matrix to have ECMp and Ligand columns
  nYm <- melt(nY, varnames=c("BWL","SE"))
  #Recreate a MEP column
  nYm$ECMp <- sub("([[:alnum:]]*_){1}","",nYm$SE)
  nYm$Ligand <- sub("([[:alnum:]]*_){2}","",nYm$BWL)
  nYm$MEP <- paste(nYm$ECMp,nYm$Ligand, sep="_")
  nYm$Barcode <- strsplit2(nYm$BWL,split="_")[,1]
  nYm$k <- k
  #Median summarize by MEP
  nYm <- data.table(nYm)
  nYm <- nYm[, BWL := NULL]
  nYm <- nYm[, SE := NULL]
  nYm$k <- k
  return(nYm)
}

RUVIIIPlate <- function(k, Y, M, cIdx){
  #browser()
  nY <- RUVIII(Y, M, cIdx, k)["newY"]
  #melt matrix to have ECMp and Ligand columns
  nYm <- melt(nY, varnames=c("Barcode","WSE"), as.is=TRUE)
  #Extract ECMp and Well values into separate columns
  nYm$ECMp <- sub("([[:alnum:]]*_){2}","",nYm$WSE)
  nYm$Well <- strsplit2(nYm$WSE,split="_")[,1]
  
  nYm <- data.table(nYm)
  nYm <- nYm[, WSE := NULL]
  nYm$k <- k
  return(nYm)
}

RZSPlate <- function(Y){
  #RZS Normalize each row of Y
  nY <- apply(Y,1, function(x){
    xCtrl <- x[grepl("A03", names(x))]
    xMedian <- median(xCtrl,na.rm=TRUE)
    #impute answers if xMAD is 0
    xMad <- mad(xCtrl,na.rm = TRUE)+.01
    xRZS <- (x-xMedian)/xMad
  })
  #melt matrix to have ECMp and Ligand columns
  nYm <- melt(nY, varnames=c("WSE","Barcode"))
  #Extract ECMp and Well values into separate columns
  nYm$ECMp <- sub("([[:alnum:]]*_){2}","",nYm$WSE)
  nYm$Well <- strsplit2(nYm$WSE,split="_")[,1]
  
  #Median summarize by MEP
  nYm <- data.table(nYm)
  nYm <- nYm[, WSE := NULL]
  return(nYm)
}

callQANorm <- function(x){
  render(paste0("Mep-LINCS_QANorm.Rmd"),
         output_file = paste0("Mep-LINCS_QANorm_",x[["CellLine"]],"_",x[["Signal"]],"_",x[["Method"]],".html"),
         output_format = "html_document") 
}

callSSQANorm <- function(x){
  render(paste0("Mep-LINCS_QANorm.Rmd"),
         output_file = paste0("Mep-LINCS_QANorm_",unique(x[["CellLine"]]),"_",unique(x[["Signal"]]),"_",unique(x[["Method"]]),".html"),
         output_format = "html_document") 
}

convert48ECMGAL <- function(filename, minEffect=.95, maxEffect=1.05){
  #Read in the spot metadata from the gal file
  spotMetadata <- readSpotMetadata(filename)
  #Relabel the column Name to ECMpAnnotID
  setnames(spotMetadata, "Name", "ECMpAnnotID")
  
  #Create a data table with the Actual-Simulated replacements
  ECMReplace <- data.table(Actual=unique(spotMetadata$ECMpAnnotID),Simulated=unique(spotMetadata$ECMpAnnotID))
  ECMReplace$Simulated[!grepl("Fiducial|PBS",ECMReplace$Simulated)]<-sprintf("ECMp%02d", 1:48)
  ECMReplace$ECMpER <- c(1,ECMpER=seq(minEffect, maxEffect,length.out = 48),0)
  
  #replace the actual ECMp names
  setkey(spotMetadata,ECMpAnnotID)
  setkey(ECMReplace,Actual)
  spotMetadata <- spotMetadata[ECMReplace]
  spotMetadata <- spotMetadata[,ECMpAnnotID :=NULL]
  setnames(spotMetadata,"Simulated","ECMpAnnotID")
  return(spotMetadata)
}
