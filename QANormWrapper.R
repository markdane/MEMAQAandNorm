
source("MEP-LINCS QANorm Functions.R")

dataFiles <- data.frame(CellLine=rep(c("PC3", "YAPC","MCF7"),each=5),
                        StainingSet = rep(c("SS2","SS3","SS3","SS1","SS2","SS3","SS3"), each=15),
                        Signal=rep(c("EdU","LineageRatioLog2","DNA2N","SCC","SCC","SCC","Ecc"), each=15),
                        Method=c("RZS","NaiveRandRUV", "NaiveReplicateRUV","RUV3", "RUV3"),
                        Unit =c("Plate","Plate","Array","Plate","ArrayWithResiduals"),
                        inputFileName=rep(c("../MEP-LINCS/PC3/SS2/AnnotatedData/PC3_SS2_Level1.txt", "../MEP-LINCS/YAPC/SS2/AnnotatedData/YAPC_SS2_Level1.txt", "../MEP-LINCS/MCF7/SS2/AnnotatedData/MCF7_SS2_Level1.txt","../MEP-LINCS/PC3/SS3/AnnotatedData/PC3_SS3_v1.1_Level1.txt", "../MEP-LINCS/YAPC/SS3/AnnotatedData/YAPC_SS3_Level1.txt", "../MEP-LINCS/MCF7/SS3/AnnotatedData/MCF7_SS3_Level1.txt","../MEP-LINCS/PC3/SS3/AnnotatedData/PC3_SS3_v1.1_Level1.txt", "../MEP-LINCS/YAPC/SS3/AnnotatedData/YAPC_SS3_Level1.txt", "../MEP-LINCS/MCF7/SS3/AnnotatedData/MCF7_SS3_Level1.txt","../MEP-LINCS/PC3/SS1/AnnotatedData/PC3_SS1_v1.3_Level1.txt", "../MEP-LINCS/YAPC/SS1/AnnotatedData/YAPC_SS1_Level1.txt", "../MEP-LINCS/MCF7/SS1/AnnotatedData/MCF7_SS1_Level1.txt","../MEP-LINCS/PC3/SS2/AnnotatedData/PC3_SS2_v1.3_Level1.txt", "../MEP-LINCS/YAPC/SS2/AnnotatedData/YAPC_SS2_Level1.txt", "../MEP-LINCS/MCF7/SS2/AnnotatedData/MCF7_SS2_Level1.txt","../MEP-LINCS/PC3/SS3/AnnotatedData/PC3_SS3_v1.1_Level1.txt", "../MEP-LINCS/YAPC/SS3/AnnotatedData/YAPC_SS3_Level1.txt", "../MEP-LINCS/MCF7/SS3/AnnotatedData/MCF7_SS3_Level1.txt","../MEP-LINCS/PC3/SS3/AnnotatedData/PC3_SS3_v1.1_Level1.txt", "../MEP-LINCS/YAPC/SS3/AnnotatedData/YAPC_SS3_Level1.txt", "../MEP-LINCS/MCF7/SS3/AnnotatedData/MCF7_SS3_Level1.txt"), each=5),stringsAsFactors = FALSE)[seq(5,105,5),]

x <- dataFiles

apply(dataFiles, 1, callQANorm)

#Setup and run with unit as staining set
ssDataFiles <- data.frame(CellLine=rep(c("PC3", "YAPC","MCF7"),each=3),
                          StainingSet = c("SS1","SS2","SS3"),
                          Signal=rep(c("SCC")),
                          Method=c("NaiveReplicateRUVSS"),
                          Unit ="StainingSet",
                          inputFileName=c("../MEP-LINCS/PC3/SS1/AnnotatedData/PC3_SS1_Level1.txt",
                                          "../MEP-LINCS/PC3/SS2/AnnotatedData/PC3_SS2_Level1.txt",
                                          "../MEP-LINCS/PC3/SS3/AnnotatedData/PC3_SS3_v1.1_Level1.txt",
                                          "../MEP-LINCS/YAPC/SS1/AnnotatedData/YAPC_SS1_Level1.txt",
                                          "../MEP-LINCS/YAPC/SS2/AnnotatedData/YAPC_SS2_Level1.txt",
                                          "../MEP-LINCS/YAPC/SS3/AnnotatedData/YAPC_SS3_Level1.txt",
                                          "../MEP-LINCS/MCF7/SS1/AnnotatedData/MCF7_SS1_Level1.txt",
                                          "../MEP-LINCS/MCF7/SS2/AnnotatedData/MCF7_SS2_Level1.txt",
                                          "../MEP-LINCS/MCF7/SS3/AnnotatedData/MCF7_SS3_Level1.txt"), stringsAsFactors = FALSE)


# for(cellLine in unique(ssDataFiles$CellLine)){
#   callSSQANorm(x=ssDataFiles[ssDataFiles$CellLine==cellLine,])
# }
