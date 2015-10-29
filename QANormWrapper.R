library("rmarkdown")
#cellLine
#ss is staining set
#s is signal
#m is method


dataFiles <- data.frame(CellLine=c("PC3"),
                        StainingSet = c("SS2"),
                        Signal=c("EdU"),
                        Method=c("NaiveReplicateRUV"),
                        inputFileName=c("../MEP-LINCS/PC3/SS2/AnnotatedData/PC3_SS2_Level1.txt"),
                        stringsAsFactors = FALSE)
x <- dataFiles

callQANorm <- function(x){
  render(paste0("Mep-LINCS_QANorm_",x[["CellLine"]],"_",x[["Signal"]],"_",x[["Method"]],".Rmd"),
         output_file = paste0("Mep-LINCS_QANorm_",x[["CellLine"]],"_",x[["Signal"]],"_",x[["Method"]],".html"),
         output_format = "html_document") 
}

apply(dataFiles, 1, callQANorm)
