library("here")
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library("optparse")
library(rlang)
library(ggplot2)
option_list = list(
  make_option(c("--WD"), type="character", 
              help="WD", metavar="character"),
  make_option(c("--BinSize"), type="character", 
              help="Binsize it should be 100, 500, 1000bp or peak", metavar="character"),
  make_option(c("--Name"), type="character", 
              help="Sample file", metavar="character"),
  make_option(c("--MinTn5"), type="character", 
              help="Tn5 cutoff", metavar="character"),
  make_option(c("--MtRatio"), type="character", 
              help="MtRatio", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
Name <- opt$Name
#Name<- "4_relk1_2"
#setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/6.CellClustering/AfterMtMapping/4_relk1_2")
#obj <- readRDS("4_relk1_2_Tn5Cut1000_Binsize500.rds")

#setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/6.CellClustering/AfterMtMapping/1_A619")
#obj <- readRDS("1_A619_Tn5Cut1000_Binsize500.rds")
#Name<-"1_A619"

setwd(paste0(opt$WD,"/",opt$Name))
obj <- readRDS(paste0(opt$Name,"_Tn5Cut",opt$MinTn5,"_Binsize",opt$BinSize,".rds"))
## Recheck QC and confirm QC with isCell
## Draw Figures

#Plotdata <- obj$meta
#Plotdata$IsCell <- "Pass"
#Plotdata[which(Plotdata$is_cell==0),]$IsCell <- "Not-Pass"

#ggplot(Plotdata, aes(y=log10(total),x=IsCell)) + 
#  geom_violin()+
#  geom_jitter(shape=16, position=position_jitter(0.2))+
#  ggtitle(paste0("is Cell passed :" ,as.character(length(which(obj$meta$is_cell==1)))))

#ggsave(paste0(Name,"_isCell_VioletPlot.pdf"), width=8, height=7)	

## 2) Re-check mitochonrial
Plotdata <- obj$meta
Plotdata$Tn5_10000 <- ">10000ofTn5"
Plotdata[which(Plotdata$total<10000),]$Tn5_10000 <- "<10000ofTn5"

ggplot(Plotdata, aes(y=pPtMt,x=Tn5_10000)) + 
  geom_violin()+ ylab("Organelle Ratio")
  #geom_jitter(shape=16, position=position_jitter(0.2))
ggsave(paste0(Name,"_MtPt_VioletOnly_Final.pdf"), width=4, height=3)	

## Filter the cells with isCell 0 and Organelle ratio less than 0.1
#opt$MtRatio <- "0.1"
#opt$Name <- "4_relk1_2"
#opt$MinTn5<-"1000"
#opt$BinSize <-"500"

SelectedCells <- rownames(obj$meta)[which(obj$meta$is_cell==1 & obj$meta$pPtMt<as.numeric(opt$MtRatio))]
obj$meta <- obj$meta[SelectedCells,]
obj$counts <- obj$counts[,SelectedCells]
obj$counts <- obj$counts[Matrix::rowSums(obj$counts) > 0,]

NewFileName <- paste0(opt$Name,"_Tn5Cut",opt$MinTn5,"_Binsize",opt$BinSize,
                      "_Mt",opt$MtRatio,".rds") 
saveRDS(obj, file=NewFileName)

