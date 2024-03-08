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
  make_option(c("--bed"), type="character",
              help="bed file from Tn5 insertion", metavar="character"),
  make_option(c("--ann"), type="character",
              help="GTF", metavar="character"),
  make_option(c("--chr"), type="character",
              help=".fai file", metavar="character"),
  make_option(c("--Name"), type="character",
              help="Sample file", metavar="character"),
  make_option(c("--MinTn5"), type="character",
              help="Sample file", metavar="character")
  );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
WD <- opt$WD
BinSize <- opt$BinSize
bed <- opt$bed
ann <- opt$ann
chr <- opt$chr
Name <- opt$Name
minimumtn5counts <- opt$MinTn5

Example <- function(){
  Name <- as.character("4_relk1")
  WD <- "/scratch/sb14489/3.scATAC/1.MaizeExample/6.CellClustering/"
  BinSize <- as.character("500")
  bed <- "/scratch/sb14489/3.scATAC/1.MaizeExample/4.Bam_FixingBarcode/Ex_4_relk1_Unique.bed"
  ann <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_MtPtAdd_Rsf.gtf"
  chr <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai"
  minimumtn5counts <- "1000"
}

### Make wd and setwd
if (!dir.exists(WD)){
  dir.create(WD)
} else {
  print("Dir already exists!")
}

WDir <- paste0(WD,"/",Name)
if (!dir.exists(WDir)){
  dir.create(WDir)
} else {
  print("Dir already exists!")
}
setwd(WDir)

obj <- loadBEDandGenomeData(bed, ann, chr)  ## Takes long
obj <- countRemoveOrganelle(obj,org_scaffolds=c("Mt","Pt"),remove_reads=TRUE) ## This function only add meta$MtPt and remove reads from bed file in Mtpt

obj <- callACRs(obj, genomesize=2.5e9,
                shift= -50,
                extsize=100,
                fdr=0.05,
                output="bulk_peaks",
                tempdir="./macs2_temp",
                verbose=T) ## It should be run to get the Ratio of Tn5s in ACR
##Q how the peaks can be generated in MtPt although we deleted all the reads

obj <- buildMetaData(obj, tss.window=2000, verbose=TRUE,organelle_scaffolds=c("Mt","Pt"))
str(obj)

## Check the cells MtPt ratio with plot
obj$meta$pPtMt <- obj$meta$ptmt/(obj$meta$total+obj$meta$ptmt) ## I think here, pOrg calculation was wrong so that's why I re calculated.
PlotData <- obj$meta
HighDepthnumber <- length(which(PlotData$total>=1000&PlotData$pPtMt>0.2))
LowDepthnumber <- length(which(PlotData$total<1000&PlotData$pPtMt>0.2))

PlotData$Tn5Number <- paste0("<1000 # Tn5\n (#Cells < 0.2 of rMtPt: ",
                             as.character(LowDepthnumber),")")
PlotData[which(PlotData$total>1000),]$Tn5Number <- paste0(">1000 # Tn5 \n( #Cells < 0.2 of rMtPt: ",
                                                          as.character(HighDepthnumber),")")
#str(obj)
ggplot(PlotData, aes(y=pPtMt,x=Tn5Number)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  ggtitle(paste0("Total barcodeNumber :" ,as.character(length(rownames(PlotData)))))

ggsave(paste0(Name,"_OrgRatio_VioletPlot.pdf"), width=8, height=7)

##########################
##Filter Cells with high MtPt Ratio
Cells_lowMtPt <- rownames(obj$meta)[which(obj$meta$pPtMt<0.2)]
#str(obj)
#head(obj$acr)
obj$meta <-obj$meta[Cells_lowMtPt,]
obj$bed <- obj$bed[obj$bed$V4%in%Cells_lowMtPt,]
#dim(obj$meta)


if (BinSize == "peak"){
  obj <- generateMatrix(obj, filtered=F,peaks=T, verbose=T,
                        organelle_scaffolds =c("Mt","Pt"))
}else{
  obj <- generateMatrix(obj, filtered=F, windows=as.numeric(BinSize),
                        peaks=F, verbose=T,organelle_scaffolds =c("Mt","Pt"))
}


## Save reds up to load data
#NewFileName <- paste0(Name,"_loadData.rds")
#saveRDS(obj, file=NewFileName)

#soc.obj <- convertSparseData(obj, verbose=T)
#isCells <- dget("/home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/5_CellClustering/isCells.R")

str(obj)

obj <- isCell(obj, num.test=16000,
                      num.tn5=NULL,
                      num.ref=5000,
                      background.cutoff=500,
                      min.pTSS=0.2,
                      min.FRiP=0.2,
                      min.pTSS.z= -2,
                      min.FRiP.z= -2,
                      verbose=F) #background.cutoff: Tn5 insertions for bad cells
#ForExData
#obj <- isCell(obj, num.test=100,
#              num.tn5=NULL,
#              num.ref=50,
#              background.cutoff=20000,
#              min.pTSS=0.1,
#              min.FRiP=0.1,
#              min.pTSS.z= -2,
#              min.FRiP.z= -2,
#              verbose=F) #background.cutoff: Maximumnumber of Tn5 insertions for bad cells

#fit <- smooth.spline(rank, depth, spar=0.1)
str(obj)
NewFileName <- paste(Name,"_Tn5Cut",minimumtn5counts,sep="")
pdf(file=paste(NewFileName,".pdf",sep=""),width=14,height=4)
obj$meta$acr <- obj$meta$acrs ##Fix the error from the variable
obj <- findCells(obj,
                 doplot=T,
                 min.cells = 100,
                 max.cells=16000,
                 set.tn5.cutoff=as.numeric(minimumtn5counts), #Override spline fitting to set minimum tn5 cout per cell.
                 min.tn5=1000, #Default:1000
                 tss.z.thresh = 3,
                 tss.min.freq = 0.2,
                 frip.min.freq = 0.35,
                 frip.z.thresh = 2,
                 filt.tss=TRUE,
                 filt.frip=TRUE)

dev.off()
str(obj)
obj <- convertSparseData(obj, verbose=T)

NewFileName <- paste0(Name,"_Tn5Cut",minimumtn5counts,"_Binsize",BinSize,".rds")
saveRDS(obj, file=NewFileName)
