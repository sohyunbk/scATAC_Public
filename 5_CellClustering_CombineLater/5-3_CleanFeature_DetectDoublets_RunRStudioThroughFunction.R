library("here")
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library("optparse")
library(rlang)
library(ggplot2)
library("RColorBrewer")

## Somehow it does not work when I submitted the job.
option_list = list(
  make_option(c("--WD"), type="character",
              help="WD", metavar="character"),
  make_option(c("--Name"), type="character",
              help="Name", metavar="character"),
  make_option(c("--PreFix"), type="character",
              help="PreFix", metavar="character"),
  make_option(c("--MinT"), type="character",
              help="MinT", metavar="character"),
  make_option(c("--MaxT"), type="character",
              help="MaxT", metavar="character"),
  make_option(c("--nPC"), type="character",
              help="nPC", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
Name <- opt$Name
PreOptions <- opt$PreFix
MinT<- opt$MinT
MaxT <- opt$MaxT
WD <- opt$WD
NumberOfPC <- opt$nPC

Ex <- function(){
  Name <- as.character("1_A619")
  MinT<- as.character(0.01)
  MaxT <- as.character(0.05)
  WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater_UMAP/"
  PreOptions <- "Tn5Cut1000_Binsize500_Mt0.05"
  #NumbeerOfWindow <- as.character(140000)
  #SVDorNMF <-as.character("SVD")
  NumberOfPC <- as.character(100)
}

Run("1_A619","0.01","0.05","/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater_UMAP/",
    "Tn5Cut1000_Binsize500_Mt0.05","100")


##### Clean First ########################################################
Run <- function(Name,MinT,MaxT,WD,PreOptions,NumberOfPC){
setwd(paste0(WD,"/",Name))
obj <- readRDS(paste0(Name,"_",PreOptions,".rds"))
str(obj)
dim(obj$counts)

cell.counts <- Matrix::colSums(obj$counts)
site.freq <- Matrix::rowMeans(obj$counts)

head(cell.counts)
head(site.freq)

tiff(file=paste(Name,"_",PreOptions,"_Distribution.tiff",sep=""),type="cairo")
layout(matrix(c(1:2), ncol=2))
par(mar=c(3,3,1,1))
plot(density(cell.counts), main="log10 cell counts", log="x")
abline(v=1000, col="red")
plot(density(site.freq), main="peak accessibility frequency", log="x")
dev.off()

### remove cells with less than 1,000 open peaks
### The distribution of average peak accessibilities doesnt show any clear (lower-tail) cut-offs,
#therefore, we will use the default thresholds (min.t=0.001, max.t=0.005). The arguments min.t and max.t set the minimum peak accessibility frequency and upper (99.5%) quantile cut-offs, respectively.
#Default: min.t=0.01, max.t=0.05

obj_AfterClean <- cleanData(obj, min.c=100,
                            min.t=as.numeric(MinT), max.t=as.numeric(MaxT), verbose=T)


## 1) Normalization-tfidf

NumberOfWindow <- as.character(140000)
SVDorNMF <-as.character("SVD")
print("Start tfidif normalization")
obj_new <- tfidf(obj_AfterClean)
print("Done with tfidf")
print("Start Reduce Dims")
#print("Done with normalization")
## 2) Reducing dimension
obj_MatrixDecomposition <- reduceDims(obj_new,method=SVDorNMF,n.pcs=as.numeric(NumberOfPC),
                                      num.var=as.numeric(NumberOfWindow),
                                      svd_slotName="MatrixDecomposition")

print("Done with ReduceDims")
str(obj_MatrixDecomposition)
obj_UMAP <- projectUMAP(obj_MatrixDecomposition, verbose=T, k.near=50, m.dist=0.01,
                        svd_slotName="MatrixDecomposition",umap_slotName="UMAP")

##PlotDrawing functions
obj_Cluster_beforeD <- callClusters(obj_UMAP,
                                    res=0.3,
                                    verbose=T,
                                    svd_slotName="MatrixDecomposition",
                                    umap_slotName="UMAP",
                                    cluster_slotName="Clusters",
                                    cleanCluster=F,
                                    e.thresh=5)

NewFileName<-paste0(Name,"_",PreOptions,"_MinT",MinT,"_MaxT",MaxT,"_PC",NumberOfPC)

ggplot(obj_Cluster_beforeD$Clusters, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
  geom_point(size=0.02) + 
  theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave(paste0(NewFileName,"_BeforeRemovingDoublets.pdf"), width=13, height=10)

## 3) Remove Doublets
obj_UMAP$meta$lib_ID <- "1"
head(obj_UMAP$meta)
str(obj_UMAP)
str(obj_UMAP$rdMethod)
obj_detectDoublets <- detectDoublets(obj_UMAP, threads=10, nTrials=5,
                                     nSample=1000, rdMethod = SVDorNMF,
                                     svd_slotName="MatrixDecomposition")
dim(obj_detectDoublets$meta)
table(obj_detectDoublets$meta$doubletscore)

obj_filterDoublets <- filterDoublets(obj=obj_detectDoublets,umap_slotname = "UMAP",
                                     embedding = "UMAP",filterRatio=1.5,
                                     removeDoublets=T, libraryVar="lib_ID",
                                     verbose=TRUE)
str(obj_filterDoublets)
table(obj_filterDoublets$meta$doubletscore)

saveRDS(obj_filterDoublets, file=paste(NewFileName,".rds",sep=""))

obj_Cluster_AfterD <- callClusters(obj_filterDoublets,
                                   res=0.3,
                                   verbose=T,
                                   svd_slotName="MatrixDecomposition",
                                   umap_slotName="UMAP",
                                   cluster_slotName="Clusters",
                                   cleanCluster=F,
                                   e.thresh=5)

ggplot(obj_Cluster_AfterD$Clusters, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
  geom_point(size=0.02) + 
  theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave(paste0(NewFileName,"_AfterRemovingDoublets.pdf"), width=13, height=10)	

## Doublets Histogram
UMAP_Table <- obj_detectDoublets$UMAP
UMAP_Table$DoubletScore <- obj_detectDoublets$meta$doubletscore
ggplot(UMAP_Table, aes(x =DoubletScore)) +
  geom_histogram()
ggsave(paste0(NewFileName,"_DoubletScoreHistogram.pdf")
       , width=12, height=10)

## Tn5 insertion score and doublet score in cluster
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(UMAP_Table$DoubletScore), max(UMAP_Table$DoubletScore)))
ggplot(UMAP_Table, aes(x=umap1, y=umap2,color=DoubletScore)) +
  geom_point(size=0.5) +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  #ylim(-0.6, 0.5)+
  #xlim(-0.6, 0.5)+
  theme_bw()+
  sc
ggsave(paste0(NewFileName,"_DoubletScoreCluster.pdf")
       , width=12, height=10)


##
UMAP_Table$NumberofTn5Insertion <- obj_detectDoublets$meta$log10nSites
head(UMAP_Table)
tail(UMAP_Table)
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(UMAP_Table$NumberofTn5Insertion),
                                      max(UMAP_Table$NumberofTn5Insertion)))
ggplot(UMAP_Table, aes(x=umap1, y=umap2,color=NumberofTn5Insertion)) +
  geom_point(size=0.5) +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  #ylim(-0.6, 0.5)+
  #xlim(-0.6, 0.5)+
  theme_bw()+
  sc

ggsave(paste0(NewFileName,"_NumberofTn5insertionCluster.pdf"),
       width=12, height=10)
}
