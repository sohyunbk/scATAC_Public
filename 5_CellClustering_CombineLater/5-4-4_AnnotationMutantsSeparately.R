library("here")
library(devtools)
library(Seurat)
library(stringr)
load_all('/home/sb14489/Socrates')
library(harmony)
library(symphony)

SampleS <- "relk1"

Ex <- function(){
  Prefix <- "Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100"
  WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater"
  NumbeerOfWindow <- as.character(0)
}


obj_All <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/CombineAll/Combined_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100.rds")
str(obj_All)
dim(subset(obj_All$meta, obj_All$meta$library=="A619_Re1"))
dim(subset(obj_All$meta, obj_All$meta$library=="A619_Re2"))
dim(subset(obj_All$meta, obj_All$meta$library=="rel2_Re1"))
dim(subset(obj_All$meta, obj_All$meta$library=="rel2_Re2"))
dim(subset(obj_All$meta, obj_All$meta$library=="bif3_Re1"))
dim(subset(obj_All$meta, obj_All$meta$library=="bif3_Re2"))
dim(subset(obj_All$meta, obj_All$meta$library=="relk1_Re1"))
dim(subset(obj_All$meta, obj_All$meta$library=="relk1_Re2"))

## Substract Sample

obj <- list()
obj$meta <- subset(obj_All$meta, obj_All$meta$SampleName==SampleS)
head(obj$meta)
head(obj_All$counts)[,c(20000:20010)]
obj$counts <- obj_All$counts[,colnames(obj_All$counts) %in% rownames(obj$meta)]
obj$counts <- obj$counts[Matrix::rowSums(obj$counts)>0,]

str(obj)

########################
## * Blacklist removal
blacklist_r <- read.table("/scratch/sb14489/3.scATAC/0.Data/BlackList/Zm.final_blaclist.Mito_Chloro_Chip.txt")
#blacklist_r <- read.table("/scratch/sb14489/3.scATAC/0.Data/BlackList/Zm.final_blaclist.Mito_Chloro_Chip_CCUpDown500.txt")
head(blacklist_r)
blacklist.gr <- GRanges(seqnames=as.character(blacklist_r$V1),
                        ranges=IRanges(start=as.numeric(blacklist_r$V2),
                                       end=as.numeric(blacklist_r$V3)),
                        names=as.character(blacklist_r$V4))

head(blacklist.gr)
chr.seq.lengths_load <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai")
chr.seq.lengths <- as.numeric(chr.seq.lengths_load$V2)
names(chr.seq.lengths) <- chr.seq.lengths_load$V1
intervals <- tileGenome(seqlengths=chr.seq.lengths, tilewidth=500, cut.last.tile.in.chrom=TRUE)
head(intervals)
intervals <- intervals[-queryHits(findOverlaps(intervals, blacklist.gr, type="any")),]
str(intervals)
regions <- as.data.frame(intervals)
head(regions)
regions <- paste(regions$seqnames, regions$start, regions$end, sep="_")
head(regions)

Temp <- obj$counts[rownames(obj$counts) %in% regions,]
dim(Temp)
dim(obj$counts)
obj$counts <- Temp


obj <- tfidf(obj, doL2=T)
str(obj)
#######################
setwd(WD)
## Make directory
if (file.exists(file.path(SampleS))){
  setwd(file.path(SampleS))
} else {
  dir.create(file.path(SampleS))
  setwd(file.path(SampleS))
}
getwd()
#######################
SVDorNMF <-as.character("SVD")
NumberOfPC <- as.character(100)
NumbeerOfWindow <- as.character(0)
###########################

out <-  paste0(SampleS,"_",Prefix,"_RemoveBLonlyMitoChloroChIP")
#out <-  paste0("Ref_",Prefix,"_RemoveMitoChloroChIP500bpCC")

saveRDS(obj, file=paste0(out,".tfidf.rds"))
#obj <- readRDS()
print("DonewithTfidf")
str(obj)
#head(obj_A619_merged$residuals)
dim(obj$residuals)
#obj_A619_merged <- readRDS(paste0(out,".tfidf.rds"))


# project with NMF -----------------------------------
obj <- reduceDims(obj,method=SVDorNMF,
                              n.pcs=100,
                              cor.max=0.7,
                              num.var=as.numeric(NumbeerOfWindow),
                              verbose=T,
                              scaleVar=T,
                              doSTD=F,
                              doL1=F,
                              doL2=T,
                              refit_residuals=F)

getwd()
head(obj$meta)


harmony_embeddings <- HarmonyMatrix(obj$PCA, meta_data=obj$meta,
                         vars_use="library", do_pca=F,
                         #theta=c(3, 2),
                         sigma=0.1,
                         nclust=30,
                         max.iter.cluster=100,
                         max.iter.harmony=30,
                         return_object=F) ##return_object should be false if I want to get pca
dim(harmony_embeddings)
harmony_embeddings[seq_len(5), seq_len(5)]

obj[['HM_EMBS']] <- harmony_embeddings
str(obj)

obj_UMAP_WithHarmony <- projectUMAP(obj, verbose=T, k.near=50, 
                                    m.dist=0.01, 
                                    svd_slotName="HM_EMBS",
                                    umap_slotName="UMAP")

K <- "50"
RES <- "0.9"

obj_Cluster_WithHarmony <- callClusters(obj_UMAP_WithHarmony, 
                                        res=as.numeric(RES),
                                        verbose=T,
                                        k.near=as.numeric(K),
                                        svd_slotName='HM_EMBS',
                                        umap_slotName="UMAP",
                                        cluster_slotName="Clusters",
                                        cleanCluster=F,
                                        e.thresh=5)




out_final <- paste0(out,"_k",K,"_res",RES)
pdf(paste0(out_final,"_WithHarmony.pdf"), width=10, height=10)
plotUMAP(obj_Cluster_WithHarmony, cluster_slotName="Clusters", cex=0.2)
dev.off()

colorr <- c("#4F96C4","#84f5d9","#DE9A89","#FDA33F","#060878","#d62744","#62a888",
            "#876b58","#800000", "#800075","#e8cf4f","#f7366d","#0bd43d",
            "#deadce","#adafde","#5703ff")
ggplot(obj_Cluster_WithHarmony$Clusters, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave(paste0(out_final,"_WithHarmonyggplotUMAP.pdf"), width=13, height=10)	

saveRDS(obj_UMAP_WithHarmony, file=paste0(out_final,".AfterHarmony.rds"))
#obj_UMAP_WithHarmony <-readRDS("bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.rds")

head(obj_Cluster_WithHarmony$meta)
write.table(obj_Cluster_WithHarmony$Clusters, paste0(out_final,".AfterHarmony.metadata.txt"), quote=F, row.names=T, col.names=T, sep="\t")
head(obj_Cluster_WithHarmony$HM_EMBS)
write.table(obj_Cluster_WithHarmony$HM_EMBS, paste0(out,".AfterHarmony.PCA.txt"), quote=F, row.names=T, col.names=T, sep="\t")

####################################################################################
## Draw the barplot by Replicates 
head(obj_Cluster_WithHarmony$Clusters)
#> head(Plotdata)
#library  Fre     Ratio                         Celltype
#1 A619_Re1  393  9.451659 BundleSheath_VascularSchrenchyma
#2 A619_Re2  667  9.976069 BundleSheath_VascularSchrenchyma
#3 bif3_Re1   88  4.853833 BundleSheath_VascularSchrenchyma
META <- obj_Cluster_WithHarmony$Clusters
ClusterLevel <- levels(factor(obj_Cluster_WithHarmony$Clusters$LouvainClusters))
RE1Meta <- META[which(META$library == paste0(SampleS,"_Re1")),]
RE2Meta <- META[which(META$library == paste0(SampleS,"_Re2")),]
head(RE1Meta)

RE1Plot <- data.frame(t(table(RE1Meta$LouvainClusters)))
RE1Plot$Var1 <- "Re1"
RE1Plot$Ratio <- (RE1Plot$Freq/sum(RE1Plot$Freq))*100
RE2Plot <- data.frame(t(table(RE2Meta$LouvainClusters)))
RE2Plot$Var1 <- "Re2"
RE2Plot$Ratio <- (RE2Plot$Freq/sum(RE2Plot$Freq))*100
PlotData <- rbind(RE1Plot,RE2Plot)
PlotData$Ratio <- round(PlotData$Ratio,2)

library(plyr)
PlotData <- ddply(PlotData, .(Var1),
                  transform, pos = cumsum(Ratio) - (0.5 * Ratio))


ggplot(PlotData, aes(fill=Var2, y=Ratio, x=Var1)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = colorr[1:length(ClusterLevel)]) +
  geom_text(aes(label = paste(round(Ratio,1),"%")), position = position_stack(vjust =  0.5))+
  theme_minimal()
ggsave(paste0(out_final,"_StackedBarplotByRe.pdf"), width=4, height=5)


