library("here")
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library("optparse")
library(rlang)
library(ggplot2)
library("RColorBrewer")
library(harmony)
library(symphony)

Ex <- function(){
  Dir <- as.character("Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100") #NewDir name
  WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/"
  NumbeerOfWindow <- as.character(0)
}

args <- commandArgs(T)
Dir <- as.character(args[1])
minimumtn5counts <- as.character(args[2])
Minimum_PeakNumberPerCell <- as.character(args[3])
MinT<- as.character(args[4])
MaxT <- as.character(args[5])
SVDorNMF <- as.character(args[6])
NumberOfPC <- as.character(args[7])
NumbeerOfWindow<- as.character(args[8])

# functions --------------------------------------------------------------

DrawFigure_Combined <- function(obj_UMAP,NewFileName) {
  ##  Clustering - labeling by sample
  UMAP_Table <- data.frame(obj_UMAP$UMAP)
  #dim(UMAP_Table)
  #length(obj_UMAP$meta$sampleID)
  UMAP_Table$Sample <- c(obj_UMAP$meta$library)
  head(UMAP_Table)
  tail(UMAP_Table)
  
  ggplot(UMAP_Table, aes(x=umap1, y=umap2,color=Sample)) +
    geom_point(size=0.5) +
    xlab("Dimension 1") +
    ylab("Dimension 2") +
    scale_x_continuous(expand=c(0.02,0)) +
    scale_y_continuous(expand=c(0.02,0)) +
    #ylim(-0.6, 0.5)+
    #xlim(-0.6, 0.5)+
    theme_bw() 
  ggsave(paste0(NewFileName,"_byReplicates.pdf"), width=12, height=10)
  
  ggplot(UMAP_Table ,aes(umap1, umap2, color = Sample)) + 
    geom_jitter(size = .1) + facet_grid(.~Sample)+ theme_bw() 
  ggsave(paste0(NewFileName,"_byReplicatesSeperate.pdf"), width=12, height=10)
  
  UMAP_Table$NumberofTn5Insertion <- obj_UMAP$meta$log10nSites
  head(UMAP_Table)
  tail(UMAP_Table)
  
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  sc <- scale_colour_gradientn(colours = myPalette(100), 
                               limits=c(min(UMAP_Table$NumberofTn5Insertion), max(UMAP_Table$NumberofTn5Insertion)))
  ggplot(UMAP_Table, aes(x=umap1, y=umap2,color=NumberofTn5Insertion)) +
    geom_point(size=0.5) +
    xlab("Dimension 1") +
    ylab("Dimension 2") +
    scale_x_continuous(expand=c(0.02,0)) +
    scale_y_continuous(expand=c(0.02,0)) +
    #ylim(-0.6, 0.5)+
    #xlim(-0.6, 0.5)+
    theme_bw()+sc
  ggsave(paste0(NewFileName,"_NumberofTn5insertionCluster.pdf")
         , width=12, height=10)
  
  }
# extract feature information for symphony integration

#obj <- obj_A619_merged
#num.var <- as.numeric(NumbeerOfWindow)

extractFeatureMetrics <- function(obj, 
                                  num.var){
  
  # internal functions
  RowVar <- function(x) {
    spm <- t(x)
    stopifnot(methods::is(spm, "dgCMatrix"))
    ans <- sapply(base::seq.int(spm@Dim[2]), function(j) {
      if (spm@p[j + 1] == spm@p[j]) {
        return(0)
      }
      mean <- base::sum(spm@x[(spm@p[j] + 1):spm@p[j +
                                                     1]])/spm@Dim[1]
      sum((spm@x[(spm@p[j] + 1):spm@p[j + 1]] - mean)^2) +
        mean^2 * (spm@Dim[1] - (spm@p[j + 1] - spm@p[j]))
    })/(spm@Dim[1] - 1)
    names(ans) <- spm@Dimnames[[2]]
    ans
  }
  
  # get sites used for clustering
  residuals_slotName <- "residuals"
  row.var <- RowVar(obj[[residuals_slotName]])
  row.means <- Matrix::rowMeans(obj[[residuals_slotName]])
  adj.row.var <- loess(row.var ~ row.means)$residuals
  
  names(adj.row.var) <- names(row.var)
  adj.row.var <- adj.row.var[order(adj.row.var, decreasing = T)]
  length(adj.row.var) #213218
  topSites <- names(adj.row.var[adj.row.var > num.var])
  
  features <- data.frame(symbol=topSites, mean=row.means[topSites], stddev=sqrt(row.var[topSites]))
  loadings <- obj$PCA_model$v
  colnames(loadings) <- paste0("PC_",seq(1:ncol(loadings)))
  loadings <- loadings[,obj$PCA_model$keep_pcs]
  rownames(loadings) <- obj$hv_sites
  
  # return
  return(list(features=features, loadings=loadings))
  
}

# new TFIDF methods
tfidf <- function(obj,
                  frequencies=T,
                  log_scale_tf=T,
                  scale_factor=10000,
                  doL2=F,
                  slotName="residuals"){
  
  # set bmat
  bmat <- obj$counts
  
  # hidden functions
  .safe_tfidf       <- function(tf, idf,  block_size=2000e6){
    result = tryCatch({
      result = tf * idf
      result
    }, error = function(e) {
      options(DelayedArray.block.size=block_size)
      DelayedArray:::set_verbose_block_processing(TRUE)
      
      tf = DelayedArray(tf)
      idf = as.matrix(idf)
      
      result = tf * idf
      result
    })
    return(result)
  }
  
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    # "term frequency" method
    tf = t(t(bmat) / Matrix::colSums(bmat))
  } else {
    # "raw count" method
    tf = bmat
  }
  
  # Either TF method can optionally be log scaled
  if (log_scale_tf) {
    if (frequencies) {
      tf@x = log1p(tf@x * scale_factor)
    } else {
      tf@x = log1p(tf@x * 1)
    }
  }
  
  # IDF w/ "inverse document frequency smooth" method
  idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  
  # TF-IDF
  tf_idf_counts = .safe_tfidf(tf, idf)
  
  # do L2?
  if(doL2){
    l2norm <- function(x){x/sqrt(sum(x^2))}
    colNorm <- sqrt(Matrix::colSums(tf_idf_counts^2))
    tf_idf_counts <- tf_idf_counts %*% Diagonal(x=1/colNorm)
  }
  
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  obj[[slotName]] <- Matrix(tf_idf_counts, sparse=T)
  obj$norm_method <- "tfidf"
  
  # return
  return(obj)
}

#############===================================================================
#obj_All <- merged.obj
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

## Substract A619.

obj_A619_merged <- list()
obj_A619_merged$meta <- subset(obj_All$meta, obj_All$meta$SampleName=="A619")
head(obj_A619_merged$meta)
obj_A619_merged$counts <- obj_All$counts[,colnames(obj_All$counts) %in% rownames(obj_A619_merged$meta)]
obj_A619_merged$counts <- obj_A619_merged$counts[Matrix::rowSums(obj_A619_merged$counts)>0,]

str(obj_A619_merged)

########################
## * Blacklist removal
blacklist_r <- read.table("/scratch/sb14489/3.scATAC/0.Data/BlackList/Zm.final_blaclist.Mito_Chloro_Chip.txt")
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
regions[c(10:20)]

Temp <- obj_A619_merged$counts[rownames(obj_A619_merged$counts) %in% regions,]
dim(Temp)
dim(obj_A619_merged$counts)
obj_A619_merged$counts <- Temp


obj_A619_merged <- tfidf(obj_A619_merged, doL2=T)
str(obj_A619_merged)
#######################
setwd(WD)
## Make directory
if (file.exists(file.path(Dir))){
  setwd(file.path(Dir))
} else {
  dir.create(file.path(Dir))
  setwd(file.path(Dir))
}
getwd()
#######################
SVDorNMF <-as.character("SVD")
NumberOfPC <- as.character(100)
NumbeerOfWindow <- as.character(0)
###########################

out <-  paste0("Ref_RemoveBLonlyMitoChloroChIP")
#out <-  paste0("Ref_",Prefix,"_RemoveMitoChloroChIP500bpCC")

saveRDS(obj_A619_merged, file=paste0(out,".tfidf.rds"))

print("DonewithTfidf")
str(obj_A619_merged)
#head(obj_A619_merged$residuals)
dim(obj_A619_merged$residuals)
#obj_A619_merged <- readRDS(paste0(out,".tfidf.rds"))


# project with NMF -----------------------------------
obj_A619_merged <- reduceDims(obj_A619_merged,method=SVDorNMF,
                              n.pcs=NumberOfPC,
                              cor.max=0.7,
                              num.var=as.numeric(NumbeerOfWindow),
                              verbose=T,
                              scaleVar=T,
                              doSTD=F,
                              doL1=F,
                              doL2=T,
                              refit_residuals=F)
                              #,svd_slotName="PCA"
                              #num.var=0



getwd()
#saveRDS(obj_A619_merged, file=paste0(out,".tfidf_ReducedD.rds"))
#obj_A619_merged <- readRDS(paste0(out,".tfidf_ReducedD.rds"))
###===== looks okay up to here ==========
str(obj_A619_merged)
# extract feature loadings and var/mean of tfidf
#saveRDS(obj_A619_merged, file=paste0(out,".tfidf_ReducedD.rds")) ## Should start from here on Oct 25th.
#obj_A619_merged <- readRDS(paste0(out,".tfidf_ReducedD.rds"))

feat.data <- extractFeatureMetrics(obj_A619_merged, 0)
str(feat.data)
ids <- rownames(obj_A619_merged$PCA)

# remove batch effects with harmony --------------------------------------
ref.obj <- HarmonyMatrix(obj_A619_merged$PCA, meta_data=obj_A619_merged$meta, 
                         vars_use="library", do_pca=F,
                         #theta=c(3, 2), 
                         sigma=0.1, 
                         nclust=30,
                         max.iter.cluster=100,
                         max.iter.harmony=30,
                         return_object=T)


# create compressed harmony reference
getwd()

sym.ref <- symphony::buildReferenceFromHarmonyObj(ref.obj,
                                                  obj_A619_merged$meta,
                                                  feat.data$features,
                                                  feat.data$loadings,               
                                                  verbose = TRUE,         
                                                  do_umap = TRUE,       
                                                  umap_min_dist = 0.01,
                                                  save_uwot_path = paste0(out,'_uwot_model'))

obj_A619_merged$PCA <- t(ref.obj$Z_corr)
colnames(obj_A619_merged$PCA) <- paste0("PC_", 2:(ncol(obj_A619_merged$PCA)+1))
rownames(obj_A619_merged$PCA) <- ids

saveRDS(sym.ref, file=paste0(out,".symphony.reference.rds"))
saveRDS(obj_A619_merged, file=paste0(out,".afterHarmony.processed.rds"))
getwd()

## Here important!!

#sym.ref<- readRDS("Ref_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100.symphony.reference.rds")
#obj_A619_merged<- readRDS("Ref_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100.afterHarmony.processed.rds")


# reduce to 2-dimensions with UMAP ---------------------------------------
obj_A619_merged$UMAP <- sym.ref$umap$embedding
rownames(obj_A619_merged$UMAP) <- rownames(obj_A619_merged$PCA)
colnames(obj_A619_merged$UMAP) <- c("umap1","umap2")
str(obj_A619_merged)
#head(obj_A619_merged$counts)



# identify clusters using neighborhood graph -----------------------------
obj_A619_merged_Cluster <- callClusters(obj_A619_merged, 
                        res=0.4,
                        k.near=10,
                        verbose=T,
                        cleanCluster=F,
                        cl.method=2,
                        e.thresh=5)
                        #m.clst=100

pdf(paste0(out,"_res0.4_knear10.pdf"), width=10, height=10)
plotUMAP(obj_A619_merged_Cluster, cluster_slotName="Clusters", cex=0.2)
dev.off()


#str(obj_A619_merged)
getwd()
# plot cluster membership on UMAP embedding ------------------------------

DrawFigure_Combined(obj_A619_merged_Cluster,out)


# output text files
meta <- obj_A619_merged_Cluster$Clusters
rd <- obj_A619_merged_Cluster$PCA

# write data
write.table(meta, file=paste0(out, ".REF_CELLs.metadata.txt"), quote=F, row.names=T, col.names=T, sep="\t")
write.table(rd, file=paste0(out, ".REF_CELLs.reduced_dimensions.txt"), quote=F, row.names=T, col.names=T, sep="\t")


############ Re-writing 
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode")
bed <- read.table("1_A619_Unique.bed")
Meta<- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/Ref/Ref_withSub_withName.metadata.txt",header=TRUE)
head(Meta)
head(bed)
Temp <- bed[bed$V4 %in% rownames(Meta),]
dim(Temp)
dim(bed)
bed <- read.table("1_A619_2_Unique.bed")
Temp2 <- bed[bed$V4 %in% rownames(Meta),]

Com <- rbind(Temp,Temp2)
write.table(Com, "1_A619_2_Combined_Filtered.bed", quote=F, row.names=F, col.names=F, sep="\t")


Tn5_Grange <-  GRanges(seqnames = obj$bed$V1,
                       ranges = IRanges(start = obj$bed$V2,
                                        end = obj$bed$V3,
                                        names = obj$bed$V4))

## load previous annotation....

## Select Features
#library("here")
#library(devtools)
#library(Seurat)
#load_all('/home/sb14489/Socrates')
#library("optparse")
#library(rlang)
#library(ggplot2)
#obj <- list()
#ann <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_MtPtAdd_Rsf.gtf"
#obj$gff <- suppressWarnings(suppressMessages(makeTxDbFromGFF(as.character(ann), format="gtf", dbxrefTag="Parent")))
#str(obj$gff)
#if (FeatureName == "gene"){
#  Ann <- genes(obj$gff)
#} else if(FeatureName == "transcript"){
#  Ann <- transcripts(obj$gff)
#} else{message("ERROR: Feature name should be 'gene' or 'transcript")}

#head(Ann)
## Broaden annotation
#nRange <- 500
#Start_new <- c((ranges(Ann)+nRange)@start)-1 ## To convert to the bed format coordinate
#Width_new <- c((ranges(Ann)+nRange)@width)+1 ## To convert to the bed format coordinate
#BroadRange_Ann <- GRanges(seqnames=Ann@seqnames,
#                          ranges= 
#                            IRanges(start=Start_new,
#                                    width=Width_new,
#                                    names=names(Ann)))
#head(BroadRange_Ann)
#obj$ann_broad <- BroadRange_Ann