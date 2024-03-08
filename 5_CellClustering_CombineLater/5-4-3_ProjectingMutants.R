library("here")
library(devtools)
library(Seurat)
library(stringr)
load_all('/home/sb14489/Socrates')
library(harmony)
library(symphony)

Ex <- function(){
  #Dir <- as.character("4_CombineAll_AfterD") #NewDir name
  minimumtn5counts <- as.character(10000)
  Minimum_PeakNumberPerCell <- as.character(100)
  MinT<- as.character(0.01)
  MaxT <- as.character(0)
  SVDorNMF <-as.character("SVD")
  NumberOfPC <- as.character(100)
  NumbeerOfWindow <- as.character(0)
  #SampleName <- "A619_Re2"
  WD <- "/scratch/sb14489/3.scATAC_flo/5.Socrates/"
  setwd(WD)
  getwd()
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

# remove high conf cells form raw data set
cleanObj <- function(raw, 
                     conf){
  
  # identify shared cells
  shared <- intersect(rownames(raw$meta), rownames(conf$meta))
  raw$counts <- raw$counts[,!colnames(raw$counts) %in% shared]
  shared.sites <- intersect(rownames(raw$counts), rownames(conf$residuals))
  raw$counts <- raw$counts[shared.sites,]
  raw$counts <- raw$counts[,Matrix::colSums(raw$counts)>0]
  raw$counts <- raw$counts[Matrix::rowSums(raw$counts)>0,]
  raw$meta <- raw$meta[colnames(raw$counts),]
  
  # return
  return(raw)
  
}

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
  
  ggplot(UMAP_Table, aes(x=umap1, y=umap2,color=NumberofTn5Insertion)) +
    geom_point(size=0.5) +
    xlab("Dimension 1") +
    ylab("Dimension 2") +
    scale_x_continuous(expand=c(0.02,0)) +
    scale_y_continuous(expand=c(0.02,0)) +
    #ylim(-0.6, 0.5)+
    #xlim(-0.6, 0.5)+
    theme_bw()
  ggsave(paste0(NewFileName,"_NumberofTn5insertionCluster.pdf")
         , width=12, height=10)
}

# extract feature information for symphony integration
extractFeatureMetrics <- function(obj, 
                                  num.var=0){
  
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

# use ref cells IDF to project new data
projectTFIDF <- function(obj, 
                         old.obj, 
                         doL2=T, 
                         slotName="residuals", 
                         scale_factor=10000){
  
  # set bmat
  bmat <- obj$counts
  old.obj$counts <- old.obj$counts[intersect(rownames(old.obj$residuals), rownames(bmat)),]
  bmat <- bmat[rownames(old.obj$counts),]
  
  # hidden functions
  .safe_tfidf <- function(tf, idf,  block_size=2000e6){
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
  
  # "term frequency" method
  tf = t(t(bmat) / Matrix::colSums(bmat))
  
  # TF log scaled
  tf@x = log1p(tf@x * scale_factor)
  
  # IDF w/ "inverse document frequency smooth" method
  idf = log(1 + ncol(old.obj$counts) / Matrix::rowSums(old.obj$counts))
  
  # TF-IDF
  tf_idf_counts = .safe_tfidf(tf, idf)
  
  # do L2?
  if(doL2){
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

# new map query
#mapQuery2(new.obj$residuals[,i],          
#          new.obj$meta[i,],      
#          sym.ref,           
#          vars = "library",         
#          do_umap = TRUE,
#          doL2 = T)

#exp_query <- new.obj$residuals[,i]
#metadata_query <- new.obj$meta[i,]
#ref_obj <- sym.ref
#vars = NULL
#do_umap = TRUE
#verbose = TRUE
#sigma = 0.1
mapQuery2 <- function(exp_query, 
                      metadata_query, 
                      ref_obj, 
                      vars = NULL, 
                      verbose = TRUE,
                      do_umap = TRUE, 
                      sigma = 0.1,
                      doSTD = F,
                      doL2 = T){
  
  # hidden functions
  cosine_normalize_cpp <- function(V, dim) {
    .Call('_symphony_cosine_normalize_cpp', PACKAGE = 'symphony', V, dim)
  }
  soft_cluster <- function(Y, Z, sigma) {
    .Call('_symphony_soft_cluster', PACKAGE = 'symphony', Y, Z, sigma)
  }
  moe_correct_ref <- function(Zq, Xq, Rq, Nr, RrZtr) {
    .Call('_symphony_moe_correct_ref', PACKAGE = 'symphony', Zq, Xq, Rq, Nr, RrZtr)
  }
  
  # verbose
  if (verbose){message("Scaling and synchronizing query gene expression")}
  
  # find shared features
  idx_shared_genes <- which(ref_obj$vargenes$symbol %in% rownames(exp_query))
  shared_genes <- ref_obj$vargenes$symbol[idx_shared_genes]
  shared_genes <- intersect(shared_genes,rownames(ref_obj$loadings))

  # verbose
  if (verbose){
    message("Found ", length(shared_genes), " out of ", length(ref_obj$vargenes$symbol),
            " reference variable genes in query dataset")
  }
  
  # align matrices
  exp_query_scaled_sync <- exp_query[shared_genes,]
  #str(exp_query_scaled_sync)
  print(dim(exp_query_scaled_sync))
  rownames(exp_query_scaled_sync) <- ref_obj$vargenes$symbol
  colnames(exp_query_scaled_sync) <- colnames(exp_query)
  dim(exp_query_scaled_sync)
  # verbose
  if (verbose){message("Project query cells using reference gene loadings")}
  dim(t(ref_obj$loadings))
  dim(exp_query_scaled_sync)
  # scale TF-IDF by feature loadings
  Z_pca_query <- t(ref_obj$loadings) %*% exp_query_scaled_sync #residual
  
  # scale PCs
  if(doSTD){
    Z_pca_query <- apply(Z_pca_query, 2, function(xx){(xx-mean(xx, na.rm=T))/sd(xx, na.rm=T)})
  }
  if(doL2){
    Z_pca_query <- apply(Z_pca_query, 2, function(xx){(xx/sqrt(sum(xx^2, na.rm=T)))})
  }
  
  # verbose
  if (verbose){message("Clustering query cells to reference centroids")}
  
  # normalize cosine distances
  Z_pca_query_cos <- cosine_normalize_cpp(Z_pca_query, 2)
  R_query <- soft_cluster(ref_obj$centroids, Z_pca_query_cos, sigma)
  
  # verbose
  if (verbose){message("Correcting query batch effects")}
  if (!is.null(vars)) {
    design <- droplevels(metadata_query)[, vars] %>% as.data.frame()
    onehot <- design %>% purrr::map(function(.x) {
      if (length(unique(.x)) == 1) {
        rep(1, length(.x))
      }
      else {
        stats::model.matrix(~0 + .x)
      }
    }) %>% purrr::reduce(cbind)
    Xq = cbind(1, intercept = onehot) %>% t()
  }
  else {
    Xq = Matrix(rbind(rep(1, ncol(Z_pca_query)), rep(1, ncol(Z_pca_query))),
                sparse = TRUE)
  }
  Zq_corr = moe_correct_ref(as.matrix(Z_pca_query), as.matrix(Xq),
                            as.matrix(R_query), as.matrix(ref_obj$cache[[1]]), as.matrix(ref_obj$cache[[2]]))
  colnames(Z_pca_query) = row.names(metadata_query)
  rownames(Z_pca_query) = paste0("PC_", seq_len(nrow(Zq_corr)))
  colnames(Zq_corr) = row.names(metadata_query)
  rownames(Zq_corr) = paste0("harmony_", seq_len(nrow(Zq_corr)))
  umap_query = NULL
  if (do_umap & !is.null(ref_obj$save_uwot_path)) {
    if (verbose)
      message("UMAP")
    ref_umap_model = uwot::load_uwot(ref_obj$save_uwot_path,
                                     verbose = FALSE)
    umap_query = uwot::umap_transform(t(Zq_corr), ref_umap_model)
    colnames(umap_query) = c("UMAP1", "UMAP2")
  }
  if (verbose)
    message("All done!")
  return(list(exp = exp_query, meta_data = metadata_query,
              Z = Zq_corr, Zq_pca = Z_pca_query, R = R_query, Xq = Xq,
              umap = umap_query))
}


#############===================================================================

sym.ref <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.symphony.reference.rds")
obj_A619_merged <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.afterHarmony.processed.rds")   

str(obj_A619_merged)
head(sym.ref$umap$embedding)


obj_A619_merged$UMAP <- sym.ref$umap$embedding
rownames(obj_A619_merged$UMAP) <- rownames(obj_A619_merged$PCA)
colnames(obj_A619_merged$UMAP) <- c("umap1","umap2")
head(obj_A619_merged$UMAP)

head(obj_A619_merged$PCA)
AllRds <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/CombineAll/Combined_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100.rds"
obj_All <- readRDS(AllRds)

str(obj_All)
str(obj_A619_merged)
new.obj <- cleanObj(obj_All, obj_A619_merged)

str(obj_A619_merged)
str(new.obj)
str(sym.ref)
## IT should be same directory with Ref as "uwot_model"!!!
Path = "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100"
setwd(Path)
getwd()
#177251


new.obj <- projectTFIDF(new.obj, obj_A619_merged, doL2=T)
str(new.obj)
head(new.obj$meta)
tail(new.obj$meta) ## new.obj --> normalized mutants samples
## Check the meta
#saveRDS(new.obj, file="new.obj.rds")
#new.obj <- readRDS("new.obj.rds")

splitter <- factor(new.obj$meta$SampleName)
ids <- split(rownames(new.obj$meta), splitter)
str(ids)
its <- 0

count(rowSums(new.obj$counts)==0)
count(rowSums(obj_A619_merged$counts)==0)

## Check Things
#i <- ids$bif3
#count(Matrix::rowSums(new.obj$residuals[,i])>0)
#count(Matrix::rowSums(new.obj$residuals[,i])==0)
#dim(new.obj$meta[i,])
#dim(new.obj$residuals[,i])
#dim(new.obj$residuals[,i])
#head(new.obj$residuals[,i])
dim(obj_A619_merged$residuals)
dim(new.obj$meta[i,])
str(sym.ref)

length(sym.ref$vargenes$symbol)
dim(new.obj$residuals[,i])
i<-NULL
str(sym.ref)
split.query <- lapply(ids, function(i){
  its <<- its + 1
  
  message(" - aligning query: ",names(ids)[its])
  message("   . query contrains ",nrow(new.obj$meta[i,]), " cells")
  
  #new.obj$residuals[,i] <- new.obj$residuals[Matrix::rowSums(new.obj$residuals[,i])>0,i]
  
  mapQuery2(new.obj$residuals[,i],          
            new.obj$meta[i,],      
            sym.ref,           
            vars = c("SampleName"),         
            do_umap = TRUE,
            doL2 = T)
  
})

#str(split.query)

# merge mapped cells meta
all.meta <- lapply(split.query, function(x){
  colnames(x$umap) <- c("umap1","umap2")
  cbind(x$meta,x$umap)
})

str(all.meta)

all.meta <- do.call(rbind, all.meta)
rownames(all.meta) <- unlist(ids)
colnames(sym.ref$umap$embedding) <- c("umap1","umap2")
ref.meta <- cbind(sym.ref$meta, sym.ref$umap$embedding)
all.meta <- rbind(ref.meta, all.meta)
head(all.meta)
#head(m.obj$UMAP)

# merge mapped cells embeddings
all.embed <- lapply(split.query, function(x){
  t(x$Z)
})
str(all.embed)
#head(all.embed)

all.embed <- do.call(rbind, all.embed)
all.embed <- rbind(t(sym.ref$Z_corr), all.embed)
colnames(all.embed) <- paste0("PC_", seq(1:ncol(all.embed)))
#head(all.embed)

# make umap
all.umap <- all.meta[,c("umap1","umap2")]
rownames(all.umap) <- rownames(all.meta)
all.meta[,c("umap1","umap2")] <- NULL
head(all.umap)

#new.obj <- cleanObj(obj_All, obj_A619_merged)
#new.obj <- cleanObj(soc.obj.raw, soc.obj)
# merge counts
shared.sites <- intersect(rownames(new.obj$residuals), rownames(obj_A619_merged$residuals))

res1 <- new.obj$residuals[shared.sites,]
res2 <- obj_A619_merged$residuals[shared.sites,]
counts <- cbind(res1, res2)

# new object
m.obj <- list(PCA=as.matrix(all.embed), 
              meta=all.meta, 
              UMAP=all.umap, 
              counts=counts)

str(m.obj)
head(m.obj$UMAP)
out<- "AnnV3_RemoveBLonlyMitoChloroChIP"
saveRDS(m.obj, file=paste0(out,".ALL_CELLs_BeforeCluster.symphony_reference.rds"))
#m.obj <- readRDS(paste0(out,".ALL_CELLs_BeforeCluster.symphony_reference.rds"))

## For 
Control_PCs <- m.obj$PCA[which(m.obj$meta$SampleName=="A619"),]
Query_PCs <- m.obj$PCA[which(m.obj$meta$SampleName!="A619"),]

head(Control_PCs)
head(Query_PCs)

## Change AnnotationName
MetaData <- read.table(paste0(Path,"/Ref_AnnV3_metadata_withNAName.txt")
  ,header=TRUE)
dim(MetaData)
dim(Control_PCs)
#rownames(MetaData)
Control_PCs_Filtered  <- Control_PCs[rownames(MetaData),]
dim(Control_PCs_Filtered)

knn_pred = class::knn(Control_PCs_Filtered, Query_PCs, MetaData$Ann_v3, k = 5, prob = TRUE)

knn_prob = attributes(knn_pred)$prob
head(MetaData)
Query_meta <- m.obj$meta[which(m.obj$meta$SampleName!="A619"),]
#Query_meta$Ann_V1 <- "NA" 
Query_meta$Ann_v3 <-  knn_pred
head(Query_meta)
dim(Query_meta)

Ref_meta <- m.obj$meta[which(m.obj$meta$SampleName=="A619"),]
Ref_meta  <- Ref_meta[rownames(MetaData),]
head(Ref_meta)
head(MetaData)
MetaData_subset <- MetaData[c("Ann_v3")]
head(MetaData_subset)

Ref_Meta <- merge(Ref_meta, MetaData_subset, by=0)
#?merge
head(Ref_Meta)
dim(Ref_Meta)
colnames(Ref_Meta)
colnames(Query_meta)
Ref_Meta <- Ref_Meta[-c(1)]

Clusters_Predicted <- rbind(Ref_Meta,Query_meta)
head(Clusters_Predicted)
ClustersList <- list()
ClustersList[["Clusters"]] <- Clusters_Predicted
ClustersList[["Clusters_Mutant"]] <- Query_meta
ClustersList[["Clusters_Control"]] <- Ref_Meta

head(ClustersList[["Clusters"]])
colorr <- c("#4F96C4","#84f5d9","#DE9A89","#FDA33F","#060878","#d62744","#62a888",
            "#876b58","#800000", "#800075","#e8cf4f","#f7366d","#0bd43d",
            "#deadce","#adafde","#5703ff")
head(Ref_Meta)

ggplot(Clusters_Predicted, aes(x=umap1, y=umap2, color=factor(Ann_v3))) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave("AnnotationV3_KnnPredictionCluster_Annotation_AllSamples.pdf", width=13, height=10)	

ggplot(Query_meta, aes(x=umap1, y=umap2, color=factor(Ann_v3))) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave("AnnotationV3_KnnPredictionCluster_Annotation_Mutant.pdf", width=13, height=10)	


ggplot(Ref_Meta, aes(x=umap1, y=umap2, color=factor(Ann_v3))) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave("AnnotationV3_KnnPredictionCluster_Annotation_Ref.pdf", width=13, height=10)	

head(Query_meta)
Bif3_Meta <- Query_meta[which(Query_meta$SampleName=="bif3"),]

head(Bif3_Meta)

ggplot(Bif3_Meta, aes(x=umap1, y=umap2, color=factor(Ann_v3))) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave("AnnotationV3_KnnPredictionCluster_Annotation_Bif3.pdf", width=13, height=10)	


# write data
write.table(Clusters_Predicted, "AnnotationV3.All_CELLs.metadata.txt", quote=F, row.names=T, col.names=T, sep="\t")
write.table(Ref_Meta, "AnnotationV3.A619.metadata.txt", quote=F, row.names=T, col.names=T, sep="\t")
write.table(Bif3_Meta, "AnnotationV3.bif3.metadata.txt", quote=F, row.names=T, col.names=T, sep="\t")

#write.table(rd, file=paste0(out, ".REF_CELLs.reduced_dimensions.txt"), quote=F, row.names=T, col.names=T, sep="\t")
levels(factor(Clusters_Predicted$Ann_v3))

