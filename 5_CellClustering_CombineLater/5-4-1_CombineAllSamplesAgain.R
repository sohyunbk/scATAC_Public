library("here")
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library("optparse")
library(rlang)
library(ggplot2)
library("RColorBrewer")

Ex <- function(){
  Dir <- as.character("CombineAll") #NewDir name
  Prefix <- "Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100"
  WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.CellClustering/Organelle5Per_CombineLater/"
}

args <- commandArgs(T)
Dir <- as.character(args[1])
minimumtn5counts <- as.character(args[2])
Minimum_PeakNumberPerCell <- as.character(args[3])
MinT<- as.character(args[4])
MaxT <- as.character(args[5])
SVDorNMF <- as.character(args[6])
NumberOfPC <- as.character(args[7])


# functions --------------------------------------------------------------

# remove high conf cells from raw data set
LoadPreviousRds <- function(Name,Prefix,libName,NameSample){
  File <- paste0(WD,"/",Name,"/",Name,"_",Prefix)
  Pastobj <- readRDS(paste(File,".rds",sep=""))
  Pastobj$meta$SampleName <- NameSample
  Pastobj$meta$library <- libName
  newobj <- list()
  newobj$meta <- Pastobj$meta 
  newobj$counts <- Pastobj$counts 
  return(newobj)
  }

NotIntersectNumber <- function(Avector,Bvector){
  print(length(union(Avector,Bvector))-length(intersect(Avector,Bvector)))}

#############===================================================================
## Combine All together in the begening!!

obj_A619_Re1 <- LoadPreviousRds("1_A619","Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100","A619_Re1","A619")
obj_A619_Re2 <- LoadPreviousRds("1_A619_2","Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100","A619_Re2","A619")
## tfidf was performed by eash sample. I will combine the replicates and then tfidf again
## Extract only meta/counts for merging! ##
#obj_A619_Re1$meta$lib_ID <- "1"
str(obj_A619_Re1)
str(obj_A619_Re2)
##Check if there are some rows in zero sum
#rownames(obj_A619_Re2$counts)[1:4]
#Temp <- obj_A619_Re2$counts[Matrix::rowSums(obj_A619_Re2$counts) > 0,]


obj_rel2_Re1 <- LoadPreviousRds("2_rel2","Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100","rel2_Re1","rel2")
obj_rel2_Re2 <- LoadPreviousRds("2_rel2_2","Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100","rel2_Re2","rel2")
str(obj_rel2_Re1)
str(obj_rel2_Re2)

obj_bif3_Re1 <- LoadPreviousRds("3_bif3","Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100","bif3_Re1","bif3")
obj_bif3_Re2 <- LoadPreviousRds("3_bif3_2","Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100","bif3_Re2","bif3")
str(obj_bif3_Re1)
str(obj_bif3_Re2)

obj_relk1_Re1 <- LoadPreviousRds("4_relk1","Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100","relk1_Re1","relk1")
obj_relk1_Re2 <- LoadPreviousRds("4_relk1_2","Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100","relk1_Re2","relk1")
str(obj_relk1_Re1)
str(obj_relk1_Re2)
#rownames(obj_relk1_Re1$counts)


NotIntersectNumber(rownames(obj_A619_Re1$counts),rownames(obj_A619_Re2$counts))
NotIntersectNumber(rownames(obj_rel2_Re1$counts),rownames(obj_rel2_Re2$counts))
NotIntersectNumber(rownames(obj_bif3_Re1$counts),rownames(obj_bif3_Re2$counts))
NotIntersectNumber(rownames(obj_relk1_Re1$counts),rownames(obj_relk1_Re2$counts))

SharedFeatures <- Reduce(intersect, list(rownames(obj_A619_Re1$counts),rownames(obj_A619_Re2$counts),
                                        rownames(obj_rel2_Re1$counts),rownames(obj_rel2_Re2$counts),
                                      rownames(obj_bif3_Re1$counts),rownames(obj_bif3_Re2$counts),
                                      rownames(obj_relk1_Re1$counts),rownames(obj_relk1_Re2$counts)))

length(SharedFeatures)
print(SharedFeatures[1:5])
#print(SharedFeatures)

files <- list(obj_A619_Re1,obj_A619_Re2, obj_rel2_Re1, obj_rel2_Re2,
              obj_bif3_Re1, obj_bif3_Re2,obj_relk1_Re1,obj_relk1_Re2 )
names(files) <- c("A619_Re1", "A619_Re2","rel2_Re1", "rel2_Re2",
                  "bi3_Re1","bi3_Re2","relk1_Re1","relk1_Re2")

merged.obj <- mergeSocratesRDS(obj.list=files)

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
str(merged.obj)
#####################
str(merged.obj)
head(merged.obj$counts[,c(1:10)])
merged.obj$counts <- merged.obj$counts[rownames(merged.obj$counts) %in% SharedFeatures,]
str(merged.obj)

Name <- paste0("Combined_",Prefix)
saveRDS(merged.obj, file=paste0(Name,".rds"))

