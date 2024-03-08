Meta_Bif <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP_k50_res0.9.AfterHarmony.metadata.txt")
head(Meta_Bif)
dim(Meta_Bif)

Meta_Bif_Ref <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/AnnotationV3.bif3.metadata.txt")
head(Meta_Bif_Ref)
dim(Meta_Bif_Ref)


Bif_Labelling <- Meta_Bif[rownames(Meta_Bif)%in%rownames(Meta_Bif_Ref),]
Bif_UMP <- Meta_Bif_Ref[rownames(Meta_Bif_Ref)%in%rownames(Meta_Bif),]
dim(Bif_Labelling)
dim(Bif_UMP)
head(Bif_UMP)
head(Bif_Labelling)
Bif_UMP$ClusterFromOriginal <- Bif_Labelling[rownames(Bif_UMP),]$LouvainClusters

colorr <- c("#4F96C4","#84f5d9","#DE9A89","#FDA33F","#060878","#d62744","#62a888",
            "#876b58","#800000", "#800075","#e8cf4f","#f7366d","#0bd43d",
            "#deadce","#adafde","#5703ff")
length(colorr)
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/")
ggplot(Bif_UMP, aes(x=umap1, y=umap2, color=factor(ClusterFromOriginal))) +
  geom_point(size=0.02) +
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave("Bif3_LabellingbyOwn_butRefProjectedUMAP.pdf", width=13, height=10)
