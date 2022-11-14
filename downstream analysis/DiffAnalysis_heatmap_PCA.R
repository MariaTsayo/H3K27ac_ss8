##------------------------------------------------------------------------------##
##  Diff analysis - heatmap - PCA                                               ##
##------------------------------------------------------------------------------##
library(tidyverse)
library(limma)
library(ggplot2)
library(ComplexHeatmap)
library(gplots)

acet <- read.csv( "H3K27ac_freq2_FDR_consensus_2022.txt", 
                  header = T, sep = "\t", check.names = FALSE)
annotation<-  #meta data
combat_edata<- read.csv( "H3K27ac_FDR_consensus_norm.Combat_2022.txt", 
                    header = T, sep = "\t", check.names = FALSE) #batch effect corrected table

######################################################

#diff analysis
design_data_ALL<- annotation %>% filter(Status == "CLL-ss8" | Status== "U-CLL") 
dfexp <-  as.matrix(combat_edata)
row.names(dfexp)<- acet$ID


list<- as.character(design_data_ALL$ID)
dfexp<- combat_edata[,list]
row.names(dfexp)<- acet$ID
pheno <- factor(design_data_ALL$Status)
phenoMat <- model.matrix(~pheno)
colnames(phenoMat) <- sub("^pheno","",colnames(phenoMat))
fit <- lmFit(object = dfexp,design = phenoMat)
gc()
set.seed(6)
fit <- eBayes(fit)
degCLLs <- topTable(fit,number =nrow(dfexp),adjust.method = "fdr",sort.by = "p")

# filtering
sign.table<- as.data.frame(degCLLs)
sign.table<- sign.table[sign.table$adj.P.Val<= 0.0001, ]

sign.table$ID<- row.names(sign.table)
list<- as.character(annotation$ID)
full.data<- combat_edata[,list]
full.data<- as.data.frame(full.data)
full.data$ID<- acet$ID
sign.table.full<- merge(sign.table,full.data,by.x="ID",by.y="ID")

#heatmap
z.mat <- t(scale(t(sign.table.full[,-c(1:7)]), center=TRUE, scale=TRUE))
design_data_ALL<-annotation
design_data_ALL$Status <- design_data_ALL$Status  %>% str_replace_all("additional_U-CLL", "U-CLL")
column_ha3 = HeatmapAnnotation(groups = design_data_ALL$Status, #No_of_reads = anno_points(design_data_ALL$reads),
                               IGHV= design_data_ALL$IGHV, tris12= design_data_ALL$tris12,
                               col = list(groups = c("GCBC" = "orange", "M-CLL" ="chartreuse","MBC"= "steelblue2", 
                                                     "NBC-B"="yellow","NBC-T"="red","PC"="purple","CLL-ss1"= "black", "CLL-ss16"="brown", 
                                                     "CLL-ss2"="beige", "CLL-ss4"="darkcyan","CLL-ss8"="pink", 
                                                     "Tris12"="mediumvioletred", "U-CLL"="grey", 
                                                      "Tris12-U-CLL"="mediumvioletred",
                                                     "IgG-CLL"= "coral",
                                                     "IgG-M-CLL"= "chocolate3", "IgG-U-CLL"= "blue"),
                                          IGHV = c("M-CLL" = "dodgerblue4", "U-CLL" = "gold", "Undetermined"="gray90"),
                                          tris12=c("tris12"='gray3', "non-tris12"="gray51")
                               ),
                               na_col="gray90"
)



cl = kmeans(z.mat, centers = 5)$cluster
Heatmap(as.matrix(z.mat), clustering_distance_columns = "euclidean",row_split = cl,
        clustering_method_columns = "complete", 
        top_annotation = column_ha3, left_annotation = 
          rowAnnotation(foo = anno_block(gp = gpar(fill = 2:8),
                                         labels_gp = gpar(col = "white", fontsize = 10))))


# merge the generated cluster with the significant regions
cut.info<- as.data.frame(cl)
cut.info$rowID<- row.names(cut.info)
sign.table.full$rowID<- row.names(sign.table.full)
sign.table.full.2<- merge(cut.info,sign.table.full, by="rowID")
acet.bed<- acet[,c(1:4)]
sign.table.full.3<- merge( acet.bed, sign.table.full.2, by="ID")

#PCA
acet<- combat_edata
data<-as.matrix(acet)
data1<-as.data.frame(data)
data1$ID<- paste("X")
data1$n<- as.character(row.names(data1))
data1$ID<- paste0(data1$ID,data1$n)
data<- as.data.frame(data)
data<-data.frame(t(data))
pcamatrix<-data[,as.vector(data1$ID)]

k=length(pcamatrix)
pcamatrix$Treatment<-annotation$Status
pcamatrix<-subset(pcamatrix,select=c(k+1,1:k))


pca<-prcomp(pcamatrix[,-1],scale. = F)
pca_data <- as.data.frame(pca$x[,c(1:10)])
pca_data$status <- annotation$Status

pca_data %>%
  arrange(PC1) %>%
  mutate(status = factor(status, levels=c("NBC-B", "NBC-T", "GCBC", "MBC", "PC",
                                      "CLL-ss1" , "CLL-ss16", 
                                      "CLL-ss2", "CLL-ss4","CLL-ss8", 
                                      "M-CLL","U-CLL", "IgG-CLL", "IgG-M-CLL", 
                                      "IgG-U-CLL","Tris12-U-CLL"))) %>%
  ggplot(aes(x=PC1, y=PC2, fill=status)) + 
  geom_point(shape=21, color="gray61", size=6) + 
  scale_fill_manual(values=c("yellow", "red", "orange", "steelblue2", "purple",
                             "black", "brown", "beige", "darkcyan", "pink",
                             "chartreuse", "grey", "coral","chocolate3",
                             "blue","mediumvioletred"
                             
  )) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  labs(x ="PC1(30%)", y = "PC2 (13%)")

