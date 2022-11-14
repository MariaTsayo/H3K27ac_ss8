##------------------------------------------------------------------------------##
##  Random resampling                                                           ##
##------------------------------------------------------------------------------##
library(limma)
library(tidyverse)  
library(boot)  
library(limma)

acet <- read.csv( "H3K27ac_freq2_FDR_consensus_2022.txt", 
                  header = T, sep = "\t", check.names = FALSE)

annotation<- #metadata

combat_edata<- read.csv("H3K27ac_FDR_consensus_norm.Combat_2022.txt", 
                         header = T, sep = "\t", check.names = FALSE)


# data manipulation
ss8<- annotation %>% filter(Status == "CLL-ss8" ) 
UCLL<- annotation %>% filter(Status == "U-CLL" ) 
list8<- as.character(ss8$ID)
listU = as.character(UCLL$ID)
dfexp8<- combat_edata[,list8]
dfexpU =combat_edata[,listU]
row.names(dfexp8)<- acet$ID
row.names(dfexpU)<- acet$ID

auto8 = colnames(dfexp8)
autoU = colnames(dfexpU)

set.seed(1)
sample8 <- sample(auto8, size=sample(5:6, 1), replace =F)
sampleU <- sample(autoU, size=sample(5:8, 1), replace =F)

df<- cbind(dfexp8[,sample8], dfexpU[,sampleU])
design_data_ALL = annotation %>% filter(ID %in% colnames(df))

pheno <- factor(design_data_ALL$Status)
phenoMat <- model.matrix(~pheno)
colnames(phenoMat) <- sub("^pheno","",colnames(phenoMat))
phenoMat;dim(phenoMat)


fit <- lmFit(object = df,design = phenoMat)
gc()
set.seed(6)
fit <- eBayes(fit)
degCLLs <- topTable(fit,number =nrow(df),adjust.method = "fdr",sort.by = "p")

sign.table<- as.data.frame(degCLLs)
sign.table<- sign.table[sign.table$adj.P.Val<= 0.0001, ]
sign.table$ID<- row.names(sign.table)
sign.table$round = paste(1)
total = sign.table

# r sample multiple times without replacement
for(i in 1:99) {
  set.seed(i)

  sample8 <- sample(auto8, size=sample(5:6, 1), replace =F)
  sampleU <- sample(autoU, size=sample(5:8, 1), replace =F)
  
  df<- cbind(dfexp8[,sample8], dfexpU[,sampleU])
  design_data_ALL = annotation %>% filter(ID %in% colnames(df))
  
  pheno <- factor(design_data_ALL$Status)
  phenoMat <- model.matrix(~pheno)
  colnames(phenoMat) <- sub("^pheno","",colnames(phenoMat))
  phenoMat;dim(phenoMat)
  

  
  colnames(df)
  fit <- lmFit(object = df,design = phenoMat)
  gc()
  set.seed(6)
  fit <- eBayes(fit)
  
  gc()
  degCLLs <- topTable(fit,number =nrow(df),adjust.method = "fdr",sort.by = "p")
  head(degCLLs)
  
  sign.table<- as.data.frame(degCLLs)
  sign.table<- sign.table[sign.table$adj.P.Val<= 0.0001, ]
  sign.table$ID<- row.names(sign.table)
  sign.table$round = paste(i)

  total = rbind(total, sign.table)
}


sum_tiss<-total %>% 
  dplyr::group_by(ID) %>% 
  dplyr::summarise(freq = n()) 

hist(sum_tiss$freq)
ggplot(sum_tiss, aes(x=freq)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 

filt=sum_tiss %>% filter(freq >=20)

#####heatmap
list<- as.character(annotation$ID)
full.data<- combat_edata[,list]
head(full.data)
full.data<- as.data.frame(full.data)
full.data$ID<- acet$ID


sign.table.full<- merge(filt,full.data,by.x="ID",by.y="ID")
colnames(sign.table.full)
sign.table.full= sign.table.full[,-2]

z.mat <- t(scale(t(sign.table.full[,-c(1)]), center=TRUE, scale=TRUE))
head(z.mat)


library(ComplexHeatmap)
library(gplots)
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




Heatmap(as.matrix(z.mat), clustering_distance_columns = "euclidean",row_split = cl,
        clustering_method_columns = "complete", 
        top_annotation = column_ha3, left_annotation = 
          rowAnnotation(foo = anno_block(gp = gpar(fill = 2:8),
                                         labels_gp = gpar(col = "white", fontsize = 10))))




#merge the cluster with the significant regions
cut.info<- as.data.frame(cl)
cut.info$rowID<- row.names(cut.info)
head(cut.info)
sign.table.full$rowID<- row.names(sign.table.full)

sign.table.full.2<- merge(cut.info,sign.table.full, by="rowID")
acet.bed<- acet[,c(1:4)]
colnames(meta_2)

sign.table.full.3<- merge( acet.bed, sign.table.full.2, by="ID")
head(sign.table.full.3)
write.table(sign.table.full.3,file="signature_RESAMPLING_ss8VSucll_2022.txt", col.names=T, row.names=F, quote=F, sep="\t")
