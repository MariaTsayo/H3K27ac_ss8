##------------------------------------------------------------------------------##
##  ComBat for removing batch effects in H3K27ac                                ##
##------------------------------------------------------------------------------##


library(sva)
library(DESeq2)
library(stringi)


acet <- read.csv( "H3K27ac_freq2_FDR_consensus_2022.txt", 
                  header = T, sep = "\t", check.names = FALSE)

annotation<-  #load metadata


###########################################################
#deseq2 normalization
design_data_ALL <- annotation
dds <- NULL 
countdata_keeps_ALL<- acet
dds <- DESeqDataSetFromMatrix(countData=countdata_keeps_ALL,
                              colData=design_data_ALL, 
                              design=~Status) # Status is CLL-ss8, CLL-ss1, UCLL, MCLL etc

dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
vsd <- vst(dds)
vvt <-assay(vsd)
res <- results(dds)

#ComBat
data<-as.matrix(vvt)
colnames(annotation)
ctype<- annotation
row.names(ctype)<-ctype$ID
modcombat <- model.matrix(~Status, data = ctype)
batch <- annotation$Read_bin

edata <- data
colnames(edata)
edata<- as.matrix(data)

combat_edata <- ComBat(dat = edata, batch = batch, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)


write.table(combat_edata, file = "H3K27ac_FDR_consensus_norm.Combat_2022.txt", quote = F, sep = "\t", row.names = F)

