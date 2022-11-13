
### required R packages
```
library("tidyverse")
library("data.table")
require("tximport")
```

### transcript to genes object
```
tx2gene = read.delim("ensembl_genes_transcripts_v100.txt.gz",
                     header = T, stringsAsFactors = F)[,c(4,2)]
colnames(tx2gene) = c("TXNAME","GENEID")
```

### import RNA-seq output from kallisto 
```
kPath = path_to_files
files = sprintf("%s/%s", kPath, list.files(kPath, pattern = ".*\\.tsv"))

rnaK = tximport(files[1], type = "kallisto", tx2gene = tx2gene,
                ignoreAfterBar = TRUE, countsFromAbundance = "scaledTPM")
data<-as.data.frame(rnaK$counts)
names(data)<- files[1]

for(i in 2:length(files))
{
  rnaK = tximport(files[i], type = "kallisto", tx2gene = tx2gene,
                  ignoreAfterBar = TRUE, countsFromAbundance = "scaledTPM")
  test2<-as.data.frame(rnaK$counts)
  names(test2)<- files2[i]
  data<- cbind(data, test2)
  
}
```

### log-transform to make numbers on scale (+1 to avoid zeroes)
```
pseudoCount = log2(data + 1) 
```

### reshape the matrix
```
df = melt(pseudoCount, variable.name = "Samples", value.name = "count")
```
