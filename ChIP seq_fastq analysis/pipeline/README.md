### trim_galore
```
path/trim_galore --path_to_cutadapt path/cutadapt-1.8.1/bin/cutadapt sample_R1.fastq > sample_trimming.output
```
### mapping
```
bwa aln -q 15 -t 16 path/genome.fa sample_R1_trimmed.fq > sample_intermediate.sai ; bwa samse path/genome.fa sample_intermediate.sai sample_R1_trimmed.fq | samtools view -bS - > sample.bam && samtools sort -@ 16 -o sample_sorted.bam sample.bam 
```
### mark duplicates
```
java -jar /work/fotis/varCall/picard.jar MarkDuplicates  I=sample_sorted.bam O=sample_duplicates.bam M=sample_dup_metrics.txt VALIDATION_STRINGENCY=LENIENT 
```
### remove duplicates
```
samtools view -@ 16 -h -F 4 -q 5 -b sample_duplicates.bam > sample_sorted.intermediate.bam 
samtools view  -@ 16 -h -F 1024 -b sample_sorted.intermediate.bam > sample_sorted.uniq.bam 
```
### mapping report
```
samtools flagstat sample_sorted.bam > sample_mapping_report.out
```
### peak calling
with Input file
```
macs2 callpeak -t sample_sorted.uniq.bam  --name=sample --gsize hs -c Input_sorted.uniq.bam --nomodel --format=BAM
```
without Inpute file
```
macs2 callpeak -t sample_sorted.uniq.bam  --name=sample --gsize hs --nomodel --format=BAM
```
### create a total table with all peaks
```
cat *peaks.narrowPeak | awk '{print $1 "," $2 ","  $3}' | csvtk csv2tab >> consensus.bed
```
### sort
```
sort -k1,1 -k2,2n consensus.bed > consensus_sorted.bed
```
### merge peaks
```
bedtools merge -i consensus_sorted.bed > consensus_merged.bed
```
### for each sample - count on the consensus
```
bedtools coverage -a consensus_merged.bed -b sample_sorted.uniq.bam -counts > consensus_sample_counts.bed  
```
### merge the files with the headers
```
(for f in *_counts.bed ; do printf "1\t2\t3\t"; printf "$f\t"; done; printf "\n") > total_consensus_count.txt
paste *_counts.bed >> total_consensus_count.txt
```
