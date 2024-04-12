About the dataset…
================
Margaux Lefebvre
2024-04-12

# The sources

The samples are drawn from several literature sources:

- [Daron *et al.* (2021)](https://doi.org/10.1126/sciadv.abc3713)
  (n=499)
- [Benavente *et al* (2021)](https://doi.org/10.1038/s41467-021-23422-3)
  (n=125)
- [Malaria Gen project *P. vivax* Genome
  Variation](https://doi.org/10.12688/wellcomeopenres.17795.1) (n=295)
- [van Dorp *et al* (2020)](https://doi.org/10.1093/molbev/msz264) for
  the Ebro sample only (n=1)

# Ebro: take care of ancient DNA samples

The study accession is PRJEB30878 in ENA. All the fastq are single ends.

|          Name          |    SRA     |
|:----------------------:|:----------:|
| CA-plasmodium-capture  | ERR3650065 |
| POS-plasmodium-capture | ERR3650068 |
| CM-plasmodium-capture  | ERR3650072 |
|         lane-8         | ERR3651363 |

In [van Dorp *et al* (2020)](https://doi.org/10.1093/molbev/msz264),
Ebro is a composite of the 4 samples.

## Filter the fastq and merge them

Version: cutadapt v2.8.

``` bash
echo "--> Start processing: $downId" #downId is the name of the sample

# remove the 3' adapters (-a), remove reads with low quality score (-q 30) and remove reads shorther than 30 bp (-m 30)
cutadapt -a AGATCGGAAGAG -j 0 -q 30 -m 30 $downId.fastq.gz -o $downId.cut.fastq.gz

# Create a composite of Ebro, by merging all the reads together
zcat CA.cut.fastq.gz CM.cut.fastq.gz POS.cut.fastq.gz lane-8.cut.fastq.gz > $downId.fastq
```

## Mapping

As I will deal with *Plasmodium falciparum* co-infection, I also create
a bam by aligning on Pf_3D7 reference genome.

Versions: bwa-mem v0.7.17, samtools v1.9.

``` bash
echo "--> Start processing: $downId" #downId is the name of the sample

# Mapping on P. vivax genome
bwa mem -t 4 P.vivax-PvP01_reference_genome.fasta $downId.cut.fastq.gz | samtools view -F 4 -b - | samtools sort - -o $downId.PV.mem.sort.bam
samtools index $downId.PV.mem.sort.bam 
# Mapping on P. falciparum genome
bwa mem -t 4 PlasmoDB-57_Pfalciparum3D7_Genome.fasta $downId.cut.fastq.gz | samtools view -F 4 -b - | samtools sort - -o $downId.PF.mem.sort.bam 
samtools index $downId.PF.mem.sort.bam
```

## Dealing with *P. falciparum* co-infection

Version: samtools v1.9, Picard tools v2.5.0.

### Extract the reads that mapped on both reference genomes

``` bash
echo "--> Start processing: $downId" #downId is the name of the sample

samtools view $downId.PV.mem.sort.bam | cut -f1 > $downId.PV_reads_mem.txt
samtools view $downId.PF.mem.sort.bam | cut -f1 > $downId.PF_reads_mem.txt

# Extract the reads that mapped on both genome
grep -Fxf $downId.PF_reads_mem.txt $downId.PV_reads_mem.txt > $downId.common_reads_mem.txt 

# For both bam, keep only the reads in common
java -Xss5m -jar /usr/local/picard-tools-2.5.0/picard.jar FilterSamReads I=$downId.PV.mem.sort.bam O=$downId.common.PV.mem.bam READ_LIST_FILE=$downId.common_reads_mem.txt FILTER=includeReadList
java -Xss5m -jar /usr/local/picard-tools-2.5.0/picard.jar FilterSamReads I=$downId.PF.mem.sort.bam O=$downId.common.PF.mem.bam READ_LIST_FILE=$downId.common_reads_mem.txt FILTER=includeReadList
```

### Extract and analyse the edit distance

The NM tag is the edit distance and it’s in the 12th column with bwa
mem.

``` bash
echo "--> Start processing: $downId" #downId is the name of the sample

samtools view $downId.common.PV.mem.bam | awk '{print $1, $12}' | sed 's/NM:i://g' > $downId.common.PV_mem.txt
samtools view $downId.common.PF.mem.bam | awk '{print $1, $12}' | sed 's/NM:i://g' > $downId.common.PF_mem.txt
```

Analyse the edit distance

``` r
library(tidyverse)
samples<-c("lane-8", "POS", "CM", "CA")
for (samp in samples){
  data_mem_PF<-read_table(paste0("./Data/",samp,".common.PF_mem.txt"), col_names =c("Reads","edit_dist_PF"))
  data_mem_PV<-read_table(paste0("./Data/",samp,".common.PV_mem.txt"), col_names = c("Reads","edit_dist_PV"))
  data_mem_tot<-inner_join(data_mem_PV, data_mem_PF)
  data_temp<-data_mem_tot[data_mem_tot$edit_dist_PF != data_mem_tot$edit_dist_PV, ]

# Plot the density for each sample
  plot_dens<-ggplot(data_temp, aes(x=edit_dist_PV, y=edit_dist_PF) ) +
 stat_density_2d(aes(fill = ..level..), geom = "polygon", bins=10) +
    scale_fill_gradient(low = "#4ECDC4", high = "#556270", name=paste0("Density of shared reads \nwith non equal distance \n(n=",nrow(data_temp),")"))+
    labs(y = "Edit distance to P. falciparum", x = "Edit distance to P. vivax")+ggtitle(paste0(samp))+
  theme_bw()
  data_temp<-data_mem_tot[data_mem_tot$edit_dist_PF <= data_mem_tot$edit_dist_PV, ]
  
# Write list of the reads that align equally or better with P. falciparum than P. vivax
write_tsv(as.data.frame(data_temp$Reads), file=paste0("./Data/",samp,".remove_reads_mem.txt"), col_names = F)
    name_plot<-paste0("plot_mem_",samp)
  assign(name_plot,plot_dens)
}

layout.matrix <- matrix(c(1,2,3,4), nrow = 2, ncol = 2)
library("gridExtra")
grid.arrange(`plot_mem_lane-8`, plot_mem_CM, plot_mem_POS, plot_mem_CA, layout_matrix=layout.matrix)
```

With the R code chunk just above, I made files called
`$downId.remove_reads_mem.txt` and it’s the list of the reads that align
equally or better with *P. falciparum* than *P. vivax*. So as they did
in [Van Dorp *et al.* (2020)](https://doi.org/10.1093/molbev/msz264), I
remove them.

### Remove the reads that map better on P. falciparum

``` bash
echo "--> Start processing: $downId" #downId is the name of the sample

java -Xss5m -jar /usr/local/picard-tools-2.5.0/picard.jar FilterSamReads I=$downId.PV.mem.sort.bam O=$downId.PV.mono.mem.bam READ_LIST_FILE=$downId.remove_reads_mem.txt FILTER=excludeReadList
```

# Produce unfiltered VCF with all samples

## Mapping

Version: cutadapt v1.18, bwa-mem v0.7.17, samtools v1.9.

``` bash
echo "--> Start processing: $downId" #downId is the name of the sample

#Rename the files in order to remain consistent with the rest of the script
mv ${downId}_1.fastq.gz ${downId}.R1.fastq.gz
mv ${downId}_2.fastq.gz ${downId}.R2.fastq.gz

  zcat ${downId}.R1.fastq.gz | sed -E 's/^((@|\+)'$downId'\.[^.]+)\.(1|2)/\1/' | bgzip > $downId.raw.1.fastq.gz
  zcat ${downId}.R2.fastq.gz | sed -E 's/^((@|\+)'$downId'\.[^.]+)\.(1|2)/\1/' | bgzip > $downId.raw.2.fastq.gz

# Remove adapters and preprocessed to eliminate low-quality reads. Reads shorter than 50 bp containing “N” were discarded.
  cutadapt -a AGATCGGAAGAGCACACGTCTGAA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 30 -m 25 --max-n 0 -o $downId.R1.fastq.gz -p $downId.R2.fastq.gz $downId.raw.1.fastq.gz $downId.raw.2.fastq.gz

#Sequenced reads were aligned to the P. vivax reference genome PVP01
  bwa mem -t 1 P.vivax-PvP01_reference_genome.fasta $downId.R1.fastq.gz $downId.R2.fastq.gz | samtools view -F 4 -b - | samtools sort - -o $downId.mapPV.sort.bam

samtools index $downId.mapPV.sort.bam
```

## Calling

Version: GATK v3.8.0, samtools v1.9, Picard tools v2.5.0.

``` bash
# Create Sequence Dictionary
java -jar /usr/local/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=P.vivax-PvP01_reference_genome.fasta O=P.vivax-PvP01_reference_genome.dict

# Mark the duplicated reads
java -Xss5m -jar /usr/local/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=$downId.mapPV.sort.bam OUTPUT=$downId.map.sort.dedup.bam METRICS_FILE=metrics.txt
samtools index $downId.map.sort.dedup.bam

#Add or Replace Read Groups
java -Xss5m -jar /usr/local/picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups I=$downId.map.sort.dedup.bam O=$downId.map.sort.dedup.rg.bam LB=LIB-$downId PL=ILLUMINA PU=H0164ALXX140820:2:1101 SM=$downId

rm $downId.map.sort.dedup.bam $downId.map.sort.dedup.bam.bai
samtools index $downId.map.sort.dedup.rg.bam

# SplitNCigarReads
java -Xss5m -jar ~/softs/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T SplitNCigarReads -R P.vivax-PvP01_reference_genome.fasta -I $downId.map.sort.dedup.rg.bam -o $downId.map.sort.dedup.rg.rq.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

rm $downId.map.sort.dedup.rg.bam $downId.map.sort.dedup.rg.bam.bai
samtools index $downId.map.sort.dedup.rg.rq.bam

# Local realignment around indels
java -Xss5m -jar ~/softs/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T RealignerTargetCreator -R P.vivax-PvP01_reference_genome.fasta -I $downId.map.sort.dedup.rg.rq.bam -o $downId.realignertargetcreator.intervals

java -Xmx8G -Djava -jar ~/softs/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T IndelRealigner -R P.vivax-PvP01_reference_genome.fasta -targetIntervals $downId.realignertargetcreator.intervals -I $downId.map.sort.dedup.rg.rq.bam -o $downId.map.sort.dedup.indelrealigner.bam
mv $downId.map.sort.dedup.indelrealigner.bam $downId.bwa.gatk.sort.bam

samtools index $downId.bwa.gatk.sort.bam

# Calling with HaplotypeCaller for ploidy 1
java -Xss5m -jar ~/softs/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T HaplotypeCaller -R P.vivax-PvP01_reference_genome.fasta -I $downId.bwa.gatk.sort.bam --genotyping_mode DISCOVERY -stand_call_conf 10 -o $downId.ploidy2.raw_variants.snp.indel.g.vcf -ERC GVCF --sample_ploidy 2
```

## Combining all the samples and keep only the nuclear genome

Version: GATK v3.8.0, bcftools v1.10.2.

``` bash
#Merge all the samples
java -Xmx8G -Djava -jar ~/softs/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T CombineGVCFs -R P.vivax-PvP01_reference_genome.fasta --variant vcfs.ploidy2.list -o vivax_temp.ploidy2.vcf
bgzip vivax_temp.ploidy2.vcf
tabix vivax_temp.ploidy2.vcf.gz

#Keep only the variants
java -Xmx8G -Djava -jar ~/softs/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T GenotypeGVCFs -R P.vivax-PvP01_reference_genome.fasta -V vivax_temp.ploidy2.vcf.gz -o  vivax_temp.raw_variants.ploidy2.vcf
bgzip vivax_temp.raw_variants.ploidy2.vcf
tabix  vivax_temp.raw_variants.ploidy2.vcf.gz

#Keep nuclear core_genome
bcftools view  vivax_temp.raw_variants.ploidy2.vcf.gz -R vivax_core_genome.bcf.txt -O z -o vivax.core.snps.ploidy2.vcf.gz
tabix vivax.core.snps.ploidy2.vcf.gz
```
