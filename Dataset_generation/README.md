About the dataset…
================
Margaux Lefebvre
2024-04-24

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

### Remove the reads that map better on *P. falciparum*

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

# Filtering of the dataset

## Quality and depth filtering

Version: bcftools v1.10.2, vcftools v0.1.16.

First, I will focus on the VCF with ploidy = 2. With this file, I can
calculate the co-infection index, essential to filter the dataset
properly.

``` bash
# Filter out a little to have only the information we need

### Keep only SNPs, regardless the number of alleles, the MAF or even the quality
vcftools --gzvcf vivax.core.snps.ploidy2.vcf.gz \
--remove-indels --non-ref-ac-any 1 \
--recode --stdout | bgzip -c > Pvivax_total_snp.ploidy2.vcf.gz

# See the info
VCF=Pvivax_total_snp.ploidy2.vcf.gz
OUT=Pvivax_total.ploidy2

vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2 &
vcftools --gzvcf $VCF --depth --out $OUT &
vcftools --gzvcf $VCF --site-mean-depth --out $OUT &
vcftools --gzvcf $VCF --site-quality --out $OUT &
vcftools --gzvcf $VCF --missing-indv --out $OUT &
vcftools --gzvcf $VCF --missing-site --out $OUT &
```

Filtering info:

- The minimum variant quality (Phred score) is 30; it’s ok.
- For the variant mean depth, we will set the minimum at 10X. The
  maximum will be set at 106X (= mean depth + twice the
  standard-deviation).
- For the individual mean depth, we will set the minimum at 10X and the
  maximum at 164X (= mean depth + twice the standard-deviation).
- We remove all the individuals with more than 50% of missing data: it
  removes 68 sample.
- Remove all the SNPs with more than 20% of missing data.
- To avoid sequencing error, we usually put the MAF at 1/number of
  samples : 1/(1133-68)=0.000938967

Applying filters to VCF

``` bash
VCF_IN=Pvivax_total_snp.ploidy2.vcf.gz
VCF_OUT=Pvivax_total_snpbi_filtered.ploidy2.vcf.gz

# set filters
MAF=0.000938967
MISS=0.8
QUAL=30
MIN_DEPTH_SNP=10
MAX_DEPTH_SNP=106
MIN_DEPTH=10
MAX_DEPTH=164

vcftools --gzvcf $VCF_IN --remove remove_miss.txt --min-alleles 2 --max-alleles 2 \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH_SNP --max-meanDP $MAX_DEPTH_SNP \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > $VCF_OUT
```

What have we done here?

- `--remove remove_miss.txt` - remove all all the individuals with more
  than 50% of missing data
- `--min-alleles 2 --max-alleles 2` - keep only the bi-allelic SNPs
- `--remove-indels` - remove all indels (SNPs only)
- `--maf`- set minor allele frequency - here 1/number of samples
- `--max-missing` - set minimum missing data. A little counter
  intuitive - 0 is totally missing, 1 is none missing. Here 0.8 means we
  will tolerate 20% missing data.
- `--minQ` - this is just the minimum quality score required for a site
  to pass our filtering threshold. Here we set it to 30.
- `--min-meanDP` - the minimum mean depth for a site.
- `--max-meanDP` - the maximum mean depth for a site.
- `--minDP` - the minimum depth allowed for a genotype - any individual
  failing this threshold is marked as having a missing genotype.
- `--maxDP` - the maximum depth allowed for a genotype - any individual
  failing this threshold is marked as having a missing genotype.

## Remove multi-clonal infections

Version: bcftools v1.10.2, vcftools v0.1.16,
[vcfdo](https://github.com/IDEELResearch/vcfdo).

> The F<sub>WS</sub> metric estimates the heterozygosity of parasites
> (HW) within an individual relative to the heterozygosity within a
> parasite population (HS) using the read count of alleles.
> F<sub>WS</sub> metric calculation for each sample was performed using
> the following equation: F<sub>WS</sub>=1− HW/HS where HW refers to the
> allele frequency of each unique allele found at specific loci of the
> parasite sequences within the individual, and HS refers to the
> corresponding allele frequencies of those unique alleles within the
> population. F<sub>WS</sub> ranges from 0 to 1; a low F<sub>WS</sub>
> value indicates low inbreeding rates within the parasite population
> and thus high within-host diversity relative to the population. An
> F<sub>WS</sub> threshold ≥ 0.95 indicates samples with clonal (single
> strain) infections, while samples with an F<sub>WS</sub> \< 0.95 are
> considered highly likely to come from mixed strain infections,
> indicating within-host diversity.

Source : [Amegashie et
al. (2020)](https://doi.org/10.1186/s12936-020-03510-3).

The F<sub>WS</sub> must be calculated by population. So we split up the
VCF by country. We only calculated F<sub>WS</sub> for the modern
samples.

``` bash
conda activate vcfdo
while read file_sample
do 
vcftools --gzvcf Pvivax_total_snpbi_filtered.ploidy2.vcf.gz --keep ./countries/$file_sample.tsv --recode --stdout | gzip -c > ./vcf_countries/$file_sample.vcf.gz #create the file by country
echo "Sample  Fws Standard_errors nb_sites" > ./fws/fws_$file_sample.txt
vcfdo wsaf -i ./vcf_countries/$file_sample.vcf.gz | vcfdo fws >> ./fws/fws_$file_sample.txt #calculate fws
done < ./list_countries.txt
```

We remove the individuals with a F<sub>WS</sub> \> 0.95 : we keep 745
individuals.

## Remove related samples

Version: bcftools v1.10.2, vcftools v0.1.16,
[hmmIBD](https://github.com/glipsnort/hmmIBD), R v4.2.

> Highly related samples and clones can generate spurious signals of
> population structure, bias estimators of population genetic variation,
> and violate the assumptions of the model-based population genetic
> approaches ([Wang 2018](https://doi.org/10.1111/1755-0998.12708)). The
> relatedness between haploid genotype pairs was measured by estimating
> the pairwise fraction of the genome identical by descent (*IBD*)
> between strains within populations.

The IBD must be calculated by countries (only with n\>=2).

We used hmmIBD that as a specific format as an input, with only haploid
information. For the sites that was heterozygote, they were marked as
missing data.

``` bash
while read country_name
do
vcftools --gzvcf Pvivax_total_snpbi_filtered.ploidy2.vcf.gz --keep ./countries/$file_sample.tsv --remove remove_fws_only.txt --recode --stdout | gzip -c > ./vcf_countries/${country_name}.vcf.gz #create the file by country and remove multi-clonal samples

# Transform in hmmIBD format
bcftools annotate -x INFO,^FORMAT/GT ./vcf_countries/${country_name}.vcf.gz | grep -v "##" |cut -d$'\t' -f1-2,10- > temp0_${country_name}

sed 's/0\/0/0/g' temp0_${country_name} > temp1_${country_name}
sed 's/1\/1/1/g' temp1_${country_name} > temp2_${country_name}
sed 's/1\/0/-1/g' temp2_${country_name} > temp3_${country_name}
sed 's/0\/1/-1/g' temp3_${country_name} > temp4_${country_name}
sed 's/0\/1/-1/g' temp4_${country_name} > temp5_${country_name}
sed 's/.\/./-1/g' temp5_${country_name} > temp6_${country_name}
sed 's/PvP01_\([0-9][0-9]*\)_v1/\1/g' temp6_${country_name} > temp7_${country_name} #Change chromosome name to chromosome number
sed 's/vPvP01_\([0-9][0-9]*\)_v1/\1/g' temp7_${country_name} > ./hmm_format/${country_name}_hmm.pf

#Calculate IBD
hmmIBD -i ./hmm_format/${country_name}_hmm.pf -o ./IBD/IBD_${country_name}
done < ./list_countries.txt
```

Isolate pairs that shared \>50% of IBD are considered highly related. In
each family of related samples, only the strain with the lowest amount
of missing data was retained:

``` r
# Read all IBD files
IBD_all <-
    list.files(path="./Data/IBD/",
               pattern = "*.hmm_fract.txt", 
               full.names = T) %>% 
    map_dfr(~read_table(.), show_col_types=F)

# Add meta-informations
metadata_vivax<-read_delim("./Data/metadata_vivax_fws.csv", 
     delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE, show_col_types = FALSE)
metadata_vivax<-metadata_vivax[,c(1,2,9,16)]
colnames(metadata_vivax)<-c("sample_ID","species","country","population")
metadata_vivax$sample1<-metadata_vivax$sample_ID
IBD_all<-inner_join(IBD_all, metadata_vivax)

total_list<-unique(c(IBD_all$sample1,IBD_all$sample2)) # List of all the samples
# Only keep pair of individuals with IBD>0.5
fam_IBD<-subset(IBD_all, IBD_all$fract_sites_IBD>0.5)

#Assign family factor by individuals
clst = data.frame(ind = c(as.character(fam_IBD$sample1[1]), as.character(fam_IBD$sample2[1])), grp = c(1,1)) # initialize data.frame
clst
for(i in 2:dim(fam_IBD)[1]){
  if(length(which(as.character(fam_IBD$sample1[i])==clst$ind))>0){
    tmp = data.frame(ind = c(as.character(fam_IBD$sample1[i]), as.character(fam_IBD$sample2[i])), grp = c(clst$grp[which(as.character(fam_IBD$sample1[i])==clst$ind)],clst$grp[which(as.character(fam_IBD$sample1[i])==clst$ind)]))
    clst = rbind(clst, tmp)
  } else if(length(which(as.character(fam_IBD$sample2[i])==clst$ind))>0){
    tmp = data.frame(ind = c(as.character(fam_IBD$sample1[i]), as.character(fam_IBD$sample2[i])), grp = c(clst$grp[which(as.character(fam_IBD$sample2[i])==clst$ind)],clst$grp[which(as.character(fam_IBD$sample2[i])==clst$ind)]))
    clst = rbind(clst, tmp)
  } else {
    tmp = data.frame(ind = c(as.character(fam_IBD$sample1[i]), as.character(fam_IBD$sample2[i])), grp = c(max(clst$grp)+1,max(clst$grp)+1))
    clst = rbind(clst, tmp)
  }
  clst = unique(clst)
}

# import the information of missing data (from vcftools, see above)
ind_miss  <- read_delim("./Data/Pvivax_total.ploidy2.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1, show_col_types = FALSE)
data_fam<-inner_join(clst, ind_miss)

### Remove with IBD only
#keep the individual in each family with the less missing data
unrelated<-data_fam %>% 
    group_by(grp) %>% 
    slice(which.min(fmiss))
```

## Create the final dataset filtered

Keep only the samples of the analysis dataset (mono-clonal and not
inbred): 622 samples.

``` bash
vcftools --gzvcf Pvivax_total_snpbi_filtered.ploidy2.vcf.gz --remove remove_fws_only.txt --remove remove_IBD_only.txt --recode --stdout | gzip -c > Pvivax_filtered_final.ploidy2.vcf.gz
```
