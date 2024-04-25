Test scenarios with ABC method
================
Margaux Lefebvre
2024-04-25

# Create the input files

Version: PLINK v2, bcftools v1.16, python v2.7.

``` bash
# Remove missing data from Ebro
bcftools view Pvivax_filtered_final.ploidy2.vcf.gz -s Ebre_composite -R core_genome.txt -o Ebro_only.vcf.gz -O z
tabix Ebro_only.vcf.gz

bcftools view -g ^miss Ebro_only.vcf.gz -Oz -o Ebro_only.no_miss.vcf.gz

bcftools query -f '%CHROM\t%POS\n' Ebro_only.no_miss.vcf.gz > SNP_position.txt

bcftools view Pvivax_filtered_final.ploidy2.vcf.gz -S subset.pop -R SNP_position.txt -Oz -o core_genome.ebro_no_miss.vcf.gz # keep only Ebro, Colombia and Mauritania
tabix core_genome.ebro_no_miss.vcf.gz

# LD pruning
plink2 --vcf core_genome.ebro_no_miss.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.5 --out Prune

sed -e ' s/:/\t/g' Prune.prune.in > bcft_Prune.prune.in

bcftools view core_genome.ebro_no_miss.vcf.gz -R bcft_Prune.prune.in -o core_genome.pruned.nomiss.vcf -c 2 -C 2
```

Create the info for the population

``` r
library(readr)
pop_file <- read_delim("./pop_file.txt", 
    delim = "\t", escape_double = FALSE, 
    col_names = FALSE, trim_ws = TRUE)

pop_file_o <- cbind(pop_file$X1, 9, pop_file$X2)
write_tsv(as.data.frame(pop_file_o), file="./Data/pop_file_DIY.txt", col_names = F)
```

Conversion (script from: <https://github.com/loire/vcf2DIYABC.snp.git>)

``` bash
python vcf2DIYABC.snp/vcf2diyabc.py core_genome.pruned.nomiss.vcf pop_file_DIY.txt #Writing outputfile as: core_genome.pruned.nomiss.vcf.DIYABC.snp
```

Check missing data

``` r
library(readr)
data_DIYABC <- read_table("./Data/core_genome.pruned.nomiss.vcf.DIYABC.snp", 
    col_names = FALSE)
data_geno<-data_DIYABC[2:48,4:5413]


# remove sites with more than 40% missing data
misssites<-colSums(data_geno =="9")*100/47
misssites[misssites>=40]<-NA
data_test<-rbind(data_geno, misssites)

new_df <- data_test[ , colSums(is.na(data_test))==0]

#remove individuals with more than 40% of missing data
new_df$missind<-(rowSums(new_df =="9")*100/5410)
remove_ind<-which(new_df$missind >=40)

data_geno_sites<-new_df[1:47,1:ncol(new_df)-1]

data_geno_info<-cbind(data_DIYABC[2:48,1:3], data_geno_sites)
data_geno_noimiss<-data_geno_info[-remove_ind,]
header<-data_DIYABC[1,1:ncol(data_geno_noimiss)]
colnames(header)<-colnames(data_geno_noimiss)
data_final<-rbind(header, data_geno_noimiss)
nrow(subset(data_final, X3=="Mauritania"))

write_tsv(data_final, "./Data/core_genome.pruned.ebro2.nomiss.haplo.clear.clean.DIYABC.txt", col_names = F)
```

# Test the scenarios

Version: abcranger v1.16.47, R v4.2.

I used the GUI interface.

``` r
library(tidyverse)
library(diyabcGUI)
diyabcGUI::diyabc()
```

I launched 10 replicates

``` bash
for i in {1..10}
do
~/Softwares/abcranger-linux-v1.16.47 -o RF_replicate$i -t 1500 -j 10
done
```

# Infer the best parameters

Version: abcranger v1.16.47.

``` bash
PARAM=(Nanc NAncEbro Nbot NCol NEbro Nghostam NMauri ra tbot tct2 tdivebraf tdivebram)
for par in ${PARAM[@]}
do
~/Softwares/abcranger-linux-v1.16.47 -o est_$par -t 1500 -j 10 --chosenscen 1 --noob 1000 --parameter $par
done
```
