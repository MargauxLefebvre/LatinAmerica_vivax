Supplementary
================
Margaux Lefebvre
2024-04-25

# Genetic map

As in [Brashear *et al.*
(2020)](https://doi.org/10.1371%2Fjournal.pntd.0008506), I will use
LDhat to infer my recombination map. For that, I will use one of the
most diverse population : Thailand (n=71).

Steps of the [LDhat pipeline](https://github.com/auton1/LDhat) and
[manual](https://raw.githubusercontent.com/auton1/LDhat/master/manual.pdf).

## Create input

Version: R v4.2, ANGSD v0.940, seqkit v2.1.0.

I kept samples with less than 5% of missing data (n=58).

Make fasta file for each sample

``` bash
echo "--> Start processing sample: $bam_name" #Deal with each chromosome separated
echo "--> Start processing chromosome: $CHR" #Deal with each chromosome separated

angsd -i $bam_name.bam -out ${bam_name}.${CHR} -r PvP01_${CHR}_v1 \
        -minMapQ 20 -minQ 20 -setMaxDepthInd 106 -doCounts 1 \
        -nThreads 8 -doFasta 2

gzip -cd ${bam_name}.${CHR}.fa.gz | sed "s/PvP01_${CHR}_v1/$bam_name/" | gzip > ${bam_name}.${CHR}.clean.fa.gz
```

Finish the input files for each chromosome

``` bash
echo "--> Start processing chromosome: $CHR" #Deal with each chromosome separated
for FILE in *.${CHR}.clean.fa.gz; do zcat $FILE >> haplo_Thai_${CHR}.fa; done

seqkit stats haplo_Thai_*.fa # get the number of sequences/genotypes for the headers

# Make the header for each chromosome
echo "58 1021664 1" > header_01.txt
echo "58 956327 1" > header_02.txt
echo "58 896704 1" > header_03.txt
echo "58 1012024 1" > header_04.txt
echo "58 1524814 1" > header_05.txt
echo "58 1042791 1" > header_06.txt
echo "58 1652210 1" > header_07.txt
echo "58 1761288 1" > header_08.txt
echo "58 2237066 1" > header_09.txt
echo "58 1548844 1" > header_10.txt
echo "58 2131221 1" > header_11.txt
echo "58 3182763 1" > header_12.txt
echo "58 2093556 1" > header_13.txt
echo "58 3153402 1" > header_14.txt

cat header_${CHR}.txt > Thai_$CHR.seq
cat haplo_Thai_${CHR}.fa >> Thai_$CHR.seq

gzip Thai_*.seq
```

## Conversion

If full sequence data is used, `convert` should be used to generate the
*sites* and *locs* files needed for subsequent analyses (these can be
renamed).

``` bash
LDhat=/path/to/LDhat
echo "--> Start processing chromosome: $CHR" #Deal with each chromosome separated

$LDhat/convert -seq Thai_${CHR}.seq -missfreqcut 0.2
mv sites.txt sites_${CHR}.txt
mv locs.txt locs_${CHR}.txt


# subset (chr 3, the one with the less nb of sites) without missing data just see the theta to use
$LDhat/convert -seq Thai_03.seq -missfreqcut 0
mv sites.txt sites_test.txt
mv locs.txt locs_test.txt
```

## Make the lookup table

``` bash
LDhat=/path/to/LDhat

# First estimation of theta
$LDhat/pairwise -seq sites_test.txt -loc locs_test.txt
# suggest Watterson estimate of 0.00011
$LDhat/complete -n 58 -rhomax 100 -n_pts 101 -theta 0.00011 -split 8
mv new_lk.txt Thai_lkt.txt
```

## Make the recombination map

``` bash
LDhat=/path/to/LDhat
echo "--> Start processing chromosome: $CHR" #Deal with each chromosome separated

$LDhat/interval -seq sites_${CHR}.txt -loc locs_${CHR}.txt -lk Thai_lkt.txt -its 1100000 -samp 100 -bpen 5
mv rates.txt rates_${CHR}.txt
mv bounds.txt bounds_${CHR}.txt
mv freqs.txt freq_${CHR}.txt

$LDhat/stat -burn 10000 -input rates_${CHR}.txt -loc locs_${CHR}.txt
mv res.txt res_${CHR}.txt
```

## Convertion of rho and genetic map

Version: R v4.2.

LDhat (interval) estimates recombination rates in the units of 4Ne r /
pb (here the position are coded in bp). Here, r is the recombination
rate in units of Morgans (M), so 0.01 is 1cM, meaning an expected number
of cross-over events per meiosis of 0.01. To convert from 4 Ne r / pb to
cM / Mb you need an estimate of Ne. According to Stairway plot, actual
median Ne of the Thai population is 6443.772 (with 95% confidence
interval of 1839.912 and 23019.76). Let’s say 6000 to simplify the
calculus. So to change *x* 4Ne r / pb in cM/Mb, I need to multiply *x*
by 10⁶ to get to units of Mb, by another hundred to go from M to cM and
divide by 24,000 to eliminate 4Ne).

For the genetic position, the formula is g= (pos - prevpos) \* 10⁻⁶ \*
rate (cM/MB) +g-1

``` r
library(readr)
ARR_CHR=c("01","02","03","04","05","06","07","08","09","10","11","12","13","14")
for (CHR in ARR_CHR){
res_recomb <- read_delim(paste0("./Data/res_",CHR,".txt"), 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE, show_col_types = FALSE)
res_recomb = res_recomb[-1,]
res_recomb$COMBINED_rate_cMMb<-res_recomb$Mean_rho*1000000*100/24000

res_recomb$Genetic_Map_cM<-NA
res_recomb$Genetic_Map_cM[1]<-0
for (i in 2:nrow(res_recomb)){
  res_recomb$Genetic_Map_cM[i]<-res_recomb$Genetic_Map_cM[i-1]+(res_recomb$Loci[i]-res_recomb$Loci[i-1])*res_recomb$COMBINED_rate_cMMb[i]/1000000
}
geneticmap<-as.data.frame(cbind(res_recomb$Loci, res_recomb$COMBINED_rate_cMMb, res_recomb$Genetic_Map_cM))
colnames(geneticmap)<-c("position","COMBINED_rate(cM/Mb)", "Genetic_Map(cM)")
write_tsv(geneticmap, file = paste0("./Data/Genetic_map_vivax_",CHR,".txt"))

# Make the map file too
geneticmap$chr<-CHR
map_file<-as.data.frame(cbind(geneticmap$chr, geneticmap$position, geneticmap$`Genetic_Map(cM)`, geneticmap$position))
write_tsv(map_file, file = paste0("./Data/mapfile_vivax_",CHR,".txt"), col_names = F)
}
```
