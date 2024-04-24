Population structure
================
Margaux Lefebvre
2024-04-24

# PCA

## LD-pruning

Version: samtools v1.10, angsd v0.940, python v3.8.12.

NB : I removed the outgroup, since it doesn’t make sense to keep *P.
vivax-like* samples.

``` bash
# Create Beagle input for the pruning
angsd -b list_bams_noPL.txt -ref P.vivax-PvP01_reference_genome.fasta -out PCA_coregenome \
        -rf core_regions.txt \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepthInd 106 -doCounts 1 \
        -GL 2 -doGlf 2 -nThreads 8 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -minMaf 0.05

# Prepare a pos file with the mafs filre
zcat PCA_coregenome.mafs.gz | cut -f 1,2 |  sed 's/:/_/g'| gzip > PCA_coregenome.pos.gz

# Run ngsLD
ngsLD \
--geno PCA_coregenome.beagle.gz \
--posH PCA_coregenome.pos.gz \
--probs \
--n_ind 620 \
--n_sites 130205 \
--max_kb_dist 1 \
--n_threads 8 \
--out PCA_coregenome.ld

python ./scripts/prune_ngsLD.py \
--input PCA_coregenome.ld  \
--max_dist 5000 \
--min_weight 0.5 \
--output PCA_coregenome.snp.unlinked.id
```

Generate an LD-pruned SNP list:

``` r
pruned_position <- as.integer(gsub(paste0("PvP01_[0-1][0-9]_v1:"), "", readLines(paste0("~/Documents/P.vivax_project/SAM_P.VIVAX/PCA/Data/PCA_coregenome.snp.unlinked.id"))))

snp_list <- read.table(paste0("./Data/PCA_coregenome.mafs.gz"), stringsAsFactors = F, header = T)[,1:4]

pruned_snp_list <- snp_list[snp_list$position %in% pruned_position, ]
  
write.table(pruned_snp_list, paste0("./Data/PCA_coregenome.snp.LDpruned.list"), col.names = F, row.names = F, quote = F, sep = "\t")
```

We keep 105,527 SNPs.

## PCA with PCAngsd

Version: angsd v0.940, PCAngsd v0.98.

``` bash
# Create input
angsd sites index PCA_coregenome.snp.LDpruned.list # mandatory

angsd -b list_bams_noPL.txt -ref P.vivax-PvP01_reference_genome.fasta -out PCA_coregenome_total \
        -rf core_regions.txt \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepthInd 106 -doCounts 1 \
        -GL 2 -doGlf 2 -nThreads 8 -doMajorMinor 1 -doMajorMinor 3 -doMAF 1 -doPost 1 -doIBS 1 -doCov 1 -makeMatrix 1 -sites PCA_coregenome.snp.LDpruned.list
        
# Calculate the PCA with PCAngsd
pcangsd -b PCA_coregenome_total.beagle.gz -o PCA_coregenome.mat
```

# Ancestry plots

Version: PCAngsd v0.98, pong v1.5.

Ancestry plots are inferred with PCAngsd and the input file is the same
as for PCA.

According to [Meisner and Albrechtsen](10.1534/genetics.118.30133), the
best K is determined by 1+ D (the optimal number of principal
components). As presented in Supplementary Fig. 4 in the paper, D would
be equal to 5 (determined by the elbow (broken-stick) method), so K=6.

``` bash
for k in {2..16}
do
 pcangsd -b PCA_coregenome_total.beagle.gz --admix --admix_K $k -o ancestry
done
```

The visualization was done with pong:

``` bash
pong -m file_map.txt -i ind2pop.txt -n country_order.txt -l color.txt 
```

# Tree of the samples (maximum likelihood)

Version: [vcf2phylip v2.0](https://doi.org/10.5281/zenodo.2540861),
python v2.7.5, iqtree v2.0.3.

I used *P. vivax-like* as an outgroup.

``` bash
# keep only the core genome
bcftools view Pvivax_filtered_final.ploidy2.vcf.gz -R core_genome.txt -o Pvivax_core.snpbi_filtered.ploidy2.vcf.gz

# change vcf to phylip file
python ./vcf2phylip/vcf2phylip.py -i PPvivax_core.snpbi_filtered.ploidy2.vcf.gz --output-prefix tree_allsamples -o p1537.PL.Cameroon

# ML tree
iqtree -s tree_allsamples.min4.phy -m MFP+ASC -T 24 --prefix tree_allsamples -o p1537.PL.Cameroon -B 1000 -alrt 1000 -st DNA #MFP+ASC = model finder for dataset with only variable sites
```
