# Unraveling the journey of Plasmodium vivax in Latin America: a genomic approach

This repository is for this paper:

**Genomic exploration of the journey of *Plasmodium vivax* in Latin America**, M. J. M. Lefebvre, F. Degrugillier, C. Arnathau, G. A. Fontecha, O. Noya, S. Houzé, C. Severini, B. Pradines, A. Berry, J-F. Trape, F. E. Sáenz, F. Prugnolle, M. C. Fontaine, V. Rougeron. *PLOS Pathogens* 21(1): e1012811. 

doi: <https://doi.org/10.1371/journal.ppat.1012811>

The languages used are mainly bash and R. For each part, software and version are specified.

# Summary

## [About the dataset...](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Dataset_generation#about-the-dataset)

-   [The sources](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Dataset_generation#the-sources)
-   [Ebro: take care of ancient DNA samples](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Dataset_generation#ebro-take-care-of-ancient-dna-samples)
-   [Produce unfiltered VCF with all samples](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Dataset_generation#produce-unfiltered-vcf-with-all-samples)
-   [Filtering of the dataset](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Dataset_generation#filtering-of-the-dataset)

## [Population structure](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Pop_structure#population-structure)

-   [PCA](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Pop_structure#pca)
-   [Ancestry plots](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Pop_structure#ancestry-plots)
-   [Tree of the samples (maximum likelihood)](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Pop_structure#tree-of-the-samples-maximum-likelihood)

## [Population graphs and *f*-statistics](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Pop_graphs_fstats#population-graphs-and-f-statistics)

-   [TreeMix](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Pop_graphs_fstats#treemix)
-   [AdmixtureBayes](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Pop_graphs_fstats#admixturebayes)
-   [*f~3~*-statistics](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Pop_graphs_fstats#f3-statistics)

## [Demographic history](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Demography#demographic-history)

-   [Nucleotide diversity (pi) and Tajima's D](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Demography#nucleotide-diversity-pi-and-tajimas-d)
-   [Coalescence with Relate and Colate](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Demography#coalescence-with-relate-and-colate)

## [Test scenarios with ABC method](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/ABC_scenarios#test-scenarios-with-abc-method)

-   [Create the input files](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/ABC_scenarios#create-the-input-files)
-   [Test the scenarios](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/ABC_scenarios#test-the-scenarios)
-   [Infer the best parameters](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/ABC_scenarios#infer-the-best-parameters)

## [Supplementary](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Supp#supplementary)

-   [Genetic map](https://github.com/MargauxLefebvre/LatinAmerica_vivax/tree/main/Supp#genetic-map)
