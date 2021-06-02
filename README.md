# MO
Code and example data for the paper "A Phylogenetic Approach to Inferring the Order in Which Mutations Arise during Cancer Progression", co-authored by Yuan Gao, Jeff Gaither, and Julia Chifman

Input Files

Genotype Matrix

The mutational profile of single cells are represented as a genotype matrix and it has to be provided as input. Each of the columns represents the mutation profile of one single cell. Each row represents mutation at one genomic site of all single cells. The entries represent the genotype of the corresponding cell for the given mutation. The genotype matrix can be binary/ternary.

a)Binary: The genotype information represents the presence/absence of a mutation. The entry at position [i, j] should be

0 if mutation i is not observed in cell j,
1 if mutation i is observed in cell j, 
3 if the genotype information is missing,or
m if the genotype information is ambiguous

b)Ternary: The genotype information represents homozygous reference, heterozygous variant and homozygous variant respectively. The entry at position [i, j] should be

0 if mutation i is not observed in cell j,
1 if mutation i is heterozygous variant in cell j,
2 if mutation i is homozygous variant in cell j, 
3 if the genotype information is missing, or
m if the genotype information is ambiguous


The Input Tree 
An inferred tree is required as an input. The computation of probabilities of ordering is based on the input tree. Many methods can be used to infer the tree topology and branch length.

Arguments
alpha and beta are set to the estimated false positive rate and estimated allelic dropout rate of the single-cell sequencing experiment. 
