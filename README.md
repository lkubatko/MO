# MO
Code and example data for the paper "A Phylogenetic Approach to Inferring the Order in Which Mutations Arise during Cancer Progression", co-authored by Yuan Gao, Jeff Gaither, and Julia Chifman

## Installation

`MO` has several dependencies:
```
install.packages(c("igraph","ape","phangorn"))
```
Once the dependencies are installed, `MO` can be installed with:
```
devtools::install_github("lkubatko/MO",ref="MO_pkg")
```
## Input Files

### Genotype Matrix

The mutational profile of single cells are represented as a genotype matrix `initial_obs` and it has to be provided as input. The first column is the gene. Starting from the second column, each of the columns represents the mutation profile of one single cell. Each row represents mutation at one genomic site of all single cells. The entries represent the genotype of the corresponding cell for the given mutation. The genotype matrix can be binary/ternary.

a)Binary: The genotype information represents the presence/absence of a mutation. The entry at position [i, j+1] should be

0 if mutation i is not observed in cell j,
1 if mutation i is observed in cell j, 
3 if the genotype information is missing,or
m if the genotype information is ambiguous

b)Ternary: The genotype information represents homozygous reference, heterozygous variant and homozygous variant respectively. The entry at position [i, j+1] should be

0 if mutation i is not observed in cell j,
1 if mutation i is heterozygous variant in cell j,
2 if mutation i is homozygous variant in cell j, 
3 if the genotype information is missing, or
m if the genotype information is ambiguous


### The Input Tree 

A phylogenetic tree `t` is required as an input. The computation of probabilities of ordering is based on the input tree. Many methods can be used to infer the tree topology and branch length.

### Arguments
`alpha` and `beta` are set to the estimated false positive rate and estimated allelic dropout rate of the single-cell sequencing experiment. `binary` is the logical parameter of whether the genotype data is binary(TRUE) or ternary(FALSE). `multinomialPrior` is the logical parameter of whether using the multinomial prior(TRUE) or not(FALSE). lambda01, lambda02 and lambda12 are rate of mutating from 0 to 1, 0 to 2, and 1 to 2, respectively. 

### Example data and Simulation data
Example data are in folder example, and simulation data are in folder simulation. For example, if you use CRC1 patient as an example, you can run 
```
eg = MO(alpha, beta, lambda01, lambda02, lambda12, initial_obs = CRC1_mutation_data, t=crc1tree, binary=TRUE, multinomialPrior=FALSE)
```

## Issues using MO?
MO is currently in beta. We expect there to be bumps in the road. If you have questions about MO usage, please submit an issue on Github or send us an email. We will do our best to respond to questions that are not otherwise answered in the documentation.
