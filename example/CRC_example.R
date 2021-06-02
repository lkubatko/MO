library(ape)
library(phangorn)

setwd("G:/Github_MO")
source("code_MO_binary_function.R")
crc1tree=read.tree("crc1.tre")

cells=crc1tree$tip.label
genotype_178_cells_gt=read.csv("CRC1_mutation_data.csv")

genotype_178_initial_obs=genotype_178_cells_gt[c(cells)]


genotype_178_normal_genotype=rep(0,16)

genotype_178_mutation_genotype=rep(1,16)


alpha = 1.52/100
beta = 7.89/100
alpha01= alpha
alpha02= alpha*beta/2
beta10= beta/2
beta12= beta/2


P1_binary_prob_matrix_0_1=c()


for (i in 1:dim(genotype_178_initial_obs)[1]){
  
  print(i)
  rd_unit_lambda <- rgamma(10000, shape = 2*10^(-0), scale = (10^(-0)/2)) 
  
  generate_prob_br_0_1_dat=data.frame(matrix(NA, nrow = 144, ncol = 10))
  
  for (indexj in 1:10){
    
    generate_prob_br <- MO_binary(alpha,beta,rd_unit_lambda[indexj],genotype_178_initial_obs[i,],crc1tree)
      
      
    generate_prob_br_0_1_single <- generate_prob_br
    
    
    generate_prob_br_0_1_dat[,indexj] = generate_prob_br_0_1_single
    
  }
  
  generate_prob_br_0_1=rowMeans(generate_prob_br_0_1_dat, na.rm = FALSE, dims = 1)
  
  P1_binary_prob_matrix_0_1 = rbind(P1_binary_prob_matrix_0_1,generate_prob_br_0_1)
  
}



#CRC2
crc2tree = read.tree("crc2.tre")

cells=crc2tree$tip.label
genotype_182_cells_gt=read.csv("CRC2_mutation_data.csv")

genotype_182_initial_obs=genotype_182_cells_gt[c(cells)]

genotype_182_normal_genotype=rep(0,36)

genotype_182_mutation_genotype=rep(1,36)



alpha = 1.74/100
beta = 12.56/100
alpha01= alpha
alpha02= alpha*beta/2
beta10= beta/2
beta12= beta/2

P1_binary_prob_matrix_0_1=c()


for (i in 1:dim(genotype_182_initial_obs)[1]){
  
  print(i)
  rd_unit_lambda <- rgamma(10, shape = 2*10^(-0), scale = (10^(-0)/2)) 
  
  generate_prob_br_0_1_dat=data.frame(matrix(NA, nrow = 172, ncol = 10))
  
  for (indexj in 1:10){
    
    generate_prob_br <- MO_binary(alpha,beta,rd_unit_lambda[indexj],genotype_182_initial_obs[i,],crc2tree)
    
    
    generate_prob_br_0_1_single <- generate_prob_br
    
    
    generate_prob_br_0_1_dat[,indexj] = generate_prob_br_0_1_single
    
  }
  
  generate_prob_br_0_1=rowMeans(generate_prob_br_0_1_dat, na.rm = FALSE, dims = 1)
  
  P1_binary_prob_matrix_0_1 = rbind(P1_binary_prob_matrix_0_1,generate_prob_br_0_1)
  
}



selected_br=rep(0,36)
max_prob=rep(0,36)
for(i in 1:dim(P1_binary_prob_matrix_0_1)[1]){
  selected_br[i]=which.max(P1_binary_prob_matrix_0_1[i,])
  max_prob[i]=max(P1_binary_prob_matrix_0_1[i,])
}

P1_binary_prob_matrix_0_1=cbind(selected_br,P1_binary_prob_matrix_0_1)
rownames(P1_binary_prob_matrix_0_1)=genotype_182_cells_gt[,1]





