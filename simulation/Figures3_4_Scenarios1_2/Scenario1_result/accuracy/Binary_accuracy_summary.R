library(ROCR)
library(pROC)
library(ape)


args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  stop("Error: No alpha and beta values are given!\n", call.=FALSE)
} else {
  alpha = as.double(args[1])
  beta = as.double(args[2])
  print(paste("alpha is ", alpha))
  print(paste("beta is ", beta))
}

if (alpha < 0.1)
{
  alpha_str = sprintf('0%s', alpha*100)
} else
{
  alpha_str = sprintf('%s', alpha*10)
}

if (beta < 0.1)
{
  beta_str = sprintf('0%s', beta*100)
} else
{
  beta_str = sprintf('%s', beta*10)
}

alpha01= alpha
alpha02= alpha*beta/2
beta10= beta/2
beta12= beta/2
gamma20 =0
gamma21 =0
sequencing_error_model=matrix(c(1-alpha01-alpha02,alpha01,alpha02,
                                beta10,1-beta10-beta12,beta12,
                                gamma20,gamma21,1-gamma21-gamma21),nrow=3,byrow = TRUE)
print(sequencing_error_model)
unit_theta = 10^(-7)
unit_gamma = 10^(-14)
unit_mu = 10 ^(-2)



multiROC_dat=data.frame(matrix(, nrow=1, ncol=0))
prob_cat_dat=data.frame(matrix(, nrow=0, ncol=5))

for (indexn in 1:1000){

  form = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/SimulateData/RandomTree/RandomTree_%s.tre', indexn)
  sampletr=read.tree(form)
  
  prob_form_0_1 = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/SimulateData/Binary_alpha0%s_beta0%s_result/binary_prob_matrix_all_0_1_out_alpha_0%s_beta_0%s_matrix%s.csv', alpha_str, beta_str, alpha_str, beta_str,indexn)
  
  
  mat_prob_form_0_1 = read.csv(prob_form_0_1)
  
  initial_obs_0_1 = as.matrix(mat_prob_form_0_1[,-c(1,2,3)])
  true_obs_0_1 = as.matrix(mat_prob_form_0_1[,c(2)])
  initial_obs_0_1[is.na(initial_obs_0_1)] <- 0
  

  
  colnames(initial_obs_0_1)=c(1:18)
  true_obs_0_1_recoded= true_obs_0_1
  
  prob_list=data.frame(prob=numeric(0))
  cat_list=data.frame(branch=numeric(0))
  indexn_list=data.frame(indexn=numeric(0))
  for (i in 1:dim(initial_obs_0_1)[1]){
    
    print(c(indexn,i))
    indexn_list[i,]=indexn
    possible_br = which.max(initial_obs_0_1[i,])
    
    if (possible_br == true_obs_0_1_recoded[i]){prob_list[i,] = max(initial_obs_0_1[i,]) 
    cat_list[i,] = 1}
    else{prob_list[i,] = max(initial_obs_0_1[i,]) 
    cat_list[i,] = 0}
    
    
  }
  
  prob_cat_list = cbind(mat_prob_form_0_1[,c(2,3)],cat_list,prob_list,indexn_list)
  colnames(prob_cat_dat)=colnames(prob_cat_list)
  prob_cat_dat=rbind(prob_cat_dat,prob_cat_list)
  #multiROC = multiclass.roc(true_obs_0_1_recoded, initial_obs_0_1)
  #multiROC_dat[,indexn] = multiROC$auc
  
  
}

#multiROC_dat_0_1_out = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/SimulateData/Binary_Summary/binary_multiROC_dat_0_1_out_alpha_0%s_beta_0%s_matrix.csv', alpha_str, beta_str)
#write.csv(multiROC_dat,file=multiROC_dat_0_1_out)
setwd('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/SimulateData/Binary_Summary/')
prob_cat_list_0_1_out = sprintf('binary_prob_cat_list_dat_0_1_out_alpha_0%s_beta_0%s_matrix.csv', alpha_str, beta_str)

write.csv(prob_cat_dat,file=prob_cat_list_0_1_out)


