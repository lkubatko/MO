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
interval_prob_cat_dat=data.frame(matrix(, nrow=0, ncol=7))

for (indexn in 1:1000){
  
  form = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/SimulateData/RandomTree/RandomTree_%s.tre', indexn)
  sampletr=read.tree(form)
  
  prob_form_0_1 = sprintf('Obs_0_1_alpha0%s_beta0%s_result/ternary_prob_matrix_all_0_1_out_alpha_0%s_beta_0%s_matrix%s.csv', alpha_str, beta_str, alpha_str, beta_str,indexn)
  
  mat_prob_form_0_1 = read.csv(prob_form_0_1)
  
  #initial_obs_0_1 = as.matrix(mat_prob_form_0_1[,-c(1,2,3)])
  true_obs_0_1 = as.matrix(mat_prob_form_0_1[,c(2)])

  true_obs_0_1_recoded= true_obs_0_1
  
  
  interval_cat_list=data.frame(cat=numeric(0))
  interval_prob_list=data.frame(interval_prob=numeric(0))
  interval_maxprob_list=data.frame(max_prob=numeric(0))
  interval_num_list=data.frame(prob_num=numeric(0))
  
  interval_indexn_list=data.frame(indexn=numeric(0))
  for (i in 1:dim(mat_prob_form_0_1[,-c(1,2,3)])[1]){
    
    print(c(indexn,i))
    interval_indexn_list[i,]=indexn
    interval_maxprob_list[i,]=max(mat_prob_form_0_1[i,-c(1,2,3)])
    possible_br = which.max(mat_prob_form_0_1[i,-c(1,2,3)])
    reduced_mat_prob_form_0_1=(mat_prob_form_0_1[i,-c(1,2,3)])[,(which(colSums(mat_prob_form_0_1[i,-c(1,2,3)]) != 0))]
    
    true_obs_0_1_recoded_form=sprintf('X.%s',true_obs_0_1_recoded[i])
    
    
    
    sorted_prob = sort(reduced_mat_prob_form_0_1,decreasing = TRUE)
    sorted_br = colnames(sorted_prob)
    sorted_sum = cumsum(t(sorted_prob))
    selected_br_num = which.min(abs(sorted_sum - 0.95))
    selected_br = as.factor(sorted_br[1:selected_br_num])
    selected_prob = sorted_sum[selected_br_num]
    selected_possible_br = as.factor(colnames(reduced_mat_prob_form_0_1)[which.max(reduced_mat_prob_form_0_1)])
    
    if (any(levels(selected_br) == (true_obs_0_1_recoded_form))){
                          interval_prob_list[i,] = selected_prob 
                          interval_cat_list[i,] = 1
                          interval_num_list[i,] = selected_br_num}
    else{interval_prob_list[i,] = selected_prob
         interval_cat_list[i,] = 0
         interval_num_list[i,] = selected_br_num}
    
  }
  
  interval_prob_cat_list = cbind(mat_prob_form_0_1[,c(2,3)],interval_cat_list,interval_prob_list,interval_maxprob_list,interval_num_list,interval_indexn_list)
  colnames(interval_prob_cat_dat)=colnames(interval_prob_cat_list)
  interval_prob_cat_dat=rbind(interval_prob_cat_dat,interval_prob_cat_list)
  #multiROC = multiclass.roc(true_obs_0_1_recoded, initial_obs_0_1)
  #multiROC_dat[,indexn] = multiROC$auc
  
  
}

#multiROC_dat_0_1_out = sprintf('/users/PAS1462/gao957/workspace/SimulateData_20_EXP10/Ternary_Summary/ternary_multiROC_dat_0_1_out_alpha_0%s_beta_0%s_matrix.csv', alpha_str, beta_str)
#write.csv(multiROC_dat,file=multiROC_dat_0_1_out)

prob_cat_list_0_1_out = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/SimulateData/Ternary_Summary/ternary_interval_prob_cat_list_dat_0_1_out_alpha_0%s_beta_0%s_matrix.csv', alpha_str, beta_str)

write.csv(interval_prob_cat_dat,file=prob_cat_list_0_1_out)


