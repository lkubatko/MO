library(plyr)
library(ape)

parameter_setting = expand.grid(alpha=c(0.05,0.1,0.2,0.4),
                                beta =c(0.05,0.1,0.2,0.4),
                                missingindex=c(0,10,20))
# columns for correct inferred pairs,true pairs,numSites,alpha,beta,missing index, genotype, Method="MO"

all_pair_result_dat= c()


for(parameterindex in 1:dim(parameter_setting)[1]){
  
  print(parameterindex)
  alpha =   parameter_setting[parameterindex,1]
  beta = parameter_setting[parameterindex,2]
  missingindex = parameter_setting[parameterindex,3]
  genotype = as.character(parameter_setting[parameterindex,4])
  genotype_lowercase = tolower(genotype)
  
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


if(missingindex==0){
  
  individual_res = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario1_result_SiFit/SiFit_SimulateData/summary/binary_ALL_order_matrix_alpha_0%s_beta_0%s.csv",alpha_str,beta_str)
  
}else{
  
  individual_res = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario1_result_SiFit/SiFit_SimulateData_%smissing/summary/binary_ALL_order_matrix_alpha_0%s_beta_0%s.csv",missingindex,alpha_str,beta_str)
  
}

  all_pair_result = read.csv(individual_res)

  all_pair_result_complete = cbind(all_pair_result[,c(2:4)],
                                   numSites = rep(18,dim(all_pair_result)[1]), 
                                   alpha_para = rep(alpha,dim(all_pair_result)[1]),
                                   beta_para = rep(beta,dim(all_pair_result)[1]),
                                   missing_index = rep(missingindex,dim(all_pair_result)[1]), 
                                   genotype_para = rep("Binary",dim(all_pair_result)[1]), 
                                   Method=rep("SiFit",dim(all_pair_result)[1]))
  colnames(all_pair_result_complete)[1:3] = c("inferred_correct","inferred_all","true")
  
  all_pair_result_dat = rbind(all_pair_result_dat,all_pair_result_complete)


}

all_pair_result_form_0_1 = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario1_result_SiFit/accuracy/scenario1_all_order_accuracy.csv')
write.csv(all_pair_result_dat,file=all_pair_result_form_0_1)
