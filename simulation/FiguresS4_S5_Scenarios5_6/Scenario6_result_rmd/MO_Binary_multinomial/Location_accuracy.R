
library(ape)
parameter_setting = expand.grid(alpha=c(0.05,0.1,0.2,0.4),
                                sites =c(20,40,80),
                                missingindex=c(0,10,20))
set.seed(1000)

####################################################################################
#for binary data location accuracy
for(parameterIndex in 1:dim(parameter_setting)[1]){
  #get the parameters in each setting
  alpha =   parameter_setting[parameterIndex,1]#fpr
  beta = parameter_setting[parameterIndex,1]#fnr
  numSites = parameter_setting[parameterIndex,2]
  missingindex = parameter_setting[parameterIndex,3]
  #number_br = 18 #number of branches in a tree
  
  #return alpha and beta string used for data input
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
 

prob_cat_dat=data.frame(matrix(, nrow=0, ncol=4))

for (indexn in 1:100){
  #form = sprintf('./Figures3_4_Scenarios1_2/Scenario1/RandomTree/RandomTree_%s.tre', indexn)
  #sampletr=read.tree(form)
  print(c(parameterIndex,indexn))
  
  if (indexn < 10){
    
    trueTree_form=sprintf("./FiguresS4_S5_Scenarios5_6/Scenario6/results_alpha_0%s_beta_0%s_%s/trees.000%s",alpha_str,alpha_str,numSites,indexn)
    
  }else if( indexn>=10 & indexn<100){
    
    trueTree_form=sprintf("./FiguresS4_S5_Scenarios5_6/Scenario6/results_alpha_0%s_beta_0%s_%s/trees.00%s",alpha_str,alpha_str,numSites,indexn)
    
    
  }else{
    
    trueTree_form=sprintf("./FiguresS4_S5_Scenarios5_6/Scenario6/results_alpha_0%s_beta_0%s_%s/trees.0%s",alpha_str,alpha_str,numSites,indexn)
    
  }
  
  
  #scanned_trueTree_form = scan(file=trueTree_form,what=character(), n = -1, sep = "")
  
  #sampletr = read.tree(text=scanned_trueTree_form)
  
  if (missingindex == 0){
    
    prob_form_0_1 = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/MO/simulation/FiguresS4_S5_Scenarios5_6/Scenario6_result_rmd/MO_Binary_multinomial/Binary_alpha0%s_beta0%s_%s_result/all_binary_prob_matrix_all_0_1_out_matrix%s_multinomial.csv",alpha_str,alpha_str,numSites,indexn)
    
  }else{
    
    prob_form_0_1 = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/MO/simulation/FiguresS4_S5_Scenarios5_6/Scenario6_result_rmd/MO_Binary_multinomial/Binary_alpha0%s_beta0%s_%s_result/all_binary_prob_matrix_all_0_1_out_matrix%s_%s_multinomial.csv",alpha_str,alpha_str,numSites,indexn,missingindex)
    
  }
  mat_prob_form_0_1 = read.csv(prob_form_0_1)
  
  #initialize a vector cat to store result: 0 means MAP estimate is different from true value, 1 means MAP estimate is the same as true value
  cat=rep(0,dim(mat_prob_form_0_1)[1])
  for (i in 1:dim(mat_prob_form_0_1)[1]){if(mat_prob_form_0_1$selected_First_branch[i] == mat_prob_form_0_1$selected_br[i]){cat[i] = 1}}
  
  prob_cat_list = cbind(mat_prob_form_0_1[,c(2,3,4)],cat)
  colnames(prob_cat_dat)=colnames(prob_cat_list)
  prob_cat_dat=rbind(prob_cat_dat,prob_cat_list)

}

prob_cat_list_0_1_out = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/MO/simulation/FiguresS4_S5_Scenarios5_6/Scenario6_result_rmd/MO_Binary_multinomial/summary/binary_location_alpha_0%s_beta_0%s_%s_%smissing.csv', alpha_str, beta_str,numSites,missingindex)
write.csv(prob_cat_dat,file=prob_cat_list_0_1_out)
}


####################################################################################
#summary of location accuracy
library(plyr)
library(ape)

parameter_setting = expand.grid(alpha=c(0.05,0.1,0.2,0.4),
                                missingindex=c(0,10,20),
                                genotype = c("Binary"),
                                sites =c(20,40,80))

all_pair_result_dat= c()


for(parameterindex in 1:dim(parameter_setting)[1]){
  
  print(parameterindex)
  alpha =   parameter_setting[parameterindex,1]
  beta = parameter_setting[parameterindex,1]
  missingindex = parameter_setting[parameterindex,2]
  genotype = as.character(parameter_setting[parameterindex,3])
  genotype_lowercase = tolower(genotype)
  numSites = parameter_setting[parameterindex,4]
  
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
    
    individual_res = sprintf("./FiguresS4_S5_Scenarios5_6/Scenario6_result_rmd/MO_Binary_multinomial/summary/%s_location_alpha_0%s_beta_0%s_%s_0missing.csv",genotype_lowercase,alpha_str,beta_str,numSites)
    
  }else if(missingindex==10){
    
    individual_res = sprintf("./FiguresS4_S5_Scenarios5_6/Scenario6_result_rmd/MO_Binary_multinomial/summary/%s_location_alpha_0%s_beta_0%s_%s_10missing.csv",genotype_lowercase,alpha_str,beta_str,numSites)
    
  }else{
    
    individual_res = sprintf("./FiguresS4_S5_Scenarios5_6/Scenario6_result_rmd/MO_Binary_multinomial/summary/%s_location_alpha_0%s_beta_0%s_%s_20missing.csv",genotype_lowercase,alpha_str,beta_str,numSites)
    
  }
  
  all_pair_result = read.csv(individual_res)
  
  all_pair_result_complete = data.frame(inferred = all_pair_result[,c("cat")],
                                        true = rep(1,dim(all_pair_result)[1]),
                                        numSitesAll = rep(numSites,dim(all_pair_result)[1]), 
                                        alpha_para = rep(alpha,dim(all_pair_result)[1]),
                                        beta_para = rep(beta,dim(all_pair_result)[1]),
                                        missing_index = rep(missingindex,dim(all_pair_result)[1]), 
                                        genotype_para = rep(genotype,dim(all_pair_result)[1]), 
                                        Method=rep("MO",dim(all_pair_result)[1]))
  
  
  colnames(all_pair_result_complete)[1:2] = c("inferred","true")
  
  all_pair_result_dat = rbind(all_pair_result_dat,all_pair_result_complete[1:100,])
  
  
}

all_pair_result_form_0_1 = sprintf('./FiguresS4_S5_Scenarios5_6/Scenario6_result_rmd/accuracy/scenario6_location_accuracy_multinomial.csv')
write.csv(all_pair_result_dat,file=all_pair_result_form_0_1)
