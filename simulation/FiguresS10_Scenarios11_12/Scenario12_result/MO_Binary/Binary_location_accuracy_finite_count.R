#library(ROCR)
#library(pROC)
library(ape)


parameter_setting = expand.grid(m10=c(0.01,0.1,1),
                                alpha=c(0.05,0.1,0.2,0.4),
                                sites =c(20,40,80))

prob_cat_dat= data.frame(matrix(, nrow=0, ncol=10))

for(paraInd in 1:dim(parameter_setting)[1]){
  
  
  m10 = parameter_setting[paraInd,1] 
  alpha =   parameter_setting[paraInd,2]
  beta = parameter_setting[paraInd,2]
  numSites= parameter_setting[paraInd,3]
  
  
  if (m10 < 0.1)
  {
    m10_str = sprintf('00%s', m10*100)
  } else if (m10 >= 0.1 & m10<1){
    m10_str = sprintf('0%s', m10*10)
  } else{
    m10_str = sprintf('%s', m10)
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
  
  
  
  prob_cat_list= data.frame(matrix(, nrow=0, ncol=10))
  for (indexn in 1:20){
    
  
    if (indexn < 10){
      
      trueTree_form=sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/trees_dir/trees.000%s",alpha_str,alpha_str,numSites,indexn)
      
    }else if( indexn>=10 & indexn<100){
      
      trueTree_form=sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/trees_dir/trees.00%s",alpha_str,alpha_str,numSites,indexn)
      
      
    }else{
      
      trueTree_form=sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/trees_dir/trees.0%s",alpha_str,alpha_str,numSites,indexn)
      
    }
    
    
    scanned_trueTree_form = scan(file=trueTree_form,what=character(), n = -1, sep = "")
    
    sampletr = read.tree(text=scanned_trueTree_form)
    
    prob_form_0_1 = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/MO_Binary/Binary_alpha0%s_beta0%s_%s_result/scaled_all_binary_prob_matrix_all_0_1_out_matrix%s_mu%s.csv",alpha_str,alpha_str,numSites,indexn,m10_str)
    formMediumTrueBr = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/true_haplotypes_dir/scaled_finite_br_collapsed_true_hap%s_%s_selectedbr.RData",alpha_str,alpha_str,numSites,indexn,m10_str)
    
    mat_prob_form_0_1 = read.csv(prob_form_0_1)
    all_selected_First_Branch_sub = get(load(formMediumTrueBr))
    
  
  initial_obs_0_1 = as.matrix(mat_prob_form_0_1[,-c(1,2)])
  true_obs_0_1 = all_selected_First_Branch_sub

  true_obs_0_1_recoded= true_obs_0_1
  
  prob_list=data.frame(prob=numeric(0))
  cat_list=data.frame(branch=numeric(0))
  indexn_list=data.frame(indexn=numeric(0))
  for (i in 1:dim(initial_obs_0_1)[1]){
    
    print(c(indexn,i))
    indexn_list[i,]=indexn
    possible_br = which.max(initial_obs_0_1[i,])
    
    if (all(possible_br == true_obs_0_1_recoded[[i]])){prob_list[i,] = max(initial_obs_0_1[i,]) 
    cat_list[i,] = 1}
    else{prob_list[i,] = max(initial_obs_0_1[i,]) 
    cat_list[i,] = 0}
    
    
  }
  
  prob_cat_list[indexn,] = cbind(sum(cat_list),dim(mat_prob_form_0_1)[1],numSites,alpha,m10,missingindex=0,Genotype="Binary",lost_rate=0,lost_br=0,Method="MO")
  
  
  }
  
  colnames(prob_cat_list)=c("inferred location","true location","sites","error","m10","missingindex","Genotype","lost_rate","lost_br","MO")
  colnames(prob_cat_dat)=c("inferred location","true location","sites","error","m10","missingindex","Genotype","lost_rate","lost_br","MO")
  prob_cat_dat=rbind(prob_cat_dat,prob_cat_list)
}


prob_cat_list_0_1_out = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/MO_Binary/binary_location_accuracy_matrix_01_0_1_out_finite.csv')

write.csv(prob_cat_dat,file=prob_cat_list_0_1_out)


