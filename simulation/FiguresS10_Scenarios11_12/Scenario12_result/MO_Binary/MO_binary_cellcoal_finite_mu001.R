library(ape)
library(phangorn)

source("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/MO_Binary/MO_binary_function.R")

parameter_setting = expand.grid(m10=c(0.01),
                                alpha=c(0.05,0.1,0.2,0.4),
                                sites =c(20,40,80))
set.seed(1000)

All_loc_dat=data.frame(matrix(ncol = 3, nrow = 0))

for(paraInd in 1:dim(parameter_setting)[1]){
  
  m10 = parameter_setting[paraInd,1] 
  alpha =   parameter_setting[paraInd,2]
  beta = parameter_setting[paraInd,2]
  numSites= parameter_setting[paraInd,3]
  
  sequencing_error_model=matrix(c(1-alpha,alpha,
                                  beta,1-beta),nrow=2,byrow = TRUE)
  print(sequencing_error_model)
  unit_theta = 1
  unit_gamma = 0#10^(-14)
  unit_mu = 0#10 ^(-2)
  number_br = 100
  number_cell = 51
  
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
  
  
  #binary_folder_form_result = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/MO_Binary/Binary_alpha0%s_beta0%s_%s_result',alpha_str, beta_str,numSites)
  #dir.create(binary_folder_form_result)
  
  location_acc = data.frame(matrix(ncol = 3, nrow = 100))
  
  for (indexn in 1:20){
    print(indexn)

    
    if (indexn < 10){trueTree_form=sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/trees_dir/trees.000%s",alpha_str,alpha_str,numSites,indexn)
    
    }else if( indexn>=10 & indexn<100){
      trueTree_form=sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/trees_dir/trees.00%s",alpha_str,alpha_str,numSites,indexn)
      
      
    }else{
      
      trueTree_form=sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/trees_dir/trees.0%s",alpha_str,alpha_str,numSites,indexn)
      
    }
    
  
    scanned_trueTree_form = scan(file=trueTree_form,what=character(), n = -1, sep = "")
    
    trueTree = read.tree(text=scanned_trueTree_form)
    
    sampletr_original= trueTree
    sampletr = sampletr_original
    sampletr$edge.length= sampletr_original$edge.length*1000
    #binary_folder_form_result_br = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/true_haplotypes_dir/br_collapsed_true_hap%s.csv",alpha_str,alpha_str,numSites,indexn)
    obs_binary_folder_form_result_br = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/scaled_finite_br_collapsed_snv_hap%s_%s.csv",alpha_str,alpha_str,numSites,indexn, m10_str)
    
    mat_obs_form_0_1 = read.csv(obs_binary_folder_form_result_br)
    
    mat_obs_form_0_1_sub=na.omit(mat_obs_form_0_1)
    
    normal_genotype_0_1 = rep(0,dim(mat_obs_form_0_1_sub)[1])
    mutation_genotype_0_1 = rep(1,dim(mat_obs_form_0_1_sub)[1])
    initial_obs_0_1 = data.frame(mat_obs_form_0_1_sub)
    
    initial_obs_0_1_recode = data.frame(matrix(ncol = dim(initial_obs_0_1)[2], nrow = dim(initial_obs_0_1)[1]))
    
    for ( nr in 1:dim(initial_obs_0_1)[1]){
      for ( nc in 1:dim(initial_obs_0_1)[2]){
        if (initial_obs_0_1[nr,nc] == "-"){initial_obs_0_1_recode[nr,nc] ="-"}
        if (initial_obs_0_1[nr,nc] == "0"){initial_obs_0_1_recode[nr,nc] ="0"}
        if (initial_obs_0_1[nr,nc] == "1"){initial_obs_0_1_recode[nr,nc] ="1"}
        if (initial_obs_0_1[nr,nc] == "2"){initial_obs_0_1_recode[nr,nc] ="1"}
      }
    }
    
    colnames(initial_obs_0_1_recode)=colnames(initial_obs_0_1)
    
    binary_prob_matrix_all_0_1=c()
    
    
    for (i in 1:dim(initial_obs_0_1_recode)[1]){
      
      print(i)
      
      #rd_unit_theta <- rbeta(10, (10^7)*unit_theta, (10^7)*(1-unit_theta))
      #rd_unit_gamma <- rbeta(10, (10^14)*unit_gamma, (10^14)*(1-unit_gamma))
      
      #rd_unit_theta =  rgamma(n = 3, shape = 100, scale = 0.01*unit_theta)
      #rd_unit_gamma = rgamma(3, shape = 100, scale = 0.01*unit_gamma)
      
      
      
      generate_prob_br <- generate_prob_binary(alpha,beta,unit_theta+unit_gamma,number_br,number_cell,
                                                          normal_genotype_0_1[i],mutation_genotype_0_1[i],initial_obs_0_1_recode[i,],sampletr)
      
      
      generate_prob_br_all_single <- c(generate_prob_br,rep(0,(number_br+1)-length(generate_prob_br)))
      
      binary_prob_matrix_all_0_1 = rbind(binary_prob_matrix_all_0_1,generate_prob_br_all_single)
      
    }
      

    
    binary_prob_matrix_all_0_1_out = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/MO_Binary/Binary_alpha0%s_beta0%s_%s_result/scaled_all_binary_prob_matrix_all_0_1_out_matrix%s_mu%s.csv",alpha_str,alpha_str,numSites,indexn,m10_str)
    
    selected_br=c()
    for(i in 1:dim(binary_prob_matrix_all_0_1)[1]){
      
      selected_br[i]=which.max(binary_prob_matrix_all_0_1[i,])
    }
    
    binary_prob_matrix_all_0_1_rownames=cbind(selected_br,binary_prob_matrix_all_0_1)
    
    location_acc[indexn,]=c(sum(binary_prob_matrix_all_0_1_rownames[,1]==binary_prob_matrix_all_0_1_rownames[,3]),dim(binary_prob_matrix_all_0_1_rownames)[1],alpha)
    write.csv(binary_prob_matrix_all_0_1_rownames,file=binary_prob_matrix_all_0_1_out)
    
    
  }
  #All_loc_dat=rbind(All_loc_dat,location_acc)
}

#write.csv(All_loc_dat,file="/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/MO_Binary/All_MO_binary_location_accuracy.csv")
