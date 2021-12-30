library(ape)
library(phangorn)

source("/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/MO/code/code_MO_ternary_function.R")

parameter_setting = expand.grid(alpha=c(0.05,0.1,0.2,0.4),
                                sites =c(20,40,80))
set.seed(20000)
All_loc_dat=data.frame(matrix(ncol = 3, nrow = 0))
for(paraInd in 1:dim(parameter_setting)[1]){
  
  alpha =   parameter_setting[paraInd,1]
  beta = parameter_setting[paraInd,1]
  numSites= parameter_setting[paraInd,2]
  

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
unit_theta = 1
unit_gamma = 0
unit_mu = 0
number_br = 100
number_cell = 51

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


ternary_folder_form_result = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/MO/simulation/FiguresS4_S5_Scenarios5_6/Scenario6_result_rmd/MO_Ternary/Ternary_alpha0%s_beta0%s_%s_result',alpha_str, beta_str,numSites)
dir.create(ternary_folder_form_result)



for (indexn in 1:100){

  print(indexn)
  
  ternary_prob_matrix_all_0_1_2_out = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/MO/simulation/FiguresS4_S5_Scenarios5_6/Scenario6_result_rmd/MO_Ternary/Ternary_alpha0%s_beta0%s_%s_result/all_ternary_prob_matrix_all_0_1_2_out_matrix%s.csv",alpha_str,alpha_str,numSites,indexn)
  
  if(file.exists(ternary_prob_matrix_all_0_1_2_out)==FALSE){
    
  if (indexn < 10){trueTree_form=sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/trees_dir/trees.000%s",alpha_str,alpha_str,numSites,indexn)
  true_dat_form = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/true_haplotypes_dir/true_hap.000%s.csv",alpha_str,alpha_str,numSites,indexn)
  obs_dat_form = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/snv_hap.000%s.csv",alpha_str,alpha_str,numSites,indexn)
  
  }else if( indexn>=10 & indexn<100){
    trueTree_form=sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/trees_dir/trees.00%s",alpha_str,alpha_str,numSites,indexn)
    true_dat_form = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/true_haplotypes_dir/true_hap.00%s.csv",alpha_str,alpha_str,numSites,indexn)
    obs_dat_form = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/snv_hap.00%s.csv",alpha_str,alpha_str,numSites,indexn)
    
    
  }else{
    
    trueTree_form=sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/trees_dir/trees.0%s",alpha_str,alpha_str,numSites,indexn)
    true_dat_form = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/true_haplotypes_dir/true_hap.0%s.csv",alpha_str,alpha_str,numSites,indexn)
    obs_dat_form = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/snv_hap.0%s.csv",alpha_str,alpha_str,numSites,indexn)
    
  }
  
  
  scanned_trueTree_form = scan(file=trueTree_form,what=character(), n = -1, sep = "")
  
  sampletr = read.tree(text=scanned_trueTree_form)
  
  
  
 
  obs_ternary_folder_form_result_br = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/br_collapsed_snv_hap%s.csv",alpha_str,alpha_str,numSites,indexn)
  mat_obs_form_0_1_2 = read.csv(obs_ternary_folder_form_result_br)

  if (dim(mat_obs_form_0_1_2)[1]>0){
    
    
    
    normal_genotype_0_1_2 = rep(0,dim(mat_obs_form_0_1_2)[1])
    mutation_genotype_0_1_2 = rep(2,dim(mat_obs_form_0_1_2)[1])
    initial_obs_0_1 = data.frame(mat_obs_form_0_1_2[,-c(1,2,3)])
    
    initial_obs_0_1_recode = data.frame(matrix(ncol = dim(initial_obs_0_1)[2], nrow = dim(initial_obs_0_1)[1]))
    
    for ( nr in 1:dim(initial_obs_0_1)[1]){
      for ( nc in 1:dim(initial_obs_0_1)[2]){
        if (initial_obs_0_1[nr,nc] == "-"){initial_obs_0_1_recode[nr,nc] ="-"}
        if (initial_obs_0_1[nr,nc] == "0"){initial_obs_0_1_recode[nr,nc] ="0"}
        if (initial_obs_0_1[nr,nc] == "1"){initial_obs_0_1_recode[nr,nc] ="1"}
        if (initial_obs_0_1[nr,nc] == "2"){initial_obs_0_1_recode[nr,nc] ="2"}
      }
    }
    
    colnames(initial_obs_0_1_recode)=colnames(initial_obs_0_1)
    
    
    
    #initialize vector to store estimated probability.Initialize vector to store MAP estimates
    ternary_prob_matrix_all_0_1=c()
    selected_br=c()
    
    for (i in 1:dim(initial_obs_0_1_recode)[1]){
      
      #compute probability that mutation i on each of the branches on the tree
      rd_unit_theta =  rgamma(n = 3, shape = 100, scale = 0.01*unit_theta)
      rd_unit_gamma = rgamma(n = 3, shape = 100, scale = 0.01*unit_gamma)
      rd_unit_mu = rgamma(n = 3, shape = 100, scale = 0.01*unit_mu)
      
      generate_prob_br_all_dat=data.frame(matrix(NA, nrow = number_br, ncol = 3))
      
      for (j in 1:3){
        
        generate_prob_br <- MO_ternary(alpha,beta,rd_unit_theta[j],rd_unit_gamma[j],rd_unit_mu[j],initial_obs_0_1_recode[i,],sampletr)
        generate_prob_br_all_single <- c(generate_prob_br,rep(0,number_br-length(generate_prob_br)))
        generate_prob_br_all_dat[,j] = generate_prob_br_all_single}
      
      probs=rowMeans(generate_prob_br_all_dat, na.rm = FALSE, dims = 1)
      #find the MAP
      selected_br[i]=which.max(probs)
      #save the ith result to the matrix
      ternary_prob_matrix_all_0_1 = rbind(ternary_prob_matrix_all_0_1,probs)
      
    }
    
    
    
  }
  
  
  ternary_prob_matrix_all_0_1_rownames=cbind(mat_obs_form_0_1_2[,c(2,3)],selected_br,ternary_prob_matrix_all_0_1)
  
  ternary_prob_matrix_all_0_1_2_out = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/MO/simulation/FiguresS4_S5_Scenarios5_6/Scenario6_result_rmd/MO_Ternary/Ternary_alpha0%s_beta0%s_%s_result/all_ternary_prob_matrix_all_0_1_2_out_matrix%s.csv",alpha_str,alpha_str,numSites,indexn)
  write.csv(ternary_prob_matrix_all_0_1_rownames,file=ternary_prob_matrix_all_0_1_2_out)
  
  }
  
  
}}




