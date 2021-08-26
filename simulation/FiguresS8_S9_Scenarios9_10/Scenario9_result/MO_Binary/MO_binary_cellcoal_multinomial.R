library(ape)
library(phangorn)

source("/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/MO_Binary/code_MO_binary_function_multinomial.R")

parameter_setting = expand.grid(alpha=c(0.05,0.1,0.2,0.4),
                                sites =c(20,40,80))
set.seed(1000)

All_loc_dat=data.frame(matrix(ncol = 3, nrow = 0))

for(paraInd in 1:dim(parameter_setting)[1]){
  
  
  alpha =   parameter_setting[paraInd,1]
  beta = parameter_setting[paraInd,1]
  numSites= parameter_setting[paraInd,2]
  
  sequencing_error_model=matrix(c(1-alpha,alpha,
                                  beta,1-beta),nrow=2,byrow = TRUE)
  print(sequencing_error_model)
  unit_theta = 10^(-7)
  unit_gamma = 10^(-9)
  unit_mu = 10 ^(-2)
  number_br = 18
  number_cell = 10
  
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
  
  
  binary_folder_form_result = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/MO_Binary/Binary_alpha0%s_beta0%s_%s_result',alpha_str, beta_str,numSites)
  dir.create(binary_folder_form_result)
  
  location_acc = data.frame(matrix(ncol = 3, nrow = 100))
  
  for (indexn in 1:100){
    
    print(indexn)
    
    
    form = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/SimulateData/RandomTree/RandomTree_%s.tre', indexn)
    
    sampletr=read.tree(form)
    
    #binary_folder_form_result_br = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/true_haplotypes_dir/br_collapsed_true_hap%s.csv",alpha_str,alpha_str,numSites,indexn)
    obs_binary_folder_form_result_br = sprintf("/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/SimulateData/Ternary_alpha0%s_beta0%s_sites%s/ternary_obs_0_1_tip_alpha_0%s_beta_0%s_matrix%s.csv", alpha_str, beta_str, numSites,alpha_str, beta_str,indexn)
    
    mat_obs_form_0_1 = read.csv(obs_binary_folder_form_result_br)
    
    mat_obs_form_0_1_sub=na.omit(mat_obs_form_0_1)
    
    normal_genotype_0_1 = rep(0,dim(mat_obs_form_0_1_sub)[1])
    mutation_genotype_0_1 = rep(1,dim(mat_obs_form_0_1_sub)[1])
    initial_obs_0_1 = data.frame(mat_obs_form_0_1_sub[,-c(1,2,3)])
    
    initial_obs_0_1_recode = data.frame(matrix(ncol = dim(initial_obs_0_1)[2], nrow = dim(initial_obs_0_1)[1]))
    
    for ( nr in 1:dim(initial_obs_0_1)[1]){
      for ( nc in 1:dim(initial_obs_0_1)[2]){
        if (initial_obs_0_1[nr,nc] == "-"){initial_obs_0_1_recode[nr,nc] ="3"}
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
      
      
      
      generate_prob_br <- MO_binary_multinomial(alpha,beta,initial_obs_0_1_recode[i,],sampletr)
        
        
      generate_prob_br_all_single <- c(generate_prob_br,rep(0,(number_br+1)-length(generate_prob_br)))
      
      binary_prob_matrix_all_0_1 = rbind(binary_prob_matrix_all_0_1,generate_prob_br_all_single)
      
    }
      

    
    binary_prob_matrix_all_0_1_out = sprintf("/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/MO_Binary/Binary_alpha0%s_beta0%s_%s_result/all_binary_prob_matrix_all_0_1_2_out_matrix%s_multinomial.csv",alpha_str,alpha_str,numSites,indexn)
    
    selected_br=c()
    for(i in 1:dim(binary_prob_matrix_all_0_1)[1]){
      
      selected_br[i]=which.max(binary_prob_matrix_all_0_1[i,])
    }
    
    binary_prob_matrix_all_0_1_rownames=cbind(mat_obs_form_0_1_sub[,c(2,3)],selected_br,binary_prob_matrix_all_0_1)
    
    location_acc[indexn,]=c(sum(binary_prob_matrix_all_0_1_rownames[,1]==binary_prob_matrix_all_0_1_rownames[,3]),dim(binary_prob_matrix_all_0_1_rownames)[1],alpha)
    write.csv(binary_prob_matrix_all_0_1_rownames,file=binary_prob_matrix_all_0_1_out)
    
    
  }
  All_loc_dat=rbind(All_loc_dat,location_acc)
}

write.csv(All_loc_dat,file="/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/MO_Binary/All_MO_binary_location_accuracy_multinomial.csv")