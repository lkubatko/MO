library(ape)
library(phangorn)

source("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/MO_Ternary/MO_ternary_function.R")

parameter_setting = expand.grid(m10=c(0.1),
                                alpha=c(0.05,0.1,0.2,0.4),
                                sites =c(20,40,80))
set.seed(1000)

All_loc_dat=data.frame(matrix(ncol = 3, nrow = 0))

for(paraInd in 6:dim(parameter_setting)[1]){
  
  m10 = parameter_setting[paraInd,1] 
  alpha =   parameter_setting[paraInd,2]
  beta = parameter_setting[paraInd,2]
  numSites= parameter_setting[paraInd,3]
  
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
  
  
  #ternary_folder_form_result = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/MO_Ternary/Ternary_alpha0%s_beta0%s_%s_result',alpha_str, beta_str,numSites)
  #dir.create(ternary_folder_form_result)
  
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
    #ternary_folder_form_result_br = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/true_haplotypes_dir/br_collapsed_true_hap%s.csv",alpha_str,alpha_str,numSites,indexn)
    obs_ternary_folder_form_result_br = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/scaled_finite_br_collapsed_snv_hap%s_%s.csv",alpha_str,alpha_str,numSites,indexn, m10_str)
    
    mat_obs_form_0_1_2 = read.csv(obs_ternary_folder_form_result_br)
    
    ternary_prob_matrix_all_0_1_all=c()
    ternary_prob_matrix_01_0_1_all=c()
    ternary_prob_matrix_02_0_1_all=c()
    ternary_prob_matrix_012_0_1_all=c()
    
    
    
    
    if (dim(mat_obs_form_0_1_2)[1]>0){
      
      ternary_prob_matrix_all_0_1=c()
      ternary_prob_matrix_01_0_1=c()
      ternary_prob_matrix_02_0_1=c()
      ternary_prob_matrix_012_0_1=c()
      
      normal_genotype_0_1_2 = rep(0,dim(mat_obs_form_0_1_2)[1])
      mutation_genotype_0_1_2 = rep(2,dim(mat_obs_form_0_1_2)[1])
      initial_obs_0_1 = data.frame(mat_obs_form_0_1_2)
      
      
      for (i in 1:dim(initial_obs_0_1)[1]){
        
        print(i)
        
        #rd_unit_theta <- rbeta(10, (10^7)*unit_theta, (10^7)*(1-unit_theta))
        #rd_unit_gamma <- rbeta(10, (10^14)*unit_gamma, (10^14)*(1-unit_gamma))
        #rd_unit_mu <- rbeta(10, 100*unit_mu, 100*(1-unit_mu))
        
        
        
        
        generate_prob_br <- generate_prob(sequencing_error_model,unit_theta,unit_gamma,unit_mu,number_br,number_cell,
                                          normal_genotype_0_1_2[i],mutation_genotype_0_1_2[i],initial_obs_0_1[i,],sampletr)
        generate_prob_br_0_1_single <- c(generate_prob_br[,1],rep(0,number_br-dim(generate_prob_br)[1]))
        generate_prob_br_0_2_single <- c(generate_prob_br[,2],rep(0,number_br-dim(generate_prob_br)[1]))
        generate_prob_br_0_1_2_single <- c(generate_prob_br[,3],rep(0,number_br-dim(generate_prob_br)[1]))
        generate_prob_br_all_single <- c(generate_prob_br[,4],rep(0,number_br-dim(generate_prob_br)[1]))
        
        
        generate_prob_br_0_1=generate_prob_br_0_1_single
        generate_prob_br_0_2=generate_prob_br_0_2_single
        generate_prob_br_0_1_2=generate_prob_br_0_1_2_single
        generate_prob_br_all=generate_prob_br_all_single
        
        
        ternary_prob_matrix_all_0_1 = rbind(ternary_prob_matrix_all_0_1,generate_prob_br_all)
        ternary_prob_matrix_01_0_1 = rbind(ternary_prob_matrix_01_0_1,generate_prob_br_0_1)
        ternary_prob_matrix_02_0_1 = rbind(ternary_prob_matrix_02_0_1,generate_prob_br_0_2)
        ternary_prob_matrix_012_0_1 = rbind(ternary_prob_matrix_012_0_1,generate_prob_br_0_1_2)
        
        
      }
      
      
      ternary_prob_matrix_all_0_1_all=rbind(ternary_prob_matrix_all_0_1_all,data.frame(ternary_prob_matrix_all_0_1))
      
      
      
    }
    
    
    
    
    ternary_prob_matrix_all_0_1_2_out = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/MO_Ternary/Ternary_alpha0%s_beta0%s_%s_result/scaled_all_ternary_prob_matrix_all_0_1_out_matrix%s_mu%s.csv",alpha_str,alpha_str,numSites,indexn,m10_str)
    
    
    
    selected_br=c()
    for(i in 1:dim(ternary_prob_matrix_all_0_1)[1]){
      
      selected_br[i]=which.max(ternary_prob_matrix_all_0_1[i,])
    }
    
    ternary_prob_matrix_all_0_1_rownames=cbind(selected_br,ternary_prob_matrix_all_0_1)
    
    location_acc[indexn,]=c(sum(ternary_prob_matrix_all_0_1_rownames[,1]==ternary_prob_matrix_all_0_1_rownames[,3]),dim(ternary_prob_matrix_all_0_1_rownames)[1],alpha)
    write.csv(ternary_prob_matrix_all_0_1_rownames,file=ternary_prob_matrix_all_0_1_2_out)
    
    
  }
  #All_loc_dat=rbind(All_loc_dat,location_acc)
}

#write.csv(All_loc_dat,file="/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_50tips/MO_Ternary/All_MO_ternary_location_accuracy.csv")
