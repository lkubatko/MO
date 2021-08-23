#ape and phangorn are dependent packages
library(ape)
library(phangorn)

source("/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario1_result/code_MO_binary_function.R")

#there are 16 settings with different error rates for complete data
parameter_setting = expand.grid(alpha=c(0.05,0.1,0.2,0.4),
                                beta =c(0.05,0.1,0.2,0.4))

#iteratively to run data in each setting
for(parameterIndex in 1:dim(parameter_setting)[1]){
  
  #get the parameters in each setting
  alpha =   parameter_setting[parameterIndex,1]#fpr
  beta = parameter_setting[parameterIndex,2]#fnr
  
  unit_theta = 10^(-7) #rate from 0 to 1
  unit_gamma = 10^(-9) #rate from 0 to 2
  unit_mu = 10 ^(-2) #rate from 1 to 2
  
  number_br = 18 #number of branches in a tree
  number_cell = 10 #number of cells in a tree
  
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
  
  binary_folder_form = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario1_result/MO_Binary_10missing/Binary_alpha0%s_beta0%s_result',alpha_str, beta_str)

  dir.create(binary_folder_form)
  
  #in each setting, there are 100 trees, iteratively to estimate mutation order based on each tree
  for (indexn in 1:1000){
    
    print(c(parameterIndex,indexn))
    
    # trueTree_form is the tree file name
    trueTree_form=sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/SimulateData/RandomTree/RandomTree/RandomTree_%s.tre', indexn)
    
    #scan and read the tree file into R 
    sampletr = read.tree(trueTree_form)
    
    #read the simulated data into R
    obs_form = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario1/SimulateData_New_10Missing/Binary_alpha0%s_beta0%s/binary_obs_0_1_tip_alpha_0%s_beta_0%s_matrix%s.csv', alpha_str, beta_str, alpha_str, beta_str,indexn)
    mat_obs_form_0_1 = read.csv(obs_form)
    
    #generate the normal,mutation, and observed data for inference
    normal_genotype_0_1 = rep(0,dim(mat_obs_form_0_1)[1])
    mutation_genotype_0_1 = rep(1,dim(mat_obs_form_0_1)[1])
    initial_obs_0_1 = data.frame(mat_obs_form_0_1[,-c(1,2,3)])
    #recode data
    initial_obs_0_1_recode = data.frame(matrix(,nrow=dim(initial_obs_0_1)[1],
                                               ncol=dim(initial_obs_0_1)[2]))
    
    
    for(rownum in 1:dim(initial_obs_0_1)[1]){
      for(colnum in 1:dim(initial_obs_0_1)[2])
        
        if(initial_obs_0_1[rownum,colnum] == "-"){
          
          initial_obs_0_1_recode[rownum,colnum] = 3 #recode missing values to 3
          
          
        }else if(initial_obs_0_1[rownum,colnum] == "0"){
          initial_obs_0_1_recode[rownum,colnum] = 0
        }else if(initial_obs_0_1[rownum,colnum] == "1"){
          initial_obs_0_1_recode[rownum,colnum] = 1
        }else{
          print("value does not exist")
        }
    }
    
    colnames(initial_obs_0_1_recode) = colnames(initial_obs_0_1)
    #initialize vector to store estimated probability
    binary_prob_matrix_all_0_1=c()
    #initialize vector to store MAP estimates
    selected_br=c()
    
    for (i in 1:dim(initial_obs_0_1_recode)[1]){
      
      #print(i)
      #compute probability that mutation i on each of the branches on the tree
      probs <- MO_binary(alpha,beta,unit_theta+unit_gamma,initial_obs_0_1_recode[i,],sampletr)
      
      #find the MAP
      selected_br[i]=which.max(probs)
      
      #save the ith result to the matrix
      binary_prob_matrix_all_0_1 = rbind(binary_prob_matrix_all_0_1,probs)
      
    }
      

    
    binary_prob_matrix_all_0_1_out = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario1_result/MO_Binary_10missing/Binary_alpha0%s_beta0%s_result/all_binary_prob_matrix_all_0_1_out_matrix%s.csv",alpha_str,beta_str,indexn)
    
    
    # combine the true mutation branch on column1 and the inferred mutation branch on column 3 for future comparison
    binary_prob_matrix_all_0_1_rownames=cbind(mat_obs_form_0_1[,c(2,3)],selected_br,binary_prob_matrix_all_0_1)
    # save inferred result for future comparison
    write.csv(binary_prob_matrix_all_0_1_rownames,file=binary_prob_matrix_all_0_1_out)
    
    
  }
  
}

