#####################################################
#####################################################
#BEGIN:REQUIRED PACKAGES
#####################################################
#####################################################
library(ape)
library(plyr)
require(xlsx)

#####################################################
#####################################################
#END:REQUIRED PACKAGES
#####################################################
#####################################################

#####################################################

#####################################################
#####################################################
#BEGIN:SIMULATE
#####################################################
#####################################################



# simulate_sample add errors to tips for mutation on each br of a tree 
simulate_sample <- function(sequencing_error_model,number_br,
                            number_cell,ts){
  
  #true_0_1_tip_matrix stores the true value if a mutation 0 to 1 occurs on a branch
  true_0_1_tip_matrix = matrix(rep(0,(length(ts$tip.label)+2)*length(ts$edge.length)),nrow=length(ts$edge.length),ncol = (length(ts$tip.label)+2))
  
  #true_0_2_tip_matrix stores the true value if a mutation 0 to 2 occurs on a branch
  true_0_2_tip_matrix = matrix(rep(0,(length(ts$tip.label)+2)*length(ts$edge.length)),nrow=length(ts$edge.length),ncol = (length(ts$tip.label)+2))
  
  #true_0_1_2_tip_matrix stores the true value if a mutation 0 to 1 then 1 to2 occurs on a branch and its descendent
  true_0_1_2_tip_matrix = list()
  
  # the observed matrix correspond to the above three matrix
  obs_0_1_tip_matrix = matrix(rep(0,(length(ts$tip.label)+2)*length(ts$edge.length)),nrow=length(ts$edge.length),ncol = (length(ts$tip.label)+2))
  obs_0_2_tip_matrix = matrix(rep(0,(length(ts$tip.label)+2)*length(ts$edge.length)),nrow=length(ts$edge.length),ncol = (length(ts$tip.label)+2))
  obs_0_1_2_tip_matrix = list()
  
  
  ####################################################################################################
  ####################################################################################################
  #extract the tree, if mutation is on one branch, then the corresponding tips will have mutation
  ####################################################################################################
  ####################################################################################################
  
  left_right <- function(edge,parent){
    child = c()
    for (i in 1:nrow(edge)) {
      if (edge[i,1] == parent) {
        child = c(child,edge[i,2])
      }
    }
    return(child)
  }
  
  build_tree <- function(edge,branch){
    child_node = left_right(edge,branch[length(branch)])
    new_branch=matrix(c(branch,child_node[1],branch,child_node[2]),nrow=2,byrow = TRUE)
    return(new_branch)
  }
  
  #####################################modify begin################################
  # find node parent
  find_ancestor <- function(edge,node){
    parent = 0
    for (i in 1:nrow(edge)) {
      if (edge[i,2] == node) {
        parent = edge[i,1]
      }
    }
    return(parent)
  }
  
  # get all unique nodes in the tree
  get_all_nodes <- function(edge)
  {
    all_nodes = integer(length(edge))
    for (i in 1:nrow(edge))
    {
      all_nodes[(i-1)*2+1] = edge[i,1]
      all_nodes[(i-1)*2+2] = edge[i,2]
    }
    all_nodes = unique(all_nodes)
    return(all_nodes)
  }
  
  # find root node
  find_root <- function(edge)
  {
    all_nodes = get_all_nodes(edge)
    
    for (i in 1:length(all_nodes))
    {
      parent = find_ancestor(edge, all_nodes[i])
      if (parent == 0)
      {
        root_node = all_nodes[i]
        break
      }
    }
  }
  
  # find two child branches and nodes if they exist. Otherwise all zeros matrix output
  find_child_branches_and_nodes <- function(edge, parent_node){
    child_branches_and_nodes = matrix(0, 2, 2)
    child_id = 1
    # first row are two nodes, second row are two branches
    for (i in 1:nrow(edge))
    {
      if (edge[i,1] == parent_node) {
        child_branches_and_nodes[1,child_id] = edge[i,2]
        child_branches_and_nodes[2,child_id] = i
        child_id = child_id + 1
      }
    }
    
    return(child_branches_and_nodes)
  }
  
  # find all child branch for current branch
  find_child_branches <- function(edge, current_edge, child_branches)
  {
    id = length(child_branches)
    right_node = edge[current_edge, 2]
    
    child_branches_and_nodes = find_child_branches_and_nodes(edge, right_node)
    
    if (child_branches_and_nodes[1,1] != 0)
    {
      # if not leaf node
      left_node = child_branches_and_nodes[1,1]
      right_node = child_branches_and_nodes[1,2]
      left_branch = child_branches_and_nodes[2,1]
      right_branch = child_branches_and_nodes[2,2]
      
      id = id + 1
      child_branches[id] = left_branch
      id = id + 1
      child_branches[id] = right_branch
      
      child_branches = find_child_branches(edge, left_branch, child_branches)
      child_branches = find_child_branches(edge, right_branch, child_branches)
      
      return(child_branches)
      
    }
    else
    {
      return(child_branches)
    }
    
  }
  
  # find all child branch for all branches
  find_all_child_branches <- function(edge){
    # get root node
    root_node = find_root(edge)
    
    all_child_branches = rep(list(list()), nrow(edge)) 
    
    for (i in 1:nrow(edge))
    {
      current_edge = i
      # iterative find all its child branches
      child_branches = integer(0)
      all_child_branches[[i]] = find_child_branches(edge, current_edge, child_branches)
      
    }
    
    return(all_child_branches)
  }
  
  find_all_tip_nodes <- function(edge)
  {
    all_parent_nodes = numeric()
    for (i in 1:nrow(edge))
    {
      all_parent_nodes = c(all_parent_nodes, edge[i, 1])
    }
    all_parent_nodes = unique(all_parent_nodes)
    
    all_nodes = get_all_nodes(edge)
    all_tip_nodes = numeric()
    for (i in 1:length(all_nodes))
    {
      if (!is.element(all_nodes[i], all_parent_nodes))
        all_tip_nodes = c(all_tip_nodes, all_nodes[i])
    }
    
    return(all_tip_nodes)
  }
  
  # find tip nodes under one edge
  find_tip_nodes_of_edge <- function(edge)
  {
    all_tip_nodes = find_all_tip_nodes(edge)
    all_child_branches = find_all_child_branches(edge)
    
    all_branch_tips = rep(list(list()), nrow(edge)) 
    
    for (i in 1:nrow(edge)) {
      child_branches = all_child_branches[[i]]
      
      tip_nodes = numeric()
      if (length(child_branches) > 0)
      {
        for (j in 1:length(child_branches)) {
          child_node = edge[child_branches[j], 2]
          if (is.element(child_node, all_tip_nodes))
          {
            tip_nodes = c(tip_nodes, child_node)
          }
        }
      }
      
      child_node = edge[i, 2]
      if (is.element(child_node, all_tip_nodes))
      {
        tip_nodes = c(tip_nodes, child_node)
      }
      
      all_branch_tips[[i]] = tip_nodes
    }
    
    return(all_branch_tips)
  }
  
  ###################################################################################
  ###################################################################################
  #find the joint prob of observation and branch_i in the subtree
  ###################################################################################
  ###################################################################################
  
  
  #find_all_possible_mutation_matrix returns the true mutation status at tips if mutation occurred on branch_i
  find_all_possible_mutation_matrix <- function(subtree, branch_i)
  {
    num_rows = nrow(subtree$edge)
    num_cols = length(subtree$tip.label)
    
    # build the branch tree sturcture from each tip to the root
    branch_trees = rep( list(list()),num_cols ) 
    num_parent = 0
    for (tip_i in 1:num_cols) {
      branch_trees[[tip_i]] = tip_i
      parent = find_ancestor(subtree$edge,tip_i)
      branch_trees[[tip_i]][num_parent+2] = parent
      num_parent=num_parent+1
      
      while (parent != num_cols+1) {
        tip_node = parent
        parent = find_ancestor(subtree$edge,tip_node)
        branch_trees[[tip_i]][num_parent+2] = parent
        num_parent=num_parent+1
        
      }
      num_parent = 0
    }
    
    # loop over all the branches, and find the possible final stage 
    # if the event occurs in that branch
    possible_true_genotype_with_1 = matrix(rep(0,num_rows*num_cols),nrow=num_rows,ncol = num_cols)
    for (j in 1:num_rows) {
      branch_edge = subtree$edge[j,]
      if (branch_edge[2] <= num_cols) {
        possible_true_genotype_with_1[j,branch_edge[2]] = 1
      }else {
        for (i in 1:num_cols) {
          list_branch = branch_trees[[i]]
          if (is.na(match(branch_edge[2],list_branch)) == FALSE) {
            possible_true_genotype_with_1[j,i] = 1
            colnames(possible_true_genotype_with_1)=subtree$tip.label
          }
        }
      }
    }
    
    descendant_branches=find_all_child_branches(subtree$edge)
    
    num_of_cases = num_rows-1
    possible_true_genotype_with_only_1 = matrix(rep(0,num_cols+2),nrow=1,ncol = num_cols+2)
    possible_true_genotype_with_only_2 = matrix(rep(0,num_cols+2),nrow=1,ncol = num_cols+2)
    possible_true_genotype_with_2 = matrix(rep(0,num_of_cases*(num_cols+2)),nrow=num_of_cases,ncol = num_cols+2)
    
    colnames(possible_true_genotype_with_only_1)= c("First_branch", "Second_branch", subtree$tip.label)
    colnames(possible_true_genotype_with_only_2)= c("First_branch", "Second_branch", subtree$tip.label)
    colnames(possible_true_genotype_with_2)= c("First_branch", "Second_branch", subtree$tip.label)
    
    # a matrix is created, where first two columns are the places for the mutation occuring.
    # if it's NA, it means that no mutation occurs to stand for situation like 0, 1, 2
    # if only this branch has one mutation
    possible_true_genotype_with_only_1[1,1] = branch_i
    possible_true_genotype_with_only_1[1,1:num_cols+2] = possible_true_genotype_with_1[branch_i, ]
    
    possible_true_genotype_with_only_2[1,1] = branch_i
    possible_true_genotype_with_only_2[1,2] = branch_i
    possible_true_genotype_with_only_2[1,1:num_cols+2] = 2*possible_true_genotype_with_1[branch_i, ]
    
    id_row = 1
    # if this branch has one mutation, and other branch has another
    for (branch_j in 1:num_rows) {
      if (branch_j == branch_i)
      {
        next
      }
      
      possible_true_genotype_with_2[id_row,1] = branch_i
      possible_true_genotype_with_2[id_row,2] = branch_j
      
      possible_true_genotype_with_2[id_row,1:num_cols+2] = possible_true_genotype_with_1[branch_i, ] + possible_true_genotype_with_1[branch_j, ]
      id_row = id_row+1
    }
    
    possible_true_genotype_with_2_sub = possible_true_genotype_with_2[ possible_true_genotype_with_2[,2] %in% descendant_branches[[branch_i]], ]
    
    return(list(possible_true_genotype_with_only_1, possible_true_genotype_with_only_2, possible_true_genotype_with_2_sub))
  }
  
  
  for (i in 1:length(ts$edge.length)){
    
    true_possible_tip_all=find_all_possible_mutation_matrix(ts,i)
    true_0_1_tip_matrix[i,]=(true_possible_tip_all[[1]])
    true_0_2_tip_matrix[i,]=(true_possible_tip_all[[2]])
    
    true_0_1_2_tip_matrix[[i]]=(true_possible_tip_all[[3]])
    
    obs_0_1_2_tip_matrix[[i]]= (true_possible_tip_all[[3]])
  }
  
  colnames(true_0_1_tip_matrix) = c("First_branch","Second_branch",ts$tip.label)
  colnames(true_0_2_tip_matrix) = c("First_branch","Second_branch",ts$tip.label)
  
  colnames(obs_0_1_tip_matrix) = c("First_branch","Second_branch",ts$tip.label)
  colnames(obs_0_2_tip_matrix) = c("First_branch","Second_branch",ts$tip.label)
  
  # sample errors at the tips when true value is 0
  sampleDist_0 = function(n) { sample(x = c(0,1,2), n, replace = T, prob = sequencing_error_model[1,])}
  
  # sample errors at the tips when true value is 1
  sampleDist_1 = function(n) { sample(x = c(0,1,2), n, replace = T, prob = sequencing_error_model[2,])}
  
  # sample errors at the tips when true value is 2
  sampleDist_2 = function(n) { sample(x = c(0,1,2), n, replace = T, prob = sequencing_error_model[3,])}
  
  ########################################################################
  # add true mutation branch in the first two columns of observed matrix
  obs_0_1_tip_matrix[,c(1,2)] = true_0_1_tip_matrix[,c(1,2)]
  # add errors independently to each tip if mutation occurs on any branch
  for (i in 1:length(ts$edge.length)){
    for (j in 1:length(ts$tip.label))
      
      if(true_0_1_tip_matrix[i,j+2]==1){
        
        obs_0_1_tip_matrix[i,j+2]=sampleDist_1(1)
        
        }else if(true_0_1_tip_matrix[i,j+2]==0){
          
        obs_0_1_tip_matrix[i,j+2]=sampleDist_0(1)
        }else{
        
        obs_0_1_tip_matrix[i,j+2]=true_0_1_tip_matrix[i,j+2]
        }
   }
  
  
  ########################################################################
  # add true mutation branch in the first two columns of observed matrix
  obs_0_2_tip_matrix[,c(1,2)] = true_0_2_tip_matrix[,c(1,2)]
  # add errors independently to each tip if mutation occurs on any branch
  for (i in 1:length(ts$edge.length)){
    for (j in 1:length(ts$tip.label))
      
      if(true_0_2_tip_matrix[i,j+2]==2){ 
        
        obs_0_2_tip_matrix[i,j+2]=sampleDist_2(1)
        
      }else if(true_0_2_tip_matrix[i,j+2]==0){
          
        obs_0_2_tip_matrix[i,j+2]=sampleDist_0(1)
      }else{
        
        obs_0_2_tip_matrix[i,j+2]=true_0_2_tip_matrix[i,j+2]
      }
  }
  
  
  
  ########################################################################
  for (k in 1:length(true_0_1_2_tip_matrix)){
    
    # add error if mutation can initiate from the br
    if (dim(true_0_1_2_tip_matrix[[k]])[1] >0){
      
      for (i in 1:dim(true_0_1_2_tip_matrix[[k]])[1]){
      
       for (j in 1:length(ts$tip.label)){
        
        if(true_0_1_2_tip_matrix[[k]][i,j+2]==2){
          
          obs_0_1_2_tip_matrix[[k]][i,j+2]=sampleDist_2(1)}else if(true_0_1_2_tip_matrix[[k]][i,j+2]==1){
            
          obs_0_1_2_tip_matrix[[k]][i,j+2]=sampleDist_1(1)}else if(true_0_1_2_tip_matrix[[k]][i,j+2]==0){
            
          obs_0_1_2_tip_matrix[[k]][i,j+2]=sampleDist_0(1)}else{
            
          obs_0_1_2_tip_matrix[[k]][i,j+2]=true_0_1_2_tip_matrix[[k]][i,j+2]}
        
      }
     }
    }else{
      
      obs_0_1_2_tip_matrix[[k]]=true_0_1_2_tip_matrix[[k]]
      
      }
  }
  
  
  
  return(list(true_0_1_tip_matrix=true_0_1_tip_matrix,
              true_0_2_tip_matrix=true_0_2_tip_matrix,
              true_0_1_2_tip_matrix=true_0_1_2_tip_matrix,
              obs_0_1_tip_matrix = obs_0_1_tip_matrix,
              obs_0_2_tip_matrix = obs_0_2_tip_matrix,
              obs_0_1_2_tip_matrix = obs_0_1_2_tip_matrix))
  
  
}


set.seed(1000)

parameter_setting = expand.grid(alpha=c(0.05,0.1,0.2,0.4),
                                beta =c(0.05,0.1,0.2,0.4))

for(parameter in 1:dim(parameter_setting)[1]){
  
alpha =   parameter_setting[parameter,1]
beta = parameter_setting[parameter,2]

if (alpha < 0.01)
{
  alpha_str = sprintf('00%s', alpha*1000)
} else if (alpha < 0.1){
  alpha_str = sprintf('0%s', alpha*100)
}else{
  alpha_str = sprintf('%s', alpha*10)
}

if (beta < 0.1)
{
  beta_str = sprintf('0%s', beta*100)
} else{
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



folder_form_0_1 = sprintf('Obs_0_1_alpha0%s_beta0%s',alpha_str, beta_str)
folder_form_0_2 = sprintf('Obs_0_2_alpha0%s_beta0%s',alpha_str, beta_str)
folder_form_0_1_2 = sprintf('Obs_0_1_2_alpha0%s_beta0%s',alpha_str, beta_str)

folder_form_0_1_result = sprintf('Obs_0_1_alpha0%s_beta0%s_result',alpha_str, beta_str)
folder_form_0_2_result = sprintf('Obs_0_2_alpha0%s_beta0%s_result',alpha_str, beta_str)
folder_form_0_1_2_result = sprintf('Obs_0_1_2_alpha0%s_beta0%s_result',alpha_str, beta_str)

setwd('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario2/SimulateData_New_EXP10')
dir.create(folder_form_0_1)
dir.create(folder_form_0_2)
dir.create(folder_form_0_1_2)

dir.create(folder_form_0_1_result)
dir.create(folder_form_0_2_result)
dir.create(folder_form_0_1_2_result)

for (indexn in 1:1000){
  
  print(c(parameter,indexn))
  form = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/SimulateData_EXP10/RandomTree/RandomTree_%s.tre', indexn)
  
  sampletr=read.tree(form)
  
  form_0_1 = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario2/SimulateData_New_EXP10/Obs_0_1_alpha0%s_beta0%s/true_0_1_tip_alpha_0%s_beta_0%s_matrix%s.Rdata', alpha_str, beta_str, alpha_str, beta_str,indexn)
  form_0_2 = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario2/SimulateData_New_EXP10/Obs_0_2_alpha0%s_beta0%s/true_0_2_tip_alpha_0%s_beta_0%s_matrix%s.Rdata', alpha_str, beta_str, alpha_str, beta_str,indexn)
  form_0_1_2= sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario2/SimulateData_New_EXP10/Obs_0_1_2_alpha0%s_beta0%s/true_0_1_2_tip_alpha_0%s_beta_0%s_matrix%s.Rdata', alpha_str, beta_str, alpha_str, beta_str,indexn)
  
  
  obs_form_0_1 = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario2/SimulateData_New_EXP10/Obs_0_1_alpha0%s_beta0%s/obs_0_1_tip_alpha_0%s_beta_0%s_matrix%s.Rdata', alpha_str, beta_str, alpha_str, beta_str,indexn)
  obs_form_0_2 = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario2/SimulateData_New_EXP10/Obs_0_2_alpha0%s_beta0%s/obs_0_2_tip_alpha_0%s_beta_0%s_matrix%s.Rdata', alpha_str, beta_str, alpha_str, beta_str,indexn)
  obs_form_0_1_2= sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario2/SimulateData_New_EXP10/Obs_0_1_2_alpha0%s_beta0%s/obs_0_1_2_tip_alpha_0%s_beta_0%s_matrix%s.Rdata', alpha_str, beta_str, alpha_str, beta_str,indexn)
  
  sim_result = simulate_sample(sequencing_error_model,length(sampletr$edge.length),length(sampletr$tip.label),sampletr)
  
  mat_form_0_1 = sim_result$true_0_1_tip_matrix
  mat_form_0_2 = sim_result$true_0_2_tip_matrix
  mat_form_0_1_2= sim_result$true_0_1_2_tip_matrix
  
  
  mat_obs_form_0_1 = sim_result$obs_0_1_tip_matrix
  mat_obs_form_0_2 = sim_result$obs_0_2_tip_matrix
  mat_obs_form_0_1_2= sim_result$obs_0_1_2_tip_matrix
  
  save(mat_form_0_1,file = form_0_1)
  
  save(mat_form_0_2,file = form_0_2)
  
  save(mat_form_0_1_2,file = form_0_1_2)
  
  
  save(mat_obs_form_0_1,file = obs_form_0_1)
  
  save(mat_obs_form_0_2,file = obs_form_0_2)
  
  save(mat_obs_form_0_1_2,file = obs_form_0_1_2)
  
}

}
############################

#####################################################
#####################################################
#END SIMULATE
#####################################################
#####################################################

