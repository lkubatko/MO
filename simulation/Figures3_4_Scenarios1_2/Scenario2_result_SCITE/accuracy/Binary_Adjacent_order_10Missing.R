library(plyr)
library(ape)
compare_near_order <- function(number_br,number_cell,ts){
  
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
  
  
  # find neighboring child branch for current branch
  find_neighboring_child_branches <- function(edge, current_edge, child_branches)
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
      
      return(child_branches)
    }
    else
    {
      return(child_branches)
    }
    
  }
  
  # find following child branch for all branches
  find_all_neighboring_child_branches <- function(edge){
    neighbor_child_branches = rep(list(list()), nrow(edge)) 
    
    for (i in 1:nrow(edge))
    {
      current_edge = i
      # iterative find all its child branches
      child_branches = integer(0)
      neighbor_child_branches[[i]] = find_neighboring_child_branches(edge, current_edge, child_branches)
      
    }
    
    return(neighbor_child_branches)
  }
  
  neighboring_child_branches = find_all_neighboring_child_branches(ts$edge)
  
  
  return(neighboring_child_branches=neighboring_child_branches)
  
  
  
  
}

parameter_setting = expand.grid(alpha=c(0.05,0.1,0.2,0.4),
                                beta =c(0.05,0.1,0.2,0.4))


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
  
  
scite_near_pair_result= data.frame(matrix(, nrow=0, ncol=3))
scite_near_pair_result_form_0_1 = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario2_result_SCITE/SCITE_Binary_10missing/summary/binary_adjacent_order_alpha_0%s_beta_0%s.csv', alpha_str, beta_str)


for (indexn in 1:1000){
  
  
  #read true tree
  form = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting/SimulateData_EXP10/RandomTree/RandomTree_%s.tre', indexn)
  sampletr=read.tree(form)
  
  #read inferred results
  scite_order_form_0_1 = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario2_result_SCITE/SCITE_Binary_10missing/SCITE_binary_alpha_0%s_beta_0%s/sub_binary_obs_0_1_tip_alpha_0%s_beta_0%s_node_order_%s.txt',alpha_str, beta_str, alpha_str, beta_str,indexn)
  mat_scite_order_form_0_1 = scan(file=scite_order_form_0_1) 
  root_tree=mat_scite_order_form_0_1[1]
  root_index=which(mat_scite_order_form_0_1%in% root_tree )
  
  mat_scite_order_form_0_1_trans=list()
  scite_near_match_matchresult= data.frame(matrix(, nrow=0, ncol=2))
  for (i in 1:length(root_index)){
    
    if(i<length(root_index)){mat_scite_order_form_0_1_trans[[i]]=mat_scite_order_form_0_1[c(root_index[i]:(root_index[i+1]-1))]}
    else{mat_scite_order_form_0_1_trans[[i]]=mat_scite_order_form_0_1[c(root_index[i]:length(mat_scite_order_form_0_1))]}
    
  }
  
  mat_scite_order_form_0_1_near_order_pair=c()
  
  for (i in 1:length(mat_scite_order_form_0_1_trans)){
    
    near_order_pair=c()
    for (j in 1:(length(mat_scite_order_form_0_1_trans[[i]])-1)){
      near_order_pair=rbind(near_order_pair,c((mat_scite_order_form_0_1_trans[[i]][j]),mat_scite_order_form_0_1_trans[[i]][j+1]))
      
    }
    mat_scite_order_form_0_1_near_order_pair=rbind(mat_scite_order_form_0_1_near_order_pair,near_order_pair)
  }
  
  exclude_node = mat_scite_order_form_0_1_trans[[i]][1]
  mat_scite_order_form_0_1_near_order_pair = as.data.frame(mat_scite_order_form_0_1_near_order_pair)
  colnames(mat_scite_order_form_0_1_near_order_pair)=c("first_br_indexn","second_br_indexn")
  mat_scite_order_form_0_1_near_order_pair = na.omit(mat_scite_order_form_0_1_near_order_pair)
  mat_scite_order_form_0_1_near_order_pair<-mat_scite_order_form_0_1_near_order_pair[!(mat_scite_order_form_0_1_near_order_pair$first_br_indexn == (exclude_node) ),]
  mat_scite_order_form_0_1_near_order_pair<-mat_scite_order_form_0_1_near_order_pair[!(mat_scite_order_form_0_1_near_order_pair$second_br_indexn == (exclude_node) ),]
  
  
  true_near_branch = compare_near_order(length(sampletr$edge.length),length(sampletr$tip.label),sampletr)
  
  all_true_near_branch = data.frame(matrix(, nrow=0, ncol=2))
  all_match_near_branch_sub = data.frame(matrix(, nrow=0, ncol=2))
  
  for (br in 1:length(true_near_branch)){
    
    if(length(true_near_branch[[br]])>0){
      for(br_sub in 1:length(true_near_branch[[br]])){
        all_true_near_branch = rbind(all_true_near_branch,c(br,(true_near_branch[[br]])[br_sub]))
      }
    }
  }
  
  all_true_near_branch_sub=na.omit(all_true_near_branch)
  colnames(all_true_near_branch_sub) = c("first_br_indexn","second_br_indexn")
  
  unique_mat_scite_order_form_0_1_near_order_pair=unique(mat_scite_order_form_0_1_near_order_pair)
  
  for(orderpair in 1:dim(unique_mat_scite_order_form_0_1_near_order_pair)[1]){
    print(c(indexn,orderpair))
    matchresult = match_df(all_true_near_branch_sub, unique_mat_scite_order_form_0_1_near_order_pair[orderpair,])
    scite_near_match_matchresult=rbind(scite_near_match_matchresult,matchresult)
  }
  
  scite_near_pair_result[indexn,]=c(dim(scite_near_match_matchresult)[1],dim(unique_mat_scite_order_form_0_1_near_order_pair)[1],dim(all_true_near_branch_sub)[1])
}
  
  colnames(scite_near_pair_result)=c("matched_paired","inferred_paired","true_pair")
  scite_near_pair_result_form_0_1 = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario2_result_SCITE/SCITE_Binary_10missing/summary/binary_adjacent_order_alpha_0%s_beta_0%s.csv', alpha_str, beta_str)
  
  write.csv(scite_near_pair_result,file=scite_near_pair_result_form_0_1)
  
}