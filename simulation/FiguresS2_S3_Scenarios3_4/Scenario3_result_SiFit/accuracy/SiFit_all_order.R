library(plyr)
library(ape)
compare_all_order <- function(number_br,number_cell,ts){
  
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
  
  
  
  all_child_branches =  find_all_child_branches(ts$edge)
  
  
  return(all_child_branches=all_child_branches)
  
  
  
  
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  stop("Error: No alpha and beta values are given!\n", call.=FALSE)
} else {
  alpha = as.double(args[1])
  beta = as.double(args[2])
  print(paste("alpha is ", alpha))
  print(paste("beta is ", beta))
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








sifit_all_pair_result= data.frame(matrix(, nrow=0, ncol=4))

for (indexn in 1:100){
  
  #form = sprintf('G:/SimulateData/RandomTree/RandomTree_%s.tre', indexn)
  form = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData/RandomTree/RandomTree_%s.tre', indexn)
  sampletr=read.tree(form)
  
  sifit_tree_form_0_1 = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SiFit_SimulateData/SiFit_binary_alpha_0%s_beta_0%s_TrueTree/True_sub_ternary_obs_0_1_tip_alpha_0%s_beta_0%s_matrix%s_mlTree.newick',alpha_str, beta_str, alpha_str, beta_str,indexn)
  mat_sifit_tree_form_0_1 = read.tree(sifit_tree_form_0_1)
  
  sifit_br_form_0_1 = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SiFit_SimulateData/SiFit_binary_alpha_0%s_beta_0%s_TrueTree/TrueTree_outcome_order_alpha_0%s_beta_0%s_matrix%s.txt.in.txt',alpha_str, beta_str, alpha_str, beta_str,indexn)
  mat_sifit_br_form_0_1 = read.csv(sifit_br_form_0_1,sep="\t",header=FALSE)
  
  sifit_order_form_0_1 = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SiFit_SimulateData/SiFit_binary_alpha_0%s_beta_0%s_TrueTree/TrueTree_outcome_order_alpha_0%s_beta_0%s_matrix%s.txt.par_child.txt',alpha_str, beta_str, alpha_str, beta_str,indexn)
  
  if(file.size(sifit_order_form_0_1)>0){
  mat_sifit_order_form_0_1 = read.csv(sifit_order_form_0_1,sep="\t",header=FALSE)
  colnames(mat_sifit_order_form_0_1)=c("parent","child","mutation")
  
  mat_sifit_br_form_0_1[,3]=match(mat_sifit_br_form_0_1[,2],mat_sifit_tree_form_0_1$edge.length)
  colnames(mat_sifit_br_form_0_1)=c("child","edgelength","truebr")
  
  mat_joint_all=merge(x=mat_sifit_order_form_0_1,y=mat_sifit_br_form_0_1,by="child",all=TRUE)
  mat_joint_all=merge(x=mat_sifit_br_form_0_1,y=mat_sifit_order_form_0_1,by="child",all=TRUE)
  
  mat_joint=na.omit(mat_joint_all)
  
  true_all_branch = compare_all_order(length(sampletr$edge.length),length(sampletr$tip.label),sampletr)
  sifit_all_branch = compare_all_order(length(mat_sifit_tree_form_0_1$edge.length),length(mat_sifit_tree_form_0_1$tip.label),mat_sifit_tree_form_0_1)
  
  
  sifit_all_match_matchresult= data.frame(matrix(, nrow=0, ncol=2))
  all_true_all_branch = data.frame(matrix(, nrow=0, ncol=2))
  all_sifit_all_branch = data.frame(matrix(, nrow=0, ncol=2))
  
  for (br in 1:length(sifit_all_branch)){
    
    if(length(sifit_all_branch[[br]])>0){
      for(br_sub in 1:length(sifit_all_branch[[br]])){
        all_sifit_all_branch = rbind(all_sifit_all_branch,c(br,(sifit_all_branch[[br]])[br_sub]))
      }
    }
  }
  
  for (br in 1:length(true_all_branch)){
    
    if(length(true_all_branch[[br]])>0){
      for(br_sub in 1:length(true_all_branch[[br]])){
        all_true_all_branch = rbind(all_true_all_branch,c(br,(true_all_branch[[br]])[br_sub]))
      }
    }
  }
  
  all_true_all_branch_sub = na.omit(all_true_all_branch)
  colnames(all_true_all_branch_sub)= c("first_br_indexn","second_br_indexn")
  
  all_sifit_all_branch_sub=data.frame(matrix(, nrow=0, ncol=2))
  for (i in 1:dim(all_sifit_all_branch)[1]){
    if((all_sifit_all_branch[i,1] %in% mat_joint$truebr)& (all_sifit_all_branch[i,2] %in% mat_joint$truebr)){
      all_sifit_all_branch_sub=rbind(all_sifit_all_branch_sub,all_sifit_all_branch[i,])
    }
  }
  colnames(all_sifit_all_branch_sub)= c("first_br_indexn","second_br_indexn")
  
  
  matched_inferred_sifit_all_branch_sub=data.frame(matrix(, nrow=0, ncol=2))
  for( i in 1:dim(all_sifit_all_branch_sub)[1]){
    
      
matched_inferred_sifit_all_branch_sub=rbind(matched_inferred_sifit_all_branch_sub,c(mat_joint$mutation[which(mat_joint$truebr %in% all_sifit_all_branch_sub[i,1])],mat_joint$mutation[which(mat_joint$truebr %in% all_sifit_all_branch_sub[i,2])]))
    
  }
  colnames(matched_inferred_sifit_all_branch_sub)= c("first_br_indexn","second_br_indexn")
  
  
  for(orderpair in 1:dim(matched_inferred_sifit_all_branch_sub)[1]){
    print(c(indexn,orderpair))
    matchresult = match_df(all_true_all_branch_sub, matched_inferred_sifit_all_branch_sub[orderpair,])
    sifit_all_match_matchresult=rbind(sifit_all_match_matchresult,matchresult)
  }
  
  sifit_all_pair_result[indexn,]=c(dim(sifit_all_match_matchresult)[1],dim(all_sifit_all_branch_sub)[1],dim(all_true_all_branch_sub)[1],dim(all_sifit_all_branch)[1])}
  if (file.size(sifit_order_form_0_1)==0){
    
    
    true_all_branch = compare_all_order(length(sampletr$edge.length),length(sampletr$tip.label),sampletr)
    sifit_all_branch = compare_all_order(length(mat_sifit_tree_form_0_1$edge.length),length(mat_sifit_tree_form_0_1$tip.label),mat_sifit_tree_form_0_1)
    
    
    all_true_all_branch = data.frame(matrix(, nrow=0, ncol=2))
    all_sifit_all_branch = data.frame(matrix(, nrow=0, ncol=2))
    
    for (br in 1:length(sifit_all_branch)){
      
      if(length(sifit_all_branch[[br]])>0){
        for(br_sub in 1:length(sifit_all_branch[[br]])){
          all_sifit_all_branch = rbind(all_sifit_all_branch,c(br,(sifit_all_branch[[br]])[br_sub]))
        }
      }
    }
    
    for (br in 1:length(true_all_branch)){
      
      if(length(true_all_branch[[br]])>0){
        for(br_sub in 1:length(true_all_branch[[br]])){
          all_true_all_branch = rbind(all_true_all_branch,c(br,(true_all_branch[[br]])[br_sub]))
        }
      }
    }
    
    all_true_all_branch_sub = na.omit(all_true_all_branch)
    colnames(all_true_all_branch_sub)= c("first_br_indexn","second_br_indexn")
    
    
    
    sifit_all_pair_result[indexn,]=c(0,0,dim(all_true_all_branch_sub)[1],dim(all_sifit_all_branch)[1])}
}
colnames(sifit_all_pair_result)=c("matched_paired","inferred_paired","true_pair","sifit_pair")
sifit_all_pair_result_form_0_1 = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SiFit_SimulateData/Binary_Summary/binary_ALL_order_matrix_01_0_1_out_alpha_0%s_beta_0%s.csv', alpha_str, beta_str)
write.csv(sifit_all_pair_result,file=sifit_all_pair_result_form_0_1)
