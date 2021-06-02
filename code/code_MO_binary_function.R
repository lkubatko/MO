library(ape)
library(phangorn)

####################################################################################################
####################################################################################################
# basic functions for computing the prior and likelihood
####################################################################################################
####################################################################################################
# find parent of the node, edge is the edge ouput from tree
find_ancestor <- function(edge,node){
  parent = 0
  for (i in 1:nrow(edge)) {
    if (edge[i,2] == node) {
      parent = edge[i,1]
    }
  }
  return(parent)
}
# get all unique nodes ordered as in the tree
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
# find and return root node of the tree
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
  return (root_node)
}
# find two child branches in row 2 and two child nodes in row 1 of parent_node if they exist. Otherwise all zeros matrix output
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

# find all child branches for each branch saved in list
find_all_child_branches <- function(edge){
  
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

####################################################################################################
####################################################################################################
#basic functions end
####################################################################################################
####################################################################################################

####################################################################################################
####################################################################################################
# posterior prob function 
####################################################################################################
####################################################################################################
# alpha is fpr and beta is fnr
# unit_lambda is transition rate from normal to mutation
# initial_obsis the observation data at the tips:  0 is normal, 1 is mutation, 3 is missing value, "m" is ambiguous value
# t is the tree
MO_binary<- function(alpha,beta,unit_lambda,initial_obs,t){
  
  ts = t
  num_rows = nrow(ts$edge) #number of edges
  num_cols = length(ts$tip.label) #number of tips
  num_node = 2*length(ts$tip.label)-1

  #initialize to store the processed data. Double check the data and find out the site with missing values, ambiguous values
  obs_genotype = initial_obs[ts$tip.label]
  tip_exclude_index = which(obs_genotype == 3)
  tip_exclude = colnames(initial_obs)[tip_exclude_index]
  
  
  # find out all descendant br with above function
  descendant_branches_all = find_all_child_branches(ts$edge) 
  # find out br with missing desendant tips
  br_exclude=c()
  for (Numnode in 1:num_node){
    Numnode_des=Descendants(ts, Numnode, type = c("tips"))
    if(all(Numnode_des[[1]] %in% tip_exclude_index)){single_br_exclude=which(ts$edge[,2] == Numnode)
    single_br_exclude_all = c(single_br_exclude,descendant_branches_all[[single_br_exclude]])
    br_exclude=c(br_exclude,single_br_exclude_all)}
  }
  
  # get the subtree with some br lengths set to 0 due to missing values
  subtree = ts
  subtree$edge.length[br_exclude] = 0
  branch_time_list = subtree$edge.length
  
  
  ####################################################
  # generate the obs matrix with ambiguous sites so that there are 2^m rows and same columns as initial data
  obs_genotype_mat = matrix(,nrow = 2^length(grep("m",obs_genotype)),ncol = dim(obs_genotype)[2])
  colnames(obs_genotype_mat) = colnames(obs_genotype)
  
  # find the index of each cell status:whether it is normal,mutation,missing or ambiguous
  ambiguity_index = which(obs_genotype=="m")
  normal_index = which(obs_genotype=="0")
  mutate_index = which(obs_genotype=="1")
  missing_index = which(obs_genotype=="3")
  
  # create all possible situations for ambiguous 
  inupt_list = rep(list(0:1), length(ambiguity_index))
  input_ambiguity = expand.grid(inupt_list)
  
  # put the possible status into the matrix for gene i, each row represent one possible situation
  obs_genotype_mat[,as.numeric(ambiguity_index)] = as.matrix(input_ambiguity)
  obs_genotype_mat[,normal_index] = 0
  obs_genotype_mat[,mutate_index] = 1
  obs_genotype_mat[,missing_index] = 3
  
  # for each of the possible situation, assign weight to them, here, we use equal weights
  ambiguity_weight = matrix(rep(1/dim(obs_genotype_mat)[1],dim(obs_genotype_mat)[1], nrow = 1))


  # build the branch tree sturcture with node labels for each lineage: from each tip to the root
  branch_trees = rep( list(list()),num_cols ) 
  num_parent = 0
  for (tip_i in 1:num_cols) {
    branch_trees[[tip_i]] = tip_i
    parent = find_ancestor(subtree$edge,tip_i)
    branch_trees[[tip_i]][num_parent+2] = parent
    num_parent = num_parent+1
    
    while (parent != num_cols+1) {
      tip_node = parent
      parent = find_ancestor(subtree$edge,tip_node)
      branch_trees[[tip_i]][num_parent+2] = parent
      num_parent = num_parent+1
      
    }
    num_parent = 0
  }
  
  
  # loop over all the branches, and find the possible final stage if only one event occurs in that branch
  possible_true_genotype = matrix(rep(0,num_rows*num_cols),nrow=num_rows,ncol = num_cols)
  for (branch_i in 1:num_rows) {
    branch_edge = subtree$edge[branch_i,]
    if (branch_edge[2] <= num_cols) {
      possible_true_genotype[branch_i,branch_edge[2]] = 1
    }else {
      for (i in 1:num_cols) {
        list_branch = branch_trees[[i]]
        if (is.na(match(branch_edge[2],list_branch)) == FALSE) {
          possible_true_genotype[branch_i,i] = 1
          colnames(possible_true_genotype) = subtree$tip.label
        }
      }
    }
  }
  

  
  
  
  ###########################################################################################
  ###########################################################################################
  #Part 1: Mutation model:find the prob of mutation on each branch, but not on other branches
  ###########################################################################################
  ###########################################################################################
  #for this locus or location, find the probability of mutation on each branch by the poisson process
  
    # prob_Bj_ind is vector storig the prob that mutation occurs on each branch
    prob_Bj_ind = 1-exp(-unit_lambda*branch_time_list)
    
    # create a matrix. In each row, diagnonal is the prob on that branch, and off-branch is the mutation prob that not on that branch
    prob_Bj_mat = matrix(rep(1-prob_Bj_ind,length(branch_time_list)),nrow=length(branch_time_list), byrow = TRUE)
    
    # diagonal is replaced with prob of mutation
    diag(prob_Bj_mat) = prob_Bj_ind
    
    # descendant branches carry the mutation so the probability is 1
    for (l in 1:length(branch_time_list)){prob_Bj_mat[l,c(descendant_branches_all[[l]])] = 1}
    
    # find the marginal prob that a mutation on each br
    prob_Bj_final_all = t(apply(data.frame(prob_Bj_mat),1,cumprod))
    prob_Bj_final = prob_Bj_final_all[,ncol(prob_Bj_final_all)]
  
 
  ###################################################################################
  ###################################################################################
  #Part 2: Error model:incorporate the errors at the tips
  ###################################################################################
  ###################################################################################
  
  # dim(obs_genotype_mat)[1] is the number of possible mutation situations in the data
  # dim(obs_genotype_mat)[2] is the number of tips or samples in the subtree
  # error_result_mat is a list of matrix, length of the list equals number of branches, and each matrix is as obs_genotype_mat, which is all possible situations
  # each matrix is for one branch
  error_result_mat = replicate(dim(possible_true_genotype)[1], 
                             matrix(rep(0,dim(obs_genotype_mat)[1]*dim(obs_genotype_mat)[2]),nrow=dim(obs_genotype_mat)[1]),
                             simplify = FALSE)
  
  # error_result_mat[[t]][[k]] is the error prob if the mutation occurs on branch k
  # error_result_mat[[t]][[k]] each line corresponds to the possible observed(ambiguity) in obs_genotype_mat
  # branch k (mutation on branch k)
  for (k in 1:dim(possible_true_genotype)[1]){
    # situation j (possible situation for ambiguous status),dim(obs_genotype_mat)[1]=2^(# of ambiguity sites)
    for (j in 1:dim(obs_genotype_mat)[1]){
      # tip for sample i 
      for (i in 1:dim(possible_true_genotype)[2]){
        # if true is 1, and observed is 1, prob is 1-beta
        if (as.matrix(possible_true_genotype)[k,i] == 1 &  as.matrix(obs_genotype_mat)[j,i] == 1){error_result_mat[[k]][j,i] = 1-beta}
        
        # if true is 0 and observed is 0, prob is 1-alpha
        else if (as.matrix(possible_true_genotype)[k,i] == 0 & as.matrix(obs_genotype_mat)[j,i] == 0){error_result_mat[[k]][j,i] = 1-alpha}
        
        # if true is 1 and observed is 0, false negative, the prob is beta
        else if (as.matrix(possible_true_genotype)[k,i] == 1 & as.matrix(obs_genotype_mat)[j,i] == 0){error_result_mat[[k]][j,i] = beta}
        
        # if true is 0 and observed is 1, false positive, the prob is alpha
        else if (as.matrix(possible_true_genotype)[k,i] == 0 & as.matrix(obs_genotype_mat)[j,i] == 1){error_result_mat[[k]][j,i] = alpha}
        
        # if missing values at the tips, just return a 1 which will not have effect on the prob
        else if (as.matrix(possible_true_genotype)[k,i] == 1 & as.matrix(obs_genotype_mat)[j,i] == 3){error_result_mat[[k]][j,i] = 1}
        else if (as.matrix(possible_true_genotype)[k,i] == 0 & as.matrix(obs_genotype_mat)[j,i] == 3){error_result_mat[[k]][j,i] = 1}
      }
    }
  }
  
  
  
  #error_prob is a matrix of ncol= number of possible situation(ambiguous), nrow= number of branches
  #error_prob, each column is one possible observed mutation, each line corresponds to error prob of the mutation branch
  error_prob_prod = matrix(unlist(lapply(error_result_mat, cumprod)), 
                           nrow = length(error_result_mat), ncol = dim(error_result_mat[[1]])[2], byrow = TRUE)
  error_prob = error_prob_prod[,ncol(error_prob_prod)]
  
  
  #weight is assigned to each possible situation(ambiguous), and the total weighted prob is calculated
  #weighted_error is the probability of observed conditioning on a mutation branch, each row represents one branch
  weighted_error = error_prob%*%ambiguity_weight
  
  #############################################################################################################
  #############################################################################################################
  #Part 3: Posterior probability that a mutation on each of the branch based on mutation model and error model 
  #############################################################################################################
  #############################################################################################################
  
  prob_Bj_S = prob_Bj_final * as.vector(weighted_error)/sum(prob_Bj_final * as.vector(weighted_error))
  
  return(prob_Bj_S)
}
####################################################################################################
####################################################################################################
# posterior prob function end
####################################################################################################
####################################################################################################
