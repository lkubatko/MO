
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

find_mutation_prob_Bj <- function(branch_time_list, unit_theta,unit_gamma,unit_mu,branch_i){
  ###############################################
  #find the prob of only 1 mutation on branch_i
  ###############################################
  # find all branches numbers in the subtree
  all_branches=c(1:dim(subtree$edge)[1])
  descendant_branches=find_all_child_branches(subtree$edge)
  # 1-1:find all descendant branches of the input branch, all these branches carry the mutation and no further mutation occurs
  Bj_with_only_1_descendant_branches=descendant_branches[[branch_i]]
  # 1-1:find the prob that all descendant branches of the input branch carry the mutation and no further mutation occurs
  prob_Bj_with_only_1_descendant_branches=exp(-unit_mu*branch_time_list[Bj_with_only_1_descendant_branches])
  
  
  # 0-0:find the complement branches of input branch and its descendant, branches in this set has no mutation
  Bj_with_only_1_no_mutation_branches=setdiff(all_branches,union(branch_i, Bj_with_only_1_descendant_branches))
  # 0-0:find the prob the complement branches of input branch and its descendant in this set has no mutation
  prob_Bj_with_only_1_no_mutation_branches=exp(-(unit_theta+unit_gamma)*(branch_time_list[Bj_with_only_1_no_mutation_branches]))
  
  
  # get the prob matrix of the only 1 mutation
  prob_Bj_with_only_1_mat= matrix(, nrow = 1, ncol = dim(subtree$edge)[1])
  prob_Bj_with_only_1_mat[,Bj_with_only_1_descendant_branches]=prob_Bj_with_only_1_descendant_branches
  prob_Bj_with_only_1_mat[,Bj_with_only_1_no_mutation_branches]=prob_Bj_with_only_1_no_mutation_branches
  prob_Bj_with_only_1_mat[,branch_i]= unit_theta*(exp(-(unit_theta+unit_gamma)*branch_time_list[branch_i])-exp(-unit_mu*branch_time_list[branch_i]))/(unit_mu-unit_theta-unit_gamma)
  prob_Bj_with_only_1=cbind(First_branch = branch_i,
                            Second_branch = 0,
                            Mutation_Prob=cumprod(prob_Bj_with_only_1_mat)[dim(subtree$edge)[1]])
  
  
  
  ###############################################
  #find the prob of only 2 mutation on branch_i
  ############################################### 
  
  #2-2:find all descendant branches of the input branch, all these branches carry the mutation and no further mutation occurs
  Bj_with_only_2_descendant_branches=descendant_branches[[branch_i]]
  #2-2:find the prob that all descendant branches of the input branch carry the 2 mutations and no further mutation occurs
  prob_Bj_with_only_2_descendant_branches=exp(-0*branch_time_list[Bj_with_only_2_descendant_branches])
  
  
  #0-0:find the complement branches of input branch and its descendant, branches in this set has no mutation
  Bj_with_only_2_no_mutation_branches=setdiff(all_branches,union(branch_i, Bj_with_only_2_descendant_branches))
  #0-0:find the prob the complement branches of input branch and its descendant in this set has no mutation
  prob_Bj_with_only_2_no_mutation_branches=exp(-(unit_theta+unit_gamma)*(branch_time_list[Bj_with_only_2_no_mutation_branches]))
  
  # get the prob matrix of the only 2 mutation
  prob_Bj_with_only_2_mat= matrix(, nrow = 1, ncol = dim(subtree$edge)[1])
  prob_Bj_with_only_2_mat[,Bj_with_only_2_descendant_branches]=prob_Bj_with_only_2_descendant_branches
  prob_Bj_with_only_2_mat[,Bj_with_only_2_no_mutation_branches]=prob_Bj_with_only_2_no_mutation_branches
  prob_Bj_with_only_2_mat[,branch_i]= ((unit_gamma-unit_mu)*(exp(-(unit_theta+unit_gamma)*branch_time_list[branch_i]))+(unit_theta)*(exp(-unit_mu*branch_time_list[branch_i])))/(unit_mu-unit_theta-unit_gamma)+1
  prob_Bj_with_only_2=cbind(First_branch = branch_i,
                            Second_branch = branch_i,
                            Mutation_Prob=cumprod(prob_Bj_with_only_2_mat)[dim(subtree$edge)[1]])
  
  
  ########################################################################
  #find the prob of 0-1 mutation on branch_i and 1-2 on a descendant branch
  ########################################################################
  Bj_with_2_descendant_branches=descendant_branches[[branch_i]]
  if(length(Bj_with_2_descendant_branches) > 0){
    # create the matrix that store the prob if there are separate mutations on two branches
    prob_Bj_with_2_mat=matrix(, nrow = length(Bj_with_2_descendant_branches), ncol = 2+dim(subtree$edge)[1])
    
    # iterate on each of the descendant branch
    for (descendant_i in 1:length(Bj_with_2_descendant_branches)){
      
      prob_Bj_with_2_mat[descendant_i,1]= branch_i
      prob_Bj_with_2_mat[descendant_i,2]= Bj_with_2_descendant_branches[descendant_i]
      
      DB_descendant_i=descendant_branches[[Bj_with_2_descendant_branches[descendant_i]]]
      
      # fin the branches without mutation
      Bj_with_2_no_mutation_branches = setdiff(all_branches,union(branch_i, Bj_with_2_descendant_branches))
      prob_Bj_with_2_no_mutation_branches = exp(-(unit_theta+unit_gamma)*(branch_time_list[Bj_with_2_no_mutation_branches]))
      
      # the branch with the first mutation occurring
      Bj_with_2_1st_mutation_branches = branch_i
      prob_Bj_with_2_1st_mutation_branches = unit_theta*(exp(-(unit_theta+unit_gamma)*branch_time_list[branch_i])-exp(-unit_mu*branch_time_list[branch_i]))/(unit_mu-unit_theta-unit_gamma)
      
      # the branch with the second mutation occurring
      Bj_with_2_2nd_mutation_branches = Bj_with_2_descendant_branches[descendant_i]
      prob_Bj_with_2_2nd_mutation_branches = 1-exp(-unit_mu*branch_time_list[Bj_with_2_2nd_mutation_branches])
      
      
      #the branches carry the first and second mutation
      Bj_with_2_carry_1st_2nd_mutation_branches = descendant_branches[[Bj_with_2_2nd_mutation_branches]]
      prob_Bj_with_2_carry_1st_2nd_mutation_branches = exp(0*branch_time_list[Bj_with_2_carry_1st_2nd_mutation_branches])
      #the branches carry the first mutation
      Bj_with_2_carry_1st_mutation_branches = setdiff(Bj_with_2_descendant_branches,union(Bj_with_2_2nd_mutation_branches,Bj_with_2_carry_1st_2nd_mutation_branches))
      prob_Bj_with_2_carry_1st_mutation_branches = exp(-unit_mu*branch_time_list[Bj_with_2_carry_1st_mutation_branches])
      
      #put these probability into the matrix
      prob_Bj_with_2_mat[descendant_i,Bj_with_2_no_mutation_branches+2] = prob_Bj_with_2_no_mutation_branches
      prob_Bj_with_2_mat[descendant_i,Bj_with_2_1st_mutation_branches+2] = prob_Bj_with_2_1st_mutation_branches
      prob_Bj_with_2_mat[descendant_i,Bj_with_2_carry_1st_mutation_branches+2] = prob_Bj_with_2_carry_1st_mutation_branches
      prob_Bj_with_2_mat[descendant_i,Bj_with_2_2nd_mutation_branches+2] = prob_Bj_with_2_2nd_mutation_branches
      prob_Bj_with_2_mat[descendant_i,Bj_with_2_carry_1st_2nd_mutation_branches+2] = prob_Bj_with_2_carry_1st_2nd_mutation_branches
      
      #find the probability of 
      prob_Bj_with_2=cbind(prob_Bj_with_2_mat[,c(1,2)],
                           Mutation_Prob=t(apply(prob_Bj_with_2_mat[,-c(1,2)],1,cumprod))[,dim(subtree$edge)[1]])
      
    }
    
  }else{prob_Bj_with_2=cbind(First_branch = branch_i,
                             Second_branch = 0,
                             Mutation_Prob = 0)}
  
  return(list(prob_0_1= prob_Bj_with_only_1, 
              prob_0_2 = prob_Bj_with_only_2, 
              prob_0_1_2 = prob_Bj_with_2))
}

parameter_setting = expand.grid(alpha=c(0.05,0.1,0.2,0.4),
                                beta =c(0.05,0.1,0.2,0.4))


for (paraind in 1:dim(parameter_setting)[1]){
  
  alpha = parameter_setting[paraind,1]
  beta = parameter_setting[paraind,2]
  
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
  unit_theta = 10^(-7)
  unit_gamma = 10^(-14)
  unit_mu = 10 ^(-2)
  number_br = 98
  number_cell = 50
  
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
  
  
  binary_folder_form_10 = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_10Missing/Binary_alpha0%s_beta0%s',alpha_str, beta_str)
  ternary_folder_form_10 = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_10Missing/Ternary_alpha0%s_beta0%s',alpha_str, beta_str)
  
  
  dir.create(binary_folder_form_10)
  dir.create(ternary_folder_form_10)
  
  binary_folder_form_20 = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_20Missing/Binary_alpha0%s_beta0%s',alpha_str, beta_str)
  ternary_folder_form_20 = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_20Missing/Ternary_alpha0%s_beta0%s',alpha_str, beta_str)
  
  
  dir.create(binary_folder_form_20)
  dir.create(ternary_folder_form_20)
  
  set.seed(100)
  for (indexn in 1:100){
    
    binary_transformed_obs_dat_out = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData/Binary_alpha0%s_beta0%s/binary_obs_0_1_tip_alpha_0%s_beta_0%s_matrix%s.csv', alpha_str, beta_str, alpha_str, beta_str,indexn)
    
    transformed_obs_dat_Binary = read.csv(file = binary_transformed_obs_dat_out)
    
    ternary_transformed_obs_dat_out = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData/Ternary_alpha0%s_beta0%s/ternary_obs_0_1_tip_alpha_0%s_beta_0%s_matrix%s.csv', alpha_str, beta_str, alpha_str, beta_str,indexn)
    
    transformed_obs_dat_Ternary = read.csv(file = ternary_transformed_obs_dat_out)
    
    
    #10 percent missing
    transformed_obs_dat_Ternary_10Missing=data.frame(matrix(, nrow=0, ncol=0))
    
    transformed_obs_dat_Binary_10Missing=data.frame(matrix(, nrow=0, ncol=0))
    
    #20 percent missing
    transformed_obs_dat_Ternary_20Missing=data.frame(matrix(, nrow=0, ncol=0))
    
    transformed_obs_dat_Binary_20Missing=data.frame(matrix(, nrow=0, ncol=0))
    
    
    for (i in 1:98){
      
      print(c(indexn,i))
      
      missing_index_20Missing = sample(1:50, 10, replace=FALSE)
      missing_index_10Missing = sample(missing_index_20Missing , 5, replace=FALSE)
      
      
      transformed_obs_Ternary = transformed_obs_dat_Ternary[i,2:53]
      transformed_obs_Binary = transformed_obs_dat_Binary[i,2:53]
      
      transformed_obs_Ternary_20Missing = transformed_obs_Ternary
      transformed_obs_Binary_20Missing = transformed_obs_Binary
      
      transformed_obs_Ternary_10Missing = transformed_obs_Ternary
      transformed_obs_Binary_10Missing = transformed_obs_Binary
      
      transformed_obs_Ternary_20Missing[1,(missing_index_20Missing+2)] = c("-")
      transformed_obs_Binary_20Missing[1,(missing_index_20Missing+2)] = c("-")
      
      transformed_obs_Ternary_10Missing[1,(missing_index_10Missing+2)] = c("-")
      transformed_obs_Binary_10Missing[1,(missing_index_10Missing+2)] = c("-")
      
      transformed_obs_dat_Ternary_20Missing=rbind(transformed_obs_dat_Ternary_20Missing,transformed_obs_Ternary_20Missing)
      transformed_obs_dat_Binary_20Missing=rbind(transformed_obs_dat_Binary_20Missing,transformed_obs_Binary_20Missing)
      
      transformed_obs_dat_Ternary_10Missing=rbind(transformed_obs_dat_Ternary_10Missing,transformed_obs_Ternary_10Missing)
      transformed_obs_dat_Binary_10Missing=rbind(transformed_obs_dat_Binary_10Missing,transformed_obs_Binary_10Missing)
      

    
      
    }
    
    #10 percent missing
    binary_transformed_obs_dat_out_10Missing = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_10Missing/Binary_alpha0%s_beta0%s/binary_obs_0_1_tip_alpha_0%s_beta_0%s_matrix%s.csv', alpha_str, beta_str, alpha_str, beta_str,indexn)
    
    write.csv(transformed_obs_dat_Binary_10Missing,file = binary_transformed_obs_dat_out_10Missing)
    
    ternary_transformed_obs_dat_out_10Missing = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_10Missing/Ternary_alpha0%s_beta0%s/ternary_obs_0_1_tip_alpha_0%s_beta_0%s_matrix%s.csv', alpha_str, beta_str, alpha_str, beta_str,indexn)
    
    write.csv(transformed_obs_dat_Ternary_10Missing,file = ternary_transformed_obs_dat_out_10Missing)
    
    #20 percent missing
    binary_transformed_obs_dat_out_20Missing = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_20Missing/Binary_alpha0%s_beta0%s/binary_obs_0_1_tip_alpha_0%s_beta_0%s_matrix%s.csv', alpha_str, beta_str, alpha_str, beta_str,indexn)
    
    write.csv(transformed_obs_dat_Binary_20Missing,file = binary_transformed_obs_dat_out_20Missing)
    
    ternary_transformed_obs_dat_out_20Missing = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_20Missing/Ternary_alpha0%s_beta0%s/ternary_obs_0_1_tip_alpha_0%s_beta_0%s_matrix%s.csv', alpha_str, beta_str, alpha_str, beta_str,indexn)
    
    write.csv(transformed_obs_dat_Ternary_20Missing,file = ternary_transformed_obs_dat_out_20Missing)
    
  }
  
  
}

