library(ape)
generate_prob_binary_Missing<- function(alpha,beta,unit_lambda,number_br,number_cell,normal_genotype,mutation_genotype,initial_obs,ts){
  ############################################
  
  # generating sub tree
  
  # for this gene or location, variable singletip_exclude are those sample equal '-'
  single_tip_exclude=c()
  # for this gene or location, variable tip_exclude are all those equal '-'
  tip_exclude=c()
  tip_exclude_index=c()
  # for this gene or location, obs_colnam are those sample names, excluding the '-'
  obs_colnam=c()
  #
  obs_genotype=c()
  
  # assign the mutation status for sample j of gene i
  for (j in c(1:(number_cell))){
    # exclude those tips with missing gene status
    if (initial_obs[j]==c("-")) { single_tip_exclude=colnames(initial_obs)[j] # single_tip_exclude=colnames(initial_obs)[1,j]
    tip_exclude=c(tip_exclude,single_tip_exclude)
    tip_exclude_index=c(tip_exclude_index,j)
    obs_genotype=c(obs_genotype,3)
    }
    # value is 0 if gene status is same as normal
    else if (as.character(initial_obs[1,j])==as.character(normal_genotype)) {
      obs_genotype=c(obs_genotype,0)
    }
    # value is 1 if gene status is same as mutant
    else if (as.character(initial_obs[1,j])==as.character(mutation_genotype)) {
      obs_genotype=c(obs_genotype,1)
    }
    # value is m if gene status is ambguity
    else {obs_genotype=c(obs_genotype,"m")}
  }
  
  # for this gene or location, exclude the sample with missing gene status
  #subtree=drop.tip(ts, tip=tip_exclude)
  
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
  
  
  #select those br with missing tips
  descendant_branches_all=find_all_child_branches(ts$edge)
  
  br_exclude=c()
  for (Numnode in (1:(number_cell+ts$Nnode))){
    Numnode_des=Descendants(ts, Numnode, type = c("tips"))
    if(all(Numnode_des[[1]] %in% tip_exclude_index)){single_br_exclude=which(ts$edge[,2] == Numnode)
    single_br_exclude_all = c(single_br_exclude,descendant_branches_all[[single_br_exclude]])
    br_exclude=c(br_exclude,single_br_exclude_all)}
  }
  
  # branch_time_list is the branch length of sub tree
  subtree=ts
  subtree$edge.length[br_exclude]=0
  branch_time_list= subtree$edge.length

  ####################################################
  #generate the obs matrix
  # obs_colnam are those with observations
  obs_colnam=(colnames(initial_obs))
  # consider the ambguity status as missing
  
  obs_genotype_mat=matrix(,nrow=2^length(grep("m",obs_genotype)),ncol=length(obs_colnam))
  colnames(obs_genotype_mat)=obs_colnam
  # find the index of each gene status
  ambiguity_index=which(obs_genotype=="m")
  normal_index=which(obs_genotype=="0")
  allele_index=which(obs_genotype=="1")
  missing_index=which(obs_genotype=="3")
  #create all possible situations for ambguity 
  inupt_list <- rep(list(0:1), length(ambiguity_index))
  input_ambiguity=expand.grid(inupt_list)
  # put the possible status into the matrix for gene i, each row represent one possible situation
  obs_genotype_mat[,as.numeric(ambiguity_index)]=as.matrix(input_ambiguity)
  obs_genotype_mat[,normal_index]=rep(0,dim(obs_genotype_mat)[1])
  obs_genotype_mat[,allele_index]=rep(1,dim(obs_genotype_mat)[1])
  obs_genotype_mat[,missing_index]=rep(3,dim(obs_genotype_mat)[1])
  # for each of the possible situation, assign weight to them, here, I use equal weights
  ambiguity_weight=matrix(rep(1/dim(obs_genotype_mat)[1],dim(obs_genotype_mat)[1],nrow=1))

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
  
  #child_branches=find_all_child_branches(ts$edge)
  #####################################modify end################################
  
  # save final stage of each tree on condition that the event occurs in one branch
  
  
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
          colnames(possible_true_genotype)=subtree$tip.label
        }
      }
    }
  }
  
  possible_true_genotype
  
  descendant_branches=find_all_child_branches(subtree$edge)
  
  num_of_cases = (num_rows+1+2)*num_rows/2
  possible_true_genotype_with_2 = matrix(rep(0,num_of_cases*(num_cols+2)),nrow=num_of_cases,ncol = num_cols+2)
  colnames(possible_true_genotype_with_2)= c("First_branch", "Second_branch", subtree$tip.label)
  # a matrix is created, where first two columns are the places for the mutation occuring.
  # if it's NA, it means that no mutation occurs to stand for situation like 0, 1, 2
  id_row = 1
  for (branch_i in 1:num_rows) {
    branch_edge = subtree$edge[branch_i,]
    # if only this branch has one mutation
    possible_true_genotype_with_2[id_row,1] = branch_i
    possible_true_genotype_with_2[id_row,1:num_cols+2] = possible_true_genotype[branch_i, ]
    id_row = id_row+1
    # if this branch has one mutation, and other branch has another
    for (branch_j in branch_i:num_rows) {
      possible_true_genotype_with_2[id_row,1] = branch_i
      possible_true_genotype_with_2[id_row,2] = branch_j
      
      possible_true_genotype_with_2[id_row,1:num_cols+2] = possible_true_genotype[branch_i, ] + possible_true_genotype[branch_j, ]
      id_row = id_row+1
    }
  }
  
  
  
  ###################################################################################
  ###################################################################################
  #Mutation model:find the prob of mutation on each branch, but not on other branches
  ###################################################################################
  ###################################################################################
  
  # branch_time_list is the branch length of sub tree
  branch_time_list= subtree$edge.length
  ##################################################
  #for gene i, find the prob on each branch
  #for this gene or location, find the probability of mutation on each branch by function
  find_mutation_prob_Bj <- function(branch_time_list, unit_lambda){
    
    prob_Bj_ind=c()
    for (l in 1:length(branch_time_list)){
      #assuming the exponential distribution on each branch,prob_Bj_ind[l] is the prob of mutation on branch l 
      prob_Bj_ind[l]=1-(exp(-unit_lambda*branch_time_list[l]))
    }
    # create a matrix. In each row, diagnonal is the prob on that branch, and off-branch is the mutation prob that not on that branch
    # probability matrix when mutation exist on only one branch, the matrix is first created with same rows, and each row is the prob of no mutation
    prob_Bj_mat=matrix(rep(1-prob_Bj_ind,length(branch_time_list)),nrow=length(branch_time_list),byrow=TRUE)
    # diagonal is replaced with prob of mutation
    diag(prob_Bj_mat) <- prob_Bj_ind
    # descendant branches carry the mutation so the probability is 1
    for (l in 1:length(branch_time_list)){ prob_Bj_mat[l,c(descendant_branches[[l]])]=1}
    #take the product of each row to find the prob of mutation on branch but not on other branches
    #first create a matrix of one column and each row represent the final prob
    prob_Bj_final=matrix(,nrow=length(branch_time_list),ncol=1)# probability of mutation on branch and other branch has no mutation
    
    for (l in 1:length(branch_time_list)){
      prob_Bj_final[l,]=tail(cumprod(prob_Bj_mat[l,]),1)
    }
    
    return(prob_Bj_final)
  }
  #use the above function to find the probability
  prob_Bj_final = find_mutation_prob_Bj(branch_time_list, unit_lambda)
  ####################################################
  
  
  # all_possible_true_genotype store the information of mutation of tips if mutation on a specific branch
  # each row number represent the branch number, and if mutation on that branch, the mutation status on the tips
  #################################################################
  #################################################################
  
  #create the error matrix
  sequencing_error_model=matrix(c(1-alpha,alpha,beta,1-beta),nrow=2,byrow = TRUE)
  
  
  # dim(all_possible_true_genotype[[t]])[1] is the number of branches in the subtree
  # dim(obs_genotype_mat)[1] is the number of possible mutation situations in the data
  # dim(obs_genotype_mat)[2] is the number of tips or samples in the subtree
  # error_result_mat is a list of matrix, in the list, the number of matrix equals number of branches, and each matrix is as obs_genotype_mat, which is all possible situations
  # each matrix is for one branch
  error_result_mat=replicate(dim(possible_true_genotype)[1], 
                             matrix(rep(0,dim(obs_genotype_mat)[1]*dim(obs_genotype_mat)[2]),nrow=dim(obs_genotype_mat)[1]),
                             simplify=FALSE)
  
  
  
  ####################################################################
  ####################################################################
  
  # error_result_mat[[t]][[k]] is the error prob if the gene mutation occurs on branch k
  # error_result_mat[[t]][[k]] each line correspond to the possible observed(ambiguity) in obs_genotype_mat
  
  # branch k (mutation on branch k)
  for (k in 1:dim(possible_true_genotype)[1]){
    # situation j (possible situation for ambguity status),dim(obs_genotype_mat)[1]=2^(# of ambiguity sites)
    for (j in 1:dim(obs_genotype_mat)[1]){
      # tip or sample i 
      for (i in 1:dim(possible_true_genotype)[2]){
        # if true is 1, and observed is 1, prob is 1-beta
        if (as.matrix(possible_true_genotype)[k,i]==1 &  as.matrix(obs_genotype_mat)[j,i]==1){error_result_mat[[k]][j,i]=1-beta}
        
        # if true is 0 and observed is 0, prob is 1-alpha
        else if (as.matrix(possible_true_genotype)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==0){error_result_mat[[k]][j,i]=1-alpha}
        
        # if true is 1 and observed is 0, false negative, the prob is beta
        else if (as.matrix(possible_true_genotype)[k,i]==1 & as.matrix(obs_genotype_mat)[j,i]==0){error_result_mat[[k]][j,i]=beta}
        
        # if true is 0 and observed is 1, false positive, the prob is alpha
        else if (as.matrix(possible_true_genotype)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==1){error_result_mat[[k]][j,i]=alpha}
        
        # if true is 1 and observed is 3, false negative, the prob is beta
        else if (as.matrix(possible_true_genotype)[k,i]==1 & as.matrix(obs_genotype_mat)[j,i]==3){error_result_mat[[k]][j,i]=1}
        
        # if true is 0 and observed is 3, false positive, the prob is alpha
        else if (as.matrix(possible_true_genotype)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==3){error_result_mat[[k]][j,i]=1}
      }
    }
  }
  
  
  length(error_result_mat)
  
  ##############################################################################
  ##############################################################################
  ##########
  # 
  
  # error_prob is a matrix of ncol= number of possible situation(ambguity), nrow= number of branches
  # error_prob, each column is one possible observed mutation, each line corresponds to the mutation branch
  error_prob=matrix(, nrow = dim(possible_true_genotype)[1], ncol = dim(obs_genotype_mat)[1])
  # branch k
  for (k in 1:dim(possible_true_genotype)[1]){
    # for situation j
    for (j in 1:dim(obs_genotype_mat)[1]){
      # find the product of the prob, error_prob[k,j] is the conditional prob of observed mutation when true mutation on branch(row) k, and ambiguity situation(column) j
      error_prob[k,j]= tail(cumprod(error_result_mat[[k]][j,]),1)
    }
  }
  
  
  
  ###########################################################################
  ###########################################################################
  
  ###########################################################################
  
  # weight is assigned to each possible situation(ambguity), and the total weighted prob is calculated
  # weighted_error is the probability of observed conditioning on a mutation branch, each row represents one branch
  weighted_error=error_prob%*%ambiguity_weight
  # all_weighted_error is the prob of the observed mutations status, conditioning on each branch(each row represent the conditional prob on that branch)
  
  ############################################################################
  ############################################################################
  
  # all_weighted_error is the prob of mutation condition on branch, allprob_Bj_final is the prob of mutation on branch and not on other branches
  # take the sum of the prob will return the prob of mutation status
  prob_S=t(weighted_error)%*%prob_Bj_final# probability of S vectors=sum of prob of S and X
  
  #############################################################################
  #############################################################################
  #for gene t
  
  prob_Bj_S=c()
  #for branch k
  for (k in 1:length(weighted_error)){
    prob_Bj_S[k]=weighted_error[k]*prob_Bj_final[k]/prob_S
  }
  
  ###############################
  #print out the branch number that has the max value
  ##############################
  
  br_index=which.max(prob_Bj_S)
  
  
  
  newdata <- c(prob_Bj_S,rep(0,number_br-length(prob_Bj_S)))
  
  return(result=newdata)
}

