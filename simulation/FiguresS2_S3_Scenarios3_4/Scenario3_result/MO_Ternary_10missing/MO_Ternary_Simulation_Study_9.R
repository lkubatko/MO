generate_prob <- function(sequencing_error_model,unit_theta,unit_gamma,unit_mu,number_br,
                          number_cell,normal_genotype,mutation_genotype,initial_obs,ts){
  
  # for this gene or location, variable singletip_exclude are those sample equal '-', variable tip_exclude are all those equal '-'
  
  single_tip_exclude = c()
  tip_exclude = c()
  
  # for this gene or location, obs_colnam are those sample names, excluding the '-'
  obs_colnam = c()
  obs_genotype = c()
  
  # for this gene or location, change the ACTG to 0/1/2, assign the mutation status for sample j 
  for (j in c(1:(number_cell))){
    # exclude those tips with missing gene status
    if (initial_obs[j]==c("-")) { #single_tip_exclude=colnames(initial_obs)[1,j] 
      single_tip_exclude=colnames(initial_obs)[j]
      tip_exclude=c(tip_exclude,single_tip_exclude)
    }
    # value is 0 if gene status is same as normal
    else if (as.character(initial_obs[1,j])==as.character(normal_genotype)) {
      obs_genotype=c(obs_genotype,0)
    }
    # value is 1 if gene status is same as mutant
    else if (as.character(initial_obs[1,j])==as.character(mutation_genotype)) {
      #obs_genotype=c(obs_genotype,1)
      obs_genotype=c(obs_genotype,2)
    }
    # value is m if gene status is ambguity
    else { #obs_genotype=c(obs_genotype,"m")
      obs_genotype=c(obs_genotype,1)}
  }
  
  # for this gene or location, exclude the sample with missing gene status
  subtree = drop.tip(ts, tip = tip_exclude)
  # branch_time_list is the branch length of sub tree
  branch_time_list = subtree$edge.length
  #generate the obs matrix
  # obs_colnam are those with observations
  obs_colnam=setdiff(colnames(initial_obs),tip_exclude)
  # consider the ambguity status as missing
  #obs_genotype_mat=matrix(,nrow=max(1,2^(count(obs_genotype=="m")[2,2])),ncol=length(obs_colnam))
  #obs_genotype_mat=matrix(,nrow=1,ncol=length(obs_colnam))
  obs_genotype_mat=matrix(,nrow=2^length(grep("m",obs_genotype)),ncol=length(obs_colnam))
  colnames(obs_genotype_mat)=obs_colnam
  # find the index of each gene status for each sample
  ambiguity_index=which(obs_genotype=="m")
  normal_index=which(obs_genotype=="0")
  single_allele_index=which(obs_genotype=="1")
  double_allele_index=which(obs_genotype=="2")
  #create all possible situations for ambguity 
  inupt_list <- rep(list(0:1), length(ambiguity_index))
  input_ambiguity=expand.grid(inupt_list)
  # put the possible status into the matrix for gene i, each row represent one possible situation
  obs_genotype_mat[,as.numeric(ambiguity_index)]=as.matrix(input_ambiguity)
  obs_genotype_mat[,normal_index]=rep(0,dim(obs_genotype_mat)[1])
  obs_genotype_mat[,single_allele_index]=rep(1,dim(obs_genotype_mat)[1])
  obs_genotype_mat[,double_allele_index]=rep(2,dim(obs_genotype_mat)[1])
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
  
  ###################################################################################
  ###################################################################################
  #Mutation model:find the prob of mutation on each branch, but not on other branches
  ###################################################################################
  ###################################################################################
  
  # results from function find_mutation_prob_Bj: prob_0_1= prob_Bj_with_only_1, prob_0_2 = prob_Bj_with_only_2, prob_0_1_2 = prob_Bj_with_2
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
                              Second_branch = 0,
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
  
  
  ###################################################################################
  ###################################################################################
  #Error model:find the prob of observation, conditioning on a branch_i
  ###################################################################################
  ###################################################################################
  
  # branch k (mutation on branch k)
  find_mutation_prob_of_obs_condition_Bj <- function(branch_time_list, sequencing_error_model,branch_i){
    
    
    possible_true_genotype_with_only_1_mutation=find_all_possible_mutation_matrix(subtree,branch_i)[1][[1]]
    possible_true_genotype_with_only_2_mutation=find_all_possible_mutation_matrix(subtree,branch_i)[2][[1]]
    possible_true_genotype_with_2_mutation=find_all_possible_mutation_matrix(subtree,branch_i)[3][[1]]
    
    
    ##################################################################
    ###find the conditional prob when only 1 mutation on branch_i
    ##################################################################
    
    possible_true_genotype_with_only_1_mutation_values = matrix(possible_true_genotype_with_only_1_mutation[,c(3:dim(possible_true_genotype_with_only_1_mutation)[2])],
                                                                nrow=nrow(possible_true_genotype_with_only_1_mutation),
                                                                byrow=T)
    colnames(possible_true_genotype_with_only_1_mutation_values)=subtree$tip.label
    #create an empty matrix to store the errors for only 1 mutation on branch_i
    error_result_mat_with_only_1_mutation=replicate(dim(possible_true_genotype_with_only_1_mutation_values)[1], 
                                                    matrix(rep(0,dim(obs_genotype_mat)[1]*(dim(obs_genotype_mat)[2])),
                                                           nrow=dim(obs_genotype_mat)[1]),
                                                    simplify=FALSE)
    
    
    #assign 
    for (k in 1:dim(possible_true_genotype_with_only_1_mutation_values)[1]){
      #assign tipnames to the matrix
      colnames(error_result_mat_with_only_1_mutation[[k]]) = subtree$tip.label
      # situation j (possible situation for ambguity status),dim(obs_genotype_mat)[1]=2^(# of ambiguity sites)
      for (j in 1:dim(obs_genotype_mat)[1]){
        # tip or sample i 
        for (i in 1:dim(possible_true_genotype_with_only_1_mutation_values)[2]){
          # if true is 1, and observed is 1, prob is 1-beta
          if (as.matrix(possible_true_genotype_with_only_1_mutation_values)[k,i]==1 &  as.matrix(obs_genotype_mat)[j,i]==1)
          {error_result_mat_with_only_1_mutation[[k]][j,i]=sequencing_error_model[2,2]}
          
          # if true is 0 and observed is 0, prob is 1-alpha
          else if (as.matrix(possible_true_genotype_with_only_1_mutation_values)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==0)
          {error_result_mat_with_only_1_mutation[[k]][j,i]=sequencing_error_model[1,1]}
          
          # if true is 1 and observed is 0, false negative, the prob is beta
          else if (as.matrix(possible_true_genotype_with_only_1_mutation_values)[k,i]==1 & as.matrix(obs_genotype_mat)[j,i]==0)
          {error_result_mat_with_only_1_mutation[[k]][j,i]=sequencing_error_model[2,1]}
          
          # if true is 0 and observed is 1, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_only_1_mutation_values)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==1)
          {error_result_mat_with_only_1_mutation[[k]][j,i]=sequencing_error_model[1,2]}
          
          # if true is 2 and observed is 2, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_only_1_mutation_values)[k,i]==2 & as.matrix(obs_genotype_mat)[j,i]==2)
          {error_result_mat_with_only_1_mutation[[k]][j,i]=sequencing_error_model[3,3]}
          
          # if true is 2 and observed is 1, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_only_1_mutation_values)[k,i]==2 & as.matrix(obs_genotype_mat)[j,i]==1)
          {error_result_mat_with_only_1_mutation[[k]][j,i]=sequencing_error_model[3,2]}
          
          # if true is 2 and observed is 0, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_only_1_mutation_values)[k,i]==2 & as.matrix(obs_genotype_mat)[j,i]==0)
          {error_result_mat_with_only_1_mutation[[k]][j,i]=sequencing_error_model[3,1]}
          
          # if true is 1 and observed is 2, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_only_1_mutation_values)[k,i]==1 & as.matrix(obs_genotype_mat)[j,i]==2)
          {error_result_mat_with_only_1_mutation[[k]][j,i]=sequencing_error_model[2,3]}
          
          # if true is 0 and observed is 2, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_only_1_mutation_values)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==2)
          {error_result_mat_with_only_1_mutation[[k]][j,i]=sequencing_error_model[1,3]}
          
        }
      }
    }
    
    ##################################################################
    ###find the conditional prob when only 2 mutation on branch_i
    ##################################################################
    
    possible_true_genotype_with_only_2_mutation_values = matrix(possible_true_genotype_with_only_2_mutation[,c(3:dim(possible_true_genotype_with_only_2_mutation)[2])],
                                                                nrow=nrow(possible_true_genotype_with_only_2_mutation),
                                                                byrow=T)
    colnames(possible_true_genotype_with_only_2_mutation_values)=subtree$tip.label
    #create an empty matrix to store the errors for only 1 mutation on branch_i
    error_result_mat_with_only_2_mutation=replicate(dim(possible_true_genotype_with_only_2_mutation_values)[1], 
                                                    matrix(rep(0,dim(obs_genotype_mat)[1]*(dim(obs_genotype_mat)[2])),
                                                           nrow=dim(obs_genotype_mat)[1]),
                                                    simplify=FALSE)
    
    
    #assign 
    for (k in 1:dim(possible_true_genotype_with_only_2_mutation_values)[1]){
      #assign tipnames to the matrix
      colnames(error_result_mat_with_only_2_mutation[[k]]) = subtree$tip.label
      # situation j (possible situation for ambguity status),dim(obs_genotype_mat)[1]=2^(# of ambiguity sites)
      for (j in 1:dim(obs_genotype_mat)[1]){
        # tip or sample i 
        for (i in 1:dim(possible_true_genotype_with_only_2_mutation_values)[2]){
          # if true is 1, and observed is 1, prob is 1-beta
          if (as.matrix(possible_true_genotype_with_only_2_mutation_values)[k,i]==1 &  as.matrix(obs_genotype_mat)[j,i]==1)
          {error_result_mat_with_only_2_mutation[[k]][j,i]=sequencing_error_model[2,2]}
          
          # if true is 0 and observed is 0, prob is 1-alpha
          else if (as.matrix(possible_true_genotype_with_only_2_mutation_values)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==0)
          {error_result_mat_with_only_2_mutation[[k]][j,i]=sequencing_error_model[1,1]}
          
          # if true is 1 and observed is 0, false negative, the prob is beta
          else if (as.matrix(possible_true_genotype_with_only_2_mutation_values)[k,i]==1 & as.matrix(obs_genotype_mat)[j,i]==0)
          {error_result_mat_with_only_2_mutation[[k]][j,i]=sequencing_error_model[2,1]}
          
          # if true is 0 and observed is 1, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_only_2_mutation_values)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==1)
          {error_result_mat_with_only_2_mutation[[k]][j,i]=sequencing_error_model[1,2]}
          
          # if true is 2 and observed is 2, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_only_2_mutation_values)[k,i]==2 & as.matrix(obs_genotype_mat)[j,i]==2)
          {error_result_mat_with_only_2_mutation[[k]][j,i]=sequencing_error_model[3,3]}
          
          # if true is 2 and observed is 1, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_only_2_mutation_values)[k,i]==2 & as.matrix(obs_genotype_mat)[j,i]==1)
          {error_result_mat_with_only_2_mutation[[k]][j,i]=sequencing_error_model[3,2]}
          
          # if true is 2 and observed is 0, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_only_2_mutation_values)[k,i]==2 & as.matrix(obs_genotype_mat)[j,i]==0)
          {error_result_mat_with_only_2_mutation[[k]][j,i]=sequencing_error_model[3,1]}
          
          # if true is 1 and observed is 2, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_only_2_mutation_values)[k,i]==1 & as.matrix(obs_genotype_mat)[j,i]==2)
          {error_result_mat_with_only_2_mutation[[k]][j,i]=sequencing_error_model[2,3]}
          
          # if true is 0 and observed is 2, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_only_2_mutation_values)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==2)
          {error_result_mat_with_only_2_mutation[[k]][j,i]=sequencing_error_model[1,3]}
          
        }
      }
    }
    
    
    
    
    #########################################################################################
    ###find the conditional prob when 1 mutation on branch_i and 1 mutation on the sub branch
    #########################################################################################
    
    possible_true_genotype_with_2_mutation_values = matrix(possible_true_genotype_with_2_mutation[,c(3:dim(possible_true_genotype_with_2_mutation)[2])],
                                                           nrow=nrow(possible_true_genotype_with_2_mutation),
                                                           byrow=T)
    if(dim(possible_true_genotype_with_2_mutation_values)[1] > 0)
    {colnames(possible_true_genotype_with_2_mutation_values)=subtree$tip.label
    #create an empty matrix to store the errors for only 1 mutation on branch_i
    error_result_mat_with_2_mutation=replicate(dim(possible_true_genotype_with_2_mutation)[1], 
                                               matrix(rep(0,dim(obs_genotype_mat)[1]*(dim(obs_genotype_mat)[2])),
                                                      nrow=dim(obs_genotype_mat)[1]),
                                               simplify=FALSE)
    prob_error_result_0_1_2 = c()
    
    #assign 
    for (k in 1:dim(possible_true_genotype_with_2_mutation_values)[1]){
      #assign tipnames to the matrix
      colnames(error_result_mat_with_2_mutation[[k]]) = subtree$tip.label
      # situation j (possible situation for ambguity status),dim(obs_genotype_mat)[1]=2^(# of ambiguity sites)
      for (j in 1:dim(obs_genotype_mat)[1]){
        # tip or sample i 
        for (i in 1:dim(possible_true_genotype_with_2_mutation_values)[2]){
          # if true is 1, and observed is 1, prob is 1-beta
          if (as.matrix(possible_true_genotype_with_2_mutation_values)[k,i]==1 &  as.matrix(obs_genotype_mat)[j,i]==1)
          {error_result_mat_with_2_mutation[[k]][j,i]=sequencing_error_model[2,2]}
          
          # if true is 0 and observed is 0, prob is 1-alpha
          else if (as.matrix(possible_true_genotype_with_2_mutation_values)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==0)
          {error_result_mat_with_2_mutation[[k]][j,i]=sequencing_error_model[1,1]}
          
          # if true is 1 and observed is 0, false negative, the prob is beta
          else if (as.matrix(possible_true_genotype_with_2_mutation_values)[k,i]==1 & as.matrix(obs_genotype_mat)[j,i]==0)
          {error_result_mat_with_2_mutation[[k]][j,i]=sequencing_error_model[2,1]}
          
          # if true is 0 and observed is 1, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_2_mutation_values)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==1)
          {error_result_mat_with_2_mutation[[k]][j,i]=sequencing_error_model[1,2]}
          
          # if true is 2 and observed is 2, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_2_mutation_values)[k,i]==2 & as.matrix(obs_genotype_mat)[j,i]==2)
          {error_result_mat_with_2_mutation[[k]][j,i]=sequencing_error_model[3,3]}
          
          # if true is 2 and observed is 1, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_2_mutation_values)[k,i]==2 & as.matrix(obs_genotype_mat)[j,i]==1)
          {error_result_mat_with_2_mutation[[k]][j,i]=sequencing_error_model[3,2]}
          
          # if true is 2 and observed is 0, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_2_mutation_values)[k,i]==2 & as.matrix(obs_genotype_mat)[j,i]==0)
          {error_result_mat_with_2_mutation[[k]][j,i]=sequencing_error_model[3,1]}
          
          # if true is 1 and observed is 2, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_2_mutation_values)[k,i]==1 & as.matrix(obs_genotype_mat)[j,i]==2)
          {error_result_mat_with_2_mutation[[k]][j,i]=sequencing_error_model[2,3]}
          
          # if true is 0 and observed is 2, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_2_mutation_values)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==2)
          {error_result_mat_with_2_mutation[[k]][j,i]=sequencing_error_model[1,3]}
          
        }
      }
      
      prob_error_result_0_1_2 [[k]] = tail(cumprod(error_result_mat_with_2_mutation[[k]]), n = 1)
    }
    } else
    {error_result_mat_with_2_mutation =0
    prob_error_result_0_1_2 = 0}
    # find the 
    
    prob_error_result_0_1 = tail(cumprod(error_result_mat_with_only_1_mutation[[1]]),n=1)
    prob_error_result_0_2 = tail(cumprod(error_result_mat_with_only_2_mutation[[1]]),n=1)
    ##################################################################
    #return the errors for the three conditions
    ##################################################################
    
    return(list(error_result_0_1=error_result_mat_with_only_1_mutation, 
                error_result_0_2=error_result_mat_with_only_2_mutation, 
                error_result_0_1_2=error_result_mat_with_2_mutation,
                prob_error_result_0_1=prob_error_result_0_1,
                prob_error_result_0_2=prob_error_result_0_2,
                prob_error_result_0_1_2=prob_error_result_0_1_2))
    
    
  }
  
  
  
  ##############################################################################
  ##############################################################################
  find_joint_obs_Bj <- function(subtree, branch_i)
  {
    #use the mutation function to find the probability
    prob_Bj_final = find_mutation_prob_Bj(branch_time_list, unit_theta,unit_gamma,unit_mu,branch_i)
    #use the error function to find the conditional probability
    prob_find_mutation_prob_of_obs_condition_Bj = find_mutation_prob_of_obs_condition_Bj(branch_time_list, sequencing_error_model,branch_i)
    
    prob_Bj_and_obs_genotype = prob_Bj_final$prob_0_1[1,3] *  prob_find_mutation_prob_of_obs_condition_Bj$prob_error_result_0_1 +
      prob_Bj_final$prob_0_2[1,3] *  prob_find_mutation_prob_of_obs_condition_Bj$prob_error_result_0_2 +
      t(as.matrix(prob_Bj_final$prob_0_1_2[,3])) %*% (as.matrix(prob_find_mutation_prob_of_obs_condition_Bj$prob_error_result_0_1_2))
    
    prob_Bj_and_obs_genotype_0_1 = prob_Bj_final$prob_0_1[1,3] *  prob_find_mutation_prob_of_obs_condition_Bj$prob_error_result_0_1
    prob_Bj_and_obs_genotype_0_2 = prob_Bj_final$prob_0_2[1,3] *  prob_find_mutation_prob_of_obs_condition_Bj$prob_error_result_0_2
    prob_Bj_and_obs_genotype_0_1_2 = t(as.matrix(prob_Bj_final$prob_0_1_2[,3])) %*% (as.matrix(prob_find_mutation_prob_of_obs_condition_Bj$prob_error_result_0_1_2))
    
    return(list(prob_Bj_and_obs_genotype = prob_Bj_and_obs_genotype, 
                prob_Bj_and_obs_genotype_0_1 = prob_Bj_and_obs_genotype_0_1, 
                prob_Bj_and_obs_genotype_0_2 = prob_Bj_and_obs_genotype_0_2,
                prob_Bj_and_obs_genotype_0_1_2 = prob_Bj_and_obs_genotype_0_1_2))
    
  }
  
  #find the joint prob of observation and branch_i on all branches of the subtree
  find_joint_obs_Bj_all_branches<-function(subtree){
    prob_joint_obs_Bj_all_branches_mat = matrix(,nrow=dim(subtree$edge)[1],ncol=4)
    colnames(prob_joint_obs_Bj_all_branches_mat) = c("prob_Bj_and_obs_genotype_0_1","prob_Bj_and_obs_genotype_0_2",
                                                     "prob_Bj_and_obs_genotype_0_1_2","prob_Bj_and_obs_genotype")
    for (br_num in 1:dim(subtree$edge)[1]){
      prob_find_joint_obs_Bj = find_joint_obs_Bj(subtree,br_num)
      
      prob_joint_obs_Bj_all_branches_mat[br_num,1] = prob_find_joint_obs_Bj$prob_Bj_and_obs_genotype_0_1
      prob_joint_obs_Bj_all_branches_mat[br_num,2] = prob_find_joint_obs_Bj$prob_Bj_and_obs_genotype_0_2
      prob_joint_obs_Bj_all_branches_mat[br_num,3] = prob_find_joint_obs_Bj$prob_Bj_and_obs_genotype_0_1_2
      prob_joint_obs_Bj_all_branches_mat[br_num,4] = prob_find_joint_obs_Bj$prob_Bj_and_obs_genotype
    }
    
    return(prob_joint_obs_Bj_all_branches_mat = prob_joint_obs_Bj_all_branches_mat)
  }
  
  find_mutation_prob_of_Bj_condition_obs<-function(subtree){
    prob_joint_obs_Bj_all_branches_mat = find_joint_obs_Bj_all_branches(subtree)
    prob_of_Bj_condition_obs_mat = matrix(,nrow=dim(subtree$edge)[1],ncol=4)
    colnames(prob_of_Bj_condition_obs_mat) = c("prob_of_Bj_condition_obs_0_1","prob_of_Bj_condition_obs_0_2",
                                               "prob_of_Bj_condition_obs_0_1_2","prob_of_Bj_condition_obs_genotype")
    for (br_num in 1:dim(subtree$edge)[1]){
      
      prob_of_Bj_condition_obs_mat[br_num,1] = prob_joint_obs_Bj_all_branches_mat[br_num,1]/sum(prob_joint_obs_Bj_all_branches_mat[,4])
      prob_of_Bj_condition_obs_mat[br_num,2] = prob_joint_obs_Bj_all_branches_mat[br_num,2]/sum(prob_joint_obs_Bj_all_branches_mat[,4])
      prob_of_Bj_condition_obs_mat[br_num,3] = prob_joint_obs_Bj_all_branches_mat[br_num,3]/sum(prob_joint_obs_Bj_all_branches_mat[,4])
      prob_of_Bj_condition_obs_mat[br_num,4] = prob_joint_obs_Bj_all_branches_mat[br_num,4]/sum(prob_joint_obs_Bj_all_branches_mat[,4])
    }
    return(prob_of_Bj_condition_obs_mat = prob_of_Bj_condition_obs_mat )
  }
  
  result = find_mutation_prob_of_Bj_condition_obs(subtree)
  
  return(result = result)
  
}


library(ape)

parameter_setting = expand.grid(alpha=c(0.05,0.1,0.2,0.4),
                                beta =c(0.05,0.1,0.2,0.4))
set.seed(1000)
for(paraInd in 9:9){
  
alpha =   parameter_setting[paraInd,1]
beta = parameter_setting[paraInd,2]
  




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


ternary_folder_form_result = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_10Missing/Ternary_alpha0%s_beta0%s_result',alpha_str, beta_str)
dir.create(ternary_folder_form_result)


for (indexn in 1:100){
  print(c(alpha,beta,indexn))
  
  form = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_10Missing/RandomTree/RandomTree_%s.tre', indexn)
  
  sampletr=read.tree(form)
  
  obs_form_0_1_2 = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_10Missing/Ternary_alpha0%s_beta0%s/ternary_obs_0_1_tip_alpha_0%s_beta_0%s_matrix%s.csv', alpha_str, beta_str, alpha_str, beta_str,indexn)
  
  
  mat_obs_form_0_1_2 = read.csv(obs_form_0_1_2)  
  
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
    initial_obs_0_1 = data.frame(mat_obs_form_0_1_2[,-c(1,2,3)])
    
    
    for (i in 1:dim(initial_obs_0_1)[1]){
      
      print(i)
      
      #rd_unit_theta <- rbeta(10, (10^7)*unit_theta, (10^7)*(1-unit_theta))
      #rd_unit_gamma <- rbeta(10, (10^14)*unit_gamma, (10^14)*(1-unit_gamma))
      #rd_unit_mu <- rbeta(10, 100*unit_mu, 100*(1-unit_mu))
      
      rd_unit_theta =  rgamma(n = 3, shape = 100, scale = 0.01*unit_theta)
      rd_unit_gamma = rgamma(3, shape = 100, scale = 0.01*unit_gamma)
      rd_unit_mu = rgamma(3, shape = 100, scale = 0.01*unit_mu)
      
      generate_prob_br_0_1_dat=data.frame(matrix(NA, nrow = number_br, ncol = 3))
      generate_prob_br_0_2_dat=data.frame(matrix(NA, nrow = number_br, ncol = 3))
      generate_prob_br_0_1_2_dat=data.frame(matrix(NA, nrow = number_br, ncol = 3))
      generate_prob_br_all_dat=data.frame(matrix(NA, nrow = number_br, ncol = 3))
      
      
      for (j in 1:3){
        
        generate_prob_br <- generate_prob(sequencing_error_model,rd_unit_theta[j],rd_unit_gamma[j],rd_unit_mu[j],number_br,number_cell,
                                          normal_genotype_0_1_2[i],mutation_genotype_0_1_2[i],initial_obs_0_1[i,],sampletr)
        generate_prob_br_0_1_single <- c(generate_prob_br[,1],rep(0,number_br-dim(generate_prob_br)[1]))
        generate_prob_br_0_2_single <- c(generate_prob_br[,2],rep(0,number_br-dim(generate_prob_br)[1]))
        generate_prob_br_0_1_2_single <- c(generate_prob_br[,3],rep(0,number_br-dim(generate_prob_br)[1]))
        generate_prob_br_all_single <- c(generate_prob_br[,4],rep(0,number_br-dim(generate_prob_br)[1]))
        
        generate_prob_br_0_1_dat[,j] = generate_prob_br_0_1_single
        generate_prob_br_0_2_dat[,j] = generate_prob_br_0_2_single
        generate_prob_br_0_1_2_dat[,j] = generate_prob_br_0_1_2_single
        generate_prob_br_all_dat[,j] = generate_prob_br_all_single
        
      }
      
      generate_prob_br_0_1=rowMeans(generate_prob_br_0_1_dat, na.rm = FALSE, dims = 1)
      generate_prob_br_0_2=rowMeans(generate_prob_br_0_2_dat, na.rm = FALSE, dims = 1)
      generate_prob_br_0_1_2=rowMeans(generate_prob_br_0_1_2_dat, na.rm = FALSE, dims = 1)
      generate_prob_br_all=rowMeans(generate_prob_br_all_dat, na.rm = FALSE, dims = 1)
      
      
      ternary_prob_matrix_all_0_1 = rbind(ternary_prob_matrix_all_0_1,generate_prob_br_all)
      ternary_prob_matrix_01_0_1 = rbind(ternary_prob_matrix_01_0_1,generate_prob_br_0_1)
      ternary_prob_matrix_02_0_1 = rbind(ternary_prob_matrix_02_0_1,generate_prob_br_0_2)
      ternary_prob_matrix_012_0_1 = rbind(ternary_prob_matrix_012_0_1,generate_prob_br_0_1_2)
      
      
    }
    
    ternary_prob_matrix_all_0_1_rownames=cbind(data.frame(mat_obs_form_0_1_2[,c(2,3)]),data.frame(ternary_prob_matrix_all_0_1))
    ternary_prob_matrix_01_0_1_rownames=cbind(data.frame(mat_obs_form_0_1_2[,c(2,3)]),data.frame(ternary_prob_matrix_01_0_1))
    ternary_prob_matrix_02_0_1_rownames=cbind(data.frame(mat_obs_form_0_1_2[,c(2,3)]),data.frame(ternary_prob_matrix_02_0_1))
    ternary_prob_matrix_012_0_1_rownames=cbind(data.frame(mat_obs_form_0_1_2[,c(2,3)]),data.frame(ternary_prob_matrix_012_0_1))
    
    
    ternary_prob_matrix_all_0_1_all=rbind(ternary_prob_matrix_all_0_1_all,ternary_prob_matrix_all_0_1_rownames)
    ternary_prob_matrix_01_0_1_all=rbind(ternary_prob_matrix_01_0_1_all,ternary_prob_matrix_01_0_1_rownames)
    ternary_prob_matrix_02_0_1_all=rbind(ternary_prob_matrix_02_0_1_all,ternary_prob_matrix_02_0_1_rownames)
    ternary_prob_matrix_012_0_1_all=rbind(ternary_prob_matrix_012_0_1_all,ternary_prob_matrix_012_0_1_rownames)
    
    
    
    
  }
  
  
  
  
  ternary_prob_matrix_all_0_1_2_out = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_10Missing/Ternary_alpha0%s_beta0%s_result/ternary_prob_matrix_all_0_1_2_out_alpha_0%s_beta_0%s_matrix%s.csv', alpha_str, beta_str, alpha_str, beta_str,indexn)
  ternary_prob_matrix_01_0_1_2_out = sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_10Missing/Ternary_alpha0%s_beta0%s_result/ternary_prob_matrix_01_0_1_2_out_alpha_0%s_beta_0%s_matrix%s.csv', alpha_str, beta_str, alpha_str, beta_str,indexn)
  ternary_prob_matrix_02_0_1_2_out= sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_10Missing/Ternary_alpha0%s_beta0%s_result/ternary_prob_matrix_02_0_1_2_out_alpha_0%s_beta_0%s_matrix%s.csv', alpha_str, beta_str, alpha_str, beta_str,indexn)
  ternary_prob_matrix_012_0_1_2_out= sprintf('/fs/project/kubatko.2-temp/gao.957/DNA_alignment/bowtie2-2.3.4.2-linux-x86_64/reference/Simulation_Setting_LargeTree/SimulateData_10Missing/Ternary_alpha0%s_beta0%s_result/ternary_prob_matrix_012_0_1_2_out_alpha_0%s_beta_0%s_matrix%s.csv', alpha_str, beta_str, alpha_str, beta_str,indexn)
  
  
  
  
  
  write.csv(ternary_prob_matrix_all_0_1_all,file=ternary_prob_matrix_all_0_1_2_out)
  write.csv(ternary_prob_matrix_01_0_1_all,file=ternary_prob_matrix_01_0_1_2_out)
  write.csv(ternary_prob_matrix_02_0_1_all,file=ternary_prob_matrix_02_0_1_2_out)
  write.csv(ternary_prob_matrix_012_0_1_all,file=ternary_prob_matrix_012_0_1_2_out)
}

}





