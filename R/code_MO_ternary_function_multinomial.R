#library(ape)
#library(phangorn)

####################################################################################################
####################################################################################################
# basic functions for computing the prior and likelihood
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

MO_ternary_multinomial <- function(alpha,beta,unit_theta,unit_gamma,unit_mu,initial_obs,t){

  ts = t
  num_rows = nrow(ts$edge) #number of edges
  num_cols = length(ts$tip.label) #number of tips
  num_node = 2*length(ts$tip.label)-1 #number of nodes


  alpha01= alpha
  alpha02= alpha*beta/2
  beta10= beta/2
  beta12= beta/2
  sequencing_error_model=matrix(c(1-alpha01-alpha02,alpha01,alpha02,
                                  beta10,1-beta10-beta12,beta12,
                                  gamma20,gamma21,1-gamma21-gamma21),nrow=3,byrow = TRUE)
  #initialize to store the processed data. Double check the data and find out the site with missing values, ambiguous values
  obs_genotype = initial_obs[ts$tip.label]
  tip_exclude_index = which(obs_genotype == 3)
  tip_exclude = colnames(initial_obs)[tip_exclude_index]

  # find out all descendant br with above function
  descendant_branches_all=find_all_child_branches(ts$edge)

  # find out br with missing desendant tips
  br_exclude=c()
  for (Numnode in 1:num_node){
    Numnode_des=Descendants(ts, Numnode, type = c("tips"))
    if(all(Numnode_des[[1]] %in% tip_exclude_index)){single_br_exclude=which(ts$edge[,2] == Numnode)
    single_br_exclude_all = c(single_br_exclude,descendant_branches_all[[single_br_exclude]])
    br_exclude=c(br_exclude,single_br_exclude_all)}
  }

  # branch_time_list is the branch length of sub tree
  subtree=ts
  subtree$edge.length[br_exclude]=0
  branch_time_list= subtree$edge.length

  #######################################################################

  #generate the obs matrix
  # obs_colnam are those with observations
  obs_colnam=colnames(initial_obs)
  # consider the ambguity status as missing

  obs_genotype_mat=matrix(,nrow=2^length(grep("m",obs_genotype)),ncol=length(obs_colnam))
  colnames(obs_genotype_mat)=obs_colnam
  # find the index of each gene status for each sample
  ambiguity_index=which(obs_genotype=="m")
  normal_index=which(obs_genotype=="0")
  single_allele_index=which(obs_genotype=="1")
  double_allele_index=which(obs_genotype=="2")
  missing_index=which(obs_genotype=="3")
  #create all possible situations for ambguity
  inupt_list <- rep(list(0:2), length(ambiguity_index))
  input_ambiguity=expand.grid(inupt_list)
  # put the possible status into the matrix for gene i, each row represent one possible situation
  obs_genotype_mat[,as.numeric(ambiguity_index)]=as.matrix(input_ambiguity)
  obs_genotype_mat[,normal_index]=0
  obs_genotype_mat[,single_allele_index]=1
  obs_genotype_mat[,double_allele_index]=2
  obs_genotype_mat[,missing_index]=3
  # for each of the possible situation, assign weight to them, here, I use equal weights
  ambiguity_weight=matrix(rep(1/dim(obs_genotype_mat)[1],dim(obs_genotype_mat)[1],nrow=1))

  ###################################################################################
  ###################################################################################
  #find the joint prob of observation and branch_i in the subtree
  ###################################################################################
  ###################################################################################

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

  ################################################################################

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

  #true genotype when there is only one mutation on a branch
  possible_true_genotype_with_only_1 = cbind("First_branch" = 1:num_rows, "Second_branch" = rep(0,num_rows), possible_true_genotype)
  #true genotype when there are only two mutations on a branch
  possible_true_genotype_with_only_2 = cbind("First_branch" = 1:num_rows, "Second_branch" = rep(0,num_rows), possible_true_genotype*2)

  #true genotype when there are only two mutations on different branches
  num_of_cases = num_rows-1
  possible_true_genotype_with_2 = matrix(rep(0,num_of_cases*(num_cols+2)),nrow=num_of_cases,ncol = num_cols+2)
  colnames(possible_true_genotype_with_2)= c("First_branch", "Second_branch", subtree$tip.label)

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


  possible_true_genotype_with_2_sub = matrix(,nrow = 0,ncol = num_cols+2)
  colnames(possible_true_genotype_with_2_sub) = c("First_branch", "Second_branch", subtree$tip.label)
  for(branch_i in 1:num_rows){

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

  possible_true_genotype_with_2_sub = rbind(possible_true_genotype_with_2_sub,
                                            possible_true_genotype_with_2[ possible_true_genotype_with_2[,2] %in% descendant_branches_all[[branch_i]], ])

  }




  ###################################################################################
  ###################################################################################
  #Mutation model:find the prob of mutation on each branch, but not on other branches
  ###################################################################################
  ###################################################################################


  ##########################################################################################
  #only one mutation on one branch
  ##########################################################################################
  # prob_Bj_ind is vector storig the prob that mutation occurs on each branch
  prob_Bj_ind_0_1 = unit_theta*(exp(-(unit_theta+unit_gamma)*branch_time_list)-exp(-unit_mu*branch_time_list))/(unit_mu-unit_theta-unit_gamma)

  # create a matrix. In each row, diagnonal is the prob on that branch, and off-branch is the mutation prob that not on that branch
  prob_Bj_mat_0_1 = matrix(rep(exp(-(unit_theta+unit_gamma)*(branch_time_list)),length(branch_time_list)),nrow=length(branch_time_list), byrow = TRUE)

  # diagonal is replaced with prob of mutation
  diag(prob_Bj_mat_0_1) = prob_Bj_ind_0_1

  # descendant branches carry the mutation so the probability is 1
  for (l in 1:length(branch_time_list)){
    prob_Bj_mat_0_1[l,c(descendant_branches_all[[l]])] = exp(-unit_mu*branch_time_list[c(descendant_branches_all[[l]])])}

  # find the marginal prob that a mutation on each br
  prob_Bj_final_all_0_1 = t(apply(data.frame(prob_Bj_mat_0_1),1,cumprod))
  prob_Bj_final_0_1 = prob_Bj_final_all_0_1[,ncol(prob_Bj_final_all_0_1)]

  prob_Bj_0_1 = cbind(First_branch = 1:num_rows,
                              Second_branch = rep(0,num_rows),
                              Mutation_Prob=prob_Bj_final_0_1)



  ##########################################################################################
  #only two mutations on one branch
  ##########################################################################################
  prob_Bj_ind_0_2 = ((unit_gamma-unit_mu)*(exp(-(unit_theta+unit_gamma)*branch_time_list))+(unit_theta)*(exp(-unit_mu*branch_time_list)))/(unit_mu-unit_theta-unit_gamma)+1
  # create a matrix. In each row, diagnonal is the prob on that branch, and off-branch is the mutation prob that not on that branch
  prob_Bj_mat_0_2 = matrix(rep(exp(-(unit_theta+unit_gamma)*(branch_time_list)),length(branch_time_list)),nrow=length(branch_time_list), byrow = TRUE)

  # diagonal is replaced with prob of mutation
  diag(prob_Bj_mat_0_2) = prob_Bj_ind_0_2

  # descendant branches carry the mutation so the probability is 1
  for (l in 1:length(branch_time_list)){
    prob_Bj_mat_0_2[l,c(descendant_branches_all[[l]])] = 1}

  # find the marginal prob that a mutation on each br
  prob_Bj_final_all_0_2 = t(apply(data.frame(prob_Bj_mat_0_2),1,cumprod))
  prob_Bj_final_0_2 = prob_Bj_final_all_0_2[,ncol(prob_Bj_final_all_0_2)]

  prob_Bj_0_2 = cbind(First_branch = 1:num_rows,
                      Second_branch = rep(0,num_rows),
                      Mutation_Prob=prob_Bj_final_0_2)

  ##########################################################################################
  #first mutation on a branch, the second mutation on its descendant branch
  ##########################################################################################
  # results from function find_mutation_prob_Bj: prob_0_1= prob_Bj_with_only_1, prob_0_2 = prob_Bj_with_only_2, prob_0_1_2 = prob_Bj_with_2
  # find all branches numbers in the subtree
  all_branches=c(1:dim(subtree$edge)[1])

  prob_Bj_0_1_2 = matrix(,nrow=0,ncol = 3)
  colnames(prob_Bj_0_1_2) = c("First_branch","Second_branch","Mutation_Prob")
  for (branch_i in 1:num_rows){
    ###############################################
    #find the prob of only 1 mutation on branch_i
    ###############################################


    ########################################################################
    #find the prob of 0-1 mutation on branch_i and 1-2 on a descendant branch
    ########################################################################
    Bj_with_2_descendant_branches=descendant_branches_all[[branch_i]]
    if(length(Bj_with_2_descendant_branches) > 0){
      # create the matrix that store the prob if there are separate mutations on two branches
      prob_Bj_with_2_mat=matrix(, nrow = length(Bj_with_2_descendant_branches), ncol = 2+dim(subtree$edge)[1])

      # iterate on each of the descendant branch
      for (descendant_i in 1:length(Bj_with_2_descendant_branches)){

        prob_Bj_with_2_mat[descendant_i,1]= branch_i
        prob_Bj_with_2_mat[descendant_i,2]= Bj_with_2_descendant_branches[descendant_i]

        DB_descendant_i=descendant_branches_all[[Bj_with_2_descendant_branches[descendant_i]]]

        # find the branches without mutation
        Bj_with_2_no_mutation_branches = setdiff(all_branches,union(branch_i, Bj_with_2_descendant_branches))
        prob_Bj_with_2_no_mutation_branches = exp(-(unit_theta+unit_gamma)*(branch_time_list[Bj_with_2_no_mutation_branches]))

        # the branch with the first mutation occurring
        Bj_with_2_1st_mutation_branches = branch_i
        prob_Bj_with_2_1st_mutation_branches = unit_theta*(exp(-(unit_theta+unit_gamma)*branch_time_list[branch_i])-exp(-unit_mu*branch_time_list[branch_i]))/(unit_mu-unit_theta-unit_gamma)

        # the branch with the second mutation occurring
        Bj_with_2_2nd_mutation_branches = Bj_with_2_descendant_branches[descendant_i]
        prob_Bj_with_2_2nd_mutation_branches = 1-exp(-unit_mu*branch_time_list[Bj_with_2_2nd_mutation_branches])

        #the branches carry the first and second mutation
        Bj_with_2_carry_1st_2nd_mutation_branches = descendant_branches_all[[Bj_with_2_2nd_mutation_branches]]
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
        prob_Bj_with_2=cbind("First_branch"=prob_Bj_with_2_mat[,1],"Second_branch"=prob_Bj_with_2_mat[,2],
                             Mutation_Prob=t(apply(prob_Bj_with_2_mat[,-c(1,2)],1,cumprod))[,dim(subtree$edge)[1]])

      }

    }else{next}

    prob_Bj_0_1_2 = rbind(prob_Bj_0_1_2,prob_Bj_with_2)
  }



  ###################################################################################
  ###################################################################################
  #Error model:find the prob of observation, conditioning on a branch_i
  ###################################################################################
  ###################################################################################


    ##################################################################
    ###find the error prob when only 1 mutation on branch_i
    ##################################################################

    possible_true_genotype_with_only_1_mutation_values = possible_true_genotype_with_only_1[,c(3:dim(possible_true_genotype_with_only_1)[2])]
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

          else if (as.matrix(possible_true_genotype_with_only_1_mutation_values)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==3)
          {error_result_mat_with_only_1_mutation[[k]][j,i]=1}

          else if (as.matrix(possible_true_genotype_with_only_1_mutation_values)[k,i]==1 & as.matrix(obs_genotype_mat)[j,i]==3)
          {error_result_mat_with_only_1_mutation[[k]][j,i]=1}

          else if (as.matrix(possible_true_genotype_with_only_1_mutation_values)[k,i]==2 & as.matrix(obs_genotype_mat)[j,i]==3)
          {error_result_mat_with_only_1_mutation[[k]][j,i]=1}
        }
      }
    }
    prob_error_result_0_1_mat = matrix(unlist(error_result_mat_with_only_1_mutation), ncol = dim(obs_genotype_mat)[2], byrow = TRUE)
    prob_error_result_0_1 = cbind(First_branch = 1:num_rows,
                                  Second_branch = rep(0,num_rows),
                                  Error_Prob = apply(prob_error_result_0_1_mat, 1, prod))
    ##################################################################
    ###find the error prob when only 2 mutation on branch_i
    ##################################################################

    possible_true_genotype_with_only_2_mutation_values = possible_true_genotype_with_only_2[,c(3:dim(possible_true_genotype_with_only_2)[2])]
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


          else if (as.matrix(possible_true_genotype_with_only_2_mutation_values)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==3)
          {error_result_mat_with_only_2_mutation[[k]][j,i]=1}

          else if (as.matrix(possible_true_genotype_with_only_2_mutation_values)[k,i]==1 & as.matrix(obs_genotype_mat)[j,i]==3)
          {error_result_mat_with_only_2_mutation[[k]][j,i]=1}

          else if (as.matrix(possible_true_genotype_with_only_2_mutation_values)[k,i]==2 & as.matrix(obs_genotype_mat)[j,i]==3)
          {error_result_mat_with_only_2_mutation[[k]][j,i]=1}
        }
      }
    }

    prob_error_result_0_2_mat = matrix(unlist(error_result_mat_with_only_2_mutation), ncol = dim(obs_genotype_mat)[2], byrow = TRUE)
    prob_error_result_0_2 = cbind(First_branch = 1:num_rows,
                                  Second_branch = rep(0,num_rows),
                                  Error_Prob = apply(prob_error_result_0_2_mat, 1, prod))



    #########################################################################################
    ###find the error prob when 1 mutation on branch_i and 1 mutation on the sub branch
    #########################################################################################

    possible_true_genotype_with_2_mutation_values = possible_true_genotype_with_2_sub[,c(3:dim(possible_true_genotype_with_2_sub)[2])]
    colnames(possible_true_genotype_with_2_mutation_values)=subtree$tip.label
    #create an empty matrix to store the errors for only 1 mutation on branch_i
    error_result_mat_with_2_mutation=replicate(dim(possible_true_genotype_with_2_sub)[1],
                                               matrix(rep(0,dim(obs_genotype_mat)[1]*(dim(obs_genotype_mat)[2])),
                                                      nrow=dim(obs_genotype_mat)[1]),
                                               simplify=FALSE)


    #assign
    for (k in 1:dim(possible_true_genotype_with_2_mutation_values)[1]){
      #assign tipnames to the matrix
      colnames(error_result_mat_with_2_mutation[[k]]) = subtree$tip.label
      # situation j (possible situation for ambguity status),dim(obs_genotype_mat)[1]=2^(# of ambiguity sites)
      for (j in 1:dim(obs_genotype_mat)[1]){
        # tip or sample i
        for (i in 1:num_cols){
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


          # if true is 0 and observed is 2, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_2_mutation_values)[k,i]==0 & as.matrix(obs_genotype_mat)[j,i]==3)
          {error_result_mat_with_2_mutation[[k]][j,i]=1}

          # if true is 0 and observed is 2, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_2_mutation_values)[k,i]==1 & as.matrix(obs_genotype_mat)[j,i]==3)
          {error_result_mat_with_2_mutation[[k]][j,i]=1}

          # if true is 0 and observed is 2, false positive, the prob is alpha
          else if (as.matrix(possible_true_genotype_with_2_mutation_values)[k,i]==2 & as.matrix(obs_genotype_mat)[j,i]==3)
          {error_result_mat_with_2_mutation[[k]][j,i]=1}
        }
      }


    }

    prob_error_result_0_1_2_mat = matrix(unlist(error_result_mat_with_2_mutation), ncol = dim(obs_genotype_mat)[2], byrow = TRUE)
    prob_error_result_0_1_2 = cbind(First_branch = possible_true_genotype_with_2_sub[,1],
                              Second_branch = possible_true_genotype_with_2_sub[,2],
                              Error_Prob = apply(prob_error_result_0_1_2_mat, 1, prod))

  ###################################################################################
  ###################################################################################
  #posterior probability
  ###################################################################################
  ###################################################################################

  if(all(prob_error_result_0_1[,1] == prob_Bj_0_1[,1]) & all(prob_error_result_0_1[,2] == prob_Bj_0_1[,2])){

    posterior_mat_1 = cbind(prob_error_result_0_1[,1:2],posterprob = prob_error_result_0_1[,3] * prob_Bj_0_1[,3])

  }


  if(all(prob_error_result_0_2[,1] == prob_Bj_0_2[,1]) & all(prob_error_result_0_2[,2] == prob_Bj_0_2[,2])){

      posterior_mat_2 = cbind(prob_error_result_0_2[,1:2],posterprob = prob_error_result_0_2[,3] * prob_Bj_0_2[,3])

  }

    prob_Bj_0_1_2_sort = prob_Bj_0_1_2[order(prob_Bj_0_1_2[,c(1)],prob_Bj_0_1_2[,c(2)]),]
    prob_error_result_0_1_2_sort = prob_error_result_0_1_2[order(prob_error_result_0_1_2[,c(1)], prob_error_result_0_1_2[,c(2)]),]

  if(all(prob_error_result_0_1_2_sort[,1] == prob_Bj_0_1_2_sort[,1]) & all(prob_error_result_0_1_2_sort[,2] == prob_Bj_0_1_2_sort[,2])){

    posterior_mat_0_1_2 = cbind(prob_error_result_0_1_2_sort[,1:2],posterprob = prob_error_result_0_1_2_sort[,3] * prob_Bj_0_1_2_sort[,3])

    }


    posterior_mat = data.frame(rbind(posterior_mat_1,posterior_mat_2,posterior_mat_0_1_2))

    posterior_mat_joint = c()
    for(singleBr in 1:num_rows){

      posterior_mat_sub = posterior_mat_0_1_2[which(posterior_mat_0_1_2[,1] == singleBr),]
      posterior_mat_joint[singleBr] = sum(posterior_mat_sub[,3])

    }
    posterior_mat_comb = data.frame(cbind(posterior_mat_1[,3],posterior_mat_2[,3],posterior_mat_joint),(posterior_mat_1[,3]+posterior_mat_2[,3]+posterior_mat_joint))

  return(result = posterior_mat_comb[,4]/sum(posterior_mat_comb[,4]))

}
