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
                                sites =c(20,40,80),
                                missingindex=c(0))

all_pair_result_dat= data.frame(matrix(, nrow=0, ncol=9))

for(paraInd in 1:dim(parameter_setting)[1]){
  
  
  alpha =   parameter_setting[paraInd,1]
  beta = parameter_setting[paraInd,1]
  numSites= parameter_setting[paraInd,2]
  missingindex = parameter_setting[paraInd,3]
  
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


all_pair_result= data.frame(matrix(, nrow=0, ncol=9))

for (indexn in 1:100){
  
  print(c(paraInd,indexn))
  
  if (indexn < 10){
    
    trueTree_form=sprintf("./FiguresS4_S5_Scenarios5_6/Scenario5/results_alpha_0%s_beta_0%s_%s/trees.000%s",alpha_str,alpha_str,numSites,indexn)
    
  }else if( indexn>=10 & indexn<100){
    
    trueTree_form=sprintf("./FiguresS4_S5_Scenarios5_6/Scenario5/results_alpha_0%s_beta_0%s_%s/trees.00%s",alpha_str,alpha_str,numSites,indexn)
    
    
  }else{
    
    trueTree_form=sprintf("./FiguresS4_S5_Scenarios5_6/Scenario5/results_alpha_0%s_beta_0%s_%s/trees.0%s",alpha_str,alpha_str,numSites,indexn)
    
  }
  
  
  scanned_trueTree_form = scan(file=trueTree_form,what=character(), n = -1, sep = "")
  
  sampletr = read.tree(text=scanned_trueTree_form)
  
  if (missingindex == 0){
    
    prob_form_0_1 = sprintf("./FiguresS4_S5_Scenarios5_6/Scenario5_result_rmd/MO_Binary/Binary_alpha0%s_beta0%s_%s_result/all_binary_prob_matrix_all_0_1_out_matrix%s.csv",alpha_str,alpha_str,numSites,indexn)
    
  }else{
    
    prob_form_0_1 = sprintf("./FiguresS4_S5_Scenarios5_6/Scenario5_result_rmd/MO_Binary/Binary_alpha0%s_beta0%s_%s_result/all_binary_prob_matrix_all_0_1_out_matrix%s_%s.csv",alpha_str,alpha_str,numSites,indexn,missingindex)
    
  }
  
  mat_prob_form_0_1 = read.csv(prob_form_0_1)
  
  true_near_branch = compare_near_order(length(sampletr$edge.length),length(sampletr$tip.label),sampletr)#all descendant branches of each branch
  
  all_true_near_branch = data.frame(matrix(, nrow=0, ncol=2))#all true possible pairs on true tree
  all_match_near_branch_sub = data.frame(matrix(, nrow=0, ncol=4))#inferred matched pair
  all_true_match_near_branch_sub = data.frame(matrix(, nrow=0, ncol=4))#true pairs
  
  #four columns include gene 1, gene 2, br of gene 1, br of gene 2
  colnames(all_match_near_branch_sub) = c("gene_ID1","gene_ID2","first_br_indexn","second_br_indexn")
  colnames(all_true_match_near_branch_sub) = c("gene_ID1","gene_ID2","first_br_indexn","second_br_indexn")
  
  for (br in 1:length(true_near_branch)){
    
    if(length(true_near_branch[[br]])>0){
      for(br_sub in 1:length(true_near_branch[[br]])){
        all_true_near_branch = rbind(all_true_near_branch,c(br,(true_near_branch[[br]])[br_sub]))
      }
    }
  }
  colnames(all_true_near_branch) = c("first_br_indexn","second_br_indexn")
  
  all_true_near_branch_sub_all = na.omit(rbind(all_true_near_branch,
                                               cbind(first_br_indexn=1:dim(sampletr$edge)[1],
                                                     second_br_indexn=1:dim(sampletr$edge)[1])))#all true possible pairs on true tree removing NA values
  
  colnames(all_true_near_branch_sub_all) = c("first_br_indexn","second_br_indexn")
  
 
  ###################################################################################################################### 
    #find matched inferred br order
    for (i in 1:(dim(mat_prob_form_0_1)[1]-1)){
      
      for (j in (i+1):dim(mat_prob_form_0_1)[1]){
        
        print(c(i,j))
        first_indexn= mat_prob_form_0_1$selected_br[i]#selected br of gene i
        second_indexn= mat_prob_form_0_1$selected_br[j]#selected br of gene j
        
      if(first_indexn == second_indexn & mat_prob_form_0_1$selected_First_branch[i] == mat_prob_form_0_1$selected_First_branch[j] & 
         mat_prob_form_0_1$selected_First_branch[i]!=0 & mat_prob_form_0_1$selected_First_branch[j]!=0){#find those on same br
        
        matchresult=data.frame(first_indexn,second_indexn)#selected mutation branch
        truematchresult = data.frame(mat_prob_form_0_1$selected_First_branch[i],mat_prob_form_0_1$selected_First_branch[j])#true mutation br
        gene_result=data.frame(i,j)
        
        colnames(matchresult)=c("first_br_indexn","second_br_indexn")
        colnames(truematchresult) = c("first_br_indexn","second_br_indexn")
        colnames(gene_result)=c("gene_ID1","gene_ID2")
        
      }else if(first_indexn != second_indexn & mat_prob_form_0_1$selected_First_branch[i] != mat_prob_form_0_1$selected_First_branch[j]){
        
       
        #selected mutation branch with gene i occurs first
        br_indexn_1=data.frame(first_indexn,second_indexn)#selected mutation br order of ij
        true_br_indexn_1=data.frame(mat_prob_form_0_1$selected_First_branch[i],mat_prob_form_0_1$selected_First_branch[j])#true mutation br order
        colnames(br_indexn_1)=c("first_br_indexn","second_br_indexn")
        colnames(true_br_indexn_1)=c("first_br_indexn","second_br_indexn")
        
        #selected mutation branch with gene j occurs first
        br_indexn_2=data.frame(second_indexn,first_indexn)#selected mutation br order of ji
        true_br_indexn_2=data.frame(mat_prob_form_0_1$selected_First_branch[j],mat_prob_form_0_1$selected_First_branch[i])#true mutation br order
        colnames(br_indexn_2)=c("first_br_indexn","second_br_indexn")
        colnames(true_br_indexn_2)=c("first_br_indexn","second_br_indexn")
        
        matchresult_1 = match_df(all_true_near_branch_sub_all, br_indexn_1)#compare selected mutation br
        truematchresult_1 = match_df(all_true_near_branch_sub_all, true_br_indexn_1)
        colnames(truematchresult_1) = c("first_br_indexn","second_br_indexn")
        
        matchresult_2 = match_df(all_true_near_branch_sub_all, br_indexn_2)#compare selected mutation br
        truematchresult_2 = match_df(all_true_near_branch_sub_all, true_br_indexn_2)
        colnames(truematchresult_2) = c("first_br_indexn","second_br_indexn")
        
        
        if(dim(matchresult_1)[1]==0 & dim(matchresult_2)[1]==1 & dim(truematchresult_1)==0 & dim(truematchresult_2)==1){
          
          matchresult=matchresult_2
          truematchresult = truematchresult_2
          gene_result=data.frame(i,j)
          colnames(gene_result)=c("gene_ID1","gene_ID2")
          colnames(matchresult)=c("first_br_indexn","second_br_indexn")
          colnames(truematchresult)=c("first_br_indexn","second_br_indexn")
          
        }else if(dim(matchresult_1)[1]==1 & dim(matchresult_2)[1]==0 & dim(truematchresult_1)==1 & dim(truematchresult_2)==0){
          
          matchresult=matchresult_1
          truematchresult = truematchresult_1
          gene_result=data.frame(i,j)
          colnames(gene_result)=c("gene_ID1","gene_ID2")
          colnames(matchresult)=c("first_br_indexn","second_br_indexn")
          colnames(truematchresult)=c("first_br_indexn","second_br_indexn")
          
        }else{
          matchresult=data.frame(matrix(, nrow=0, ncol=2))
          truematchresult = data.frame(matrix(, nrow=0, ncol=2))
          gene_result=data.frame(matrix(, nrow=0, ncol=2))
          colnames(gene_result)=c("gene_ID1","gene_ID2")
          colnames(matchresult)=c("first_br_indexn","second_br_indexn")
          colnames(truematchresult)=c("first_br_indexn","second_br_indexn")
          
        }
        
        
      }else{
        matchresult=data.frame(matrix(, nrow=0, ncol=2))
        truematchresult=data.frame(matrix(, nrow=0, ncol=2))
        gene_result=data.frame(matrix(, nrow=0, ncol=2))
        colnames(gene_result)=c("gene_ID1","gene_ID2")
        colnames(matchresult)=c("first_br_indexn","second_br_indexn")
        colnames(truematchresult)=c("first_br_indexn","second_br_indexn")
        
        
      }
        
      
      
    all_match_near_branch_sub=rbind(all_match_near_branch_sub, cbind(gene_result,matchresult))
      
      }
      
      }
    #find true br order
    for (i in 1:(dim(mat_prob_form_0_1)[1]-1)){
    
      for (j in (i+1):dim(mat_prob_form_0_1)[1]){
      
      gene_result=data.frame(i,j)

      first_br_indexn= mat_prob_form_0_1$selected_First_branch[i]#true br of gene i
      second_br_indexn= mat_prob_form_0_1$selected_First_branch[j]#true br of gene j
      if(first_br_indexn == second_br_indexn & first_br_indexn!=0 & second_br_indexn!=0){
        gene_result=data.frame(i,j)
        truematchresult=data.frame(first_br_indexn,second_br_indexn)
        
      }else{
        
        
        br_indexn_1=data.frame(first_br_indexn,second_br_indexn)
        br_indexn_2=data.frame(second_br_indexn,first_br_indexn)
        colnames(br_indexn_2) = colnames(br_indexn_1)
        
        matchresult_1 = match_df(all_true_near_branch_sub_all, br_indexn_1)
        matchresult_2 = match_df(all_true_near_branch_sub_all, br_indexn_2)
        
        if(dim(matchresult_1)[1]==1 & dim(matchresult_2)[1]==0){
          truematchresult = matchresult_1
          gene_result=data.frame(i,j)
        }else if(dim(matchresult_1)[1]==0 & dim(matchresult_2)[1]==1){
          truematchresult = matchresult_2 
          gene_result=data.frame(i,j)
        }else{
          truematchresult=data.frame(matrix(, nrow=0, ncol=2))
          gene_result=data.frame(matrix(, nrow=0, ncol=2))
          
        }
      }
      colnames(gene_result)=c("gene_ID1","gene_ID2")
      
      if(dim(truematchresult)[1]>0){all_true_match_near_branch_sub=rbind(all_true_match_near_branch_sub, cbind(gene_result,truematchresult))}
      
      }
    
    }
  
  unique_all_match_near_branch_sub=unique(all_match_near_branch_sub)
  unique_all_true_match_near_branch_sub=unique(all_true_match_near_branch_sub)
  
  all_pair_result[indexn,]=c(dim(unique_all_match_near_branch_sub)[1],dim(unique_all_true_match_near_branch_sub)[1],numSites,alpha,missingindex,Genotype="Binary",lost_rate=0,lost_br=0,Method="MO")
  

}

all_pair_result_dat = rbind(all_pair_result_dat,all_pair_result)
colnames(all_pair_result)=c("inferred_pair","true_pair","sites","error","missing_index","Genotype","lost_rate","lost_br","Method")
all_pair_result_form_0_1_sub = sprintf('./FiguresS4_S5_Scenarios5_6/Scenario5_result_rmd/MO_Binary/summary/binary_adjacent_order_alpha_0%s_beta_0%s_%s.csv',alpha_str,alpha_str,numSites)
write.csv(all_pair_result,file=all_pair_result_form_0_1_sub)

}
