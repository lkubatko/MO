
library(ape)
library(geiger)
library(ips)
library(phangorn)
library(TreeSearch)
#library(expm)
library(doParallel)
#library(truncnorm)
library(data.tree)
library(gtools)
library(dplyr)
library(phytools)


   
    
    for(iter in 1:1000){
     
    
    true_tree_form = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario1/RandomTree/RandomTree_%s.tre",iter)
    
    tree_form = sprintf('/fs/project/kubatko.2-temp/gao.957/workspace/github_MO/Figures3_4_Scenarios1_2/Scenario1/RandomTree/Scaled_RandomTree_%s.tre',iter)
    
    true_tree = read.tree(true_tree_form)
    
    converted_true_tree = true_tree
    
    converted_true_tree$edge.length = (true_tree$edge.length) * (10^(-7))
    
    write.tree(converted_true_tree, file = tree_form)
    
    
    
    }

  