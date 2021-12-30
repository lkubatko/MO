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

parameter_setting = expand.grid(alpha=c(0.05,0.1,0.2,0.4),
                                sites =c(20,40,80))


##########
for(paraInd in 1:dim(parameter_setting)[1]){
  
  
  
  alpha =   parameter_setting[paraInd,1]
  beta = parameter_setting[paraInd,1]
  numSites= parameter_setting[paraInd,2]
  
  sequencing_error_model=matrix(c(1-alpha,alpha,
                                  beta,1-beta),nrow=2,byrow = TRUE)
  print(sequencing_error_model)
  unit_theta = 1
  unit_gamma = 0#10^(-14)
  unit_mu = 0#10 ^(-2)
  number_br = 20
  number_cell = 11
  
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
  
  
  for (indexn in 1:100){
    print(indexn)
    obs_dat_form = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/br_collapsed_snv_hap%s.csv",alpha_str,beta_str,numSites,indexn)
    obs_dat = read.csv(obs_dat_form)
    
    obs_dat_form_10 = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/br_collapsed_snv_hap%s_10.csv",alpha_str,beta_str,numSites,indexn)
    obs_dat_10 = read.csv(obs_dat_form_10)
    
    obs_dat_form_20 = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/br_collapsed_snv_hap%s_20.csv",alpha_str,beta_str,numSites,indexn)
    obs_dat_20 = read.csv(obs_dat_form_20)
    
    obs_dat_recode = data.frame(matrix(ncol = dim(obs_dat)[2], nrow = dim(obs_dat)[1]))
    obs_dat_recode[,1:3] = obs_dat[,1:3]
    obs_dat_recode_10 = data.frame(matrix(ncol = dim(obs_dat_10)[2], nrow = dim(obs_dat_10)[1]))
    obs_dat_recode_10[,1:3] = obs_dat_10[,1:3]
    obs_dat_recode_20 = data.frame(matrix(ncol = dim(obs_dat_20)[2], nrow = dim(obs_dat_20)[1]))
    obs_dat_recode_20[,1:3] = obs_dat_20[,1:3]
    
    for ( nr in 1:dim(obs_dat)[1]){
      for ( nc in 4:dim(obs_dat)[2]){
        if (obs_dat[nr,nc] == "-"){obs_dat_recode[nr,nc] ="-"}
        if (obs_dat[nr,nc] == "0"){obs_dat_recode[nr,nc] ="0"}
        if (obs_dat[nr,nc] == "1"){obs_dat_recode[nr,nc] ="1"}
        if (obs_dat[nr,nc] == "2"){obs_dat_recode[nr,nc] ="1"}
      }
    }
    
    for ( nr in 1:dim(obs_dat_10)[1]){
      for ( nc in 4:dim(obs_dat_10)[2]){
        if (obs_dat_10[nr,nc] == "-"){obs_dat_recode_10[nr,nc] ="-"}
        if (obs_dat_10[nr,nc] == "0"){obs_dat_recode_10[nr,nc] ="0"}
        if (obs_dat_10[nr,nc] == "1"){obs_dat_recode_10[nr,nc] ="1"}
        if (obs_dat_10[nr,nc] == "2"){obs_dat_recode_10[nr,nc] ="1"}
      }
    }
    
    for ( nr in 1:dim(obs_dat_20)[1]){
      for ( nc in 4:dim(obs_dat_20)[2]){
        if (obs_dat_20[nr,nc] == "-"){obs_dat_recode_20[nr,nc] ="-"}
        if (obs_dat_20[nr,nc] == "0"){obs_dat_recode_20[nr,nc] ="0"}
        if (obs_dat_20[nr,nc] == "1"){obs_dat_recode_20[nr,nc] ="1"}
        if (obs_dat_20[nr,nc] == "2"){obs_dat_recode_20[nr,nc] ="1"}
      }
    }
    
    colnames(obs_dat_recode)=colnames(obs_dat)
    colnames(obs_dat_recode_10)=colnames(obs_dat_10)
    colnames(obs_dat_recode_20)=colnames(obs_dat_20)
    
    gene_name= c(1:dim(obs_dat)[1])
    obs_dat_gene_name=cbind(gene_name,obs_dat_recode)
    obs_dat_gene_name_10=cbind(gene_name,obs_dat_recode_10)
    obs_dat_gene_name_20=cbind(gene_name,obs_dat_recode_20)
    
    
    
    obs_dat_form_gene_name = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/indexed_br_collapsed_snv_hap%s.csv",alpha_str,alpha_str,numSites,indexn)
    
    obs_dat_form_gene_name_10 = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/indexed_br_collapsed_snv_hap%s_10.csv",alpha_str,alpha_str,numSites,indexn)
    
    obs_dat_form_gene_name_20 = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/indexed_br_collapsed_snv_hap%s_20.csv",alpha_str,alpha_str,numSites,indexn)
    
    
    write.csv(obs_dat_gene_name, file = obs_dat_form_gene_name,row.names=FALSE)
    write.csv(obs_dat_gene_name_10, file = obs_dat_form_gene_name_10,row.names=FALSE)
    write.csv(obs_dat_gene_name_20, file = obs_dat_form_gene_name_20,row.names=FALSE)
    
    
  }
  
}

#########index
for(paraInd in 1:dim(parameter_setting)[1]){
  
  print(paraInd)
  
  alpha =   parameter_setting[paraInd,1]
  beta = parameter_setting[paraInd,1]
  numSites= parameter_setting[paraInd,2]
  
  sequencing_error_model=matrix(c(1-alpha,alpha,
                                  beta,1-beta),nrow=2,byrow = TRUE)
  print(sequencing_error_model)
  unit_theta = 1
  unit_gamma = 0#10^(-14)
  unit_mu = 0#10 ^(-2)
  number_br = 20
  number_cell = 11
  
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
  
  
  for (indexn in 1:100){
    print(indexn)
    
    
    obs_dat_form_10_k1 = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/br_collapsed_snv_hap%s_lostrate01_k1.csv",alpha_str,beta_str,numSites,indexn)
    obs_dat_10_k1 = read.csv(obs_dat_form_10_k1)
    
    obs_dat_form_20_k1= sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/br_collapsed_snv_hap%s_lostrate02_k1.csv",alpha_str,beta_str,numSites,indexn)
    obs_dat_20_k1 = read.csv(obs_dat_form_20_k1)
    
    
    obs_dat_form_10_k2 = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/br_collapsed_snv_hap%s_lostrate01_k2.csv",alpha_str,beta_str,numSites,indexn)
    obs_dat_10_k2 = read.csv(obs_dat_form_10_k2)
    
    obs_dat_form_20_k2= sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/br_collapsed_snv_hap%s_lostrate02_k2.csv",alpha_str,beta_str,numSites,indexn)
    obs_dat_20_k2 = read.csv(obs_dat_form_20_k2)
    
    
    
    
    obs_dat_recode_10_k1 = data.frame(matrix(ncol = dim(obs_dat_10_k1)[2], nrow = dim(obs_dat_10_k1)[1]))
    obs_dat_recode_10_k1[,1:3] = obs_dat_10_k1[,1:3]
    
    obs_dat_recode_20_k1 = data.frame(matrix(ncol = dim(obs_dat_20_k1)[2], nrow = dim(obs_dat_20_k1)[1]))
    obs_dat_recode_20_k1[,1:3] = obs_dat_20_k1[,1:3]
    
    obs_dat_recode_10_k2 = data.frame(matrix(ncol = dim(obs_dat_10_k2)[2], nrow = dim(obs_dat_10_k2)[1]))
    obs_dat_recode_10_k2[,1:3] = obs_dat_10_k2[,1:3]
    
    obs_dat_recode_20_k2 = data.frame(matrix(ncol = dim(obs_dat_20_k2)[2], nrow = dim(obs_dat_20_k2)[1]))
    obs_dat_recode_20_k2[,1:3] = obs_dat_20_k2[,1:3]
    
    for ( nr in 1:dim(obs_dat_10_k1)[1]){
      for ( nc in 4:dim(obs_dat_10_k1)[2]){
        if (obs_dat_10_k1[nr,nc] == "-"){obs_dat_recode_10_k1[nr,nc] ="-"}
        if (obs_dat_10_k1[nr,nc] == "0"){obs_dat_recode_10_k1[nr,nc] ="0"}
        if (obs_dat_10_k1[nr,nc] == "1"){obs_dat_recode_10_k1[nr,nc] ="1"}
        if (obs_dat_10_k1[nr,nc] == "2"){obs_dat_recode_10_k1[nr,nc] ="1"}
      }
    }
    
    for ( nr in 1:dim(obs_dat_20_k1)[1]){
      for ( nc in 4:dim(obs_dat_20_k1)[2]){
        if (obs_dat_20_k1[nr,nc] == "-"){obs_dat_recode_20_k1[nr,nc] ="-"}
        if (obs_dat_20_k1[nr,nc] == "0"){obs_dat_recode_20_k1[nr,nc] ="0"}
        if (obs_dat_20_k1[nr,nc] == "1"){obs_dat_recode_20_k1[nr,nc] ="1"}
        if (obs_dat_20_k1[nr,nc] == "2"){obs_dat_recode_20_k1[nr,nc] ="1"}
      }
    }
    
    for ( nr in 1:dim(obs_dat_10_k2)[1]){
      for ( nc in 4:dim(obs_dat_10_k2)[2]){
        if (obs_dat_10_k2[nr,nc] == "-"){obs_dat_recode_10_k2[nr,nc] ="-"}
        if (obs_dat_10_k2[nr,nc] == "0"){obs_dat_recode_10_k2[nr,nc] ="0"}
        if (obs_dat_10_k2[nr,nc] == "1"){obs_dat_recode_10_k2[nr,nc] ="1"}
        if (obs_dat_10_k2[nr,nc] == "2"){obs_dat_recode_10_k2[nr,nc] ="1"}
      }
    }
    
    for ( nr in 1:dim(obs_dat_20_k2)[1]){
      for ( nc in 4:dim(obs_dat_20_k2)[2]){
        if (obs_dat_20_k2[nr,nc] == "-"){obs_dat_recode_20_k2[nr,nc] ="-"}
        if (obs_dat_20_k2[nr,nc] == "0"){obs_dat_recode_20_k2[nr,nc] ="0"}
        if (obs_dat_20_k2[nr,nc] == "1"){obs_dat_recode_20_k2[nr,nc] ="1"}
        if (obs_dat_20_k2[nr,nc] == "2"){obs_dat_recode_20_k2[nr,nc] ="1"}
      }
    }
    
    colnames(obs_dat_recode_10_k1)=colnames(obs_dat_10_k1)
    colnames(obs_dat_recode_20_k1)=colnames(obs_dat_20_k1)
    colnames(obs_dat_recode_10_k2)=colnames(obs_dat_10_k2)
    colnames(obs_dat_recode_20_k2)=colnames(obs_dat_20_k2)
    
    gene_name= c(1:dim(obs_dat_10_k1)[1])
    
    obs_dat_gene_name_10_k1=cbind(gene_name,obs_dat_recode_10_k1)
    obs_dat_gene_name_20_k1=cbind(gene_name,obs_dat_recode_20_k1)
    obs_dat_gene_name_10_k2=cbind(gene_name,obs_dat_recode_10_k2)
    obs_dat_gene_name_20_k2=cbind(gene_name,obs_dat_recode_20_k2)
    
    
    
    obs_dat_form_gene_name_10_k1 = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/indexed_br_collapsed_snv_hap%s_lostrate01_k1.csv",alpha_str,alpha_str,numSites,indexn)
    
    obs_dat_form_gene_name_20_k1 = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/indexed_br_collapsed_snv_hap%s_lostrate02_k1.csv",alpha_str,alpha_str,numSites,indexn)
    
    obs_dat_form_gene_name_10_k2 = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/indexed_br_collapsed_snv_hap%s_lostrate01_k2.csv",alpha_str,alpha_str,numSites,indexn)
    
    obs_dat_form_gene_name_20_k2 = sprintf("/fs/project/kubatko.2-temp/gao.957/workspace/cellcoal-master/cellcoal_10tips/results_alpha_0%s_beta_0%s_%s/snv_haplotypes_dir/indexed_br_collapsed_snv_hap%s_lostrate02_k2.csv",alpha_str,alpha_str,numSites,indexn)
    
    
    
    write.csv(obs_dat_gene_name_10_k1, file = obs_dat_form_gene_name_10_k1,row.names=FALSE)
    write.csv(obs_dat_gene_name_20_k1, file = obs_dat_form_gene_name_20_k1,row.names=FALSE)
    write.csv(obs_dat_gene_name_10_k2, file = obs_dat_form_gene_name_10_k2,row.names=FALSE)
    write.csv(obs_dat_gene_name_20_k2, file = obs_dat_form_gene_name_20_k2,row.names=FALSE)
    
  }
  
}
