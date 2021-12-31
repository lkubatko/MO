#' Compute the probability that a mutation ocurrs on branches of a phylogeny
#'
#' This function computes the probability that a mutation occurs on each branch of a tree. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param alpha false positive error rate/probability
#' @param beta false negative error rate/probability
#' @param lambda01 unit_theta=lambda01 rate from 0 to 1
#' @param lambda02 unit_gamma=lambda02 rate from 0 to 2
#' @param lambda12 unit_mu=lambda12 rate from 1 to 2
#' @param initial_obs observed mutation data matrix of all cells in the tree, each row is one mutation and each column is one cell
#' @param t phylogenetic tree
#' @param binary a TRUE/FALSE of whether the data is binary or ternary
#' @param multinomialPrior a TRUE/FALSE of whether the prior of mutations are multinomial
#' @return A probability matrix of the mutations on each branch
#' @export
MO <- function(alpha,beta,lambda01,lambda02,lambda12,initial_obs,t,binary=c(TRUE,FALSE),multinomialPrior=c(FALSE,TRUE),seed=1000){
   
  alpha = alpha
  beta = beta
  unit_theta = lambda01
  unit_gamma = lambda02
  unit_mu = lambda12
  initial_obs = initial_obs
  t = t
  binary = binary
  multinomialPrior = multinomialPrior
  brNum = 2*(length(t$tip.label)-1)
  set.seed(seed)
  if(multinomialPrior==FALSE){
    if(binary==TRUE){
      
      P1_binary_prob_matrix_0_1 = c()
      selected_br = c()
      for (i in 1:dim(initial_obs)[1]){
        
        print(i)
        #rd_unit_lambda <- rgamma(10, shape = 2*unit_gamma, scale = (unit_gamma/2)) 
        rd_unit_theta =  rgamma(n = 10, shape = 2, scale = 0.5*unit_theta)
        rd_unit_gamma = rgamma(n = 10, shape = 2, scale = 0.5*unit_gamma)
        
        
        generate_prob_br_0_1_dat=data.frame(matrix(NA, nrow = brNum, ncol = 10))
        
        for (indexj in 1:10){
          
          generate_prob_br <- MO_binary(alpha,beta,rd_unit_theta[indexj] + rd_unit_gamma[indexj],initial_obs[i,],t)
          
         
          generate_prob_br_0_1_single <- generate_prob_br
          
          
          generate_prob_br_0_1_dat[,indexj] = generate_prob_br_0_1_single
          
        }
        
        generate_prob_br_0_1=rowMeans(generate_prob_br_0_1_dat, na.rm = FALSE, dims = 1)
        
        P1_binary_prob_matrix_0_1 = rbind(P1_binary_prob_matrix_0_1,generate_prob_br_0_1)
        selected_br[i] = which.max(generate_prob_br_0_1)
      }
      
      prob = cbind.data.frame(initial_obs[,1],selected_br,P1_binary_prob_matrix_0_1)
      colnames(prob) = c("gene","MutationBr",paste0("br",1:brNum))
      rownames(prob) = NULL
      }else if(binary==FALSE){
        
        #initialize vector to store estimated probability.Initialize vector to store MAP estimates
        ternary_prob_matrix_all_0_1=c()
        selected_br=c()
        
        for (i in 1:dim(initial_obs_0_1_recode)[1]){
          
          #compute probability that mutation i on each of the branches on the tree
          rd_unit_theta =  rgamma(n = 10, shape = 100, scale = 0.01*unit_theta)
          rd_unit_gamma = rgamma(n = 10, shape = 100, scale = 0.01*unit_gamma)
          rd_unit_mu = rgamma(n = 10, shape = 100, scale = 0.01*unit_mu)
          
          generate_prob_br_all_dat=data.frame(matrix(NA, nrow = brNum, ncol = 10))
          
          for (j in 1:10){
            
            generate_prob_br <- MO_ternary(alpha,beta,rd_unit_theta[j],rd_unit_gamma[j],rd_unit_mu[j],initial_obs_0_1_recode[i,],sampletr)
            generate_prob_br_all_single <- c(generate_prob_br,rep(0,brNum-length(generate_prob_br)))
            generate_prob_br_all_dat[,j] = generate_prob_br_all_single}
          
          probs=rowMeans(generate_prob_br_all_dat, na.rm = FALSE, dims = 1)
          #find the MAP
          selected_br[i]=which.max(probs)
          #save the ith result to the matrix
          ternary_prob_matrix_all_0_1 = rbind(ternary_prob_matrix_all_0_1,probs)
          
        }
        
        # combine the true mutation branch and the inferred mutation branch for future comparison
        #ternary_prob_matrix_all_0_1_rownames=cbind(selected_br,ternary_prob_matrix_all_0_1)
       
        prob = cbind.data.frame(initial_obs[,1],selected_br,ternary_prob_matrix_all_0_1)
        colnames(prob) = c("gene","MutationBr",paste0("br",1:brNum))
        rownames(prob) = NULL
      }
      
      
          }else{
            
      P1_binary_prob_matrix_0_1=c()
      selected_br = c()
      for (i in 1:dim(initial_obs)[1]){
        
        print(i)
        
        generate_prob_br_0_1_dat=data.frame(matrix(NA, nrow = brNum, ncol = 10))
        
       
        generate_prob_br_0_1 <- MO_binary_multinomial(alpha,beta,initial_obs[i,],t)

        P1_binary_prob_matrix_0_1 = rbind(P1_binary_prob_matrix_0_1,generate_prob_br_0_1)
        selected_br[i] = which.max(generate_prob_br_0_1)
      }
      
      prob = cbind.data.frame(initial_obs[,1],selected_br,P1_binary_prob_matrix_0_1)
      colnames(prob) = c("gene","MutationBr",paste0("br",1:brNum))
      rownames(prob) = NULL
    }
  
  
  return(prob)
}
