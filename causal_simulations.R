simulate_from_regions(region_list,h2,n)

simulate_betas_w_ld = function(causal_number, ld_matrix, h2, n, idx=NULL){
  if(!is.null(idx)){
    betas = rnorm(n=length(idx),mean=0,sd=sqrt(h2/length(idx)))
    beta_vec = rep(0,ncol(ld_matrix))
    beta_vec[idx] = betas
    noises = rnorm(n=length(idx),mean=0,sd=sqrt((1-h2)/(causal_number *n)))
    noise_vec = rep(0,ncol(ld_matrix))
    noise_vec[idx] = noises
  }else{
    betas = rnorm(n=causal_number,mean=0,sd=sqrt(h2/causal_number))
    noises = rnorm(n=causal_number,mean=0,sd=sqrt((1-h2)/(causal_number *n)))
    beta_idx = sample(1:ncol(ld_matrix),size = causal_number,replace = F)
    beta_vec = rep(0,ncol(ld_matrix))
    beta_vec[beta_idx] = betas
    noise_vec = rep(0,ncol(ld_matrix))
    noise_vec[beta_idx] = noises
  }
  #beta_noises = rnorm(n=ncol(ld_matrix),0,1/n)
  #beta_c = rnorm(n = 100, mean = 0, sd=(sqrt(1/n)))
  beta_hat = ld_matrix %*% beta_vec
  beta_noises = ld_matrix %*% noise_vec
  beta_hat = beta_hat + beta_noises
  ses = ((1-h2)/(n)) * ld_matrix
  z = beta_hat/sqrt(diag(ses))
  return(list(beta_gwas=beta_hat, true_beta=beta_vec, ses=ses,z=z,idx=beta_idx))
}

recapitulate_results = function(sim_list, ld_matrix_inv, h2,n){
  beta_hat =  ld_matrix_inv %*% res.sim$beta_gwas
  ses = ((1-h2)/(n)) * ld_matrix_inv
  z = beta_hat/sqrt(diag(ses))
  return(list(beta_hat = beta_hat,z_hat=z, se_hat=sqrt(diag(ses))))
}

# simulate_unscaled_betas = function(effect, idx, ld_matrix){
#   n_snps = dim(ld_matrix)[1]
#   effects = rep(0, n_snps)
#   effects[idx] = effect
#   noises = rnorm(n_snps,0,1)
#   chol_ld_matrix = invert_ld_matrix(ld_matrix)
#   C = chol_ld_matrix
#   D_I <- chol2inv(chol_ld_matrix)
#   betas_ld <- ld_matrix %*% effects
#   noises_ld <- t(C) %*% noises
#   betas = betas_ld + noises_ld
#   return(list(true_betas=betas_ld, betas=betas,C=C,D_I=D_I))
# }
# 
# simulate_unscaled_betas_w_ld = function(effect, idx, ld_matrix, C, D_I, n){
#   n_snps = dim(ld_matrix)[1]
#   effect = effect/sqrt(n)
#   effects = rep(0, n_snps)
#   effects[idx] = effect
#   noises = rnorm(n_snps,0,1/n)
#   betas_ld <- ld_matrix %*% effects
#   noises_ld <- t(C) %*% noises
#   betas = betas_ld + noises_ld
#   return(list(true_betas=betas_ld, betas=betas))
#   
# }

# pick_snp = function(maf, snpdat, ldscore=NULL, sample_proportion=0.1){
#   #  if(ldscore!=NULL){ ; }
#   snpdat$FREQ1[snpdat$FREQ1>0.5] = 1- snpdat$FREQ1[snpdat$FREQ1>0.5]
#   snpdat$diffs = abs(maf - snpdat$FREQ1)
#   snpdat_2choose = order(snpdat$diffs)[1:max(1,dim(snpdat)[1]*sample_proportion)]
#   #  hist(snpdat$FREQ1[snpdat_2choose], xlim=c(0,0.5))
#   idx = sample(snpdat_2choose, size=1)
#   return(idx)
# }