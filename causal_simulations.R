#simulate_from_regions(region_list,h2,n)

r_squared_one = function(idx, ld_matrix,betas){
  i = 1
  idxs = c()
  betas_out = c()
  for(j in idx){
    idx_which = which((ld_matrix[j,c(1:ncol(ld_matrix))]^2 >=.995))
    if(length(idx_which) <=1){
      idxs = c(idxs,j)
      betas_out = c(betas_out, betas[i])
      next
    }
    tmp_betas = betas[i]/length(idx_which)
    tmp_betas = ifelse(ld_matrix[j,idx_which] < 0,-tmp_betas,tmp_betas)
    betas_out = c(betas_out,tmp_betas)
    i = i + 1
    idxs = c(idxs,idx_which)
  }
  return(cbind(idxs,betas_out))
}

sim_partition = function(c1,c2,ld_matrix, category,h2_cat,h2_out,n,chol_mat){
  # H2_in and H2_out.
 # assign_enrichments=rbinom(n=causal_number,size=2,prob=0.5)
  # Out of category
  beta_vec = rep(0,ncol(ld_matrix))
  noise_vec = rep(0,ncol(ld_matrix))
  #no_enrichments = sum(assign_enrichments[assign_enrichments==1])
  betas = rnorm(n=c1,mean=0,sd=sqrt(h2_cat/c1))
  print(sqrt(h2_cat/c2*var(betas)))
 # betas = betas * sqrt(h2_cat/(c1*var(betas)))
  
  samp_idx = (1:ncol(ld_matrix))[!category]
  beta_idx = sample(samp_idx,size = c1,replace = F)
  #res = r_squared_one(beta_idx, ld_matrix, betas)
  #beta_idx = res[,1]
  #betas = res[,2]
  noises = rnorm(n=length(betas),mean=0,sd=1)
  beta_vec[beta_idx] = betas
  noise_vec[beta_idx] = noises
  # In category
  
  samp_idx = (1:ncol(ld_matrix))[category]
  #no_enrichments = sum(assign_enrichments[assign_enrichments==2])
  betas = rnorm(n=c2,mean=0,sd=sqrt(h2_out/c2))
  #print(sqrt((h2_cat/(c1*var(betas)))))
  #betas = betas * sqrt(h2_out/(c2*var(betas)))
  noises = rnorm(n=c2,mean=0,sd=(1))
  beta_idx = sample(samp_idx,size = c2,replace = F)
 # res = r_squared_one(beta_idx, ld_matrix, betas)
  #beta_idx = res[,1]
  #betas = res[,2]
  noises = rnorm(n=length(betas),mean=0,sd=(1-h2))
  beta_vec[beta_idx] = betas
  noise_vec[beta_idx] = noises  
  
  beta_hat_w_nois = beta_vec+ noise_vec
  beta_hat = ld_matrix %*% beta_vec
  beta_noises = ld_matrix %*% noise_vec
  #beta_hat = beta_hat + beta_noises
  beta_noise2 = sqrt(1/n) * chol_mat
  beta_noise2 = beta_noise2 %*% noise_vec
  print(sd(beta_noise2))
  print(sd(beta_noises))
  print(summary(beta_noise2))
  beta_hat = beta_hat + beta_noise2
  
  
  ses = (((1-h2_cat)/(n)) * ld_matrix)
  #print(beta_hat)
  
  chol_mat=NULL
  sim_results = list(beta_gwas=beta_hat,true_beta=beta_vec,beta_noise=beta_hat_w_nois,beta_hat = betas + noises)
  hess_estimate = hess_estimator(sim_results,ld_matrix,n)
 # message(paste("H^2 hess =",hess_estimate))
  #message(paste("H^2 ld ="), sum(beta_vec^2))
  return(list(beta_gwas=beta_hat,true_beta=beta_vec,ses=ses,noises=noises,sim_results=sim_results,h2_ld_score=sum(beta_vec^2),h2_hess=hess_estimate)) 
    
}



simulate_betas_w_ld = function(causal_number, ld_matrix, h2, n, idx=NULL,betas=NULL, only_these_idxs=NULL,chol_mat=NULL){
  if(!is.null(idx)){
    if(is.null(betas)){
      betas = rnorm(n=length(idx),mean=0,sd=sqrt(h2/length(idx)))
      }
    beta_vec = rep(0,ncol(ld_matrix))
    print(length(beta_vec))
    beta_vec[idx] =  betas * sqrt(h2_out/(c2*var(betas)))
    noises = rnorm(n=length(idx),mean=0,sd=sqrt((1-h2)/(n)))
    noise_vec = rep(0,ncol(ld_matrix))
    noise_vec[idx] = noises
    beta_idx = idx
  }else{
    betas = rnorm(n=causal_number,mean=0,sd=sqrt(h2/causal_number))
    #plot(betas)
    #print(var(betas))
    #print(causal_number*var(betas))
    #print(sqrt(h2/(causal_number*var(betas))))
    betas = betas * sqrt(h2/(causal_number*var(betas)))
   # plot(betas)
    noises = rnorm(n=causal_number,mean=0,sd=(1-h2))
    if(!is.null(only_these_idxs)){
      beta_idx = sample(1:ncol(ld_matrix),size = causal_number,replace = F)
    }else{
      beta_idx = sample(1:ncol(ld_matrix),size = causal_number,replace = F)
    }
    beta_vec = rep(0,ncol(ld_matrix))
    beta_vec[beta_idx] = betas
    noise_vec = rep(0,ncol(ld_matrix))
    noise_vec[beta_idx] = noises
  }
  #beta_noises = rnorm(n=ncol(ld_matrix),0,1/n)
  #beta_c = rnorm(n = 100, mean = 0, sd=(sqrt(1/n)))
  beta_hat_w_nois = beta_vec+ noise_vec
  beta_hat =  ld_matrix %*% beta_vec
  beta_noise2 = sqrt((1-h2)/n) * chol_mat
  beta_noise2 = t(chol_mat) %*% noise_vec
  #noise_vec = sqrt(1/n) %*% ld_matrix
 # noise_vec = noise_vec *sqrt((1-h2)/n)
#  beta_noises = ld_matrix %*% noise_vec
  beta_noise3 = sqrt((1-h2)/n) * ld_matrix %*% noise_vec
  
  
  #print(sd(beta_noises))
  #print(sd(beta_noise2))
  #print(sd(beta_noise3))
  #print(summary(beta_noises - beta_noise2))
  
  ses = (((1-h2)/(n)) * ld_matrix)
  z = beta_hat/sqrt(diag(ses))
  
 # beta_effective = ld_matrix %*% betas_v
 # noise_effective = ld_matrix %*% noises
  #beta_effective_noise = beta_effective + noise_effective
  a = list(beta_gwas=beta_hat, true_beta=beta_vec, ses=ses,z=z,idx=beta_idx,beta_noise=beta_hat_w_nois,beta_hat = betas + noises)
  hess_estimate = hess_estimator(a,ld_matrix,n)
  return(list(beta_gwas=beta_hat, true_beta=beta_vec, ses=ses,z=z,idx=beta_idx,beta_noise=beta_hat_w_nois,beta_hat = betas + noises,h2_hess=hess_estimate,h2_ld_score=sum(a$true_beta^2)))
}

hess_estimator = function(res.sim,ld_matrix,n,rank=502){
 #print(n* t(res.sim$true_beta) %*% ld_matrix %*% res.sim$true_beta)
  x = (t(res.sim$true_beta) %*% ld_matrix %*% res.sim$true_beta)
  x2 = (n*t(res.sim$true_beta) %*% ld_matrix %*% res.sim$true_beta - rank)/(n-rank)
  #print(x2)
  
  return(x)
}

recapitulate_results = function(res.sim, ld_matrix_inv, h2,n,k){
  beta_hat =  ld_matrix_inv %*% res.sim$beta_gwas
  ses = abs(((1-h2)/(n)) * ld_matrix_inv)
  #print(is.na(ses))
  z = beta_hat/sqrt(diag(ses))
  #with trace 
 # tr()
  print( t(res.sim$beta_gwas) %*% ld_matrix_inv %*% res.sim$beta_gwas)
# try_this = (1- t(res.sim$beta_gwas) %*% ld_matrix_inv %*% res.sim$beta_gwas)/n * k
 #print(try_this)
  return(list(beta_hat = beta_hat,z_hat=z, se_hat=(diag(ses))))
}

simulate_w_stats = function(h2,n,causals){
  
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