invert_ld_matrix_corcorp = function(ld_matrix){
  #Use default
  ld_shrink = cor.shrink(ld_matrix)
  ld_inv = chol2inv(chol(ld_shrink))
  return(list(ld_matrix_inv=ld_inv,lambda=attr(ld_shrink, "lambda")))
}

invert_ld_matrix_pca = function(ld_matrix, k){
  # invert LD matrix with principal components.
  ld_matrix_inv = svd(ld_matrix)
  test_d = 1/s_vd$d
  test_d[(k):length(test_d)] = 0
  ld_invert = ((s_vd$v) %*% diag(test_d) %*% t(s_vd$u)) 
  return(list(ld_matrix_inv=ld_invert, k=k))
}

# TODO Look at using fast.SVD from the corpcor R package.
#
#
invert_ld_matrix = function(ld_matrix,max_iter =100,lambda = NULL,method=c("pca"),k=NULL){
  # Check if the method is one of the valid ones.
  if(!(method %in% c("pca","ridge","pca_ridge","corpcop","random"))){
    stop("Method must be one of pca, ridge, pca_ridge, corpcor, random")
  }
  if(method =="pca"){
    ld_invert = invert_ld_matrix_pca(ld_matrix, k)  
  }else if(method == "ridge"){
    ld_invert = invert_ld_matrix_ridge(ld_matrix, lambda=lambda)
  }else if(method == "corcor"){
    ld_invert = invert_ld_matrix_corpcor(ld_matrix)
  }else if(method == "random"){
    # Randamised PCA
    stop("Not implemented ")
  }else if(method == "pca_ridge"){
    # New method for regularising matrices.
    stop("Not implemented")
  }
  return(ld_invert)
}  

invert_ld_matrix_ridge = function(ld_matrix,max_iter =100,lambda = NULL ){
  invert_successful = F
  num_iterations = 1
  reg_factor = 0.001
  if(!is.null(lambda)){
    diag(ld_matrix) = diag(ld_matrix) + lambda
    chol_ld_matrix <- chol(ld_matrix)
    return(chol2inv(chol_ld_matrix))
  }
  while(num_iterations < max_iter){
    if(invert_successful){
      break
    }
    if(num_iterations > 1){
      ld_matrix = res.chol$ld_matrix
      reg_factor = res.chol$reg_factor
    }
    res.chol = tryCatch({
      num_iterations = num_iterations + 1
      chol_ld_matrix <- chol(ld_matrix)
      invert_successful = T
    },error=function(e){
      message(num_iterations)
      diag(ld_matrix) = diag(ld_matrix)+reg_factor
      reg_factor = reg_factor * 1.25
      return(list(ld_matrix=ld_matrix,reg_factor=reg_factor))
    })
  }
  if(!invert_successful){
    stop("Could not invert LD matrix.")
  }
  return(list(ld_matrix=chol2inv(chol_ld_matrix),reg_factor=reg_factor))
}