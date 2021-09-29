get_regularization <- function(){
  
  
  if(is.null(list_structured[[2]]))
  {
    if(!is.null(lambda_lasso) & is.null(lambda_ridge)){
      
      reg = function(x) tf$keras$regularizers$l1(l=lambda_lasso)(
        model$trainable_weights[[mono_layer_ind]]
      )
      
    }else if(!is.null(lambda_ridge) & is.null(lambda_lasso)){
      
      reg = function(x) tf$keras$regularizers$l2(l=lambda_ridge)(
        model$trainable_weights[[mono_layer_ind]]
      )
      
      
    }else if(!is.null(lambda_ridge) & !is.null(lambda_lasso)){
      
      reg = function(x) tf$keras$regularizers$l1_l2(
        l1=lambda_lasso,
        l2=lambda_ridge)(model$trainable_weights[[mono_layer_ind]])
      
    }else{
      
      reg = NULL # no penalty
      
    }
    
    reg2 = NULL
    
  }else{
    
    if(penalize_bsp)
      bspP <- secondOrderPenBSP(order_bsp, order_diff = order_bsp_penalty)
    bigP <- list_structured[[2]]
    if(!is.null(deep_part_ia))
    {
      
      bigP <- bdiag(list(bigP, diag(rep(1, ncol(deep_part_ia)[[1]]))))
      
    }
    if(length(bigP@x)==0 & penalize_bsp==0){
      reg = NULL
    }else if(penalize_bsp==0){
      
      reg = function(x) k_mean(k_batch_dot(model$trainable_weights[[mono_layer_ind]], k_dot(
        # tf$constant(
        sparse_mat_to_tensor(as(kronecker(diag(rep(1, ncol(input_theta_y)[[1]])),bigP),
                                "CsparseMatrix")),
        # dtype = "float32"),
        model$trainable_weights[[mono_layer_ind]]),
        axes=2) # 1-based
      )
      
    }else if(length(bigP@x)==0)
    {
      
      reg = function(x) k_mean(k_batch_dot(model$trainable_weights[[mono_layer_ind]],
                                           k_dot(
                                             # tf$constant(
                                             sparse_mat_to_tensor(as(kronecker(penalize_bsp*bspP,
                                                                               diag(rep(1, ncol(interact_pred)[[1]]))),
                                                                     "CsparseMatrix")),
                                             # dtype = "float32"),
                                             model$trainable_weights[[mono_layer_ind]]),
                                           axes=2) # 1-based
      )
      
    }else{
      
      reg = function(x) k_mean(k_batch_dot(model$trainable_weights[[mono_layer_ind]],
                                           k_dot(
                                             # tf$constant(
                                             sparse_mat_to_tensor(
                                               as(kronecker(penalize_bsp*bspP, diag(rep(1, ncol(interact_pred)[[1]])))
                                                  + kronecker(diag(rep(1, ncol(input_theta_y)[[1]])),bigP), "CsparseMatrix")),
                                             # dtype = "float32"),
                                             model$trainable_weights[[mono_layer_ind]]),
                                           axes=2) # 1-based
      )
      
    }

  }
  
  return(reg)
  
}