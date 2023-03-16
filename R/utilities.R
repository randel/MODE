############## filterzerovar ##############
filterzerovar <- function(mat){

  if(class(mat)[[1]] == "dgCMatrix"){
    mat <- mat[!rowVars(mat) == 0,]
  }else{
    mat <- mat[!apply(mat, 1, function(x) var(x) == 0),]
  }
  return(mat)
}

############## BetaToMvalue ##############
BetaToMvalue = function(Beta){
  Beta [which(Beta ==0)] <- 0.01
  Beta [which(Beta <0.01)] <- 0.01
  Beta [which(Beta ==1)] <- 0.99
  Beta [which(Beta >0.99)] <- 0.99
  Mvalue <- log2(Beta/(1 - Beta))

  return(Mvalue)
}
