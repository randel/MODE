############## filterzerovar ##############3
filterzerovar <- function(mat){

  if(class(mat)[[1]] == "dgCMatrix"){
    mat <- mat[!rowVars(mat) == 0,]
  }else{
    mat <- mat[!apply(mat, 1, function(x) var(x) == 0),]
  }
  return(mat)
}

