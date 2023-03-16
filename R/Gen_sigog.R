#' Gen_sigog: Generate standardizing factors
#'
#' @param omic omic indicator
#' @param Datalist
#' @param reflist
#' @param cp
#' @param refind
#' @param genewise
#' @param bulk
#' @param sig
#' @param ens_ad
#'
#' @return
#' @export
#'
#' @examples
Gen_sigog <-   function(omic,Datalist =FHS_ct4_list ,reflist = ref_list_log,cp = NULL,refind,genewise = T, bulk = NULL,sig = NULL,ens_ad){
  if(is.null(bulk) & is.null(sig)){
    bulk = Datalist[[omic]]$count_bulk
    sig = reflist[[omic]][[refind]]$ref_matrix
  }


  mrk = Reduce(intersect, list(rownames(sig), rownames(bulk)))
  sig <- sig[mrk,]
  bulk <- bulk[mrk,]

  if(is.null(cp)){
    cp = ens_ad$ensemble_p
  }


  sbjind = intersect(rownames(cp),colnames(bulk))
  bulk = bulk[,sbjind]
  cp = cp[sbjind,]

  ct_name = colnames(cp)
  mrk = intersect(rownames(sig),rownames(bulk))

  sig = sig[mrk,ct_name]
  bulk = bulk[mrk,]
  prod <-  sig%*%t(cp)
  tmp <-  bulk - prod
  if(genewise){
    ss <- apply(tmp, 1, function(x) sqrt(mean(x^2)))
  }else{
    ss <- sqrt(mean(tmp^2))
  }

  new_bulk <-  bulk/ss
  new_sig = sig/ss
  return(list(bulk = new_bulk ,sig = new_sig, ss= ss))
}
