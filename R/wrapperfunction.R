#' wrapperfunction is a wrapper for MODE
#'
#' @param Datalist
#' @param reflist
#' @param outpath
#' @param cp
#' @param CP_list
#' @param omicparam
#' @param Marker.Method
#' @param nmrk
#' @param genewise
#'
#' @return
#' @export
#'
#' @examples
wrapperfunction = function(Datalist =FHS_ct4_list ,
                               reflist = ref_list_log,outpath = "D:/Manqi/Multi-omics",cp = NULL,CP_list,omicparam = NULL
                               ,Marker.Method = "none",nmrk = 50,genewise = T){

  if(is.null(omicparam)){
    omicparam <- data.frame(omics = c("DNAm","RNA"),
                            ref = c("K1000C2","Darmanis"))
    omicparam$names = paste0(omicparam$omics,omicparam$ref)
  }
  if(is.null(cp)){
    ens_ad <- CP_list$EnsDeconv

    new_data <-  lapply(names(Datalist), Gen_sigog,ens_ad = ens_ad,Datalist =Datalist ,reflist = reflist,genewise=genewise)
  }else{

    new_data <-  lapply(1:nrow(omicparam),function(i){
      res = Gen_sigog(omic = omicparam$omics[i],cp = cp,Datalist =Datalist ,reflist = reflist,
                      refind = omicparam$ref[i],genewise = genewise)
    } )
  }

  names(new_data) = omicparam$names

  nomic =nrow(omicparam)/length(Datalist)
  Dataall = list(new_data[[i]],new_data[[j]])

  res_con =multi(Dataall = Dataall,outpath = outpath,Marker.Method = Marker.Method,nmrk = nmrk)

  return(res_con)
}
