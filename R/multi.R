#' Title
#'
#' @param Dataall
#' @param outpath
#' @param Marker.Method
#' @param nmrk
#'
#' @return
#' @export
#'
#' @examples
multi = function(Dataall,outpath,Marker.Method = "none",nmrk = 50){

  # subj <-  overlap[[ind]]
  # Data_sub <- Dataall[strsplit(ind, "_")[[1]]]
  sig <-  do.call(rbind, lapply(Dataall, function(res){
    sig = res$sig
    return(sig)
  }))
  bulk <-  do.call(rbind, lapply(Dataall, function(res){
    bulk = res$bulk
    return(bulk)
  }))

  mrk = Reduce(intersect, list(rownames(sig), rownames(bulk)))
  sig <- sig[mrk,]
  bulk <- bulk[mrk,]
  #EnsDeconv with marker.method == "none", no markers number is needed
  mettt = data.frame(deconv_clust = colnames(sig), SamplesName=colnames(sig))
  ref_list = list()
  ref_list$"ref" = list()
  ref_list$"ref"$ref_matrix = sig
  ref_list$"ref"$meta_ref = mettt

  params = get_params(data_type = "singlecell-rna", data_name = "ref", n_markers = nmrk, Marker.Method = Marker.Method,
                      TNormalization = "none", CNormalization = "none",
                      dmeths = c("CIBERSORT","EPIC","FARDEEP","DCQ"),Scale = "linear")

  allgene_res = gen_all_res_list(count_bulk = bulk, ref_list = ref_list,
                                 params = params,parallel_comp = T,ncore = 6,outpath =paste0(outpath,"/ind"))

  ind = sapply(allgene_res, function(x) {
    length(x[["a"]][["p_hat"]][[1]])
  })
  allgene_res = allgene_res[which(ind == 1)]

  return(allgene_res)
}
