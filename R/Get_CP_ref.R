#' Get_CP_ref This fuction is to get CP estimate
#'
#' @param omic Character
#' @param Datalist List of bulk data
#' @param reflist List of reference data
#' @param refind Character
#' @param Marker.Method Marker selection methods
#' @param nmrk Number of markers
#' @param ref_list Specific ref_list, can be null if omic and reflist is provided.
#' @param count_bulk Specific bulk data, can be null if omic and Datalist is provided.
#'
#' @return
#' @export
#'
#' @examples
Get_CP_ref = function(omic = "RNA/Mic",Datalist,reflist,refind, Marker.Method = "none",nmrk = 50,ref_list = NULL,count_bulk = NULL){
  params = get_params(data_type = "singlecell-rna", data_name = refind, n_markers = nmrk, Marker.Method = Marker.Method,
                      TNormalization = "none", CNormalization = "none",
                      dmeths = c("CIBERSORT","EPIC","FARDEEP","DCQ"),Scale = "linear")
  if(is.null(ref_list) & is.null(count_bulk)){
    ref_list = reflist[[omic]]
    count_bulk = as.matrix(Datalist[[omic]]$count_bulk)
  }


  mrk = Reduce(intersect, list(rownames(ref_list[[refind]]$ref_matrix), rownames(count_bulk)))
  ref_list[[refind]]$ref_matrix = ref_list[[refind]]$ref_matrix[mrk,]
  count_bulk = count_bulk[mrk,]

  res = EnsDeconv(count_bulk = count_bulk, ref_list = ref_list,
                        params = params,parallel_comp = T,ncore = 6)

  return(res)

}
