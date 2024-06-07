#' Extract trajectory
#'
#' This function extracts all corrdinates for a trajectory created in monocle3. .
#' @param monocle_object A monocle3 cell_data_set object.
#' @param reduction The dimensionality reduction to use.
#' @param is_3D Whether trajectory has been drawn for a 3D dimension plot.
#' @param spread Effective scale of embedded points. Determines, along with min.dist, how clumped embedded points are.
#' @param bg_col The background color. Currently not supported.
#' @param file The file name to use for the GIF.
#' @param fps Frames per secund; controls the speed at which the plot will rotate.
#' @param trajectory_coords A data frame with coordinates for plotting a trajectory.
#' @param trajectory_col The color for the trajectory, if coordinates are provided.
#' @export
#' @examples
#' KG_get_trajectory(cds)
KG_get_trajectory = function(monocle_object,
                             reduction = "UMAP",
                             is_3D = FALSE){
  suppressPackageStartupMessages(c(require(dplyr),
                                   require(igraph),
                                   require(monocle3)))

  # Extract data point coordinates
  if(is_3D){
    data_coords = t(monocle_object@principal_graph_aux[[reduction]]$dp_mst) %>%
      as.data.frame() %>%
      dplyr::select(UMAP_1, UMAP_2, UMAP_3) %>%
      dplyr::mutate(sample_name = rownames(.),
                    sample_state = rownames(.))

    # Extract trajectory coordinates
    graph = monocle_object@principal_graph[[reduction]]

    edges = graph %>%
      igraph::as_data_frame() %>%
      dplyr::select(source = "from", target = "to") %>%
      dplyr::left_join(data_coords %>%
                         dplyr::select(source = "sample_name",
                                       source_dim1 = "UMAP_1",
                                       source_dim2 = "UMAP_2",
                                       source_dim3 = "UMAP_3"),
                       by = "source") %>%
      dplyr::left_join(data_coords %>%
                         dplyr::select(target = "sample_name",
                                       target_dim1 = "UMAP_1",
                                       target_dim2 = "UMAP_2",
                                       target_dim3 = "UMAP_3"),
                       by = "target")
  } else{
    data_coords = t(monocle_object@principal_graph_aux[[reduction]]$dp_mst) %>%
      as.data.frame() %>%
      dplyr::select(UMAP_1, UMAP_2) %>%
      dplyr::mutate(sample_name = rownames(.),
                    sample_state = rownames(.))

    # Extract trajectory coordinates
    graph = monocle_object@principal_graph[[reduction]]

    edges = graph %>%
      igraph::as_data_frame() %>%
      dplyr::select(source = "from", target = "to") %>%
      dplyr::left_join(data_coords %>%
                         dplyr::select(source = "sample_name",
                                       source_dim1 = "UMAP_1",
                                       source_dim2 = "UMAP_2"),
                       by = "source") %>%
      dplyr::left_join(data_coords %>%
                         dplyr::select(target = "sample_name",
                                       target_dim1 = "UMAP_1",
                                       target_dim2 = "UMAP_2"),
                       by = "target")
  }
  return(edges)
  gc()
}
