#' Create a 3D plot as a GIF
#'
#' This function generates 3D dimension plots as GIFs.
#' @param seurat_object A Seurat object
#' @param dims A vector of dimensions to use.
#' @param pt.size Size of data points.
#' @param cols A vector of colors to use for the shading of data points.
#' @param group.by Group identifier to use for shading.
#' @param run_umap Whether or not to run UMAP from scratch.
#' @param reduction The dimensionality reduction to use.
#' @param min.dist Controls how tightly embeddings are allowed to compress points together. Larger values result in more even distribution. Smaller values optimize more accurately the local structure.
#' @param spread Effective scale of embedded points. Determines, along with min.dist, how clumped embedded points are.
#' @param bg_col The background color. Currently not supported.
#' @param file The file name to use for the GIF.
#' @param fps Frames per secund; controls the speed at which the plot will rotate.
#' @param trajectory_coords A data frame with coordinates for plotting a trajectory.
#' @param trajectory_col The color for the trajectory, if coordinates are provided.
#' @export
#' @examples
#' KG_3D_to_gif(B_cell_dataset, file = "3D_plot.gif")
#'
KG_3D_to_gif = function(seurat_object,
                        dims = 1:30,
                        pt.size = 2,
                        cols = NULL,
                        group.by = NULL,
                        run_umap = TRUE,
                        reduction = "pca",
                        min.dist = 0.3,
                        spread = 1,
                        bg_col = "white", # currently not supported
                        file = "plot.gif",
                        fps = 5,
                        trajectory_coords = NULL,
                        trajectory_col = NULL){ # currently not supported

  suppressPackageStartupMessages({
    require(SeuratObject)
    require(Seurat)
    require(scales)
    require(rgl)
    require(magick)
    require(stats)})

  # Set the scene
  seurat_object[['ident']] = factor(SeuratObject::Idents(seurat_object))
  group.by = group.by %||% 'ident'
  cols = cols %||% scales::hue_pal()(length(unique(seurat_object[['ident', drop = T]])))
  trajectory_col = trajectory_col %||% "black"

  # Extract cell embeddings
  if(run_umap){
    temp = Seurat::RunUMAP(seurat_object,
                           dims = dims,
                           reduction = reduction,
                           min.dist = min.dist,
                           spread = spread,
                           n.components = 3)
    df = data.frame(UMAP_1 = temp@reductions$umap@cell.embeddings[,1],
                    UMAP_2 = temp@reductions$umap@cell.embeddings[,2],
                    UMAP_3 = temp@reductions$umap@cell.embeddings[,3],
                    Idents = temp[[group.by]])
    names(df)[4] = "Idents"
    rm(temp)
    gc()
  } else{
    df = data.frame(UMAP_1 = seurat_object@reductions$umap@cell.embeddings[,1],
                    UMAP_2 = seurat_object@reductions$umap@cell.embeddings[,2],
                    UMAP_3 = seurat_object@reductions$umap@cell.embeddings[,3],
                    Idents = seurat_object[[group.by]])
    names(df)[4] = "Idents"
  }
  # Map colors to IDs
  coloramp = stats::setNames(cols, unique(df$Idents))

  # Set color scheme for the plot based on background color
  material_color = if(bg_col == "black") "white" else "black"

  # Plot 3D data
  rgl::open3d(windowRect = c(1920, 143, 2550, 681)) +
    #material3d(color = "black") +
    rgl::bg3d(bg_col) +
    rgl::plot3d(x = df$UMAP_1,
                y = df$UMAP_2,
                z = df$UMAP_3,
                xlab = names(df)[1],
                ylab = names(df)[2],
                zlab = names(df)[3],
                col = coloramp[df$Idents],
                size = pt.size)

  # Plot trajectory if coordinates have been provided
  if(!is.null(trajectory_coords)){
    for(i in 1:nrow(trajectory_coords)){
      rgl::lines3d(x = as.vector(t(trajectory_coords[i, c("source_dim1", "target_dim1")])),
                   y = as.vector(t(trajectory_coords[i, c("source_dim2", "target_dim2")])),
                   z = as.vector(t(trajectory_coords[i, c("source_dim3", "target_dim3")])),
                   color = trajectory_col,
                   lwd = 2)


      #rgl::lines3d(x = as.vector(t(trajectory_coords$source_dim1[i], trajectory_coords$target_dim1[i])),
      #             y = as.vector(t(trajectory_coords$source_dim2[i], trajectory_coords$target_dim2[i])),
      #             z = as.vector(t(trajectory_coords$source_dim3[i], trajectory_coords$target_dim3[i])),
      #             color = trajectory_col,
      #             lwd = 2)
    }
  }

  # Keep track of snapshot file names
  frame_files = vector()

  # Take snapshots
  for(theta in seq(0, 360, by = 5)){
    rgl::view3d(theta = theta, phi = 30, zoom = 1)
    Sys.sleep(0.1)
    # Capture each frame
    frame_name = sprintf("frame%03d.png", theta)
    rgl::rgl.snapshot(frame_name)
    frame_files = c(frame_files, frame_name)
  }

  # Convert to GIF
  img_list = lapply(list.files(pattern = "frame", full.names = TRUE),
                    magick::image_read)
  img = magick::image_join(img_list)
  img_animated = magick::image_animate(img, fps = fps)
  magick::image_write(img_animated, file)
  # Clean up
  rgl::close3d(#dev = cur3d(),
    silent = TRUE)
  file.remove(frame_files)
  gc()
}
