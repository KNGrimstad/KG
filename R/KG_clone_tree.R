#' Construct and plot clonal trees
#'
#' This function uses the dowser and ggtree packages to construct custom clonal trees.
#' @param clones An airrClone object from dowser::formatClones.
#' @param tips The same as in plotTrees. Note: Does not accept factors.
#' @param tip_palette A vector of colors and names specifying the colors for each state.
#' @param build The program used to construct trees (pratchet, pml, dnappars, dnaml, igphyml).
#' @param tipsize The size of tree tips.
#' @param scale The width of branch length scale bar.
#' @param layout The layout of the trees (rectangular, circular, slanted, fan, inward_circular, radial, unrooted, equal_angle, daylight, dendrogram, ape, ellipse, roundrect).
#' @param x_max The maximum value of the x-axis. Default is the maximum of the SHM for each specific clone.
#' @param breaks The breaks for the x-axis.
#' @param common_scale Define if branches should have the same scale across trees.
#' @param stroke Thickness of tip border.
#' @export
#' @examples
#' KG_clone_tree(B_cell_clones)
KG_clone_tree = function(clones,
                         tips = NULL,
                         tip_palette = NULL,
                         build = "pml",
                         tipsize = 7,
                         scale = 0,
                         layout = "rectangular",
                         x_max = NULL,
                         breaks = NULL,
                         common_scale = FALSE,
                         stroke = 0.5){

  suppressPackageStartupMessages({
    require(ggtree)
    require(dowser)
    require(ggplot2)})

  # Initiate empty list to store plots in
  tree_list = list()

  # Specify colors for tips, if NULL
  if(is.null(tip_palette)){

    # For main isotypes
    if(tips == "isotype"){
      tip_palette = c(
        "IgM" = "cornflowerblue",
        "IgD" = "wheat2",
        "IgG" = "brown4",
        "IgA" = "tan2",
        "IgE" = "lightgreen",
        "Germline" = "lightgrey"
      )

    } else if(tips == "c_call"){
      # For Ig subtypes
      tip_palette = c(
        "IGHM"    ="cornflowerblue",
        "IGHD"    ="wheat2",
        "IGHG3"   ="sienna2",
        "IGHG1"   ="brown4",
        "IGHA1"   ="tan1",
        "IGHG2"   ="indianred3",
        "IGHG4"   ="sienna4",
        "IGHE"    ="lightgreen",
        "IGHA2"   ="tan3",
        "Germline"="lightgrey")
    } else{
      tip_palette = NULL
    }
  }

  # Construct trees
  trees = dowser::getTrees(clones,
                           build = build,
                           nproc = 1)

  clone_trees = dowser::plotTrees(trees,
                                  tips = tips,
                                  tip_palette = tip_palette,
                                  tipsize = tipsize,
                                  scale = scale,
                                  layout = layout,
                                  common_scale = common_scale)

  # Plot clone trees
  for(i in seq(nrow(trees))){

    if(is.null(breaks)){
      breaks = as.numeric(formatC(seq(min(trees[[5]][[i]]$edge.length), max(trees[[5]][[i]]$edge.length) + 0.02, 0.02), 1))
    }
    if(layout == "rectangular"){
      tree_list[[i]] = clone_trees[[i]] +
        ggtree::geom_tippoint(shape = 21, size = tipsize, stroke = stroke) +
        ggtree::theme_tree2() +
        ggtree::xlim_tree(x_max) +
        ggplot2::scale_x_continuous(breaks = breaks) +
        ggplot2::ggtitle(paste("Clonotype ", clone_trees[[i]]$labels$title, sep = "")) +
        ggplot2::xlab("SHM frequency") +
        ggplot2::theme(axis.title.x = ggplot2::element_text(hjust = 0.5))
    } else{
      tree_list[[i]] = clone_trees[[i]] +
        ggtree::geom_tippoint(shape = 21, size = tipsize, stroke = stroke) +
        ggtree:theme_tree() +
        ggplot2::ggtitle(paste("Clonotype ", trees[[5]][[i]]$name, sep = "")) +
        ggplot2::xlab("SHM frequency") +
        ggplot2::theme(axis.title.x = ggplot2::element_text(hjust = 0.5))
    }
  }
  return(tree_list)
}
