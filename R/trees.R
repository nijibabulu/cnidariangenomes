# this computes branch lengths for a tree, which is required for adding a
# root node. since compute.brlen uses a the proportional depth (method = 1L,
# last argument to node_depth), I adapted this to the "equal length" method
compute.cg_brlen <- function(phy, power=1, scale.lengths=T) {
  tr <- reorder(phy, "postorder")
  Ntip <- length(phy$tip.label)
  Nedge <- dim(phy$edge)[1]
  Nnode <- phy$Nnode
  xx <- .C(ape::node_depth, as.integer(Ntip),
           as.integer(tr$edge[, 1]), as.integer(tr$edge[, 2]), as.integer(Nedge),
           double(Ntip + Nnode), 2L)[[5]] - 1
  m <- if(scale.lengths) { Ntip - 1 } else { 1 }
  phy$edge.length <- (xx[phy$edge[, 1]]/m)^power - (xx[phy$edge[,  2]]/m)^power
  phy
}

clade_label_node <- function(a, b, .tree) {
  if(missing(b) || a == b)
    tidytree::nodeid(.tree, a)
  else
    tidytree::MRCA(.tree, a, b)
}

# tree must be a phylo object to add root edge
plot_phylogenetic_tree <- function(tree, clade_labels = list(), rotations=list(), node_labels = NULL,
                                   node_label_nudge_x = 0,  node_label_nudge_y = 0, treecolor = "grey15",
                                   treelayout = "rectangular", tree_rotation = NULL,
                                   tiplab_fontface = "italic", tiplab_size=3.8,
                                   nodelab_size = tiplab_size,
                                   clade_label_offset=5, tip_label_offset = 0,
                                   right_expansion = .5, rootedge = NULL, group_colors = NULL, tag = "", ...) {

  if(!is.null(rootedge)) {
    if(is.null(tree$edge.length)) {
      tree <- compute.cg_brlen(tree, scale.lengths = F)
    }
  }

  if(is.null(group_colors)) {
    p <- ggtree::ggtree(tree, ggplot2::aes(x, y), color = treecolor, layout = treelayout, ...)
  } else {
    p <- ggtree::ggtree(tree, ggplot2::aes(x, y, color = group), layout = treelayout, ...)
  }


  p <- p +
    ggtree::geom_tiplab(fontface = tiplab_fontface, size = tiplab_size) +
    #ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult = c(.1,right_expansion))) +
    ggplot2::scale_x_continuous(expand = c(0.1, 0, right_expansion, 0)) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::labs(tag = tag) +
    ggtree::theme_tree2(axis.text.x = ggplot2::element_blank(),
                        axis.ticks.x = ggplot2::element_blank(),
                        axis.line.x = ggplot2::element_blank(),
                        legend.position = "none")

  #p <- ggplot2::ggplot(tree, aesthetic) +
    #ggtree::geom_tree() +
    #ggtree::geom_tiplab(fontface = "italic", offset = tip_label_offset, hjust = tip_label_hjust)  +
    #ggplot2::scale_x_continuous(expand = c(0.1, 0, right_expansion, 0)) +
    #ggplot2::scale_color_manual(values = group_colors) +
    #ggtree::theme_tree2(axis.text.x = ggplot2::element_blank(),
                        #axis.ticks.x = ggplot2::element_blank(),
                        #axis.line.x = ggplot2::element_blank())

  if(!is.null(rootedge)) {
    p <- p + ggtree::geom_rootedge(rootedge = rootedge, color = treecolor)
  }

  # apply clade labels
  p <- purrr::map(clade_labels, purrr::lift(clade_label_node), .tree=tree) %>%
    purrr::map2(names(clade_labels), ~ggtree::geom_cladelabel(node=.x, label=.y, offset=clade_label_offset)) %>%
    purrr::reduce(`+`, .init=p)

  if(!is.null(tree_rotation)) {
    p <- ggtree::rotate_tree(p, tree_rotation)
  }

  # add node labels
  if(!is.null(node_labels)) {
    p$data <- p$data %>% dplyr::left_join(
      tibble::tibble(
        node = purrr::map_int(node_labels, purrr::lift(tidytree::MRCA), .data=tree) ,
        node_label = names(node_labels)
      ), by="node")

    p <- p + ggtree::geom_text(ggplot2::aes(label=node_label),
                               nudge_x = node_label_nudge_x, nudge_y = node_label_nudge_y,
                               size = nodelab_size, hjust = 1)
    }

  # apply rotations
  p <- purrr::map(rotations, purrr::lift(tidytree::MRCA), .data=tree) %>%
    purrr::reduce(ggtree::rotate, .init=p)

  p
}

plot_macrosynteny_tree <- function(tree, ingroup, outgroup, witness, alg_species=NULL, ...) {
  if(!is.null(alg_species)) {
    alg_node <- tidytree::MRCA(tree, alg_species[1], alg_species[2])
  } else {
    alg_node <- tidytree::MRCA(tree, ingroup, outgroup)
  }
  tree <- ggtree::groupOTU(tree, .node=c(ingroup, outgroup, witness))

  plot_phylogenetic_tree(tree, group_colors = c("black", "red"), ...) +
    ggtree::geom_point2(ggplot2::aes(subset=node == alg_node), color = "black")

}
