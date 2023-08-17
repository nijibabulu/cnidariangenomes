
low.tree <- treeio::read.newick("data/clock/low.nwk")
high.tree <- treeio::read.newick("data/clock/high.nwk")

mean.lengths <- (low.tree$edge.length + high.tree$edge.length)/2
range.tree <- low.tree
range.tree$edge.length <- mean.lengths

data <- low.tree %>% treeio::as.treedata() %>% tibble::as_tibble() %>%
  dplyr::full_join(high.tree %>% treeio::as.treedata() %>% tibble::as_tibble(), by=c("parent", "node", "label"),
                   suffix = c(".low", ".high")) %>%
  dplyr::rename(labels=label) %>%
  dplyr::mutate(length = (branch.length.low + branch.length.high)/2,
                length_range = purrr::transpose(list(branch.length.low, branch.length.high)))

parent_data <- data %>% dplyr::select(node=parent, node_range=length_range) %>% dplyr::left_join(data) %>%
  dplyr::group_by(node) %>% dplyr::mutate(node_rangea=purrr::map_dbl(node_range, purrr::pluck, 1),
                                   node_rangeb=purrr::map_dbl(node_range, purrr::pluck, 2)) %>%
  dplyr::filter(node_rangea-node_rangeb == max(node_rangea-node_rangeb, na.rm = T))
tree <- treeio::treedata(phylo = range.tree, data = tibble::as_tibble(parent_data))

cb_split <- tidytree::MRCA(tree, treeio::nodeid(tree, "re"), treeio::nodeid(tree, "hs"))
cb_range <- c(595.7, 688.3)
cb_delta <- mean(cb_range)-cb_range
cb_branch_length <- tree@data[tree@data$node == cb_split,]$length
tree@data[tree@data$node == cb_split,]$node_range <- list(as.list(cb_branch_length-cb_delta))

names <- tibble::tribble(~label, ~name,
                         "el", "E. lineata",
                         "ec", "E. carnea",
                         "xs", "Xenia sp.",
                         "sc", "S. callimorphus",
                         "nv", "N. vectensis",
                         "ms", "M. senile",
                         "ep", "E. pallida",
                         "ad", "A. digitifera",
                         "am", "A. millepora",
                         "aa", "A. aurita",
                         "re", "R. esculentum",
                         "hm", "H. magnipapillata",
                         "ch", "C. hemisphaera",
                         "bf", "B. floridae",
                         "lo", "L. oculatus",
                         "dm", "D. melanogaster",
                         "ce", "C. elegans",
                         "em", "E. muelleri",
                         "hs", "H. sapiens")
tree@phylo$tip.label <-
  tibble::tibble(label=tree@phylo$tip.label) %>%
  dplyr::left_join(names, by="label") %>%
  dplyr::pull(name)

length_to_mya <- function(x) {
  round(max(tree@phylo$edge.length, na.rm = T)-x, digits=2)
}

p <- ggtree::ggtree(tree, right=T) +
  ggtree::geom_tiplab(align=TRUE, linetype='dashed', linesize=.3) +
  ggtree::geom_range("node_range", color='red', size=2, alpha=.5) +
  ggtree::theme_tree2() +
  ggplot2::scale_x_continuous(
    breaks = seq(36.78-100,843.2,100),
    labels = length_to_mya,
    expand = c(0.1, 0, .4, 0)) +
  ggplot2::labs(x="Age, Mya")


save_fig(p, "EDF1_clock", width = 6, height = 4)
ggsave("figures/EDF1_clock.pdf", width=6, height=4)
ggsave("figures/EDF1_clock.png", width=6, height=4)

