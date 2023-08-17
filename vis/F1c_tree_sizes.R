size_plot <- function(tree, genome_sizes, node_labels = NULL, rotations = list(),
                      right_expansion = .25, tip_label_offset = 2.5,
                      clade_label_offset = 5, treecolor = "grey60",
                      bar_expand_mult = c(0.01, .05), rootedge = .1,
                      tiplab_size = 3.0, color_evidence = T, ...) {
  p.tree <- plot_phylogenetic_tree(tree,
    rotations = rotations,
    node_labels = node_labels,
    right_expansion = right_expansion,
    tip_label_offset = tip_label_offset,
    tiplab_size = tiplab_size,
    clade_label_offset = clade_label_offset,
    treecolor = treecolor,
    rootedge = rootedge, ...
  ) +
    ggplot2::theme(plot.margin = ggplot2::margin(5.5, 0, 5.5, 5.5))


  tree_order <- p.tree$data$label[!is.na(p.tree$data$label)]
  tree_order <- p.tree$data %>%
    dplyr::arrange(y) %>%
    dplyr::filter(!is.na(label)) %>%
    dplyr::pull(label)
  genome_sizes <- genome_sizes %>%
    dplyr::right_join(p.tree$data %>% dplyr::filter(!is.na(label)), by = c(Species = "label"))

  p.tree$data <- p.tree$data %>% dplyr::mutate(label = stringr::str_replace_all(label, "_", " "))
  # genome_sizes <- genome_sizes %>%
  # dplyr::arrange(match(Species, tree_order)) %>%
  # dplyr::mutate(Species = factor(Species, levels=Species))
  bar_aes <- ggplot2::aes(x = y, y = Length)
  if (color_evidence) bar_aes <- utils::modifyList(bar_aes, ggplot2::aes(fill = Evidence))
  p.bar <-
    ggplot2::ggplot(genome_sizes, bar_aes) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = bar_expand_mult)) +
    ggthemes::scale_fill_ptol() +
    ggplot2::labs(y = "Genome Size (MB)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 0)
    )

  patchwork::wrap_plots(p.tree, p.bar, nrow = 1, guides = "collect", design = "AAB") &
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2)) &
    ggplot2::theme(legend.position = "bottom")
}

pelagiidae <- "Sanderia_malayensis"
ulmaridae <- "Aurelia_aurita"
rhizostomeae <- "(Nemopilema_nomurai, Rhopilema_esculentum)"
scyphozoa <- stringr::str_glue("({pelagiidae}, ({ulmaridae}, {rhizostomeae}))")
# C. xamachana, A.aurita relationship based on
# https://en.wikipedia.org/wiki/Scyphozoa: Ulmaridae is a separate order,
# Cassiopeidae and Rhizostomeae are sibling orders.

cubozoa <- "Alatina_alata"
staurozoa <- "Calvadosia_cruxmelitensis"
hydrozoa <- "(Clytia_hemisphaerica, (Hydra_veridissima, Hydra_vulgaris))"
hydrozoa <- "(Clytia_hemisphaerica, Hydra_vulgaris)"
medusozoa <- stringr::str_glue("({hydrozoa}, {scyphozoa})")
actiniaria <- "((Scolanthus_callimorphus, Nematostella_vectensis), Exaiptasia_pallida)"
robusta <- "(Pocillopora_damicornis, Stylophora_pistillata)"
stony_corals <- stringr::str_glue("((Acropora_digitifera, Acropora_millepora), {robusta})")
octocorallia <- "((Dendronephthya_gigantea, Xenia_sp.), Renilla_muelleri)"
hexacorallia <- stringr::str_glue("({stony_corals},{actiniaria})")
anthozoa <- stringr::str_glue("({octocorallia},{hexacorallia})")
nv.nwk <- stringr::str_glue("({medusozoa},{anthozoa});")
nv.tree <- ape::read.tree(text = nv.nwk)

clades <- list(
  Actinaria = c("Nematostella_vectensis", "Exaiptasia_pallida"),
  Scleractinia = c("Acropora_millepora", "Stylophora_pistillata"),
  Medusozoa = c("Hydra_vulgaris", "Aurelia_aurita"),
  Ecdysozoa = c("Drosophila_melanogaster", "Caenorhabditis_elegans"),
  Chordata = c("Branchiostoma_floridae", "Homo_sapiens")
)

sizes <- readxl::read_xlsx("data/tables/genomesizes.xlsx")





rotations <- list(
  c("Hydra_vulgaris", "Clytia_hemisphaerica"),
  c("Nematostella_vectensis", "Stylophora_pistillata"),
  c("Nematostella_vectensis", "Exaiptasia_pallida")
)

node_labels <- list(
  Actiniaria = c("Nematostella_vectensis", "Exaiptasia_pallida"),
  Scleractinia = c("Acropora_millepora", "Stylophora_pistillata"),
  Anthozoa = c("Nematostella_vectensis", "Dendronephthya_gigantea"),
  Medusozoa = c("Hydra_vulgaris", "Aurelia_aurita"),
  Hexacorallia = c("Stylophora_pistillata", "Exaiptasia_pallida"),
  Octocorallia = c("Renilla_muelleri", "Dendronephthya_gigantea"),
  Hydrozoa = c("Hydra_vulgaris", "Clytia_hemisphaerica"),
  Scyphozoa = c("Aurelia_aurita", "Sanderia_malayensis")
)
p <-
  size_plot(nv.tree, sizes,
    node_labels = node_labels, rotations = rotations,
    node_label_nudge_y = .4, node_label_nudge_x = -.01, right_expansion = .54, clade_label_offset = 0,
    rootedge = 3, bar_expand_mult = c(0.015, .05)
  )


save_fig(f1, prefix = "F1c_tree_sizes", width = 6, height = 7)
