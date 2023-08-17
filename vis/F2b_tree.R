# compute or retrieve cache of computed ALGs
# NB: if not computed yet, this takes a long time and a lot of resources!
#     We need this only for the colors, so if not interested in this pick a
#     number for the alg_ogs$M$alg value below
infer_cached_algs(
  "((Em,Aq),(Ta,((((((Nv,Ep),Am),Xs),(Aa,Re)),(Ch,Hv)),((Tc,Cr),(Py,(Bf,Lo))))));",
  critical_nodes = list(B="Tc,Cr,Py,Bf,Lo",
                        C="Nv,Ep,Am,Xs,Aa,Re,Ch,Hv",
                        M="Em,Aq,Ta,Nv,Ep,Am,Xs,Aa,Re,Ch,Hv,Tc,Cr,Py,Bf,Lo"),
  ortholog_source = "data/orthologs/OrthologousGroups.fixed.txt")

C <- "((Hydra_vulgaris, (Aurelia_aurita,Rhopilema_esculentum)),(Xenia_sp.,(Acropora_millepora,(Exaiptasia_pallida,(Scolanthus_callimorphus,Nematostella_vectensis)))))"
B <- "(((Homo_sapiens, Lepisosteus_oculatus),Branchiostoma_floridae),((Patinopectin_yessoensis,Streblospio_benedicti),(Carcinoscorpius_rotundicauda,(Trigoniulus_corallinus,Drosophila_melanogaster))))"

tree <- treeio::read.newick(text = stringr::str_glue("
                                      ((({B},{C}),
                                      Trichoplax_adhaerens),
                                      (Ephydatia_muelleri,
                                      Amphimedon_queenslandica));") %>%
                              stringr::str_remove_all("[[:space:]]"))
clade_nodes <-
  c("Anthozoa"=treeio::MRCA(tree, "Nematostella_vectensis", "Xenia_sp."),
    "Medusozoa"=treeio::MRCA(tree, "Aurelia_aurita", "Hydra_vulgaris"),
    "Ecdysozoa"=treeio::MRCA(tree, "Carcinoscorpius_rotundicauda", "Drosophila_melanogaster"),
    "Spiralia"=treeio::MRCA(tree, "Patinopectin_yessoensis", "Streblospio_benedicti"),
    "Chordata"=treeio::MRCA(tree, "Branchiostoma_floridae", "Homo_sapiens"))

tree <- tidytree::groupClade(tree,  clade_nodes)
attr(tree, "group") <- factor(attr(tree, "group") , c(NA, names(clade_nodes)))
tree$tip.label <- stringr::str_replace_all(tree$tip.label, "_", " ")

# the colors used are complementary to the ALG colors in the main figure
n_clades <- 5
n_cols_main <- dplyr::n_distinct(alg_ogs$M$alg) + n_clades
cols_main <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(n_cols_main)
cols_tree_inds <- seq(4, n_cols_main, length.out = n_clades) %>% round()
cols_tree <- cols_main[cols_tree_inds]

p <-
  ggtree::ggtree(tree, ggplot2::aes(color=group), ladderize = F) +
  ggtree::geom_tiplab( hjust=0, offset = -7, align = F, linetype = "dashed", show.legend = F) +
  ggplot2::scale_color_manual(breaks = levels(attr(tree, "group")), values = cols_tree, na.value = "black") +
  ggplot2::scale_x_reverse() +
  ggtree::theme_tree(legend.position = "bottom", legend.title = ggplot2::element_blank())
save_fig(p, prefix="F2b_tree", height=5, width=5)
