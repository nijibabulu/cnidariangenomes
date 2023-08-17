
# compute or retrieve cache of computed ALGs
# NB: if not computed yet, this takes a long time and a lot of resources!
infer_cached_algs(
  "((Em,Aq),(Ta,((((((Nv,Ep),Am),Xs),(Aa,Re)),(Ch,Hv)),((Tc,Cr),(Py,(Bf,Lo))))));",
  critical_nodes = list(B="Tc,Cr,Py,Bf,Lo",
                        C="Nv,Ep,Am,Xs,Aa,Re,Ch,Hv",
                        M="Em,Aq,Ta,Nv,Ep,Am,Xs,Aa,Re,Ch,Hv,Tc,Cr,Py,Bf,Lo"),
  ortholog_source = "data/orthologs/OrthologousGroups.fixed.txt")

# convert the alg orthogroups for a list for use in the matrix-based functions below.
alg_ogs_list <- alg_ogs %>%
  purrr::map(~dplyr::group_by(.x, name = alg) %>%
               dplyr::summarize(value = list(group)) %>%
               tibble::deframe())

# also create extant genomes lists for the graph comparison
extant_genome_ogs_list <- rlang::set_names(species_big, species) %>%
  purrr::map(msynt_species_to_groups, msynts$M)

# the colors used are synced complementary to the ALG colors in figure 2a
n_clades <- 5
n_cols_main <- dplyr::n_distinct(alg_ogs$M$alg) + n_clades
cols_main <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(n_cols_main)
cols_tree_inds <- seq(4, n_cols_main, length.out = n_clades) %>% round()
cols_M <- cols_main[-cols_tree_inds]

# generate ALG composition figure
p_composition <- patchwork::wrap_plots(
  plot_ancestral_genome_composition("C", "M", alg_ogs_list, min_overlap = 50, cols=cols_M, font_size=13),
  plot_ancestral_genome_composition("M", "M", alg_ogs_list, min_overlap = 50, cols=cols_M, font_size=13, type = "float", space = c(.5,.5)),
  plot_ancestral_genome_composition("B", "M", alg_ogs_list, min_overlap = 50, cols=cols_M, font_size=13, type = "invert", space = c(.5,.5)),
  ncol = 1
) & ggplot2::theme(plot.margin=grid::unit(c(10,0,10,0),"pt"))
save_fig(p_composition, "F2c_composition",  width=6, height=6)

# generate a list of comparisons to use in the graph figure
c_chroms <- c("nv", "hv", "xs")
p_chroms <- c("em")
b_chroms <- c("py", "bf", "tc")

genome_comparisons <-
  tibble::tribble(~from, ~to,
                  "M", tibble::tibble(to=c("B","C", "em")),
                  "C", tibble::tibble(to=c("nv", "hv")),
                  "B", tibble::tibble(to="bf")) %>%
    tidyr::unnest()

# analyze the overlap of the OGs from each node
overlap_mats <-
  purrr::map2(genome_comparisons$from, genome_comparisons$to,
              compute_group_overlap_mat, c(alg_ogs_list, extant_genome_ogs_list),
              return_normalized = T, normalization = "column") %>%
  rlang::set_names(stringr::str_glue_data(genome_comparisons, "{from}_{to}"))


p_graph_composition <-
  plot_genomes_graph(overlap_mats, lower_thresh = .1, font_sizes = c(4,4.5), cols = cols_M)
save_fig(p_graph_composition, "F2c_graph_composition", types = "pdf", width=8, height=8)
