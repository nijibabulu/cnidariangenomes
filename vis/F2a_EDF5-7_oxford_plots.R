main_figure <- function(ms, ...) {
  a <- plot_ms("Em", "Nv", ms, a_chroms = stringr::str_glue("em.{seq(23)}"), xtitle = "Ephydatia muelleri", ...)
  b <- plot_ms("Sc", "Nv", ms, a_chroms = stringr::str_glue("sc.{seq(15)}"), xtitle = "Scolanthus callimorphus", ...)
  c <- plot_ms("Xs", "Nv", ms, a_chroms = stringr::str_glue("xs.{seq(15)}"), xtitle = "Xenia sp.", ...)
  d <- plot_ms("Re", "Nv", ms, xtitle = "Rhopilema esculentum", ...)
  e <- plot_ms("Sb", "Nv", ms,  a_chroms = stringr::str_glue("sb.{seq(11)}"), xtitle = "Streblospio benedicti", ...)
  f <- plot_ms("Tc", "Nv", ms, min_scaffold_size = 50, xtitle = "Trigoniulus corallinus", ...)
  g <- plot_ms("Dm", "Nv", ms, min_scaffold_size = 50, xtitle = "Drosophila melanogaster", ...)
  h <- plot_ms("Hv", "Nv", ms, xtitle = "Hydra vulgaris", ...)
  i <- plot_ms("Py", "Nv", ms, xtitle = "Patinopecten yessoensis", ...)
  j <- plot_ms("Bf", "Nv", ms, xtitle = "Branchiostoma floridae", ...)
  k <- plot_ms("Lo", "Nv", ms,  xtitle = "Lepisosteus oculatus", min_scaffold_size = 50, ...)
  l <- plot_ms("Hs", "Nv", ms,  xtitle = "Homo sapiens", a_chroms=  c("X", stringr::str_glue("hs.{seq(22)}")), min_scaffold_size = 50, ...)
  patchwork::wrap_plots(a,b,c,d,e,f,g,h,i,j,k,l, ncol=4)
}

sup_figure_M <- function(ms, anc = "M", ...) {
  a <- plot_ms("Em", anc, ms, a_chroms = stringr::str_glue("em.{seq(23)}"), xtitle = "Ephydatia muelleri",  ...)
  if(anc == "Nv") {
    b <- plot_ms("Sc", anc, ms, a_chroms = stringr::str_glue("sc.{seq(15)}"), xtitle = "Scolanthus callimorphus", ...)
  } else {
    b <- plot_ms("Nv", anc, ms, a_chroms = stringr::str_glue("nv.{seq(15)}"), xtitle = "Nematostella vectensis", ...)
  }
  c <- plot_ms("Xs", anc, ms, a_chroms = stringr::str_glue("xs.{seq(15)}"), xtitle = "Xenia sp.", ...)
  d <- plot_ms("Re", anc, ms, xtitle = "Rhopilema esculentum", ...)
  e <- plot_ms("Sb", anc, ms,  a_chroms = stringr::str_glue("sb.{seq(11)}"), xtitle = "Streblospio benedicti", ...)
  f <- plot_ms("Cr", anc, ms,  a_chroms = stringr::str_glue("cr.{seq(16)}"), xtitle = "Carcinoscorpius rotundicauda", ...)
  g <- plot_ms("Tc", anc, ms, min_scaffold_size = 50, xtitle = "Trigoniulus corallinus", ...)
  h <- plot_ms("Dm", anc, ms, min_scaffold_size = 50, xtitle = "Drosophila melanogaster", ...)
  i <- plot_ms("Py", anc, ms, xtitle = "Patinopecten yessoensis", ...)
  j <- plot_ms("Bf", anc, ms, xtitle = "Branchiostoma floridae", ...)
  k <- plot_ms("Lo", anc, ms,  xtitle = "Lepisosteus oculatus", min_scaffold_size = 50, ...)
  l <- plot_ms("Hs", anc, ms,  xtitle = "Homo sapiens", a_chroms=  c("X", stringr::str_glue("hs.{seq(22)}")),
               min_scaffold_size = 50, show_legend = T,  ...)
  patchwork::wrap_plots(a,b,c,d,e,f,g,h,i,j,k,l, ncol=4, guides = "collect") & ggplot2::theme(legend.position = "bottom")
}

sup_figure_C <- function(ms, anc = "C", ...) {
  cm <- seriate_macrosynteny_cluster_method(method = "TSP")
  if(anc == "Nv") {
    a <- plot_ms("Sc", anc, ms, a_chroms = stringr::str_glue("sc.{seq(15)}"), xtitle = "Scolanthus callimorphus",  ...)
  } else {
    a <- plot_ms("Nv", anc, ms, a_chroms = stringr::str_glue("nv.{seq(15)}"), xtitle = "Nematostella vectensis",  ...)
  }
  b <- plot_ms("Ep", anc, ms, xtitle = "Exaiptasia pallida", show_chrom_names = c(F, T), linesize = .1, clust_method = cm,  ...)
  c <- plot_ms("Am", anc, ms, xtitle = "Acropora millepora", show_chrom_names = c(F, T), linesize = .1, clust_method = cm, ...)
  d <- plot_ms("Xs", anc, ms, a_chroms = stringr::str_glue("xs.{seq(15)}"), xtitle = "Xenia sp.",  ...)
  e <- plot_ms("Re", anc, ms, xtitle = "Rhopilema esculentum",  ...)
  f <- plot_ms("Aa", anc, ms, xtitle = "Aurelia aurita", show_chrom_names = c(F, T), linesize = .1, clust_method = cm,  ...)
  h <- plot_ms("Hv", anc, ms, xtitle = "Hydra vulgaris", show_legend = T,   linesize = .1, ...)
  patchwork::wrap_plots(a,b,c,d,e,f,h, ncol=4, guides = "collect") & ggplot2::theme(legend.position = "bottom")
}

sup_figure_B <- function(ms, anc = "B", ...) {
  a <- plot_ms("Cr", anc, ms,  a_chroms = stringr::str_glue("cr.{seq(16)}"), xtitle = "Carcinoscorpius rotundicauda",  ...)
  b <- plot_ms("Tc", anc, ms, min_scaffold_size = 50, xtitle = "Trigoniulus corallinus",  ...)
  c <- plot_ms("Dm", anc, ms, min_scaffold_size = 50, xtitle = "Drosophila melanogaster",  ...)
  d <- plot_ms("Ce", anc, ms,  xtitle = "Caenorhabditis elegans",  ...)
  e <- plot_ms("Sb", anc, ms,  a_chroms = stringr::str_glue("sb.{seq(11)}"), xtitle = "Streblospio benedicti",  ...)
  f <- plot_ms("Py", anc, ms, xtitle = "Patinopecten yessoensis",  ...)
  g <- plot_ms("Bf", anc, ms, xtitle = "Branchiostoma floridae",  ...)
  h <- plot_ms("Lo", anc, ms,  xtitle = "Lepisosteus oculatus", min_scaffold_size = 50,  ...)
  i <- plot_ms("Hs", anc, ms,  xtitle = "Homo sapiens", a_chroms=  c("X", stringr::str_glue("hs.{seq(22)}")), min_scaffold_size = 50, show_legend = T,  ...)

  patchwork::wrap_plots(a,b,c,d,e,f,g,h,i,  guides = "collect",
                        design = "abcd
                                  ef##
                                  ghi#") &
    ggplot2::theme(legend.position = "bottom")
}

# compute or retrieve cache of computed ALGs
infer_cached_algs(
  "((Em,Aq),(Ta,((((((Nv,Ep),Am),Xs),(Aa,Re)),(Ch,Hv)),((Tc,Cr),(Py,(Bf,Lo))))));",
  critical_nodes = list(B="Tc,Cr,Py,Bf,Lo",
                        C="Nv,Ep,Am,Xs,Aa,Re,Ch,Hv",
                        M="Em,Aq,Ta,Nv,Ep,Am,Xs,Aa,Re,Ch,Hv,Tc,Cr,Py,Bf,Lo"),
  ortholog_source = "data/orthologs/OrthologousGroups.fixed.txt")

# determine the colors for the main figure and supplemental macrosynteny figure
# total colors is the number of algs plus 5 clades colored in the tree and text
n_clades <- 5
n_cols_main <- dplyr::n_distinct(alg_ogs$M$alg) + n_clades
cols_main <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(n_cols_main)
cols_tree_inds <- seq(4, n_cols_main, length.out = n_clades) %>% round()
cols_tree <- cols_main[cols_tree_inds]
cols_M <- cols_main[-cols_tree_inds]

# colors for the supplemental figures do not include the clades
n_cols_C <- dplyr::n_distinct(alg_ogs$C$alg)
cols_C <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(n_cols_C)

n_cols_B <- dplyr::n_distinct(alg_ogs$B$alg)
cols_B <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(n_cols_B)

# add additional orthologs to groups not directly used in ALG inference for figures
M_msynt <- msynts$M %>%
  msynt_append(pos, oma, "Dm", "Nv") %>%
  msynt_append(pos, oma, "Sc", "Nv") %>%
  msynt_append(pos, oma, "Hs", "Nv") %>%
  msynt_append(pos, oma, "Sb", "Nv") %>%
  msynt_append(pos, oma, "Ce", "Nv")

B_msynt <- msynts$B %>%
  msynt_append(pos, oma, "Dm", "Bf") %>%
  msynt_append(pos, oma, "Sc", "Bf") %>%
  msynt_append(pos, oma, "Hs", "Bf") %>%
  msynt_append(pos, oma, "Sb", "Bf") %>%
  msynt_append(pos, oma, "Ce", "Bf")

C_msynt <- msynts$C %>%
  msynt_append(pos, oma, "Dm", "Nv") %>%
  msynt_append(pos, oma, "Sc", "Nv") %>%
  msynt_append(pos, oma, "Hs", "Nv")  %>%
  msynt_append(pos, oma, "Sb", "Nv")

# generate oxford plot figures for supplement
p_supM <- sup_figure_M(M_msynt, cols = cols_M)
save_fig(p_supM, prefix="EDF7_M_oxford_plots", height=8.5, width=10.5)

p_supB <- sup_figure_B(B_msynt, cols = cols_B)
save_fig(p_supB, prefix="EDF6_B_oxford_plots", height=8.5, width=10.5)

p_supC <- sup_figure_C(C_msynt, cols = cols_C)
save_fig(p_supC, prefix="EDF5_C_oxford_plots", height=6.5, width=10.5)

# generate main oxford plot figure
p_main <- main_figure(M_msynt, cols = cols_M)
save_fig(p_main, prefix="F2a_Nv_oxford_plots", height=8, width=10.5)
