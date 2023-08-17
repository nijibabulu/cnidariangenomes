subset_msynt <- function(ms, sp, include_chroms) {
  ms_flt <-
    ms %>% dplyr::filter(species == sp)
  if (is.null(include_chroms)) {
    ms_flt
  } else {
    ms_flt %>% dplyr::filter(chrom %in% include_chroms)
  }
}

analyze_macrosynteny <- function(a, b, data, adjust_method = "BH", a_chroms = NULL, b_chroms = NULL) {
  ms_a <- subset_msynt(data, a, a_chroms)
  ms_b <- subset_msynt(data, b, b_chroms)

  alg_genes <- auto_detect_alg_genes.msynt(ms_a, ms_b, adjust_method = adjust_method)

  list(ms_a, ms_b, alg_genes)
}

# use_alg uses the alg column in the data rather than computing them
plot_ms <- function(a, b, data, use_alg = T, a_chroms = NULL, b_chroms = NULL,
                    adjust_method = "BH", remove_chrom_prefix = T,
                    xtitle = a, ytitle = b, min_scaffold_size = 0,
                    show_chrom_names = T, drop_genes_missing_algs = T,
                    clust_method = macrosynteny_cluster, ...) {
  if(drop_genes_missing_algs) data <- dplyr::filter(data, !is.na(alg))

  c(ms_a, ms_b, alg_genes) %<-% analyze_macrosynteny(a, b, data, adjust_method, a_chroms, b_chroms)

  if(remove_chrom_prefix) {
    ms_a$chrom <- stringr::str_remove(ms_a$chrom, "^.*[.]")
    ms_b$chrom <- stringr::str_remove(ms_b$chrom, "^.*[.]")
  }

  ranks <- rank_genes(ms_a, ms_b,  clust_method = clust_method)
  big_chroms.x <- dplyr::count(ranks, chrom.x) %>% dplyr::filter(n > min_scaffold_size) %>% dplyr::pull(chrom.x)
  big_chroms.y <- dplyr::count(ranks, chrom.y) %>% dplyr::filter(n > min_scaffold_size) %>% dplyr::pull(chrom.y)
  ranks <- dplyr::filter(ranks, chrom.x %in% big_chroms.x & chrom.y %in% big_chroms.y)

  if(use_alg) alg_genes <- data

  plot_macrosynteny_dots(ranks, alg_genes, show_chrom_names = show_chrom_names,
                         xtitle = xtitle, ytitle = ytitle, ...) +
    ggplot2::theme(axis.title = ggplot2::element_text(size=12, face = "italic"),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size=7),
                   plot.tag = ggplot2::element_text(size=ggplot2::rel(2)),
                   plot.margin = ggplot2::unit(c(10,5.5,10,5.5), "points"))
}
compute_alg_links <- function(a, b, data, p_adjust_method = "BH", c_adjust_method = "ci_mm") {
  c(ms_a, ms_b, alg_genes) %<-% analyze_macrosynteny(a,b,data,p_adjust_method)
  cmat <- compute_cvalues(ms_a, ms_b, adjust = c_adjust_method, alg_genes = alg_genes)
  cmat_ci <- compute_cvalues(ms_a, ms_b, adjust = "ci", alg_genes = alg_genes)
  cmat_raw <- compute_cvalues(ms_a, ms_b, adjust = "none", alg_genes = alg_genes)

  alg_genes %>% dplyr::count(chrom.x, chrom.y) %>%
    dplyr::mutate(cvalue = purrr::pmap_dbl(., function(chrom.x, chrom.y, n) cmat[chrom.x, chrom.y]),
                  cvalue_ci = purrr::pmap_dbl(., function(chrom.x, chrom.y, n) cmat_ci[chrom.x, chrom.y]),
                  cvalue_raw =purrr::pmap_dbl(., function(chrom.x, chrom.y, n) cmat_raw[chrom.x, chrom.y]) )
}

collapse_codes <- function(codes_df) {
  codes_counts <- dplyr::count(codes_df, chr_list) %>% dplyr::arrange(dplyr::desc(n))
  collapsed_codes_list <-
    purrr::reduce2(codes_counts$chr_list, codes_counts$n,
                 function(acc, chr_list, n) {
                   if(all(is.na(chr_list))) {
                     acc
                   } else {
                     is <- which(purrr::map_lgl(acc, ~all(na.omit(.x[[1]] == chr_list))))
                     if(length(is) == 0)
                       c(acc, list(list(chr_list, n)))
                     else {
                       i <- is[1]
                       c(cur_list,cur_n) %<-% acc[[i]]
                       if(length(na.omit(chr_list)) > length(na.omit(cur_list)))
                         cur_list <- chr_list
                       acc[[i]] <- list(cur_list, cur_n + n)
                       acc
                     }
                   }
                 }, .init = list())
  collapsed_codes_list %>%
    purrr::transpose() %>%
    purrr::set_names("code_obj", "count") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(count = unlist(count))
}
compute_mutual_connectivity <- function(alg_graph) {
  mutual_connectivity <- igraph::E(alg_graph) %>%
    purrr::map_int(~igraph::ends(alg_graph, .x) %>%
                     igraph::ego(graph = alg_graph, nodes = .) %>%
                     purrr::lift(igraph::graph.intersection)() %>%
                     length())

  alg_graph <- igraph::set_edge_attr(alg_graph, "mutual_connectivity", value = mutual_connectivity) #
}


compute_alg_graph <- function(alg_links, exclude_vertices=c()) {
  alg_graph <- alg_links %>%
    igraph::graph_from_data_frame(directed = F) %>%
    igraph::set_edge_attr("edge_betweenness", value = igraph::edge_betweenness(.)) %>%
    igraph::set_edge_attr("edge_connectivity", value = igraph::edge_connectivity(.))

  # add the "mutual connectivity" of each edge; remove edges which join seaparte communities
  alg_graph <- compute_mutual_connectivity(alg_graph)

  compute_species_cliques(alg_graph)
}
compute_gene_edges <- function(sp, tree, ms, verbose=F) {
  if(verbose) cat("computing edges for", sp[[1]], sp[[2]], fill = T)
  depth <- ape::node.depth(tree)[treeio::MRCA(tree, sp[1], sp[2])]
  ms %>%
    dplyr::filter(species %in% sp) %>%
    dplyr::select(chrom, species, group) %>%
    dplyr::mutate(species = dplyr::case_when(species == sp[1] ~ "A", species == sp[2] ~ "B")) %>%
    tidyr::pivot_wider(group, names_from="species", values_from="chrom") %>%
    tidyr::drop_na() %>%
    dplyr::group_split(A,B) %>%
    purrr::map_dfr(~tidyr::expand_grid(from=.x$group, to=.x$group) %>% dplyr::filter(from < to)) %>%
    dplyr::mutate(depth = depth, species=stringr::str_c(sp, collapse=','))
}
add_node_label <- function(graph, node_label, name = "label") {
   tidygraph::as_tbl_graph(graph) %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate("{name}" := node_label)
}

add_edge_label <- function(tbg, what) {
  nodes <- tidygraph::activate(tbg, nodes) %>% tibble::as_tibble() %>%
    tibble::rowid_to_column("nodeid") %>%
    dplyr::select(nodeid, dplyr::all_of(what))
  name_prepended <- function(name, prefix) rlang::set_names(name, paste0(prefix, name))
  tbg <- tidygraph::activate(tbg, edges) %>%
    dplyr::left_join(dplyr::rename(nodes, dplyr::all_of(name_prepended(what, "to_"))),
                     by=c(to="nodeid")) %>%
    dplyr::left_join(dplyr::rename(nodes, dplyr::all_of(name_prepended(what, "from_"))),
                     by=c(from="nodeid"))
  tbg
}

# add a a characteristics to all nodes in a graph, and apply that
# to the edges if an edge connects two nodes with the same label
add_graph_label <- function(graph, node_label, name="label", fill_value=0) {
  ec <- stringr::str_glue("{c('from','to')}_{name}")
  add_node_label(graph, node_label, name) %>%
    add_edge_label(name) %>%
    tidygraph::activate(edges) %>%
    dplyr::mutate(leiden_group = ifelse(.data[[ec[1]]] == .data[[ec[2]]], .data[[ec[1]]], 0))
}


# create a groups tbl by adding the label from groups to the nodes and then
# converting the names to integers for comparison to msynts
as_groups_tbl <- function(graph, groups, col_name="alg", name_col="name") {
  add_node_label(graph, groups, col_name) %>%
    tibble::as_tibble() %>%
    dplyr::mutate("{name_col}" := as.integer(.data[[name_col]]))
}

filter_iqr_outliers <- function(groups, which=c("low", "high", "both")) {
  which <- match.arg(which)
  szs <- purrr::map_int(groups, length)
  bps <- boxplot.stats(szs)
  out <- c()
  if(which != "low") out <- c(out, bps$out[bps$out > bps$stats[5]])
  if(which != "high") out <- c(out, bps$out[bps$out < bps$stats[1]])
  groups[!names(groups) %in% names(out)]
}

compute_group_overlap_mat <- function(from, to, genomes_og,
                              seriation_method="TSP", normalization=c("row", "column", "outer"),
                              min_scaffold = 0, remove_self_edges = T, return_normalized = F) {
  normalization <- match.arg(normalization)
  a <- genomes_og[[from]]
  b <- genomes_og[[to]]
  prefix_a <- stringr::str_split(names(a), '[.]', simplify = T)[,1][1]
  prefix_b <- stringr::str_split(names(b), '[.]', simplify = T)[,1][1]
  c <- c(a, b)
  m <- c %>% purrr::cross2(.,.) %>% purrr::map_int(purrr::compose(length, purrr::lift(intersect))) %>%
    matrix(nrow=length(a)+length(b), dimnames=list(names(c), names(c)))
  m <- m[rowSums(m) > min_scaffold, colSums(m) > min_scaffold]
  mn <- if(normalization == "row") m/rowSums(m)
  else if(normalization == "column") t(t(m)/colSums(m))
  else if(normalization == "outer") m/(outer(rowSums(m), colSums(m), FUN = "+")-m)
  o <- seriation::seriate(dist(mn), method=seriation_method)[[1]] %>% as.integer()
  m <- m[o,o]
  if(remove_self_edges)
    m <- m[stringr::str_detect(rownames(m), prefix_a), stringr::str_detect(colnames(m), prefix_b)]
  if(return_normalized) {
    if(normalization == "row") m/rowSums(m)
    else if(normalization == "column") t(t(m)/colSums(m))
    else if(normalization == "outer") m/(outer(rowSums(m), colSums(m), FUN = "+")-m)
  } else
    m
}
compute_weighted_graph <- function(subtree, msynt, verbose=F, min_weight = 1, weight = c("n_links", "n_species", "total_depth")) {
  if(verbose) cat(paste0(subtree$tip.label, collapse = ","), fill = T)
  weight <- match.arg(weight)
  edge_tbl <-
    subtree$tip.label %>%
    purrr::cross2(.,.,`>=`) %>%
    purrr::map_dfr(compute_gene_edges, subtree, msynt) %>%
    tidyr::separate(species, c("a", "b"), ",") %>%
    tidyr::pivot_longer(c(a,b))  %>%
    dplyr::group_by(from,to)

  if(weight == "n_species") {
    weighted_edge_tbl <-  dplyr::summarize(edge_tbl, weight=dplyr::n_distinct(value), .groups="drop")
  }
  else if(weight == "n_links")  {
    weighted_edge_tbl <-  dplyr::summarize(edge_tbl, weight=dplyr::n(), .groups = "drop")
  }
  else if(weight == "total_depth") {
    weighted_edge_tbl <-  dplyr::summarize(edge_tbl, weight=sum(depth), .groups="drop")
  }

  weighted_edge_tbl %>%
      dplyr::filter(weight >= min_weight) %>%
      igraph::graph_from_data_frame(directed = F)

}
# retrieve the largest component from a graph
largest_component <- function(graph) {
  subcomponents <- igraph::decompose.graph(graph)
  sizes <- purrr::map_int(subcomponents, igraph::vcount)
  subcomponents[order(sizes, decreasing = T)][[1]]
}

# convert a species chromosomes to a list of chromosomes
# represented by their linkage groups
msynt_species_to_groups <- function( species, ms) {
  dplyr::filter(ms, species == !!species) %>%
    dplyr::select(chrom, group)  %>%
    dplyr::group_by(chrom) %>%
    dplyr::group_map(~rlang::set_names(list(as.character(.x$group)), .y)) %>%
    purrr::flatten()
}
overlap_mat_to_edge_df <- function(m, min_weight=-Inf, max_weight=Inf) {
  m %>%
    as.data.frame() %>%
    tibble::rownames_to_column("from") %>%
    tidyr::pivot_longer(-from, names_to="to", values_to="weight") %>%
    dplyr::filter(weight >= min_weight & weight <= max_weight) %>%
    dplyr::mutate(to=factor(to, levels = colnames(m)))
}

plot_genomes_graph <- function(overlap_mats, upper_thresh=.8, lower_thresh=.2, font_sizes = c(3,4), cols=NULL) {
  g <- purrr::map_dfr(overlap_mats, overlap_mat_to_edge_df) %>%
    #g <- purrr::map_dfr(list(me, mc, cn,cr, bp, mb, bb), overlap_mat_to_edge_df) %>%
    #g <- purrr::map_dfr(list(me, mn, mp, mr, mb), overlap_mat_to_edge_df) %>%
    dplyr::filter(weight > lower_thresh) %>%
    dplyr::mutate(alpha=ifelse(weight > upper_thresh, "one-to-one", "split")) %>%
    igraph::graph_from_data_frame() %>%
    tidygraph::as_tbl_graph()
  g <- g %>% tidygraph::activate(nodes) %>%
    dplyr::mutate(in_incidence=igraph::incident_edges(g, igraph::V(g), mode = "in") %>% purrr::map_int(length)) %>%
    dplyr::mutate(root=tidygraph::map_bfs_chr(tidygraph::node_is_source(),
                                              .f = function(node, path, rank, parent,  ...)  {
                                                #cat(stringr::str_glue("{str(path$result[[1]])} {parent}"), fill=T)
                                                if(nrow(path) == 0) tidygraph::.N()$name[node]
                                                else {path$result[[1]]} }))

  g <- add_edge_label(g, "root")
  g <- add_edge_label(g, "in_incidence")

  g <- g %>% tidygraph::activate(edges) %>%
    dplyr::mutate(highlight=to_in_incidence > 1 | alpha == "split")


  if(is.null(cols)) {
    n_roots <- g %>% tidygraph::activate(edges) %>%
      dplyr::pull(from_root) %>% dplyr::n_distinct()
    cols <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(n_roots)
  }

  roots <- g %>% tidygraph::activate(edges) %>% dplyr::pull(from_root) %>% as.factor() %>%
    forcats::fct_reorder( stringr::str_remove(., "^..") %>% as.integer())
  g <- g %>% tidygraph::activate(nodes) %>%
    dplyr::mutate(ancestral = stringr::str_detect(name, "^[[:upper:]]"),
                  fontface = ifelse(ancestral, "bold", "plain"),
                  size = ifelse(ancestral, 4, 3.8))

  #p <- ggraph::create_layout(sg, "focus", focus = 1) %>%
  p <- ggraph::create_layout(g, "nicely") %>%
    #dplyr::mutate(y=ifelse(stringr::str_detect(name, "em"), 1.5, y)) %>%
    dplyr::rename(x="y", y="x") %>%
    ggraph::ggraph() +
    ggraph::geom_edge_link(ggplot2::aes(width=weight, alpha=highlight, color=factor(from_root))) +#
    ggraph::geom_node_text(ggplot2::aes(label=name, fontface = fontface, size = factor(ancestral))) +
    #ggplot2::scale_color_manual(limits = roots, values = cols) +
    ggraph::scale_edge_color_manual(limits = levels(roots), values = cols) +
    ggraph::scale_edge_alpha_discrete(range = c(0.3, 1.0), guide = F) +
    ggraph::scale_edge_width_continuous(range = c(1,3)) +
    ggplot2::scale_size_manual(values=font_sizes) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
                   legend.position = "none")
  p
}
plot_ancestral_genome_composition <- function(genome_name, anc_genome_name, genome_ogs, cols = NULL, min_overlap = 10, alpha = 1, font_size=15, type = c("normal", "float", "invert"), space = c(0,0)) {
  type <- match.arg(type)
  if(genome_name == anc_genome_name) {
    m <- diag(lengths(genome_ogs[[genome_name]]))
    dimnames(m) <- list(names(genome_ogs[[genome_name]]), names(genome_ogs[[genome_name]]))
  } else
    m <- compute_group_overlap_mat(genome_name, anc_genome_name, genome_ogs)

  data <-
    m %>%
    as.data.frame() %>%
    tibble::rownames_to_column("c.from") %>%
    tidyr::pivot_longer(-c.from, names_to="c.to", values_to="count") %>%
    # for self comparisons, we get rid of links that are there twice
    dplyr::distinct(c.from, c.to, .keep_all = T) %>%
    dplyr::filter(count >= min_overlap) %>%
    dplyr::mutate(dplyr::across(starts_with("c."),  ~as.integer(stringr::str_remove(.x, "^.[.]"))))

  # order by the ancestral genome chromosomal order
  genome_order <- dplyr::arrange(data, c.to, dplyr::desc(count)) %>%
    dplyr::distinct(c.from, .keep_all = T) %>%
    dplyr::pull(c.from)
  data$c.from <- factor(data$c.from, genome_order)
  data$c.to <- factor(data$c.to)

  alpha_hex = format(as.hexmode(as.integer(alpha*255)), upper.case=T)
  if(is.null(cols))  {
    cols <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(genome_ogs[[anc_genome_name]]))
  }
  if(all(stringr::str_length(cols) == 7)) {
    alpha_hex = format(as.hexmode(as.integer(alpha*255)), upper.case=T)
    cols <- stringr::str_c(cols, alpha_hex, sep = "")
  }

  plot_data <- data %>%
    dplyr::arrange(c.to, dplyr::desc(count)) %>%
    dplyr::group_by(c.from) %>%
    dplyr::mutate(xmin = as.numeric(c.from), xmax=xmin+.9, ymax=cumsum(count), ymin=dplyr::lag(ymax, default=0))
  if(type == "float")
    plot_data <- plot_data %>% dplyr::mutate(ymin = -ymax/2, ymax = ymax/2)
  if(type == "invert")
    plot_data <- plot_data %>% dplyr::mutate(ymin = -ymin, ymax = -ymax)

  ggplot2::ggplot(plot_data, ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=factor(c.to))) +
    ggplot2::geom_rect() +
    ggplot2::geom_text(ggplot2::aes(x=(xmin+xmax)/2,y=(ymin+ymax)/2,label=c.to)) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_y_continuous(expand = c(0,0)) + #, limits = c(min(plot_data$ymin), max(plot_data$ymax))) +
    ggplot2::xlim(-space[1], max(plot_data$xmax)+space[2]) +
    ggplot2::theme_void() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), legend.position = "none")
}

plot_cluster_internals <- function(clustering_info, title, tags=NULL, gamma_range=c(-Inf, Inf)) {
  df <- clustering_info$modularity_data
  pgamma <- df %>%
    dplyr::filter(gamma >= gamma_range[1], gamma <= gamma_range[2]) %>%
    dplyr::mutate(gamma = sprintf("%.2f", gamma)) %>%
    ggplot2::ggplot(ggplot2::aes(x=factor(gamma), y = modularity)) +
    #ggplot2::geom_violin(scale="count") +
    ggplot2::stat_summary(fun.data=ggplot2::mean_sdl, fun.args = list(mult=1), geom="errorbar") +
    ggplot2::stat_summary(fun.y=mean, geom="point", color="red") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = stringr::str_glue("{title} Modularity by Resolution Parameter"), x = "Resolution Parameter", y = "Modularity",
                  tag = if(is.null(tags)) ggplot2::waiver() else tags[1]) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0))
  pp <- preplot(clustering_info$gamma_fit)
  pptbl <- tibble::tibble(`p, Frequency`=pp$fit, `Resolution Parameter`=pp$xev[[1]])
  htbl <- tibble::tibble(gamma_samples=clustering_info$gamma_samples)
  psample <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x=`Resolution Parameter`,y=`p, Frequency`), data = pptbl) +
    ggplot2::geom_histogram(ggplot2::aes(x=gamma_samples, y=ggplot2::after_stat(count)/sum(count)),
                            data=htbl, override.aes = T, fill = "red", alpha=.7, bins = 30) +
    ggplot2::labs(title = stringr::str_glue("{title} Resolution Parameter Sampling"),
                  tag = if(is.null(tags)) ggplot2::waiver() else tags[2]) +
    ggplot2::theme_bw()

  pfix <- tibble::tibble(`Number of Clusters`=factor(clustering_info$cluster_sizes),
                         Modularity=clustering_info$cluster_modularity) %>%
    ggplot2::ggplot(ggplot2::aes(x=`Number of Clusters`, y=Modularity)) +
    #ggplot2::geom_violin(scale="count") +
    ggplot2::stat_summary(fun.data=ggplot2::mean_sdl, fun.args = list(mult=1), geom="errorbar") +
    #ggplot2::stat_summary(fun.y=mean, geom="point", color="red") +
    ggplot2::stat_summary(fun.y=mean, geom="point", color="red") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = stringr::str_glue("{title} Partition Modularity vs No. of ALGs"),
                  tag = if(is.null(tags)) ggplot2::waiver() else tags[3])

  patchwork::wrap_plots(pgamma, psample, pfix, ncol=1)
}

plot_orthologs_summary <- function(ms) {
 p1 <- ms %>% dplyr::group_by(`Number of Species`=factor(group_size) %>% forcats::fct_rev()) %>%
   dplyr::summarize(`Number of Orthologous Groups` = dplyr::n()) %>%
    ggplot2::ggplot(ggplot2::aes(x=`Number of Species`, y=`Number of Orthologous Groups`)) +
    ggplot2::geom_col()  +
   ggplot2::coord_flip() +
    ggplot2::theme_bw()
    species_names = c(Aa="A. aurita", Am="A. millepora", Aq = "A. queenslandica", Bf="B. floridae",
                      Ch="C. hemisphaerica", Cr= "C. rotundicauda", Em="E. muelleri",
                      Ep= "E. pallida", Hm="H. magnipapillata", Lo="L. oculatus",
                      Nv="N. vectensis", Py = "P. yessoensis", Re="R. esculentum",
                      Ta = "T. adhaerens", Tc="T. corallinus", Xs="Xenia sp.",
                      Hv="H. vulgaris")
  p2 <- ms %>% dplyr::count(sp=factor(species, levels=c("Em", "Aq", "Ta", "Nv", "Ep", "Am", "Xs", "Re", "Aa", "Ch", "Hm", "Tc", "Cr", "Py", "Bf", "Lo", "Hv"))) %>%
    dplyr::mutate(Species = factor(species_names[as.character(sp)], levels = species_names[as.character(sp)])) %>%
     ggplot2::ggplot(ggplot2::aes(x=Species, y = n)) +
     ggplot2::geom_col() +
     ggplot2::theme_bw() +
    ggplot2::labs(y = "Number of Orthologous Groups") +
     ggplot2::coord_flip()
  patchwork::wrap_plots(p1,p2, ncol=1) + patchwork::plot_annotation(tag_levels = "a")
}
msynt_append <- function(ms, pos, orthogroups, sp, link_sp,
                         og_fmt = c("oma", "rbh"), max_group_size = 2) {
  og_fmt <- match.arg(og_fmt)
  if(og_fmt == "oma") {
    orthogroups <- orthogroups %>%
      dplyr::filter(species %in% c(sp, link_sp)) %>%
      dplyr::select(-species)
  }
  og_flt <- orthogroups %>%
    dplyr::group_by(group) %>%
    dplyr::filter(dplyr::n() <= 2) %>%
    dplyr::ungroup()

  mapped_ogs <-
    dplyr::inner_join(ms, dplyr::rename(og_flt, rbh_group=group), by="Geneid") %>%
    dplyr::select(group,rbh_group,alg) %>%
    dplyr::inner_join(og_flt, by=c(rbh_group="group")) %>%
    dplyr::anti_join(ms, by="Geneid") %>%
    dplyr::select(-rbh_group)

  create_msynt(pos[[sp]], dplyr::select(mapped_ogs, -alg), mapped_ogs) %>%
    dplyr::mutate(species=sp) %>%
    dplyr::bind_rows(ms)
}
