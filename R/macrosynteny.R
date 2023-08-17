#' create a msynt [tibble][tibble::tibble-package] for tracking macrosynteny
#'
#' This function creates the main per-species [tibble][tibble::tibble-package]
#' needed to track ancestral linkage groups. This includes 5 columns:
#' \itemize{
#'   \item{chrom}{The chromosome on which the gene sits}
#'   \item{Geneid}{The name of the gene}
#'   \item{position}{The position of the gene on the chromosome}
#'   \item{group}{The orthology group of the gene. This column will need to correspond to
#'                the `group` column of another msynt table}
#'   \item{group_size}{the size of the orthology group. This is useful for filtering
#'                     suspiciously large groups.}
#'   \item{alg}{optionally an ancestral linkage group if provided}
#' }
#'
#' @param gene_positions a [tibble][tibble::tibble-package] of gene positions with the
#'   columns `chrom`, `Geneid` and `position` as documented above
#' @param orthogroups a [tibble][tibble::tibble-package] of gene orthogroups with
#'   columns `Geneid` and `group` as documented above.
#' @param algs a [tibble][tibble::tibble-package] of ancestral linkage groups with the
#'   columns `Geneid` and `alg`.
#' @param max_group_size filter orthogroups whose size are greater than this
#' @param include_chrom a character vector of genes to select the genes on
#' @param min_scaffold_size the minimum number of genes on a scaffold/chromosome to include the
#'    scaffold/chromosome. Overridden by `include_chrom`.
#' @param max_positions restrict the number of times a gene may occur in the genome.
#'    usually this should be 1.
#'
#' @export
#'
create_msynt <- function(gene_positions, orthogroups, algs = NULL,
                             max_group_size = NULL, include_chrom = NULL,
                             min_scaffold_size = NULL, max_positions = 1) {
  group_size <- orthogroups %>% dplyr::count(group, name = "group_size")
  ms <- gene_positions %>%
    dplyr::group_by(Geneid) %>%
    dplyr::filter(dplyr::n() == max_positions) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(orthogroups, by="Geneid") %>%
    dplyr::left_join(group_size, by="group")
  if(!is.null(max_group_size)) {
    ms <- ms %>% dplyr::filter(group_size <= max_group_size)
  }
  if(!is.null(include_chrom)) {
    ms <- ms %>% dplyr::filter(chrom %in% include_chrom)
  }
  if(!is.null(min_scaffold_size)) {
    big_chroms <- ms %>% dplyr::count(chrom) %>% dplyr::filter(n > min_scaffold_size) %>% dplyr::pull(chrom)
    ms <- ms %>% dplyr::filter(chrom %in% big_chroms)
  }
  if(!is.null(algs)) {
    ms <- ms %>% dplyr::left_join(algs %>% dplyr::select(Geneid, alg), by="Geneid")
  } else {
    ms <- ms %>% dplyr::mutate(alg=NA)
  }
  ms
}


#' convert a  [tibble][tibble::tibble-package] of various formats to position format
#'
#' These functions convert a tibble of another shape to a position-shaped [tibble][tibble::tibble-package]
#' to make it compatible with [create_msynt()]
#'
#' @export
psl_to_position <- function(psl) dplyr::select(psl,  chrom="tName", Geneid="qName", position="tStart")

#' @rdname psl_to_position
#' @export
msynt_to_position <- function(msynt) dplyr::select(msynt, chrom, Geneid, position)

#' @rdname psl_to_position
#' @export
exonerate_gff_to_position <- function(gff) dplyr::select(gff, chrom, Geneid="sequence", position="start")

#' @rdname gff_to_position
#' @export
gff_to_position <- function(gff) dplyr::select(gff, chrom, Geneid="geneid", position="start")

#' @rdname gtf_to_position
#' @export
gtf_to_position <- function(gff) dplyr::select(gff, chrom="seqname", Geneid="gene_id", position="start")

#' convert a  [tibble][tibble::tibble-package] of various formats to agora gene list format
#'
#' These functions convert a tibble of another shape to a position-shaped [tibble][tibble::tibble-package]
#' to make it compatible with the AGORA pipeline
#'
#' @export
gff_to_agora <- function(gff) dplyr::select(gff, chrom, position="start", end, strand, Geneid="geneid")

#' @rdname gff_to_agora
#' @export
exonerate_gff_to_agora <- function(gff) dplyr::select(gff, chrom, position="start", end, strand, Geneid="sequence")

#' @rdname gff_to_agora
#' @export
psl_to_agora <- function(psl) dplyr::select(psl, chrom="tName", position="tStart", end="tEnd", strand, Geneid="qName")


#' Write an AGORA [tibble][tibble::tibble-package] in to a path
#'
#' @param agora a [tibble][tibble::tibble-package] containing the columns chrom, Geneid, position, end and strand
#' @param path the output path of the file
#'
#' @export
write_agora <- function(agora, path) {
  agora %>%
    dplyr::mutate(strand = ifelse(strand == "-", -1, 1)) %>%
    readr::write_tsv(col_names = F, file = path)
}

#' Change the chromosome names to the the standard format: <SP>.<CHROM>
#' where <SP> is the two letter canonical species name (e.g. hs for human)
#' and the chromosome number <CHROM>. Chromosome numbers are in descending order of
#' length.
#'
#' @param pos a position-formatted [tibble][tibble::tibble-package]
#' @param prefix the two-letter species prefix
#' @param fasta optionally a fasta file to determine actual chromosome length.
#'
#' @export

make_chrom_names <- function(pos, prefix, fasta=NULL) {
  if(is.null(fasta)) {
    tbl <- pos %>%
      dplyr::group_by(chrom) %>%
      dplyr::summarise(length=max(position)) %>%
      dplyr::arrange(dplyr::desc(length)) %>%
      tibble::rowid_to_column("seqnum") %>%
      dplyr::mutate(new_name=stringr::str_glue("{prefix}.{seqnum}")) %>%
      dplyr::select(name=chrom, new_name)
  } else {
    tbl <- fasta_summary(fasta, compute_stats = F) %>%
      dplyr::mutate(new_name=stringr::str_glue("{prefix}.{seqnum}")) %>%
      dplyr::select(name, new_name)
  }
  pos %>%
    dplyr::left_join(tbl, by=c(chrom="name")) %>%
    dplyr::mutate(chrom=new_name) %>%
    dplyr::select(-new_name) #%>%
    #dplyr::select(chrom, Geneid, dplyr::everything())
}

#' create a matrix containing the number of genes shared in all pairs of chromosomes
#'
#' @param ms.x a msynt [tibble][tibble::tibble-package]
#' @param ms.y a msynt [tibble][tibble::tibble-package]
#' @param metric a column in the msynt (usually either `group` or `alg`) by which to make the matrix
#'
#' @export
macrosynteny_matrix <- function(ms.x, ms.y, metric="group") {
  m <- dplyr::inner_join(ms.x, ms.y, by=metric) %>%
    dplyr::filter(!is.na(!!rlang::ensym(metric))) %>%
    dplyr::count(chrom.x, chrom.y) %>%
    tidyr::pivot_wider(names_from = chrom.y, values_from = n, values_fill = list(n = 0)) %>%
    tibble::column_to_rownames("chrom.x") %>%
    as.matrix()

  m
}

#' plot heatmap of macrosynteny
#'
#' Plots a [ComplexHeatmap::Heatmap] corresponding to the number of genes
#' shared by each pair of chromosomes
#'
#' @param ms.x a msynt [tibble][tibble::tibble-package]
#' @param ms.y a msynt [tibble][tibble::tibble-package]
#' @param metric a column in msynt (usually either `group` or `alg`) by which the pairs of chromosomes are measured
#' @param scale whether to scale the initial matrix. This works better for comparing to scaffolded genomes.
#' @param ... arguments passed to [ComplexHeatmap::Heatmap()]
#'
  #' @export
plot_macrosynteny_heatmap <- function(ms.x, ms.y, metric="group", clust_method = macrosynteny_cluster, scale=F, ...) {
  scaled_chrom_matrix <- macrosynteny_matrix(ms.x, ms.y, metric=metric, scale = scale)
  hm <- ComplexHeatmap::Heatmap(scaled_chrom_matrix,
                                cluster_rows = clust_method,
                                cluster_columns = clust_method,
                                ...)

  grid::grid.grabExpr(ComplexHeatmap::draw(hm))
}

#' mark genes with the top `ngroups` macrosyntenic pairs
#'
#' Mark the genes in `ms.x` with macrosyntenic group ID based on ancestral
#' linkage groups (ALG). An ALG in this case based on pairs of chromosomes
#' that have a high number of shared genes. Genes on the highest
#' `ngroups` pairs of chromosomes are marked as having an ALG, the rest
#' are given `NA`. The resulting [tibble][tibble::tibble-package] will have
#' the columns
#' \itemize{
#'   \item{Geneid} The name of the gene
#'   \item{alg} Ancestral linkage group
#'   }
#'
#' @param ms.x a msynt [tibble][tibble::tibble-package]. Genes from this
#' @param ms.y a msynt [tibble][tibble::tibble-package]
#' @param ngroups the number of groups to detect. The biggest pairs
#'
#' @export
detect_alg_genes.msynt <- function(ms.x, ms.y, ngroups, max_size=NULL) {
  joined <- dplyr::inner_join(ms.x, ms.y, by="group")
  groups <- joined %>%
    dplyr::count(chrom.x, chrom.y) %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::slice(1:ngroups) %>%
    tibble::rowid_to_column("alg")

  dplyr::inner_join(joined,groups, by=c("chrom.x", "chrom.y")) %>%
    dplyr::select(Geneid.x, Geneid.y, alg, group, chrom.x, chrom.y) %>%
    tidyr::pivot_longer(-c(chrom.x, chrom.y, alg, group), values_to="Geneid") %>%
    dplyr::select(Geneid, alg, chrom.x, chrom.y, group)
}

detect_alg_genes <- function(pos_a, pos_b, ogs, n_algs, max_group_size = 2, min_scaffold_size=50) {
  ms_ab <- create_msynt(pos_a, ogs, max_group_size = max_group_size, min_scaffold_size = min_scaffold_size)
  ms_ba <- create_msynt(pos_b, ogs, max_group_size = max_group_size, min_scaffold_size = min_scaffold_size)
  detect_alg_genes.msynt(ms_ab, ms_ba, ngroups = n_algs)
}

quantile_detect_alg_genes.msynt <- function(ms.x, ms.y, q, max_size=NULL) {
  joined <- dplyr::inner_join(ms.x, ms.y, by="group")
  groups <- joined %>%
    dplyr::count(chrom.x, chrom.y)

  threshold <- quantile(groups$n, q)

  groups <- groups %>%
    dplyr::filter(n >= threshold) %>%
    tibble::rowid_to_column("alg")

  dplyr::inner_join(joined,groups, by=c("chrom.x", "chrom.y")) %>%
    dplyr::select(Geneid.x, Geneid.y, alg, group, chrom.x, chrom.y) %>%
    tidyr::pivot_longer(-c(chrom.x, chrom.y, alg, group), values_to="Geneid") %>%
    dplyr::select(Geneid, alg, chrom.x, chrom.y, group)
}


quantile_detect_alg_genes <- function(pos_a, pos_b, ogs, q, max_group_size = 2, min_scaffold_size=50) {
  ms_ab <- create_msynt(pos_a, ogs, max_group_size = max_group_size, min_scaffold_size = min_scaffold_size)
  ms_ba <- create_msynt(pos_b, ogs, max_group_size = max_group_size, min_scaffold_size = min_scaffold_size)
  quantile_detect_alg_genes.msynt(ms_ab, ms_ba, q = q)
}

auto_detect_alg_genes.msynt <- function(ms.x, ms.y, max_group_size = 2, min_scaffold_size=50, adjust_method="none") {
  all_algs <- detect_alg_genes.msynt(ms.x, ms.y, ngroups = dplyr::n())

  within_alg <- all_algs %>% dplyr::group_by(chrom.x, chrom.y) %>% dplyr::summarize(n=dplyr::n_distinct(group), .groups="drop") %>%
    tidyr::pivot_wider(names_from=chrom.y, values_from=n) %>% tibble::column_to_rownames("chrom.x") %>% as.matrix()

  outside_alg <- purrr::cross2(seq(nrow(within_alg)), seq(ncol(within_alg))) %>%
    purrr::map_int(~sum(within_alg[-.x[[1]],-.x[[2]]], na.rm = T)) %>%
    matrix(nrow=nrow(within_alg), dimnames = dimnames(within_alg))

  inside_chr.x <- -1*sweep(within_alg, 1, rowSums(within_alg, na.rm = T))
  inside_chr.y <- -1*sweep(within_alg, 2, colSums(within_alg, na.rm = T))

  n_tests <- sum(!is.na(within_alg))

  pval <- purrr::cross_df(list(x=rownames(within_alg), y=colnames(within_alg))) %>% purrr::pmap_dbl(
    function(x,y)  {
      if(is.na(within_alg[x,y])) { 1 } else {
        matrix(c(within_alg[x,y], inside_chr.x[x,y], inside_chr.y[x,y], outside_alg[x,y]), nrow=2) %>%
          fisher.test(alternative="greater") %>% `[[`("p.value")
      }
    }
  ) %>% matrix(nrow=nrow(within_alg), dimnames = dimnames(within_alg))
  #qval <- pval * n_tests
  qval <- as.numeric(pval) %>%
    p.adjust(method = adjust_method) %>%
    matrix(nrow=nrow(pval), ncol=ncol(pval), dimnames=dimnames(pval))

  groups <- purrr::cross_df(list(chrom.x=rownames(within_alg), chrom.y=colnames(within_alg))) %>%
    dplyr::mutate(q=purrr::map2_dbl(chrom.x, chrom.y, ~qval[.x,.y])) %>%
    dplyr::arrange(q) %>%
    dplyr::filter(q < 0.05) %>%
    tibble::rowid_to_column("alg")

  joined <- dplyr::inner_join(ms.x, ms.y, by="group")

  dplyr::inner_join(joined,groups, by=c("chrom.x", "chrom.y")) %>%
    dplyr::select(Geneid.x, Geneid.y, alg, group, chrom.x, chrom.y) %>%
    tidyr::pivot_longer(-c(chrom.x, chrom.y, alg, group), values_to="Geneid") %>%
    dplyr::select(Geneid, alg, chrom.x, chrom.y, group)
}

#' compute chromosomal c-values
#'
#' Construct a matrix of c-values, which is the number of genes that are
#' conserved in each pair of chromosomes, as compared to the total number
#' of genes in each of the chromosomes. Optionally, adjust c-values to
#' correct for pairs of genomes with poorly conserved chromosomes.
#'
#' @param ms.x a msynt [tibble][tibble::tibble-package]. Genes from this
#' @param ms.y a msynt [tibble][tibble::tibble-package]
#' @param adjust one of
#'     max - divide all c-values by the maximum c-value
#'     minmax - divide c-values by the lower maximum c-value assigned to the
#'              chromosome of either species. This is useful if the links
#'              where a normal c-value would penalize chromosomal fusions
#'     ci - divide by the conservation index, as defined by the number of
#'          genes placed in ancestral linkage groups divided by the total
#'          number of orthologous genes placed on the genome.
#'     ci_mm - adjust the minmax c-values by the conservation
#'
#' @export
compute_cvalues <- function(ms.x, ms.y, adjust = c("max", "minmax", "ci", "ci_mm", "none"), alg_genes = NULL) {
  adjust <- match.arg(adjust)
  msmat <- macrosynteny_matrix(ms.x, ms.y)

  # note we subtract msmat from the denominator so as not to count each cell twice
  cmat <- msmat/(outer(rowSums(msmat), colSums(msmat), FUN = "+")-msmat)


  if(stringr::str_detect(adjust, "ci")) {
    if(is.null(alg_genes)) stop("Need alg genes from *detect_alg_genes() in order to compute ci")
    ci <- dplyr::n_distinct(alg_genes$group)/sum(msmat)
  }

  if(adjust == "max") {
    cmat <- cmat/max(cmat)
  } else if(adjust == "minmax" || adjust == "ci_mm") {
    cmat <- cmat/outer(apply(cmat, 1, max), apply(cmat, 2, max), FUN = "pmin")
    if(adjust == "ci_mm") {
      cmat <- cmat * ci
    }
  } else if(adjust == "ci") {
    cmat <- cmat/ci
  }

  cmat
}

auto_detect_alg_genes <- function(pos_a, pos_b, ogs, max_group_size = 2, min_scaffold_size=50) {
  ms_ab <- create_msynt(pos_a, ogs, max_group_size = max_group_size, min_scaffold_size = min_scaffold_size)
  ms_ba <- create_msynt(pos_b, ogs, max_group_size = max_group_size, min_scaffold_size = min_scaffold_size)
  auto_detect_alg_genes.msynt(ms_ab, ms_ba, max_group_size = max_group_size, min_scaffold_size = min_scaffold_size)
}

#' Default clustering method for macrosynteny analysis
#'
#' @param m the matrix of counts common to scaffolds or chromosomes
#' @param priority the priority of each oto be reordered when all other things are equal. passed to reorder.dendrogram()
#'
#' @export
macrosynteny_cluster <- function(m, priority, ...) {
  d <- dist(m, method = "euclidean")
  h <- hclust(d, method = "ward.D2")
  dend <- as.dendrogram(h)
  if(!is.null(priority)) {
    dend <- reorder(dend, priority)
  }
  order.dendrogram(dend)
}

#' Generate a seriation based method
#'
#' This generates a method like macrosynteny_cluster but based on a chosen method in from the seriation package.
#'
#' @param method the method to be used. see the seriation package for details
#' @param dist whether to convert the matrix to a distance before passing to seriat()
#'
#' @export
seriate_macrosynteny_cluster_method <- function(method="TSP", dist = T, ...) {
  function(m, priority) {
    m.s <- t(scale(t(m)))
    d <- if(dist) dist(m, method = "euclidean") else m
    s <- seriation::seriate(d, method = method, ...)
    seriation::get_order(s)
  }
}

# rank the genes in the order of chromosomes that have the most in common, as they would be organized on a heatmap
# order.x and order.y may be a return value of rank_genes. this will force the chromosome order of x or y
# the order previously found
#' Rank genes across the genome based on macrosynnteny and chromosome position
#'
#' Assign a position in the genome by discretizing chromosomes to one gene per
#' position from 1 to the number of genes in the genome. The genes will be ordered
#' based on both similarity in number of shared genes per chromosome pair, then by
#' position in the chromosome
#'
#' @param  ms.x a msynt [tibble][tibble::tibble-package].
#' @param  ms.y a msynt [tibble][tibble::tibble-package].
#' @param metric a column in msynt (usually either `group` or `alg`) by which the genes are ranked
#' @param clust_method a function accepting the matrix and an additional priority parameter (which may be ignored) for dendrogram-based methods
#' @param order.x pass a customized ordering of genes for the x axis
#' @param order.y pass a customized ordering of genes for the y axis
# @param scale whether to scale the initial matrix. This works better for comparing to scaffolded genomes.
#'
#' @export
rank_genes <- function(ms.x, ms.y, metric="group", order.x=F, order.y=F,
                       clust_method = macrosynteny_cluster) {
  chrmat <- macrosynteny_matrix(ms.x, ms.y, metric = metric)

  .get_order <- function(ranks, matnames, axis) {
    chrom_rank_col <- stringr::str_c("chrom_rank.", axis)
    chrom_col <- stringr::str_c("chrom.", axis)
    names <- dplyr::arrange(ranks, !!rlang::ensym(chrom_rank_col)) %>%
      dplyr::pull(!!rlang::ensym(chrom_col)) %>%
      unique() %>% rev()
    match(names, matnames)
  }
  x_order <- if(order.x!=F) .get_order(order.x, rownames(chrmat), "x") else clust_method(chrmat, rowMeans(chrmat))
  y_order <- if(order.y!=F) .get_order(order.y, colnames(chrmat), "y") else clust_method(t(chrmat), colMeans(chrmat))

  .add_ranks <- function(ms, names, chrom_order)  {
    ms %>%
      dplyr::full_join(tibble::tibble(chrom=names[rev(chrom_order)]) %>%
                       tibble::rowid_to_column("chrom_rank"), by="chrom") %>%
      dplyr::arrange(chrom_rank, position) %>%
      #tibble::rowid_to_column("rank") %>%
      dplyr::mutate(rank=seq(nrow(.))) %>%
      dplyr::select(chrom, chrom_rank, Geneid, group, rank, position)
  }
  gene_ranks.x <- .add_ranks(ms.x, rownames(chrmat), x_order)
  gene_ranks.y <- .add_ranks(ms.y, colnames(chrmat), y_order)

  dplyr::inner_join(gene_ranks.x, gene_ranks.y, by="group") %>%
    dplyr::mutate(rank.x = rank(rank.x), rank.y=rank(rank.y))
}

#' Align ranks between two ranked genomes
#'
#' @param gene_ranks.xy gene ranks returned from [gene_ranks()]
#' @param gene_ranks.zy gene ranks returned from [gene_ranks()]
#'
#' @return A list of two gene ranks [tibble][tibble::tibble-package]s with the y ranks of the y of both aligned.
#'
#' @export
join_ranks <- function(gene_ranks.xy, gene_ranks.zy) {
  aligned_ranks <-
    dplyr::rename_with(dplyr::arrange(gene_ranks.zy, rank.y), ~stringr::str_replace(.x, "x$", "z"), ends_with("x")) %>%
    dplyr::rename(group.z = "group") %>%
    dplyr::full_join(dplyr::rename(dplyr::arrange(gene_ranks.xy, rank.y), group.x="group"),
                     by=c("chrom.y", "chrom_rank.y", "Geneid.y", "position.y")) %>%
    dplyr::arrange(chrom_rank.y, position.y) %>%
    tibble::rowid_to_column(var = "rank.y") %>%
    dplyr::select(-rank.y.x, -rank.y.y)

  aligned_ranks.xy <- aligned_ranks %>%
    dplyr::select(ends_with("x")|ends_with("y"))

  aligned_ranks.zy <- aligned_ranks %>%
    dplyr::select(ends_with("z")|ends_with("y")) %>%
    dplyr::rename_with(~stringr::str_replace(.x, "z$", "x"), ends_with("z"))

  purrr::map(list(aligned_ranks.xy, aligned_ranks.zy),
             ~dplyr::rename_with(.x, ~"group", starts_with("group")) %>%
              dplyr::select(chrom.x, chrom_rank.x, Geneid.x, group, position.x, rank.x,
                            chrom.y, chrom_rank.y, Geneid.y, position.y, rank.y))
}

#' Make a macrosynteny dotplot
#'
#' @param gene_ranks the gene ranks [tibble][tibble::tibble-package] from [gene_ranks()]
#' @param alg_genes optional tibble of alg genes to mark linkage groups
#' @param alg_genes_only only plot genes with an alg group
#' @param show_chrom_names show the chromosome names on the axes. May be T, F or a vector
#'                         of 2 logicals incidating whether to show the x and y chrom names, respectively.
#' @param xtitle a title for the x axis (usually the species)
#' @param ytitle a title for the y axis
#' @param point_size the size of the points as passed to [ggplot2::geom_point]
#' @param shape the shape of points as passed to [ggplot2::geom_point]
#' @param annotate_algs show names of the ancestral linkage groups by color
#' @param nonorthologous_color genes which do not fit into an ancestral linkage group
#'
#' @export
#'
#'
# TODO the null value of alg_genes is handled twice, figure out which to keep
plot_macrosynteny_dots <- function (gene_ranks, alg_genes = NULL, alg_genes_only = F, show_chrom_names = F, xtitle = "", ytitle = "",
                                    point_size = .2, shape=16, linetype="solid", linesize=.2, linecolor="darkgrey",
                                    annotate_algs = F, show_legend = F, show_na = T, cols = NULL, nonorthologous_color="#777777FF") {
  if(!is.null(alg_genes))
    gene_ranks <- dplyr::full_join(gene_ranks, alg_genes %>% dplyr::select(Geneid, alg), by=c(Geneid.y = "Geneid")) %>% dplyr::arrange(dplyr::desc(is.na(alg)))
  else
    gene_ranks <- dplyr::mutate(gene_ranks, alg=NA)

  if(alg_genes_only)
    gene_ranks <- gene_ranks %>% dplyr::filter(!is.na(alg))

  chroms.x <- gene_ranks %>% dplyr::group_by(chrom.x) %>% dplyr::summarize(start = min(rank.x), end = max(rank.x), .groups = "drop") %>% dplyr::slice(1:(dplyr::n() -1))
  chroms.y <- gene_ranks %>% dplyr::group_by(chrom.y) %>% dplyr::summarize(start = min(rank.y), end = max(rank.y), .groups = "drop") %>% dplyr::slice(1:(dplyr::n() -1))

  if(!is.null(alg_genes))  {
    p <- ggplot2::ggplot(gene_ranks, ggplot2::aes(x=rank.x, y = rank.y, color = factor(alg), fill = factor(alg))) +
      ggplot2::geom_point(size=point_size, shape=shape, key_glyph = "rect")
  } else {
    p <- ggplot2::ggplot(gene_ranks, ggplot2::aes(x=rank.x, y = rank.y)) +
      ggplot2::geom_point(size=point_size, shape=shape, color = nonorthologous_color, key_glyph = "rect")
  }

  if(length(show_chrom_names) == 1) { show_chrom_names <- rep(show_chrom_names, 2) }
  if(show_chrom_names[1]) {
    p <- p + ggplot2::scale_x_continuous(expand = c(0,0), breaks = purrr::map2_dbl(chroms.x$start, chroms.x$end, ~(.x+.y)/2), labels = chroms.x$chrom.x)
  } else {
    p <- p + ggplot2::scale_x_continuous(expand = c(0,0))
  }
  if(show_chrom_names[2]) {
    p <- p + ggplot2::scale_y_continuous(expand = c(0,0), breaks = purrr::map2_dbl(chroms.y$start, chroms.y$end, ~(.x+.y)/2), labels = chroms.y$chrom.y)
  } else {
    p <- p + ggplot2::scale_y_continuous(expand = c(0,0))
  }

  legend_pos <- if(show_legend) "bottom" else "none"

  if(is.null(cols))
    p <- p +  ggplot2::scale_color_discrete(na.value = nonorthologous_color, guide = F, na.translate = show_na) +
    ggplot2::scale_fill_discrete(na.value = nonorthologous_color, na.translate = show_na)
  else
    p <- p +  ggplot2::scale_color_manual(values = cols, na.value = nonorthologous_color, guide = F, na.translate = show_na) +
    ggplot2::scale_fill_manual(values = cols, na.value = nonorthologous_color, na.translate  = show_na)

  p <- p +
    ggplot2::geom_hline(yintercept = chroms.y$end, linetype = linetype, lwd = linesize, color = linecolor) +
    ggplot2::geom_vline(xintercept = chroms.x$end, linetype = linetype, lwd = linesize, color = linecolor) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = xtitle, y = ytitle, color = "ALG", fill = "ALG") +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = legend_pos)

  if(stringr::str_length(xtitle) == 0)
    p <- p + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  if(stringr::str_length(ytitle) == 0)
    p <- p + ggplot2::theme(axis.title.y = ggplot2::element_blank())

  if(annotate_algs) {
    centers <- gene_ranks %>%
      dplyr::group_by(alg) %>%
      dplyr::summarize(center.x = mean(rank.x, na.rm=T), center.y = mean(rank.y, na.rm=T))
    p <- p + ggplot2::geom_text(data = centers, mapping = ggplot2::aes(x=center.x, y=center.y, label=alg),
                                inherit.aes = F)
  }

  if(show_chrom_names[2] == F)
    p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  if(show_chrom_names[1])
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = .5))
  else
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank())

  p
}

