best_hits <- function(blast6, column="bitscore", keyf=max) {
  blast6 %>% dplyr::group_by(qaccver) %>% dplyr::filter(!!rlang::ensym(column) == keyf(!!rlang::ensym(column))) %>% dplyr::ungroup()
}

reciprocal_hits <- function(blast6x, blast6y, ...) {
  dplyr::inner_join(blast6x, blast6y, by=c(qaccver="saccver", saccver="qaccver"), ...)
}

rbb_to_orthogroups <- function(reciprocals) {
  tidygraph::as_tbl_graph(reciprocals %>% dplyr::select(saccver, qaccver)) %>%
    tidygraph::activate(nodes) %>%
    tidygraph::mutate(group = tidygraph::group_components()) %>%
    tidygraph::as_tibble() %>%
    dplyr::rename(Geneid=name)
}

read_pair_hits <- function(name_a, name_b, dir=".", glue="{x}-{y}") {
  tibble::tibble(x=c(name_a, name_b), y=c(name_b, name_a)) %>%
    stringr::str_glue_data(glue) %>%
    purrr::map(~fs::path(dir, .x)) %>%
    purrr::map(~find_file_alternatives(.x, suffixes = c("", ".gz", ".bz2"))) %>%
    purrr::map(read_blast6)
}

pair_reciprocal_hits <- function(name_a, name_b, dir=".", glue="{x}-{y}", subset_a=NULL, subset_b=NULL, pep_gene_tbl_a=NULL, pep_gene_tbl_b=NULL, ...) {
  hits <- read_pair_hits(name_a, name_b, dir, glue)

  filter_subset <- function(x,y,set) dplyr::filter(x, !!rlang::ensym(y) %in% set)
  if(!is.null(subset_a))
    hits <- purrr::map2(hits, c("qaccver", "saccver"), filter_subset, set = subset_a)
  if(!is.null(subset_b))
    hits <- purrr::map2(hits, c("saccver", "qaccver"), filter_subset, set = subset_b)

  convert_pep <- function(x,y,tbl) {
    dplyr::full_join(x, dplyr::rename(tbl, !!rlang::ensym(y):= "pep"), by=y) %>%
    dplyr::select(-!!rlang::ensym(y), !!rlang::ensym(y) := "gene")
  }
  if(!is.null(pep_gene_tbl_a))
    hits <- purrr::map2(hits, c("qaccver", "saccver"), convert_pep, tbl = pep_gene_tbl_a)
  if(!is.null(pep_gene_tbl_b))
    hits <- purrr::map2(hits, c("saccver", "qaccver"), convert_pep, tbl = pep_gene_tbl_b)

  hits %>%
    purrr::map(best_hits, ...) %>%
    purrr::lift(reciprocal_hits)()
}

pair_orthogroups <- function(...) {
  pair_reciprocal_hits(...) %>%
    rbb_to_orthogroups()
}
