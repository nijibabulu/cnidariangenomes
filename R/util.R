save_fig <- function(plot, prefix, dir = "figures/", types = c("pdf", "png"), ...) {
  if(!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }
  purrr::walk(types, ~ggplot2::ggsave(plot = plot,
                                      file = fs::path_join(c(dir, stringr::str_glue("{prefix}.{.x}"))),
                                      ...),
              ...)
}

find_file_alternatives <- function(file, suffixes = "") {
  files <- stringr::str_glue("{file}{suffixes}")
  found_file <- purrr::detect(files, fs::file_exists)
  if(is.null(found_file)) {
    stop(stringr::str_glue("No file found among {stringr::str_c(files, collapse = ', ')}"))
  }
  found_file
}

#expects chrom, start, end, outputs Chrom, Start, End, Value
#' Calculate the number of features binned by a given size
#'
#' @param features a [tibble][tibble::tibble-package] containing the columns chrom, start and end
#' @param genome a [tibble][tibble::tibble-package] containing the columns name and length, describing each chromosome
#' @param binsize the size of the bins to calculate to
#'
#' @return a [tibble][tibble::tibble-package] containing the columns Chr, Start, End, describing the bin,
#'   and Value, number of features found in that bin
#'
#' @export
calculate_feature_density <- function(features, genome, binsize = 100000) {
  empty <- purrr::map2_dfr(genome$name, genome$length,
                           ~tibble::tibble(Chr=.x,
                                           Start=seq(0, .y, binsize),
                                           End=pmin(Start+binsize, .y),
                                           Value=0, Length=0)) %>%
    dplyr::filter(Chr %in% unique(features$chrom))
  features %>%
    dplyr::mutate(length = end-start) %>%
    dplyr::left_join(genome %>% dplyr::select(chrom="name", chr_length="length"), by = "chrom") %>%
    dplyr::mutate(bin = cut(start, breaks = c(seq(0, max(chr_length), binsize), Inf))) %>%
    dplyr::mutate(bin = stringr::str_sub(bin, 2, -2)) %>%
    tidyr::separate(bin, c("bin_start", "bin_end"), sep = ",", convert = TRUE) %>%
    dplyr::mutate(bin_end = pmin(bin_end, chr_length)) %>%
    dplyr::group_by(chrom, bin_start, bin_end) %>%
    dplyr::summarize(n=dplyr::n(), l=sum(length)) %>%
    dplyr::ungroup() %>%
    dplyr::select(Chr = "chrom", Start = "bin_start", End = "bin_end", Value = "n", Length = "l") %>%
    dplyr::bind_rows(dplyr::anti_join(empty, ., by=c("Chr", "Start", "End"))) %>%
    dplyr::mutate(ScaledValue = scale(Value/(End-Start))[,1])
}

#expects chrom, start, end, outputs Chrom, Start, End, Value
#' Calculate the number of features binned by a given size
#'
#' @param features a [tibble][tibble::tibble-package] containing the columns chrom, start and end
#' @param genome a [tibble][tibble::tibble-package] containing the columns name and length, describing each chromosome
#' @param binsize the size of the bins to calculate to
#'
#' @return a [tibble][tibble::tibble-package] containing the columns Chr, Start, End, describing the bin,
#'   Value and ScaledValue, the total bases in that bin covered by one or more features.
#'
#' @export
calculate_feature_genomic_coverage <- function(features, genome, binsize = 100000) {
  GenomicRanges::GRanges(features$chrom, IRanges::IRanges(features$start, features$end)) %>%
    GenomicRanges::reduce() %>%
    tibble::as_tibble() %>%
    dplyr::rename(chrom="seqnames") %>%
    dplyr::left_join(genome %>% dplyr::select(chrom="name", chr_length="length"), by = "chrom") %>%
    dplyr::mutate(bin = cut(start, breaks = c(seq(0, max(chr_length)+binsize, binsize)))) %>%
    dplyr::mutate(bin = stringr::str_sub(bin, 2, -2)) %>%
    tidyr::separate(bin, c("bin_start", "bin_end"), sep = ",", convert = TRUE) %>%
    dplyr::mutate(bin_end = pmin(chr_length, bin_end)) %>%
    dplyr::mutate(start_in_bin = pmax(start, bin_start),
                  end_in_bin = pmin(end,bin_end),
                  length = end_in_bin-start_in_bin) %>%
    dplyr::group_by(chrom, bin_start, bin_end) %>%
    dplyr::summarize(Value = sum(length)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ScaledValue = scale(Value/(bin_end-bin_start))[,1]) %>%
    dplyr::rename(Chr = "chrom", Start = "bin_start", End = "bin_end")
}


computeFPKM <- function(mat, meta) {
  lengths_m <- meta %>% dplyr::select(Length) %>% purrr::flatten_dbl()
  10**9*t(t(mat/lengths_m)/colSums(mat))
}

setsummary <- function(x,y) {
  xnm <- deparse(substitute(x))
  ynm <- deparse(substitute(y))
  xynm <- stringr::str_c(xnm, ynm, sep = "_")
  res <- list(setdiff(x,y), intersect(x,y), setdiff(y,x))
  rlang::set_names(res, c(xnm, xynm, ynm))
}

#' Return the next index of `x` which satisfies a predicate `p()`
#'
#' The advantage of using this function over `which()` is that it does not
#' calculate for all values in the vector `x`, offering a very big speedup.
#'
#' @param x a vector
#' @param p a predicate for which the next index should be true
#' @param cur the starting index
#'
#' @return
#'   the first index `i` after `cur` for which `p(x[i])` is true
#' @export
#'
#' @examples
#'   which_next(seq(10), function(x) x > 5)
#'   which_next(seq(10), function(x) x > 0)
which_next <- function(x, p, cur = 1) {
  i <- cur
  while(i < length(x) && !p(x[i])) {
    i <- i + 1
  }
  i
}
