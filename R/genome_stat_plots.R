

fasta_summaries <- function(..., genome_size = NULL) {
  fastas <- list(...)
  purrr::map2_dfr(fastas, names(fastas), fasta_summary, genome_size = genome_size)
}

read_fai <- function(fai) {
  readr::read_tsv(fai, col_types = "ciiii",
                  col_names=c("name", "length", "offset", "linebases", "linewidth"))
}

compute_nZ <- function(z, cumlength, length) {
  length[
    purrr::accumulate(z, ~which_next(cumlength, function(x) x >= .y, .x), .init = 1)[-1]
  ]
}

fasta_summary <- function(fasta_file, genome_size=NULL, assembly=NULL, compute_stats=TRUE) {
  index = stringr::str_glue("{fasta_file}.fai")
  already_indexed = fs::file_exists(index)
  if(!already_indexed) {
    Rsamtools::indexFa(fasta_file)
  }
  faidx <- read_fai(index)
  if(!already_indexed) {
    fs::file_delete(index)
  }
  genome_size <- if(is.null(genome_size)) { 1 } else { genome_size }
  if(is.null(assembly)) {
    assembly <-basename(fasta_file)
  }

  nx <- function(data, x, sz) {
    purrr::map_int(x, ~data[data$cumlength >= .x*sz,]$length[1])
  }
  ngx <- function(data, x) {
    purrr::map_int(x, ~data[data$cumlength >= .x*genome_size,]$length[1])
  }

  assy_size <- sum(dplyr::pull(faidx, length))
  summary <-
    faidx %>%
    dplyr::select(name, length) %>%
    dplyr::arrange(desc(length)) %>%
    tibble::rowid_to_column("seqnum") %>%
    dplyr::mutate(assembly = assembly)
  if(compute_stats)  {
    summary <-
      summary %>%
      dplyr::mutate(
        x=seqnum/dplyr::n(),
        cumlength=as.integer(cumsum(length)),
        nx = compute_nZ(x*assy_size, cumlength, length),
        ngx = compute_nZ(x*genome_size, cumlength, length)
      )
  }
  summary
}

summarize_assemblies <- function(seq_lengths) {
 seq_lengths %>% dplyr::group_by(assembly) %>%
    dplyr::summarize(
              median_length=median(length),
              total_length=sum(length),
              nseqs=dplyr::n(),
              n50=nx[abs(x-.5) == min(abs(x-.5))][1],
              ng50=ngx[abs(x-.5) == min(abs(x-.5))][1]) %>%
    dplyr::mutate(size_label=as.character(stringr::str_glue("{assembly}\nn={nseqs} L={round(total_length/1e6)}Mb")),
           n50_label=stringr::str_glue("{assembly} N50={round(n50/1e6, digits=1)}Mb"),
           ng50_label=stringr::str_glue("{assembly} NG50={round(ng50/1e6, digits=1)}Mb"))
}

assembly_size_plot <- function(seq_lengths, assembly_summaries, log_x_axis = F, upper_expansion = 0.2, label = T) {
  p <- ggplot2::ggplot(seq_lengths, ggplot2::aes(x=seqnum, y=cumlength, color=assembly)) +
    ggplot2::geom_line()

  if(label)  {
    p <- p + ggplot2::geom_point(data=assembly_summaries, ggplot2::aes(x=nseqs,y=total_length), color="black") +
    ggrepel::geom_text_repel(data=assembly_summaries,
                    ggplot2::aes(x=nseqs,y=total_length,label=size_label),
                    color="black",
                    force=100)
  }

  if(log_x_axis)
    p <- p + ggplot2::scale_x_log10()

  p <- p +
    ggthemes::scale_color_ptol() +
    ggplot2::labs(x="Number of Sequences", y="Cumulative Length") +
    ggplot2::ggtitle("Assembly Size") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult=c(0,upper_expansion))) +
    #cowplot::theme_cowplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(
                   plot.title = ggplot2::element_text(hjust=0.5, face = "bold"),
                   panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90))
  if(label)
    p <- p + ggplot2::theme(legend.position = "none")

  p

}


assembly_nx_plot <- function(seq_lengths, assembly_summaries, genome_size = NULL, label = T) {
    if(is.null(genome_size)) {
      p <- ggplot2::ggplot(seq_lengths, ggplot2::aes(x=x, y=nx, color=assembly)) +
          ggplot2::geom_line() +
          ggplot2::labs(x="x", y="Nx")
      if(label) {
        p <- p +
          ggplot2::geom_point(data=assembly_summaries, ggplot2::aes(x=.5,y=n50), color="black") +
          ggrepel::geom_text_repel(data=assembly_summaries,
                          ggplot2::aes(x=.5,y=n50,label=n50_label),
                          color="black")
      }
    } else {
        p <- p + ggplot2::ggplot(seq_lengths, ggplot2::aes(x=x, y=ngx, color=assembly)) +
          ggplot2::geom_line() +
          ggplot2::labs(x="x", y="NGx") +
        if(label) {
          p <- p + ggplot2::geom_point(data=seq_length_summaries, ggplot2::aes(x=.5,y=ng50), color="black") +
            ggrepel::geom_text_repel(data=seq_length_summaries,
                                     ggplot2::aes(x=.5,y=ng50,label=ng50_label),
                                     color="black")
        }
    }

    p <- p + ggthemes::scale_color_ptol() +
          ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                 labels = scales::trans_format("log10", scales::math_format(10^.x))) +
          ggplot2::ggtitle("Contiguity") +
          #cowplot::theme_cowplot() +
          ggplot2::theme_bw() +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face = "bold"),
                         panel.grid = ggplot2::element_blank(),
                         axis.text.x = ggplot2::element_text(angle = 90))
  if(label)
    p <- p + ggplot2::theme(legend.position = "none")

  p
}

seq_length_distrib_plot <- function(seq_lengths, assembly_summaries, title = "Length Distribution", bins = 50, label = T) {
  p <- ggplot2::ggplot(seq_lengths, ggplot2::aes(length, color=assembly, label=assembly)) +
    ggplot2::geom_freqpoly(bins=bins) +
    ggplot2::scale_x_log10() +
    ggthemes::scale_color_ptol() +
    #cowplot::theme_cowplot() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(title) +
    ggplot2::theme(legend.position="none",
                   plot.title = ggplot2::element_text(hjust=0.5, face = "bold"),
                   panel.grid = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle = 90))

  pb <- ggplot2::ggplot_build(p)

  length_summary <- seq_lengths %>% dplyr::group_by(assembly) %>%
    dplyr::summarize(med = median(length)) %>%
    dplyr::full_join(pb$data[[1]], by=c("assembly"="label")) %>%
    dplyr::mutate(diff = abs(10**x-med)) %>%
    dplyr::group_by(assembly) %>%
    dplyr::slice(which.min(diff)) %>%
    dplyr::mutate(plot_label = stringr::str_glue("{assembly}\nmedian = {med}"))

  if(label) {
    p <- p  +
      ggrepel::geom_text_repel(data=length_summary, ggplot2::aes(x=10**x, y=y, label=plot_label), color='black') +
      ggplot2::geom_point(data=length_summary, ggplot2::aes(x=10**x, y=y), inherit.aes = F)
  }

  p <- p + ggplot2::scale_y_log10()
  if(label)
    p <- p + ggplot2::theme(legend.position = "none")

  p
}


