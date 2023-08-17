resnv56 <- genomescope("data/kmer_histograms/nvmalem56.histo", k = 56, readlength = 75)
ressc18 <- genomescope("data/kmer_histograms/scm18.histo", k = 18, readlength = 100)

genomescope_plot56 <-
  patchwork::wrap_plots(plot_genomescope(resnv56, title = "N. vectensis") +
                          ggplot2::theme(plot.margin = ggplot2::unit(c(5.5,6,5.5,5.5), "points")),
                        plot_genomescope(ressc18, title = "S. callimorphus"),
                        guides = "collect") +
  patchwork::plot_annotation(tag_levels="a") &
  ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold"))

save_fig(genomescope_plot56, prefix = "EDF2a_size_estimates", width = 7, height = 3)

