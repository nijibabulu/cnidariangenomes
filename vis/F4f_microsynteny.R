ms <- readr::read_table("data/tables/microsynteny.tsv", col_types = "cicf") |>
  dplyr::mutate(Clade = factor(Clade, levels = c("Deuterostome", "Spiralia", "Cnidaria")))

pos <- ggplot2::position_dodge(width = 0.6)
p <-
  ggplot2::ggplot(ms, ggplot2::aes(x=Clade, fill=Clade, y=Counts)) +
  ggplot2::theme_light() +
  ggplot2::facet_wrap(~NumGenes, strip.position = "bottom") +
  ggplot2::geom_dotplot(binaxis="y", stackdir="centerwhole", position = "dodge", dotsize=0.5) +
  ggplot2::stat_summary(fun.min = median, position = pos, geom = "errorbar",
                        linewidth=1.5, ggplot2::aes(col=Clade)) +
  ggplot2::scale_y_continuous(trans=scales::pseudo_log_trans(),
                              breaks = c(0, 10, 100, 1000)) +
  ggsci::scale_fill_d3() +
  ggsci::scale_color_d3() +
  ggplot2::labs(x="Minimum syntenic block length", y="Microsyntenic blocks", fill="", col="") +
  ggplot2::theme(legend.position="top",
                 axis.text.x = ggplot2::element_blank())

save_fig(p, "F4f_microsynteny", width = 4, height = 6.4)
