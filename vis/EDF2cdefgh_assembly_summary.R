
# Note to the reader -- Nemve1.ctg.fasta is generated using the AGOUTI shred
# tool with default settings. For more information about the tool see the
# https://github.com/svm-zhang/AGOUTI
# The input is the original Nematostella vectensis assembly downloadable from
# https://mycocosm.jgi.doe.gov/Nemve1/Nemve1.home.html

nv_seq_lengths <- fasta_summaries(nv1contigs="data/assemblies/nv_seqs/Nemve1/Nemve1.ctg.fasta",
                                  nv1scaffolds="data/assemblies/nv_seqs/Nemve1/Nemve1.fasta",
                                  nv2contigs="data/assemblies/nv_seqs/simr_r0p3_c0p045.contigs.fasta",
                                  nv2purged="data/assemblies/nv_seqs/simr_r0p3_c0p045.contigs.ph.fasta",
                                  nv2dovetail="data/assemblies/nv_seqs/nematostella_vectensis_07Nov2019_hUzXf.fasta") %>%
                        dplyr::bind_rows(nv1_contigs %>% dplyr::mutate(assembly="nv1contigs"))
nv_seq_lengths <- nv_seq_lengths %>% dplyr::mutate(assembly=factor(assembly, levels = c("nv1contigs", "nv1scaffolds", "nv2contigs",  "nv2purged", "nv2dovetail")))
nv_assembly_summaries <- summarize_assemblies(nv_seq_lengths)
nv_summary <-
  patchwork::wrap_plots(assembly_size_plot(nv_seq_lengths, nv_assembly_summaries, label = F, log_x_axis = T),
                        assembly_nx_plot(nv_seq_lengths, nv_assembly_summaries, label = F),
                        seq_length_distrib_plot(nv_seq_lengths, nv_assembly_summaries, label = F, bins = 10),
                        nrow = 1, guides = "collect") & ggplot2::theme(legend.title = ggplot2::element_blank())


sc_seq_lengths <- fasta_summaries(sc1contigs="data/assemblies/sc_seqs/scol3PB_r0p3_c0p045.contigs.fasta",
                                  sc1purged="data/assemblies/sc_seqs/sc_purged_contigs.fasta")
sc_assembly_summaries <- summarize_assemblies(sc_seq_lengths)
sc_summary <-
  patchwork::wrap_plots(assembly_size_plot(sc_seq_lengths, sc_assembly_summaries, label = F, upper_expansion = 0),
                        assembly_nx_plot(sc_seq_lengths, sc_assembly_summaries, label = F),
                        seq_length_distrib_plot(sc_seq_lengths, sc_assembly_summaries, label = F, bins = 10),
                        nrow = 1, guides = "collect") & ggplot2::theme(legend.title = ggplot2::element_blank())


combined_seq_lengths <- dplyr::bind_rows(nv_seq_lengths, sc_seq_lengths)
combined_assembly_summaries <- dplyr::bind_rows(nv_assembly_summaries, sc_assembly_summaries)
combined_seq_lengths <- combined_seq_lengths %>% dplyr::mutate(assembly=factor(assembly, levels = c("nv1contigs", "nv1scaffolds", "nv2contigs",  "nv2purged", "nv2dovetail", "sc1contigs", "sc1purged")))
combined_summary <-
    patchwork::wrap_plots(assembly_size_plot(combined_seq_lengths, combined_assembly_summaries, label = F, log_x_axis = T),
                        assembly_nx_plot(combined_seq_lengths, combined_assembly_summaries, label = F),
                        seq_length_distrib_plot(combined_seq_lengths, combined_assembly_summaries, label = F, bins = 10),
                        nrow = 1, guides = "collect") & ggplot2::theme(legend.title = ggplot2::element_blank())


finalized_summaries <- fasta_summaries(sc2chromosomes="data/reviewed_assemblies/sc/sc.fasta",
                                       nv2chromosomes="data/reviewed_assemblies/nv/nv.fasta") %>%
  summarize_assemblies() %>%
  dplyr::mutate(n50=NA)



reapr_stats <- purrr::map_dfr(c("nemVec1.fmt/unknown_out", "nemVec1.fmt/male_out", "nemVec1.fmt/female_out",
                 "nv_dovetail_scaffolds/unknown_out", "nv_dovetail_scaffolds/male_out", "nv_dovetail_scaffolds/female_out"),
               ~readr::read_tsv(stringr::str_glue("data/reapr/{.x}/05.summary.report.tsv"),
                                col_types = readr::cols())) %>%
  dplyr::mutate(assembly=c(rep("nv1",3), rep("nv2",3)),
                sample=rep(c("unknown","male","female"),2) %>% as.character()) %>%
  dplyr::select(assembly, sample, bases, N50,gaps, gaps_bases, bases_br, sequences_br, N50_br, gaps_br, error_free, FCD, FCD_gap, frag_cov) %>%
  dplyr::mutate(error_free_frac = error_free/bases, broken_frac=N50_br/N50)

mod_reapr_stats <- reapr_stats %>%
  dplyr::filter(sample != "unknown") %>%
  dplyr::mutate(sample = factor(sample,
                                levels=c("male", "female"),
                                labels=c("Animal 1 Corrected", "Animal 2 Corrected")))

broken_n50_plot <-
  mod_reapr_stats %>% dplyr::select(assembly, sample,  N50=N50_br) %>%
  dplyr::bind_rows(reapr_stats %>% dplyr::distinct(assembly, N50) %>% dplyr::select(assembly, N50) %>% dplyr::mutate(sample = "Original")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels=c("Original", "Animal 1 Corrected", "Animal 2 Corrected")))  %>%
  ggplot2::ggplot(ggplot2::aes(x=assembly, y=N50, fill=sample)) +
  ggplot2::geom_col(position="dodge") +
  ggthemes::scale_fill_ptol() +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.5, 0.5))) +
  ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", scales::math_format(10^.x)),
                         expand = ggplot2::expansion(mult = c(0, 0.01))) +
  ggplot2::theme_bw()  +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol=1)) +
  ggplot2::theme(legend.title = ggplot2::element_blank(),
                 legend.position = "bottom",
                 axis.title.x = ggplot2::element_blank(),
                 legend.key.size = ggplot2::unit(0.2, "cm")
)


error_free_plot <- mod_reapr_stats %>% dplyr::select(assembly, sample, error_free_frac) %>%
  ggplot2::ggplot(ggplot2::aes(x = assembly, y = error_free_frac, fill = sample)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::scale_fill_manual(values = ggthemes::ptol_pal()(3)[2:3]) +
  ggplot2::scale_y_continuous(limits = c(0,1)) +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.5, 0.5))) +
  ggplot2::labs(y = "Error Free Fraction") +
  ggplot2::theme_bw() +
  ggplot2::guides(fill = "none") +
  ggplot2::theme(legend.title = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_blank()  )

reapr_plot <- patchwork::wrap_plots(broken_n50_plot, error_free_plot, nrow = 1, guides = "collect") &
  patchwork::plot_annotation(title = "c", subtitle = "REAPR Assessment") &
  ggplot2::theme(legend.position = "bottom")


assembly_summary_table <- dplyr::bind_rows(nv_assembly_summaries, sc_assembly_summaries, finalized_summaries) %>%
  dplyr::select(Assembly = assembly,  Length = total_length, `No. of Sequences`=nseqs, `Median Length` = median_length, N50=n50) %>%
  tidyr::separate(Assembly, c("Assembly", "Level"), sep=3) %>%
  dplyr::mutate(Level = factor(Level, levels=c("contigs", "scaffolds", "purged", "dovetail", "chromosomes"))) %>%
  dplyr::arrange(Assembly, Level) %>%
  dplyr::mutate(Level=as.character(Level))

summary_table_flexraster <-
  assembly_summary_table %>%
  flextable::flextable() %>%
  flextable::merge_v(j = "Assembly") %>%
  flextable::bold(part = "header") %>%
  flextable::colformat_int( big.mark = ",") %>%
  flextable::theme_vanilla() %>%
  flextable::as_raster()

summary_table_plot <-
  summary_table_flexraster %>%
  grid::rasterGrob(width = ggplot2::unit(4, "in"), height = ggplot2::unit(2, "in"))

assembly_summary <- patchwork::wrap_plots(combined_summary, reapr_plot, summary_table_plot,
                      design="AAA
                              BCC" ) +
  patchwork::plot_annotation(tag_level = list(letters[3:8])) &
  ggplot2::theme(plot.title = ggplot2::element_text(size = ggplot2::rel(.9)),
                 axis.title = ggplot2::element_text(size = ggplot2::rel(.8)),
                 axis.text = ggplot2::element_text(size = ggplot2::rel(.7)),
                 legend.text = ggplot2::element_text(size = ggplot2::rel(.9)),
                 plot.tag = ggplot2::element_text(face = "bold"))

save_fig(assembly_summary, "EDF2cdefgh_assembly_summary", width = 7, height = 5)
