buscos <- parse_buscos(nv1scaffolds="data/busco_summaries/short_summary_nemVec1.txt",
                       nv2contigs="data/busco_summaries/short_summary_simr_r0p3_c0p045.txt",
                       nv2purged="data/busco_summaries/short_summary_simr_r0p3_c0p045.contigs.ph.txt",
                       nv2dovetail="data/busco_summaries/short_summary_nematostella_vectensis_07Nov2019_hUzXf.txt",
                       nv2chroms="data/busco_summaries/short_summary_nv_bwa_lowmasked_gapped_chroms.txt",
                       nv2chromsonly="data/busco_summaries/short_summary_nv_dovetail_4_gapped_chroms_only.txt",
                       sc1contigs="data/busco_summaries/short_summary_scol3PB_r0p3_c0p045.contigs.fasta.txt",
                       sc1purged="data/busco_summaries/short_summary_sc_markII.2purge_haplotigs_contigs.txt",
                       sc1chroms="data/busco_summaries/short_summary_sc_markII.4purged_remove_tigs_misplaced_genes.markII.review_gapped_chroms.final.txt",
                       sc1chromsonly="data/busco_summaries/short_summary_sc_markII.4purged_remove_tigs_misplaced_genes.markII.review_gapped_chroms.final_chromsonly.txt") %>%
  dplyr::mutate(names = factor(names, levels=c("nv1scaffolds", "nv2contigs", "nv2purged", "nv2dovetail", "nv2chroms", "nv2chromsonly", "sc1contigs", "sc1purged", "sc1chroms", "sc1chromsonly")) %>%
                  forcats::fct_rev())
busco_plot <- plot_buscos(buscos)
save_fig(busco_plot, "EDF2i_busco", width=6, height=4)
