nv_trans_table <- readr::read_tsv("data/assemblies/nv/trans_table.txt",
                                  col_names = c("old_name", "new_name"), col_types = "cc") %>%
  dplyr::mutate(old_name=stringr::str_remove(old_name, "Chromosome"))
sc_trans_table <- readr::read_tsv("data/assemblies/sc/trans_table.txt",
                                  col_names = c("old_name", "new_name"), col_types = "cc") %>%
  dplyr::mutate(old_name=stringr::str_remove(old_name, "Chromosome"))


nv_asm <- read_assembly("data/assemblies/nv", trans_table = nv_trans_table)
sc_asm <- read_assembly("data/assemblies/sc", trans_table = sc_trans_table)

nv_p <- plot_global_contact_map(nv_asm,  title = "N. vectensis", limit_quantile = .996,
                                guide_height = 0.5, triangular = T)
sc_p <- plot_global_contact_map(sc_asm,  title = "S. callimorphus", limit_quantile = .995,
                                triangular = T, guide_breaks = c(0,1,2), guide_height = 0.5)

p <- patchwork::wrap_plots(nv_p + ggplot2::labs(tag = "d"),
                           sc_p + ggplot2::labs(tag = "e"), nrow = 1)
save_fig(p, "F1de_contact_maps", width = 8, height = 4)
