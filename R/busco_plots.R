


parse_busco <- function(file) {
  readr::read_lines(file) %>%
    stringr::str_remove("^[[:space:]]") %>%
    purrr::keep(stringr::str_detect(.,"C:")) %>%
    tibble::tibble(Summary=., file=file) %>%
    tidyr::extract(Summary,
                   regex="C:([0-9.]+)%.S:([0-9.]+)%,D:([0-9.]+)%.,F:([0-9.]+)%,M:([0-9.]+)%",
                   into=c("Complete", "Single", "Duplicated", "Fragmented" , "Missing"),
                   remove=FALSE,
                   convert=TRUE) %>%
    dplyr::mutate(score = 3*Single + 2*Duplicated + Fragmented - Missing)
}

parse_buscos <- function(...) {
  files <- list(...)
  names <- names(files)
  tibble::tibble(file=as.character(files),
         names = if(is.null(names)) fs::path_file(file) else names) %>%
  dplyr::full_join(purrr::map_dfr(.$file, parse_busco), by="file")
}

plot_buscos <- function(busco_summaries, size_ratio=1, font_family="sans") {
  busco_colors = c("#56B4E9", "#3492C7", "#F0E442", "#F04442")

  busco_theme <- ggplot2::theme(plot.title = ggplot2::element_text(family=font_family, colour = "black", size = ggplot2::rel(2.2)*size_ratio, face = "bold")) +
    ggplot2::theme(legend.position="top",legend.title = ggplot2::element_blank()) +
    ggplot2::theme(legend.text = ggplot2::element_text(family=font_family, size = ggplot2::rel(1.2)*size_ratio)) +
    ggplot2::theme(panel.background = ggplot2::element_rect(color="#FFFFFF", fill="white")) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank()) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(family=font_family, colour = "black", size = ggplot2::rel(1.66)*size_ratio)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(family=font_family, colour = "black", size = ggplot2::rel(1.66)*size_ratio)) +
    ggplot2::theme(axis.line = ggplot2::element_line(size=1*size_ratio, colour = "black")) +
    ggplot2::theme(axis.ticks.length = ggplot2::unit(.85, "cm")) +
    ggplot2::theme(axis.ticks.y = ggplot2::element_line(colour="white", size = 0)) +
    ggplot2::theme(axis.ticks.x = ggplot2::element_line(colour="#222222")) +
    ggplot2::theme(axis.ticks.length = ggplot2::unit(0.4, "cm")) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(family=font_family, size=ggplot2::rel(1.2)*size_ratio))

  busco_summaries %>%
    tidyr::pivot_longer( c(Missing, Fragmented, Duplicated, Single), values_to = "Percent") %>%
    dplyr::mutate(name=factor(name, levels=c("Single", "Duplicated", "Fragmented", "Missing"))) %>%
    ggplot2::ggplot() +
    ggplot2::geom_col(ggplot2::aes(x=names, y=Percent, fill=name)) +
    ggplot2::coord_flip() +
    ggplot2::theme_gray(base_size = 8) +
    ggplot2::scale_y_continuous(labels = c("0","20","40","60","80","100"), breaks = c(0,20,40,60,80,100)) +
    ggplot2::scale_fill_manual(values = busco_colors,
                      labels =c(" Complete (C) and single-copy (S)  ",
                                " Complete (C) and duplicated (D)",
                                " Fragmented (F)  ",
                                " Missing (M)"))  +
    ggplot2::geom_text(ggplot2::aes(x = names, y = 3, label=Summary), hjust = 0) +
    busco_theme +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
        legend.position = "top",
        legend.title = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow=2,byrow=T))



}

if(0) {
g_summaries <- tibble(file=c(dir_ls("assembly_summaries")) %>% as.character()) %>%
  mutate(names = path_file(file)) %>%
  full_join(map_dfr(.$file, parse_busco), by="file")

p_summaries <- tibble(file=dir_ls("polish_summaries") %>% as.character()) %>%
  mutate(names = path_file(file)) %>%
  full_join(map_dfr(.$file, parse_busco), by="file") %>%
  mutate(polish=case_when(
                    str_detect(names, "female_uknown") ~ "Female+Unknown",
                    str_detect(names, "female_") ~ "Female",
                    str_detect(names, "male_") ~ "Male"),
         genome=if_else(str_detect(names, "dovetailreview"), "Dovetail", "NoDovetail")
  )

sc_sizes <- tibble(file=dir_ls("sc_sizes/") %>% as.character()) %>%
  mutate(names = path_file(file),
         size= map_int(file, compose(as.integer, read_lines)))
sc_summaries <- tibble(file=dir_ls("sc_busco_summaries", glob = "*") %>% as.character()) %>%
  mutate(names= path_file(file)) %>%
  full_join(map_dfr(.$file, parse_busco), by="file") %>%
  mutate(base =str_remove(.$names, ".filtered.*")) %>%
  filter(str_detect(names, "run_joined$", negate=T)) %>%
  filter(str_detect(names, "run_joined[^_]", negate=T)) %>%
  filter(names != "run_joined_opt.unfilt") %>%
  mutate(filter=case_when(
           str_detect(names, "filtered$") ~ "StandardFilter",
           str_detect(names, "filtered_rm$") ~ "AggressiveFilter",
           TRUE ~ "NoFilter"),
         prot=case_when(
           str_detect(base, "_$") ~ "NoProt",
           str_detect(base, "ane$") ~ "AnemoneProts",
           str_detect(base, "cni$") ~ "CnidarianProts",
           str_detect(base, "edw$") ~ "EdwardsiidProts",
           str_detect(base, "joined") ~ "Joined",
           TRUE ~ "ERROR"),
         mask=case_when(
           str_detect(base, "run_rb") ~ "StandardMask",
           str_detect(base, "run_rm") ~ "AggressiveMask",
           str_detect(base, "joined") ~ "JoinedMask",
           TRUE ~ "ERROR"),
         ) %>%
  unite(GeneSet, prot, mask, filter, sep="", remove = F) %>%
  left_join(sc_sizes %>% select(-file), by="names")




  sc_summaries %>% pivot_longer(c(Missing, Fragmented, Duplicated, Single), values_to = "Percent") %>%
    mutate(name=factor(name, levels=c("Single", "Duplicated", "Fragmented", "Missing"))) %>%
    ggplot() +
    geom_col(aes(x=GeneSet, y=Percent, fill=name)) +
    coord_flip() +
    theme_gray(base_size = 8) +
    scale_y_continuous(labels = c("0","20","40","60","80","100"), breaks = c(0,20,40,60,80,100)) +
    scale_fill_manual(values = busco_colors,
                      labels =c(" Complete (C) and single-copy (S)  ",
                                " Complete (C) and duplicated (D)",
                                " Fragmented (F)  ",
                                " Missing (M)"))  +
    geom_text(aes(x = GeneSet, y = 3, label=Summary), hjust = 0) +
    busco_theme +
    facet_grid(prot~filter, scales="free")


}
points_plot <- function(df, labeled=c(""), ...) {
 df %>%
  mutate(label=if_else(GeneSet %in% labeled, str_remove(Summary, ",n:978"), "")) %>%
  ggplot() +
  geom_point(aes(x=size, y=log(score), ...), size=2) +
  geom_text_repel(aes(x=size, y=log(score), label=label), color = "black") +
  theme_bw() +
  theme(legend.title = ggplot2::element_blank())
}

if(0) {
labeled <- c("AnemoneProtsStandardMaskStandardFilter",
             "AnemoneProtsAggressiveMaskStandardFilter",
             "CnidarianProtsStandardMaskAggressiveFilter",
             "JoinedJoinedMaskNoFilter")

sc_summaries %>% points_plot(color=filter)
ggsave(file = "01ScBuscoOverview.pdf", height=4, width=6)

sc_summaries %>% points_plot(
  color=filter, labeled=arrange(., score) %>% slice(1,n()) %>% pull(GeneSet))
ggsave(file = "011ScBuscoOverview.pdf", height=4, width=6)

sc_summaries %>% filter(filter == "NoFilter", mask %in% c("StandardMask", "JoinedMask")) %>%
  points_plot(shape=prot)
ggsave(file = "02ScBuscoProtFocus.pdf", height=4, width=6)

sc_summaries %>% filter(filter == "NoFilter", mask %in% c("StandardMask", "JoinedMask")) %>%
  points_plot(shape=prot,
              labeled=arrange(., score) %>% slice(1,n()) %>% pull(GeneSet))
ggsave(file = "03ScBuscoProtFocusLabel.pdf", height=4, width=6)

sc_summaries %>% filter(filter != "NoFilter") %>% points_plot(color=mask, shape=prot, labeled = labeled)
ggsave(file = "04ScBuscoFilterMask.pdf", height=4, width=6)

sc_summaries %>% filter(filter == "NoFilter", mask %in% c("StandardMask", "JoinedMask")) %>%
  points_plot(shape=prot)

sc_summaries %>% points_plot(color=mask, shape=prot)
sc_summaries %>% points_plot(color=mask, shape=prot,
                             labeled=c("JoinedJoinedMaskNoFilter",
                                       "NoProtStandardMaskAggressiveFilter"))
sc_summaries %>% points_plot(color=filter)

sc_summaries %>%
  filter(filter == "NoFilter") %>%
  points_plot(color = "mask")

sc_summaries %>%
  points_plot(labeled = "JoinedJoinedMaskNoFilter")
  ggplot() +
  geom_point(aes(x=size, y=log(score), color=filter, shape = prot), size=2) +
  geom_text_repel(aes(x=size, y=log(score), label=label), color = "black") +
  theme_bw() +
  theme(legend.title = element_blank())

ggsave(file = "ScolanthusBusco.pdf", width=6, height=4)

sc_summaries %>%
  mutate(label=if_else(GeneSet %in% labeled, str_c(GeneSet, str_remove(Summary, ",n:978"), sep="\n"), "")) %>%
  ggplot() +
  geom_point(aes(x=size, y=log(score), color=filter, shape = prot), size=2) +
  theme_bw() +
  theme(legend.title = element_blank())

ggsave(file = "ScolanthusBuscoNoLabel.pdf", width=6, height=4)

sc_summaries %>%
  mutate(label=if_else(GeneSet %in% labeled, str_c(GeneSet, str_remove(Summary, ",n:978"), sep="\n"), "")) %>%
  filter(filter != "NoFilter") %>%
  ggplot() +
  geom_point(aes(x=size, y=log(score), color=filter, shape = prot), size=2) +
  geom_text_repel(aes(x=size, y=log(score), label=label), color = "black") +
  theme_bw() +
  theme(legend.title = element_blank())

ggsave(file = "ScolanthusBuscoNoUnfiltered.pdf", width=6, height=4)



genomes_plot <- c(nemVec1="nemVec1",
                  nemVecDovetailChroms="nv_bwa_lowmasked_gapped_chroms",
                  nemVecPacBioOnlyChroms="nv_rm_ph_gapped_chroms",
                  scoCalChroms="sc_purged_chroms")
g_summaries %>%
  filter(names %in% genomes_plot) %>%
  mutate(names = factor(names, levels = genomes_plot, labels=names(genomes_plot))) %>%
  mutate(names = fct_rev(names)) %>%
  pivot_longer(c(Missing, Fragmented, Duplicated, Single), values_to = "Percent")  %>%
  mutate(name=factor(name, levels=c("Single", "Duplicated", "Fragmented", "Missing"))) %>%
  ggplot() +
  geom_col(aes(x=names, y=Percent, fill=name)) +
  coord_flip() +
  theme_gray(base_size = 8) +
  scale_y_continuous(labels = c("0","20","40","60","80","100"), breaks = c(0,20,40,60,80,100)) +
  scale_fill_manual(values = busco_colors,
                    labels =c(" Complete (C) and single-copy (S)  ",
                              " Complete (C) and duplicated (D)",
                              " Fragmented (F)  ",
                              " Missing (M)"))  +
  geom_text(aes(x = names, y = 3, label=Summary), hjust = 0) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "top",
        legend.title = element_blank()) +
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave(file = "InitialSummaries.pdf", height = 4, width = 7)

unpolished_genomes_plot <- c(
                  nemVec1="nemVec1",
                  nemVecDovetailChroms="nv_bwa_lowmasked_gapped_chroms")
p_summaries %>%
  filter(genome == "Dovetail") %>%
  unite(label,genome,polish,sep="") %>%
  mutate(label=factor(label)) %>%
  bind_rows(g_summaries %>%
              filter(names %in% unpolished_genomes_plot) %>%
              mutate(label=factor(names, levels=unpolished_genomes_plot, labels=names(unpolished_genomes_plot)))) %>%
  mutate(label=factor(label, levels=c(
    "nemVec1", "nemVecDovetailChroms", "DovetailMale", "DovetailFemale", "DovetailFemale+Unknown"
  )) %>% fct_rev()) %>%
    pivot_longer(c(Missing, Fragmented, Duplicated, Single), values_to = "Percent")  %>%
  mutate(name=factor(name, levels=c("Single", "Duplicated", "Fragmented", "Missing"))) %>%
  ggplot() +
  geom_col(aes(x=label, y=Percent, fill=name)) +
  coord_flip() +
  theme_gray(base_size = 8) +
  scale_y_continuous(labels = c("0","20","40","60","80","100"), breaks = c(0,20,40,60,80,100)) +
  scale_fill_manual(values = busco_colors,
                    labels =c(" Complete (C) and single-copy (S)  ",
                              " Complete (C) and duplicated (D)",
                              " Fragmented (F)  ",
                              " Missing (M)"))  +
  geom_text(aes(x = label, y = 3, label=Summary), hjust = 0) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "top",
        legend.title = element_blank()) +
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave(file="PolishedSummaries.pdf", height = 4, width = 7)

}
