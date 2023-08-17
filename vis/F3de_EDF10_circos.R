library(circlize)

pal <- scales::hue_pal()




# -- READ IN GENOMES

nv <- fasta_summary("data/reviewed_assemblies/nv/nv.fasta") %>%
  dplyr::mutate(name = stringr::str_c("nv.", stringr::str_remove(name, "chr")),
                color = "white") %>%
  dplyr::arrange(desc(seqnum))
sc <- fasta_summary("data/reviewed_assemblies/sc/sc.fasta") %>%
  dplyr::mutate(name = stringr::str_c("sc.", stringr::str_remove(name, "chr")),
                color = "lightgrey")
bf <- fasta_summary("data/genomes/branchiostoma/GCA_000003815.2_Bfl_VNyyK_genomic.fna") %>%
  dplyr::mutate(ncbi_name=name, name = stringr::str_c("bf.", dplyr::row_number()),
                color = "lightgrey") %>%
  dplyr::slice(1:19)

em <- fasta_summary("data/genomes/ephydatia/Emu_genome_v1.fa") %>%
  dplyr::mutate(emu_name=name, name=stringr::str_c("em.", dplyr::row_number()),
                color="lightgrey") %>%
  dplyr::slice(1:23)


# -- CONSTRUCT GENOME PAIRS FOR PLOTTING
nvsc <- dplyr::bind_rows(nv,sc) %>% dplyr::filter(!stringr::str_detect(name, "Un"))
nvbf <- dplyr::bind_rows(nv,bf) %>% dplyr::filter(!stringr::str_detect(name, "Un"))
nvem <- dplyr::bind_rows(nv,em) %>% dplyr::filter(!stringr::str_detect(name, "Un"))


# -- TRANSLATION TABLES FOR GENOMES
nv_trans_table <- readr::read_tsv("data/assemblies/nv/trans_table.txt",
                               col_names = c("old_name", "new_name"), col_types = "cc")
sc_trans_table <- readr::read_tsv("data/assemblies/sc/trans_table.txt",
                               col_names = c("old_name", "new_name"), col_types = "cc")

# -- READ IN REPEAT LOCATIONS
nv_repeats <- read_repeatmasker_gff3("data/reviewed_assemblies/nv/rm.gff3") %>%
  dplyr::left_join(nv_trans_table, by=c(chrom="old_name")) %>%
  dplyr::select(-chrom, chrom="new_name") %>%
  dplyr::filter(!stringr::str_detect(chrom, "chrUn")) %>%
  dplyr::mutate(chrom = stringr::str_c("nv.", stringr::str_remove(chrom, "chr")))

sc_repeats <- read_repeatmasker_gff3("data/reviewed_assemblies/sc/rm.gff3") %>%
  dplyr::left_join(sc_trans_table, by=c(chrom="old_name")) %>%
  dplyr::select(-chrom, chrom="new_name") %>%
  dplyr::filter(!stringr::str_detect(chrom, "chrUn")) %>%
  dplyr::mutate(chrom = stringr::str_c("sc.", stringr::str_remove(chrom, "chr")))

bf_repeats <- read_repeatmasker_out("data/genomes/branchiostoma/GCA_000003815.2_Bfl_VNyyK_rm.out.gz") %>%
  dplyr::rename(ncbi_name=chrom, start=chrom_begin, end=chrom_end) %>%
  dplyr::inner_join(bf %>% dplyr::select(ncbi_name, chrom=name), by="ncbi_name")

# -- READ IN POSITIONS OF GENES

nv_pos <- read_psl("data/reviewed_assemblies/nv/nve.psl") %>%
  dplyr::select(chrom="tName", Geneid="qName", Start="tStart", End="tEnd") %>%
  dplyr::left_join(nv_trans_table, by=c(chrom="old_name")) %>%
  dplyr::select(-chrom, chrom="new_name") %>%
  dplyr::filter(!stringr::str_detect(chrom, "chrUn")) %>%
  dplyr::mutate(chrom = stringr::str_c("nv.", stringr::str_remove(chrom, "chr")))

sc_pos <-  read_psl("data/reviewed_assemblies/sc/rb.psl") %>%
  dplyr::select(chrom="tName", Geneid="qName", Start="tStart", End="tEnd", Strand="strand") %>%
  dplyr::filter(!stringr::str_detect(chrom, "chrUn")) %>%
  dplyr::mutate(chrom = stringr::str_c("sc.", stringr::str_remove(chrom, "chr")))

bf_pos <- read_exonerate_gff("data/genomes/branchiostoma/Braflo-amphioxus_7u5tJ.genes.gff") %>%
  dplyr::mutate(sequence = stringr::str_replace(sequence, "1$", ".1"), position=start) %>%
  dplyr::rename(Start=start, End=end, Geneid=sequence, Strand=strand) %>%
  dplyr::select(-phase, -geneid, -gene_orientation, -source)

bf_chroms <- bf_pos %>% dplyr::group_by(chrom) %>%
  dplyr::summarize(length=max(position) + 1000, .groups="drop") %>%
  dplyr::arrange(dplyr::desc(length)) %>%
  dplyr::mutate(name=stringr::str_glue("bf.{dplyr::row_number()}"), color="grey") %>% head(n=19)

bf_pos <- bf_pos %>%
  dplyr::inner_join(bf_chroms, by="chrom") %>% dplyr::mutate(chrom=name)

em_pos <- read_exon_gff("data/genomes/ephydatia/Emu_genome_v1.gff", feature="mRNA", geneid_attribute = "Parent")  %>%
  dplyr::select(chrom, Geneid="geneid", Start="start", End="end", Strand="strand") %>%
  dplyr::mutate (chrom = stringr::str_glue("em.{as.integer(stringr::str_remove(chrom, 'scaffold_'))}")) %>%
  dplyr::filter(chrom %in% stringr::str_glue("em.{seq(23)}"))

nv_sc_ogs <- pair_orthogroups("rb", "nv", dir="data/blast_hits")
nv_bf_ogs <- pair_orthogroups("nv", "bf", dir="data/blast_hits")
nv_em_ogs <- pair_orthogroups("nv", "em", dir="data/blast_hits")

bf_nv_algs <- auto_detect_alg_genes(nv_pos, bf_pos, nv_bf_ogs)
nv_sc_algs <- auto_detect_alg_genes(nv_pos, sc_pos, nv_sc_ogs)
nv_em_algs <- auto_detect_alg_genes(nv_pos, em_pos, nv_em_ogs)


# -- CONSTRUCT LINKS BETWEEN POSITIONS IN THE GENOMES
# note that chrom.x must always be the anchor genome
chord_links <- function(algs, genomes) {
  alg_counts <- dplyr::count(algs, chrom.x, chrom.y, alg)
  norm_alg_counts.x <- alg_counts %>%
    dplyr::group_by(chrom.x) %>% dplyr::summarize(total.x=sum(n), .groups="drop") %>%
    dplyr::inner_join(alg_counts, by="chrom.x") %>%
    dplyr::group_by(chrom.x) %>%
    dplyr::arrange(alg) %>%
    dplyr::mutate(End.x = cumsum(n/total.x),
                  Start.x = c(0,End.x) %>% head(-1)) %>%
    dplyr::ungroup()
  norm_alg_counts <- alg_counts %>%
    dplyr::group_by(chrom.y) %>% dplyr::summarize(total.y=sum(n), .groups="drop") %>%
    dplyr::inner_join(norm_alg_counts.x, by="chrom.y") %>%
    dplyr::arrange(alg) %>%
    dplyr::group_by(chrom.y) %>%
    dplyr::mutate(End.y = cumsum(n/total.y),
                  Start.y = c(0,End.y) %>% head(-1)) %>%
    dplyr::ungroup()

  norm_alg_counts %>%
      dplyr::inner_join(genomes %>% dplyr::select(chrom.x=name,length.x=length), by="chrom.x") %>%
      dplyr::mutate_at(c("Start.x", "End.x") , ~.x*length.x) %>%
      dplyr::inner_join(genomes %>% dplyr::select(chrom.y=name,length.y=length), by="chrom.y") %>%
      dplyr::mutate_at(c("Start.y", "End.y"), ~.x*length.y) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(color = factor(chrom.x, labels = stringr::str_c(pal(length(unique(algs$chrom.x))), "50")))
}

bf_nv_chord_links <- chord_links(bf_nv_algs, nvbf)
nv_sc_chord_links <- chord_links(nv_sc_algs, nvsc)
nv_em_chord_links <- chord_links(nv_em_algs, nvem)



# -- READ IN GENES TO PLOT
nv_hox <- readxl::read_xlsx("data/tables/nv_hox_genes.xlsx", sheet="filtered",
                                   col_names=c("gr_name", "nve", "nv1_pos", "nv2_pos", "orientation", "comments", "short_name", "name_comment", "bob_comment")) %>%
  tidyr::separate(nv2_pos, c("pos", "length"), sep=" ") %>%
  tidyr::separate(pos, c("chrom", "range"), sep=":") %>%
  tidyr::separate(range, c("start", "end"), sep="[.][.]", convert=T )

bf_hox <- readxl::read_xlsx("data/tables/bf_hox_genes.xlsx", sheet="filtered")

bf_aa <- Biostrings::readAAStringSet("data/genomes/branchiostoma/Bfl.fastp")
names(bf_aa) <- stringr::str_remove(names(bf_aa), " .*$")
bf_hox_seq <- bf_aa[bf_hox %>% dplyr::filter(!stringr::str_detect(gene_name, "Hox")) %>% dplyr::pull(saccver)]
names(bf_hox_seq) <-
  tibble::tibble(name=names(bf_hox_seq)) %>%
  dplyr::left_join(bf_hox, by=c(name="saccver")) %>%
  dplyr::mutate(combined_name = stringr::str_glue("{gene_name} {name}")) %>%
  dplyr::pull(combined_name)

em_nv_confirmed_mbh <-
  nv_em_ogs %>% dplyr::filter(nv_em_ogs$Geneid %in% (nv_hox$nve %>% toupper())) %>%
  dplyr::inner_join(nv_em_ogs, by="group") %>%
  dplyr::filter(stringr::str_starts(Geneid.y, "Em")) %>%
  dplyr::inner_join(nv_hox %>%
                     dplyr::mutate(nve=toupper(nve)) %>%
                     dplyr::select(nve, short_name),
                   by=c(Geneid.x="nve"))




nv_plot_genes <- nv_hox %>%
  dplyr::select(chrom, start, end, label=short_name) %>%
  dplyr::mutate(chrom=stringr::str_glue("nv.{stringr::str_remove(chrom, 'chr')}"))


sc_plot_genes <- nv_hox %>%
  dplyr::mutate(nve=toupper(nve)) %>%
  dplyr::inner_join(nv_sc_ogs, by=c(nve="Geneid")) %>%
  dplyr::left_join(nv_sc_ogs, by="group") %>%
  dplyr::filter(!stringr::str_starts(Geneid, "NVE")) %>%
  dplyr::inner_join(sc_pos, by="Geneid") %>%
  dplyr::select(chrom="chrom.y", start="Start", end="End", label="short_name") %>%
  dplyr::filter(!is.na(chrom)) %>%
  dplyr::bind_rows(tibble::tibble(chrom=c("sc.4", "sc.10"), start=c(12475674,584017), end=c(12499588,599151), label=c("Engrailed", "Xlox/Cdx")))

bf_plot_genes <-
  bf_hox %>%
  dplyr::inner_join(bf_pos, by=c(saccver="Geneid")) %>%
  dplyr::select(chrom=name, start=position, end=position, label=gene_name)

  #dplyr::mutate(start = start/1e7, end=end/1e7)

em_plot_genes <-
  em_nv_confirmed_mbh %>%
  dplyr::inner_join(em_pos, by=c(Geneid.y="Geneid")) %>%
  dplyr::select(chrom, start=Start, end=End, label=short_name) %>%
  dplyr::bind_rows(
    tibble::tribble(~chrom, ~start, ~end, ~label,
                    "em.13", 5015528, 5015701, "(AmqHex)",
                    "em.13", 5117533, 5117712, "AmqMsx",
                    "em.13", 5038085, 5038255, "AmqNK5/6/7/A",
                    #"em.13", 5821610, 5821774, "(AmqNK5/6/7/B)",
                    "em.16", 7320520, 7320699, "AmqBarH")
  )
  #dplyr::mutate(start = start/1e7, end=end/1e7)



nv_em_plot_genes <- nv_plot_genes %>% dplyr::filter(label %in% c(em_plot_genes$label, "Msx", "Msx21", "Hex",
                                                                 stringr::str_c("Nk", c(3,5,6,7))))

combine_plot_genes <- function(tbl, genes, sep="-") {
  genes_rows <- tibble::enframe(genes, value="label") %>%
    dplyr::inner_join(tbl, by="label") %>%
    dplyr::arrange(start)

  if(is.null(names(genes))) {
    genes_rows <- genes_rows %>% dplyr::mutate(display = label)
  } else {
    genes_rows <- genes_rows %>%
      dplyr::mutate(display=ifelse(stringr::str_length(name) > 0, stringr::str_trim(name), label))
  }
  combined_row <- genes_rows %>%
    dplyr::group_by(chrom) %>%
    dplyr::summarize(chrom=chrom[1],
                     start=min(start), end=max(end),
                     label=stringr::str_c(display %>% purrr::keep(~stringr::str_length(.x) > 0), collapse=sep),
                     .groups="drop")
  tbl %>% dplyr::filter(!label %in% genes) %>%
    dplyr::bind_rows(combined_row)
}






nv_linked_plot_genes <-
  purrr::reduce(
    list(
       c("Nk2.2E", stringr::str_glue("Nk2.2{c('D','C','A','B')}")) %>% purrr::set_names(c("Nk2 (x5)", rep(" ", 4))),
       c("HoxC", D="HoxDa", " "="HoxDb", "Evx"),
       c("Hbn", "Rx", "Otp"),
       #c(`Six1/2`="Six1/2-1", "Goosecoid", `Six1/2`="Six1/2-2"),
       c("Mnx", "Rough"),
       c(Msx="Msx21", "HoxE", R="Anthox9"),
       c(OtxC="OtxC", B="OtxB", A="OtxA"),
       c(MoxA="MoxA", C="MoxC", D="MoxD", B="MoxB"),
       c(`Emx (x2)`="Emx", ` `="Emx2", "Hex"),
       c("Nk6", `7`="Nk7"),
       c("Nk3", `4`="Nk4", Msx="Msx2"),
       c("Nk5", `1`="Nk1"),
       #c("Nk7", "Nk1"),
       c("Dlx", `En-like`="Engrailed") ),
  combine_plot_genes, .init = nv_plot_genes)




nv_linked_plot_genes <-
  purrr::reduce(list(c("Six3/6", "OtxC"),
                     c("HoxB", "Gbx"),
                     c(Tlx="TLX3", Lbx="Ladybird"),
                     #c("HoxC-D-Evx-HoxA", "Mnx-Rough"),
                     c("Six3/6", "OtxC-B-A")),
                combine_plot_genes,
                .init=nv_linked_plot_genes,
                sep = ",")

sc_linked_plot_genes <-
  purrr::reduce(
    list(
      c(stringr::str_glue("Nk2.2{c('D','C','A','B')}")) %>% purrr::set_names(c("Nk2 (x4)", rep(" ", 3))),
      c(HoxD="HoxDa"),
      #c("HoxC", D="HoxDa", " "="HoxDb", "Evx", "HoxA"),
      #c(Alx="Alx1"),
      c("Hbn", "Rx", "Otp"),
      c(`Six1/2 (x2)`="Six1/2-1", ` `="Six1/2-2"),
      c(Mox="MoxB"),
      c(Emx="Emx2"),
      c("Mnx", "Rough"),
      c(Msx="Msx21"), # "HoxE", R="Anthox9"),
      c(`Otx (x3)`="OtxC", ` `="OtxB", ` `="OtxA"),
      c(`Mox (x4)`="MoxC", ` `="MoxA", ` `="MoxD"),
      #c(`Emx (x2)`="Emx", ` `="Emx2", "Hhex"),
      #c("Nk6", Hmx="Hmx3A"),
      #c("Nk3", `4`="Nk4"),
      c(Msx="Msx2", Hmx="Hmx3", "Nk1"),
      c("Dlx", `En-like`="Engrailed"),
      c(Lbx="Ladybird")
      #c("Xlox/Cdx", "Gsx")
      ),
    combine_plot_genes, .init = sc_plot_genes)




  #c("Lhx9", "B-H2-like"),
       #c("Nk1", "Hmx3", "Msx2", "Nk3"),
       #c("Six4/5", "Msx2l"),
       #c("Ro", "HoxA", "Evx", D="HoxDa", "HoxC"),
       #c("B-H2-Like","Lhx9")
bf_linked_plot_genes <-
  purrr::reduce(
    list(c(stringr::str_glue("Hox{seq(15)}"), "EvxA", B="EvxB") %>%
           purrr::set_names(c("B", "EvxA", "15", rep(" ", 13), "Hox1")),
         c("Dlx", "Engrailed"),
         c(`3`="Nk3", `4`="Nk4", "Nk5"),
         c("Hex", "Msx"),
         c("Mnx", "Rough"),
         c("Otp", "Rx"),
         c("Nk2-1", `2`="Nk2-2"),
         c("Six3", `1/2`="Six1/2", `4/5`="Six4/5"),
         c("Otx", "Goosecoid"),
         c("Nk6", `7`="Nk7", Lbx="Ladybird", "Tlx"),
         c(b="Nk1a", Nk1a="Nk1b"),
         c("Gsx", "Xlox", "Cdx")),
    combine_plot_genes, .init = bf_plot_genes
  )

nvsc_genes <- dplyr::bind_rows(sc_linked_plot_genes, nv_linked_plot_genes) %>% as.data.frame()
nvbf_genes <- dplyr::bind_rows(nv_linked_plot_genes, bf_linked_plot_genes) %>% as.data.frame()
nvem_genes <- dplyr::bind_rows(nv_linked_plot_genes, em_plot_genes) %>% as.data.frame()

# nvsc
heat_binsize <- 5e5
nvsc_repeat_dens <- list(nv_repeats, sc_repeats) %>%
  purrr::map_dfr(calculate_feature_genomic_coverage, genome = nvsc, binsize = heat_binsize)
nvsc_gene_dens <-  list(nv_pos, sc_pos) %>%
  purrr::map(dplyr::select, Geneid, start="Start", end = "End", chrom) %>%
  purrr::map_dfr(calculate_feature_density, genome = nvsc, binsize = heat_binsize)

# nvbf
nvbf_repeat_dens <- list(nv_repeats, bf_repeats) %>%
  purrr::map_dfr(calculate_feature_genomic_coverage, genome = nvbf, binsize = heat_binsize)
nvbf_gene_dens <-  list(nv_pos, bf_pos) %>%
  purrr::map(dplyr::select, Geneid, start="Start", end = "End", chrom) %>%
  purrr::map_dfr(calculate_feature_density, genome = nvbf, binsize = heat_binsize)

trunc_to <- function(x, to = 2) {
  f <- 10**(floor(log10(x))-(to-1))
  trunc(x/f) * f
}


plot_circos <- function(outname, pair, links, genes, repeat_dens = NULL, gene_dens = NULL,
                        width = 6, height = 6) {
  fct <- factor(pair$name, levels = pair$name)
  colors <- dplyr::select(pair, name, color) %>% tibble::deframe()

  small.gap <- 1
  big.gap <- 10
  nameTrackHeight <- convert_height(5, "mm")
  gridTrackHeight <- convert_height(2, "mm")
  heatTrackHeight <- convert_height(9, "mm")
  hoxNamesCex <- 0.75

  lens <- rle(pair$assembly)$lengths

  gap.after <-  c(rep(small.gap, lens[1]-1), big.gap, rep(small.gap, lens[2]-1), big.gap)

  circos.par(gap.after = gap.after)
  circos.par(cell.padding = c(0, 0, 0, 0))
  circos.par(track.margin = c(0,0))

  pdf(stringr::str_glue("figures/{outname}.pdf"), width=width, height=height)
  circos.initialize(factors = fct,
                    xlim = cbind(rlang::rep_along(pair$length, 0), pair$length))
  circos.genomicLabels(genes, labels.column = 4, cex = hoxNamesCex, side = "outside", col = "black")

  #axis
  circos.trackPlotRegion(ylim = c(0, .9), bg.border = NA,
                         panel.fun = function(x, y) {
                           xlim = get.cell.meta.data("xlim")
                           current.sector.index = get.cell.meta.data("sector.index")
                           col = colors[current.sector.index]
                           circos.axis("top", major.at = c(xlim[1], trunc_to(xlim[2], to=2)),
                                       labels =round(xlim/1e6),
                                       labels.cex = 0.5,
                                       direction = "inside")
                         }, track.height  = gridTrackHeight)


  # handle naming
  circos.trackPlotRegion(ylim = c(0, 3), bg.border = NA,
                         panel.fun = function(x, y) {
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           current.sector.index = get.cell.meta.data("sector.index")
                           i = get.cell.meta.data("sector.numeric.index")
                           color = colors[current.sector.index]
                           circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], col = color,  border = color)
                           circos.text(mean(xlim), 0, labels = stringr::str_sub(current.sector.index, 4),
                                       cex = 0.8, facing = "inside", #niceFacing = TRUE,
                                       adj = c(0.5,0))
                         }, track.height = nameTrackHeight, track.margin  = c(0.01,0.08))

  col_fun <- colorRamp2(c(-2,0,2), c("#619CFF", "white", "#F8766D"))
  if(!is.null(repeat_dens)) {
    circos.trackPlotRegion( repeat_dens$Chr, repeat_dens$Start, repeat_dens$ScaledValue, #/1e7
                            ylim = c(0,1),track.height = heatTrackHeight, track.margin  = c(0,0.02),
                            panel.fun = function(x,y)  {
                              xend <- pmin(x+heat_binsize, CELL_META$cell.xlim[2])
                              circos.rect(x, 0, xend, 1, border = NA, col = col_fun(y))
                            })
  }

  if(!is.null(gene_dens)) {
    circos.trackPlotRegion(gene_dens$Chr, gene_dens$Start, gene_dens$ScaledValue, #/1e7
                         ylim = c(0,1), track.height = heatTrackHeight, track.margin  = c(0,0.02),
                         panel.fun = function(x,y)  {
                           xend <- pmin(x+heat_binsize, CELL_META$xlim[2])
                           circos.rect(x, 0, xend, 1, border = NA, col = col_fun(y))
                         })
  }

  region1 <- links %>% dplyr::select(chrom=chrom.x, start=Start.x, end=End.x) %>% as.data.frame()
  region2 <- links %>% dplyr::select(chrom=chrom.y, start=Start.y, end=End.y) %>% as.data.frame()
  color <-  links %>% dplyr::pull(color) %>% as.character()
  circos.genomicLink(region1, region2, col=color)
  dev.off()
}

plot_circos("F3d_nvsc", nvsc, nv_sc_chord_links, nvsc_genes, nvsc_repeat_dens, nvsc_gene_dens)
plot_circos("F3e_nvbf", nvbf, bf_nv_chord_links, nvbf_genes, nvbf_repeat_dens, nvbf_gene_dens)
plot_circos("EDF10_nvem", nvem, nv_em_chord_links, nvem_genes, width = 8, height = 8)
