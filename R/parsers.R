read_blast6 <- function(file) {
  cols <- c("qaccver","saccver","pident","length","mismatch",
            "gapopen","qstart","qend","sstart","send","evalue","bitscore")
  readr::read_tsv(file, col_names = cols, col_types = "ccdiiiiiiidd")
}

read_megablast <- function(file) {
  cols <- c("qseqid", "sseqid", "pident", "length",
            "qstart", "qend", "sstart", "send", "qseq", "sseq")
  readr::read_tsv(file, col_types = "ccddddddcc",  col_names=cols)
}

read_domtbl <- function(file) {
  col_names <- c("target_name","accession","tlen","query_name","query_accession",
                 "qlen","full_E-value","full_score","full_bias",
                 "n","of","c-Evalue","i-Evalue","score","bias",
                 "hmm_from","hmm_to","ali_from","ali_to","env_from","env_to",
                 "acc","descr")
  readr::read_lines("data/rb_.domtbl") %>%
    purrr::discard(stringr::str_starts, pattern="#") %>%
    stringr::str_split(" +", n=23) %>%
    unlist() %>%
    matrix(ncol=23, dimnames=list(NULL, col_names), byrow=T) %>%
    tibble::as_tibble() %>%
    purrr::map_dfc(type.convert)
}


read_psl <- function(file, top_matches_only = T) {
  readr::read_delim(file, delim = "\t", skip = 5, trim_ws = T,
             col_names = c("matches", "misMatches", "repMatches", "nCount",
                           "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert",
                           "strand", "qName", "qSize", "qStart", "qEnd",
                           "tName", "tSize", "tStart", "tEnd",
                           "blockCount", "blockSizes", "qStarts", "tStarts"),
             col_types = "iiiiiiiicciiiciiicccc") %>%
    dplyr::group_by(qName) %>%
    dplyr::filter(matches == max(matches)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Superscaffold_ID = tName %>% stringr::str_remove("PGA_scaffold_") %>% stringr::str_remove("_.*") %>% as.integer())
}

read_gtf <- function(file, ...) {
  cols <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
  readr::read_tsv(file, col_names = cols, col_types="ccciicccc", comment = "#")
}

read_augustus_gtf <- function(file) { read_gtf(file) }

read_simple_gtf <- function(file, geneid_attribute="gene_id", collapse_to = NULL, ...) {
  geneid_regex <- stringr::str_glue('{geneid_attribute} "([^"]*)"')
  gtf <- read_gtf(file)  %>%
    dplyr::filter(stringr::str_detect(attributes, geneid_attribute )) %>%
    dplyr::mutate(gene_id=stringr::str_match(attributes, geneid_regex)[,2])
  if(!is.null(collapse_to)) {
    gtf <- dplyr::group_by(gtf, gene_id) %>%
      dplyr::summarize(seqname=seqname[1], source=source[1], feature=collapse_to,
                       start=min(start), end=max(end), score='.',
                       strand=strand[1], frame='.') %>%
      dplyr::relocate(gene_id, .after=dplyr::last_col())
  }
  gtf

}

read_gtf_with_attributes <- function(file, include_raw_attributes = F) {
  gtf <- read_gtf(file)
  enhanced_gtf <- gtf %>%
    dplyr::mutate(attribute_types=stringr::str_remove_all(attributes, '"[^"]+";') %>% stringr::str_trim() %>% stringr::str_split(" +"),
                  attribute_values=stringr::str_match_all(gtf$attributes,  '"([^"]+)";') %>% purrr::map(~.x[,2]))
  attribute_types <- purrr::flatten_chr(enhanced_gtf$attribute_types) %>% unique()
  attribute_values <- purrr::map2(enhanced_gtf$attribute_types, enhanced_gtf$attribute_values,
                                      ~.y[match(attribute_types, .x)])
  attribute_values.matrix <- do.call(rbind, attribute_values)
  colnames(attribute_values.matrix) <- attribute_types
  if(!include_raw_attributes)
    gtf <- gtf %>% dplyr::select(-attributes)
  dplyr::bind_cols(gtf, as.data.frame(attribute_values.matrix))
}


write_gtf_with_attributes <- function(gtf_with_attributes, path) {
  attributes <- dplyr::select(gtf_with_attributes, -c(1:8), -matches("attributes"))
  attribute_strs <- purrr::map2(colnames(attributes), attributes,
                                    ~ifelse(is.na(.y), "", stringr::str_glue('{.x} "{.y}";'))) %>%
    purrr::pmap_chr(stringr::str_c,  sep = " ")
  gtf <- gtf_with_attributes %>% dplyr::select(c(1:8)) %>% dplyr::mutate(attributes=attribute_strs)
  purrr::pmap(gtf, stringr::str_c, sep="\t") %>%  readr::write_lines(path=path)
}

read_gff <- function(file, ...) {
  readr::read_delim(file, delim = "\t", trim_ws = T, comment = "#",
                    col_names = c("chrom", "source", "feature_type", "start", "end",
                                  "score", "strand", "phase", "attributes"),
                    col_types = "ccciicccc", ...)
}

read_exonerate_gff <- function(file, collapse = T, feature = "gene") {
  gff <- read_gff(file) %>%
    tidyr::separate(attributes, c("geneid", "sequence", "gene_orientation"), " ; ") %>%
    dplyr::mutate_at(c("geneid", "sequence", "gene_orientation"), ~stringr::word(., start = 2, end = 2))

  if(collapse) {
    gff %>%
      dplyr::filter(feature_type == feature) %>%
      dplyr::group_by(sequence, strand, chrom) %>%
      dplyr::mutate(start = min(start), end = max(end)) %>%
      dplyr::distinct() %>%
      dplyr::ungroup()
  }
  else {
    gff
  }
}

read_exon_gtf <- function(file, feature_type = "CDS", geneid_attribute = "protein_id") {
  .mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  read_gtf_with_attributes(file) %>%
    dplyr::filter(feature == feature_type) %>%
    dplyr::grouped_df(geneid_attribute) %>%
    dplyr::summarize(seqname = seqname[1],
                     source = source[1],
                     feature = feature[1],
                     start = min(start),
                     end = max(end),
                     score = ".",
                     strand = .mode(strand),
                     phase = ".") %>%
    dplyr::rename(gene_id = dplyr::all_of(geneid_attribute))
}

read_exon_gff <- function(file, feature = "CDS", geneid_attribute = "Name") {
  .mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  read_gff(file) %>%
    dplyr::filter(feature_type == feature) %>%
    dplyr::mutate(attributes = isolate_gff_attribute(attributes, geneid_attribute)) %>%
    dplyr::rename(geneid=attributes) %>%
    dplyr::group_by(geneid) %>%
    dplyr::summarize(chrom = unique(chrom),
                     source = unique(source),
                     feature_type = unique(feature_type),
                     start = min(start),
                     end = max(end),
                     score = ".",
                     strand = .mode(strand),
                     phase = ".")
}

isolate_gff_attribute <- function(attributes, attribute) {
  stringr::str_remove(attributes, stringr::str_glue(".*{attribute}=")) %>% stringr::str_remove(";.*$")
}

read_simple_gff <- function(file, feature = "mRNA", geneid_attribute = "Name", ...) {
  read_gff(file, ...) %>%
    dplyr::filter(feature_type == feature) %>%
    dplyr::mutate(attributes = isolate_gff_attribute(attributes, geneid_attribute)) %>%
    dplyr::rename(geneid = "attributes")
}

read_wb_gff <- function(file) {
  read_simple_gff(file) %>%
    dplyr::mutate(geneid = strip_geneid_version(geneid))
}

strip_geneid_version <- function(geneids) stringr::str_remove(geneids, "[.]\\d*$")

ens_pep_to_gene_tx_tbl <- function(file) {
  readr::read_lines(file) %>%
    purrr::keep(~ stringr::str_starts(.x, ">")) %>%
    stringr::str_split(" ") %>%
    purrr::map_dfr(~tibble::tibble(pep=stringr::str_remove(.x[1], ">"),
                                   gene=strip_geneid_version(.x[4])))
}

read_repeatmasker_gff3 <- function(file) {
  read_gff(file) %>%
    tidyr::separate(attributes, c("target", "subjStart", "subjEnd"), sep = " ", convert=TRUE) %>%
    dplyr::mutate(target = stringr::str_remove(target, "Target="))
}

read_repeatmasker_out <- function(file) {
  readr::read_table(file, skip = 3,
                    col_types = "idddciicccccccc",
                    col_names = c("SW_Score", "perc_div", "perc_del", "perc_ins",
                                  "chrom", "chrom_begin", "chrom_end", "chrom_left",
                                  "strand", "repeat_type", "repeat_class",
                                  "repeat_begin", "repeat_end", "repeat_left",
                                  "id", "star"))
}

read_bedtools_nuc <- function(file) {
  readr::read_tsv(file, col_types = "ciiddiiiiiii",
                  comment = "#",
                  col_names = c("chrom", "start", "end", "pct_at", "pct_gc",
                                "num_a", "num_c", "num_g", "num_t", "num_n", "num_oth", "length"))
}


read_repeatmasker_tbl <- function(file) {
  lines <- readr::read_lines(file)
  start <- which(stringr::str_detect(lines, "-----"))
  lines <- lines[-seq(1,start)]
  end <- which(stringr::str_detect(lines, "======"))
  lines <- lines[1:end]
  #lines <- lines[-seq(30,length(lines))]
  #data <- stringr::str_split(lines, ":")

  parse_line <- function(line) {
    sline <- stringr::str_squish(line)
    nspaces <- stringr::str_match(line, "^ *")[,1] %>% stringr::str_length()
    if(nspaces > 6) return(tibble::tibble())
    col <- dplyr::case_when(
      nspaces == 0 ~ "class",
      nspaces == 3 ~ "order",
      nspaces == 4 ~ "superfamily",
      nspaces == 5 ~ "family"
    )
    if(stringr::str_detect(sline, "[(:]")) {
      c(first,rest) %<-% stringr::str_split(sline, "[:]|\\(.*,")[[1]]
      fields <- c(first, stringr::str_split(stringr::str_squish(rest), " +")[[1]])
    } else {
      fields <- stringr::str_split(sline, " +")[[1]]
    }
    if(fields[1] == "Total interspersed repeats") {fields <- c(fields[1], "", fields[2:6])}
    if(length(fields) < 6) return(tibble::tibble())
    r <- tibble::tibble(class="", order="", superfamily="", family="", n=fields[2], bp=fields[3], pct=fields[5])
    r[,col] = fields[1]
    r
  }

  purrr::map_dfr(lines, parse_line) %>%
    purrr::map_dfc(utils::type.convert, as.is = T)
}

read_msynt <- function(file) {
  readr::read_tsv(file, col_types = "cciii",
                  col_names = c("chrom", "Geneid", "position", "group", "group_size"))
}

read_show_coords_tab <- function(file) {
  readr::read_tsv(file, skip=4, col_types = "iiiiiidcc",
                  col_names=c("rstart", "rend", "qstart", "qend", "rlen", "qlen", "pct_idy", "ref", "query"))
}

read_emapper <- function(file) {
  readr::read_tsv(file, skip = 4,
                  col_names = c(
                    "tx_id", "seed", "seed_evalue", "seed_score",  "pred_tax_group",
                    "pred_name", "go_terms",   "ec_number", "kegg_ko", "kegg_pathway",
                    "kegg_module", "kegg_reaction", "kegg_rclass", "brite", "kegg_tc",
                    "cazy", "bigg_reaction", "eggnog_tax_scope", "eggnog_ogs",
                    "bestog", "cog_cat", "desc"),
                  col_types = "ccddcccccccccccccccccc",
                  comment = "#")
}

# strip_geneids takes only the first word as the gene id.
read_oma <- function(oma_file, strip_geneids = T, unite_species = T, sep = "_") {
  readLines(oma_file) %>%
    purrr::discard(~"#" == stringr::str_trim(.x) %>% stringr::str_sub(1,1)) %>%
    purrr::map_dfr(~{
      og_gene_sp <- stringr::str_split(.x, "\t")[[1]]
      og <- og_gene_sp[1]
      gene_sp <- stringr::str_split(og_gene_sp[-1], ":", simplify = T)
      genes <- gene_sp[,2]
      sp <- gene_sp[,1]

      if(strip_geneids) {
        genes <- stringr::word(genes)
      }
      if(unite_species) {
        genes <- stringr::str_c(sp, genes, sep = sep)
      }

      tibble::tibble(Geneid = genes, species = sp, group = og)
    })
}
