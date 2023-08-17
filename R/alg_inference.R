
#' @title Load gene positions from annotation files
#'
#' @description
#' Load gene positions from annotation files
#'
#' @param cache_name Name of cache file
#' @param dirty If TRUE, recompute the cache
#'
#' @return A list of data frames with columns chrom, Geneid, position
#' @export
#'
load_cached_gene_positions <- function(cache_name = "pos", dirty = F) {
  cache_file <- stringr::str_glue("cache/{cache_name}.Rdata")
  if(dirty || !fs::file_exists(cache_file)) {
    nv_trans_table <- readr::read_tsv("data/assemblies/nv/trans_table.txt",
                                      col_names = c("old_name", "new_name"), col_types = "cc")
    nv_pos <- read_psl("data/reviewed_assemblies/nv/nve.psl") %>% psl_to_position() %>%
      dplyr::left_join(nv_trans_table, by=c(chrom="old_name")) %>% dplyr::rename(old_name = "chrom", chrom="new_name") %>% dplyr::select( -old_name) %>%
      dplyr::mutate(chrom=stringr::str_replace(chrom, "chr", "nv."))
    sc50 <- readr::read_lines("data/reviewed_assemblies/sc/rb_filtered_rm_50pct_tx_ids.txt")
    sc_pos <- read_psl("data/reviewed_assemblies/sc/rb.psl") %>% psl_to_position() %>% dplyr::filter(Geneid %in% sc50) %>%
      dplyr::mutate(chrom=stringr::str_replace(chrom, "chr", "sc."))

    hs_gene_pep <- ens_pep_to_gene_tx_tbl("data/genomes/human/Homo_sapiens.GRCh38.pep.all.flt.fa.bz2")
    hs_pos <- read_simple_gff("data/genomes/human/Homo_sapiens.GRCh38.100.chr.gff3", feature = "gene",
                              geneid_attribute = "ID") %>%
      gff_to_position() %>%
      dplyr::full_join(hs_gene_pep, by=c(Geneid="gene")) %>%
      dplyr::select(-Geneid) %>%
      dplyr::select(chrom, Geneid="pep",  position) %>%
      dplyr::mutate(chrom = stringr::str_c("hs.", chrom))
    dm_names <- c(`NT_037436.4`="dm.3L", `NT_033779.5`="dm.2L", `NT_033778.4`="dm.2R", `NT_033777.3`="dm.3R",
                  `NC_004353.4`="dm.4", `NC_004354.4`="dm.X", `NW_007931083.1`="dm.Un8", `NC_024512.1`="dm.Y",
                  `NW_001844935.1`="dm.3CEN", `NW_001846712.1`="dm.Y2", `NW_001845284.1`="dm.X2", `NC_024511.2`="dm.MT")
    dm_pos <- read_exon_gff("data/genomes/drosophila/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff") %>% gff_to_position() %>%
      dplyr::mutate(chrom=dm_names[chrom] %>% unname())

    ce_pos <- read_wb_gff("data/genomes/celegans/c_elegans.PRJNA13758.WS276.annotations.mRNA_only.gff3") %>% gff_to_position() %>%
      dplyr::group_by(chrom, Geneid) %>% dplyr::summarize(position = min(position), .groups = "drop_last") %>% dplyr::ungroup() %>%
      dplyr::mutate(chrom=stringr::str_c("ce.", chrom))
    em_pos <- read_exon_gff("data/genomes/ephydatia/Emu_genome_v1.gff", geneid_attribute = "Parent", feature = "mRNA") %>% gff_to_position() %>%
      dplyr::mutate(chrom=stringr::str_glue("em.{as.integer(stringr::str_remove(chrom, 'scaffold_'))}"))

    sb_pos <- read_exon_gff("data/genomes/streblspio/Sbenedicti_v2.maker.gff.gz", geneid_attribute = "ID", feature = "mRNA") %>% gff_to_position() %>%
      make_chrom_names("sb")

    aq_pos <- read_simple_gff("data/genomes/amphimedon/Amphimedon_queenslandica.Aqu1.39.gff3", feature="mRNA" ) %>% gff_to_position() %>%
      dplyr::mutate(Geneid=stringr::str_remove(Geneid, "ID=transcript:")) %>% make_chrom_names("aq")
    py_pos <- read_simple_gff("data/genomes/yesso/PY-genome/PYgene-strucanno.gffread.cds.chr.gff3", feature = "gene") %>%
      gff_to_position() %>% make_chrom_names("py", fasta = "data/genomes/yesso/PY-genome/PY.chr.genome.fa")
    re_pos <- read_exon_gff("data/genomes/rhopilema//Rhopilema_esculentum.gff3", feature="mRNA", geneid_attribute = "ID") %>% gff_to_position() %>%
      dplyr::inner_join(readr::read_tsv("data/genomes/rhopilema/trans_table.txt", col_types="cc", col_names=c("re_name", "chrom")), by="chrom") %>%
      dplyr::select(-chrom) %>% dplyr::rename(chrom="re_name")
    bf_pos <- read_exonerate_gff("data/genomes/branchiostoma/Braflo-amphioxus_7u5tJ.genes.gff") %>%
      dplyr::mutate(sequence = stringr::str_replace(sequence, "1$", ".1")) %>%
      exonerate_gff_to_position()
    bf_chroms <- bf_pos %>%
      dplyr::group_by(chrom) %>%
      dplyr::summarize(length=max(position) + 1000, .groups="drop") %>%
      dplyr::arrange(dplyr::desc(length)) %>%
      dplyr::mutate(name=stringr::str_glue("bf.{dplyr::row_number()}")) %>%
      head(n=19)
    bf_pos <- bf_pos %>%
      dplyr::inner_join(bf_chroms, by="chrom") %>%
      dplyr::select(-chrom, -length) %>% dplyr::rename(chrom="name")

    xs_trans_table <- fasta_summary("data/genomes/xenia/xenSp1.scaffolds.fa") %>%
      dplyr::mutate(new_name=stringr::str_glue("xs.{seqnum}")) %>%
      dplyr::select(old_name=name, new_name)
    xs_pos <- read_simple_gff("data/genomes/xenia/xenSp1.gff3", geneid_attribute = "ID") %>% gff_to_position() %>%
      dplyr::left_join(xs_trans_table, by=c(chrom="old_name")) %>%
      dplyr::rename(old_name = "chrom", chrom="new_name") %>%
      dplyr::select(-old_name) %>%
      dplyr::filter()
    lo_pos <- read_exon_gff("data/genomes/spottedgar/Lepisosteus_oculatus.LepOcu1.100.chr.gff3", geneid_attribute = "ID") %>% gff_to_position() %>%
      dplyr::mutate(Geneid = stringr::str_remove(Geneid, "CDS:") %>%
      stringr::str_c(".1")) %>%
      dplyr::mutate(chrom=stringr::str_c("lo.", stringr::str_remove(chrom, "LG")))
    tc_pos <- read_simple_gff("data/genomes/millipede/Trigoniulus_corallinus_hic.gff3.gz", geneid_attribute = "ID") %>%
      gff_to_position() %>% make_chrom_names("tc", fasta = "data/genomes/millipede/Trigoniulus_corallinus_GCA_013389805.1_ASM1338980v1_genomic.scafname.fna")


    aa_pos <- read_simple_gff("data/genomes/aurelia/atlantic/AUR21_r04_wref.gff3", feature = "transcript") %>%
      gff_to_position() %>% make_chrom_names("aa")
    am_pos <- read_simple_gff("data/genomes/millepora/amilgenomev1.1/amil_1.1.maker_006.genes.gff3") %>% gff_to_position() %>%
      dplyr::mutate(Geneid=stringr::str_c("amil.", Geneid)) %>% make_chrom_names("am")
    ch_pos <- read_simple_gff("data/genomes/clytia/merged_transcript_models.gff3", geneid_attribute = "ID", feature = "transcript") %>%
      dplyr::mutate(geneid = stringr::str_c(geneid, "-protein")) %>% gff_to_position() %>% make_chrom_names("ch")

    cr_pos <- read_simple_gff("data/genomes/horseshoe/Annotation/MHSC_MAKER_annotation.gff", feature = "gene", geneid_attribute = "ID") %>%
      gff_to_position() %>% make_chrom_names("cr")
    ep_pos <- read_exon_gff("data/genomes/exaiptasia/GCF_001417965.1_Aiptasia_genome_1.1_genomic.gff") %>% gff_to_position() %>%
      make_chrom_names("ep")
    hm_pos <- read_simple_gff("data/genomes/hydra/all.hydra2.0_genemodels.gff3") %>% gff_to_position() %>%
      make_chrom_names("hm")
    hv_pos <- read_exon_gtf("data/genomes/hydrav/GCF_022113875.1_Hydra_105_v3_genomic.gtf") %>%
      gtf_to_position() %>% make_chrom_names("hv")
    ta_pos <- read_simple_gtf("data/genomes/trichoplax/Triad1_best_genes.Tadh_P.ids.gff", geneid_attribute = "proteinId", collapse_to="gene") %>%
      gtf_to_position() %>% make_chrom_names("ta")

    pos <- list(Aq=aq_pos, Ch=ch_pos, Cr=cr_pos, Ep=ep_pos, Hm=hm_pos, Hv=hv_pos, Ta=ta_pos, Aa=aa_pos, Am=am_pos, Nv=nv_pos, Em=em_pos, Py=py_pos, Sb=sb_pos, Re=re_pos, Bf=bf_pos, Xs=xs_pos, Lo=lo_pos, Tc=tc_pos, Hs=hs_pos, Dm=dm_pos, Sc=sc_pos, Ce=ce_pos) %>%
      purrr::lmap(~list(dplyr::mutate(.x[[1]], Geneid=stringr::str_glue("{names(.x)}_{Geneid}"))) %>%
                    purrr::set_names(names(.x)))
    fs::dir_create(fs::path_dir(cache_file))
    save(pos, file = cache_file)
  }
  load(file = cache_file, envir = .GlobalEnv)
}

#' Perform ancestral linkage group inference
#'
#' This function performs ancestral linkage group inference. The results are
#' cached in a file so that they can be reused.
#' The first step is to determine a weighted graph for each of the critical ndoes
#' that are indicated (e.g., "Bilaterian", "Cnidarian", "Metazoan"). The weights are
#' determined by the number of orthologs co-occurring on the same scaffold. Next, perform
#' optimized leiden clustering. The main objective of the optimization is to avoid
#' pre-determining the number of clusters. The optimization is performed by trying
#' different gamma values. The gamma value that results in the highest modularity
#' is used.
#'
#' @param tree_text A newick tree string
#' @param critical_nodes A vector of node names that must be present in the tree
#' @param ortholog_source Output table from OMA standalone
#' @param dirty If TRUE, the cache will be ignored and the inference will be run
#' from scratch
#' @param cache_name Name of the cache file to use
#' @param min_links Minimum number of orthologs required to infer a link
#' @param min_scaffold_size Minimum number of genes included in a scaffold for use in inference
#' @param gamma_opt_tries Number of random gamma values to try when optimizing gamma
#' @param cluster_opt_tries Number of random clusterings to try when optimizing clustering
#' @param opt_beta Beta value to use when optimizing gamma (see Leiden algorithm)
#' @param gammas A vector of gamma values to try when optimizing gamma
#' @param weight How to calculate edge weight when computing the weighted graph.
#' If "n_links", the number of orthologs is used. If "n_species", the number of species is used.
#' If "total_depth", the sum of the total tree depth for each link to the gene is used.
#' @param use_gamma_prior If TRUE, a prior distribution will be calculated for gamma which
#' is then sampled. The prior distribution is determined based on the modularity of the graphs
#' that are inferred when using gamma values.
#' @param n_iterations Number of iterations to run the inference algorithm.
#' @param min_alg_group_size Minimum orthologous group size to include in the
#' list of groups on the inferred consensus ALGs
infer_cached_algs <- function(tree_text, critical_nodes, ortholog_source,
                              dirty=F, cache_name="algs", min_links = 1,
                              min_scaffold_size = 50, gamma_opt_tries = 100,
                              cluster_opt_tries = 1000, opt_beta = 0.01,
                              gammas =  seq(.5, 1.5, .05),
                              weight = c("n_links", "n_species", "total_depth"),
                              use_gamma_prior = TRUE, n_iterations = 0,
                              min_alg_group_size = 3) {
  load_cached_gene_positions()
  cache_file <- stringr::str_glue("cache/{cache_name}.RData")

  if(dirty || !fs::file_exists(cache_file)) {
    # pull out parameters to use
    weight <- match.arg(weight)
    tree <- treeio::read.newick(text = tree_text)
    species_big <- tree$tip.label
    species <- stringr::str_to_lower(tree$tip.label)

    # check if the critical nodes exist before loading to avoid a delayed error
    subtree_names <- ape::subtrees(tree) %>%
      purrr::keep(~ape::Ntip(.x) > 2) %>%
      purrr::map_chr(~paste0(.x$tip.label, collapse = ","))
    missing <- setdiff(critical_nodes, subtree_names)
    if(length(missing) > 0)
      stop(stringr::str_glue("Nodes {stringr::str_c(missing, collapse=', ')} are missing from the tree {tree_text}."))

    # convert oma input to msynt format
    oma <- read_oma(ortholog_source)
    ogs <- purrr::map_dfr(pos[species_big], create_msynt, oma, min_scaffold_size = min_scaffold_size)

    # compute weighted graph for each critical node
    weighted_subgraphs <-
      ape::subtrees(tree) %>%
      purrr::keep(~ape::Ntip(.x) > 2) %>%
      rlang::set_names(purrr::map_chr(., ~paste0(.x$tip.label, collapse = ","))) %>%
      purrr::map(compute_weighted_graph, ogs, weight = weight, min_weight = min_links)

    # perform clustering for each critical node
    critical_subgraphs <- purrr::map(critical_nodes, ~weighted_subgraphs[[.x]])
    clusters <- critical_subgraphs %>%
      purrr::map(largest_component) %>%
      purrr::map(optimized_leiden_clustering, gamma_opt_tries = gamma_opt_tries,
                 cluster_opt_tries = cluster_opt_tries,
                 opt_beta = opt_beta,
                 n_iterations = n_iterations,
                 use_gamma_prior = use_gamma_prior,
                 gammas = gammas)


    og_szs <- dplyr::distinct(ogs, group, group_size) %>% tibble::deframe()

    # determine ALGs from consensus clusters and give them names based on clade and size
    alg_ogs <- clusters %>%
      purrr::map("gene_support") %>%
      purrr::map(dplyr::select, group, alg=consensus_chrom) %>%
      purrr::map(dplyr::filter, og_szs[group] >= min_alg_group_size) %>%
      purrr::map(dplyr::add_count, alg) %>%
      purrr::map(dplyr::mutate, alg = match(-n, sort(unique(-n))), .keep = "unused") %>%
      purrr::imap(~dplyr::mutate(.x, alg = factor(glue::glue("{.y}.{alg}"))))

    ancestral_msynts <- alg_ogs %>%
      purrr::map(tibble::rowid_to_column, "position") %>%
      purrr::imap(~dplyr::mutate(.x, Geneid = group, species = .y,
                                 chrom = as.character(alg)))

    msynts <- alg_ogs %>%
      purrr::map2(ancestral_msynts,
                  ~dplyr::left_join(dplyr::select(ogs, -alg), .x, by = "group") %>%
                    dplyr::bind_rows(.y))

    save(tree, species, species_big, pos, oma, weighted_subgraphs, clusters,
         ogs, alg_ogs, msynts, file = cache_file)
  }

  load(file = stringr::str_glue(cache_file), envir = .GlobalEnv)
}

#' Optimized clustering algorithm
#'
#' This function performs optimized clustering using the Leiden algorithm.
#' Instead of requiring a resolution (gamma) value, the algorithm tries
#' several gamma values and chooses the one that results in the highest
#' modularity. Optionally, a distribution of gamma values can be inferred
#' for the final `cluster_opt_tries` number of clusterings. This distribution
#' can then be used to sample gamma values.
#'
#' @param g An igraph object
#' @param n_iterations Number of iterations to run the Leiden algorithm
#' @param gammas A vector of gamma values to try
#' @param gamma_opt_tries Number of times to try each gamma value
#' @param cluster_opt_tries Number of times to try each gamma value
#' @param use_gamma_prior If TRUE, use the gamma prior to sample gamma values
#' @param opt_beta The beta value to use for the Leiden algorithm
#'
#' @return A list with the following elements:
#' \itemize{
#'  \item modularity_data: A data frame with the modularity values for each gamma value
#'  \item cluster_sizes: A vector of cluster sizes for each clustering
#'  \item cluster_modularity: A vector of modularity values for each clustering
#'  \item gamma_priors: A data frame with the gamma priors
#'  \item gamma_fit: The gamma prior fit
#'  \item gamma_samples: The gamma values sampled
#'  \item gene_support: The gene support values
#'  \item gene_support_se: The gene support standard errors
#'  \item optimized_gamma: The optimized gamma value
#'  \item optimized_k: The optimized number of clusters
#'  \item optimized_clusters_pool: The optimized cluster pool
#' }
#'
#' @export
optimized_leiden_clustering <- function(g, n_iterations=0, gammas = seq(.5, 1.5, .05),
                                        gamma_opt_tries = 100, cluster_opt_tries = 1000,
                                        use_gamma_prior = TRUE, opt_beta = 0.01) {
  cluster_stats <- function(i, gamma, pb = NULL) {
    if(!is.null(pb)) pb$tick()
    cl <- igraph::cluster_leiden(g,
                                 n_iterations=n_iterations,
                                 objective_function = "modularity",
                                 resolution_parameter = gamma)
    m <- igraph::membership(cl)
    mod <- igraph::modularity(g, m)
    tibble::tibble(gamma = gamma, iter = i, n_clusters = length(cl), modularity = mod)
  }

  # compute the modularity of the clusterings uniformly over the gamma range
  call_df <- purrr::cross_df(list(i = seq(gamma_opt_tries), gamma = gammas))
  pb <- progress::progress_bar$new(total = nrow(call_df),
                                   format = "gamma optimization [:bar] :percent ETA: :eta Elapsed: :elapsedfull",
                                   clear = FALSE)
  mod_df <- purrr::pmap_dfr(call_df, cluster_stats, pb = pb) %>%
    dplyr::arrange(modularity) %>%
    tibble::rowid_to_column("rank") %>%
    dplyr::mutate(p = rank / sum(rank))

  opt_gamma <-  mod_df %>%
    dplyr::group_by(gamma) %>%
    dplyr::summarize(m = mean(modularity)) %>%
    dplyr::filter(m == max(m)) %>%
    dplyr::pull(gamma) %>%
    purrr::pluck(1)

  # get the priors based on their ranks
  gamma_priors <- mod_df %>% dplyr::group_by(gamma) %>%
    dplyr::summarize(prior = sum(p))

  # create a distribution based on a local fit
  fit <- locfit::locfit(prior ~ gamma, data = gamma_priors)
  dist_obj <- distr::AbscontDistribution(d = function(.x) predict(fit, .x),  up1 = max(gammas), low1 = min(gammas))
  gammadist <- distr::r(dist_obj)

  gamma_samples <-
    if(use_gamma_prior)  gammadist(cluster_opt_tries)
    else  rep(opt_gamma, cluster_opt_tries)

  # cluster with the maximized gamma value for n_tries times and select the
  # number of clusters which occurs most frequently
  pb <- progress::progress_bar$new(total = cluster_opt_tries,
                                   format = "cluster optimization [:bar] :percent ETA: :eta Elapsed: :elapsedfull",
                                   clear = F)
  cl <- purrr::map(gamma_samples,
                   ~{pb$tick()
                     igraph::cluster_leiden(g, n_iterations = n_iterations, beta = opt_beta, objective_function = "modularity", resolution_parameter = .x)})

  # weight gamma choice by rank in upper quartile, modularity in upper quartile
  tab <- table(lengths(cl))

  # among the ones with the number of clusters we selected,
  # pick the one that maximizes modularity
  cl_mod <- purrr::map_dbl(cl, ~igraph::modularity(g, igraph::membership(.x)))

  # optimum k is the k for which the mean modularity among the cluster pools is maximal
  optimum_k <- split(cl_mod, lengths(cl)) %>%
    purrr::map(mean) %>%
    which.max() %>%
    names() %>%
    as.numeric()

  gv_support <- compute_gene_support(cl, optimum_k)
  se_support <- compute_gene_support(cl, optimum_k, method = "SE")




  list(modularity_data=mod_df,
       cluster_sizes = lengths(cl),
       cluster_modularity = cl_mod,
       gamma_priors = gamma_priors,
       gamma_fit = fit,
       gamma_samples = gamma_samples,
       gene_support = gv_support,
       gene_support_se = se_support,
       optimized_gamma = opt_gamma,
       optimized_k = optimum_k,
       optimized_clusters_pool = cl)
}

#' Compute gene support for a list of clusterings
#'
#' @param pool A list of clusterings
#' @param k The number of clusters
#' @param method The method to use for computing gene support
#' @param nruns The number of runs to use for computing gene support
#' @param verbose Whether to print progress
#' @return A tibble with the gene support values
#' @export
compute_gene_support <- function(pool, k, method = "GV1", nruns = 5, verbose = F) {
  ens <- purrr::map(pool, ~clue::as.cl_partition(.x$membership)) %>%
    clue::as.cl_ensemble()
  cons <- clue::cl_consensus(ens, method = method, control = list(k = k, nruns= nruns, verbose = verbose))
  consensus_chroms <- apply(as.matrix(cons$.Data), 1, which.max)
  support <- apply(as.matrix(cons$.Data), 1, max)
  tibble::tibble(group = pool[[1]]$names,
                 consensus_chrom = consensus_chroms,
                 support_freq = support,
                 support = support * length(pool))
}
