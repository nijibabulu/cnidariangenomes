read_assembly_file <- function(file) {
  lines <- readr::read_lines(file)
  seq_lines <- purrr::keep(lines, stringr::str_starts, ">")
  scaf_lines <-  purrr::discard(lines, stringr::str_starts, ">")
  list(seqs=readr::read_delim(I(seq_lines), delim=" ",
                              col_names=c("chr", "Superscaffold_ID", "length"), col_types="ccc"),
       scafs=purrr::map_chr(scaf_lines, stringr::str_split, " ", simplify=T))
}

read_assembly <- function(dir, trans_table = NULL, use_assembly_ids=F) {
  mat_files <- fs::dir_ls(dir, glob =  "*.mat")
  res <- stringr::str_match(mat_files, "([[:digit:]]+)k.mat")[,2]
  names <- stringr::str_glue("mat{res}k")
  mats <- purrr::map(mat_files, readr::read_tsv,
                     col_types = "iid", col_names = c("x", "y", "z")) %>%
    purrr::set_names(names)
  ctgs <- readr::read_delim(fs::dir_ls(dir, glob = "*.scaffold_track.txt")[1], delim="\t",
                     col_types="ciiciicciiii")
  chrs <- readr::read_delim(fs::dir_ls(dir, glob = "*.superscaf_track.txt")[1], delim="\t",
                     col_types="ciiciicciiii")
  if(!is.null(trans_table)) {
    chrs <- dplyr::left_join(chrs, trans_table, by=c("Superscaffold_ID"="old_name")) %>%
      dplyr::select(-Superscaffold_ID) %>%
      dplyr::rename(Superscaffold_ID="new_name")
  } else if(use_assembly_ids) {
    chrs <-
      read_assembly_file(fs::dir_ls(dir, glob="*assembly")[1])$seqs %>%
      dplyr::mutate(chr=stringr::str_remove(chr, ">")) %>%
      dplyr::right_join(chrs, by="Superscaffold_ID") %>%
      dplyr::select(-Superscaffold_ID, -length, Superscaffold_ID=chr)
  }
  c(mats, list(chrs=chrs, ctgs=ctgs,  id=dir)) #mat10k=mat10k,
}

# sort an assembly by chromosome size. note that scaffold info is lost here
# currently just designed for global contact maps
# this does not work because of rounding differences; best to simply re-render the matrices
# after manually reordering
sort_assembly_by_chromosome_size <- function(asm) {
  sorted_chrs <- asm$chrs %>%
    dplyr::select(old_start = x1, old_end = x2) %>%
    dplyr::mutate(length=old_end-old_start) %>%
    dplyr::arrange(dplyr::desc(length)) %>%
    dplyr::mutate(new_start = c(0, head(length, -1) %>% cumsum()), new_end = new_start + length)

  move_matrix_coords <- function(old_start, old_end, length, new_start, new_end, v, mat) {
    s <- rlang::ensym(v)
    max_coord <- max(mat$x)
    max_digit <- log10(max_coord) %>% ceiling()
    rnd_level <- purrr::map_lgl(seq(max_digit), ~ round(max_coord/10**.x) * 10**.x == max_coord) %>% sum()

    #rnd_old_start <- mat %>% dplyr::filter(!!s <= old_start) %>% dplyr::pull(!!s) %>% max()
    #rnd_old_end <- mat %>% dplyr::filter(!!s >= old_end) %>% dplyr::pull(!!s) %>% min()
    #rnd_new_start <- mat %>% dplyr::filter(!!s >= new_start) %>% dplyr::pull(!!s) %>% min()
    rnd_old_start <- round(old_start, -rnd_level)
    rnd_old_end <- round(old_end, -rnd_level)
    rnd_new_start <- round(new_start, -rnd_level)
    cat(stringr::str_glue("{rnd_old_start} {rnd_old_end} {rnd_new_start}\n"))
    mat %>% dplyr::filter(!!s >= rnd_old_start & !!s < rnd_old_end) %>%
      dplyr::mutate(!!s := !!s - rnd_old_start + rnd_new_start)
  }

  translate_matrix_coords <- function(mat) {
    max_coord <- max(mat$x)
    max_digit <- log10(max_coord) %>% ceiling()
    rnd_level <- purrr::map_lgl(seq(max_digit), ~ round(max_coord/10**.x) * 10**.x == max_coord) %>% sum()

    transl_tbl <- dplyr::select(sorted_chrs, -length) %>%
      dplyr::mutate_all(~floor(./10**rnd_level)*10**rnd_level) %>%
      purrr::pmap_dfr(function(old_start, old_end, new_start, new_end) {
        tibble::tibble(old=seq(old_start, old_end, by = 10**rnd_level),
               new=seq(new_start, new_end, by = 10**rnd_level))
      })

    mat %>% dplyr::full_join(transl_tbl, by=c(x="old")) %>% dplyr::select(-x) %>% dplyr::rename(x=new) %>%
      dplyr::full_join(transl_tbl, by=c(y="old")) %>% dplyr::select(-y) %>% dplyr::rename(y=new)
  }
  orig_mats <- asm[stringr::str_detect(names(asm), "mat")]
  corr_mats <- purrr::map(orig_mats, translate_matrix_coords)
  #ycor_mats <- purrr::map(xcor_mats, ~purrr::pmap_dfr(sorted_chrs, move_matrix_coords, v="y", mat = .x))

  #xcor_mats <- purrr::map(orig_mats, ~purrr::pmap_dfr(sorted_chrs, move_matrix_coords, v="x", mat = .x))
  #ycor_mats <- purrr::map(xcor_mats, ~purrr::pmap_dfr(sorted_chrs, move_matrix_coords, v="y", mat = .x))
  chrs <- asm$chrs %>%
    dplyr::full_join(sorted_chrs, by=c(x1="old_start", x2="old_end")) %>%
    dplyr::mutate_at(dplyr::vars(tidyr::ends_with("1")), ~new_start) %>%
    dplyr::mutate_at(dplyr::vars(tidyr::ends_with("2")), ~new_end) %>%
    dplyr::select(-new_start, -new_end, -length)

  c(list(chrs=new_chrs), corr_mats)
}


rotate_mat <- function(df, degree, x = "x", y = "y") {
    xs <- df[[x]]
    ys <- df[[y]]
    dfr <- df
    degree <- pi * degree / 180
    l <- sqrt(xs^2 + ys^2)
    teta <- atan(xs / ys)
    teta[is.nan(teta)] <- 0
    dfr[,x] <- l * cos(teta - degree)
    dfr[,y] <- l * sin(teta - degree)
    return(dfr)
  }

prep_matrix <- function(matrix, sym=T) {
  max_z = matrix %>% dplyr::filter(x != y) %>% dplyr::pull(z) %>% max(na.rm=T)
  mat_nodiag <- matrix %>% dplyr::mutate(z=dplyr::if_else(x == y, max_z, z))
  if(sym)
    mat_nodiag <- mat_nodiag %>% dplyr::bind_rows(
      mat_nodiag %>%
        dplyr::filter(x != y) %>%
        dplyr::mutate(tmp=x) %>%
        dplyr::mutate(x=y, y=tmp) %>%
        dplyr::select(-tmp))
  mat_nodiag
}

# plot_genes should be a tibble with the genes to plot by name in a column "name", and a column indicating the family with "family"
plot_global_contact_map <- function(asm, matrix = "mat1000k", n_chrs = 15, gene_positions = NULL, plot_genes = tibble::tibble(), sym = T,
                                    output_dir = NULL, output_format = "pdf", output_height = 7, output_width = 7, limit_quantile = .9999,
                                    output_file = if(is.null(output_dir)) NULL else path(output_dir, stringr::str_glue("{asm$id}_global.{output_format}")),
                                    gradient_colors =  c("white", "white", "white", "orange", "red", "purple", "black"),
                                    family_colors = c("black", "grey25", "grey42"), title ="",
                                    sort_chrs = F, triangular=F, guide_breaks = ggplot2::waiver(), guide_width = NULL, guide_height = NULL,
                                    guide_ticks_color = "black")  {
  if(sort_chrs) {
    asm <- sort_assembly_by_chromosome_size(asm)
  }
  u = if(triangular) sqrt(.5) else 1
  mat_poly <- prep_matrix(asm[[matrix]], sym=sym) %>%
    dplyr::mutate(x = u*x/1e6, y = u*y/1e6) %>%
    dplyr::mutate(ul.x=x,ul.y=y,ur.x=x+u,ur.y=y,lr.x=x+u,lr.y=y+u,ll.x=x,ll.y=y+u) %>%
    tidyr::unite("id", c(x,y), sep=",") %>%
    tidyr::pivot_longer(-c(id,z)) %>%
    tidyr::separate(name, c("corner","coord")) %>%
    tidyr::pivot_wider(names_from=c(coord), id_cols=c(id,corner,z))
  chrs <- asm$chrs[1:n_chrs,] %>%
    dplyr::mutate(ul.x = u*x1, ul.y = u*y1,
                  ur.x = u*x2, ur.y = u*y1,
                  ll.x = u*x1, ll.y = u*y2,
                  lr.x = u*x2, lr.y = u*y2)
  if(triangular) {
    a = 45
    mat_poly <- rotate_mat(mat_poly, a)
    chrs <- chrs %>%
      rotate_mat(a, x="ul.x", y="ul.y") %>% rotate_mat(a, x="ur.x", y="ur.y") %>%
      rotate_mat(a, x="ll.x", y="ll.y") %>% rotate_mat(a, x="lr.x", y="lr.y")
  }
  #mat_poly = mat_poly %>% dplyr::filter(x <= 10, y <= 10)
  p <- ggplot2::ggplot(mat_poly, ggplot2::aes(x=x,y=y,fill=log10(z),group=id)) +
    ggplot2::geom_polygon() +
    #ggplot2::geom_rect(ggplot2::aes(xmin = x1/1e6, xmax = x2/1e6, ymin = y1/1e6, ymax = y2/1e6),
                       #alpha = 0, color = "grey50", data = chrs, inherit.aes = F) +
   ggplot2::geom_segment(ggplot2::aes(x = ul.x/1e6, y = ul.y/1e6, xend = ur.x/1e6, yend = ur.y/1e6),
                          color = "grey50", data = chrs, inherit.aes = F) +
    ggplot2::geom_segment(ggplot2::aes(x = ur.x/1e6, y = ur.y/1e6, xend = lr.x/1e6, yend = lr.y/1e6),
                          color = "grey50", data = chrs, inherit.aes = F) +
    ggplot2::geom_segment(ggplot2::aes(x = lr.x/1e6, y = lr.y/1e6, xend = ll.x/1e6, yend = ll.y/1e6),
                          color = "grey50", data = chrs, inherit.aes = F) +
    ggplot2::geom_segment(ggplot2::aes(x = ll.x/1e6, y = ll.y/1e6, xend = ul.x/1e6, yend = ul.y/1e6),
                          color = "grey50", data = chrs, inherit.aes = F) +
    ggplot2::scale_fill_gradientn(colors = gradient_colors,
                                           name = "log10 Norm. Contact Intensity",
                                  limits=c(0,quantile(log10(mat_poly$z), limit_quantile, na.rm = T)),
                                  breaks = guide_breaks) +
    ggplot2::scale_y_reverse() +
    #ggplot2::scale_x_continuous(position = "top") +
    ggplot2::labs(y = "Genomic Position (MB)") +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(title) +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_text(face = "italic", hjust = 0.5),
                   legend.title = ggplot2::element_text(size = 8)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5,
                                                   barwidth = guide_width,
                                                   barheight = guide_height,
                                                   ticks.colour = guide_ticks_color))
  if(triangular) {
    chrom_centers <- (chrs$x1 + chrs$x2)/2
    #chrom_nrs <- rank(chrs$x1-chrs$x2)
    #chrom_labels <- stringr::str_glue("chr{chrom_nrs}")
    chrom_labels <- chrs$Superscaffold_ID
    p <- p +
      ggplot2::scale_x_continuous(breaks = chrom_centers/1e6,
                                  labels = chrom_labels,
                                  minor_breaks = seq(0, max(chrs$x2), by = 25)) +
      ggplot2::coord_cartesian(ylim = c(0, NA), expand = F) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                     axis.text.y = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     panel.grid.minor.y = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_line())
  } else {
    p <- p +
      ggplot2::coord_equal() +
      ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }
  if(!is.null(gene_positions)) {
    labels <- gene_positions %>%
      dplyr::right_join(plot_genes, by="Geneid") %>%
      dplyr::inner_join(asm$chrs, by=c(chrom="Superscaffold_ID")) %>%
      dplyr::mutate(position = position+x1)
    p <- p + ggplot2::geom_point(ggplot2::aes(x = position/1e6,  y = position/1e6, color = family),
                        size = 2,
                        data = labels,
                        inherit.aes = F) +
      ggrepel::geom_text_repel(ggplot2::aes(x = position/1e6,  y = position/1e6, label = name, color = family),
                      inherit.aes = F,
                      nudge_x = 2,
                      nudge_y = 2,
                      data = labels) +
      ggplot2::scale_color_manual(values = family_colors) +
      ggplot2::guides(color = "none")

  }

  if(!is.null(output_file)) {
    ggplot2::ggsave(plot = p, filename = output_file, height = output_height, width = output_width)
  } else {
    p
  }
}

plot_chromosome_contact_map <- function(asm, id, lims=NULL, matrix = "mat100k", output_dir = NULL, output_format = "pdf", output_filename_extra="", output_height = 7, output_width = 7,
                                        gene_positions = NULL, plot_genes = tibble::tibble(), gene_nudge_x = 0, gene_nudge_y = 0, plot_ctgs = T, sym=T,
                                        title = "", subtitle = "",limit_quantile = .9999, limits = NULL, segments = list(), segment_offset = 1,
                                        title_face = "plain", subtitle_face=title_face, colors = c("white", "white", "orange", "red", "purple", "black"), relative_coords = T, show_guide = T,
                                        output_file = if(is.null(output_dir)) NULL else path(output_dir, stringr::str_glue("{asm$id}_{id}{output_filename_extra}.{output_format}")))  {

  chr <- asm$chrs %>% dplyr::filter(Superscaffold_ID == id) %>% as.list()
  coord_ratio <- chr$sx2/chr$x2
  if(coord_ratio != 1) {
    warning("adjusting matrix coordinates")
  }

  mat_filt <- prep_matrix(asm[[matrix]], sym=sym) %>%
    dplyr::mutate(x = x/coord_ratio, y = y/coord_ratio) %>%
    dplyr::filter(x >= chr$x1, x <= chr$x2, y >= chr$x1, y <= chr$x2) %>%
    dplyr::mutate(x = x - chr$x1, y = y - chr$x1)
  ctgs <- asm$ctgs %>%
    dplyr::filter(x1 >= chr$x1, x2 <= chr$x2, y1 >= chr$x1, y2 <= chr$x2) %>%
    dplyr::mutate(x1 = x1 - chr$x1, x2 = x2 - chr$x1, y1 = y1 - chr$x1, y2 = y2 - chr$x1)

  if(!is.null(lims)) mat_filt <- mat_filt %>% dplyr::filter(x >= lims[1], x <= lims[2], y >= lims[1], y <= lims[2])
  zlimits <- if(!is.null(limits)) limits else c(0,quantile(log10(mat_filt$z), limit_quantile, na.rm = T))
  legend_position <-if(show_guide) { "bottom" } else { "none" }
  p <- ggplot2::ggplot(mat_filt, ggplot2::aes(x = x/1e6, y = y/1e6, fill = log10(z))) +
    ggplot2::geom_raster()

  if(length(segments)) {
    offs <- segment_offset/sqrt(2)
    seg_data <- purrr::map_dfr(segments, ~tibble::tibble(x = .x[1]/1e6+offs, y=.x[1]/1e6-offs, xend=.x[2]/1e6+offs, yend=.x[2]/1e6-offs))
    p <- p + ggplot2::geom_segment(data = seg_data, ggplot2::aes(x=x,y=y,xend=xend,yend=yend), inherit.aes = F)
  }

  if(plot_ctgs) {
    p <- p + ggplot2::geom_rect(ggplot2::aes(xmin = x1/1e6, xmax = x2/1e6, ymin = y1/1e6, ymax = y2/1e6), alpha = 0, color = "grey50", data = ctgs, inherit.aes = F)
  }

  p <- p + ggplot2::scale_y_reverse(expand = c(0,0)) +
    ggplot2::scale_x_continuous(position = "top", expand = c(0,0)) +
    ggplot2::scale_fill_gradientn(colors = colors,
                                  name = "log10 Norm. Contact Intensity",
                                  na.value = "white",
                                  limits=zlimits) +

    #ggplot2::scale_fill_gradientn(colors = colors,
    #name = "Log10 Normalized Contact Intensity") +
    ggplot2::labs(y = "Chromosmal Position (MB)") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::coord_equal() +
    ggplot2::theme(legend.position = legend_position,
                   axis.title.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = title_face, hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(face = subtitle_face, hjust = 0.5),
                   legend.title = ggplot2::element_text(size = 8),
                   panel.grid = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_line()) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5))
  if(!is.null(gene_positions)) {
    labels <- gene_positions %>%
      dplyr::right_join(plot_genes, by="Geneid") %>%
      dplyr::inner_join(asm$chrs, by=c(chrom="Superscaffold_ID")) %>%
      dplyr::mutate(position = position+x1) %>%
      dplyr::filter(position >= lims[1] & position <= lims[2])
    if(nrow(labels) > 0) {
      p <- p + ggplot2::geom_point(ggplot2::aes(x = position/1e6,  y = position/1e6),
                                   size = 2,
                                   data = labels,
                                   inherit.aes = F) +
        ggrepel::geom_text_repel(ggplot2::aes(x = position/1e6,  y = position/1e6, label = name),
                                 inherit.aes = F,
                                 nudge_x = gene_nudge_x,
                                 nudge_y = gene_nudge_y,
                                 data = labels)
    }
  }
  if(!is.null(output_file)) {
    ggplot2::ggsave(filename = output_file, plot = p, height = output_height, width = output_width)
  } else {
    p
  }
}

plot_chromosome_contact_map_new <- function(asm, id, lims=NULL,  scale_kb=100, matrix = stringr::str_glue("mat{scale_kb}k"), output_dir = NULL, output_format = "pdf", output_filename_extra="", output_height = 7, output_width = 7,
                            gene_positions = NULL, plot_genes = tibble::tibble(), gene_nudge_x = 0, gene_nudge_y = 0, plot_ctgs = T,
                            title = "", subtitle = "",limit_quantile = .9999, limits = NULL, triangular=F, sym=T,
                            title_face = "plain", subtitle_face=title_face, colors = c("white", "white", "orange", "red", "purple", "black"), relative_coords = T, show_guide = T,
                            output_file = if(is.null(output_dir)) NULL else path(output_dir, stringr::str_glue("{asm$id}_{id}{output_filename_extra}.{output_format}")))  {
  #mat_sym <- prep_matrix(asm[[matrix]])
  chr <- asm$chrs %>% dplyr::filter(Superscaffold_ID == id) %>% as.list()
  asm_mat_filt <- asm[[matrix]] %>%
    dplyr::filter(x >= chr$x1, x <= chr$x2, y >= chr$x1, y <= chr$x2) %>%
    dplyr::mutate(x = x - chr$x1, y = y - chr$x1)
  if(!is.null(lims)) asm_mat_filt <- asm_mat_filt %>% dplyr::filter(x >= lims[1], x <= lims[2], y >= lims[1], y <= lims[2])
  u = if(triangular) sqrt(.5) else 1
  scale_b <- scale_kb * 1000
  mat_poly <- prep_matrix(asm_mat_filt, sym=sym) %>%
    dplyr::mutate(x = u*x/scale_b, y = u*y/scale_b) %>%
    dplyr::mutate(ul.x=x,ul.y=y,ur.x=x+u,ur.y=y,lr.x=x+u,lr.y=y+u,ll.x=x,ll.y=y+u) %>%
    tidyr::unite("id", c(x,y), sep=",") %>%
    tidyr::pivot_longer(-c(id,z)) %>%
    tidyr::separate(name, c("corner","coord")) %>%
    tidyr::pivot_wider(names_from=c(coord), id_cols=c(id,corner,z))

  ctgs <- asm$ctgs %>% dplyr::filter(x1 >= lims[1], x2 <= lims[2], y1 >= lims[1], y2 <= lims[2])
  zlimits <- if(!is.null(limits)) limits else c(0,quantile(log10(mat_poly$z), limit_quantile, na.rm = T))
  legend_position <-if(show_guide) { "bottom" } else { "none" }
  p <- ggplot2::ggplot(mat_poly, ggplot2::aes(x = x, y = y, fill = log10(z))) +
    ggplot2::geom_polygon()

  if(plot_ctgs) {
    p <- p + ggplot2::geom_rect(ggplot2::aes(xmin = x1/1e6, xmax = x2/1e6, ymin = y1/1e6, ymax = y2/1e6), alpha = 0, color = "grey50", data = ctgs, inherit.aes = F)
  }

  p <- p + ggplot2::scale_y_reverse(expand = c(0,0)) +
    ggplot2::scale_x_continuous(position = "top", expand = c(0,0)) +
    ggplot2::scale_fill_gradientn(colors = colors,
                                  name = "log10 Norm. Contact Intensity",
                                  na.value = "white",
                                  limits=zlimits) +

    #ggplot2::scale_fill_gradientn(colors = colors,
                                  #name = "Log10 Normalized Contact Intensity") +
    ggplot2::labs(y = "Position (MB)") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, subtitle = subtitle) #+

  if(triangular) {
    chrom_centers <- (chrs$x1 + chrs$x2)/2
    #chrom_nrs <- rank(chrs$x1-chrs$x2)
    #chrom_labels <- stringr::str_glue("chr{chrom_nrs}")
    chrom_labels <- chrs$Superscaffold_ID
    p <- p +
      ggplot2::scale_x_continuous(breaks = chrom_centers/1e6,
                                  labels = chrom_labels,
                                  minor_breaks = seq(0, max(chrs$x2), by = 25)) +
      ggplot2::coord_cartesian(ylim = c(0, NA), expand = F) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                     axis.text.y = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     panel.grid.minor.y = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_line())
  } else {
    p <- p +
      ggplot2::coord_equal() +
      ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }
  p <- p +
    ggplot2::theme(legend.position = legend_position,
                   plot.title = ggplot2::element_text(face = title_face, hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(face = subtitle_face, hjust = 0.5),
                   legend.title = ggplot2::element_text(size = 8),
                   panel.grid = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_line()) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5))
  if(!is.null(gene_positions)) {
    labels <- gene_positions %>%
      dplyr::right_join(plot_genes, by="Geneid") %>%
      dplyr::inner_join(asm$chrs, by=c(chrom="Superscaffold_ID")) %>%
      dplyr::mutate(position = position+x1) %>%
      dplyr::filter(position >= lims[1] & position <= lims[2])
    if(nrow(labels) > 0) {
      p <- p + ggplot2::geom_point(ggplot2::aes(x = position/1e6,  y = position/1e6),
                            size = 2,
                            data = labels,
                            inherit.aes = F) +
              ggrepel::geom_text_repel(ggplot2::aes(x = position/1e6,  y = position/1e6, label = name),
                               inherit.aes = F,
                               nudge_x = gene_nudge_x,
                               nudge_y = gene_nudge_y,
                               data = labels)
    }
  }
  if(!is.null(output_file)) {
    ggplot2::ggsave(filename = output_file, plot = p, height = output_height, width = output_width)
  } else {
    p
  }
}

plot_assembly <- function(dir, n_chroms, matrix_to_plot = "mat25k",
                          gradient_colors = c("white", "orange", "red", "purple", "black"), ...) {
  asm <- read_assembly(dir)
  plot_global_contact_map(asm, output_dir = dir)
  purrr::walk(1:n_chroms, ~plot_chromosome(asm = asm, id = .x, output_dir = dir, matrix = matrix_to_plot,
                                           colors = gradient_colors))
}

if(F) {
  plot_assembly("nv_bwa_lowmasked", 15, T)
  plot_assembly("nv_ph_bwa_lowmasked", 15, T)
  plot_assembly("nv_rm_ph", 15, T)

  plot_assembly("sc_purged", 15, F)
  plot_assembly("sc_purged_rm", 15, F)

}
