circos.genomicLabels.fix <-
  function (bed, labels = NULL, labels.column = NULL, facing = "clockwise",
          niceFacing = TRUE, col = par("col"), cex = 0.8, font = par("font"),
          padding = 0.4, connection_height = convert_height(5, "mm"),
          line_col = par("col"), line_lwd = par("lwd"), line_lty = par("lty"),
          labels_height = min(c(convert_height(1.5, "cm"), max(strwidth(labels,
                                                                        cex = cex, font = font)))), side = c("inside", "outside"),
          track.margin = circos.par("track.margin"))
{
  bed = circlize:::validate_data_frame(bed)
  circlize:::validate_region(bed)
  if (is.null(labels)) {
    labels = bed[[labels.column]]
  }
  bed[[4]] = as.vector(labels)
  bed[[1]] = factor(as.vector(bed[[1]]), levels = get.all.sector.index())
  od = order(bed[[1]], bed[[2]])
  bed = bed[od, , drop = FALSE]
  bed[[1]] = as.vector(bed[[1]])
  if (length(col) == 1)
    col = rep(col, nrow(bed))
  if (length(cex) == 1)
    cex = rep(cex, nrow(bed))
  if (length(font) == 1)
    font = rep(font, nrow(bed))
  if (length(line_col) == 1)
    line_col = rep(line_col, nrow(bed))
  if (length(line_lwd) == 1)
    line_lwd = rep(line_lwd, nrow(bed))
  if (length(line_lty) == 1)
    line_lty = rep(line_lty, nrow(bed))
  col = col[od]
  cex = cex[od]
  font = font[od]
  line_col = line_col[od]
  line_lwd = line_lwd[od]
  line_lty = line_lty[od]
  chr = get.all.sector.index()[1]
  sector_data = circlize:::get.sector.data(chr)
  chr_width = sector_data["start.degree"] - sector_data["end.degree"]
  extend = (360 - chr_width)/chr_width
  extend = c(0, extend)
  all_chr = unique(bed[, 1])
  bed2 = NULL
  for (cr in all_chr) {
    sub_bed = bed[bed[, 1] == cr, ]
    if (cr != chr) {
      x1 = reverse.circlize(circlize(sub_bed[, 2], y = rep(1,
                                                           nrow(sub_bed)), sector.index = cr), sector.index = chr)[,
                                                                                                                   1]
      x2 = reverse.circlize(circlize(sub_bed[, 3], y = rep(1,
                                                           nrow(sub_bed)), sector.index = cr), sector.index = chr)[,
                                                                                                                   1]
      sub_bed[, 2:3] = data.frame(start = x1, end = x2)
    }
    bed2 = rbind(bed2, sub_bed)
  }
  bed2[, 1] = chr
  if (!facing %in% c("clockwise", "reverse.clockwise")) {
    stop_wrap("facing can only be 'clockwise' or `reverse.clockwise`.")
  }
  op = circos.par("points.overflow.warning")
  circos.par(points.overflow.warning = FALSE)
  if (side == "inside") {
    circos.genomicTrackPlotRegion(bed2, ylim = c(0, 1),
                                  track.margin = c(convert_height(0.5, "mm"), track.margin[2]),
                                  cell.padding = c(0, 0, 0, 0), track.height = connection_height,
                                  bg.border = NA)
    i_track = get.cell.meta.data("track.index")
    circos.genomicTrackPlotRegion(bed2, ylim = c(0, 1),
                                  panel.fun = function(region, value, ...) {
                                    l = bed2[[1]] == CELL_META$sector.index
                                    circos.genomicText(region, value, y = 1, labels.column = 1,
                                                       facing = facing, adj = c((facing == "clockwise") +
                                                                                  0, 0.5), posTransform = posTransform.text,
                                                       col = col[l], cex = cex[l], font = font[l],
                                                       niceFacing = niceFacing, padding = padding,
                                                       extend = extend)
                                  }, track.height = labels_height, bg.border = NA,
                                  track.margin = c(track.margin[1], 0), cell.padding = c(0,
                                                                                         0, 0, 0))
    circos.genomicPosTransformLines(bed2, posTransform = function(region,
                                                                  value) {
      l = bed2[[1]] == get.cell.meta.data("sector.index",
                                          track.index = i_track + 1)
      posTransform.text(region, y = 1, labels = value[[1]],
                        cex = cex[l], font = font[l], track.index = i_track +
                          1, padding = padding, extend = extend)
    }, direction = "inside", horizontalLine = "top", track.index = i_track,
    cell.padding = c(0, 0, 0, 0), col = line_col, lwd = line_lwd,
    lty = line_lty)
  }
  else {
    circos.genomicTrackPlotRegion(bed2, ylim = c(0, 1),
                                  panel.fun = function(region, value, ...) {
                                    l = bed2[[1]] == CELL_META$sector.index
                                    circos.genomicText(region, value, y = 0, labels.column = 1,
                                                       facing = facing, adj = c((facing == "reverse.clockwise") +
                                                                                  0, 0.5), posTransform = posTransform.text,
                                                       col = col[l], cex = cex[l], font = font[l],
                                                       niceFacing = niceFacing, padding = padding,
                                                       extend = extend)
                                  }, track.height = labels_height, bg.border = NA,
                                  track.margin = c(0, track.margin[2]), cell.padding = c(0,
                                                                                         0, 0, 0))
    i_track = get.cell.meta.data("track.index")
    circos.genomicPosTransformLines(bed2, posTransform = function(region,
                                                                  value) {
      l = bed2[[1]] == get.cell.meta.data("sector.index",
                                          track.index = i_track)
      posTransform.text(region, y = 0, labels = value[[1]],
                        cex = cex[l], font = font[l], track.index = i_track,
                        padding = padding, extend = extend)
    }, direction = "outside", horizontalLine = "bottom",
    col = line_col, lwd = line_lwd, lty = line_lty,
    track.height = connection_height, track.margin = c(track.margin[1],
                                                       convert_height(0.5, "mm")), cell.padding = c(0,
                                                                                                    0, 0, 0))
  }
  circos.par(points.overflow.warning = op)
}
