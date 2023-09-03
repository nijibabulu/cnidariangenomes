## Author: Juan Montenegro Cabrera
##
## Email: jdmontenegroc@gmail.com
## Github: https://github.com/jdmontenegro

### this script uses

### libraries
library(ggplot2)  # easy plotting
library(ggpmisc)  # add R^2 values to linear regression
library(data.table) # easy load zipped tables
library(scales) # scaling options for ggplot
library(ggpubr) # to calculate regression values on the fly for ggplots
library(extrafont) # we need to import "Arial" font type for the plots
library(GenomicFeatures)  # easy load of gene annotations
library(ChIPseeker)  # some useful functions to get gene coordinates
library(dplyr)  # to summarise and group datasets in a tibble
library(readr) # to read in tsv files
library(ggrepel) # to add labels to points
library(ggsci) # to use d3 color palette
library(ggstar) # to use better points figures than standard geom_point

# The following are not required in the description. To install them, please use
# BiocManager::install(c("TxDb.Mmusculus.UCSC.mm10.knownGene",
#                        "TxDb.Dmelanogaster.UCSC.dm6.ensGene",
#                        "TxDb.Celegans.UCSC.ce11.refGene",
#                        "TxDb.Hsapiens.UCSC.hg38.knownGene"))

library("TxDb.Mmusculus.UCSC.mm10.knownGene")		# load mouse annotations
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")			# drosophila dm6 annotations
library("TxDb.Celegans.UCSC.ce11.refGene")				# c elegans ce11 annotations
library("TxDb.Hsapiens.UCSC.hg38.knownGene")				# h sapiens hg38 annotations

### create tibble with Species and genome sizes
species=tibble(Species=factor(c("Nv", "Hv", "Ce", "Dm", "Of", "My", "Lv", "Mm", "Hs"), levels=c("Nv", "Hv", "Ce", "Dm", "Of", "My", "Lv", "Mm", "Hs")),
               sizes=c(252278125, 900935055, 100286401, 143726002, 499544331, 987588634, 869598128, 2730855475, 3099734149))

### load annotations
gtf_nv<-"data/refs/tcs2_internal.cds.full.gtf.gz"
gff_hv<-"data/refs/HVAEP.GeneModels.gff3.gz"
gff_of<-"data/refs/GCA_903813345.2_Owenia_chromosome_genomic.gff.gz"
gff_lv<-"data/refs/GCF_018143015.1_Lvar_3.0_genomic.gff.gz"
gff_my<-"data/refs/GCF_002113885.1_ASM211388v2_genomic.gff.gz"

### create txdb for each species
txdbs=list(
  Nv=makeTxDbFromGFF(gtf_nv),
  Hv=makeTxDbFromGFF(gff_hv),
  Mm=TxDb.Mmusculus.UCSC.mm10.knownGene,
  Hs=TxDb.Hsapiens.UCSC.hg38.knownGene,
  Ce=TxDb.Celegans.UCSC.ce11.refGene,
  Dm=TxDb.Dmelanogaster.UCSC.dm6.ensGene,
  Of=makeTxDbFromGFF(gff_of),
  Lv=makeTxDbFromGFF(gff_lv),
  My=makeTxDbFromGFF(gff_my)
)

atacseqFiles=tibble(Species=species$Species, files=c(
  "data/atac/gastrula_peaks.rm.narrowPeak",
  "data/atac/consensusAEP_noRep.bed",
  "data/atac/merged_ATC.Emb.05.AllAg.AllCell_Ce.bed",
  "data/atac/merged_dm6_embryo_q20.ATACseq.bedGraph",
  "data/atac/ofusiformis_final_peaks.rm.tsv",
  "data/atac/myessoensis_final_peaks.rm.tsv",
  "data/atac/lvariegatus_final_peaks.rm.tsv",
  "data/atac/merged_embr_e11.5.pooled_peaks.narrowPeak",
  "data/atac/merged_ATC.Emb.05.AllAg.AllCell_Hs.bed"
  )
)

### extract genes coordinates

geneGRs<- lapply(
  txdbs, function(x) genes(x))

### load ATACseq peaks
load_peaks<-function(file, sp, tab){
  tmp<-read_tsv(file, col_names=F, skip=1)
  tmp<-tmp[, c(1:3)]
  names(tmp)<-c("chr", "start", "end")
  tmp$Species<-sp
  tab<-bind_rows(tab, tmp)
  return(tab)
}

peaks <-tibble()

for (i in 1:nrow(atacseqFiles)){
  peaks<-load_peaks(atacseqFiles$files[i], atacseqFiles$Species[i], peaks)
}

### function to create peaks GenomicRange list
peakGRs<-makeGRangesFromDataFrame(
    peaks,
	keep.extra.columns=TRUE,
	ignore.strand=TRUE,
	seqinfo=NULL,
	seqnames.field=c("chr"),
	start.field=c("start"),
	end.field=c("end"),
	starts.in.df.are.0based=FALSE
)

### find nearest
minDist<-tibble()
for (i in seq_along(geneGRs)){
  tmp<-tibble()
  sp<-names(geneGRs)[i]
  x<-peakGRs[peakGRs$Species==sp, ]
  y<-geneGRs[[i]]
  tmp<-as_tibble(distanceToNearest(x, y, ignore.strand = TRUE))
#  print(c(sp, length(x), length(y), dim(tmp)))
  tmp$normDist<-tmp$distance*10000000/species$sizes[species$Species==sp]
  tmp$Species<-rep(sp, nrow(tmp))
  minDist<-bind_rows(minDist, tmp[, c(3:5)])
}

# extracting only intergenic
intergenDist<-minDist[minDist$distance>0, ]

### remove outliers
quants <- tapply(
  intergenDist$distance,
  intergenDist$Species,
  quantile
)
Q1s <- sapply(
  1:9,
  function(i) quants[[i]][2]
)
Q3s <- sapply(
  1:9,
  function(i) quants[[i]][4]
)

IQRs <- tapply(
  intergenDist$distance,
  intergenDist$Species,
  IQR
)

Lowers <- Q1s - 1.5*IQRs
Uppers <- Q3s + 1.5*IQRs

datas <- split(intergenDist, intergenDist$Species)

intergenDistClean <- NULL
for (i in 1:9){
  out <- subset(
    datas[[i]],
    datas[[i]]$distance > Lowers[i] & datas[[i]]$distance < Uppers[i]
  )
  intergenDistClean <- tibble(rbind(intergenDistClean, out))
}

intergenDistClean <- left_join(intergenDistClean, species, by="Species")
intergenDistClean$Species<-factor(intergenDistClean$Species, levels=c("Nv", "Hv", "Ce", "Dm", "Of", "My", "Lv", "Mm", "Hs"))

### plotting
# boxplot plot
ggplot(intergenDistClean, aes(x=Species, y=normDist, col=Species)) +
  geom_boxplot(alpha=0.5, aes(fill=Species)) +
  theme_bw() +
  theme(text=element_text(family="Arial", size=12),
        axis.title=element_text(size=14),
        axis.title.y=element_text(margin=margin(r=10)),
        axis.title.x=element_text(margin=margin(t=10)),
        legend.title=element_text(size=14)) +
  labs(x="Species", y="Distance of intergenic ATACseq peaks (bp)")

# violin plot + boxplot
ggplot(intergenDistClean, aes(y=Species, x=distance, col=Species)) +
  geom_violin(trim=FALSE, alpha=0.5, aes(fill=Species)) +
  geom_boxplot(width=0.2, col="black") +
  scale_x_continuous(trans='log10',
                     breaks=c(10,100,1000,10000,100000,1000000),
                     labels=c("10", "100", "1Kb", "10Kb", "100Kb", "1M")) +
  theme_bw() +
  theme(text=element_text(family="Arial", size=14),
        axis.title=element_text(size=16),
        axis.text.y=element_blank(),
        legend.position = "None") +
  labs(y="", x="Distance of genes to open chromatin (bp)")

ggsave("figures/F4d_intergen_dist.png", width=5, height = 8)
ggsave("figures/F4d_intergen_dist.pdf", width=5, height = 8)

# density plots
ggplot(intergenDistClean, aes(x=distance, , col=Species, fill=Species)) +
  geom_density(alpha=0.3, linetype="dashed") +
  theme_bw() +
  theme(text=element_text(family="Arial", size=14),
        axis.title=element_text(size=18),
        legend.title=element_text(size=18),
        axis.title.y=element_text(margin=margin(r=10)),
        axis.title.x=element_text(margin=margin(t=10))) +
  labs(x="Distance to closest gene (bp)", y="Frequency of ATACseq peaks") +
  xlim(0, 5e4)

ggplot(intergenDistClean, aes(x=log(distance+1), col=Species, group=Species, fill=Species)) +
  geom_density(alpha=0.1, linetype="dashed") +
  theme_light() +
  theme(text=element_text(family="Arial", size=14),
        axis.title=element_text(size=18),
        legend.title=element_text(size=18),
        axis.title.y=element_text(margin=margin(r=10)),
        axis.title.x=element_text(margin=margin(t=10))
        ) +
  labs(x="Log of distance to closest gene (bp)", y="Frequency of ATACseq peaks")

# generate a means and median table per Species
distStatsTab<-intergenDistClean %>%
  group_by(Species) %>%
  summarise(
    distMean=mean(distance),
    normDistMean=mean(normDist),
    distMedian=median(distance),
    normDistMeadian=median(normDist),
    genome=mean(sizes)
  )

distStatsTab$bodyOrg<-c("non-bilaterian", "non-bilaterian", rep("bilaterian", 7))
distStatsTab$tads<-c(rep("no TADs", 3), rep("TADs", 6))

# plotting mean distance
ggplot(distStatsTab, aes(x=genome, y=distMean)) +
  theme_bw() +
  stat_poly_line(linewidth=0.2, fill="lightgray", linetype="dashed") +
  stat_poly_eq() +
  geom_point(size=5, aes(col=tads, shape=bodyOrg)) +
  scale_fill_d3() + scale_color_d3() +
  labs(x="Genome size",
       y="Distance from genes to nearest ATAC peak (bp)",
       col="TAD presence",
       shape="Body organization"
      ) +
  theme(text=element_text(family="Arial", size=14),
        axis.title=element_text(size=18),
        axis.title.y=element_text(margin=margin(r=10)),
        axis.title.x=element_text(margin=margin(t=10)),
       ) +
  scale_y_continuous(breaks=seq(0,50000,10000),
                     labels=paste(seq(0,50,10), "Kb")) +
  scale_x_continuous(breaks=seq(0,3e9,5e8),
                     labels=paste(seq(0,3,0.5), "Gb")) +
  geom_text_repel(point.padding=0.1, aes(label=Species))

ggsave("figures/F4e_dist_stats.png", width=6.5, height = 7)
ggsave("figures/F4d_dist_stats.pdf", width=6.5, height = 7)

ggplot(distStatsTab, aes(x=genome, y=distMedian)) +
  theme_bw() +
  stat_poly_line(linewidth=0.2, fill="lightgray", linetype="dashed") +
  stat_poly_eq() +
  geom_point(size=5, aes(col=tads, shape=bodyOrg)) +
  scale_fill_d3() + scale_color_d3() +
  labs(x="Genome size",
       y="Distance from genes to nearest ATAC peak",
       col="TAD presence",
       shape="Body organization"
      ) +
  theme(text=element_text(family="Arial", size=14),
        axis.title=element_text(size=18),
        axis.title.y=element_text(margin=margin(r=10)),
        axis.title.x=element_text(margin=margin(t=10)),
       ) +
  scale_y_continuous(breaks=seq(0,50000,10000),
                     labels=paste(seq(0,50,10), "Kb")) +
  scale_x_continuous(breaks=seq(0,3e9,5e8),
                     labels=paste(seq(0,3,0.5), "Gb")) +
  geom_text_repel(point.padding=0.1, aes(label=Species))

