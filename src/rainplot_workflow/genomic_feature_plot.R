#!/usr/bin/env Rscript

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
if (length(script_arg) == 1) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
} else {
  script_dir <- getwd()
}
activate_path <- file.path(script_dir, "renv", "activate.R")
if (file.exists(activate_path)) {
  source(activate_path)
}

suppressWarnings(suppressMessages({
  library(Gviz)
  library(rtracklayer)
  library(GenomicRanges)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    "Usage: genomic_feature_plot.R <gff_file> <chrom> <start> <end> <output_png>",
    call. = FALSE
  )
}

gff_file <- args[[1]]
chrom_input <- args[[2]]
region_start <- as.integer(args[[3]])
region_end <- as.integer(args[[4]])
output_png <- args[[5]]

seq_aliases <- c(
  "CM007964.1" = "NC_001133.9",
  "CM007965.1" = "NC_001134.8",
  "CM007966.1" = "NC_001135.5",
  "CM007967.1" = "NC_001136.10",
  "CM007968.1" = "NC_001137.3",
  "CM007969.1" = "NC_001138.5",
  "CM007970.1" = "NC_001139.9",
  "CM007971.1" = "NC_001140.6",
  "CM007972.1" = "NC_001141.2",
  "CM007973.1" = "NC_001142.9",
  "CM007974.1" = "NC_001143.9",
  "CM007975.1" = "NC_001144.5",
  "CM007976.1" = "NC_001145.3",
  "CM007977.1" = "NC_001146.8",
  "CM007978.1" = "NC_001147.6",
  "CM007979.1" = "NC_001148.4",
  "CM007980.1" = "NC_001224.1",
  "CM007981.1" = "NC_001224.1",
  "1" = "NC_001133.9",
  "2" = "NC_001134.8",
  "3" = "NC_001135.5",
  "4" = "NC_001136.10",
  "5" = "NC_001137.3",
  "6" = "NC_001138.5",
  "7" = "NC_001139.9",
  "8" = "NC_001140.6",
  "9" = "NC_001141.2",
  "10" = "NC_001142.9",
  "11" = "NC_001143.9",
  "12" = "NC_001144.5",
  "13" = "NC_001145.3",
  "14" = "NC_001146.8",
  "15" = "NC_001147.6",
  "16" = "NC_001148.4",
  "chr1" = "NC_001133.9",
  "chr2" = "NC_001134.8",
  "chr3" = "NC_001135.5",
  "chr4" = "NC_001136.10",
  "chr5" = "NC_001137.3",
  "chr6" = "NC_001138.5",
  "chr7" = "NC_001139.9",
  "chr8" = "NC_001140.6",
  "chr9" = "NC_001141.2",
  "chr10" = "NC_001142.9",
  "chr11" = "NC_001143.9",
  "chr12" = "NC_001144.5",
  "chr13" = "NC_001145.3",
  "chr14" = "NC_001146.8",
  "chr15" = "NC_001147.6",
  "chr16" = "NC_001148.4",
  "chrI" = "NC_001133.9",
  "chrII" = "NC_001134.8",
  "chrIII" = "NC_001135.5",
  "chrIV" = "NC_001136.10",
  "chrV" = "NC_001137.3",
  "chrVI" = "NC_001138.5",
  "chrVII" = "NC_001139.9",
  "chrVIII" = "NC_001140.6",
  "chrIX" = "NC_001141.2",
  "chrX" = "NC_001142.9",
  "chrXI" = "NC_001143.9",
  "chrXII" = "NC_001144.5",
  "chrXIII" = "NC_001145.3",
  "chrXIV" = "NC_001146.8",
  "chrXV" = "NC_001147.6",
  "chrXVI" = "NC_001148.4",
  "MT" = "NC_001224.1",
  "chrM" = "NC_001224.1",
  "chrMT" = "NC_001224.1"
)

resolve_seqname <- function(chrom) {
  resolved <- unname(seq_aliases[chrom])
  if (length(resolved) == 1 && !is.na(resolved)) {
    return(resolved)
  }
  chrom
}

make_placeholder_plot <- function(message_text) {
  png(output_png, width = 3200, height = 500, res = 200)
  par(mar = c(1, 1, 2, 1))
  plot.new()
  title(main = sprintf("Genomic features: %s:%s-%s", chrom_input, region_start, region_end))
  text(0.5, 0.5, message_text, cex = 1)
  dev.off()
}

resolved_seqname <- resolve_seqname(chrom_input)
gff <- import(gff_file)

child_feature_types <- c(
  "mRNA", "transcript", "exon", "CDS", "five_prime_UTR", "three_prime_UTR",
  "pseudogenic_transcript", "pseudogenic_exon", "ncRNA", "lnc_RNA"
)

region_gr <- GRanges(
  seqnames = resolved_seqname,
  ranges = IRanges(start = region_start, end = region_end)
)

features <- subsetByOverlaps(gff, region_gr, ignore.strand = TRUE)
features <- features[as.character(seqnames(features)) == resolved_seqname]
features <- features[!(tolower(as.character(mcols(features)$type)) %in% tolower(child_feature_types))]
features <- features[width(features) > 0]

if (length(features) == 0) {
  make_placeholder_plot("No genomic features found in the requested region.")
  quit(save = "no")
}

feature_type <- as.character(mcols(features)$type)
feature_type[is.na(feature_type) | feature_type == ""] <- "feature"

feature_label <- as.character(mcols(features)$gene)
missing_label <- is.na(feature_label) | feature_label == ""
feature_label[missing_label] <- as.character(mcols(features)$Name)[missing_label]
missing_label <- is.na(feature_label) | feature_label == ""
feature_label[missing_label] <- as.character(mcols(features)$locus_tag)[missing_label]
missing_label <- is.na(feature_label) | feature_label == ""
feature_label[missing_label] <- as.character(mcols(features)$ID)[missing_label]
missing_label <- is.na(feature_label) | feature_label == ""
feature_label[missing_label] <- feature_type[missing_label]

strand_values <- strand(features)
strand_values[strand_values == "*"] <- "+"
strand(features) <- strand_values

mcols(features)$feature <- feature_type
mcols(features)$id <- feature_label
mcols(features)$group <- feature_label

plot_height <- max(3, min(12, 2.5 + (length(features) * 0.18)))

feature_track <- AnnotationTrack(
  range = features,
  genome = "sacCer3",
  chromosome = resolved_seqname,
  name = "Features",
  id = feature_label,
  feature = feature_type,
  group = feature_label,
  stacking = "squish",
  shape = "arrow"
)

displayPars(feature_track) <- list(
  featureAnnotation = "id",
  cex = 0.65,
  fontsize = 8,
  just.group = "above",
  col = "black",
  fill = "#6AA84F",
  fontcolor = "black",
  background.panel = "white",
  col.line = "black"
)

axis_track <- GenomeAxisTrack(
  name = "Position",
  add53 = TRUE,
  littleTicks = TRUE,
  labelPos = "below"
)

png(output_png, width = 3200, height = plot_height * 200, res = 200)
plotTracks(
  list(axis_track, feature_track),
  from = region_start,
  to = region_end,
  chromosome = resolved_seqname,
  main = sprintf("Genomic features: %s:%s-%s", chrom_input, region_start, region_end),
  sizes = c(0.22, 0.78),
  background.title = "white",
  col.title = "black",
  cex.title = 0.9
)
dev.off()
