#!/usr/bin/env Rscript

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
if (length(script_arg) == 1) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
} else {
  script_dir <- getwd()
}

workflow_root <- dirname(dirname(script_dir))
g4_bed_file <- file.path(workflow_root, "data", "bed", "g4.motifs.bed")

local_r_library <- Sys.getenv("R_LIBS_USER", unset = "")
if (!nzchar(local_r_library)) {
  local_r_library <- file.path(script_dir, ".r_library")
}

.libPaths(c(local_r_library, .libPaths()))

required_pkgs <- c("Gviz", "GenomicRanges")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE, lib.loc = local_r_library)]

if (length(missing_pkgs) > 0) {
  stop(
    sprintf(
      "Missing required R packages in local R library (%s): %s\nRun ensure_r_environment.R first.",
      local_r_library,
      paste(missing_pkgs, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressWarnings(suppressMessages({
  library(Gviz)
  library(GenomicRanges)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    "Usage: genomic_feature_plot.R <chrom> <start> <end> <output_png> <ncbi_gff>",
    call. = FALSE
  )
}

chrom_input <- args[[1]]
region_start <- as.integer(args[[2]])
region_end <- as.integer(args[[3]])
output_png <- args[[4]]
gff_file <- args[[5]]

if (is.na(region_start) || is.na(region_end) || region_start > region_end) {
  stop("Start and end coordinates must be valid integers with start <= end.", call. = FALSE)
}

if (!file.exists(gff_file)) {
  stop(sprintf("NCBI GFF file not found: %s", gff_file), call. = FALSE)
}

seq_aliases <- c(
  "CM007964.1" = "chrI",
  "CM007965.1" = "chrII",
  "CM007966.1" = "chrIII",
  "CM007967.1" = "chrIV",
  "CM007968.1" = "chrV",
  "CM007969.1" = "chrVI",
  "CM007970.1" = "chrVII",
  "CM007971.1" = "chrVIII",
  "CM007972.1" = "chrIX",
  "CM007973.1" = "chrX",
  "CM007974.1" = "chrXI",
  "CM007975.1" = "chrXII",
  "CM007976.1" = "chrXIII",
  "CM007977.1" = "chrXIV",
  "CM007978.1" = "chrXV",
  "CM007979.1" = "chrXVI",
  "CM007980.1" = "chrM",
  "CM007981.1" = "chrM",
  "NC_001133.9" = "chrI",
  "NC_001134.8" = "chrII",
  "NC_001135.5" = "chrIII",
  "NC_001136.10" = "chrIV",
  "NC_001137.3" = "chrV",
  "NC_001138.5" = "chrVI",
  "NC_001139.9" = "chrVII",
  "NC_001140.6" = "chrVIII",
  "NC_001141.2" = "chrIX",
  "NC_001142.9" = "chrX",
  "NC_001143.9" = "chrXI",
  "NC_001144.5" = "chrXII",
  "NC_001145.3" = "chrXIII",
  "NC_001146.8" = "chrXIV",
  "NC_001147.6" = "chrXV",
  "NC_001148.4" = "chrXVI",
  "NC_001224.1" = "chrM",
  "1" = "chrI",
  "2" = "chrII",
  "3" = "chrIII",
  "4" = "chrIV",
  "5" = "chrV",
  "6" = "chrVI",
  "7" = "chrVII",
  "8" = "chrVIII",
  "9" = "chrIX",
  "10" = "chrX",
  "11" = "chrXI",
  "12" = "chrXII",
  "13" = "chrXIII",
  "14" = "chrXIV",
  "15" = "chrXV",
  "16" = "chrXVI",
  "chr1" = "chrI",
  "chr2" = "chrII",
  "chr3" = "chrIII",
  "chr4" = "chrIV",
  "chr5" = "chrV",
  "chr6" = "chrVI",
  "chr7" = "chrVII",
  "chr8" = "chrVIII",
  "chr9" = "chrIX",
  "chr10" = "chrX",
  "chr11" = "chrXI",
  "chr12" = "chrXII",
  "chr13" = "chrXIII",
  "chr14" = "chrXIV",
  "chr15" = "chrXV",
  "chr16" = "chrXVI",
  "chri" = "chrI",
  "chrii" = "chrII",
  "chriii" = "chrIII",
  "chriv" = "chrIV",
  "chrv" = "chrV",
  "chrvi" = "chrVI",
  "chrvii" = "chrVII",
  "chrviii" = "chrVIII",
  "chrix" = "chrIX",
  "chrx" = "chrX",
  "chrxi" = "chrXI",
  "chrxii" = "chrXII",
  "chrxiii" = "chrXIII",
  "chrxiv" = "chrXIV",
  "chrxv" = "chrXV",
  "chrxvi" = "chrXVI",
  "MT" = "chrM",
  "chrMT" = "chrM",
  "chrMito" = "chrM",
  "Mito" = "chrM"
)

resolve_chromosome <- function(chrom) {
  resolved <- unname(seq_aliases[chrom])
  if (length(resolved) > 0 && !is.na(resolved[[1]])) {
    return(resolved[[1]])
  }

  chrom_lower <- tolower(chrom)
  resolved <- unname(seq_aliases[chrom_lower])
  if (length(resolved) > 0 && !is.na(resolved[[1]])) {
    return(resolved[[1]])
  }

  if (grepl("^chr", chrom, ignore.case = TRUE)) {
    suffix <- sub("^chr", "", chrom, ignore.case = TRUE)
    return(paste0("chr", toupper(suffix)))
  }

  paste0("chr", chrom)
}

parse_attributes <- function(attr_string) {
  attrs <- list()
  if (!nzchar(attr_string)) {
    return(attrs)
  }

  fields <- strsplit(attr_string, ";", fixed = TRUE)[[1]]
  for (field in fields) {
    field <- trimws(field)
    if (!nzchar(field)) {
      next
    }

    if (grepl("=", field, fixed = TRUE)) {
      parts <- strsplit(field, "=", fixed = TRUE)[[1]]
      key <- trimws(parts[[1]])
      value <- trimws(paste(parts[-1], collapse = "="))
      attrs[[key]] <- utils::URLdecode(value)
    } else if (grepl(" ", field, fixed = TRUE)) {
      parts <- strsplit(field, " ", fixed = TRUE)[[1]]
      key <- trimws(parts[[1]])
      value <- trimws(gsub('^"|"$', "", paste(parts[-1], collapse = " ")))
      attrs[[key]] <- utils::URLdecode(value)
    }
  }

  attrs
}

first_non_empty <- function(values) {
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0) {
    return(NA_character_)
  }
  values[[1]]
}

note_token <- function(attrs) {
  note_value <- first_non_empty(c(attrs[["Note"]], attrs[["description"]]))
  if (is.na(note_value)) {
    return(NA_character_)
  }

  trimws(strsplit(note_value, ";", fixed = TRUE)[[1]][[1]])
}

normalize_feature_type <- function(feature_type) {
  if (feature_type == "long_terminal_repeat") {
    return("repeat_region")
  }
  feature_type
}

choose_feature_name <- function(feature_type, attrs) {
  if (feature_type == "gene") {
    return(first_non_empty(c(attrs[["Name"]], attrs[["gene"]], attrs[["locus_tag"]], attrs[["ID"]], "Gene")))
  }

  if (feature_type %in% c("tRNA", "rRNA", "ncRNA", "snoRNA", "snRNA")) {
    return(first_non_empty(c(attrs[["gene"]], attrs[["product"]], attrs[["Name"]], attrs[["locus_tag"]], note_token(attrs), attrs[["ID"]], feature_type)))
  }

  if (feature_type %in% c("telomere", "centromere", "origin_of_replication", "repeat_region", "mobile_genetic_element")) {
    return(first_non_empty(c(note_token(attrs), attrs[["Name"]], attrs[["gene"]], attrs[["locus_tag"]], attrs[["ID"]], feature_type)))
  }

  first_non_empty(c(attrs[["Name"]], attrs[["gene"]], attrs[["product"]], attrs[["locus_tag"]], attrs[["ID"]], feature_type))
}

keep_feature_row <- function(feature_type, attrs) {
  wanted <- c(
    "gene",
    "tRNA",
    "rRNA",
    "ncRNA",
    "snoRNA",
    "snRNA",
    "centromere",
    "telomere",
    "repeat_region",
    "origin_of_replication",
    "mobile_genetic_element"
  )

  if (!feature_type %in% wanted) {
    return(FALSE)
  }

  if (feature_type == "gene") {
    gene_biotype <- first_non_empty(c(attrs[["gene_biotype"]], ""))
    if (gene_biotype %in% c("tRNA", "rRNA", "ncRNA", "snoRNA", "snRNA")) {
      return(FALSE)
    }
  }

  if (feature_type == "centromere") {
    note_value <- first_non_empty(c(attrs[["Note"]], ""))
    if (grepl("CDEI|CDEII|CDEIII", note_value)) {
      return(FALSE)
    }
  }

  TRUE
}

make_placeholder_plot <- function(message_text, plot_start, plot_end, outpath) {
  png(outpath, width = 3200, height = 900, res = 200)
  par(mar = c(1, 1, 2, 1))
  plot.new()
  title(main = sprintf("W303 NCBI annotations: %s:%s-%s", chrom_input, plot_start, plot_end))
  text(0.5, 0.55, message_text, cex = 1.05)
  text(0.5, 0.42, sprintf("Annotation source: %s", basename(gff_file)), cex = 0.9)
  dev.off()
}

read_features_from_gff <- function(path, requested_chrom, region_start, region_end) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)

  rows <- list()
  row_index <- 1L

  repeat {
    lines <- readLines(con, n = 50000, warn = FALSE)
    if (length(lines) == 0) {
      break
    }

    for (line in lines) {
      if (!nzchar(line) || startsWith(line, "#")) {
        next
      }

      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(fields) != 9) {
        next
      }

      seqid <- fields[[1]]
      feature_type <- normalize_feature_type(fields[[3]])

      if (resolve_chromosome(seqid) != requested_chrom) {
        next
      }

      attrs <- parse_attributes(fields[[9]])
      if (!keep_feature_row(feature_type, attrs)) {
        next
      }

      start_1based <- suppressWarnings(as.integer(fields[[4]]))
      end_1based <- suppressWarnings(as.integer(fields[[5]]))
      if (is.na(start_1based) || is.na(end_1based)) {
        next
      }

      feature_start <- start_1based - 1L
      feature_end <- end_1based
      if (feature_end <= region_start || feature_start >= region_end) {
        next
      }

      feature_name <- choose_feature_name(feature_type, attrs)
      strand_value <- fields[[7]]
      if (!strand_value %in% c("+", "-")) {
        strand_value <- "*"
      }

      rows[[row_index]] <- data.frame(
        chrom = requested_chrom,
        start = max(region_start, feature_start),
        end = min(region_end, feature_end),
        raw_start = feature_start,
        raw_end = feature_end,
        name = feature_name,
        strand = strand_value,
        feature_type = feature_type,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }

  if (length(rows) == 0) {
    return(data.frame())
  }

  features <- do.call(rbind, rows)
  features <- unique(features)
  rownames(features) <- NULL
  features
}

feature_rank <- c(
  "gene" = 1,
  "origin_of_replication" = 2,
  "centromere" = 3,
  "telomere" = 4,
  "tRNA" = 5,
  "rRNA" = 6,
  "snoRNA" = 7,
  "snRNA" = 8,
  "ncRNA" = 9,
  "repeat_region" = 10,
  "mobile_genetic_element" = 11,
  "g4_motif" = 12
)

feature_color <- function(feature_type) {
  switch(
    feature_type,
    "gene" = "#1F4E79",
    "origin_of_replication" = "#38761D",
    "centromere" = "#7F6000",
    "telomere" = "#CC0000",
    "tRNA" = "#8E24AA",
    "rRNA" = "#6A1B9A",
    "snoRNA" = "#AB47BC",
    "snRNA" = "#AD1457",
    "ncRNA" = "#5E35B1",
    "repeat_region" = "#F39C12",
    "mobile_genetic_element" = "#D35400",
    "g4_motif" = "#008080",
    "#4A4A4A"
  )
}

feature_type_label <- function(feature_type) {
  switch(
    feature_type,
    "gene" = "Genes",
    "origin_of_replication" = "Origins",
    "centromere" = "Centromeres",
    "telomere" = "Telomeres",
    "tRNA" = "tRNAs",
    "rRNA" = "rRNAs",
    "snoRNA" = "snoRNAs",
    "snRNA" = "snRNAs",
    "ncRNA" = "ncRNAs",
    "repeat_region" = "Repeat regions",
    "mobile_genetic_element" = "Mobile genetic elements",
    "g4_motif" = "G4 motifs",
    feature_type
  )
}

is_directional_feature <- function(feature_type, strand_values) {
  has_strand <- any(strand_values %in% c("+", "-"))
  feature_type %in% c(
    "gene",
    "tRNA",
    "rRNA",
    "ncRNA",
    "snoRNA",
    "snRNA",
    "mobile_genetic_element"
  ) && has_strand
}

draw_direction_backdrop <- function(region_start, region_end, y, strand_value, color, region_span) {
  if (!strand_value %in% c("+", "-")) {
    return(invisible(NULL))
  }

  step <- max(400L, ceiling(region_span / 32))
  arrow_len <- step * 0.7
  centers <- seq(region_start + step / 2, region_end - step / 2, by = step)
  if (length(centers) == 0) {
    return(invisible(NULL))
  }

  backdrop_color <- adjustcolor(color, alpha.f = 0.22)
  if (strand_value == "+") {
    x0 <- centers - arrow_len / 2
    x1 <- centers + arrow_len / 2
  } else {
    x0 <- centers + arrow_len / 2
    x1 <- centers - arrow_len / 2
  }

  arrows(
    x0 = x0,
    y0 = rep(y, length(centers)),
    x1 = x1,
    y1 = rep(y, length(centers)),
    length = 0.06,
    angle = 22,
    code = 2,
    col = backdrop_color,
    lwd = 1
  )
}

draw_feature_ruler <- function(feature_row, y, region_start, region_end, region_span) {
  feature_type <- feature_row$feature_type[[1]]
  strand_value <- feature_row$strand[[1]]
  track_color <- feature_color(feature_type)
  ruler_color <- adjustcolor("#4A4A4A", alpha.f = 0.65)
  label_y <- y + 0.21
  tick_half_height <- 0.16

  segments(region_start, y, region_end, y, col = ruler_color, lwd = 1)

  if (is_directional_feature(feature_type, strand_value)) {
    draw_direction_backdrop(region_start, region_end, y, strand_value, track_color, region_span)
  }

  feature_start <- feature_row$start[[1]]
  feature_end <- feature_row$end[[1]]
  feature_width <- max(1L, feature_end - feature_start)

  segments(feature_start, y, feature_end, y, col = track_color, lwd = 2.6)
  segments(
    c(feature_start, feature_end),
    c(y - tick_half_height, y - tick_half_height),
    c(feature_start, feature_end),
    c(y + tick_half_height, y + tick_half_height),
    col = track_color,
    lwd = 2.2
  )

  if (feature_width < max(2L, ceiling(region_span * 0.0015))) {
    center_x <- (feature_start + feature_end) / 2
    segments(center_x, y - 0.22, center_x, y + 0.22, col = track_color, lwd = 2.4)
  }

  label_x <- min(region_end, max(region_start, (feature_start + feature_end) / 2))
  text(
    x = label_x,
    y = label_y,
    labels = feature_row$name[[1]],
    cex = 0.58,
    font = 2,
    col = "black"
  )
}

read_g4_from_bed <- function(path, requested_chrom, region_start, region_end) {
  if (!file.exists(path)) {
    return(data.frame())
  }

  g4_df <- tryCatch(
    read.delim(path, header = FALSE, stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )

  if (nrow(g4_df) == 0) {
    return(data.frame())
  }

  if (ncol(g4_df) < 3) {
    return(data.frame())
  }

  colnames(g4_df)[1:3] <- c("chrom", "start", "end")
  if (ncol(g4_df) >= 4) {
    colnames(g4_df)[4] <- "name"
  } else {
    g4_df$name <- paste0("G4_", seq_len(nrow(g4_df)))
  }

  g4_df$start <- suppressWarnings(as.integer(g4_df$start))
  g4_df$end <- suppressWarnings(as.integer(g4_df$end))
  g4_df <- g4_df[!is.na(g4_df$start) & !is.na(g4_df$end), , drop = FALSE]
  if (nrow(g4_df) == 0) {
    return(data.frame())
  }

  keep <- vapply(g4_df$chrom, resolve_chromosome, character(1)) == requested_chrom
  g4_df <- g4_df[keep, , drop = FALSE]
  if (nrow(g4_df) == 0) {
    return(data.frame())
  }

  g4_df <- g4_df[g4_df$end > region_start & g4_df$start < region_end, , drop = FALSE]
  if (nrow(g4_df) == 0) {
    return(data.frame())
  }

  data.frame(
    chrom = requested_chrom,
    start = pmax(region_start, g4_df$start),
    end = pmin(region_end, g4_df$end),
    name = g4_df$name,
    feature_type = "g4_motif",
    stringsAsFactors = FALSE
  )
}

draw_g4_ruler <- function(g4_df, y, region_start, region_end) {
  g4_color <- "#008080"
  ruler_color <- adjustcolor("#4A4A4A", alpha.f = 0.65)
  tick_half_height <- 0.18

  segments(region_start, y, region_end, y, col = ruler_color, lwd = 1)
  text(
    x = region_start,
    y = y + 0.21,
    labels = "G4 motifs",
    pos = 4,
    offset = 0.2,
    cex = 0.6,
    font = 2,
    col = "black"
  )

  invisible(lapply(seq_len(nrow(g4_df)), function(i) {
    motif_start <- g4_df$start[[i]]
    motif_end <- g4_df$end[[i]]
    motif_center <- (motif_start + motif_end) / 2

    segments(motif_start, y, motif_end, y, col = g4_color, lwd = 2)
    segments(
      c(motif_start, motif_end),
      c(y - tick_half_height, y - tick_half_height),
      c(motif_start, motif_end),
      c(y + tick_half_height, y + tick_half_height),
      col = g4_color,
      lwd = 2
    )

    if ((motif_end - motif_start) < 2) {
      segments(motif_center, y - 0.22, motif_center, y + 0.22, col = g4_color, lwd = 2.2)
    }
  }))
}

render_feature_plot <- function(plot_start, plot_end, outpath, plot_chrom = chrom_input) {
  requested_chrom <- resolve_chromosome(plot_chrom)
  feature_df <- read_features_from_gff(gff_file, requested_chrom, plot_start, plot_end)
  g4_df <- read_g4_from_bed(g4_bed_file, requested_chrom, plot_start, plot_end)

  if (nrow(feature_df) == 0 && nrow(g4_df) == 0) {
    make_placeholder_plot(
      "No genes, genomic features, or G4 motifs were found in the requested W303 region.",
      plot_start,
      plot_end,
      outpath
    )
    return(invisible(NULL))
  }

  feature_order <- unname(feature_rank[feature_df$feature_type])
  feature_order[is.na(feature_order)] <- 999L

  feature_df <- feature_df[order(feature_order, feature_df$raw_start, feature_df$raw_end, feature_df$name), ]
  rownames(feature_df) <- NULL

  region_span <- max(1L, plot_end - plot_start + 1L)
  extra_rows <- if (nrow(g4_df) > 0) 1L else 0L
  total_rows <- nrow(feature_df) + extra_rows
  plot_height_inches <- max(5.9, 2.1 + 0.32 * total_rows)

  png(outpath, width = 3200, height = plot_height_inches * 200, res = 200)
  par(
    mar = c(5.8, 2.2, 3.8, 1.2),
    xaxs = "i",
    yaxs = "i"
  )

  plot(
    NA,
    xlim = c(plot_start, plot_end),
    ylim = c(0.5, total_rows + 1.2),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    bty = "n"
  )

  title(
    main = sprintf(
      "W303 NCBI annotations: %s:%s-%s",
      plot_chrom,
      plot_start,
      plot_end
    ),
    cex.main = 1.15,
    font.main = 2
  )

  axis_ticks <- pretty(c(plot_start, plot_end), n = 6)
  axis_ticks <- axis_ticks[axis_ticks >= plot_start & axis_ticks <= plot_end]
  axis(
    side = 3,
    at = axis_ticks,
    labels = sprintf("%s kb", format(round(axis_ticks / 1000), trim = TRUE)),
    cex.axis = 0.72,
    tick = TRUE,
    line = -0.4
  )

  legend_types <- unique(feature_df$feature_type)
  if (nrow(g4_df) > 0) {
    legend_types <- unique(c(legend_types, "g4_motif"))
  }
  legend_order <- unname(feature_rank[legend_types])
  legend_order[is.na(legend_order)] <- 999L
  legend_types <- legend_types[order(legend_order, legend_types)]

  abline(v = axis_ticks, col = "#EAEAEA", lwd = 0.8)

  if (nrow(feature_df) > 0) {
    feature_y_positions <- rev(seq(from = 1L + extra_rows, length.out = nrow(feature_df)))
    invisible(lapply(seq_len(nrow(feature_df)), function(i) {
      draw_feature_ruler(
        feature_df[i, , drop = FALSE],
        y = feature_y_positions[[i]],
        region_start = plot_start,
        region_end = plot_end,
        region_span = region_span
      )
    }))
  }

  if (nrow(g4_df) > 0) {
    draw_g4_ruler(
      g4_df,
      y = 1,
      region_start = plot_start,
      region_end = plot_end
    )
  }

  legend(
    "bottom",
    inset = c(0, -0.12),
    xpd = NA,
    horiz = TRUE,
    legend = vapply(legend_types, feature_type_label, character(1)),
    col = vapply(legend_types, feature_color, character(1)),
    lwd = 3,
    seg.len = 1.7,
    bty = "n",
    cex = 0.78,
    title = "NCBI Annotations Key"
  )

  dev.off()
  invisible(NULL)
}

metadata_path <- file.path(dirname(output_png), "rainplot_regions.tsv")

if (file.exists(metadata_path)) {
  region_df <- tryCatch(
    read.delim(metadata_path, stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )

  required_cols <- c("rainplot_path", "chrom", "read_start", "read_end")
  if (nrow(region_df) > 0 && all(required_cols %in% colnames(region_df))) {
    for (i in seq_len(nrow(region_df))) {
      rainplot_path <- region_df$rainplot_path[[i]]
      annotation_path <- file.path(
        dirname(rainplot_path),
        sub("^rainplot_", "annotation_", basename(rainplot_path))
      )
      plot_chrom <- first_non_empty(c(region_df$chrom[[i]], chrom_input))
      plot_start <- as.integer(region_df$read_start[[i]])
      plot_end <- as.integer(region_df$read_end[[i]])

      if (is.na(plot_start) || is.na(plot_end)) {
        next
      }
      if (plot_end <= plot_start) {
        plot_end <- plot_start + 1L
      }

      render_feature_plot(plot_start, plot_end, annotation_path, plot_chrom)

      if (i == 1L) {
        render_feature_plot(plot_start, plot_end, output_png, plot_chrom)
      }
    }

    if (file.exists(output_png) && file.info(output_png)$size > 0) {
      quit(save = "no")
    }
  }
}

render_feature_plot(region_start, region_end, output_png, chrom_input)
