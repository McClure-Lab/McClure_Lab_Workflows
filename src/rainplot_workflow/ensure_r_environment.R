#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: ensure_r_environment.R <project_dir>", call. = FALSE)
}

project_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
setwd(project_dir)

options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(
  RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
  RENV_CONSENT = "yes"
)

required_packages <- c("BiocManager", "Gviz", "rtracklayer", "GenomicRanges")
lockfile_path <- file.path(project_dir, "renv.lock")
activate_path <- file.path(project_dir, "renv", "activate.R")

if (!requireNamespace("renv", quietly = TRUE)) {
  message("[INFO] Installing renv...")
  install.packages("renv", repos = "https://cloud.r-project.org")
}

renv::consent(provided = TRUE)

if (!file.exists(activate_path)) {
  message("[INFO] Initializing renv project...")
  renv::init(bare = TRUE, restart = FALSE)
} else {
  message("[INFO] renv project already initialized.")
}

if (file.exists(lockfile_path)) {
  message("[INFO] Restoring R packages from renv.lock...")
  renv::restore(prompt = FALSE)
} else {
  message("[INFO] No renv.lock found yet. A new one will be created after package installation.")
}

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  message(sprintf(
    "[INFO] Installing missing renv packages: %s",
    paste(missing_packages, collapse = ", ")
  ))
  renv::install(missing_packages)
} else {
  message("[INFO] Required R packages already available in renv.")
}

message("[INFO] Snapshotting renv.lock...")
renv::snapshot(prompt = FALSE)

message(sprintf("[INFO] renv ready in: %s", project_dir))
