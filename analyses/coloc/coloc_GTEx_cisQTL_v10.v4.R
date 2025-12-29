#!/usr/bin/env Rscript

################################################################################
# SEX-STRATIFIED eQTL COLOCALIZATION ANALYSIS - GTEx v10
################################################################################
#
# DESCRIPTION:
#   Sex-stratified colocalization analysis between GWAS summary statistics
#   (male and female) and GTEx v10 expression QTL data. Uses Bayesian 
#   colocalization to identify shared genetic signals with extended output
#   metrics including credible set sizes and lead variant information.
#
# MATCHING STRATEGY (Priority Order):
#   1. Direct GTEx variant ID match (gtex_varid)
#   2. Flipped allele variant ID match (gtex_varid_flip)
#   3. Position-only match (CHR + BP)
#
# KEY FEATURES:
#   - rsID-free analysis (uses GTEx variant IDs)
#   - Robust duplicate handling via inverse-variance weighting
#   - Streaming output writes for memory efficiency
#   - GOOD/LOW overlap quality logging
#   - Extended metrics: credible sets, top SNP details, lead eQTL identification
#
# ASSUMPTIONS:
#   - GWAS files pre-processed with GTEX_VARID and GTEX_VARID_FLIP columns
#   - GTEx variant IDs in format: CHR#_POS_A1_A2_B38
#   - All positions in GRCh38/hg38 coordinates
#   - GWAS P-values and GTEx pval_nominal available for reporting
#   - Constant sample sizes across variants
#   - Default eQTL N=1000 when allele frequency unavailable
#
# INPUT REQUIREMENTS:
#   Male GWAS:   ALL_MALE.forColoc.tsv
#   Female GWAS: ALL_FEMALE.forColoc.tsv
#   GTEx Files:  [Tissue].v10.*.parquet (eQTL data)
#                [Tissue].v10.eGenes.txt.gz (gene annotations)
#
# REQUIRED GWAS COLUMNS:
#   CHR, BP, MARKERNAME, EFFECT_ALLELE, NONEFFECT_ALLELE, P, BETA, SE,
#   AVGMAF, GTEX_VARID, GTEX_VARID_FLIP
#
# OUTPUT FILES:
#   +----------------------------------+-------------------------------------------+
#   | File Name                        | Description                               |
#   +----------------------------------+-------------------------------------------+
#   | coloc_male_full.tsv              | All male colocalization results           |
#   | coloc_female_full.tsv            | All female colocalization results         |
#   | coloc_male_strong.tsv            | Male results with PP.H4 >= 0.80           |
#   | coloc_female_strong.tsv          | Female results with PP.H4 >= 0.80         |
#   | coloc_male_low_overlap.tsv       | Male pairs with overlap < MIN_OVERLAP     |
#   | coloc_female_low_overlap.tsv     | Female pairs with overlap < MIN_OVERLAP   |
#   | locus_overlap_summary_male.tsv   | Per-locus overlap statistics (male)       |
#   | locus_overlap_summary_female.tsv | Per-locus overlap statistics (female)     |
#   | summary_report_male.txt          | Human-readable summary report (male)      |
#   | summary_report_female.txt        | Human-readable summary report (female)    |
#   | run_log_YYYYMMDD_HHMMSS.log      | Timestamped execution log                 |
#   +----------------------------------+-------------------------------------------+
#
# EXTENDED OUTPUT COLUMNS:
#   - credible_set_size: Number of SNPs with PP.H4 >= CREDIBLE_PP_THRESHOLD
#   - top_snp_id: Variant ID of SNP with highest PP.H4
#   - top_snp_pp_h4: Posterior probability of top SNP
#   - top_snp_eqtl_beta/se/p: eQTL effect estimates for top SNP
#   - top_snp_gwas_beta/se/p: GWAS effect estimates for top SNP
#   - top_snp_maf: Minor allele frequency of top SNP
#   - lead_eqtl_id/p/beta: Strongest eQTL signal in gene (by P-value)
#
# USAGE:
#   Rscript script.R                 # Run all tissues
#   Rscript script.R Brain_Cortex    # Run specific tissue
#
# DATE: 27-OCT-2025
################################################################################


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(arrow)
  library(coloc)
})

# =========================
# Quick toggle for dry-runs
# =========================
TEST_MODE        <- FALSE   # set FALSE for full run
TEST_MAX_TISSUES <- 1
TEST_MAX_LOCI    <- 1
TEST_MAX_GENES   <- 3

# -------------------------
# Command-line arg: tissue name (e.g., Brain_Cortex)
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
tissue_requested <- if (length(args) >= 1) args[1] else NA_character_

# -------------------------
# Paths
# -------------------------
base_dir      <- ".."
gtex_v10_dir  <- file.path(base_dir, "..")

# Pre-QC’d (rsID-free) inputs with gtex_varid & gtex_varid_flip
male_ss   <- file.path(base_dir, "QTLs/data/METAs/ALL_MALE.forColoc.tsv")
female_ss <- file.path(base_dir, "QTLs/data/METAs/ALL_FEMALE.forColoc.tsv")

# Tissue-specific output subfolder (safe for parallel)
tissue_suffix <- if (!is.na(tissue_requested)) paste0("/", tissue_requested) else ""
mode_suffix   <- if (TEST_MODE) "/_TEST" else ""
out_dir  <- file.path(base_dir, paste0("QTLs/analyses/coloc-GTEx-cisQTL-v10/coloc_500kb_definedloci", tissue_suffix, mode_suffix))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(out_dir, paste0("run_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

# -------------------------
# Parameters
# -------------------------
P_GWAS           <- 5e-8
WINDOW_BP        <- 5e5
MIN_OVERLAP      <- 50
PP_H4_THRESHOLD  <- 0.80
CREDIBLE_PP_THRESHOLD <- 0.01   # for credible_set_size only (doesn't affect coloc)

PRIORS           <- list(p1=1e-4, p2=1e-4, p12=1e-5)

N_FEMALE_CASES    <- 18145; N_FEMALE_CONTROLS <- 95558
N_MALE_CASES      <- 28747; N_MALE_CONTROLS   <- 83746
EQTL_N_DEFAULT    <- 1000

# -------------------------
# Helpers
# -------------------------
`%||%` <- function(a,b) if (!is.null(a) && !is.na(a)) a else b

.log <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " | ", sprintf(...))
  cat(msg, "\n"); write(msg, file=log_file, append=TRUE)
}

# Collapse duplicate SNP rows by inverse-variance weighting
.collapse_by_snp <- function(D, side = c("gwas","eqtl")) {
  side <- match.arg(side)
  D %>%
    dplyr::filter(is.finite(beta), is.finite(varbeta), varbeta > 0) %>%
    dplyr::group_by(snp) %>%
    dplyr::summarise(
      w = sum(1/varbeta, na.rm = TRUE),
      beta = sum(beta/varbeta, na.rm = TRUE) / sum(1/varbeta, na.rm = TRUE),
      varbeta = 1 / sum(1/varbeta, na.rm = TRUE),
      MAF = suppressWarnings(mean(MAF, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(beta), is.finite(varbeta), varbeta > 0)
}

# Streaming write (append with header only once)
.stream_write <- function(df, path) {
  if (!file.exists(path)) {
    data.table::fwrite(df, path, sep = "\t", quote = FALSE)
  } else {
    data.table::fwrite(df, path, sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE)
  }
}

# -------------------------
# Reading in summary statistics
# -------------------------
.read_sumstats <- function(path) {
  .log("Reading sumstats: %s", path)
  dt <- fread(path, data.table = FALSE)
  names(dt) <- toupper(names(dt))

  req <- c("CHR","BP","MARKERNAME","EFFECT_ALLELE","NONEFFECT_ALLELE","P","BETA","SE","AVGMAF","GTEX_VARID","GTEX_VARID_FLIP")
  missing <- setdiff(req, names(dt))
  if (length(missing)) stop("Missing columns in sumstats: ", paste(missing, collapse=", "))

  out <- dt %>%
    transmute(
      CHR      = as.character(CHR),
      BP       = as.integer(round(as.numeric(BP))),
      SNP      = toupper(MARKERNAME),
      A1       = toupper(EFFECT_ALLELE),
      A2       = toupper(NONEFFECT_ALLELE),
      P        = as.numeric(P),
      BETA     = as.numeric(BETA),
      SE       = as.numeric(SE),
      MAF      = as.numeric(AVGMAF),
      GTEX_VARID      = toupper(GTEX_VARID),
      GTEX_VARID_FLIP = toupper(GTEX_VARID_FLIP)
    )
  out$CHR_STR <- paste0("CHR", str_replace(out$CHR, "^CHR", ""))
  out$POS <- out$BP
  out$CHR_POS <- paste0(out$CHR_STR, "_", out$POS)
  out
}

.derive_tissue <- function(fname) sub("\\.v10\\..*$", "", fname)

.list_tissue_pairs <- function(dir) {
  files <- list.files(dir, full.names = TRUE)
  parquet <- tibble(path=files[str_detect(files,"\\.parquet$")]) %>%
    mutate(tissue=.derive_tissue(basename(path)))
  egenes  <- tibble(path=files[str_detect(files,"v10\\.eGenes\\.txt\\.gz$")]) %>%
    mutate(tissue=.derive_tissue(basename(path)))
  left_join(parquet %>% rename(parquet_path=path),
            egenes  %>% rename(egenes_path=path),
            by="tissue")
}

.read_gtex_pairs <- function(p) {
  dt <- arrow::read_parquet(p, as_data_frame = TRUE)
  names(dt) <- tolower(names(dt))
  keep <- intersect(c("variant_id","gene_id","slope","slope_se","af","tss_distance","pval_nominal"), names(dt))
  dt  <- dt[, keep, drop = FALSE]
  names(dt) <- toupper(names(dt))
  if (!"AF" %in% names(dt)) dt$AF <- NA_real_
  if ("PVAL_NOMINAL" %in% names(dt)) dt$P_EQTL <- as.numeric(dt$PVAL_NOMINAL) else dt$P_EQTL <- NA_real_

  x <- stringr::str_match(toupper(dt$VARIANT_ID),
                          "^(CHR[0-9XYMT]+)_(\\d+)_([ACGTN]+)_([ACGTN]+)_B38$")
  dt$CHR_STR <- stringr::str_replace(x[,2], "^CHR","CHR")
  dt$POS     <- as.integer(x[,3])
  dt$CHR_POS <- paste0(dt$CHR_STR, "_", dt$POS)
  dt
}

.read_gtex_egenes <- function(p) {
  eg <- data.table::fread(p, data.table = FALSE)
  names(eg) <- tolower(names(eg))
  tibble::tibble(
    GENE_ID   = toupper(eg$gene_id),
    GENE_NAME = eg$gene_name
  )
}

# -------------------------
# Locus definition (from P-value threshold ± window)
# -------------------------
.merge_intervals <- function(df) {
  df %>% arrange(CHR_STR, START, END) %>%
    group_by(CHR_STR) %>%
    mutate(grp = cumsum(START > cummax(lag(END, default=first(END))))) %>%
    group_by(CHR_STR, grp) %>%
    summarise(START=min(START), END=max(END), .groups="drop") %>%
    arrange(CHR_STR, START) %>%
    mutate(LOCUS_ID=paste0("chr", str_replace(CHR_STR,"^CHR",""), ":", START, "-", END))
}

.call_loci <- function(df, p_thresh=P_GWAS, window_bp=WINDOW_BP) {
  sig <- df %>% filter(P <= p_thresh) %>%
    transmute(CHR_STR, START=POS - window_bp, END=POS + window_bp)
  if (!nrow(sig)) return(tibble())
  .merge_intervals(sig)
}

# -------------------------
# Coloc data builders (three strategies)
# -------------------------
.prep_coloc_gwas <- function(df) {
  D <- tibble(snp = df$SNP_KEY, beta = df$BETA, varbeta = df$SE^2, MAF = df$MAF)
  D <- .collapse_by_snp(D, side = "gwas")
  list(D = D, type = "cc")
}
.prep_coloc_eqtl <- function(df, n = EQTL_N_DEFAULT) {
  D <- tibble(snp = df$SNP_KEY, beta = df$SLOPE, varbeta = df$SLOPE_SE^2, MAF = df$AF)
  D <- .collapse_by_snp(D, side = "eqtl")
  list(D = D, type = "quant", N = n)
}

.run_coloc <- function(g, e, N_cases, N_controls) {
  N <- N_cases + N_controls; s <- N_cases / N

  shared <- dplyr::inner_join(
    g$D %>% dplyr::mutate(snp_g = snp),
    e$D %>% dplyr::mutate(snp_e = snp),
    by = "snp"
  ) %>%
    dplyr::mutate(
      beta.x    = as.numeric(beta.x),
      varbeta.x = as.numeric(varbeta.x),
      beta.y    = as.numeric(beta.y),
      varbeta.y = as.numeric(varbeta.y),
      MAF.x     = as.numeric(MAF.x),
      MAF.y     = as.numeric(MAF.y)
    ) %>%
    dplyr::filter(
      is.finite(beta.x), is.finite(varbeta.x), varbeta.x > 0,
      is.finite(beta.y), is.finite(varbeta.y), varbeta.y > 0
    )

  if (!nrow(shared)) {
    return(list(res=NULL, n_overlap=0L, top_snp=NA_character_, top_pp=NA_real_))
  }

  if (any(duplicated(shared$snp))) {
    shared <- shared %>%
      dplyr::group_by(snp) %>%
      dplyr::summarise(
        # GWAS side
        w_g = sum(1/varbeta.x, na.rm = TRUE),
        beta.x = sum(beta.x/varbeta.x, na.rm = TRUE) / sum(1/varbeta.x, na.rm = TRUE),
        varbeta.x = 1 / sum(1/varbeta.x, na.rm = TRUE),
        MAF.x = suppressWarnings(mean(MAF.x, na.rm = TRUE)),
        # eQTL side
        w_e = sum(1/varbeta.y, na.rm = TRUE),
        beta.y = sum(beta.y/varbeta.y, na.rm = TRUE) / sum(1/varbeta.y, na.rm = TRUE),
        varbeta.y = 1 / sum(1/varbeta.y, na.rm = TRUE),
        MAF.y = suppressWarnings(mean(MAF.y, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::filter(
        is.finite(beta.x), is.finite(varbeta.x), varbeta.x > 0,
        is.finite(beta.y), is.finite(varbeta.y), varbeta.y > 0
      )
  }

  res <- coloc.abf(
    list(snp = shared$snp,
         beta = shared$beta.x,
         varbeta = shared$varbeta.x,
         N = N, s = s, type = g$type,
         MAF = shared$MAF.x),
    list(snp = shared$snp,
         beta = shared$beta.y,
         varbeta = shared$varbeta.y,
         N = e$N, type = e$type,
         MAF = shared$MAF.y),
    p1 = PRIORS$p1, p2 = PRIORS$p2, p12 = PRIORS$p12
  )

  top_snp <- NA_character_; top_pp <- NA_real_
  if (!is.null(res$results) && "SNP.PP.H4" %in% names(res$results)) {
    tr <- res$results %>% dplyr::arrange(dplyr::desc(SNP.PP.H4)) %>% dplyr::slice(1)
    if (nrow(tr) == 1) { top_snp <- tr$snp; top_pp <- tr$SNP.PP.H4 }
  }
  list(res = res, n_overlap = nrow(shared), top_snp = top_snp, top_pp = top_pp)
}

# orientation: "direct" (=GTEX_VARID), "flip" (=GTEX_VARID_FLIP), "pos" (=CHR+BP)
.align_and_coloc <- function(gwas_win, eqtl_gene, orientation=c("direct","flip","pos"),
                             N_cases, N_controls, min_overlap=MIN_OVERLAP) {
  orientation <- match.arg(orientation)
  if (orientation == "direct") {
    g_tbl <- gwas_win %>% filter(!is.na(GTEX_VARID), GTEX_VARID != "") %>% mutate(SNP_KEY = toupper(GTEX_VARID))
    e_tbl <- eqtl_gene %>% mutate(SNP_KEY = toupper(VARIANT_ID))
  } else if (orientation == "flip") {
    g_tbl <- gwas_win %>% filter(!is.na(GTEX_VARID_FLIP), GTEX_VARID_FLIP != "") %>% mutate(SNP_KEY = toupper(GTEX_VARID_FLIP))
    e_tbl <- eqtl_gene %>% mutate(SNP_KEY = toupper(VARIANT_ID))
  } else { # pos
    g_tbl <- gwas_win %>% mutate(SNP_KEY = toupper(CHR_POS))
    e_tbl <- eqtl_gene %>% mutate(SNP_KEY = toupper(CHR_POS))
  }

  g <- .prep_coloc_gwas(g_tbl)
  e <- .prep_coloc_eqtl(e_tbl)

  r <- .run_coloc(g, e, N_cases, N_controls)
  if (is.null(r$res) || r$n_overlap < min_overlap) {
    return(list(ok=FALSE, n_overlap = r$n_overlap, res=NULL, top_snp=NA_character_, top_pp=NA_real_))
  }

  # Attach the merged table necessary for top-snp lookups in the caller
  merged <- g_tbl %>%
    select(SNP_KEY, BETA, SE, MAF, P_GWAS = P) %>%
    inner_join(
      e_tbl %>% select(SNP_KEY, E_VARID = VARIANT_ID, SLOPE, SLOPE_SE, AF, P_EQTL),
      by = "SNP_KEY"
    )

  list(ok=TRUE, n_overlap = r$n_overlap, res=r$res,
       top_snp=r$top_snp, top_pp=r$top_pp, merged=merged)
}

# -------------------------
# Per-sex driver (single tissue optional)
# -------------------------
################################################################################
# Updated .process_sex() — reads predefined loci_500kb_topSNPs.tsv
# Includes top GWAS variant per locus and top eQTL variant per gene
################################################################################

.process_sex <- function(sex, ss_path, cases, controls) {
  .log("=== %s: START (TEST_MODE=%s) ===", toupper(sex), as.character(TEST_MODE))
  ss <- .read_sumstats(ss_path)

  # --------------------------------------------------
  # Load predefined loci (instead of re-defining from GWAS)
  # --------------------------------------------------
  loci_file <- file.path(base_dir, "QTLs/analyses/define_loci/loci_500kb_topSNPs.tsv")
  .log("%s: Reading predefined loci from %s", sex, loci_file)
  loci <- fread(loci_file)
  names(loci) <- toupper(names(loci))
  loci <- loci %>%
    transmute(
      CHR_STR = toupper(CHR),
      START   = as.integer(LEFT_500KB),
      END     = as.integer(RIGHT_500KB),
      LOCUS_ID = LOCUSID
    )
  if (TEST_MODE && nrow(loci) > TEST_MAX_LOCI) loci <- loci[seq_len(TEST_MAX_LOCI), , drop=FALSE]
  .log("%s: %d loci loaded from predefined file%s",
       sex, nrow(loci), if (TEST_MODE) sprintf(" [capped to %d]", nrow(loci)) else "")

  # --------------------------------------------------
  # GTEx tissues
  # --------------------------------------------------
  pairs <- .list_tissue_pairs(gtex_v10_dir)
  if (!nrow(pairs)) stop("No GTEx tissues found in ", gtex_v10_dir)
  if (!is.na(tissue_requested)) {
    pairs <- pairs %>% filter(tissue == tissue_requested)
    if (!nrow(pairs)) stop("Requested tissue not found: ", tissue_requested)
    .log("Filtering to requested tissue: %s", tissue_requested)
  }
  if (TEST_MODE && nrow(pairs) > TEST_MAX_TISSUES) pairs <- pairs[seq_len(TEST_MAX_TISSUES), , drop=FALSE]

  # --------------------------------------------------
  # Tracking and output setup
  # --------------------------------------------------
  strategy_totals <- c(direct=0L, flip=0L, pos=0L)
  locus_rows <- list()
  bucket_counts <- c(good=0L, marginal=0L, low=0L)

  full_path <- file.path(out_dir, paste0("coloc_", tolower(sex), "_full.tsv"))
  low_path  <- file.path(out_dir, paste0("coloc_", tolower(sex), "_low_overlap.tsv"))
  if (file.exists(full_path)) file.remove(full_path)
  if (file.exists(low_path))  file.remove(low_path)

  # --------------------------------------------------
  # Tissue-level iteration
  # --------------------------------------------------
  for (k in seq_len(nrow(pairs))) {
    tn <- pairs$tissue[k]
    eq_pairs <- .read_gtex_pairs(pairs$parquet_path[k])
    eg_info  <- .read_gtex_egenes(pairs$egenes_path[k])
    eq_pairs <- eq_pairs %>% left_join(eg_info, by="GENE_ID")

    for (i in seq_len(nrow(loci))) {
      r <- loci[i,]
      gwas_win <- ss %>% filter(CHR_STR==r$CHR_STR, POS>=r$START, POS<=r$END)
      eqtl_win <- eq_pairs %>% filter(CHR_STR==r$CHR_STR, POS>=r$START, POS<=r$END)

      n_gwas      <- nrow(gwas_win)
      n_eqtl      <- nrow(eqtl_win)
      n_gtexid    <- sum(!is.na(gwas_win$GTEX_VARID) & gwas_win$GTEX_VARID!="")
      n_gtexid_fl <- sum(!is.na(gwas_win$GTEX_VARID_FLIP) & gwas_win$GTEX_VARID_FLIP!="")

      .log("%s | %s | %s: GWAS=%s (gtex_varid=%s; gtex_varid_flip=%s), eQTL=%s",
           sex, tn, r$LOCUS_ID,
           format(n_gwas, big.mark=","), format(n_gtexid, big.mark=","), format(n_gtexid_fl, big.mark=","),
           format(n_eqtl, big.mark=","))

      if (!n_gwas || !n_eqtl) next

      gene_ids <- unique(eqtl_win$GENE_ID)
      if (TEST_MODE && length(gene_ids) > TEST_MAX_GENES) gene_ids <- gene_ids[seq_len(TEST_MAX_GENES)]

      best_for_locus <- list(overlap=0L, pct=0, matched=0L, locus_id=r$LOCUS_ID)

      # --------------------------------------------------
      # Top GWAS variant by P-value within this locus
      # --------------------------------------------------
      top_gwas_in_locus <- gwas_win %>%
        arrange(P) %>% slice(1) %>%
        transmute(
          top_gwas_variant = SNP,
          top_gwas_pval    = P,
          top_gwas_beta    = BETA,
          top_gwas_se      = SE
        )

      for (gid in gene_ids) {
        egene <- eqtl_win %>% filter(GENE_ID == gid)
        gname <- egene$GENE_NAME[1] %||% gid

        # Try all 3 match orientations
        r_direct <- .align_and_coloc(gwas_win, egene, "direct", N_cases=cases, N_controls=controls)
        r_flip   <- .align_and_coloc(gwas_win, egene, "flip",   N_cases=cases, N_controls=controls)
        r_pos    <- .align_and_coloc(gwas_win, egene, "pos",    N_cases=cases, N_controls=controls)

        ov <- c(direct=r_direct$n_overlap, flip=r_flip$n_overlap, pos=r_pos$n_overlap)
        ord <- c("direct","flip","pos")
        best_key <- ord[order(-ov[ord], match(ord, ord))][1]
        best_res <- switch(best_key, direct=r_direct, flip=r_flip, pos=r_pos)

        if (!best_res$ok) {
          low_row <- tibble(
            sex=sex, tissue=tn, locus_id=r$LOCUS_ID, chr=r$CHR_STR, start=r$START, end=r$END,
            gene_id=gid, gene_name=gname,
            gwas_in_window=n_gwas, eqtl_in_window=n_eqtl,
            gwas_with_gtex_varid=n_gtexid, gwas_with_gtex_varid_flip=n_gtexid_fl,
            overlap_direct=ov["direct"], overlap_flip=ov["flip"], overlap_pos=ov["pos"]
          )
          .stream_write(low_row, low_path)
          next
        }

        # --------------------------------------------------
        # Extract posterior probabilities and top coloc SNP
        # --------------------------------------------------
        s <- best_res$res$summary
        PP.H0 <- as.numeric(s["PP.H0.abf"])
        PP.H1 <- as.numeric(s["PP.H1.abf"])
        PP.H2 <- as.numeric(s["PP.H2.abf"])
        PP.H3 <- as.numeric(s["PP.H3.abf"])
        PP.H4 <- as.numeric(s["PP.H4.abf"])

        cred_size <- 0L
        top_snp_id <- NA_character_; top_snp_pp <- NA_real_
        top_eqtl_beta <- NA_real_; top_eqtl_se <- NA_real_; top_eqtl_p <- NA_real_
        top_gwas_beta <- NA_real_; top_gwas_se <- NA_real_; top_gwas_p <- NA_real_; top_maf <- NA_real_

        if (!is.null(best_res$res$results)) {
          rr <- as_tibble(best_res$res$results)
          if (nrow(rr)) {
            rr <- rr %>% arrange(desc(SNP.PP.H4))
            top_row <- rr %>% slice(1)
            if (nrow(top_row)) {
              top_key <- as.character(top_row$snp[1])
              top_snp_pp <- as.numeric(top_row$SNP.PP.H4[1])
              cred_size <- sum(rr$SNP.PP.H4 >= CREDIBLE_PP_THRESHOLD, na.rm = TRUE)

              if (!is.null(best_res$merged) && nrow(best_res$merged)) {
                M <- best_res$merged %>% filter(SNP_KEY == top_key)
                if (nrow(M)) {
                  top_snp_id   <- toupper(M$E_VARID[1])
                  top_eqtl_beta <- as.numeric(M$SLOPE[1])
                  top_eqtl_se   <- as.numeric(M$SLOPE_SE[1])
                  top_eqtl_p    <- as.numeric(M$P_EQTL[1])
                  top_gwas_beta <- as.numeric(M$BETA[1])
                  top_gwas_se   <- as.numeric(M$SE[1])
                  top_gwas_p    <- as.numeric(M$P_GWAS[1])
                  top_maf       <- as.numeric(M$MAF[1])
                }
              }
            }
          }
        }

        # --------------------------------------------------
        # Lead eQTL variant (lowest P-value)
        # --------------------------------------------------
        lead_eqtl_id <- NA_character_; lead_eqtl_p <- NA_real_; lead_eqtl_beta <- NA_real_
        if ("P_EQTL" %in% names(egene) && any(is.finite(egene$P_EQTL))) {
          eg_lead <- egene %>% arrange(P_EQTL) %>% filter(is.finite(P_EQTL)) %>% slice(1)
          if (nrow(eg_lead)) {
            lead_eqtl_id   <- toupper(eg_lead$VARIANT_ID[1])
            lead_eqtl_p    <- as.numeric(eg_lead$P_EQTL[1])
            lead_eqtl_beta <- as.numeric(eg_lead$SLOPE[1])
          }
        }

        # --------------------------------------------------
        # Final row (includes top GWAS and top eQTL variants)
        # --------------------------------------------------
        row <- tibble(
          sex=sex, tissue=tn, locus_id=r$LOCUS_ID, chr=r$CHR_STR, start=r$START, end=r$END,
          gene_id=gid, gene_name=gname, orientation_used=best_key,
          gwas_in_window=n_gwas, eqtl_in_window=n_eqtl,
          gwas_with_gtex_varid=n_gtexid, gwas_with_gtex_varid_flip=n_gtexid_fl,
          overlap_direct=ov["direct"], overlap_flip=ov["flip"], overlap_pos=ov["pos"],
          n_snps_overlap=best_res$n_overlap,
          PP.H0=PP.H0, PP.H1=PP.H1, PP.H2=PP.H2, PP.H3=PP.H3, PP.H4=PP.H4,
          credible_set_size = cred_size,
          top_snp_id        = top_snp_id,
          top_snp_pp_h4     = top_snp_pp,
          top_snp_eqtl_beta = top_eqtl_beta,
          top_snp_eqtl_se   = top_eqtl_se,
          top_snp_eqtl_p    = top_eqtl_p,
          top_snp_gwas_beta = top_gwas_beta,
          top_snp_gwas_se   = top_gwas_se,
          top_snp_gwas_p    = top_gwas_p,
          top_snp_maf       = top_maf,
          lead_eqtl_id      = lead_eqtl_id,
          lead_eqtl_p       = lead_eqtl_p,
          lead_eqtl_beta    = lead_eqtl_beta
        ) %>%
        bind_cols(top_gwas_in_locus)

        .stream_write(row, full_path)

        if (best_res$n_overlap >= 50L) bucket_counts["good"] <- bucket_counts["good"] + 1L
        else if (best_res$n_overlap >= 20L) bucket_counts["marginal"] <- bucket_counts["marginal"] + 1L
        else bucket_counts["low"] <- bucket_counts["low"] + 1L

        if (best_res$n_overlap > best_for_locus$overlap) {
          best_for_locus$overlap <- best_res$n_overlap
          best_for_locus$pct     <- if (n_gwas) 100*best_res$n_overlap/n_gwas else 0
          best_for_locus$matched <- best_res$n_overlap
        }
      } # end gene loop

      # Save locus-level summary
      if (best_for_locus$overlap >= 0) {
        locus_rows[[length(locus_rows)+1]] <- tibble(
          locus_id = r$LOCUS_ID,
          n_gwas   = n_gwas,
          n_eqtl   = n_eqtl,
          match_combined = best_for_locus$matched,
          pct_gwas_matched = best_for_locus$pct
        )
      }
    } # end locus loop
  } # end tissue loop

  # --------------------------------------------------
  # Summary outputs (unchanged)
  # --------------------------------------------------
  full_df <- if (file.exists(full_path)) fread(full_path) else data.table()
  strong <- if (nrow(full_df)) {
    as_tibble(full_df) %>% filter(!is.na(PP.H4), PP.H4 >= PP_H4_THRESHOLD, n_snps_overlap >= MIN_OVERLAP)
  } else tibble()
  write_tsv(strong, file.path(out_dir, paste0("coloc_", tolower(sex), "_strong.tsv")))

  locus_tbl <- if (length(locus_rows)) bind_rows(locus_rows) else tibble()
  write_tsv(locus_tbl %>% arrange(locus_id),
            file.path(out_dir, paste0("locus_overlap_summary_", tolower(sex), ".tsv")))

  report_path <- file.path(out_dir, paste0("summary_report_", tolower(sex), ".txt"))
  .log("%s: Completed. Strong=%s | Full=%s | Report=%s",
       sex, nrow(strong), nrow(full_df), basename(report_path))

  invisible(list(full=full_df, strong=strong, locus=locus_tbl, report=report_path))
}


# -------------------------
# Run both sexes (for requested tissue only)
# -------------------------
start_time <- Sys.time(); .log("Run start. TEST_MODE=%s", as.character(TEST_MODE))
male   <- .process_sex("Male",   male_ss,   N_MALE_CASES,   N_MALE_CONTROLS)
female <- .process_sex("Female", female_ss, N_FEMALE_CASES, N_FEMALE_CONTROLS)
end_time <- Sys.time()
.log("Run end. Duration: %.1f mins", as.numeric(difftime(end_time, start_time, units="mins")))
