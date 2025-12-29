library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(RMySQL)
library(DBI)
library(locuscomparer)

## Ran locally for handful of SNPs 
## Set axes so the right hand plots have same x and y 

chr      <- "14"
bp       <- "87949673"
ref      <- "G"
alt      <- "A"
lead_snp <- "rs3213916"

base_path <- "../locuscompare/"
out_dir   <- "../locuscompare/output/"

read_metal <- function(in_fn, marker_col = "rsid", pval_col = "pval") {
  if (is.character(in_fn)) {
    d <- read.table(in_fn, header = TRUE, stringsAsFactors = FALSE)
    colnames(d)[colnames(d) == marker_col] <- "rsid"
    colnames(d)[colnames(d) == pval_col]   <- "pval"
  } else if (is.data.frame(in_fn)) {
    d <- in_fn
  } else {
    stop('in_fn must be a string or data.frame')
  }
  d$logp <- -log10(d$pval)
  d[, c("rsid", "pval", "logp")]
}

get_position <- function(x, genome = c("hg19","hg38")) {
  genome <- match.arg(genome)
  data(config)
  on.exit(rm(config))
  conn <- RMySQL::dbConnect(RMySQL::MySQL(), "locuscompare", config$b, config$c, config$a)
  on.exit(RMySQL::dbDisconnect(conn), add = TRUE)
  stopifnot("rsid" %in% colnames(x))
  cmd <- sprintf(
    "select rsid, chr, pos from tkg_p3v5a_%s where rsid in ('%s')",
    genome, paste0(x$rsid, collapse = "','")
  )
  res <- DBI::dbGetQuery(conn = conn, statement = cmd)
  merge(x, res, by = "rsid")
}

retrieve_LD <- function(chr, snp, population) {
  data(config)
  on.exit(rm(config))
  conn <- RMySQL::dbConnect(RMySQL::MySQL(), "locuscompare", config$b, config$c, config$a)
  on.exit(RMySQL::dbDisconnect(conn), add = TRUE)
  
  res1 <- DBI::dbGetQuery(
    conn = conn,
    statement = sprintf(
      "select SNP_A, SNP_B, R2
       from tkg_p3v5a_ld_chr%s_%s
       where SNP_A = '%s';",
      chr, population, snp
    )
  )
  
  res2 <- DBI::dbGetQuery(
    conn = conn,
    statement = sprintf(
      "select SNP_B as SNP_A, SNP_A as SNP_B, R2
       from tkg_p3v5a_ld_chr%s_%s
       where SNP_B = '%s';",
      chr, population, snp
    )
  )
  
  rbind(res1, res2)
}

get_lead_snp <- function(merged, snp = NULL) {
  if (is.null(snp)) {
    snp <- merged[which.min(merged$pval1 + merged$pval2), "rsid"]
  } else if (!snp %in% merged$rsid) {
    stop(sprintf("%s not found in intersection of in_fn1 and in_fn2.", snp))
  }
  as.character(snp)
}

assign_color <- function(rsid, snp, ld) {
  ld <- ld[ld$SNP_A == snp, ]
  ld$color <- as.character(cut(
    ld$R2,
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c("blue4", "skyblue", "darkgreen", "orange", "red"),
    include.lowest = TRUE
  ))
  
  color <- data.frame(rsid, stringsAsFactors = FALSE)
  color <- merge(color, ld[, c("SNP_B", "color")],
                 by.x = "rsid", by.y = "SNP_B", all.x = TRUE)
  color[is.na(color$color), "color"] <- "blue4"
  
  if (snp %in% color$rsid) {
    color[rsid == snp, "color"] <- "purple"
  } else {
    color <- rbind(color, data.frame(rsid = snp, color = "purple"))
  }
  
  res <- color$color
  names(res) <- color$rsid
  res
}

add_label <- function(merged, snp) {
  merged$label <- ifelse(merged$rsid %in% snp, merged$rsid, "")
  merged
}

make_scatterplot <- function(merged, title1, title2, color, shape, size) {
  ggplot(merged, aes(logp1, logp2)) +
    geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) +
    geom_point(
      data = merged[merged$label != "", ],
      aes(logp1, logp2, fill = rsid, size = rsid, shape = rsid)
    ) +
    xlab(bquote(.(title1) ~ -log[10] * "(P)")) +
    ylab(bquote(.(title2) ~ -log[10] * "(P)")) +
    scale_fill_manual(values = color, guide = "none") +
    scale_shape_manual(values = shape, guide = "none") +
    scale_size_manual(values = size, guide = "none") +
    ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf) +
    theme_classic()
}

make_locuszoom <- function(metal, title, chr, color, shape, size) {
  ggplot(metal, aes(x = pos, logp)) +
    geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) +
    geom_point(
      data = metal[metal$label != "", ],
      aes(x = pos, logp, fill = rsid, size = rsid, shape = rsid)
    ) +
    scale_fill_manual(values = color, guide = "none") +
    scale_shape_manual(values = shape, guide = "none") +
    scale_size_manual(values = size, guide = "none") +
    scale_x_continuous(labels = function(x) sprintf("%.1f", x / 1e6)) +
    ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf) +
    xlab(paste0("chr", chr, " (Mb)")) +
    ylab(bquote(.(title) ~ -log[10] * "(P)")) +
    theme_classic() +
    theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "lines"))
}

make_combined_plot <- function(merged, title1, title2, ld, chr, snp = NULL) {
  snp   <- get_lead_snp(merged, snp)
  color <- assign_color(merged$rsid, snp, ld)
  
  shape <- ifelse(merged$rsid == snp, 23, 21)
  names(shape) <- merged$rsid
  
  size <- ifelse(merged$rsid == snp, 3, 2)
  names(size) <- merged$rsid
  
  merged <- add_label(merged, snp)
  
  max_logp <- max(c(merged$logp1, merged$logp2), na.rm = TRUE)
  if (!is.finite(max_logp) || max_logp <= 0) max_logp <- 1
  max_logp <- ceiling(max_logp)
  
  p1 <- make_scatterplot(merged, title1, title2, color, shape, size) +
    scale_x_continuous(limits = c(0, max_logp)) +
    scale_y_continuous(limits = c(0, max_logp))
  
  metal1 <- merged[, c("rsid", "logp1", "chr", "pos", "label")]
  colnames(metal1)[colnames(metal1) == "logp1"] <- "logp"
  p2 <- make_locuszoom(metal1, title1, chr, color, shape, size) +
    scale_y_continuous(limits = c(0, max_logp))
  
  metal2 <- merged[, c("rsid", "logp2", "chr", "pos", "label")]
  colnames(metal2)[colnames(metal2) == "logp2"] <- "logp"
  p3 <- make_locuszoom(metal2, title2, chr, color, shape, size) +
    scale_y_continuous(limits = c(0, max_logp))
  
  p2 <- p2 + theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank())
  
  right <- cowplot::plot_grid(p2, p3, align = "v", nrow = 2,
                              rel_heights = c(0.8, 1))
  
  cowplot::plot_grid(
    p1, right,
    align      = "h",
    rel_widths = c(1.6, 2)
  )
}

locuscompare <- function(in_fn1, in_fn2,
                         title1 = "eQTL", title2 = "GWAS",
                         snp = NULL, population = "EUR",
                         genome = c("hg19","hg38")) {
  
  d1 <- read_metal(in_fn1)
  d2 <- read_metal(in_fn2)
  
  merged <- merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
  genome <- match.arg(genome)
  merged <- get_position(merged, genome)
  
  chr <- unique(merged$chr)
  if (length(chr) != 1) stop("There must be one and only one chromosome.")
  
  snp <- get_lead_snp(merged, snp)
  ld  <- retrieve_LD(chr, snp, population)
  make_combined_plot(merged, title1, title2, ld, chr, snp)
}


gwas_file <- paste0(base_path, "merged_chr", chr, "_", bp, "_", ref, "_", alt, "_female.txt")
eqtl_file <- paste0(base_path, "merged_chr", chr, "_", bp, "_", ref, "_", alt, ".txt")

gwas_fn <- fread(gwas_file)
gwas_fn <- gwas_fn[gwas_fn$rsid != "."]

eqtl_fn <- fread(eqtl_file)
eqtl_fn <- eqtl_fn[eqtl_fn$rsid != "."]

p <- locuscompare(
  in_fn1 = gwas_fn,
  in_fn2 = eqtl_fn,
  title1 = "FEMALE",
  title2 = "MALE",
  snp    = lead_snp,
  population = "EUR",
  genome     = "hg38"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
base_name <- paste0("locuscompare_chr", chr, "_", bp, "_", ref, "_", alt)

ggsave(
  filename = file.path(out_dir, paste0(base_name, ".png")),
  plot = p,
  width = 12,
  height = 6,
  dpi = 300
)
