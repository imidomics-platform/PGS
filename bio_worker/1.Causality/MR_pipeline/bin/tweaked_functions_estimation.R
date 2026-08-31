my_ld_matrix_local <- function (variants, refPanel, bfile, plink_bin){
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", 
                  "sh")
  message("p1")
  fn <- tempfile()
  fwrite(data.frame(variants), file = fn, row.names = F, 
         col.names = F, quote = F)
  bim <- refPanel[match(variants, refPanel$V2),]
  fun2 <- paste0(shQuote(plink_bin, type = shell), 
                 " --bfile ", shQuote(bfile, type = shell), 
                 " --extract ", shQuote(fn, type = shell), 
                 " --r square ", 
                 " --keep-allele-order ", 
                 " --out ", shQuote(fn, type = shell))
  system(fun2, ignore.stdout = T, ignore.stderr = T)
  res <- fread(paste0(fn, ".ld"), header = FALSE)
  system(paste0("rm ", fn, "*"))
  
  if(length(colnames(res))==length(rownames(res)) & length(colnames(res)) == length(paste(bim$V2, bim$V5, bim$V6, sep = "_"))){
    res <- as.matrix(res)
    rownames(res) <- colnames(res) <- paste(bim$V2, bim$V5, bim$V6, sep = "_")
    cat("Returning ok")
    return(res)
  } else{
    cat("Dimension disagreement: ", dim(res), " vs ", length(length(paste(bim$V2, bim$V5, bim$V6, sep = "_"))), "\n")
    cat("Duplicated variants? ", any(duplicated(variants)), "\n")
    return(NULL)
  }
}

ld_matrix <- function (variants, with_alleles = TRUE, pop = "EUR", bfile = NULL, plink_bin = NULL){
  if (length(variants) > 500 & is.null(bfile)) {
    stop("SNP list must be smaller than 500. Try running locally by providing local ld reference with bfile argument. See vignettes for a guide on how to do this.")
  }
  if (is.null(bfile)) {
    message("Please look at vignettes for options on running this locally if you need to run many instances of this command.")
  }
  if (!is.null(bfile)) {
    message("p2")
    return(my_ld_matrix_local(variants, refPanel = refPanel, bfile = bfile, plink_bin = plink_bin))
  }
  res <- api_query("ld/matrix", query = list(rsid = variants, 
                                             pop = pop), access_token = NULL) %>% get_query_content()
  if (all(is.na(res))) 
    stop("None of the requested variants were found")
  variants2 <- res$snplist
  res <- res$matrix
  res <- matrix(as.numeric(res), nrow(res), ncol(res))
  variants3 <- do.call(rbind, strsplit(variants2, split = "_"))
  if (with_alleles) {
    rownames(res) <- variants2
    colnames(res) <- variants2
  }
  else {
    rownames(res) <- variants3[, 1]
    colnames(res) <- variants3[, 1]
  }
  missing <- variants[!variants %in% variants3[, 1]]
  if (length(missing) > 0) {
    warning("The following variants are not present in the LD reference panel\n", 
            paste(missing, collapse = "\n"))
  }
  ord <- match(variants3[, 1], variants)
  res <- res[order(ord), order(ord)]
  return(res)
}

dat_to_MRInput <- function(dat, get_correlations=FALSE, pop="EUR"){
  out <- plyr::dlply(dat, c("exposure", "outcome"), function(x)
  {
    x <- plyr::mutate(x)
    message("Converting:")
    message(" - exposure: ", x$exposure[1])
    message(" - outcome: ", x$outcome[1])
    if(get_correlations & length(unique(x$SNP[x$SNP %in% refPanel$V2]))>1)
    {
      key <- paste0(gene, "_", is)
      ld_file <- paste0(ld_dir, "/", key, ".rds")
      if(file.exists(ld_file)){
        message("p4")
        message(" - reading LD matrix")
        ld_obj <- readRDS(ld_file)
        snps <- x$SNP
        idx <- match(snps, ld_obj$snps)
        ld <- ld_obj$ld[idx, idx]
        bim <- refPanel[match(snps, refPanel$V2),]
        rownames(ld) <- colnames(ld) <- paste(bim$V2, bim$V5, bim$V6, sep = "_")
      } else{
        message("p3")
        message(" - obtaining LD matrix")
        ld <- ld_matrix(unique(x$SNP), plink_bin = plink_bin, bfile = plink_ref)
      }
      
      # Rearrange so that order is same as in x!
      snpnames <- do.call(rbind, strsplit(rownames(ld), split = "_"))[,1]
      x <- x[match(snpnames, x$SNP),]
      # ld <- ld_matrix(unique(x$SNP), pop=pop)
      out <- harmonise_ld_dat(x, ld)
      if(is.null(out))
      {
        return(NULL)
      }
      x <- out$x
      ld <- out$ld
      
      if(is.matrix(ld)){
        MendelianRandomization::mr_input(
          bx = x$beta.exposure,
          bxse = x$se.exposure,
          by = x$beta.outcome,
          byse = x$se.outcome,
          exposure = x$exposure[1],
          outcome = x$outcome[1],
          snps = x$SNP,
          effect_allele=x$effect_allele.exposure,
          other_allele=x$other_allele.exposure,
          eaf = x$eaf.exposure,
          correlation = ld
        )
      } else{
        MendelianRandomization::mr_input(
          bx = x$beta.exposure,
          bxse = x$se.exposure,
          by = x$beta.outcome,
          byse = x$se.outcome,
          exposure = x$exposure[1],
          outcome = x$outcome[1],
          snps = x$SNP,
          effect_allele=x$effect_allele.exposure,
          other_allele=x$other_allele.exposure,
          eaf = x$eaf.exposure
        )
      }
      
    } else {
      MendelianRandomization::mr_input(
        bx = x$beta.exposure,
        bxse = x$se.exposure,
        by = x$beta.outcome,
        byse = x$se.outcome,
        exposure = x$exposure[1],
        outcome = x$outcome[1],
        snps = x$SNP,
        effect_allele=x$effect_allele.exposure,
        other_allele=x$other_allele.exposure,
        eaf = x$eaf.exposure
      )
    }
  })
  return(out)
}

# mr_raps <- function (b_exp, b_out, se_exp, se_out, parameters = default_parameters()) {
#   cpg <- requireNamespace("mr.raps", quietly = TRUE)
#   if (!cpg) {
#     stop("Please install the mr.raps package using devtools::install_github('qingyuanzhao/mr.raps')")
#   }
#   data <- data.frame(beta.exposure = b_exp, beta.outcome = b_out, 
#                      se.exposure = se_exp, se.outcome = se_out)
#   out <- suppressMessages(mr.raps::mr.raps(b_exp = b_exp, b_out = b_out, se_exp = se_exp, se_out = se_out, diagnosis = FALSE))
#   list(b = out$beta.hat, se = out$beta.se, pval = pnorm(-abs(out$beta.hat/out$beta.se)) * 
#          2, nsnp = length(b_exp))
# }
# 
# assignInNamespace("mr_raps", mr_raps, ns = "TwoSampleMR")

fit.mixture.model <- function (z, n = 2, ntry = 20, force.mu.zero = TRUE, diagnostics = FALSE) 
{
  loglike <- function(param, z) {
    n <- length(param)/3
    p <- param[1:n]
    p <- p/sum(p)
    mu <- param[n + 1:n]
    sigma <- param[2 * n + 1:n]
    l <- matrix(0, length(z), length(p))
    for (i in 1:length(p)) {
      l[, i] <- p[i] * dnorm(z, mu[i], sqrt(sigma[i]^2 + 
                                              1))
    }
    row_sums <- rowSums(l)
    row_sums[row_sums < 1e-300] <- 1e-300
    -sum(log(row_sums))
  }
  get.random.init <- function(n = 2) {
    p <- rgamma(n, 1)
    p <- p/sum(p)
    mu <- rnorm(n)
    sigma <- c(0.1, rexp(n - 1) * 2)
    c(p, mu, sigma)
  }
  if (force.mu.zero) {
    mu.lower <- -1e-08
    mu.upper <- 1e-08
  }
  else {
    mu.lower <- -Inf
    mu.upper <- Inf
  }
  res <- list()
  for (try in 1:ntry) {
    try(res[[try]] <- optim(get.random.init(n), function(param) loglike(param, 
                                                                        z), method = "L-BFGS-B", lower = c(rep(0.01, n), 
                                                                                                           rep(mu.lower, n), rep(0.01, n)), upper = c(rep(0.99, 
                                                                                                                                                          n), rep(mu.upper, n), rep(Inf, n))), silent = TRUE)
  }
  i <- which.min(sapply(1:length(res), function(i) {
    tmp <- res[[i]]$value
    if (is.null(tmp)) {
      Inf
    }
    else {
      tmp
    }
  }))
  if (length(i) == 0) {
    stop("Failed to initialize, increase ntry")
  }
  res <- res[[i]]
  param <- res$par
  p <- param[1:n]
  p <- p/sum(p)
  mu <- param[n + 1:n]
  mu[abs(mu) < 1e-06] <- 0
  sigma <- param[2 * n + 1:n]
  if (diagnostics) {
    print("Estimated mixture model:")
    print(paste0("p = ", signif(p, 3), ", ", "mu = ", signif(mu, 
                                                             3), ", ", "sigma = ", signif(sigma, 3)))
    print("Generating diagnostic plots...")
    t <- seq(-10, 10, 0.1)
    fitted.cdf <- sapply(1:n, function(i) p[i] * pnorm(t, 
                                                       mu[i], sqrt(sigma[i]^2 + 1)))
    fitted.cdf <- rowSums(fitted.cdf)
    emp.cdf <- sapply(t, function(t) mean(z < t))
    par(mfrow = c(1, 2))
    plot(t, fitted.cdf, col = "red")
    lines(t, emp.cdf)
    plot(fitted.cdf, emp.cdf)
    abline(0, 1)
  }
  list(p = p, mu = mu, sigma = sigma)
}
assignInNamespace("fit.mixture.model", fit.mixture.model, ns = "mr.raps")

# b_exp <- x$beta.exposure
# b_out <- x$beta.outcome
# se_exp <- x$se.exposure
# se_out <- x$se.outcome
# 
#   if (!(requireNamespace("mr.raps", quietly = TRUE))) {
#     stop("You can install mr.raps with install.packages('mr.raps', repos = c('https://mrcieu.r-universe.dev', 'https://cloud.r-project.org'))")
#   }
#   data <- data.frame(beta.exposure = b_exp, beta.outcome = b_out, 
#                      se.exposure = se_exp, se.outcome = se_out)
#   out <- suppressMessages(mr.raps::mr.raps(data, diagnostics = FALSE, 
#                                            over.dispersion = parameters$over.dispersion, loss.function = parameters$loss.function, 
#                                            shrinkage = parameters$shrinkage))
#   list(b = out$beta.hat, se = out$beta.se, pval = stats::pnorm(-abs(out$beta.hat/out$beta.se)) * 
#          2, nsnp = length(b_exp))

