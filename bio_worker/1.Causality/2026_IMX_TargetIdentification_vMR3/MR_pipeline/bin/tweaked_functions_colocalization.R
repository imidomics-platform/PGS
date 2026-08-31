my_ld_matrix_local <- function (variants, refPanel, bfile, plink_bin){
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", 
                  "sh")
  fn <- tempfile()
  fwrite(data.frame(variants), file = fn, row.names = F, 
         col.names = F, quote = F)
  bim <- refPanel[match(variants, refPanel$V2),]
  fun2 <- paste0(shQuote(plink_bin, type = shell), " --bfile ", 
                 shQuote(bfile, type = shell), " --extract ", shQuote(fn, 
                                                                      type = shell), " --r square ", " --keep-allele-order ", 
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

read_outcome <- function(study, rsids, chrom.pos){
  
  if(study=="RA_Ha2020"){
    out.dat <- ha[ha$SNP %in% chrom.pos,]
    out.dat$SNP <- rsids[match(out.dat$SNP, chrom.pos)]
  }
  if(study=="RA_Ishi2022"){
    out.dat <- ishi[ishi$SNP %in% chrom.pos,]
    out.dat$SNP <- rsids[match(out.dat$SNP, chrom.pos)]
  }
  if(study=="RA_Okada2014"){
    out.dat <- extract_outcome_data(snps=rsids, outcomes="ieu-a-833")
  }
  if(study=="PS_Stuart2022"){
    out.dat <- ps.stuart[ps.stuart$SNP %in% chrom.pos,]
    out.dat$SNP <- rsids[match(out.dat$SNP, chrom.pos)]
  }
  if(study=="CD_deLange2017"){
    out.dat <- extract_outcome_data(snps=rsids, outcomes="ebi-a-GCST004132")
  }
  if(study=="UC_deLange2017"){
    out.dat <- extract_outcome_data(snps=rsids, outcomes="ebi-a-GCST004133")
  }
  if(study=="SLE_Bentham2015"){
    out.dat <- extract_outcome_data(snps=rsids, outcomes="ebi-a-GCST003156")
  }
  if(study=="PSA_Stuart2015"){
    out.dat <- psa.stuart[psa.stuart$SNP %in% chrom.pos,]
    out.dat$SNP <- rsids[match(out.dat$SNP, chrom.pos)]
  } 
  if(study=="PSA_Finn2021"){
    out.dat <- extract_outcome_data(snps=rsids, outcomes="finn-b-L12_PSORI_ARTHRO", opengwas_jwt = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJzZXJnaW8ubWFydGluZXpAaW1pZG9taWNzLmNvbSIsImlhdCI6MTc1MDIzNTQ2MywiZXhwIjoxNzUxNDQ1MDYzfQ.QJdkY8J_sOCZyddcf3t9dkuzJbi0ltselrtUGtZoOadO5vJvyWb-ZXpuWRsN-FLX4BMdBKkbYqdTuge_xFUO-VWv1EqvnDDcOktTZy1UbFxGnaXCA83n0v0OQG2ULdqUdSgRVtZUW11PUoGeDqejnCNI35WjraiT0BG0HC6gRD27JPRW_4qwJwKz0hiUiTfhOT1zJ8h5d3U8WqYQ6ZI_QPS0F3_FC2awu39FSZoOViXuoYCN_T5j5wbfc3qIW_xt2GodoVY9iC7dFleeyFd4veQbmCOKoDG79nHARMkkzgdW-Eie1-OAp5zOH4mWdvVuCkNBj5Bh7giYku4ufYX8MA")
    out.dat$samplesize.outcome <- 1637+212242
  }
  
  out.dat <- out.dat[!duplicated(out.dat$SNP),]
  return(out.dat)
}

my_check_dataset <- function (d, suffix = "", req = c("snp"), warn.minp = 1e-06){
  if (!is.list(d)) 
    stop("dataset ", suffix, ": is not a list")
  nd <- names(d)
  n <- 0
  for (v in nd) {
    if (v %in% req && !(v %in% nd)) 
      stop("dataset ", suffix, ": missing required element ", 
           v)
    if (any(is.na(d[[v]]))) 
      stop("dataset ", suffix, ": ", v, " contains missing values")
  }
  if ("snp" %in% nd && any(duplicated(d$snp))) 
    stop("dataset ", suffix, ": duplicated snps found")
  if ("snp" %in% nd && is.factor(d$snp)) 
    stop("dataset ", suffix, ": snp should be a character vector but is a factor")
  if ("MAF" %in% nd && (!is.numeric(d$MAF) || any(is.na(d$MAF)) || 
                        any(d$MAF <= 0) || any(d$MAF >= 1))) 
    stop("dataset ", suffix, ": MAF should be a numeric, strictly >0 & <1")
  l <- -1
  shouldmatch <- c("pvalues", "MAF", "beta", "varbeta", "snp", 
                   "position")
  for (v in shouldmatch) if (v %in% nd) 
    if (l < 0) {
      l <- length(d[[v]])
    }
  else {
    if (length(d[[v]]) != l) {
      stop("dataset ", suffix, ": lengths of inputs don't match: ")
      print(intersect(nd, shouldmatch))
    }
  }
  if (("N" %in% req) && (!("N" %in% nd) || is.null(d$N) || 
                         is.na(d$N))) 
    stop("dataset ", suffix, ": sample size N not set")
  if (!("type" %in% nd)) 
    stop("dataset ", suffix, ": variable type not set")
  if (any(!(d$type %in% c("quant", "cc"))))
    stop("dataset ", suffix, ": ", "type must be quant or cc")
  if (("s" %in% nd) && (!is.numeric(d$s) || d$s <= 0 || d$s >= 1)) 
    stop("dataset ", suffix, ": ", "s must be between 0 and 1")
  if (!("beta" %in% nd) || !("varbeta" %in% nd)) {
    if (!("pvalues" %in% nd) || !("MAF" %in% nd)) 
      stop("dataset ", suffix, ": ", "require p values and MAF if beta, varbeta are unavailable")
    if (all(d$type == "cc") && !("s" %in% nd)) 
      stop("dataset ", suffix, ": ", "require, s, proportion of samples who are cases, if beta, varbeta are unavailable")
    p = d$pvalues
  }
  else {
    p = pnorm(-abs(d$beta/sqrt(d$varbeta))) * 2
  }
  if (min(p) > warn.minp) 
    warning("minimum p value is: ", format.pval(min(p)), 
            "\nIf this is what you expected, this is not a problem.\nIf this is not as small as you expected, please check the 02_data vignette.")
  if (all(d$type == "quant") && !("sdY" %in% nd)) 
    if (!("MAF" %in% nd && "N" %in% nd)) 
      stop("dataset ", suffix, ": ", "must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")
  if ("LD" %in% nd) {
    if (nrow(d$LD) != ncol(d$LD)) 
      stop("LD not square")
    if (!identical(colnames(d$LD), rownames(d$LD))) 
      stop("LD rownames != colnames")
    if (length(setdiff(d$snp, colnames(d$LD)))) 
      stop("colnames in LD do not contain all SNPs")
  }
  NULL
}

my.coloc.abf <- function (dataset1, dataset2, MAF = NULL, p1 = 1e-04, p2 = 1e-04,  p12 = 1e-05) {
  if (!("MAF" %in% names(dataset1)) & !is.null(MAF)) 
    dataset1$MAF <- MAF
  if (!("MAF" %in% names(dataset2)) & !is.null(MAF)) 
    dataset2$MAF <- MAF
  my_check_dataset(d = dataset1, 1)
  my_check_dataset(d = dataset2, 2)
  df1 <- my.process.dataset(d = dataset1, suffix = "df1")
  df2 <- my.process.dataset(d = dataset2, suffix = "df2")
  merged.df <- merge(df1, df2)
  if (!nrow(merged.df)) 
    stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")
  merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + 
                                        lABF.df2)
  my.denom.log.abf <- coloc:::logsum(merged.df$internal.sum.lABF)
  merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - 
                               my.denom.log.abf)
  pp.abf <- coloc:::combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, 
                        p1, p2, p12)
  common.snps <- nrow(merged.df)
  results <- c(nsnps = common.snps, pp.abf)
  output <- list(summary = results, results = merged.df, priors = c(p1 = p1, 
                                                                    p2 = p2, p12 = p12))
  class(output) <- c("coloc_abf", class(output))
  return(output)
}

my.process.dataset <- function (d, suffix){
  nd <- names(d)
  if (!"type" %in% nd) 
    stop("dataset ", suffix, ": ", "The variable type must be set, otherwise the Bayes factors cannot be computed")
  if (any(!(d$type %in% c("quant", "cc"))))
    stop("dataset ", suffix, ": ", "type must be quant or cc")
  if (all(d$type == "cc") & "pvalues" %in% nd) {
    if (!("s" %in% nd)) 
      stop("dataset ", suffix, ": ", "please give s, proportion of samples who are cases, if using p values")
    if (!("MAF" %in% nd)) 
      stop("dataset ", suffix, ": ", "please give MAF if using p values")
    if (any(d$s <= 0) || any(d$s >= 1))
      stop("dataset ", suffix, ": ", "s must be between 0 and 1")
  }
  if (all(d$type == "quant")) {
    if (!("sdY" %in% nd || ("MAF" %in% nd && "N" %in% nd))) 
      stop("dataset ", suffix, ": ", "must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")
  }
  if ("beta" %in% nd && "varbeta" %in% nd) {
    if (length(d$beta) != length(d$varbeta)) 
      stop("dataset ", suffix, ": ", "Length of the beta vectors and variance vectors must match")
    if (!("snp" %in% nd)) 
      d$snp <- sprintf("SNP.%s", 1:length(d$beta))
    if (length(d$snp) != length(d$beta)) 
      stop("dataset ", suffix, ": ", "Length of snp names and beta vectors must match")
    if (all(d$type == "quant") && !("sdY" %in% nd)) 
      d$sdY <- coloc:::sdY.est(d$varbeta, d$MAF, d$N)
    df <- my.approx.bf.estimates(z = d$beta/sqrt(d$varbeta), 
                              V = d$varbeta, type = d$type, suffix = suffix, sdY = d$sdY)
    df$snp <- as.character(d$snp)
    if ("position" %in% nd) 
      df <- cbind(df, position = d$position)
    return(df)
  }
  if ("pvalues" %in% nd & "MAF" %in% nd & "N" %in% nd) {
    if (length(d$pvalues) != length(d$MAF)) 
      stop("Length of the P-value vectors and MAF vector must match")
    if (!("snp" %in% nd)) 
      d$snp <- sprintf("SNP.%s", 1:length(d$pvalues))
    df <- data.frame(pvalues = d$pvalues, MAF = d$MAF, N = d$N, 
                     snp = as.character(d$snp))
    snp.index <- which(colnames(df) == "snp")
    colnames(df)[-snp.index] <- paste(colnames(df)[-snp.index], 
                                      suffix, sep = ".")
    keep <- which(df$MAF > 0 & df$pvalues > 0)
    df <- df[keep, ]
    abf <- my.approx.bf.p(p = df$pvalues, f = df$MAF, type = d$type, 
                       N = df$N, s = d$s, suffix = suffix)
    df <- cbind(df, abf)
    if ("position" %in% nd) 
      df <- cbind(df, position = d$position[keep])
    return(df)
  }
  stop("Must give, as a minimum, one of:\n(beta, varbeta, type, sdY)\n(beta, varbeta, type, MAF)\n(pvalues, MAF, N, type)")
}

my.approx.bf.estimates <- function (z, V, type, suffix = NULL, sdY = 1) {
  sd.prior <- if(all(type == "quant")) {
    0.15 * sdY
  }
  else {
    0.2
  }
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  ret <- data.frame(V, z, r, lABF)
  if (!is.null(suffix)) 
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}

my.approx.bf.p <- function (p, f, type, N, s, suffix = NULL){
  if (all(type == "quant")) {
    sd.prior <- 0.15
    V <- Var.data(f, N)
  }
  else {
    sd.prior <- 0.2
    V <- Var.data.cc(f, N, s)
  }
  z <- qnorm(0.5 * p, lower.tail = FALSE)
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  ret <- data.frame(V, z, r, lABF)
  if (!is.null(suffix)) 
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}