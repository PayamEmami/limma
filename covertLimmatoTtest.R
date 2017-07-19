convertLimmaToregularTtest<-function(fit)
{
  require(limma)
  ebfit <- eBayes(fit)
  ebfit$metabolites <- ebfit$genes
  ebfit$genes <- ebfit$rank <- ebfit$assign <- NULL
  ebfit$qr <- ebfit$qraux <- ebfit$pivot <- ebfit$tol <- NULL
  ebfit$cov.coefficients <- ebfit$pivot <- ebfit$lods <- NULL
  beta <- ebfit$coeff
  tstat <- sweep(data.matrix(ebfit$coef/ebfit$stdev.unscaled), 
                 1, ebfit$sigma, "/")
  tpval <- 2 * pt(-abs(tstat), df = ebfit$df.residual)
  ebfit$t <- tstat
  df <- ebfit$df.residual
  fstat <- classifyTestsF(ebfit, df = df, fstat.only = TRUE)
  Fstat <- as.vector(fstat)
  df1 <- attr(fstat, "df1")
  df2 <- attr(fstat, "df2")
  if (df2[1] > 1e+06) {
    Fpval <- pchisq(df1 * Fstat, df1, lower.tail = FALSE)
  }
  else {
    Fpval <- pf(Fstat, df1, df2, lower.tail = FALSE)
  }
  se <- beta/tstat
  ebfit$F <- Fstat
  ebfit$F.p.value <- Fpval
  ebfit$p.value <- tpval
  ebfit$df.total <- ebfit$s2.post <- ebfit$stdev.unscaled <- NULL
  ebfit$var.prior <- ebfit$proportion <- ebfit$s2.prior <- NULL
  ebfit$df.prior <- NULL
  ebfit$std.error <- se
  ebfit$t <- tstat
  ebfit$df <- df
  
  return(ebfit)
}