library(dplyr)
library(readr)
library(parallel)
library(corpora)

memCheck <- function(show=FALSE) {res <- gc(); gc(reset=TRUE); if (!show) return(res); print(res); invisible(res)}

# !! we still use tidiverse to read and prepare data because base R can't do efficient group-by operations
# !! (so even getting the global joint and marginal frequencies is more than painful)
# !! but we'll try to get by with base R functionality inside the bootstrapping algorithm
bnc <- read_tsv("bnc_adj_n.tsv", quote="", col_names=qw("f text l1 l2"))
dim(bnc)
print(object.size(bnc), units="MiB")

# determine marginal frequencies and sample size for each text --> input for bootstrap.test()
bnc |> summarise(N = sum(f), .by=text) -> sample.size
bnc |> summarise(f1 = sum(f), .by=c(text, l1)) -> margin.adj
bnc |> summarise(f2 = sum(f), .by=c(text, l2)) -> margin.noun
memCheck(TRUE)

bnc |> select(text, l1, l2, f) -> bnc

# compute AM as efficiently as possible (avoiding intermediate terms if we can)
AM <- function (measure=qw("MI2 simple-ll"), f, f1, f2, N) {
  measure <- match.arg(measure)
  stopifnot(length(f) == length(f1))
  stopifnot(length(f) == length(f2))
  stopifnot(length(N) == 1 || length(f) == length(N))
  if (measure == "MI2") {
    stopifnot(all(f > 0))
    O11 <- f
    E11 <- f1 * f2 / N
    log2(O11 * O11 / E11)
  }
  else if (measure == "simple-ll") {
    O11 <- f
    E11 <- f1 * f2 / N
    res <- O11 * log2(O11 * O11 / E11)
    res[O11 == 0] <- 0
    res
  }
  else {
    stop("not implemented")
  }
}


# test bootstrapping algorithms using base R functionality only
#   - mem=TRUE:     run GC and print memory usage at intermediate steps (reduces RAM overhead)
#   - verbose=TRUE: show progress bar (should be used in interactive mode only)
#   - check=FALSE:  disable consistency checks for slightly faster performance
bootstrap.test <- function (cooc, margin1, margin2, size, measure, min.tf=NULL, min.df=NULL, replicates=10, details=FALSE,
                            parallel=1L, mem=FALSE, check=TRUE, verbose=TRUE) {
  cooc <- as.data.frame(cooc) # make sure we don't accidentally rely on special features of tibbles
  margin1 <- as.data.frame(margin1)
  margin2 <- as.data.frame(margin2)
  size <- as.data.frame(size)
  if (check) {
    stopifnot(all(qw("text l1 l2 f") %in% colnames(cooc)))
    stopifnot(all(qw("text l1 f1") %in% colnames(margin1)))
    stopifnot(all(qw("text l2 f2") %in% colnames(margin2)))
    stopifnot(all(qw("text N") %in% colnames(size)))
  }
  if (check) {
    stopifnot(!any(duplicated(size$text)))
    stopifnot(all(cooc$text %in% size$text))
    stopifnot(all(margin1$text %in% size$text))
    stopifnot(all(margin2$text %in% size$text))
    stopifnot(all(cooc$l1 %in% margin1$l1))
    stopifnot(all(cooc$l2 %in% margin2$l2))
  }
  if (parallel > 1) verbose <- FALSE # progress bar doesn't work well with parallelisation

  if (verbose) pb <- txtProgressBar(min=0, max=replicates + 1, style=3)

  cooc <- transform(cooc, pair=factor(paste(l1, l2, sep="\t"))) # code pair types as factor

  # apply global tf / df threshold on pair types (unless done in advance on cooc, which is recommended)
  total_f <- rowsum(
    cbind(cooc$f, rep(1, nrow(cooc))),
    cooc$pair, reorder=TRUE)
  pairs <- data.frame(pair=rownames(total_f), tf=total_f[, 1], df=total_f[, 2], row.names=NULL)
  if (!is.null(min.tf) || !is.null(min.df)) {
    if (!is.null(min.tf)) pairs <- subset(pairs, tf >= min.tf)
    if (!is.null(min.df)) pairs <- subset(pairs, df >= min.df)
    cooc <- droplevels(subset(cooc, pair %in% pairs$pair)) # filter raw cooc data
    if (mem) memCheck(TRUE)
  }
  # split filtered pairs into l1 and l2 (strsplit() creates list, which is probably less efficient)
  pairs <- transform(
    pairs,
    l1 = sub("\t.*$", "", pair, perl=TRUE),
    l2 = sub("^.*\t", "", pair, perl=TRUE))
  pairs$pair <- factor(pairs$pair) # convert to factor with same ordering as in cooc
  if (check) stopifnot(all.equal( levels(pairs$pair), levels(cooc$pair) ))
  if (mem) memCheck(TRUE)

  text.id <- unique(size$text)
  n.texts <- length(text.id)

  if (!(parallel > 1)) {
    # -- standard iterative algorithm
    boot.res <- matrix(nrow=nrow(pairs), ncol=replicates) # per-allocate boostrap data table
    for (i in seq_len(replicates)) {
      if (verbose) setTxtProgressBar(pb, i)

      # obtain bootstrapped per-text counts
      text.cnt <- structure(as.vector(rmultinom(1, n.texts, rep(1, n.texts))), names=text.id)
      cooc <- transform(cooc, f.boot = f * text.cnt[text])
      margin1 <- transform(margin1, f1.boot = f1 * text.cnt[text])
      margin2 <- transform(margin2, f2.boot = f2 * text.cnt[text])
      size <- transform(size, N.boot = N * text.cnt[text])
      if (mem) memCheck(TRUE)

      # aggregate bootstrapped counts in vectors with rowsum
      f <- rowsum(cooc$f.boot, cooc$pair, reorder=TRUE)[, 1] # extract single column vector with labels preserved
      f1 <- rowsum(margin1$f1.boot, margin1$l1)[, 1]
      f2 <- rowsum(margin2$f2.boot, margin2$l2)[, 1]
      N <- sum(size$N.boot)

      # look up marginal frequencies for full frequency signatures
      if (check) stopifnot(all.equal( names(f), as.character(pairs$pair) ))
      f1 <- f1[ pairs$l1 ]
      f2 <- f2[ pairs$l2 ]
      if (mem) memCheck(TRUE)

      boot.res[, i] <- AM(measure, f, f1, f2, N)
      rm(f, f1, f2, N)
      if (mem) memCheck(TRUE)
    }

    if (verbose) setTxtProgressBar(pb, replicates + 1)
  }
  else {
    # -- parallel processing (here with fork(), using a cluster will be much worse due to data transfer)
    .worker <- function (n) {
      text.cnt <- structure(as.vector(rmultinom(1, n.texts, rep(1, n.texts))), names=text.id)
      cooc <- transform(cooc, f.boot = f * text.cnt[text])
      margin1 <- transform(margin1, f1.boot = f1 * text.cnt[text])
      margin2 <- transform(margin2, f2.boot = f2 * text.cnt[text])
      size <- transform(size, N.boot = N * text.cnt[text])

      # aggregate bootstrapped counts in vectors with rowsum
      f <- rowsum(cooc$f.boot, cooc$pair, reorder=TRUE)[, 1] # extract single column vector with labels preserved
      f1 <- rowsum(margin1$f1.boot, margin1$l1)[, 1]
      f2 <- rowsum(margin2$f2.boot, margin2$l2)[, 1]
      N <- sum(size$N.boot)

      # look up marginal frequencies for full frequency signatures
      if (check) stopifnot(all.equal( names(f), as.character(pairs$pair) ))
      f1 <- f1[ pairs$l1 ]
      f2 <- f2[ pairs$l2 ]

      AM(measure, f, f1, f2, N)
    }
    boot.res <- mclapply(seq_len(replicates), .worker, mc.cores=parallel)
    boot.res <- do.call(cbind, boot.res) # convert into expected matrix format
  }

  if (details) {
    # return full bootstrap data for detailed analysis
    pairs <- cbind(pairs[, qw("l1 l2")], boot.res)
  }
  else {
    # return summary statistics (to be improved)
    pairs$mean <- rowMeans(boot.res)
    pairs$sd <- apply(boot.res, 1, sd)
    pairs$med <- apply(boot.res, 1, median)
  }

  if (verbose) close(pb)
  pairs[, -1]
}

## -- these lines are for interactive testing and exploration
if (FALSE) {
  system.time(res <- bootstrap.test(bnc, margin.adj, margin.noun, sample.size, "simple-ll", replicates=10,
                                    min.tf=50, parallel=5, mem=FALSE, details=FALSE, verbose=TRUE, check=TRUE)); memCheck()
  head(res, 10) # already sorted by pair (l1, l2)
}

## benchmark different algorithms when running as script
##  - consistency checks are disabled for optimal performance
##  - print first line of bootstrapped table to verify functionality
cat("\n\n-- base R algorithm\n")
system.time(res <- bootstrap.test(bnc, margin.adj, margin.noun, sample.size, "simple-ll", replicates=50,
                                  min.tf=5, mem=FALSE, details=FALSE, verbose=FALSE, check=FALSE))
memCheck()
print(head(res, 1))

cat("\n\n-- parallel on 4 cores (watch Activity Monitor for true RAM usage)\n")
system.time(res <- bootstrap.test(bnc, margin.adj, margin.noun, sample.size, "simple-ll", replicates=50,
                                  min.tf=5, parallel=4, mem=FALSE, details=FALSE, verbose=FALSE, check=FALSE))
memCheck()
print(head(res, 1))

cat("\n\n-- parallel on 8 cores (watch Activity Monitor for true RAM usage)\n")
system.time(res <- bootstrap.test(bnc, margin.adj, margin.noun, sample.size, "simple-ll", replicates=50,
                                  min.tf=5, parallel=8, mem=FALSE, details=FALSE, verbose=FALSE, check=FALSE))
memCheck()
print(head(res, 1))
