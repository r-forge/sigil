library(data.table)
library(corpora)

memCheck <- function(show=FALSE) {res <- gc(); gc(reset=TRUE); if (!show) return(res); print(res); invisible(res)}

# per-text adj-noun frequency counts from BNC
bnc <- fread("bnc_adj_n.tsv", sep="\t", header=FALSE, col.names=qw("f text l1 l2"), encoding="UTF-8")
bnc[, f := as.numeric(f)] # make sure we use floating-point for all calculations
dim(bnc)
print(object.size(bnc), units="MiB")

# aggregate global counts for entire corpus, then annotate marginal frequencies
global <- bnc[, .(f = sum(f)), by=.(l1, l2)]
global.adj <- global[, .(f1 = sum(f)), by=l1]
global.noun <- global[, .(f2 = sum(f)), by=l2]
global <- global |> merge(global.adj, by="l1") |> merge(global.noun, by="l2")
global[, N := sum(f)]
setkey(global, l1, l2)
setcolorder(global, qw("l1 l2 f f1 f2 N"))
print(object.size(global), units="MiB") # ~ 55 MiB, quite small
memCheck(TRUE)

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

# global association scores
system.time(global[, MI2 := AM("MI2", f, f1, f2, N)])
system.time(global[, LL := AM("simple-ll", f, f1, f2, N)])
memCheck(TRUE)

# determine marginal frequencies and sample size for each text --> simple input for association.measure.boot() ?
sample.size <- bnc[, .(N = sum(f)), by=text]
margin.adj <- bnc[, .(f1 = sum(f)), by=.(text, l1)]
margin.noun <- bnc[, .(f2 = sum(f)), by=.(text, l2)]
memCheck(TRUE)

# merge into main table
# bnc <- merge(bnc, sample.size, by="text")
# bnc <- merge(bnc, margin.adj, by=qw("text adj"))
# bnc <- merge(bnc, margin.noun, by=qw("text noun"))
# dim(bnc)
# print(object.size(bnc), units="MiB")

# aggregate overall pair frequency and df for frequency thresholds
# pairs <- bnc[, .(tf=sum(f), df=.N), by=.(adj, noun)]
# bnc <- merge(bnc, pairs, by=qw("adj noun"))
# dim(bnc)
# print(object.size(bnc), units="MiB") # ~ 226 MiB for ~ 4.6M rows (1.4M pair types across 4045 texts)

# restore preferred column order
# setcolorder(bnc, qw("text adj noun f f1 f2 N tf df"))
setcolorder(bnc, qw("text l1 l2 f"))

# apply frequency threshold (bootstrapping makes little sense for a hapax pair type)
# nrow(pairs) # 1.46M pair types
# table(pairs$tf)[1:10] # threshold tf >= 5 leaves less than 200k pair types
# table(pairs$df)[1:10] # similar for df; but cases with df==1 vs. tf>=5 are particularly interesting
# bnc5 <- bnc[tf >= 5]
# dim(bnc5)
# print(object.size(bnc5), units="MiB") # ~ 136 MiB for ~ 2.9M rows


# test different bootstrapping algorithms
#   - j1=TRUE:  use data.table joins on keyed tables
#   - j1=FALSE: use standard R vector lookup
#   - j2=TRUE:  merge by index joins with temporary tables
#   - j2=FALSE: merge with merge() for data.tables
#   - mem=TRUE: run GC and print memory usage at intermediate steps (reduces RAM overhead)
#   - verbose=TRUE: show progress bar (should be used in interactive mode only)
#   - check=FALSE:  disable consistency checks for slightly faster performance
bootstrap.test <- function (cooc, margin1, margin2, size, measure, min.tf=NULL, min.df=NULL, replicates=10, details=FALSE,
                            j1=FALSE, j2=TRUE, mem=FALSE, check=TRUE, verbose=TRUE) {
  if (check) {
    stopifnot(all(qw("text l1 l2 f") %in% colnames(cooc)))
    stopifnot(all(qw("text l1 f1") %in% colnames(margin1)))
    stopifnot(all(qw("text l2 f2") %in% colnames(margin2)))
    stopifnot(all(qw("text N") %in% colnames(size)))
  }
  setkey(cooc, l1, l2)
  setkey(margin1, l1)
  setkey(margin2, l2)

  if (check) {
    stopifnot(!any(duplicated(size$text)))
    stopifnot(all(cooc$text %in% size$text))
    stopifnot(all(margin1$text %in% size$text))
    stopifnot(all(margin2$text %in% size$text))
    stopifnot(all(cooc$l1 %in% margin1$l1))
    stopifnot(all(cooc$l2 %in% margin2$l2))
  }

  if (verbose) pb <- txtProgressBar(min=0, max=replicates + 1, style=3)
  pairs <- cooc[, .(tf=sum(f), df=.N), keyby=.(l1, l2)]

  # apply global tf / df threshold on pairs (unless done in advance on cooc, which is recommended)
  if (!is.null(min.tf) || !is.null(min.df)) {
    if (!is.null(min.tf)) pairs <- pairs[tf >= min.tf]
    if (!is.null(min.df)) pairs <- pairs[df >= min.df]
    cooc <- cooc[pairs[, .(l1, l2)], on=qw("l1 l2"), nomatch=NULL]
    if (mem) print(memCheck())
  }

  text.id <- unique(size$text)
  n.texts <- length(text.id)
  boot.res <- matrix(nrow=nrow(pairs), ncol=replicates) # per-allocate boostrap data table
  if (j1) bootcnt <- data.table(text=text.id, cnt=0, key="text") # pre-allocate data.table for j1 method
  for (i in seq_len(replicates)) {
    if (verbose) setTxtProgressBar(pb, i)
    if (j1) {
      # use data.table joins on keyed table (j1=TRUE)
      bootcnt[, cnt := as.vector(rmultinom(1, n.texts, rep(1, n.texts)))]
      cooc[bootcnt, f.boot := f * i.cnt, on="text"]
      margin1[bootcnt, f1.boot := f1 * i.cnt, on="text"]
      margin2[bootcnt, f2.boot := f2 * i.cnt, on="text"]
      size[bootcnt, N.boot := N * i.cnt, on="text"]
    }
    else {
      # use standard R vector lookup (j1=FALSE)
      text.cnt <- structure(as.vector(rmultinom(1, n.texts, rep(1, n.texts))), names=text.id)
      cooc[, f.boot := f * text.cnt[text]]
      margin1[, f1.boot := f1 * text.cnt[text]]
      margin2[, f2.boot := f2 * text.cnt[text]]
      size[, N.boot := N * text.cnt[text]]
    }
    if (mem) print(memCheck())
    if (j2) {
      # merge by index joins with temporary tables (j2=TRUE)
      boot <- cooc[, .(f = sum(f.boot)), keyby=.(l1, l2)]
      tmp <- margin1[, .(f1 = sum(f1.boot)), keyby=l1]
      boot[tmp, f1 := i.f1, on="l1"]
      tmp <- margin2[, .(f2 = sum(f2.boot)), keyby=l2]
      boot[tmp, f2 := i.f2, on="l2"]
    }
    else {
      # merge with merge.data.table (j2=FALSE)
      boot <- cooc[, .(f = sum(f.boot)), keyby=.(l1, l2)] |>
        merge(margin1[, .(f1 = sum(f1.boot)), by=l1], by="l1") |>
        merge(margin2[, .(f2 = sum(f2.boot)), by=l2], by="l2")
      setkey(boot, l1, l2)
    }
    N <- sum(size$N.boot)
    boot[, am := AM(measure, f, f1, f2, N)]
    if (check) stopifnot(all(pairs$l1 == boot$l1) && all(pairs$l2 == boot$l2)) # in debug mode
    boot.res[, i] <- boot$am
    if (mem) print(memCheck())
  }

  if (verbose) setTxtProgressBar(pb, replicates + 1)

  if (details) {
    # return full bootstrap data for detailed analysis
    pairs <- cbind(pairs[, .(l1, l2)], boot.res)
  }
  else {
    # return summary statistics (to be improved)
    pairs$mean <- rowMeans(boot.res)
    pairs$sd <- apply(boot.res, 1, sd)
    pairs$med <- apply(boot.res, 1, median)
  }

  if (verbose) close(pb)
  pairs
}

# -- these lines are for interactive testing and exploration
if (FALSE) {
  system.time(res <- bootstrap.test(bnc, margin.adj, margin.noun, sample.size, "simple-ll", replicates=50,
                                    min.tf=5, j1=FALSE, j2=TRUE, mem=FALSE, details=FALSE, verbose=TRUE, check=TRUE)); memCheck()
  res
}

# benchmark different algorithms when running as script
#  - consistency checks are disabled for optimal performance
#  - print first line of bootstrapped table to verify functionality
cat("\n\n-- standard R vector lookup + merge()\n")
system.time(res <- bootstrap.test(bnc, margin.adj, margin.noun, sample.size, "simple-ll", replicates=50,
                                  min.tf=5, j1=FALSE, j2=FALSE, mem=FALSE, details=FALSE, verbose=FALSE, check=FALSE))
memCheck()
print(head(res, 1))

cat("\n\n-- join on keyed data.tables + merge()\n")
system.time(res <- bootstrap.test(bnc, margin.adj, margin.noun, sample.size, "simple-ll", replicates=50,
                                  min.tf=5, j1=TRUE, j2=FALSE, mem=FALSE, details=FALSE, verbose=FALSE, check=FALSE))
memCheck()
print(head(res, 1))

cat("\n\n-- standard R vector lookup + index join with temporary data.table\n")
system.time(res <- bootstrap.test(bnc, margin.adj, margin.noun, sample.size, "simple-ll", replicates=50,
                                  min.tf=5, j1=FALSE, j2=TRUE, mem=FALSE, details=FALSE, verbose=FALSE, check=FALSE))
memCheck()
print(head(res, 1))

cat("\n\n-- join on keyed data.tables + index join with temporary data.table\n")
system.time(res <- bootstrap.test(bnc, margin.adj, margin.noun, sample.size, "simple-ll", replicates=50,
                                  min.tf=5, j1=TRUE, j2=TRUE, mem=FALSE, details=FALSE, verbose=FALSE, check=FALSE))
memCheck()
print(head(res, 1))





