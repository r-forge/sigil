library(dplyr)
library(readr)
library(corpora)

memCheck <- function(show=FALSE) {res <- gc(); gc(reset=TRUE); if (!show) return(res); print(res); invisible(res)}

# per-text adj-noun frequency counts from BNC
bnc <- read_tsv("bnc_adj_n.tsv", quote="", col_names=qw("f text l1 l2"))
dim(bnc)
print(object.size(bnc), units="MiB")

# aggregate global counts for entire corpus, then annotate marginal frequencies
bnc |> summarise(f = sum(f), .by=c(l1, l2)) -> global
global |> summarise(f1 = sum(f), .by=l1) -> global.adj
global |> summarise(f2 = sum(f), .by=l2) -> global.noun
global |>
  inner_join(global.adj, by="l1") |>
  inner_join(global.noun, by="l2") |>
  mutate(N = sum(f)) ->
  global
print(object.size(global), units="MiB") # ~ 77 MiB
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
system.time(global |> mutate(MI2 = AM("MI2", f, f1, f2, N)) -> global)
system.time(global |> mutate(LL = AM("simple-ll", f, f1, f2, N)) -> global)
memCheck(TRUE)

# determine marginal frequencies and sample size for each text --> input for bootstrap.test()
bnc |> summarise(N = sum(f), .by=text) -> sample.size
bnc |> summarise(f1 = sum(f), .by=c(text, l1)) -> margin.adj
bnc |> summarise(f2 = sum(f), .by=c(text, l2)) -> margin.noun
memCheck(TRUE)

bnc |> select(text, l1, l2, f) -> bnc

# test bootstrapping algorithms using tidyverse() functionality
#   - method="basic":     direct translation of the simplest data.table() implementation
#   - method="optimised": try avoiding intermediate data and allocations
#   - sort=TRUE:    pre-sort frequency tables alphabetically by l1 and l2
#   - mem=TRUE:     run GC and print memory usage at intermediate steps (reduces RAM overhead)
#   - verbose=TRUE: show progress bar (should be used in interactive mode only)
#   - check=FALSE:  disable consistency checks for slightly faster performance
bootstrap.test <- function (cooc, margin1, margin2, size, measure, min.tf=NULL, min.df=NULL, replicates=10, details=FALSE,
                            method=c("basic", "optimised"), sort=FALSE, mem=FALSE, check=TRUE, verbose=TRUE) {
  method <- match.arg(method)
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

  if (verbose) pb <- txtProgressBar(min=0, max=replicates + 1, style=3)

  if (sort) {
    cooc |> arrange(l1, l2, text) -> cooc
    margin1 |> arrange(l1, text) -> margin1
    margin2 |> arrange(l2, text) -> margin2
  }

  # apply global tf / df threshold on pair types (unless done in advance on cooc, which is recommended)
  cooc |> summarise(tf = sum(f), df = n(), .by=c(l1, l2)) -> pairs
  if (!is.null(min.tf) || !is.null(min.df)) {
    if (!is.null(min.tf)) pairs |> filter(tf >= min.tf) -> pairs
    if (!is.null(min.df)) pairs |> filter(df >= min.df) -> pairs
    cooc |> semi_join(pairs, by=c("l1", "l2")) -> cooc
    if (mem) memCheck(TRUE)
  }

  text.id <- unique(size$text)
  n.texts <- length(text.id)
  boot.res <- matrix(nrow=nrow(pairs), ncol=replicates) # per-allocate boostrap data table
  for (i in seq_len(replicates)) {
    if (verbose) setTxtProgressBar(pb, i)

    if (method == "basic") {
      # obtain bootstrapped per-text counts by standard R vector lookup
      text.cnt <- structure(as.vector(rmultinom(1, n.texts, rep(1, n.texts))), names=text.id)
      cooc |> mutate(f.boot = f * text.cnt[text]) -> cooc
      margin1 |> mutate(f1.boot = f1 * text.cnt[text]) -> margin1
      margin2 |> mutate(f2.boot = f2 * text.cnt[text]) -> margin2
      size |> mutate(N.boot = N * text.cnt[text]) -> size
      if (mem) memCheck(TRUE)
      # aggregate bootstrapped counts in temporary tables
      cooc |> summarise(f = sum(f.boot), .by=c(l1, l2)) -> boot
      margin1 |> summarise(f1 = sum(f1.boot), .by=l1) -> tmp1
      margin2 |> summarise(f2 = sum(f2.boot), .by=l2) -> tmp2
      # merge frequency signatures with inner joins
      boot |>
        inner_join(tmp1, by="l1", relationship="many-to-one") |>
        inner_join(tmp2, by="l2", relationship="many-to-one") ->
        boot
      N <- sum(size$N.boot)
      rm(tmp1, tmp2, text.cnt)
    }
    else if (method == "optimised") {
      # obtain bootstrapped per-text counts by standard R vector lookup,
      # but aggregate them directly into temporary tables
      # NB: vector lookup within the summarise() turned out to be extremely inefficient
      text.cnt <- structure(as.vector(rmultinom(1, n.texts, rep(1, n.texts))), names=text.id)
      cooc |> mutate(f.boot = f * text.cnt[text]) |> summarise(f = sum(f.boot), .by=c(l1, l2)) -> boot
      margin1 |> mutate(f1.boot = f1 * text.cnt[text]) |> summarise(f1 = sum(f1.boot), .by=l1) -> tmp1
      margin2 |> mutate(f2.boot = f2 * text.cnt[text]) |> summarise(f2 = sum(f2.boot), .by=l2) -> tmp2
      if (mem) memCheck(TRUE)
      # merge frequency signatures with inner joins
      boot |>
        inner_join(tmp1, by="l1") |>
        inner_join(tmp2, by="l2") ->
        boot
      size |> summarise(N = sum(N * text.cnt[text])) |> pull(N) -> N
      rm(tmp1, tmp2, text.cnt)
    }
    if (mem) memCheck(TRUE)

    if (check) stopifnot(all(pairs$l1 == boot$l1) && all(pairs$l2 == boot$l2)) # in debug mode
    boot.res[, i] <- AM(measure, boot$f, boot$f1, boot$f2, N)
    if (mem) memCheck(TRUE)
  }

  if (verbose) setTxtProgressBar(pb, replicates + 1)

  if (details) {
    # return full bootstrap data for detailed analysis
    pairs |> select(l1, l2) |> cbind(boot.res) -> pairs
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

## -- these lines are for interactive testing and exploration
if (FALSE) {
  system.time(res <- bootstrap.test(bnc, margin.adj, margin.noun, sample.size, "simple-ll", replicates=10,
                                    min.tf=10, method="opt", sort=TRUE, mem=FALSE, details=FALSE, verbose=TRUE, check=TRUE)); memCheck()
  arrange(res, l1, l2)
}

## benchmark different algorithms when running as script
##  - consistency checks are disabled for optimal performance
##  - print first line of bootstrapped table to verify functionality
cat("\n\n-- basic algorithm\n")
system.time(res <- bootstrap.test(bnc, margin.adj, margin.noun, sample.size, "simple-ll", replicates=50,
                                  min.tf=5, method="basic", sort=FALSE, mem=FALSE, details=FALSE, verbose=FALSE, check=FALSE))
memCheck()
print(arrange(res, l1, l2), n=1)

cat("\n\n-- optimised algorithm\n")
system.time(res <- bootstrap.test(bnc, margin.adj, margin.noun, sample.size, "simple-ll", replicates=50,
                                  min.tf=5, method="optimised", sort=FALSE, mem=FALSE, details=FALSE, verbose=FALSE, check=FALSE))
memCheck()
print(arrange(res, l1, l2), n=1)

cat("\n\n-- optimised algorithm + sorted tables\n")
system.time(res <- bootstrap.test(bnc, margin.adj, margin.noun, sample.size, "simple-ll", replicates=50,
                                  min.tf=5, method="optimised", sort=TRUE, mem=FALSE, details=FALSE, verbose=FALSE, check=FALSE))
memCheck()
print(arrange(res, l1, l2), n=1)
