##
## test am.score() against UCS reference implementation (Evert 2004)
##
library(corpora)

## KrennPPV data set includes both frequency signatures and reference scores
PPV <- KrennPPV[, 1:9] # without AM scores

## original data set has z-score and chi-squared without Yates' correction
my.z.score <- function(O, E, ...) (O - E) / sqrt(E)
my.chisq <- function(O11, O12, O21, O22, E11, E12, E21, E22, ...) {
  X2 <- (O11 - E11)^2 / E11 + (O12 - E12)^2 / E12 + (O21 - E21)^2 / E21 + (O22 - E22)^2 / E22
  sign(O11 - E11) * X2
}

## UCS implementation doesn't properly handle corner cases for Yates' correction,
## so provide a matching implementation here
my.chisq.corr <- function (E11, O11, O12, O21, O22, R1, R2, C1, C2, N, ...) {
  term <- abs(O11 * O22 - O12 * O21) - N/2
  X2 <- (N * term^2) / (R1 * R2 * C1 * C2)
  sign(O11 - E11) * X2
}

## now recompute AM scores
PPV <- transform(
  PPV,
  MI = am.score(f=freq, f1=f.PP, f2=f.verb, N=N, measure="MI") * log10(2),
  Dice = am.score(f=freq, f1=f.PP, f2=f.verb, N=N, measure="Dice"),
  z.score = am.score(f=freq, f1=f.PP, f2=f.verb, N=N, measure=my.z.score),
  t.score = am.score(f=freq, f1=f.PP, f2=f.verb, N=N, measure="t"),
  chisq = am.score(f=freq, f1=f.PP, f2=f.verb, N=N, measure=my.chisq),
  chisq.corr = am.score(f=freq, f1=f.PP, f2=f.verb, N=N, measure=my.chisq.corr),
  log.like = am.score(f=freq, f1=f.PP, f2=f.verb, N=N, measure="G2"),
  Fisher = am.score(f=freq, f1=f.PP, f2=f.verb, N=N, measure="Fisher.pv", p.adjust=FALSE)
)

## compare with reference implementation (should be highly accurate after all adjustements)
stopifnot(all.equal(PPV, KrennPPV, tol=1e-12))

## standard application of am.score() with lookup of marginal frequencies in separate tables
idx <- !duplicated(PPV$PP)
f1.tbl <- structure(PPV$f.PP[idx], names=PPV$PP[idx])
idx <- !duplicated(PPV$verb)
f2.tbl <- structure(PPV$f.verb[idx], names=PPV$verb[idx])
sample.size <- max(PPV$N)

P2 <- KrennPPV[, 1:9] # recompute AM scores with lookup in marginal tables
P2 <- transform(
  P2,
  MI = am.score(PP, verb, freq, f1.tbl, f2.tbl, sample.size, measure="MI") * log10(2),
  Dice = am.score(PP, verb, freq, f1.tbl, f2.tbl, sample.size, measure="Dice"),
  z.score = am.score(PP, verb, freq, f1.tbl, f2.tbl, sample.size, measure=my.z.score),
  t.score = am.score(PP, verb, freq, f1.tbl, f2.tbl, sample.size, measure="t"),
  chisq = am.score(PP, verb, freq, f1.tbl, f2.tbl, sample.size, measure=my.chisq),
  chisq.corr = am.score(PP, verb, freq, f1.tbl, f2.tbl, sample.size, measure=my.chisq.corr),
  log.like = am.score(PP, verb, freq, f1.tbl, f2.tbl, sample.size, measure="G2"),
  Fisher = am.score(PP, verb, freq, f1.tbl, f2.tbl, sample.size, measure="Fisher.pv", p.adjust=FALSE)
)

stopifnot(all.equal(P2, PPV, tol=1e-14)) # should be completely identical
