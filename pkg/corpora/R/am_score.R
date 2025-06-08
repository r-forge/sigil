am.score <- function(w1, w2, f, f1, f2, N, measure,
                    span.size=1, log=FALSE, labels=FALSE, conf.level=.95, p.adjust=TRUE, param=list()) {
  if (log && (missing(w1) || missing(w2))) stop("w1 and w2 must be specified if log=TRUE")
  if (!is.function(measure)) {
    measure <- match.arg(measure, names(builtin.am))
    measure <- builtin.am[[measure]]
  }
  
  l <- length(f) # ensure that all vectors have the same common length (determined by f)
  if (!missing(w1)) .match.len("w1", len=l, adjust=TRUE, check.numeric=FALSE)
  if (!missing(w2)) .match.len("w2", len=l, adjust=TRUE, check.numeric=FALSE)
  .match.len("N", len=l, adjust=TRUE)
  
  if (is.null(names(f1))) {
    .match.len("f1", len=l)
  } else {
    f1 <- f1[w1]
    if (any(is.na(f1))) stop("f1 must provide marginal frequencies for all distinct strings in w1")
  }

  if (is.null(names(f2))) {
    .match.len("f2", len=l)
  } else {
    f2 <- f2[w2]
    if (any(is.na(f2))) stop("f2 must provide marginal frequencies for all distinct strings in w2")
  }
  
  if (length(conf.level) != 1) stop("conf.level must be a scalar") # validate other arguments
  if (any(conf.level <= 0) || any(conf.level > 1)) stop("conf.level must be in range [0,1]")
  
  ## Bonferroni correction factor (1 for p.adjust=FALSE)
  if (isFALSE(p.adjust)) {
    p.adjust <- 1
  }
  else if (isTRUE(p.adjust)) {
    p.adjust <- l
  }
  else {
    if (!(is.numeric(p.adjust) && length(p.adjust) == 1)) stop("p.adjust must either be TRUE/FALSE or a number specifying the family size")
  }

  ## ensure that all frequency data are floating-point (double), in particular guarding against bit64::integer64 vectors
  .ensure.double(c("f", "f1", "f2", "N"))
  
  ## apply span size adjustment (guarding against R1 > N and inconvenient R2 = 0)
  if (length(span.size) != 1 || any(span.size < 1)) stop("span.size must be an integer >= 1")
  if (span.size != 1) f1 <- pmin(f1 * span.size, N - 1)
  
  ## add implicit parameters
  if (!is.list(param) || is.null(names(param))) stop("param must be a named list")
  param$conf.level <- conf.level
  param$p.adjust <- p.adjust # effective family size m
  
  ## wrapper function to compute observed and expected frequencies as needed
  compute.AM <- function (
    AM, f, f1, f2, N, param,
    O=f, E=f1*f2/N,
    R1=f1, R2=N-f1, C1=f2, C2=N-f2,
    O12=f1-f, O21=f2-f, O22=N-f1-f2+f,
    E12=f1*C2/N, E21=R2*f2/N, E22=R2*C2/N) {
    AM(f=f, f1=f1, f2=f2, N=N, param=param,
       O=O, E=E, R1=R1, R2=R2, C1=C1, C2=C2,
       O11=O, O12=O12, O21=O21, O22=O22, E11=E, E12=E12, E21=E21, E22=E22)
  }
  
  ## compute AM scores and apply optional log transform
  scores <- compute.AM(measure, f, f1, f2, N, param)
  if (log) scores <- sign(scores) * log2(abs(scores) + 1)
  
  if (labels) names(scores) <- paste(w1, w2)
  scores
}


######################################################################
## implementations of built-in association measures (exported for reference)

builtin.am <- list(
  MI = function (O, E, ...) {
    log2(O / E)
  },
  MI.k = function (O, E, param, ...) {
    k <- if ("k" %in% names(param)) param$k else 2
    log2(O^k / E)
  },
  G2 = function (O11, O12, O21, O22, E11, E12, E21, E22, ...) {
    .G2.term <- function (o, e) {
      res <- o * log(o / e)
      res[o <= 0] <- 0
      res
    }

    G2 <- .G2.term(O11, E11) + .G2.term(O12, E12) + .G2.term(O21, E21) + .G2.term(O22, E22)
    sign(O11 - E11) * 2 * G2
  },
  G2.pv = function (O11, O12, O21, O22, E11, E12, E21, E22, param, ...) {
    .G2.term <- function (o, e) {
      res <- o * log(o / e)
      res[o <= 0] <- 0
      res
    }
    m <- param$p.adjust # family size for Bonferroni correction
    
    G2 <- .G2.term(O11, E11) + .G2.term(O12, E12) + .G2.term(O21, E21) + .G2.term(O22, E22)
    pv <- pchisq(2 * G2, df=1, lower.tail=FALSE, log.p=TRUE)
    if (m > 1) pv <- pmin(pv + log(m), 0) # Bonferroni: p = p * m 
    sign(O11 - E11) * (-pv / log(10)) # convert to -log10(p), add sign for negative assocation
  },
  simple.ll = function (O, E, ...) {
    term <- O * log(O / E)
    term[O <= 0] <- 0
    sign(O - E) * 2 * (term + (O - E))
  },
  t = function (O, E, ...) {
    (O - E) / sqrt(O)
  },
  X2 = function (O11, O12, O21, O22, R1, R2, C1, C2, N, ...) {
    ## common form for homogeneity test with Yates' correction (Evert 2004: 82)
    term <- abs(O11 * O22 - O12 * O21)
    term <- pmax(term - N/2, 0) # Yates' correction
    X2 <- (N * term^2) / (R1 * R2 * C1 * C2)
    sign(O11 - E11) * X2
  },
  z = function (O, E, ...) {
    ## z-score with Yates' correction (and smooth transition to O-E = 0)
    .yates.corr <- function (x) {
      x.abs <- abs(x)
      sign(x) * ifelse(x.abs >= 1, x.abs - 0.5, x.abs / 2)
    }
    .yates.corr(O - E) / sqrt(E)
  },
  Dice = function (O11, R1, C1, ...) {
    2 * O11 / (R1 + C1)
  },
  DP = function (O11, O21, R1, R2, ...) {
    O11 / R1 - O21 / R2
  },
  LRC = function (O11, O21, R1, R2, param, ...) {
    keyness(O11, R1, O21, R2, "PositiveLRC", conf.level=param$conf.level, p.adjust=param$p.adjust)
  }
)
