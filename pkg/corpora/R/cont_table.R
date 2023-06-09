cont.table <- function (k1, n1, k2, n2, as.list=NA) {
  l <- .match.len(c("k1", "n1", "k2", "n2"), adjust=TRUE) # ensure that all vectors have the same length
  if (missing(as.list)) as.list <- if (l > 1) TRUE else FALSE
  if (l > 1 && !as.list) stop("k1, n1, k2, n2 must be single numbers (as.list=FALSE)")

  # sanity checks
  if (any(k1 < 0) || any(k1 > n1) || any(n1 <= 0)) stop("k1 and n1 must be integers with 0 <= k1 <= n1")
  if (any(k2 < 0) || any(k2 > n2) || any(n2 <= 0)) stop("k2 and n2 must be integers with 0 <= k2 <= n2")
  if (any(k1 + k2 <= 0)) stop("either k1 or k2 must be non-zero")

  # construct list of 2x2 contingency tables
  table.list <- lapply(1:l, function (i) matrix(c(k1[i], n1[i]-k1[i], k2[i], n2[i]-k2[i]), nrow=2, byrow=FALSE))
  
  if (as.list) table.list else table.list[[1]]
}
