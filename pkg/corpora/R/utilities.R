# validate that vector arguments have same length (or are scalars) 
#  - if adjust=TRUE replicate scalars to common length (NB: directly modifies the specified variables in the calling frame)
#  - checks that all vectors are numeric (unlesss check.numeric=FALSE)
#  - len= can be explicitly specified, otherwise inferred from arguments
#  - invisibly returns len
.match.len <- function (vars, len=NULL, adjust=FALSE, check.numeric=TRUE, envir=parent.frame()) {
  vecs <- setNames(lapply(vars, get, envir=envir), vars)
  ok <- sapply(vecs, is.numeric)
  if (check.numeric && any(!ok)) stop("argument(s) ", paste(vars[!ok], collapse=", "), " must be numeric vector(s)")
  if (is.null(len)) len <- max(sapply(vecs, length))
  for (v in vars) {
    if (length(vecs[[v]]) == 1) {
      if (adjust) assign(v, rep(vecs[[v]], len), envir=envir)
    }
    else if (length(vecs[[v]]) != len) {
      stop(sprintf("argument %s should be of length %d or a scalar (%s must have same length)", v, len, paste(vars, collapse=", ")))
    }
  }
  invisible(len)
}

# enforce double-precision floating-point format for vector (throwing error if it isn't numeric)
#  - operates on variables in parent environment like .match.len() does, for convenience and to avoid copies
#  - tries to detect bit64::integer64 and other edge cases
#  - dimensions and names are not preserved since conversion is performed with as.do
.ensure.double <- function (vars, envir=parent.frame()) {
  vecs <- setNames(lapply(vars, get, envir=envir), vars)
  ok <- sapply(vecs, is.numeric)
  if (any(!ok)) stop("argument(s) ", paste(vars[!ok], collapse=", "), " must be numeric vector(s)")
  for (v in vars) {
    vec <- vecs[[v]]
    if (storage.mode(vec) != "double") {
      storage.mode(vec) <- "double"
      assign(v, vec, envir=envir)
    } 
    else if ("integer64" %in% class(vec)) {
      assign(v, as.double(vec), envir=envir)
    }
  }
}
