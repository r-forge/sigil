##    Author: Stefan Evert
##   Purpose: Utility functions used by several plotting scripts
##   Created: Sun Aug  6 18:42:17 2006
##  Modified: Fri Jun 27 12:07:48 2014 (evert)    

##
## compute type probabilities pi_k from LNRE model
##
lnre.pi.k <- function (M, k) {
  plnre(M, tqlnre(M, k-1)) - plnre(M, tqlnre(M, k))
}

##
## ZM or fZM model with specified parameters for the Zipf-Mandelbrot law
##
make.zm.model <- function (a, b, S=Inf) {
  alpha <- 1/a
  ## for the fZM, this is an approximation, but it works well as long as A is very small
  B <- (1 - alpha) / (b * alpha) 
  if (is.infinite(S)) {
    lnre("zm", alpha=alpha, B=B)
  }
  else {
    ## choose A so that vocabulary size S has specified value
    root.S <- function (A, S) {
      ((1-alpha)/alpha) * (A^(-alpha) - B^(-alpha)) / (B^(1-alpha) - A^(1-alpha)) - S
    }
    .res <- uniroot(root.S, lower=1e-30, upper=B*(1-1e-12), tol=1e-50, S=S)
     if (abs(.res$f.root) > 0.5) stop("parameter estimation failed")
    lnre("fzm", alpha=alpha, A=(.res$root), B=B)
  }
}

##
## plot Zipf ranking of ZM/fZM model with specified ZM parameters (showing paramter values)
##
test.zm.model <- function (a, b, S=Inf, log=FALSE) {
  model <- make.zm.model(a=a, b=b, S=S)
  print(model)
  m.values <- if (log) 1:500 else 1:50
  zipf.plot(lnre.pi.k(model, m.values), p.max=.1, log=log, legend=zm.legend(model))
  print("PLOT")
}

##
## compute parameters of Zipf-Mandelbrot law from LNRE model (formatted as expression for legends) 
##
zm.legend <- function (model, single.line=FALSE) {
  if (! (inherits(model, "lnre.zm") || inherits(model, "lnre.fzm")))
    stop("'model' must be a ZM or fZM LNRE model")
  alpha <- model$param$alpha
  B <- model$param$B
  a <- 1 / alpha
  if (inherits(model, "lnre.zm")) {     # ZM model
    b <- (1 - alpha) / (B * alpha)
    if (single.line) {
     bquote( list(a == .(round(a,2)), b == .(round(b,2))) )
    }
    else {
      eval(bquote(expression(
          a == .(round(a,2)), b == .(round(b,2))
          )))
    }
  }
  else {                                # fZM model
    S <- model$S
    C <- model$param2$C
    b <- C / (alpha * B^alpha)
    if (single.line) {
      bquote( list(a == .(round(a,2)), b == .(round(b,2)), S == .(round(S,1))) )
    }
    else {
      eval(bquote(expression(
          a == .(round(a,2)), b == .(round(b,2)), S == .(round(S,1))
          )))
    }
  }
}

##
## Zipf plot of type probabilites (predicted by LNRE models)
##
zipf.plot <- function (p, k.max=if (log) 500 else 50,
                       p.max=max(p)*1.1, p.min=1e-4,
                       xlab=expression(k), ylab=expression(pi[k]), main="",
                       log=FALSE, legend=NULL, file=NULL) {
  zipfR.begin.plot(filename=file)
  .xlim <- if (log) c(1, k.max) else c(0, k.max)
  .ylim <- if (log) c(p.min, p.max) else c(0, p.max)
  .logtype <- if (log) "xy" else ""
  plot(1:length(p), p, type="o", pch=20, cex=1.2,
       log=.logtype, xlim=.xlim, xaxs="i", ylim=.ylim, yaxs="i",
       xlab=xlab, ylab=ylab, main=main)
  if (!missing(legend)) legend("topright", inset=.05, legend=legend, cex=1.4)
  zipfR.end.plot()
}

##
## fake histogram plot, i.e boxes centered around equidistant data points
##
fake.hist <- function (x, y, width=mean(diff(x))/2, add=FALSE,
                       xlim=NULL, ylim=NULL, areas=FALSE,
                       col="black", fill=NA, lwd=1, ...)
{
  if (areas) y <- y / width
  if (missing(xlim)) xlim <- c(min(x)-width, max(x)+width)
  if (missing(ylim)) ylim <- c(0, max(y)*1.1)


  if (!add) {                           # set up new plot (axes, labels, titles)
    plot(0, 0, type="n", xlim=xlim, ylim=ylim, ...)
  }
  rect(x - width, 0, x + width, y, border=col, col=fill, lwd=lwd)
  box()    # make sure that histogram doesn't overpaint box around plot region
}

##
## -- compare two distributions by plotting kernel densities (with optional rug plots for small data sets) --
##  xlim, ylim: if not set, will automatically determine display window to fit data
##  rug:        if TRUE, add colour-coded rug plots to the density curves
##  legend:     optional legend (character vector of length 2)
##  paired:     if set, perform paired (TRUE) or unpaired (FALSE) t-test for location difference
##  info.cex:   text scaling factor for t-test information
##
compare.distributions <- function (x1, x2=numeric(0), rug=FALSE, paired=NA, info.cex=1, lwd=3, xlim=NULL, ylim=NULL, legend=NULL, main="", xaxs="i", yaxs="i", ...) {
	compare <- !missing(x2)
	if (missing(xlim)) xlim <- extendrange(c(x1, x2))
	d1 <- density(x1, from=xlim[1], to=xlim[2])
	if (compare) d2 <- density(x2, from=xlim[1], to=xlim[2]) else d2 <- d1
	if (missing(ylim)) ylim <- c(0, max(c(d1$y, d2$y)) * 1.05)
	plot(d1, type="l", col="red", lwd=lwd, main=main, xlim=xlim, ylim=ylim, xaxs=xaxs, yaxs=yaxs, ...)
	if (compare) lines(d2, col="blue", lwd=lwd, ...)
	if (rug) {
		rug(x1, ticksize=.03, lwd=.5, col="red")
		if (compare) rug(x2, ticksize=-.03, lwd=.5, col="blue")
	}
	if (!is.null(legend)) legend("topright", inset=.02, bg="white", legend=legend, lwd=lwd+1, col=c("red", "blue"))
	if (compare && !missing(paired)) {
		result <- t.test(x1, x2, paired=paired)
		print(result)
		legend("topleft", inset=c(0, .02), bty="n", cex=info.cex, legend=c(if (paired) "paired t-test" else "two-sample t-test", sprintf("p = %.5f", result$p.value), sprintf("t = %.1f", result$statistic)))
	}
}

##
## -- print legend with info about ZM/fZM model fit --
##   position:   "left", "top", ...
##   param:      if TRUE, print parameters alpha and S
##   gof:        if TRUE, print goodness-of-fit information
##   cex:        additional font scaling (for readability)
##
lnre.info <- function(position, model, param=FALSE, gof=TRUE, cex=1) {
	info <- expression()
	if (param) {
		info <- append(info, bquote(alpha == .(round(model$param$alpha,5))))
		info <- append(info, bquote(S == .(round(model$S))))
	}
	if (param && gof) info <- append(info, bquote(""))
	if (gof) {
		info <- append(info, bquote(X^2 == .(round(model$gof$X2, 2))))
		info <- append(info, bquote(df == .(model$gof$d)))
		info <- append(info, bquote(p == .(signif(model$gof$p, 4))))
	}
	legend(position, inset=.02, legend=info, cex=cex, bty="n")
}
