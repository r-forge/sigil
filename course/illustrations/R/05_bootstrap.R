##
## Confidence intervals for LNRE models from parametric bootstrapping
##

library(zipfR)
zipfR.pick.device() # call with command-line option --pdf to generate PDF files (or pass "--pdf" here)

zipfR.par(init.par=list(mar=c(4, 4, 1, 1)+.1, cex=1.3), width=8, height=5)

## additional utility and plotting functions
source("R/05_utilities.R")

## bootstrapping validation for prefix ultra-
ItaUltra.spc # small sample size of 3467 tokens --> parameter estimates unreliable
N0 <- N(ItaUltra.spc)

## estimate fZM model
model <- lnre("fzm", ItaUltra.spc)
model # reasonable goodness of fit

## parametric bootstrap with 100 replicates
extract.stats <- function (m) data.frame(alpha=m$param$alpha, A=m$param$A, B=m$param$B, S=m$S, X2=m$gof$X2)
runs <- lnre.bootstrap(model, N0, lnre, extract.stats, type="fzm")
runs <- do.call(rbind, runs) # combine into data frame

zipfR.begin.plot(file="../latex-slides/img/05-bootstrap-ultra-alpha")
compare.distributions(runs$alpha, rug=TRUE, xlim=c(0, 1), xlab=expression(alpha))
abline(v=model$param$alpha, lwd=2)

zipfR.begin.plot(file="../latex-slides/img/05-bootstrap-ultra-A")
compare.distributions(runs$A, rug=TRUE, xlim=c(0, 2e-5), xlab=expression(A))
abline(v=model$param$A, lwd=2)

zipfR.begin.plot(file="../latex-slides/img/05-bootstrap-ultra-B")
compare.distributions(runs$B, rug=TRUE, xlim=c(0, 1), xlab=expression(B))
abline(v=model$param$B, lwd=2)

zipfR.begin.plot(file="../latex-slides/img/05-bootstrap-ultra-X2")
compare.distributions(runs$X2, rug=TRUE, xlim=c(0, 30), xlab=expression(X^2))
critical <- qchisq(.05, df=model$gof$df, lower=FALSE)
abline(v=critical, lwd=2, lty="dashed", col="blue")
text(critical, 0, expression(p < .05), srt=90, adj=c(-0.2, 1.4), col="blue")

zipfR.begin.plot(file="../latex-slides/img/05-bootstrap-ultra-S")
compare.distributions(runs$S, rug=TRUE, xlab=expression(S))

zipfR.begin.plot(file="../latex-slides/img/05-bootstrap-ultra-S-zoomed")
compare.distributions(runs$S, rug=TRUE, xlim=c(0, 100e3), xlab=expression(S))

zipfR.end.plot(shutdown=TRUE)

