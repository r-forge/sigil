##
## Illustration of hypothesis test and effect size measures for frequency comparison
##

screen.device <- quartz # change this if not working on a Mac

##
## 1) visualize effect size measures with contour plots
eff.dp <- function (pi1, pi2, ...) {
  pi1 - pi2 # difference of proportions
}
eff.r <- function (pi1, pi2, log=FALSE, ...) {
  if (log) log2(pi1 / pi2) else pi1 / pi2 # relative risk
}
eff.or <- function (pi1, pi2, log=FALSE, ...) {
  or <- (pi1 / (1 - pi1)) / (pi2 / (1 - pi2))
  if (log) log2(or) else or # relative risk
}
eff.phi <- function (pi1, pi2, r=c(1, 1), abs=FALSE, ...) {
  r <- r / sum(r) # normalize relative sample sizes
  pi <- pi1 * r[1] + pi2 * r[2]
  denom <- sqrt( pi * (1 - pi) / (r[1] * r[2]) )
  numer <- pi1 * (1 - pi2) - pi2 * (1 - pi1)
  if (abs) abs(numer) / denom else numer / denom
}

## test that population based phi-coefficient agrees with definition based on X2 statistic
test.phi <- function (ct = NULL) {
  if (is.null(ct)) ct <- matrix(round(runif(4, 1, 50)), 2, 2)
  n <- colSums(ct) # sample sizes
  V <- sqrt(chisq.test(ct, correct=FALSE)$statistic / sum(n))
  phi <- eff.phi(ct[1,1] / n[1], ct[1,2] / n[2], r=n)
  print(ct)
  cat(sprintf("V = %+7.4f   phi = %+7.4f\n", V, phi))
  invisible(ct)
}
for (i in 1:5) suppressWarnings(test.phi())

## visualize effect size measure with a contour plot
effect.contour <- function (FUN, pi.max=1, grid=100, zlim=NULL, nlevels=20, main="", ...) {
  pi.vec <- pi.max * (1:(grid-1)) / grid # exclude extremes to avoid invalid effect sizes
  eff <- outer(pi.vec, pi.vec, FUN, ...) # evaluate effect size measure on grid
  str(pi.vec)
  str(eff)
  col.range <- c(hsv(h=2/3, v=1, s=(100:1)/100), rgb(1,1,1), hsv(h=0, v=1, s=(1:100/100))) 
  if (is.null(zlim)) zlim <- c(-1,1) * max(abs(eff))
  image(pi.vec, pi.vec, t(eff), xlim=c(0, pi.max), ylim=c(0, pi.max), zlim=zlim,
        col=cor.colors, xlab=expression(pi[2]), ylab=expression(pi[1]), main=main)
  contour(pi.vec, pi.vec, t(eff), nlevels=nlevels, labcex=1, add=TRUE)
}

screen.device(width=7, height=7, bg="white")
par(cex=1.2, mar=c(4,4,2,2)+.2)

## visualize effect size measures for full range (0, 1) of proportions 
effect.contour(eff.dp, zlim=c(-1,1), main=expression(pi[1] - pi[2]))
dev.copy2pdf(file="../keynote-slides/img/effect_size_difference.pdf")

effect.contour(eff.r, log=TRUE, main=expression(paste("relative risk: ", log[2](r))))
dev.copy2pdf(file="../keynote-slides/img/effect_size_rrisk.pdf")

effect.contour(eff.or, log=TRUE, main=expression(paste("odds ratio: ", log[2](theta))))
dev.copy2pdf(file="../keynote-slides/img/effect_size_odds.pdf")

effect.contour(eff.phi, zlim=c(-1,1), main=expression(paste(phi, "-coefficient (equal samples)")))
dev.copy2pdf(file="../keynote-slides/img/effect_size_phi_1_1.pdf")

effect.contour(eff.phi, abs=TRUE, nlevels=10, zlim=c(-1,1), main=expression(paste(sqrt(X^2 / n), " (equal samples)")))
dev.copy2pdf(file="../keynote-slides/img/effect_size_phiabs_1_1.pdf")

effect.contour(eff.phi, r=c(10,1), zlim=c(-1,1), main=expression(paste(phi, "-coefficient (sample sizes 10 : 1)")))
dev.copy2pdf(file="../keynote-slides/img/effect_size_phi_10_1.pdf")

effect.contour(eff.phi, r=c(1,10), zlim=c(-1,1), main=expression(paste(phi, "-coefficient (sample sizes 1 : 10)")))
dev.copy2pdf(file="../keynote-slides/img/effect_size_phi_1_10.pdf")

## zoom in on small probabilities (0, .01)
effect.contour(eff.dp, pi.max=.01, zlim=c(-1, 1), main=expression(pi[1] - pi[2]))
dev.copy2pdf(file="../keynote-slides/img/effect_size_zoom_difference.pdf")

effect.contour(eff.r, log=TRUE, pi.max=.01, main=expression(paste("relative risk: ", log[2](r))))
dev.copy2pdf(file="../keynote-slides/img/effect_size_zoom_rrisk.pdf")

effect.contour(eff.or, log=TRUE, pi.max=.01, main=expression(paste("odds ratio: ", log[2](theta))))
dev.copy2pdf(file="../keynote-slides/img/effect_size_zoom_odds.pdf")

effect.contour(eff.phi, pi.max=.01, zlim=c(-1,1), main=expression(paste(phi, "-coefficient (equal samples)")))
dev.copy2pdf(file="../keynote-slides/img/effect_size_zoom_phi_1_1.pdf")

effect.contour(eff.phi, r=c(10,1), pi.max=.01, zlim=c(-1,1), main=expression(paste(phi, "-coefficient (sample sizes 10 : 1)")))
dev.copy2pdf(file="../keynote-slides/img/effect_size_zoom_phi_10_1.pdf")

effect.contour(eff.phi, r=c(1,10), pi.max=.01, zlim=c(-1,1), main=expression(paste(phi, "-coefficient (sample sizes 1 : 10)")))
dev.copy2pdf(file="../keynote-slides/img/effect_size_zoom_phi_1_10.pdf")



