##
## Illustration of theoretical sampling distributions
##

screen.device <- quartz # change this if not working on a Mac
blue1 <- "#4477AA"
blue2 <- "#B0D4FF"
red1 <- "#DD4433"
red2 <- "#FF9980"
green1 <- "#229966"

library(corpora)

##
## 1) Sampling distribution for sample size n=100 and null proportion pi=15% of passive sentences
k <- 0:30
n <- 100
p <- .15
perc.k <- round(100 * dbinom(k, n, p), digits=1)
bars.1 <- rbind(perc.k * (k != 19), perc.k * (k == 19))  # observed value k=19
bars.2 <- rbind(perc.k * (k < 19), perc.k * (k >= 19))   # upper tail (k >= 19)
bars.3 <- rbind(perc.k * (k > 10 & k < 19), perc.k * (k <= 10 | k >= 19))   # both tails (k >= 19 | k <= 10)

screen.device(width=11, height=6, bg="white")
par(cex=1.5, mar=c(4,4,1,0)+.5)

mids <- barplot(bars.1, names=ifelse(k%%2==1, k, NA), ylim=c(0,12),
                col=c("grey", "black"), density=c(NA, 10), angle=45, 
                xlab="value k of observed frequency X", ylab="percentage of samples with X=k")
text(mids, perc.k + .5, perc.k, cex=.7)
dev.copy2pdf(file="../keynote-slides/img/binomial_sampling_dist_1.pdf", onefile=FALSE) # should adjust width/height automatically

mids <- barplot(bars.2, names=ifelse(k%%2==1, k, NA), ylim=c(0,12),
                col=c("grey", "black"), density=c(NA, 10), angle=45, 
                xlab="value k of observed frequency X", ylab="percentage of samples with X=k")
text(mids, perc.k + .5, perc.k, cex=.7)
dev.copy2pdf(file="../keynote-slides/img/binomial_sampling_dist_2.pdf", onefile=FALSE)

mids <- barplot(bars.3, names=ifelse(k%%2==1, k, NA), ylim=c(0,12),
                col=c("grey", "black"), density=c(NA, 10), angle=45, 
                xlab="value k of observed frequency X", ylab="percentage of samples with X=k")
text(mids, perc.k + .5, perc.k, cex=.7)
dev.copy2pdf(file="../keynote-slides/img/binomial_sampling_dist_3.pdf", onefile=FALSE)

## colour version for presentation
mids <- barplot(bars.1, names=ifelse(k%%2==1, k, NA), ylim=c(0,12), col=c("grey", blue1), 
                xlab="value k of observed frequency X", ylab="percentage of samples with X=k")
text(mids, perc.k + .5, perc.k, cex=.7)
dev.copy2pdf(file="../keynote-slides/img/binomial_sampling_dist_c1.pdf", onefile=FALSE) # should adjust width/height automatically

mids <- barplot(bars.2, names=ifelse(k%%2==1, k, NA), ylim=c(0,12), col=c("grey", blue1), 
                xlab="value k of observed frequency X", ylab="percentage of samples with X=k")
text(mids, perc.k + .5, perc.k, cex=.7)
dev.copy2pdf(file="../keynote-slides/img/binomial_sampling_dist_c2.pdf", onefile=FALSE)

mids <- barplot(bars.3, names=ifelse(k%%2==1, k, NA), ylim=c(0,12), col=c("grey", blue1), 
                xlab="value k of observed frequency X", ylab="percentage of samples with X=k")
text(mids, perc.k + .5, perc.k, cex=.7)
dev.copy2pdf(file="../keynote-slides/img/binomial_sampling_dist_c3.pdf", onefile=FALSE)

dev.off()


##
## 2) Effect size and power: comparison of barplots (for one-sided test)
two.binom <- function(n, p1, p2=0, range=c(0, 40), y.max=NA, threshold=NA, effect=FALSE) {
  e1 <- n * p1
  e2 <- n * p2
  k <- round(seq(range[1], range[2], 1))
  d1 <- dbinom(k, n, p1)
  d2 <- dbinom(k, n, p2)
  perc1 <- 100 * d1
  perc2 <- 100 * d2
  if (is.na(y.max)) y.max <- max(perc1, perc2)
  if (is.na(threshold)) {
    cols <- c(blue2, red2)
  } else {
    cols1 <- ifelse(k >= threshold, blue1, blue2)
    cols2 <- ifelse(k < threshold, red1, red2)
    cols <- rbind(cols1, cols2)
  }
  l <- diff(range)
  step <- if (l >= 300) 50 else if (l >= 100) 25 else if (l >= 30) 5 else 2
  labels <- ifelse(k %% step == 0, k, NA)
  mids <- barplot(rbind(perc1, perc2), beside=TRUE, names.arg=labels, col=cols, border=cols, space=c(0, .5), ylim=c(0, y.max),
                  xlab="value k of observed frequency X", ylab="percentage of samples with X=k")
  x2scr <- function (x) (x - range[1]) * diff(range(mids)) / diff(range) + mids[1]  # convert x-axis to screen coordinates 
  if (effect) {
    abline(v=x2scr(e1), col=blue1, lwd=5, lty="22")
    text(x2scr(e1), y.max, sprintf("%.0f", e1), adj=c(1, -0.5), srt=90, cex=.8, font=2, col=blue1)
    if (p2 > 0) {
      abline(v=x2scr(e2), col=red1, lwd=5, lty="22")
      text(x2scr(e2), y.max, sprintf("%.0f", e2), adj=c(1, 1.5), srt=90, cex=.8, font=2, col=red1)
      arrows(x2scr(e1), y.max * .90, x2scr(e2), y.max * .90, code=3, col=green1, lwd=2, length=.15)
      text(x2scr((e1 + e2) / 2), y.max * .90, bquote(delta == .(sprintf("%.2f", p2 - p1))), pos=3, col=green1)
    }
  }
  if (!is.na(threshold)) {
    risk1 <- 100 * pbinom(threshold - 1, n, p1, lower.tail=FALSE)
    text(max(mids), y.max, sprintf("type I risk  = %4.1f%%", risk1), adj=c(1, 1.2), family="mono", col=blue1)
    if (p2 > 0) {
      risk2 <- 100 * pbinom(threshold - 1, n, p2, lower.tail=TRUE)
      text(max(mids), y.max, sprintf("type II risk = %4.1f%%", risk2), adj=c(1, 2.8), family="mono", col=red1)
    }
  }
}

screen.device(width=11, height=6, bg="white")
par(cex=1.4, mar=c(4,4,1,0)+.5)

two.binom(100, .15, range=c(5, 40), y.max=16)
dev.copy2pdf(file="../keynote-slides/img/binomial_effect_power_100_c1.pdf", onefile=FALSE)

two.binom(100, .15, range=c(5, 40), y.max=16, threshold=22)
dev.copy2pdf(file="../keynote-slides/img/binomial_effect_power_100_c2.pdf", onefile=FALSE)

two.binom(100, .15, .25, range=c(5, 40), y.max=16, threshold=22)
dev.copy2pdf(file="../keynote-slides/img/binomial_effect_power_100_c3.pdf", onefile=FALSE)

two.binom(100, .15, .25, range=c(5, 40), y.max=16, threshold=22, effect=TRUE)
dev.copy2pdf(file="../keynote-slides/img/binomial_effect_power_100_c4.pdf", onefile=FALSE)

two.binom(100, .15, .20, range=c(5, 40), y.max=16, threshold=22, effect=TRUE)
dev.copy2pdf(file="../keynote-slides/img/binomial_effect_power_100_c5.pdf", onefile=FALSE)

two.binom(100, .15, .30, range=c(5, 40), y.max=16, threshold=22, effect=TRUE)
dev.copy2pdf(file="../keynote-slides/img/binomial_effect_power_100_c6.pdf", onefile=FALSE)

two.binom(1000, .15, .20, range=c(50, 400), y.max=4, threshold=170, effect=TRUE)
dev.copy2pdf(file="../keynote-slides/img/binomial_effect_power_1000_c5a.pdf", onefile=FALSE)

two.binom(1000, .15, .20, range=c(100, 300), y.max=4, threshold=170, effect=TRUE)
dev.copy2pdf(file="../keynote-slides/img/binomial_effect_power_1000_c5b.pdf", onefile=FALSE)

two.binom(1000, .15, .25, range=c(50, 400), y.max=4, threshold=170, effect=TRUE)
dev.copy2pdf(file="../keynote-slides/img/binomial_effect_power_1000_c4a.pdf", onefile=FALSE)

two.binom(1000, .15, .25, range=c(100, 300), y.max=4, threshold=170, effect=TRUE)
dev.copy2pdf(file="../keynote-slides/img/binomial_effect_power_1000_c4b.pdf", onefile=FALSE)

dev.off()


##
## 3) Power analysis of binomial test
x.max <- 0.20
delta <- seq(0, x.max, .001)

power.curve <- function (delta, n, alpha=.05, p0=.15, ...) {
  threshold <- qbinom(alpha, n, p0, lower.tail=FALSE) # reject if k > threshold
  risk <- 100 * pbinom(threshold, n, p0 + delta)
  lines(delta, risk, ...)
}

n <- rep(c(100, 1000, 1e5), each=3)
alpha <- rep(c(.05, .01, .001), 3)
col <- rep(c(red1, blue1, green1), each=3)
lty <- rep(c("solid", "solid", "dashed"), 3)
lwd <- rep(c(4, 2, 2), 3)
legend.txt <- expression(
  list(n == 100, alpha == .05),
  list(n == 100, alpha == .01),
  list(n == 100, alpha == .001),
  list(n == 1000, alpha == .05),
  list(n == 1000, alpha == .01),
  list(n == 1000, alpha == .001),
  list(n == "100k", alpha == .05),
  list(n == "100k", alpha == .01),
  list(n == "100k", alpha == .001)
  )

power.plot <- function (delta, k) {
  plot(0, 0, type="n", xlim=c(0, x.max), ylim=c(0, 100), xlab=expression("effect size " * (delta)), ylab="risk of type II error (%)\nfor one-sided binomial test", main="")
  abline(h=seq(0, 100, 10), lwd=1, col="grey")
  for (i in 1:k) {
    power.curve(delta, n=n[i], alpha=alpha[i], col=col[i], lty=lty[i], lwd=lwd[i])
  }
  legend("topright", inset=.03, legend=legend.txt[1:k], col=col[1:k], lwd=lwd[1:k]+1, lty=lty[1:k])
}
  
screen.device(width=11, height=6, bg="white")
par(cex=1.4, mar=c(4,5,1,1)+.5, xaxs="i", yaxs="i")

power.plot(delta, 1)
dev.copy2pdf(file="../keynote-slides/img/binomial_power_curve_c1.pdf", onefile=FALSE)

power.plot(delta, 2)
dev.copy2pdf(file="../keynote-slides/img/binomial_power_curve_c2.pdf", onefile=FALSE)

power.plot(delta, 3)
dev.copy2pdf(file="../keynote-slides/img/binomial_power_curve_c3.pdf", onefile=FALSE)

power.plot(delta, 6)
dev.copy2pdf(file="../keynote-slides/img/binomial_power_curve_c4.pdf", onefile=FALSE)

power.plot(delta, 9)
dev.copy2pdf(file="../keynote-slides/img/binomial_power_curve_c5.pdf", onefile=FALSE)

dev.off()


##
## 4) Illustration of procedure for computing confidence intervals (by inverted hypothesis tests)
visual.conf.plot <- function (p, n=1000, f=190, xlim=c(150,250), ylim=NULL, alpha=.05, mark=NULL) {
  e <- p * n                            # expected frequency
  k <- 0:n                              # range of possible values for binomial variable
  prob.binom <- dbinom(k, n, p)         # vector of binomial probabilities
                                        # idx <- order((k - e)^2)               # re-order possible observations k according to their chi-squared value
  ## above carries out an exact chi-squared test; but binom.test() uses likelihood criterion:
  idx <- order(prob.binom, decreasing=TRUE) # order possible observations according to their likelihood
  tail.binom <- numeric(length(k))
  tail.binom[idx] <- rev( cumsum( prob.binom[rev(idx)] ) ) # accurate two-sided chi-squared tail probabilities
  perc.binom <- 100 * prob.binom                           # show percentages on y-axis of plot
  if (missing(ylim)) ylim <- c(0, max(perc.binom) * 1.1)   # calculate y-axis limits automatically unless specified
  reject <- tail.binom[f+1] < alpha
  text.reject <- if (reject) "rejected" else "plausible" # whether H0 can be rejected (shown in title of plot)
  hist.col <- ifelse(tail.binom < alpha, blue1, blue2) # colours for histogram bars (two-sided alpha tail shown in black)
  plot(0, 0, type="n",
       xlim=xlim, ylim=ylim, yaxs="i",
       xlab="observed frequency k", ylab="percentage of samples",
       main=bquote(H[0]: mu == .(100 * p)*"%" %->% .(text.reject)))
  if (!missing(mark)) {
    mark.pval <- sapply(mark, function (mu) binom.test(f, n, mu)$p.value)
    mark.col <- ifelse(mark.pval < alpha, red1, green1)
    abline(v = n * mark, lwd=1, col=mark.col)
  }
  points(k, perc.binom, type="h", lwd=5, col=hist.col)
  col.reject <- if (reject) red1 else green1
  abline(v=f, lwd=4, lty="32", col=col.reject)
  text(f, ylim[2] * .95, pos=2, cex=1.2, col=col.reject, labels=bquote(f == .(f)), srt=90)
}
visual.conf.plot(.216, mark=c(.16, .165, .17))

screen.device(width=11, height=6, bg="white")
par(cex=1.4, mar=c(4,4,2,0)+.5)
# par(cex=1.3, mar=c(4,4,2,0)+.1)

mu.vals <- c(.16, .165, .166, .167, .17, .18, .19, .20, .21, .215, .216, .22, .23, .24)
for (mu in mu.vals) {
  visual.conf.plot(mu, mark=seq(.08, mu, .0005), ylim=c(0, 5))
  dev.copy2pdf(file=sprintf("../keynote-slides/img/confidence_interval_%03d.pdf", 1000 * mu), onefile=FALSE, bg="white")
}

dev.off()


##
## 5) Choosing sample size based on graph of confidence intervals
draw <- function (K, col="black", lty="solid", lwd=3) {
  p <- 100 * (0:K) / K
  lines(p, 100 * prop.cint(0:K, K)$upper, lwd=lwd, col=col, lty=lty)
  lines(p, 100 * prop.cint(0:K, K)$lower, lwd=lwd, col=col, lty=lty)
}
draw.ppm <- function (K, col="black", lty="solid", lwd=3) {
  p <- seq(0, 20, .1)
  n <- K * 1e6
  O <- K * p
  cint <- prop.cint(O, n)
  lines(p, 1e6 * cint$upper, lwd=lwd, col=col, lty=lty)
  lines(p, 1e6 * cint$lower, lwd=lwd, col=col, lty=lty)
}

screen.device(width=8, height=8, bg="white")
par(cex=1.4, mar=c(4,4,1,1)+.1, xaxs="i", yaxs="i")

plot(c(0,100), c(0,100), type="l", lwd=3, xlab="sample: k / n (%)", ylab=expression("population: " * pi * " (%)"), main="")
for (x in 0:10 * 10) { abline(h=x); abline(v=x) }
draw(500, col="orange")
draw(200, col="darkgreen")
draw(100, col="red")
draw(50, col="blue")
draw(20, col="grey")
legend("topleft", inset=.05, bg="white",
       legend=c("MLE", "n = 500", "n = 200", "n = 100", "n = 50", "n = 20"),
       lwd=4, col=c("black", "orange", "darkgreen", "red", "blue", "grey"))
dev.copy2pdf(file="../keynote-slides/img/sample_size_1.pdf", onefile=FALSE, bg="white")

plot(c(0,25), c(0,25), type="l", lwd=3, xlab="sample: k / n (%)", ylab=expression("population: " * pi * " (%)"), main="")
for (x in 0:10 * 10) { abline(h=x); abline(v=x) }
for (x in 0:4 * 5) { abline(h=x,lty="dashed"); abline(v=x,lty="dashed") }
draw(500, col="orange")
draw(200, col="darkgreen")
draw(100, col="red")
draw(50, col="blue")
draw(20, col="grey")
legend("topleft", inset=.05, bg="white",
       legend=c("MLE", "n = 500", "n = 200", "n = 100", "n = 50", "n = 20"),
       lwd=4, col=c("black", "orange", "darkgreen", "red", "blue", "grey"))
dev.copy2pdf(file="../keynote-slides/img/sample_size_2.pdf", onefile=FALSE, bg="white")

plot(c(0,20), c(0,20), type="l", lwd=3, xlab="sample: k / n (ppm)", ylab=expression("population: " * pi * " (ppm)"), main="")
for (x in 0:10 * 10) { abline(h=x); abline(v=x) }
for (x in 0:4 * 5) { abline(h=x,lty="dashed"); abline(v=x,lty="dashed") }
draw.ppm(100, col="orange")
draw.ppm(10, col="darkgreen")
draw.ppm(5, col="red")
draw.ppm(1, col="blue")
draw.ppm(.5, col="grey")
legend("topright", inset=.05, bg="white",
       legend=c("MLE", "n = 100M", "n = 10M", "n = 5M", "n = 1M", "n = 500k"),
       lwd=4, col=c("black", "orange", "darkgreen", "red", "blue", "grey"))
dev.copy2pdf(file="../keynote-slides/img/sample_size_3.pdf", onefile=FALSE, bg="white")
