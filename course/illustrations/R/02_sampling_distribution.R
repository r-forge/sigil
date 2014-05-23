##
## Illustration of theoretical sampling distributions
##

screen.device <- quartz # change this if not working on a Mac
blue1 <- "#4477AA"
blue2 <- "#B0D4FF"
red1 <- "#DD4433"
red2 <- "#FF9980"
green1 <- "#229966"


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
par(cex=1.5, mar=c(4,4,1,0)+.5);

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
## 2) Effect size and power: comparison of barplots
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
par(cex=1.5, mar=c(4,4,1,0)+.5);

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
## 2) Greatly increased sampling variation if data from larger sampling units are pooled
##    (here: two books from sections of library with pi=10% and pi=40%, vs. random sample with pi=25%)
k <- 0:50
n1 <- 50
n2 <- 50
n <- n1 + n2
p1 <- .1
p2 <- .4
p <- (n1 * p1 + n2 * p2) / n

d.random <- dbinom(k, n, p) # binomial sampling distribution for random sample of size n
tail.random <- sapply(k, function (k) binom.test(k, n, p)$p.value < .05)
## for pooled data, we have sample of 100 items from section 1 in 25% of cases, sample of 100 items from section 2 in another 25% of cases
d.p11 <- .25 * dbinom(k, n, p1)
d.p22 <- .25 * dbinom(k, n, p2)
## in the remaining 50% of cases, we pool a sample of 50 items from section 1 and a sample of 50 items from section 2
## (this is close to a random sample of 100 items from the full library, but we derive the precise distribution by explicit summation)
k.grid <- expand.grid(k1=0:n1, k2=0:n2)
k.grid <- transform(k.grid, d1=dbinom(k1, n1, p1), d2=dbinom(k2, n2, p2))
k.grid <- transform(k.grid, d=d1*d2, k=k1+k2)
d.p12 <- .50 * as.vector(rowsum(k.grid$d, factor(k.grid$k, levels=0:n)))
d.pooled <- d.p11 + d.p12[seq(along=k)] + d.p22

screen.device(width=12, height=5, bg="white")
par(cex=1.5, mar=c(4,4,0,0)+.5, xaxs="i", yaxs="i");

barplot(100 * d.random, names.arg=ifelse(k %% 5 == 0, k, NA), ylim=c(0,10),
        col=ifelse(tail.random, "black", "grey"), density=ifelse(tail.random, 10, NA), angle=45, 
        xlab="value k of observed frequency X", ylab="percentage of samples with X=k")
dev.copy2pdf(file="../keynote-slides/img/pooled_data_100_tokens.pdf", onefile=FALSE)

barplot(100 * d.pooled, names.arg=ifelse(k %% 5 == 0, k, NA), ylim=c(0,10), col="grey",
        xlab="value k of observed frequency X", ylab="percentage of samples with X=k")
dev.copy2pdf(file="../keynote-slides/img/pooled_data_2_books.pdf", onefile=FALSE)

dev.off()


## 3) Comparing binomial and Gaussian sampling distributions
##    (to illustrate that Gaussian distribution can model overdispersion by controlling spread of the distribution)
compare.binomial.gaussian <- function (p, n=100, sigma=sqrt(n*p*(1-p)), k.range=0:n, y.max=NULL) {
  mu <- n * p
  d.binom <- 100 * dbinom(k.range, n, p)
  d.gauss <- 100 * dnorm(k.range, mu, sigma)
  if (missing(y.max)) y.max <- max(d.binom, d.gauss)
  plot(range(k.range), c(0, y.max), type="n", xaxs="i", yaxs="i", xlab="", ylab="")
  lines(k.range, d.binom, lwd=3, col="blue")
  lines(k.range, d.gauss, lwd=3, col="red")
  legend.txt <- eval(bquote(expression(B(.(n),.(p)), N(.(mu),.(round(sigma,2))))))
  legend("topright", inset=.02, legend=legend.txt, lwd=4, col=c("blue", "red"))
}

screen.device(width=12, height=8, bg="white")
par(mar=c(2,2,2,2)+.5, mfrow=c(2,3), cex=.8)

compare.binomial.gaussian(.2, sigma=2.5, k.range=0:40, y.max=20)
compare.binomial.gaussian(.2, k.range=0:40, y.max=20)
compare.binomial.gaussian(.2, sigma=8, k.range=0:40, y.max=20)

compare.binomial.gaussian(.1, sigma=2, k.range=0:40, y.max=20)
compare.binomial.gaussian(.1, k.range=0:40, y.max=20)
compare.binomial.gaussian(.1, sigma=6, k.range=0:40, y.max=20)

dev.copy2pdf(file="../keynote-slides/img/binomial_vs_gaussian.pdf", onefile=FALSE)
dev.off()
