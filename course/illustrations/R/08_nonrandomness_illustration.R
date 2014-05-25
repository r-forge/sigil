##
##  Non-randomness for passives in the Brown and LOB corpora
##

## *** TODO: revise this codes and re-create images ***

##
## 1) Greatly increased sampling variation if data from larger sampling units are pooled
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


## 2) Comparing binomial and Gaussian sampling distributions
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
