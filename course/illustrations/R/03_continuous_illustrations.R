##
##  Illustrations and examples for Unit 3: Descriptive and Inferential Statistics for Continuous Data
##

library(corpora)

# generate simulated population data
FakeCensus <- simulated.census()
WackypediaStats <- simulated.wikipedia()

source("functions.R")
dev.save <- function (basename) {
	dev.copy2pdf(file=paste("../img/",basename,".pdf",sep=""), bg="white", onefile=FALSE)
}

##
## Introduction: examples and characteristic measures
##

# samples of census and wikipedia stats tables (rounded) 
tmp <- head(FakeCensus, 20)
tmp <- transform(tmp, height=round(height, 2), weight=round(weight, 2))
print(tmp)

tmp2 <- head(WackypediaStats, 20)
tmp2 <- transform(tmp2, ttr=round(ttr, 3), avglen=round(avglen, 3))
print(tmp2)

# we can also calculate some derived quantities
FakeCensus <- transform(FakeCensus, bmi=weight/(height/100)^2)

# mean, variance, etc. for census data
m <- nrow(FakeCensus)
print(m)

mean(FakeCensus$height)
mu.w <- mean(FakeCensus$weight)
print(mu.w)
mean(FakeCensus$shoe.size)

# sample and population variance are virtually indistinguishable
var(FakeCensus$weight) # sample (denominator m-1)
sigma2.w <- sum((FakeCensus$weight - mu.w)^2) / m # population
print(sigma2.w)

# mean, variance, s.d. for different attributes
mvs <- function (x) c(mean(x), var(x), sd(x))
mvs(FakeCensus$height)
mvs(FakeCensus$weight)
mvs(FakeCensus$shoe.size)

##
## Some abbreviations
##

FC <- FakeCensus
WS <- WackypediaStats

m.FC <- nrow(FC)
m.WS <- nrow(WS)

##
## Histogram, density, types of distributions
##

square.device(main=FALSE)

# discrete distribution can be tabulated
plot(table(FC$shoe.size) / m.FC, lwd=3, xlab="Shoe size", ylab="Proportion of population")
dev.save("ingary_shoesize_table")

# continuous distribution must be cut into bins: histogram
hmin <- 110; hmax <- 230
hist(FC$height, col=fillblue, labels=TRUE, main="", xlab="body height", xlim=c(120,220), ylim=c(0,65000), breaks=seq(hmin, hmax, 5))
dev.save("ingary_hist_freq_1")
hist(FC$height, col=fillblue, main="", xlab="body height", xlim=c(120,220), ylim=c(0,65000), breaks=seq(hmin, hmax, 2.5))
dev.save("ingary_hist_freq_2")

# density scale: area of bar = relative frequency / proportion in population
hist(FC$height, col=fillblue, freq=FALSE, main="", xlab="body height", xlim=c(120,220), ylim=c(0,.025), breaks=seq(hmin, hmax, 5))
dev.save("ingary_hist_1")
hist(FC$height, col=fillblue, freq=FALSE, main="", xlab="body height", xlim=c(120,220), ylim=c(0,.025), breaks=seq(hmin, hmax, 2.5))
dev.save("ingary_hist_2")
hist(FC$height, col=fillblue, freq=FALSE, main="", xlab="body height", xlim=c(120,220), ylim=c(0,.025), breaks=seq(hmin, hmax, 1))
dev.save("ingary_hist_3")
hist(FC$height, col=fillblue, border=NA, freq=FALSE, main="", xlab="body height", xlim=c(120,220), ylim=c(0,.025), breaks=seq(hmin, hmax, .4))
dev.save("ingary_hist_4")
lines(density(FC$height, from=hmin, to=hmax), col=myred, lwd=6)
dev.save("ingary_hist_4_curve")

# area under density function = range probability
par.save <- par(mar=c(2,0,0,0)+.1, cex=2.5)
x <- seq(0, 11, .1)
y <- dchisq(x, 3)
plot(x, y, type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,10), ylim=c(0,.25))
abline(h=0)
axis(1, c(3, 6), labels=expression(a, b))
idx <- x >= 3 & x <= 6
polygon(c(x[idx], rev(x[idx])), c(y[idx], rep(0, sum(idx))), col=fillblue)
segments(c(3,6), c(0,0), c(3,6), dchisq(c(3,6), 3), col="black", lwd=3)
lines(x, y, lwd=6, col=myred)
dev.save("area_under_density")
par(par.save)

# different "typical" kinds of distributions
wide.device()
mean.sd.plot(WackypediaStats$avglen, xlim=c(3.2,5.5)) # almost Gaussian
dev.save("disttype_symmetric")
mean.sd.plot(FakeCensus$height, xlim=c(125,215)) # similar to Gaussian, but more bulgy
dev.save("disttype_bulgy")
x <- exp(rnorm(1e6,0,1)) # truncated lognormal
mean.sd.plot(x[x<=10], xlim=c(-2,6), median=TRUE)
dev.save("disttype_skewed")
mean.sd.plot(FakeCensus$weight, xlim=c(25,140), median=TRUE)
dev.save("disttype_complicated")
mean.sd.plot(FakeCensus$shoe.size, median=TRUE)
dev.save("disttype_bimodal")
dev.off()

##
## Parameters of the normal distribution
##

k <- 0:30
n <- 100
p <- .15
mu <- n * p
sigma <- sqrt(n * p * (1-p))
x.norm <- seq(min(k)-1, max(k)+1, .05)   # x values of bell curve
y.norm <- 100 * dnorm(x.norm, mean=mu, sd=sigma) # y values of bell curve, scaled to percentages

dev.new(width=12, height=7, bg="white")
par(xaxt="n", yaxt="n", xaxs="i", yaxs="i", cex=1.5, mar=c(2,2,.5,.5)+.1)

plot(x.norm, y.norm, type="n", xlim=c(0,30), ylim=range(y.norm)+c(-.05,.2), xlab="", ylab="")
title(xlab=expression(t), ylab=expression(g(t)), line=1)
abline(h=0)                             # show x-axis
abline(v=mu, lwd=2, col=myblue)         # mean = mode of the distribution
abline(v=mu+sigma, lwd=1, col=myred)    # plus/minus 1 standard deviation
abline(v=mu-sigma, lwd=1, col=myred)
abline(v=mu+2*sigma, lwd=1, col=myred)  # plus/minus 2 standard deviations
abline(v=mu-2*sigma, lwd=1, col=myred)
y0 <- 100 * dnorm(mu+sigma, mean=mu, sd=sigma)

idx.tail <- x.norm >= mu + 2*sigma # shade two-sided 5% tails of distribution (beyond 2*sigma)
x.tail <- x.norm[idx.tail]
y.tail <- y.norm[idx.tail]
x.poly <- c(x.tail, rev(x.tail))
y.poly <- c(y.tail, rep(0, length=length(y.tail)))
polygon(x.poly, y.poly, col=fillblue, border="black")
idx.tail <- x.norm <= mu - 2*sigma # left tail
x.tail <- x.norm[idx.tail]
y.tail <- y.norm[idx.tail]
x.poly <- c(x.tail, rev(x.tail))        # boundary of polygon shape for tail
y.poly <- c(y.tail, rep(0, length=length(y.tail)))
polygon(x.poly, y.poly, col=fillblue, border="black")

lines(x.norm, y.norm, type="l", lwd=6, col="black") # draw density function here to cover background material
text(mu, .25, pos=4, labels=expression(mu), cex=1.6, col=myblue) # label for mu
# arrows for plus/minus 1 and 2 standard deviations
arrows(mu, y0, c(mu+sigma, mu-sigma), y0, lwd=3, col=myred, angle=20, length=.2) # uses magic replication of single ("scalar") values
text(c(mu+sigma/2, mu-sigma/2), y0, pos=3, labels=expression(sigma), cex=1.4, col=myred)
y1 <- 100 * dnorm(mu+2*sigma, mean=mu, sd=sigma)
arrows(mu, y1, c(mu+2*sigma, mu-2*sigma), y1, lwd=3, col=myred, angle=20, length=.2)
text(c(mu+sigma/2, mu-sigma/2), y1, pos=3, labels=expression(2*sigma), cex=1.4, col=myred)

dev.save("gaussian_parameters")
dev.off()

##
## Assessing normality
##

set.seed(10)
x <- sample(FakeCensus$weight, 200) # 200 items from "complicated" distribution
y <- sample(FakeCensus$height, 200) # 200 items from bulgy distribution (but reasonably close to normality)

par.save <- par(cex=1.4)

normality.plot(y, ylim=c(0,.03), xlim=c(100,220), breaks=seq(100,220,10), show.mu=TRUE, show.normal=FALSE)
dev.save("normality_bulgy_1")
normality.plot(y, ylim=c(0,.03), xlim=c(100,220), breaks=seq(100,220,10), show.mu=TRUE, show.normal=TRUE)
dev.save("normality_bulgy_2")
normality.plot(x, ylim=c(0,.03), xlim=c(0,150), breaks=seq(0,150,10), show.mu=TRUE, show.normal=TRUE)
dev.save("normality_complicated")

qqnorm(y, pch=20, main="")
qqline(y, col="red", lwd=3)
dev.save("qqnorm_bulgy")
qqnorm(x, pch=20, main="")
qqline(x, col="red", lwd=3)
dev.save("qqnorm_complicated")

par(par.save)


##
## -- draft material --
##


mean.sd.plot(FakeCensus$height, median=TRUE) # similar to Gaussian, but more bulgy: well-described
mean.sd.plot(FakeCensus$weight, median=TRUE) # skewed
mean.sd.plot(FakeCensus$bmi, median=TRUE)    # somewhat skewed, but still well-described
mean.sd.plot(FakeCensus$shoe.size, median=TRUE)  # bimodal

mean.sd.plot(log10(WackypediaStats$tokens), median=TRUE) # these are all very well-behaved
mean.sd.plot(log10(WackypediaStats$types), median=TRUE)
mean.sd.plot(WackypediaStats$ttr, median=TRUE)
mean.sd.plot(WackypediaStats$avglen, median=TRUE)
