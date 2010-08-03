##
##  Support functions for illustrations and example
##

# square-format quartz device
square.device <- function (width=8, cex=1.2, main=FALSE, ...) {
	quartz(width=width, height=width, bg="white")
	margins <- if (main) c(4,4,3,1) else c(4,4,1,1)
	par(cex=cex, mar=margins+.1, ...)
}

# wide format quartz device (for fullscreen image on slide)
wide.device <- function (width=10, aspect=4/3, cex=1.2, main=FALSE, ...) {
	quartz(width=width, height=width/aspect, bg="white")
	margins <- if (main) c(4,4,3,1) else c(4,4,1,1)
	par(cex=cex, mar=margins+.1, ...)
}

# some nice colours
myblue <- rgb(.1,.2,.8)
myred <- rgb(1,.1,.25)
mygreen <- rgb(.2,.6,.1)
fillblue <- rgb(.2,.4,1)
fillgrey <- rgb(.6,.6,.6)

# histogram vs. normal approximation (NB: pass breaks=, xlim=, ylim=, ... directly to hist() function)
normality.plot <- function (x, density.estimate=FALSE, xlab="", ...) {
	hist(x, col=fillgrey, freq=FALSE, xlab=xlab, ...)
	x.min <- par("usr")[1]
	x.max <- par("usr")[2]
	if (density.estimate) lines(density(x, from=x.min, to=x.max), lwd=6, col=myblue)
	mu <- mean(x)
	sigma <- sd(x)
	y <- seq(x.min, x.max, length.out=200)
	lines(y, dnorm(y, mu, sigma), lwd=6, col=myred)
}

# visualise mean +/- 1 or 2 s.d. for various shapes of distribution
mean.sd.plot <- function (x, xlim=range(x), xlab="", main="", median=FALSE, ...) {
	plot(density(x), type="l", lwd=6, col="black", xlim=xlim, xlab=xlab, main=main, ...)
	abline(h=0)
	mu <- mean(x)
	sigma <- sd(x)
	y <- sum(c(0,0,.95,.05) * par("usr")) # clever, eh?
	abline(v = mu, lwd=4, col=myblue)
	text(mu, y, expression(mu), adj=c(0,1.5), srt=90, col=myblue, cex=1.3)
	abline(v = mu + sigma, lwd=3, col=myred)
	text(mu + sigma, y, expression(mu + sigma), adj=c(0,1.5), srt=90, col=myred, cex=1.3)
	abline(v = mu - sigma, lwd=3, col=myred)
	text(mu - sigma, y, expression(mu - sigma), adj=c(0,-.5), srt=90, col=myred, cex=1.3)
	abline(v = mu + 2*sigma, lwd=2, col=myred)
	text(mu + 2*sigma, y, expression(mu + 2*sigma), adj=c(0,1.5), srt=90, col=myred, cex=1.3)
	abline(v = mu - 2*sigma, lwd=2, col=myred)
	text(mu - 2*sigma, y, expression(mu - 2*sigma), adj=c(0,-.5), srt=90, col=myred, cex=1.3)
	if (median) {
		med <- median(x)
		y <- sum(c(0,0,.05,.95) * par("usr"))
		abline(v = med, lwd=4, lty="dashed", col=mygreen)
		text(med, y, "median", adj=c(1,1.5), srt=90, col=mygreen, cex=1.2)
	}
}













