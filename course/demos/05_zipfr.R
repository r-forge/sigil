## SIGIL: Statistical Inference - a Gentle Introduction for Linguists
## Unit 5: Word Frequency Distributions
## -- code examples --

## load library & check for documentation (really good in this case)
library(zipfR)
?zipfR
data(package="zipfR")
help.start() # you'll find the zipfR tutorial here

## importing data from package or external files
data(ItaRi.spc) # frequency distribution of Italian prefix ri-
data(ItaRi.emp.vgc)

## here is how you load data from disk files (see help pages for file formats)
# my.spc <- read.spc("my.spc.txt") 
# my.vgc <- read.vgc("my.vgc.txt")
# my.tfl <- read.tfl("my.tfl.txt")
# my.spc <- tfl2spc(my.tfl)

## printing and plotting frequency spectra
summary(ItaRi.spc)
ItaRi.spc

N(ItaRi.spc)
V(ItaRi.spc)
Vm(ItaRi.spc,1)
Vm(ItaRi.spc,1:5)

Vm(ItaRi.spc,1) / N(ItaRi.spc) # Baayen's productivity index P

plot(ItaRi.spc)
plot(ItaRi.spc, log="x")

## vocabulary growth curves (VGC)
summary(ItaRi.emp.vgc)
ItaRi.emp.vgc
N(ItaRi.emp.vgc)

plot(ItaRi.emp.vgc, add.m=1)

## binomial interpolation of VGCs: test sampling assumptions
ItaRi.bin.vgc <- vgc.interp(ItaRi.spc, N(ItaRi.emp.vgc), m.max=1)
summary(ItaRi.bin.vgc)

plot(ItaRi.emp.vgc, ItaRi.bin.vgc, legend=c("observed", "interpolated"))

## comparison with Italian prefix ultra- (mini-exercise!)
data(ItaUltra.spc)
data(ItaUltra.emp.vgc)

plot(ItaRi.spc, ItaUltra.spc, legend=c("ri-", "ultra-")) # comparison of frequency spectra: widely different sample sizes
plot(ItaRi.emp.vgc, ItaUltra.emp.vgc, add.m=1, legend=c("ri-", "ultra-")) # comparison of vocabulary growth

ItaUltra.bin.vgc <- vgc.interp(ItaUltra.spc, N(ItaUltra.emp.vgc), m.max=1)
plot(ItaRi.bin.vgc, ItaUltra.bin.vgc, add.m=1, legend=c("ri-", "ultra-")) # interpolated vocabulary growth

## estimating LNRE models
ItaUltra.fzm <- lnre("fzm", ItaUltra.spc) # also try "zm" and "gigp"
summary(ItaUltra.fzm)

## expected frequency spectrum
ItaUltra.fzm.spc <- lnre.spc(ItaUltra.fzm, N(ItaUltra.fzm))
ItaUltra.fzm.spc

plot(ItaUltra.spc, ItaUltra.fzm.spc, legend=c("observed", "fzm"))
plot(ItaUltra.spc, ItaUltra.fzm.spc, legend=c("observed", "fzm"), m.max=10) # first 10 only

## extrapolation of ultra- VGC to sample size of ri- data
ItaUltra.ext.vgc <- lnre.vgc(ItaUltra.fzm, N(ItaRi.emp.vgc))

## compare the two growth curves
plot(ItaRi.bin.vgc, ItaUltra.ext.vgc, N0=N(ItaUltra.fzm), legend=c("ri-", "ultra-"))

plot(ItaRi.bin.vgc, ItaUltra.ext.vgc, N0=N(ItaUltra.fzm), legend=c("ri-", "ultra-"), xlim=c(0, 100e3)) # zoom in

## how reliable is the model fit for the small ItaUltra data set?
extract.stats <- function (m) data.frame(alpha=m$param$alpha, A=m$param$A, B=m$param$B, S=m$S, X2=m$gof$X2)
runs <- lnre.bootstrap(ItaUltra.fzm, N(ItaUltra.fzm), lnre, extract.stats, type="fzm")
runs <- do.call(rbind, runs)

hist(runs$alpha, freq=FALSE, xlim=c(0, 1))
lines(density(runs$alpha), lwd=2, col="red")
abline(v=ItaUltra.fzm$param$alpha, lwd=2, col="blue")

hist(runs$X2, freq=FALSE)
lines(density(runs$X2), lwd=2, col="red")

hist(runs$S) # what is wrong here?
