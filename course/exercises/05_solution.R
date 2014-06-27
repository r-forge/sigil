##
## Our solution to Exercise #5(type richness and vocabulary growth)
##

## Problem 1: German NP and PP constructions

library(zipfR)  # start by loading the zipfR library

## - Use the German NP and PP data sets included in zipfR
N(TigerNP.spc)  # check which data set is larger
N(TigerPP.spc)

N1 <- N(TigerNP.spc) # NP data set is larger, so we will plot curves up to N1 = N(TigerNP.spc)
N.values <- seq(0, N1, length.out=100) # sequence with specified number of steps between endpoints

## - For the data set with more tokens, compute a binomially interpolated vocabulary growth curve (VGC).
NP.vgc <- vgc.interp(TigerNP.spc, N.values) # binomial interpolation 

## - For the smaller one, estimate an LNRE model ...
PP.gigp <- lnre("gigp", spc=TigerPP.spc)    # estimate GIGP model (for a change) on PP data
summary(PP.gigp)  # note the excellent goodness-of-fit (fZM is even better)

## - ... and use it to compute an expected VGC up to the size of the larger data set.
PP.vgc <- lnre.vgc(PP.gigp, N.values)

## - Plot the interpolated and extrapolated VGCs, and determine which of the two constructions appears to be more productive
## vertical line (N0 parameter) indicates size of smaller sample (PP data), i.e. where extrapolation with LNRE model starts
plot(NP.vgc, PP.vgc, legend=c("NP","PP"), N0=N(TigerPP.spc))

## Conclusion: both curves have similar shape, but the NP curve produces consistenly more types, so the NP rules seem to be more productive than PP rules.

## An interesting extension is to compare the interpolated / model-based growth curves with actual vocabulary growth (for NP and PP rules) in the Tiger corpus
plot(NP.vgc, PP.vgc, TigerNP.emp.vgc, TigerPP.emp.vgc, legend=c("NP (interpolated)", "PP (GIGP)", "NP (actual)", "PP (actual)"))

## Actual PP growth seems to be somewhat steeper than predicted by the model, possibly because of non-randomness in the data.  As a consequence, the true divergence of vocabulary growth curves may not be as large as suggested by our analysis (in fact, PP and NP curves almost seem to run parallel in the range from N=60,000 to N=90,000).


## Problem 2: Data lost with cut-off points

## - First, load the file bigrams.100k.tfl (containing bigrams extracted from the first 100,000 tokens of the Brown corpus)
bigrams.tfl <- read.tfl("bigrams.100k.tfl") # don't forget to set the working directory first
bigrams.spc <- tfl2spc(bigrams.tfl)         # compute frequency spectrum

## - What proportion of bigram types occur only once? What proportion of types occur 4 times or less?
print(bigrams.spc)  # you could also calculate the proportions "by hand" from this printout
plot(bigrams.spc)   # visualise the frequency spectrum

Vm(bigrams.spc, 1) / V(bigrams.spc)  # proportion of bigram types occurring once
sum(Vm(bigrams.spc, 1:4)) / V(bigrams.spc) # note how Vm() takes a vector to return V1, V2, V3, V4

## - Now, suppose that we want to use our data to estimate the proportion of bigram types that would be lost by using the same two frequency cut-offs if we had a 1 million word corpus.  How do you go about it?
bigrams.zm <- lnre("zm", spc=bigrams.spc)  # estimate LNRE model
## Baroni & Evert (2007) suggest that ZM model might give the best extrapolation accuracy if there is non-randomness, even though GIGP has much better goodness-of-fit

## now we can use the model to predict relevant information for a sample of 1M tokens
N <- 10 * N(bigrams.spc)  # number of bigram tokens is smaller than number of words, so we can't just set N = 1e6 in order to make predictions for corpus of 1M tokens; instead, we estimate the number of bigram tokens in the 1M word corpus as 10 times the number of bigram tokens in the 100k corpus
EVm(bigrams.zm, 1, N) / EV(bigrams.zm, N)
sum(EVm(bigrams.zm, 1:4, N)) / EV(bigrams.zm, N)
# note that the type proportions are somewhat smaller than for 100k tokens

# Final remark: comparison with the actual frequency spectrum for bigrams from the full 1M word Brown corpus shows that the predictions made by ZM are indeed much better than those of GIGP, especially for the proportion of hapax legomena (true value: 74.7%, ZM predicts 80.5%, GIGP predicts 65.4%).


## Problem 3: Reliability of the fitted model

TigerNP.fzm <- lnre("fzm", TigerNP.spc)
extract.stats <- function (m) data.frame(alpha=m$param$alpha, S=m$S, X2=m$gof$X2)

## helper function: histogram with rug plot and density curve
show.hist <- function (x, true.val=NA, ...) {
  hist(x, freq=FALSE, ...)
  rug(x)
  if (!missing(true.val)) abline(v=true.val, lwd=2, lty="dashed")
  lines(density(x), lwd=3, col="blue")
}

## reliability for full Tiger treebank (N = 109k) 
runs1 <- lnre.bootstrap(TigerNP.fzm, N(TigerNP.fzm), lnre, extract.stats, type="fzm") 
runs1 <- do.call(rbind, runs1)

summary(runs1) # get an overview

show.hist(runs1$alpha, xlim=c(0, 1), xlab=expression(alpha))
show.hist(runs1$S, xlab="S") # usually, there are some clear outliers
show.hist(runs1$X2, xlab=expression(X^2))
summary(TigerNP.fzm) # X2 of the original model is in the typical range

## reliability for 10% of the Tiger treebank (N = 11k)
runs2 <- lnre.bootstrap(TigerNP.fzm, 11e3, lnre, extract.stats, type="fzm") 
runs2 <- do.call(rbind, runs2)

summary(runs2) # already hints that there may be problems ...

show.hist(runs2$alpha, xlim=c(0, 1), xlab=expression(alpha)) # somewhat more variation
show.hist(runs2$S, xlab="S") # now we get very bad outliers
S2 <- runs2$S[runs2$S < 100e3] # discard outliers (just setting xlim= doesn't help in this case)
show.hist(S2, xlab="S", xlim=c(0, 100e3)) # estimation of S from small samples seems inherently unreliable 
show.hist(runs2$X2, xlab=expression(X^2)) # even though there are no clear parameter estimation failures 

## reliability of extrapolation from 10% sample to original data size
expected.vals <- function (m) data.frame(V=EV(m, N(TigerNP.spc)), V1=EVm(m, 1, N(TigerNP.spc)))
runs3 <- lnre.bootstrap(TigerNP.fzm, 11e3, lnre, expected.vals, type="fzm") 
runs3 <- do.call(rbind, runs3)

summary(runs3)

show.hist(runs3$V, true.val=V(TigerNP.spc), xlim=c(0, 5000)) # some variability, but reasonable predictions
show.hist(runs3$V1, true.val=Vm(TigerNP.spc, 1), xlim=c(0, 5000)) # relative variability much larger for V1

