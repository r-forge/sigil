## 
## Solution to Exercise #3 (by-text frequency comparison)
## 

library(SIGIL)
Pass <- PassiveBrownFam
nrow(Pass)
head(Pass) # take a first look at the data

## -- Question 1 --

## Tables of text counts by corpus and genre 
table(Pass$corpus)
table(Pass$genre)

## Is the genre distribution the same for each corpus?
xtabs(~ genre + corpus, data=Pass) # cross-tabulation
tb <- table(Pass$genre, Pass$corpus) # alternative (equivalent)
print(tb)

## Optional: visualization (with some nice display options)
par(mar=c(7,4,2,1)+.1) # better margins for plots below
barplot(t(tb), beside=TRUE, col=1:5, border=NA, las=3)
legend("topright", legend=colnames(tb), fill=1:5)

## -- Question 2 --

## Obtain separate data sets for each of the four target corpora
Brown <- subset(Pass, corpus == "Brown")
nrow(Brown) # plausibility check is a good idea
table(Brown$corpus)
LOB <- subset(Pass, corpus == "LOB")
Frown <- subset(Pass, corpus == "Frown")
FLOB <- subset(Pass, corpus == "FLOB")

## -- Question 3 --

## Compute mean and standard deviation of passive percentages for each corpus
mean(Brown$p.pass) # mean for Brown corpus
sd(Brown$p.pass)   # s.d. for Brown corpus
## same for the other 3 corpora

cat(sprintf("Brown: μ = %.2f%%, σ = %.2f%%\n", mean(Brown$p.pass), sd(Brown$p.pass))) # nicer formatting
cat(sprintf("LOB:   μ = %.2f%%, σ = %.2f%%\n", mean(LOB$p.pass),   sd(LOB$p.pass)))
cat(sprintf("Frown: μ = %.2f%%, σ = %.2f%%\n", mean(Frown$p.pass), sd(Frown$p.pass)))
cat(sprintf("FLOB:  μ = %.2f%%, σ = %.2f%%\n", mean(FLOB$p.pass),  sd(FLOB$p.pass)))

## Which other summary statistics might be useful?
summary(Brown$p.pass) # median, min, max, quartiles
library(e1071) # skewness and kurtosis not in standard R, but in package e1071
skewness(Brown$p.pass)
kurtosis(Brown$p.pass)

## -- Question 4 --

## Plot a histogram for each corpus
hist(Brown$p.pass, xlim=c(0,100), ylim=c(0,150), xlab="passives (%)", ylab="number of texts", main="Brown") # specify xlim and ylim explicitly so plots can be compared
## same for other three corpora

## Add contour line showing the estimated density function
hist(Brown$p.pass, freq=FALSE, xlim=c(0,100), ylim=c(0,0.07), xlab="passives (%)", ylab="density", main="Brown") # needs density scaling (note differen ylim)
lines(density(Brown$p.pass), col="blue", lwd=3)
rug(Brown$p.pass, col="blue")

## Optional: can you combine all four histograms / density curves?
d.Brown <- density(Brown$p.pass) # compute density curves
d.LOB <- density(LOB$p.pass)
d.Frown <- density(Frown$p.pass)
d.FLOB <- density(FLOB$p.pass)
plot(d.Brown, col=1, lwd=3, xlim=c(0,100), ylim=c(0,0.07), xlab="passives (%)", ylab="density", main="")
lines(d.LOB, col=2, lwd=3)
lines(d.Frown, col=3, lwd=3)
lines(d.FLOB, col=4, lwd=3)
legend("topright", inset=.02, legend=c("Brown", "LOB", "Frown", "FLOB"), col=1:4, lwd=3)
## it's difficult to combine histograms in a meaningful way, but you could show the four plots
## in a single window or PDF file (see below)

## -- Question 5 --

## Add suitable Gaussian density curves
hist(Brown$p.pass, freq=FALSE, xlim=c(0,100), ylim=c(0,0.07), xlab="passives (%)", ylab="density", main="Brown") # needs density scaling (note differen ylim)
lines(density(Brown$p.pass), col="black", lwd=3) # same as above
xC <- seq(0, 100, .5) # equidistant x-coordinates for Gaussian density
yC <- dnorm(xC, mean=mean(Brown$p.pass), sd=sd(Brown$p.pass)) # Gaussian density function fitted to observed data, i.e. use observed mean and s.d. as parameters
lines(xC, yC, col="red", lwd=3)
## same for other three corpora

## Use quantile-quantile plots to assess normality (combined in single window)
par(mfrow=c(2,2)) # 4 plots in rectangular 2x2 layout, filled by row
qqnorm(Brown$p.pass, main="Brown") # top left
qqline(Brown$p.pass, col="blue")
qqnorm(LOB$p.pass, main="LOB")     # top right
qqline(LOB$p.pass, col="blue")
qqnorm(Frown$p.pass, main="Frown") # bottom left
qqline(Frown$p.pass, col="blue")
qqnorm(FLOB$p.pass, main="FLOB")   # bottom right
qqline(FLOB$p.pass, col="blue")
par(mfrow=c(1,1)) # reset to normal plots

## conclusion: all four corpora show a substantial deviation from the Gaussian distribution,
## so some care needs to be taken when applying t-tests and other methods that make normality
## assumptions; it doesn't necessarily mean that t-tests cannot be used, though

## -- Question 6 --

## Meaningful comparisons are possible for
##   Brown vs. LOB     .. AmE/BrE in the 1960s
##   Frown vs. FLOB    .. AmE/BrE in the 1990s
##   Brown vs. Frown   .. AmE in the 1960s vs. 1990s
##   LOB vs. FLOB      .. BrE in the 1960s vs. 1990s

## Which version of the t-test?
##   t-test for two independent samples

## Specify the precise null hypothesis
##   H0: mu_1 = mu_2
##   where mu_1 is the average passive frequency in texts from variety 1 (and analogously for mu_2)

## Report and interpret full results of the significance tests
t.test(Brown$p.pass, LOB$p.pass) ## p = .223 n.s.
## interpretation: in the 1960s, there is no significant evidence for a difference between American and British English (with respect to frequency of passive use)

t.test(Brown$p.pass, Frown$p.pass) ## p = .000006 ***
## interpretation: there is a significant difference between AmE in the 1960s and 1990s
## the effect size is at least 1.45 percent points, i.e. percentage of passive VPs is at least 1.45 higher in the 1960s (= Brown corpus) than in the 1990s (= Frown corpus)

t.test(LOB$p.pass, FLOB$p.pass) ## p = .023 *
## interpretation: there is a significant difference between BrE in the 1960s and 1990s
## evidence isn't very strong, though, and we cannot be sure that the effect size is larger than 0.18 percent points

t.test(Frown$p.pass, FLOB$p.pass) ## p = .000061 ***
## interpretation: there is a significant difference between AmE and BrE in the 1990s
## the effect size is at least 1.03 percent points (AmE lower than BrE)

## Optional: non-parametric Mann-Whitney tests
wilcox.test(Brown$p.pass, LOB$p.pass) # p = .0495 * (marginally significant)
## interpretation: Brown vs. LOB is not significant according to t-test, but marginally significant according to Mann-Whitney (p just below 0.05)
## this discrepancy may be due to the non-Gaussian distribution we observed in Q5, since the t-test assumes Gaussian populations

wilcox.test(LOB$p.pass, FLOB$p.pass) # p = .091 n.s.
## in contrast to the t-test, Mann-Whitney doesn't find a significant difference between BrE in the 1960s and 1990s

## same for the other two comparisons (both are highly significant)

## -- Question 7 --

## Is there a significant difference between AmE and BrE?
##  - no significant difference in 1960s (Brown vs. LOB), but significant in 1990s (Frown vs. FLOB)
##  - effect size in 1990s seems large enough (> 1% point) for a linguistic interpretation

## Do your data corroborate the claim of a decline in passive use?
##  - significance of difference between LOB and FLOB is doubtful -> no clear evidence for decline in BrE
##  - but massive decline in AmE with effect size of at least 1.45% points (Brown vs. Frown)

## NB: there is an entirely different interpretation of the test results
##  - we only have obtained clearly significant results for Frown vs. Brown and Frown vs. FLOB
##  - so perhaps Frown is "the odd one out", and there are no substantial differences between the other three corpora
##  - this implies that something has changed in AmE in the 1990s, but doesn't support the conclusion of a general difference between AmE and BrE, or a decline of passive use in English
##  - it is also possible that the lower passive frequency is an artefact of the corpus design (e.g. slightly different genres)

## Linguistic interpretation
##  - in contrast to Leech et al. (2009), the data do not support a decline in passive use in the English language led by American English
##  - however, this does not prove that there is no such decline (and the observed data suggest there may be)
##  - all we can state is that the available corpora do not provide sufficient evidence for the decline
##  - and we can put an upper limit on the possible effect size (viz. the "larger" end of the confidence intervals)

## -- Question 8 --

## Display side-by-side boxplots of passive frequencies in different genres
boxplot(p.pass ~ genre, data=Brown, las=3)
abline(h=mean(Brown$p.pass), col="red")
## same for other three corpora

## Which significance test can be used to assess differences between genres?
##   ANOVA (analysis of variance)
## with null hypothesis
##   H0: mu_1 = mu_2 = ... = mu_15
res <- aov(p.pass ~ genre, data=Brown)
summary(res) # highly significant difference

## Why are pairwise comparisons between all genres using t-tests problematic?
n <- nlevels(Brown$genre) # number of genres
choose(15, 2) # there are 105 pairwise tests
## for each individual test, risk of type I error is 5% (at nominal significance level p < .05)
## but among 105 tests, we have to expect 5 false positives (assuming H0 is always false)
## family-wise error rate = risk of at least one type I error in entire group of tests
dbinom(5, 105, p=0.05) # probability of exactly 5 false positives (assuming independent tests)
pbinom(0, 105, p=0.05, lower=FALSE) # FWER = probability of at least one type I error

## -- Question 9 --

## see Unit #8: Non-randomness
