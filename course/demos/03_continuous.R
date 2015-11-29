## SIGIL: Statistical Inference - a Gentle Introduction for Linguists
## Unit 3: Descriptive and inferential statistics for continuous data
## -- code examples --

library(SIGIL)
library(corpora) # for sample.df

##
## Part A: Descriptive Statistics
##

## generate simulated population data
Census <- simulated.census()
Wiki <- simulated.wikipedia()

## histograms visualize distribution of variable
hist(Census$height, col="lightblue", breaks=50, freq=FALSE)
lines(density(Census$height), lwd=4, col="red")

hist(Census$shoe.size, col="lightblue", breaks=seq(29.75,55.25,.5), freq=FALSE)
lines(density(Census$shoe.size), lwd=4, col="red")

## histogram for discrete data: use table() rather than hist()
plot(2 * table(Census$shoe.size) / nrow(Census))
lines(density(Census$shoe.size), lwd=4, col="red")

## compare histogram with Gaussian distribution
hist(Census$weight, col="lightblue", breaks=50, freq=FALSE, xlim=c(30,150), ylim=c(0,.03))
lines(density(Census$weight), lwd=4, col="red")
mu <- mean(Census$weight)
sigma <- sd(Census$weight)
data.frame(mu=mu, sigma=sigma)
abline(v=mu, lwd=3, lty="dashed")
abline(v=c(mu-2*sigma, mu+2*sigma), lwd=2)

## normal (Gaussian) approximation: best parameters are empirical mu and sigma
t <- seq(35, 145, length.out=200)
g.t <- dnorm(t, mean=mu, sd=sigma)
lines(t, g.t, lwd=4, col="blue")

## histograms for Wackypedia 
hist(Wiki$ttr, col="lightblue", breaks=50, freq=FALSE)
lines(density(Wiki$ttr), lwd=4, col="red")

hist(Wiki$avglen, col="lightblue", breaks=50, freq=FALSE)
lines(density(Wiki$avglen), lwd=4, col="red")

## a skewed distribution
x <- Wiki$types # cheat: remove outliers for reasonable plot
sum(x >= 700)
x <- x[x < 700]
hist(x, col="lightblue", breaks=50, freq=FALSE)
lines(density(x), lwd=4, col="red")

## log-transformation of skewed distribution often has a Gaussian shape
x <- Wiki$types # no need to remove outliers here
hist(log10(x), col="lightblue", breaks=50, freq=FALSE) # log-transformed counts
lines(density(log10(x)), lwd=4, col="red")

hist(log10(x), col="lightblue", breaks=200, freq=FALSE) # type counts are really discrete

## Now take a small random sample (i.e. a survey)
Survey <- sample.df(Census, 100, sort=TRUE)

## Assessing normality: histogram & density function
x <- Survey$weight

hist(x, freq=FALSE) # example code from slides
lines(density(x))
xG <- seq(min(x),max(x),len=100)
yG <- dnorm(xG,mean(x),sd(x))
lines(xG,yG,col="red")

qqnorm(x)           # example code from slides
qqline(x,col="red")

## histograms and Q-Q plots for some other variables
hist(Survey$height, breaks=20, freq=FALSE) # is normality assumption plausible?
lines(density(Survey$height), lwd=2, col="red")

qqnorm(Survey$height) # q-q plot better for small samples
qqline(Survey$height, lwd=2, col="blue")

qqnorm(Survey$weight) # skewed distribution in q-q plot
qqline(Survey$weight, lwd=2, col="blue")

qqnorm(Survey$shoe.size) # bimodal discrete distribution in q-q plot
qqline(Survey$shoe.size, lwd=2, col="blue")


##
## Part B: Inferential Statistics
##

## statistical distributions in R: random numbers, density, tail probability
x <- rnorm(50, mean=100, sd=15) # random sample of 50 IQ scores
hist(x, freq=FALSE, breaks=seq(45,155,10)) # histogram in density scale

xG <- seq(45, 155, 1) # plot theoretical density in steps of 1 IQ point
yG <- dnorm(xG, mean=100, sd=15)
lines(xG, yG, col="blue", lwd=2)

## What is the probability of an IQ score above 150?
pnorm(150, mean=100, sd=15, lower.tail=FALSE) # upper tail probabiltiy

## What does it mean to be among the bottom 25% of the population?
qnorm(.25, mean=100, sd=15) # inverse (lower) tail probability

## now do the same for a chi-squared distribution with 5 degrees of freedom
## (hint: the parameter you're looking for is df=5)
xC <- seq(0, 10, .1)
yC <- dchisq(xC, df=5)
plot(xC, yC, type="l", col="blue", lwd=2)

pchisq(10, df=5, lower.tail=FALSE) # tail probability for Z2 >= 10

## what is the appropriate rejection criterion for a variance test with alpha = 0.05?
qchisq(0.025, df=5, lower.tail=FALSE) # two-sided test for V > n * (sigma_0)^2
qchisq(0.025, df=5, lower.tail=TRUE)  # two-sided test for V < n * (sigma_0)^2

## For the following examples, let us take a reproducible sample of the population
Survey <- Census[1:100, ]

## statistical tests for variance (chi-squared) and mean (t)
x <- Survey$height # sample data
n <- length(x)

## chi-squared test for a hypothesis about the s.d. (with unknown mean)
## H0: sigma = 12 -- one-sided test against sigma > sigma0
sigma0 <- 12
S2 <- sum((x - mean(x))^2) / (n-1) # unbiased estimator of sigma^2
X2 <- (n-1) * S2 / sigma0^2        # has chi-squared distribution under H0
pchisq(X2, df=n-1, lower.tail=FALSE)

## here's a trick to carry out an approximate two-sided test (try e.g. with sigma0=20)
alt.higher <- S2 > sigma0^2
2 * pchisq(X2, df=n-1, lower.tail=!alt.higher)

## Student's t-test for a hypothesis about the mean (with unknown s.d.)
## H0: mu = 165 cm
mu0 <- 165
x.bar <- mean(x) # sample mean
s2 <- var(x)     # sample variance
t.score <- (x.bar - mu0) / (sqrt(s2 / n)) # t statistic
print(t.score) # positive indicates mu > mu_0, negative mu < mu_0
-qt(0.05 / 2, n-1) # two-sided rejection threshold for |t| at alpha = .05
2 * pt(abs(t.score), n-1, lower=FALSE) # two-sided p-value

t.test(x, mu=165) # agrees with our "manual" t-test + confidence interval

## split sample into men and women
Male <- subset(Survey, sex == "m")
Female <- subset(Survey, sex == "f")

t.test(Male$height, Female$height)
t.test(Male$height, Female$height, var.equal=TRUE) # assuming equal variances

t.test(height ~ sex, data=Survey) # "formula" interface
boxplot(height ~ sex, data=Survey)

s2.male <- var(Male$height)
s2.female <- var(Female$height)
f <- s2.male / s2.female
f

var.test(Male$height, Female$height) # F test in R

before <- runif(100, 1, 6)
after <- before - 0.5 + rnorm(100)
after <- pmax(1, pmin(after, 6))

data.frame(before=before, after=after)[1:10,]

boxplot(before, after)
t.test(before, after)

plot(before, after)
abline(0, 1)

data.frame(before=before, after=after, diff=after-before)[1:10,]

t.test(after - before, mu=0)
boxplot(before, after, after - before)

t.test(after, before, paired=TRUE) # paired t-test

## test 20 different drugs, two are significant at 5% level
dbinom(2, 20, .05) # this happens with a probability of 18.9%
pbinom(1, 20, .05, lower.tail=FALSE) # Pr( >= 2 type I errors )

## Sidak correction
pbinom(1, 20, .01, lower.tail=FALSE) # individual tests at 1% level
pbinom(1, 20, .001, lower.tail=FALSE) # individual tests at 0.1% level

prop.cint(1, 20, alternative="greater") # Bonferroni correction

summary(chickwts)
boxplot(weight ~ feed, data=chickwts)

# TukeyHSD(weight ~ feed, data=chickwts)

?aov
res <- aov(weight ~ feed, data=chickwts) # ANOVA
res
summary(res)
TukeyHSD(res)

library(languageR)

data(ratings) # size, weight and familiarity ratings
data(english) # visual lexical decision and naming latencies

# - t-test for size, familiarity, frequency of plants vs. animals
# - paired t-test for size/weight ratings

