## SIGIL: Statistical Inference - a Gentle Introduction for Linguists
## Unit 3: Descriptive and inferential statistics for continuous data
## -- code examples --

library(corpora)

# generate simulated population data
Census <- simulated.census()
Wiki <- simulated.wikipedia()

# histograms
hist(Census$height, col="lightblue", breaks=50, freq=FALSE)
lines(density(Census$height), lwd=4, col="red")

hist(Census$shoe.size, col="lightblue", breaks=seq(29.75,55.25,.5), freq=FALSE)
lines(density(Census$shoe.size), lwd=4, col="red")

# discrete data!
plot(2 * table(Census$shoe.size) / nrow(Census))
lines(density(Census$shoe.size), lwd=4, col="red")

hist(Census$weight, col="lightblue", breaks=50, freq=FALSE, xlim=c(30,150), ylim=c(0,.03))
lines(density(Census$weight), lwd=4, col="red")
mu <- mean(Census$weight)
sigma <- sd(Census$weight)
data.frame(mu=mu, sigma=sigma)
abline(v=mu, lwd=3, lty="dashed")
abline(v=c(mu-2*sigma, mu+2*sigma), lwd=2)
# compare with normal approximation
t <- seq(35, 145, length.out=200)
g.t <- dnorm(t, mean=mu, sd=sigma)
lines(t, g.t, lwd=4, col="blue")

# histograms for Wackypedia 
hist(Wiki$ttr, col="lightblue", breaks=50, freq=FALSE)
lines(density(Wiki$ttr), lwd=4, col="red")

hist(Wiki$avglen, col="lightblue", breaks=50, freq=FALSE)
lines(density(Wiki$avglen), lwd=4, col="red")

x <- Wiki$types # cheat: remove outliers
sum(x >= 700)
x <- x[x < 700]
hist(x, col="lightblue", breaks=50, freq=FALSE)
lines(density(x), lwd=4, col="red")

hist(log10(x), col="lightblue", breaks=50, freq=FALSE) # log-transformed counts
lines(density(log10(x)), lwd=4, col="red")

hist(log10(x), col="lightblue", breaks=200, freq=FALSE) # type counts are really discrete

# one-sample hypothesis tests
Survey <- sample.df(Census, 100, sort=TRUE)

hist(Survey$height, breaks=20, freq=FALSE) # is normality assumption plausible?
lines(density(Survey$height), lwd=2, col="red")

qqnorm(Survey$height) # q-q plot better for small samples
qqline(Survey$height, lwd=2, col="blue")

qqnorm(Survey$weight) # skewed distribution in q-q plot
qqline(Survey$weight, lwd=2, col="blue")

qqnorm(Survey$shoe.size) # bimodal discrete distribution in q-q plot
qqline(Survey$shoe.size, lwd=2, col="blue")

# H0: mu = 165 cm
x.bar <- mean(Survey$height) # sample mean
s <- sd(Survey$height) # sample s.d.
n <- nrow(Survey) # sample size
t.score <- (x.bar - 165) / (sqrt(s^2 / n)) # t statistic
t.score
-qt(0.05 / 2, n - 1) # rejection threshold at alpha = .05

t.test(Survey$height, mu=165) # agrees with our t-test + confidence interval

# split sample into men and women
Male <- subset(Survey, sex == "m")
Female <- subset(Survey, sex == "f")

t.test(Male$height, Female$height)
t.test(Male$height, Female$height, var.equal=TRUE) # assuming equal variances

var.test(Male$height, Female$height)



library(languageR)

data(ratings) # size, weight and familiarity ratings
data(english) # visual lexical decision and naming latencies

# paired t-test for size/weight ratings
