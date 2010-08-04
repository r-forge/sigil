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

# Now take a small random sample (i.e. a survey)
Survey <- sample.df(Census, 100, sort=TRUE)

# Assessing normality: histogram & density function
x <- Survey$weight

hist(x, freq=FALSE) # example code from slides
lines(density(x))
xG <- seq(min(x),max(x),len=100)
yG <- dnorm(xG,mean(x),sd(x))
lines(xG,yG,col="red")

qqnorm(x)
qqline(x,col="red")

# histograms and Q-Q plots for some other variables
hist(Survey$height, breaks=20, freq=FALSE) # is normality assumption plausible?
lines(density(Survey$height), lwd=2, col="red")

qqnorm(Survey$height) # q-q plot better for small samples
qqline(Survey$height, lwd=2, col="blue")

qqnorm(Survey$weight) # skewed distribution in q-q plot
qqline(Survey$weight, lwd=2, col="blue")

qqnorm(Survey$shoe.size) # bimodal discrete distribution in q-q plot
qqline(Survey$shoe.size, lwd=2, col="blue")

# statistical distributions in R: random numbers, density, tail probability
x <- rnorm(50, mean=100, sd=15) # random sample of 50 IQ scores
hist(x, freq=FALSE, breaks=seq(45,155,10)) # histogram in density scale
xG <- seq(45, 155, 1) # plot theoretical density in steps of 1 IQ point
yG <- dnorm(xG, mean=100, sd=15)
lines(xG, yG, col="blue", lwd=2)

# now do the same for a chi-squared distribution with 5 degrees of freedom
# (hint: the parameter you're looking for is df=5)

pchisq(10, df=5, lower.tail=FALSE) # tail probability for >= 5

# what is the appropriate rejection criterion for a variance test with alpha = 0.05?
qchisq(0.05, df=5, lower.tail=FALSE)

# statistical tests for variance (chi-squared) and mean (t)
x <- Survey$height # sample data
n <- length(x)

# chi-squared test for a hypothesis about the s.d. (with unknown mean)
# H0: sigma = 13 -- one-sided test against sigma > sigma0
sigma0 <- 13
S2 <- sum((x - mean(x))^2) / (n-1) # unbiased estimator of sigma^2
X2 <- (n-1) * S2 / sigma0^2        # has chi-squared distribution under H0
pchisq(X2, df=n-1, lower.tail=FALSE)

# here's a trick to carry out an approximate two-sided test (try e.g. with sigma0=20)
alt.higher <- S2 > sigma0^2
2 * pchisq(X2, df=n-1, lower.tail=!alt.higher)

# Student's t-test for a hypothesis about the mean (with unknown s.d.)
# H0: mu = 165 cm
mu0 <- 165
x.bar <- mean(x) # sample mean
s2 <- sum((x - x.bar)^2) / (n-1) # sample variance
s2 <- sd(x)^2 # easier with built-in function (check that values are the same!)
t.score <- (x.bar - mu0) / (sqrt(s^2 / n)) # t statistic
print(t.score) # positive indicates mu > mu_0, negative mu < mu_0
-qt(0.05 / 2, n-1) # two-sided rejection threshold for |t| at alpha = .05

t.test(x, mu=165) # agrees with our "manual" t-test + confidence interval

# split sample into men and women
Male <- subset(Survey, sex == "m")
Female <- subset(Survey, sex == "f")

t.test(Male$height, Female$height)
t.test(Male$height, Female$height, var.equal=TRUE) # assuming equal variances

var.test(Male$height, Female$height)



library(languageR)

data(ratings) # size, weight and familiarity ratings
data(english) # visual lexical decision and naming latencies

# - t-test for size, familiarity, frequency of plants vs. animals
# - paired t-test for size/weight ratings
