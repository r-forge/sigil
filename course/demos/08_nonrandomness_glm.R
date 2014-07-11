## SIGIL: Statistical Inference - a Gentle Introduction for Linguists
## Unit 8: The non-randomness of corpus data & generalised linear models
## -- code examples --

## read example data (passive counts on all texts from 9 genre categories in Brown/LOB)
Passives <- read.delim("passives_by_text.tbl")

## data set is also included in the SIGIL package
library(SIGIL)
Passives <- BrownLOBPassives

## display 10 random rows to get an idea of the table layout
Passives[sample(nrow(Passives), 10), ]

## add relative frequency of passives in each file (percentage)
Passives <- transform(Passives, relfreq = 100 * passive / n_s)

## split into Brown and LOB sentences
Brown <- subset(Passives, lang == "AmE")
LOB <- subset(Passives, lang == "BrE")

n.B <- nrow(Brown) # 310 texts
n.L <- nrow(LOB)   # 312 texts (two more in skills/hobbies)

n.sent <- 100      # assume each text consists of 100 sentencs (simplification)
p.B <- sum(Brown$passive) / sum(Brown$n_s)  # overall proportion of passives in Brown
p.L <- sum(LOB$passive) / sum(LOB$n_s)      # and LOB

#### plot empirical distribution of passive counts across Brown texts
## (NB: assumes that all texts contain the same number of sentences, which is not true ...)
my.k <- 0:60 # range of passive counts (for comparison against binomial)
my.breaks <- seq(from=0.5, to=60.5, by=1) # bin "breaks" are between integer counts

hist(Brown$passive, breaks=my.breaks, col="grey", ylim=c(0,30), 
	xlab="passive count", ylab="number of files", main="Brown")
my.expected <- n.B * dbinom(my.k, size=n.sent, prob=p.B)
lines(my.k, my.expected, col="red", lwd=2)

hist(LOB$passive, breaks=my.breaks, col="grey", ylim=c(0,30),
	xlab="passive count", ylab="number of files", main="LOB")
my.expected <- n.L * dbinom(my.k, size=n.sent, prob=p.L)
lines(my.k, my.expected, col="red", lwd=2)

#### perform chi-squared test for pooled data
passives.B <- sum(Brown$passive) # pool Brown data
n_s.B <- sum(Brown$n_s)
passives.L <- sum(LOB$passive)   # pool LOB data
n_s.L <- sum(LOB$n_s)
prop.test(c(passives.L, passives.B), c(n_s.L, n_s.B)) # executes chi-squared test

#### perform t-test for per-text data (measurement: relative frequency of passives)
t.test(LOB$relfreq, Brown$relfreq)

## alternative: use formula interface to t.test (relfreq against language variety)
t.test(relfreq ~ lang, data=Passives) # this puts AmE first! (why?)

## visualize the underlying data with boxplot
boxplot(relfreq ~ lang, notch=TRUE, data=Passives, ylim=c(0,100))

#### side-by-side boxplots for each text type can be created with the "lattice" package
library(lattice)

bwplot(relfreq ~ lang | genre, data=Passives) # bw = "Box and Whiskers"
## NB: it's much trickier to customise lattice plots than standard plots

### linear model based on genre and AmE/BrE (without interaction)
LM <- lm(relfreq ~ genre + lang, data=Passives) # response ~ predictor1 + predictor2 + …
anova(LM) # analysis of variance: do factors significantly improve the model?
## see ?anova.lm for details (?anova is not very informative)

summary(LM) # shows coefficients for each genre / language variety
confint(LM)
## NB: difficult to interpret, because arbitrary factor level is chosen as baseline

## more intuitive effects: obtain model predictions for each genre in each language variety
## we need to construct a dummy data set for which predictions are made
Predictions <- unique(Passives[, c("genre", "lang")]) # all unique combinations of genre + AmE/BrE
Predictions <- Predictions[order(Predictions$genre, Predictions$lang), ]
transform(Predictions, predicted=predict(LM, newdata=Predictions))
cbind(Predictions, predict(LM, newdata=Predictions, interval="confidence"))
cbind(Predictions, predict(LM, newdata=Predictions, interval="prediction"))

## alternative: partial effects (using add-on package "effects")
library(effects)
eff <- Effect(c("genre", "lang"), LM)
print(eff)
plot(eff, multiline=TRUE, rotx=30)

topR <- list(x=0.95, y=0.95, corner=c(1,1)) # some display options
plot(eff, multiline=TRUE, rotx=30, ci.style="bars", key.args=topR)

plot(Effect("genre", LM), rotx=30, ylim=c(0,50))
plot(Effect("lang", LM), rotx=30, ylim=c(0,50))

#### include genre/variety interaction term: not significant
LM <- lm(relfreq ~ genre + lang + genre:lang, data=Passives)
anova(LM)

#### should always look at diagnostic plot (4 panels)
par(mfrow=c(2,2))
plot(LM)
par(mfrow=c(1,1))


### GLM is more appropriate: predictions transformed to [0,1] range, binomial sampling component
## note that we need some additional options below to get the intended results!

## responses are now pairs of passive/active counts (k, n-k) = "successes" / "failures"
## we collect them in a response matrix where each row corresponds to one text
response.matrix <- cbind(Passives$passive, Passives$n_s - Passives$passive)

## genre * lang is a shorthand for "main effects" + interactions = genre + lang + genre:lang
GLM <- glm(response.matrix ~ genre * lang, family="binomial", data=Passives) # include interaction term

anova(GLM, test="Chisq") # interaction is also significant in this case

summary(GLM) # even more difficult to interpret
confint(GLM)

#### show predictions for each genre and language variety
transform(Predictions, predicted=100*predict(GLM, type="response", newdata=Predictions))
res <- predict(GLM, type="response", newdata=Predictions, se.fit=TRUE) # predictions with standard errors
transform(Predictions, predicted=100*res$fit, lwr=100*(res$fit-1.96*res$se.fit), upr=100*(res$fit+1.96*res$se.fit))
## predictions look better than for LM, narrower confidence intervals

#### We can't compute prediction intervals for new texts -- why?

## can also plot partial effects (note the significant interaction here!)
eff <- Effect(c("genre", "lang"), GLM)
print(eff)
plot(eff, multiline=TRUE, rotx=30, ci.style="bars", key.args=topR)

plot(Effect("genre", GLM), rotx=30, ylim=c(-3,0))
plot(Effect("lang", GLM), rotx=30, ylim=c(-3,0))

## don't forget to look at model diagnostics
par(mfrow=c(2,2))
plot(GLM) # much nicer than LM diagnostics
par(mfrow=c(1,1))






