##
##  Non-randomness for passives in the Brown and LOB corpora
##
library(lattice)

screen.device <- quartz   # change this if not on Mac

## Load metadata and frequency tables for Brown and LOB
MetaBrown <- read.delim("data/brown_meta.tbl")
MetaLOB <- read.delim("data/lob_meta.tbl")

PassBrown <- read.delim("data/brown_passives.tbl")
PassLOB <- read.delim("data/lob_passives.tbl")

PassBrown <- merge(PassBrown, MetaBrown, by="id")
PassLOB <- merge(PassLOB, MetaLOB, by="id")
stopifnot(nrow(PassBrown) == 500 && nrow(PassLOB) == 500) # check that IDs have matched properly

## Load latent dimensions from POS-based multivariate analysis and DSM document vectors
LatentDim <- read.delim("data/brown_lob_svd5.tbl")
LatentDim2 <- read.delim("data/brown_lob_dsm8.tbl")
colnames(LatentDim2) <- c("corpus", "id", "cat", paste("dsm", 1:8, sep=""))
stopifnot(all(LatentDim$corpus == LatentDim2$corpus))
stopifnot(all(LatentDim$id == LatentDim2$id))
LatentDim <- cbind(LatentDim, LatentDim2[, -(1:3)])
LatentDim <- transform(LatentDim, corpus = ifelse(corpus == "brown", "AmE", "BrE")) # adjust corpus labels

## use only subsets of the corpora (to simplify illustrations)
select.cats <- c("A", "B", "E", "H", "J", "K", "M", "N", "P")
Brown2 <- subset(PassBrown, cat %in% select.cats)
LOB2 <- subset(PassLOB, cat %in% select.cats)
size2.brown <- nrow(Brown2)  # sample sizes are different now
size2.lob <- nrow(LOB2)
print(c(Brown=size2.brown, LOB=size2.lob)) # 310 vs. 312

##
## 1) Observed proportions (relative frequencies) of passives in entire corpus
##

p2.brown <- sum(Brown2$passive) / sum(Brown2$n_s)
p2.lob <- sum(LOB2$passive) / sum(LOB2$n_s)
print(round(100 * data.frame("AmE"=p2.brown, "BrE"=p2.lob, row.names="passives (%)"), 2))

##
## 2) Frequency comparison of passives: binomial test vs. t-test
##

n2.brown <- sum(Brown2$n_s)
n2.lob <- sum(LOB2$n_s)
k2.brown <- sum(Brown2$passive)
k2.lob <- sum(LOB2$passive)
p2.brown <- k2.brown / n2.brown
p2.lob <- k2.lob / n2.lob
print(data.frame(passives=c(k2.brown,k2.lob), sample_size=c(n2.brown,n2.lob), percentage=round(100 * c(p2.brown,p2.lob),1), row.names=c("AmE", "BrE")))

cont.table <- cbind(c(k2.lob, n2.lob-k2.lob), c(k2.brown, n2.brown-k2.brown))
colnames(cont.table) <- c("BrE", "AmE")
rownames(cont.table) <- c("passive", "active")
print(cont.table)

print(prop.test(c(k2.lob, k2.brown), c(n2.lob, n2.brown)))
print(chisq.test(cont.table))
print(fisher.test(cont.table))

print(t.test(LOB2$passive/LOB2$n_s, Brown2$passive/Brown2$n_s))

##
## 3) Combined data frame for both corpora + latent dimensions
##

brown.names <- c("press reportage", "press editorial", "press reviews", "religion", "skills / hobbies", "popular lore", "belles lettres", "miscellaneous", "learned", "general fiction", "detective", "science fiction", "adventure", "romance", "humour")
brown.cats <- data.frame(cat=c("A", "B", "C", "D", "E", "F", "G", "H", "J", "K", "L", "M", "N", "P", "R"),
                         name=factor(brown.names, levels=brown.names)) # make sure that factor levels for names are in correct order
                         
brown2.names <- as.character(subset(brown.cats, cat %in% select.cats)$name)
brown2.cats <- data.frame(cat=select.cats, name=factor(brown2.names, levels=brown2.names))
Brown2$corpus <- "AmE"
LOB2$corpus <- "BrE"
PassBoth2 <- rbind(Brown2, LOB2)
PassBoth2 <- transform(PassBoth2, corpus=factor(corpus), relative_frequency=100*passive/n_s)
PassBoth2 <- merge(PassBoth2, brown2.cats, by="cat")
imaginative.prose <- c("K","L","M","N","P","R")
PassBoth2 <- transform(PassBoth2, major=ifelse(cat %in% imaginative.prose, "imaginative", "informative"))

# add latent dimensions to PassBoth2 data frame
PassBoth2 <- merge(PassBoth2, LatentDim[,-3], by=c("corpus", "id"))


##
## 4) Analyse register variation with linear model
##

screen.device(width=8, height=8, bg="white")
par(cex=1.3, mar=c(4,4,2,1)+.5, oma=c(0,0,0,0), mfrow=c(2,2))

## get possible value combinations of predictor variables and observed mean relative frequencies
Predictors2 <- aggregate(data.frame(observed_mean=PassBoth2$relative_frequency), list(corpus=PassBoth2$corpus, name=PassBoth2$name), mean)

## linear models are really inappropriate, but give a nice illustration of the basic approach
LM <- lm(relative_frequency ~ name + corpus, data=PassBoth2)
print(anova(LM))  # AmE/BrE has (weakly: *) significant effect once register variation is taken into account
print(summary(LM)) # effect size and effects of individual registers (model estimates)
print(confint(LM)) # 95% confidence interval for effect size of AmE/BrE
plot(LM) # model diagnostics (residuals increase with relative frequency!)
dev.copy2pdf(file="../keynote-slides/img/AmBrE_LM_diagnostics.pdf", onefile=FALSE)

print(colSums(anova(LM))) # total variance without any predictors = 189861 (df = 621)
anova(lm(relative_frequency ~ 1, data=PassBoth2)) # another way to compute total variance

LM.predict <- data.frame(Predictors2, predicted.LM=predict(LM, newdata=Predictors2))
print(LM.predict) # compare observed mean frequencies and mean frequencies predicted by model

print(anova(lm(relative_frequency ~ corpus * name, data=PassBoth2))) # no evidence for interaction between register and language

## when is a model fit good enough? -- binomial sampling variation as baseline
binom.p <- PassBoth2$relative_frequency / 100 # assume each text is sampled from population with its MLE proportion
binom.n <- 100 # assume each test is a sample of n = 100 sentences (so passive count == relative frequency in %)
binom.var <- binom.n * binom.p * (1 - binom.p) # variance of sampling distribution for each text
print(sum(binom.var)) # we're still far from this residual variance (10200) ...

## another way to get the "baseline" residual variance from binomial sampling
PassRandom <- transform(PassBoth2, passive=rbinom(1:nrow(PassBoth2), n_s, relative_frequency/100), predictor=relative_frequency)
PassRandom <- transform(PassRandom, relative_frequency=100 * passive / n_s)
LM.random <- lm(relative_frequency ~ predictor, data=PassRandom) # this should give a perfect model (for true Exp[Y]) if coefficient == 1
print(summary(LM.random)) # coefficient is very close to 1, so residual variance is purely binomial variation
print(anova(LM.random))
## -- simulated residual variance is ca. 13000 +/- 500 (larger than theoretical approx. above)

## use latent dimensions from multivariate analysis as predictors (generalisation of "genre")
LM.svd <- lm(relative_frequency ~ name + corpus + svd1 + svd2 + svd3 + svd4 + svd5, data=PassBoth2)
print(anova(LM.svd))  # all latent dimensions are significant, AmE/BrE has become even more significant (**)
print(summary(LM.svd)) # effect sizes of genres and latent dimensions (difficult to interpret)
print(confint(LM.svd)) # effect size of AmE/BrE now at least 1.4 percent points!
plot(LM.svd) # model diagnostics look much nicer now (heteroscedasticity from binomial distribution now obvious)
dev.copy2pdf(file="../keynote-slides/img/AmBrE_LM_svd_diagnostics.pdf", onefile=FALSE)
## -- residual variance down to 47970 from 77061 (vs. 10000 from binomial sampling variation)
## -- R2 up to 74% from 59% without latent dimensions

LM.dsm <- lm(relative_frequency ~ name + corpus + dsm1 + dsm3 + dsm4 + dsm5 + dsm7 + dsm8, data=PassBoth2)
print(anova(LM.dsm)) # topical genres are worse than syntactic registers, but also considerably better than metadata only

LM.svdsm <- lm(relative_frequency ~ name + corpus + svd1 + svd2 + svd3 + svd4 + svd5 + dsm1 + dsm4 + dsm6 + dsm7 + dsm8, data=PassBoth2)
print(anova(LM.svdsm)) # only moderate improvement over latent registers, but now down to 42900 from 47970

## some significant interactions, but no really substantial improvement in (adjusted) model fit
LM.svdX <- lm(relative_frequency ~ name + corpus + svd1 + svd2 + svd3 + svd4 + svd5 + svd1 * svd3, data=PassBoth2)
print(anova(LM.svdX))  # all latent dimensions are significant, AmE/BrE has become even more significant (**)
print(summary(LM.svdX)) # effect sizes of genres and latent dimensions (difficult to interpret)
plot(LM.svdX) # model diagnostics look much nicer now (heteroscedasticity from binomial distribution now obvious)

## compute AIC for incremental linear models: are we overfitting the data?
AIC.inc <- AIC(
  lm(relative_frequency ~ corpus, data=PassBoth2),
  lm(relative_frequency ~ corpus + name, data=PassBoth2),
  lm(relative_frequency ~ corpus + name + svd1, data=PassBoth2),
  lm(relative_frequency ~ corpus + name + svd1 + svd2, data=PassBoth2),
  lm(relative_frequency ~ corpus + name + svd1 + svd2 + svd3, data=PassBoth2),
  lm(relative_frequency ~ corpus + name + svd1 + svd2 + svd3 + svd4, data=PassBoth2),
  lm(relative_frequency ~ corpus + name + svd1 + svd2 + svd3 + svd4 + svd5, data=PassBoth2),
  lm(relative_frequency ~ corpus + name + svd1 + svd2 + svd3 + svd4 + svd5 + corpus * name, data=PassBoth2),
  lm(relative_frequency ~ corpus + name + svd1 + svd2 + svd3 + svd4 + svd5 + corpus * name + svd1 * svd3, data=PassBoth2),
  lm(relative_frequency ~ corpus + name + (svd1 + svd2 + svd3 + svd4 + svd5)^2, data=PassBoth2),
  lm(relative_frequency ~ corpus * name + (svd1 + svd2 + svd3 + svd4 + svd5)^2, data=PassBoth2),
  lm(relative_frequency ~ corpus * name + svd1 * svd2 * svd3 * svd4 * svd5, data=PassBoth2),
  lm(relative_frequency ~ (corpus + name + svd1 + svd2 + svd3 + svd4 + svd5)^2, data=PassBoth2),
  lm(relative_frequency ~ corpus + name + svd1 + svd2 + svd3 + svd4 + svd5 + dsm1 + dsm2 + dsm3 + dsm4 + dsm5 + dsm6 + dsm7 + dsm8, data=PassBoth2),
  lm(relative_frequency ~ corpus * name + (svd1 + svd2 + svd3 + svd4 + svd5 + dsm1 + dsm2 + dsm3 + dsm4 + dsm5 + dsm6 + dsm7 + dsm8)^2, data=PassBoth2),
  lm(relative_frequency ~ (corpus + name + svd1 + svd2 + svd3 + svd4 + svd5 + dsm1 + dsm2 + dsm3 + dsm4 + dsm5 + dsm6 + dsm7 + dsm8)^2, data=PassBoth2)
)
rownames(AIC.inc) <- c("AmE/BrE", "AmE/BrE + Genre", "AmE/BrE + Genre + 1 SVD", "AmE/BrE + Genre + 2 SVD", "AmE/BrE + Genre + 3 SVD", "AmE/BrE + Genre + 4 SVD", "AmE/BrE + Genre + 5 SVD", "AmE/BrE * Genre + 5 SVD", "AmE/BrE * Genre + 5 SVD (1i)", "AmE/BrE + Genre + (5 SVD)^2", "AmE/BrE * Genre + (5 SVD)^2", "AmE/BrE * Genre + 5 SVD**", "(AmE/BrE + Genre + 5 SVD)^2", "AmE/BrE + Genre + 5 SVD + 8 DSM", "AmE/BrE * Genre + (5 SVD + 8 DSM)^2", "(AmE/BrE + Genre + 5 SVD + 8 DSM)^2")
print(AIC.inc) 
## -- all terms improve AIC, though very small for svd4 and corpus:name (due to large number of parameters: still explains some variance)
## -- including too many interaction terms increases AIC again -> overfitted
## -- best model has all pairwise SVD and DSM interactions, but it seems more reasonable to include just 2 important interactions
## -- latent genre dimensions and interactions added by manual exploration, no systematic AIC analysis

LM.final <- lm(relative_frequency ~ name + corpus + svd1 + svd2 + svd3 + svd4 + svd5 + dsm1 + dsm4 + dsm6 + dsm7 + dsm8 + svd1 * svd3 + dsm4 * dsm6, data=PassBoth2)
print(anova(LM.final))  # all latent dimensions are significant, AmE/BrE has become even more significant (**)
print(summary(LM.final)) # effect sizes of genres and latent dimensions (difficult to interpret)
print(confint(LM.final)) # always look at confidence intervals for parameter estimates!
plot(LM.final) # model diagnostics look much nicer now (heteroscedasticity from binomial distribution now obvious)
## -- combination of latent registers and genres brings residual variance down to 40500, R2 ca. 78%, good diagnostics
dev.copy2pdf(file="../keynote-slides/img/AmBrE_LM_final_diagnostics.pdf", onefile=FALSE)

tmp <- anova(LM.final) # determine residual variance after each step
totalVar <- sum(tmp$`Sum Sq`) 
explainedVar <- cumsum(c(0, tmp$`Sum Sq`))
print(totalVar - explainedVar) # residual variance
print(explainedVar / totalVar) # R2 (unadjusted)

print(AIC( LM, LM.svd, LM.svdsm, LM.svdX, LM.final )) # final model has good AIC with small set of parameters!

## only one text type shows a significantly different effect (fewer passives in BrE)
LM.misc <- lm(relative_frequency ~ name + corpus + svd1 + svd2 + svd3 + svd4 + svd5 + dsm1 + dsm4 + dsm6 + dsm7 + dsm8 + svd1 * svd3 + dsm4 * dsm6 + corpus * (name == "miscellaneous"), data=PassBoth2)
print(anova(LM.misc))
print(confint(LM.misc))	

LM.alltt <- lm(relative_frequency ~ name + corpus + svd1 + svd2 + svd3 + svd4 + svd5 + dsm1 + dsm4 + dsm6 + dsm7 + dsm8 + svd1 * svd3 + dsm4 * dsm6 + corpus * (name == "miscellaneous") + corpus*name, data=PassBoth2)
print(anova(LM.final, LM.misc, LM.alltt)) # no significant effects of other text types (p = .212)



##
## 5) Generalized linear models (binomial family) are more appropriate
##

## appropriate modelling of frequency data with generalized linear model
response.matrix <- cbind(PassBoth2$passive, PassBoth2$n_s - PassBoth2$passive) # response matrix with "successes" and "failures"

GLM <- glm(response.matrix ~ name * corpus, family=binomial, data=PassBoth2) # here, we get a significant interaction effect
print(anova(GLM, test="Chisq")) # both AmE/BrE and interaction are highly significant
print(summary(GLM)) # no reliable estimate for AmE/BrE effect; residual deviance is 3143.9 (df=604)
## -- significant effects for registers with unusually large or small differences between AmE and BrE
## -- negative interaction effects show registers that "go against the trend"

plot(GLM) # model diagnostics look quite reassuring now
dev.copy2pdf(file="../keynote-slides/img/AmBrE_GLM_diagnostics.pdf", onefile=FALSE)

## GLM <- glm(response.matrix ~ name + corpus : name, family=binomial, data=PassBoth2) # study AmE/BrE differences within each register
## print(anova(GLM, test="Chisq")) # overall, AmE/BrE is highly significant
## print(summary(GLM)) # significant effects in some registers, both positive and negative; residual deviance still high
## -- NB: this is essentially the same model as "corpus * name", just with a different parameterisation

## estimate "baseline" deviance for binomial sampling distribution based on simulated data
PassRandom <- transform(PassRandom, log.odds=log(predictor / (100 - predictor))) # linear component predicts log odds (for logit link)
random.response.matrix <- cbind(PassRandom$passive, PassRandom$n_s - PassRandom$passive) 
GLM.random <- glm(random.response.matrix ~ log.odds, family=binomial, data=PassRandom)
print(summary(GLM.random)) # coefficient close to 1, as it should be
plot(GLM.random) # diagnostic plot looks very good
## -- simulated residual deviance ca. 630 (df=620) with much less random variation than LM residual variance

GLM.predict <- data.frame(LM.predict, predicted.GLM=100 * predict(GLM, newdata=Predictors2, type="response"))
print(GLM.predict) # it is not obvious that GLM predictions are "better" than those of LM (but model is much more appropriate)

## include same SVD dimensions and interactions that gave the best LM fit (+ corpus-genre interaction)
GLM.svdX <- glm(response.matrix ~ corpus * name + svd1 + svd2 + svd3 + svd5 + svd1 * svd3, family=binomial, data=PassBoth2)
print(anova(GLM.svdX, test="Chisq")) # both AmE/BrE and interaction are highly significant
print(summary(GLM.svdX)) # significant estimate for AmE/BrE effect; residual deviance is 1887.9 (df=599)
plot(GLM.svdX) # model diagnostics look still good, with a few outliers

## include same SVD dimensions and interactions that gave the best LM fit (+ corpus-genre interaction)
GLM.final <- glm(response.matrix ~ name * corpus + svd1 + svd2 + svd3 + svd5 + dsm1 + dsm2 + dsm4 + dsm6 + dsm7 + dsm8 + svd1 * svd3 + dsm2 * dsm4 + dsm2 * dsm6 + dsm1 * dsm4 + dsm4 * dsm8 + dsm4 * dsm6, family=binomial, data=PassBoth2)
print(anova(GLM.final, test="Chisq")) # both AmE/BrE and interaction are highly significant
print(summary(GLM.final)) # significant estimate for AmE/BrE effect; residual deviance is 1602 (df=588); big AIC improvement

plot(GLM.final) # model diagnostics look even better than for GLM.svdX
dev.copy2pdf(file="../keynote-slides/img/AmBrE_GLM_final_diagnostics.pdf", onefile=FALSE)

dev.off()


##
## 6) Illustrate how linear models explain variation (by plotting predictors and residuals)
##

screen.device(width=12, height=8, bg="white")
par(cex=1.3, mar=c(1,4,2,1)+.5, oma=c(0,0,0,0), mfrow=c(2,1))

latent.predictor <- predict(lm(relative_frequency ~ svd1 + svd2 + svd3 + svd5 + svd1*svd3, data=PassBoth2))
PassSorted2 <- PassBoth2[ order(PassBoth2$name, PassBoth2$corpus, latent.predictor, PassBoth2$id), ]

show.predictors <- function (formula, data, show.formula=deparse(substitute(formula))) {
  LM <- lm(formula, data=data)
  x <- 1:nrow(data)
  y <- data$relative_frequency
  y.lm <- predict(LM)
  plot.sym <- ifelse(data$corpus == "AmE", 3, 20)
  plot(x, y, pch=plot.sym, cex=0.7, ylim=c(0,60), ylab="relative frequency (%)",
       main=paste("Linear model predictions (", show.formula, ")", sep=""),
       xaxs="i", yaxs="i", xaxt="n") # observed relative frequencies as points
  breaks <- tapply(x, data$name, min)
  abline(v=breaks, col="#666666") # indicate boundaries between sections of corpora
  text(breaks+3, 59, adj=c(1,1), srt=90, labels=levels(data$name), col="#008800", cex=.9, font=2) # and show register names
  line.col <- ifelse(data$corpus == "AmE", "red", "blue")
  segments(x-.5, y.lm, x+.5, y.lm, col=line.col, lwd=3) # predicted relative frequencies as coloured lines
  legend(992, 42, xjust=1, yjust=1, legend=c("AmE","BrE"), pch=c(3,20), col=c("red","blue"), lwd=3, bg="white")
  R2 <- summary(LM)$r.squared # print some goodness-of-fit information
  .tmp <- anova(LM)
  SS.resid <- .tmp$"Sum Sq"[length(.tmp$"Sum Sq")]
  print(data.frame(R2, SS.resid, row.names=show.formula))
}

show.residuals <- function (formula, data) {
  LM <- lm(formula, data=data)
  x <- 1:nrow(data)
  y <- data$relative_frequency
  y.lm <- predict(LM)
  y.resid <- y - y.lm
  plot.sym <- ifelse(data$corpus == "AmE", 3, 20)
  plot(x, y.resid, pch=plot.sym, cex=0.7, ylim=c(-30,30), ylab="residuals (%)", main=paste("Unexplained residuals of linear model"),
       xaxs="i", yaxs="i", xaxt="n") # residuals as points with "stems"
  segments(x, y.resid, x, 0, lwd=.5, col="black")
  breaks <- tapply(x, data$name, min)
  abline(v=breaks, col="#666666") # indicate boundaries between sections of corpora
  text(breaks+3, 29, adj=c(1,1), srt=90, labels=levels(data$name), col="#008800", cex=.9, font=2) # and show register names
  abline(h=0, lwd=3, col="black") # zero line for residuals
  legend(992, -28, xjust=1, yjust=0, legend=c("AmE","BrE"), pch=c(3,20), col="black", bg="white")
}

show.predictors(relative_frequency ~ 1, PassSorted2, show.formula="p ~ 1")
show.residuals(relative_frequency ~ 1, PassSorted2)
dev.copy2pdf(file="../keynote-slides/img/AmBrE_LM_visual_1.pdf", onefile=FALSE)

show.predictors(relative_frequency ~ 1 + name, PassSorted2, show.formula="p ~ 1 + genre")
show.residuals(relative_frequency ~ 1 + name, PassSorted2)
dev.copy2pdf(file="../keynote-slides/img/AmBrE_LM_visual_2.pdf", onefile=FALSE)

show.predictors(relative_frequency ~ 1 + name + corpus, PassSorted2, show.formula="p ~ 1 + genre + Am/Br")
show.residuals(relative_frequency ~ 1 + name + corpus, PassSorted2)
dev.copy2pdf(file="../keynote-slides/img/AmBrE_LM_visual_3.pdf", onefile=FALSE)

show.predictors(relative_frequency ~ 1 + name + corpus + svd1 + svd2 + svd3 + svd5 + svd1 * svd3, PassSorted2, show.formula="p ~ 1 + genre + Am/Br + latent registers")
show.residuals(relative_frequency ~ 1 + name + corpus + svd1 + svd2 + svd3 + svd5 + svd1 * svd3, PassSorted2)
dev.copy2pdf(file="../keynote-slides/img/AmBrE_LM_visual_4.pdf", onefile=FALSE)

show.predictors(relative_frequency ~ 1 + name + corpus + svd1 + svd2 + svd3 + svd4 + svd5 + dsm1 + dsm4 + dsm6 + dsm7 + dsm8 + svd1 * svd3 + dsm4 * dsm6, PassSorted2, show.formula="p ~ 1 + genre + Am/Br + latent registers + latent topics")
show.residuals(relative_frequency ~ 1 + name + corpus + svd1 + svd2 + svd3 + svd4 + svd5 + dsm1 + dsm4 + dsm6 + dsm7 + dsm8 + svd1 * svd3 + dsm4 * dsm6, PassSorted2)
dev.copy2pdf(file="../keynote-slides/img/AmBrE_LM_visual_5.pdf", onefile=FALSE)

dev.off()
