## 
## Solution to Exercise #2 (frequency comparison in BNC)
## 

## -- Question 2 --

## Let us look at split infinitives ("to boldly go"), which can be identified with the simple query
##   to _{ADV} _V?I
## resulting in 4,351 hits.

## The distribution across (derived) text types is as follows (see ?BNCmeta for additional information and references)
text_types <- c("academic", "newspaper", "prose", "fiction", "misc_published", "unpublished", "spoken_other", "spoken_conversation") # ordered by assumed level of formality
f.split_inf <- c(426, 388, 712, 390, 896, 325, 976, 238)

SplitInf <- data.frame(text_type=text_types, f=f.split_inf) # represent as data frame

## -- Question 3 --

## Are word tokens a sensible unit of measurement?
##  - no, because the frequency of split infinitives per 1000 tokens would correlate strongly with the overall frequency of verbs in each text type
##  - in other words, the frequency of verbs would be a confounding factor
##  - the most appropriate unit of measurement are occurrences of infinitives with adverbs, because the speaker has a choice whether to use the split infinitive or not for each of these tokens
##  - simpler alternative: verb infinitives (with or without adverb), which are easier to identify; this leaves frequency of adverbial modification as a confounding factor

## Use a second query to obtain suitable per-category totals.
##  - for verb infinitives as unit of measurement, a suitable query is
##      to ( _{ADV} )? _V?I
##  - for verb infinitives with adverbial modifiers as unit of measurement, one might be tempted to use the following query
##      ( _{ADV} to _V?I | to _{ADV} _V?I | to _V?I _{ADV} )
##  - however, a quick inspection of the concordances shows a large number of false positives (where an initial _{ADV} is either mistagged or does not modify the infinitive)
##  - we therefore opt for the simple unit of measurement; the confounding effect of adverbial modification is accepted

## The simple query "to ( _{ADV} )? _V?I" results in 1,561,969 hits, already indicating that split infinitives are very infrequent.
## (NB: if you are not allowed to compute a distribution analysis over such a large number of hits, it's acceptable to take a random sample of 10% and multiply the observed frequency counts by 10.)
f.inf <- c(242541, 156577, 397672, 259808, 286053, 75372, 99793, 44153) # distribution across (derived) text types

SplitInf$n <- f.inf # sample size = # units of measurement for each text type

## compute relative frequency per 1000 infinitives to get an impression of the distribution
SplitInf <- transform(SplitInf, p = 1000 * f / n)
SplitInf # low rate of split infinitives in formal text types, high rate in spoke language

## optional visualization with a barplot
par(mar=c(10,4,1,1)+.1) # make room for text type labels
barplot(SplitInf$p, names=SplitInf$text_type, ylab="frequency per 1000 infinitives", las=3) # vertical orientation of labels on x-axis

## -- Question 4 --

## Perform frequency comparison tests for pairs of text types. Which differences are significant, and does the effect size make them linguistically relevant?
##  - we use prop.test() in order to obtain easily interpretable effect sizes

## newspaper vs. spoken conversation (print SplitInf to check row indices)
prop.test(SplitInf$f[c(2,8)], SplitInf$n[c(2,8)]) # ***
## effect size >= .00217 = 2.17 per 1000 infinitives is highly relevant:
## almost twice as many split infinitives in spoken conversation than in newspapers

## spoken conversations vs. other spoken material
prop.test(SplitInf$f[c(8,7)], SplitInf$n[c(8,7)]) # ***
## effect size >= .00346 = 3.46 per 1000 infinitives is also relevant

## prose vs. fiction
prop.test(SplitInf$f[c(3,4)], SplitInf$n[c(3,4)]) # **
## still significant, but effect size >= 8.76e-5  = 0.0876 per 100 infinitives is irrelevant

## academiv vs. prose
prop.test(SplitInf$f[c(1,3)], SplitInf$n[c(1,3)]) # n.s.

## NB: confidence intervals for effect sizes in the first two tests are negative because there are fewer split infinitives in the category listed first;
## always use the confidence interval closest to 0 as a conservative estimate

## -- Question 5 --

## What would be a fundamental problem of carrying out all 28 pair-wise comparisons?
##  - assume all H0 are true, i.e. there are no frequency differences between text types
##  - if we reject H0 at the customary significance level p < .05, there's still a 5% risk of a type I error for each test
##  - for 28 pairwise comparisons, we have to expect 5% * 28 = 1.4 type I error in total, i.e. a high risk that there will be at least one false positive
##  - the problem of repeated tests will be discussed in more detail in Unit 3b

## -- Question 6 --

## Construct a suitable 2 x n contingency table in order to compare all text types at once.
ct <- rbind(SplitInf$f, SplitInf$n - SplitInf$f) # recall that 2nd row = # of non-split inf rather than sample size
rownames(ct) <- c("split", "normal") # contingency tables should always be labelled
colnames(ct) <- SplitInf$text_type
ct

## What exactly is the null hypothesis of this test?
##   H0: pi_1 = pi_2 = ... = pi_8
## i.e. the relative frequency of split infinitives is equal in all 8 text types (but we don't know its precise value)

## Is there a significant difference between the categories?
chisq.test(ct)   # ***, p < 2e-16

## exact Fisher test does not work for such a large sample size,
## so we need to use a Monte Carlo simulation (B = number of simulation runs)
fisher.test(ct, simulate.p.value=TRUE, B=1e6) # ***
## note that p = 1e-6 = 1 / B, which indicates that the true p-value is even smaller (cf. chisq.test() above)

## -- Question 7 --

## The phrase "{click/V} on" is significantly more frequent in "other published material" than in any other text type. Can you think of a possible explanation?
##  - note low dispersion: 107 hits in "other published material" are concentrated in merely 6 out of 710 individual texts
##  - most hits are in HAC (59x), FT8 (21x) and CTX (13x), accounting for 93 of 107 instances
##  - quick inspection of concordance and metadata shows these texts are software manuals, giving instructions how to perform various actions with mouse clicks
