##    Author: Stefan Evert
##   Purpose: R code examples from collocations session (part 1)
##   Created: Tue Jul 22 21:59:36 2008
##  Modified: Tue Jul 22 22:30:12 2008 (severt)   

A <- matrix(c(10,47,82,956), nrow=2, ncol=2, byrow=TRUE)
print(A)

# construct matrix from row (or column) vectors
A <- rbind(c(10,47), c(82,956))

# chi-squared test is the standard independence test
chisq.test(A)

# => use test statistic as association score, p-value for interpretation
# Is there significant evidence for a collocation?

# Fisher's exact test works better for small samples and skewed tables
fisher.test(A)

# chi-squared statistic X^2 as association score
chisq.test(A)$statistic

# p-value of Fisher's test and corresponding association score
fisher.test(A)$p.value
-log10(fisher.test(A)$p.value)

# NB: chi-squared and Fisher scores are not on same scale!

# log odds ratio and conservative estimate
log(fisher.test(A)$estimate)
log(fisher.test(A)$conf.int[1])

str(fisher.test(A))  # or read help page carefully

# define two further (invented) contingency tables
B1 <- rbind(c(16,84), c(84,816))
B2 <- rbind(c(1,99), c(99,801))

# calculate chi-squared and Fisher scores for the two tables,
# as well as estimates for their log odds ratios
chisq.test(B1)$statistic
chisq.test(B2)$statistic # higher scores despite lower cooccurrence frequency and same marginals?

fisher.test(B1)$p.value
fisher.test(B2)$p.value

log(fisher.test(B1)$conf.int)
log(fisher.test(B2)$conf.int)

# Do the results look plausible to you? What is wrong?

fisher.test(B1, alternative="greater")    # high scores (significance and log odds ratio)
fisher.test(B2, alternative="greater")    # low scores (significance and log odds ratio)


Brown <- read.delim("brown_bigrams.tbl")

# Now select a number of bigrams (e.g. low and high cooccurrence frequency, or
# specific part-of-speech combinations), construct the corresponding
# contingency tables in matrix form, and calculate the different association
# scores you know.  Can you find a bigram with strong negative association?

# NB: You can use the same tests for corpus frequency comparisons.  Assume
# that a certain expression occurs 50 times in the 100,000 tokens of corpus A,
# and twice in the 1,000 tokens of corpus B.  What is an appropriate
# contingency table for these data, and what results do you obtain from the
# chi-squared and Fisher test?

