## SIGIL: Statistical Inference - a Gentle Introduction for Linguists
## Unit 8: The non-randomness of corpus data & generalised linear models
## -- code examples --

# read example data (passive counts on all texts from 9 genre categories in Brown/LOB)
Brown <- read.delim("passives_by_text.brown.tbl")
LOB <- read.delim("passives_by_text.lob.tbl")

n.B <- nrow(Brown) # 310 texts
n.L <- nrow(LOB)   # 312 texts (two more in skills/hobbies)

n.sent <- 100      # assume each text consists of 100 sentencs (simplification)
p.B <- sum(Brown$passive) / sum(Brown$n_s)  # overall proportion of passives in Brown
p.L <- sum(LOB$passive) / sum(LOB$n_s)      # and LOB

# plot empirical distribution of passive counts across Brown texts
# (NB: assumes that all texts contain the same number of sentences, which is not true ...)
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

# perform chi-squared test for pooled data
passives.B <- sum(Brown$passive) # pool Brown data
n_s.B <- sum(Brown$n_s)
passives.L <- sum(LOB$passive)   # pool LOB data
n_s.L <- sum(LOB$n_s)
prop.test(c(passives.L, passives.B), c(n_s.L, n_s.B)) # executes chi-squared test

# perform t-test for per-text data (measurement: relative frequency of passives)
relfreq.L <- LOB$passive / LOB$n_s
relfreq.B <- Brown$passive / Brown$n_s
t.test(relfreq.L, relfreq.B)

