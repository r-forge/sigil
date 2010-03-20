##
## Hypothesis tests for frequency comparison and simple script automation.
##

# relevant information from output of binom.test() is the p-value
binom.test(19, 100, p=.15) # two-sided
binom.test(19, 100, p=.15, alternative="greater") # one-sided

binom.test(19, 100, p=.15)$p.value
binom.test(23, 100, p=.15)$p.value
binom.test(190, 1000, p=.15)$p.value

# focus on confidence interval now (no null hypothesis specified)
binom.test(190, 1000, conf.level=.99)

library(corpora) # package 'corpora' can be installed from CRAN
prop.cint(190, 1000, conf.level=.99)
?prop.cint # "conf. intervals for proportions"

# plot confidence interval curves (e.g. for choosing sample size)
plot(0:1, 0:1, type="n", xaxs="i", yaxs="i", xlab="sample proportion p", ylab=expression("population  value " * pi), main="99% confidence intervals")
abline(h=seq(0,1,.1), v=seq(0,1,.1), col="grey") # background grid
abline(0, 1, lwd=2, col="black") # maximum-likelihood estimate = MLE
n <- 50 # sample size n = 50
k <- 0:n
conf <- prop.cint(k, n, conf.level=.99) # ignore spurious warnings
lines(k/n, conf$lower, lwd=2, col="blue") # curve of lower limits
lines(k/n, conf$upper, lwd=2, col="blue") # curve of upper limits
n <- 200 # sample size n = 200
k <- 0:n
conf <- prop.cint(k, n, conf.level=.99)
lines(k/n, conf$lower, lwd=2, col="red")
lines(k/n, conf$upper, lwd=2, col="red")
legend("topleft", inset=.05, bg="white", lwd=3, col=c("black","red","blue"), legend=c("MLE", "n - 200", "n = 50"))

# frequency comparison: proportions test, chi-squared, Fisher's test
prop.test(c(19,25), c(100,200)) # notice conf. int. for difference of proportions

ct <- cbind(c(19,81), c(25,175)) # construct contingency table by columns
ct
chisq.test(ct)  # chi-squared test (no confidence interval)
fisher.test(ct) # Fisher's exact test (conf. int. for odds ratio)

ct <- cont.table(19, 100, 25, 200) # lazy solution from 'corpora' library

## Now, perform frequency comparison for number of passives in Brown (AmE)
## and LOB (BrE) corpora (both pooled data and for each section / genre).

Brown <- read.csv("passives.brown.csv") # load CSV data files
LOB <- read.csv("passives.lob.csv")

Brown # table format should be self-explanatory
LOB

Brown.all <- colSums(Brown[, 2:4]) # pooled data = column sums
LOB.all <- colSums(LOB[, 2:4])

Brown.all # now copy relevant data into contingency table
LOB.all

ct <- cbind(c(10123, 49576-10123), c(10934, 49742-10934)) # Brown, LOB
ct
fisher.test(ct)

prop.test(c(10123, 10934), c(49576, 49742)) # effect size = difference of proportions

# we can do the same (manually) for every genre in the Brown / LOB corpus;
# convenient: write our own function that shows relevant information concisely
do.test <- function (k1, n1, k2, n2) {
	res <- prop.test(c(k1, k2), c(n1, n2)) # result of proportions test
	# data frames are a nice way to display summary tables
	fmt <- data.frame(p=res$p.value, lower=res$conf.int[1], upper=res$conf.int[2])
	fmt # return value of function = last expression
}

do.test(10123, 49576, 10934, 49742) # for the pooled data again
do.test(146, 975, 134, 947) # humour genre (category R)

# we can also extract these values directly from the data frames
do.test(Brown$passive[15], Brown$n_s[15], LOB$passive[15], LOB$n_s[15])

# nicer version of do.test() function labels result with genre code
do.test <- function (k1, n1, k2, n2, cat="") {
	res <- prop.test(c(k1, k2), c(n1, n2)) # result of proportions test
	# data frames are a nice way to display summary tables
	fmt <- data.frame(p=res$p.value, lower=res$conf.int[1], upper=res$conf.int[2])
	rownames(fmt) <- cat
	fmt # return value of function = last expression
}

do.test(Brown$passive[15], Brown$n_s[15], LOB$passive[15], LOB$n_s[15], cat=Brown$cat[15])

all(Brown$cat == LOB$cat) # our code relies on same ordering of genres in Brown and LOB

# we can carry out tests for all genres with a 'for' loop
for (i in 1:15) {
	print(do.test(Brown$passive[i], Brown$n_s[i], LOB$passive[i], LOB$n_s[i], cat=Brown$cat[i]))
}

# it would be nice to collect all these rows into a single table;
# you can either iteratively build this table, or use advanced R skills
result.list <- lapply(1:15, function (i) {
	do.test(Brown$passive[i], Brown$n_s[i], LOB$passive[i], LOB$n_s[i], cat=Brown$name[i])
})
result <- do.call(rbind, result.list) # think of this as an idiom you have to remember ...

result
round(result, 5) # non-engineering notation is easier to read

# which differences are significant? are the effect sizes linguistically relevant?

# Homework: Extend the do.test() function so that the two sample proportions are
# also shown in the table. Do you need to modify any of the other code as well?


