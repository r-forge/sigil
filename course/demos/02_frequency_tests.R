## SIGIL: Statistical Inference - a Gentle Introduction for Linguists
## Unit 2: Corpus frequency data & statistical inference
## -- code examples --

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

# nicer version of user function with genre category labels do.test <- function (k1, n1, k2, n2, cat="") {
	res <- prop.test(c(k1, k2), c(n1, n2))
	data.frame(
		p=res$p.value,
		lower=100*res$conf.int[1], # scaled to % points
		upper=100*res$conf.int[2],
		row.names=cat   # add genre as row label
	) # return data frame directly without local variable fmt
}

do.test(	Brown$passive[15], Brown$n_s[15],
			LOB$passive[15], LOB$n_s[15], cat=Brown$name[15])

# ad-hoc convenience function to reduce typing/editing
# (works only if global Brown/LOB variables are set correctly!)
quick.test <- function (i) {
		do.test(	k1=Brown$passive[i], n1=Brown$n_s[i],
					k2=LOB$passive[i], n2=LOB$n_s[i],
					cat=Brown$name[i] )
}

quick.test(15)
quick.test(9)

# loop over all 15 rows (or use 1:nrow(Brown) to count rows automatically)
for (i in 1:15) {
	print( quick.test(i) )
}

# our code relies on same ordering of genres in Brown and LOB
all(Brown$cat == LOB$cat) 

# we need some R wizardry to collect the 15 genre results into a single table
#  - lapply() applies a function to each element of a list or vector
#  - returns list with result values (here: table rows)
#  - use rbind() to combine rows into single table
#  - rows have to be passed as individual arguments using do.call()
res.list <- lapply(1:15, quick.test)

res <- do.call(rbind, res.list) # as if each list element were an argument

res
round(res, 3) # easier to read after rounding

# Which differences are significant? Are the effect sizes linguistically relevant?

# Use cont.table() from corpora package to create list of contingency tables,
# then lapply() or sapply() Fisher's test and chi-squared directly to each table.

