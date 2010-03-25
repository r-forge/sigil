## SIGIL: Statistical Inference - a Gentle Introduction for Linguists
## Unit 4: Collocations & contingency tables (Part 2)
## -- code examples --

Brown <- read.delim("brown_bigrams.tbl")

## demonstration: looping over large data sets is extremely slow
attach(Brown) # just this once :-) ... and don't forget to detach below

# readers with basic Java/C/Perl/... knowledge should understand this code
result <- numeric(nrow(Brown))
for (i in 1:nrow(Brown)) {
    if ((i %% 100) == 0) cat(i, " bigrams done\n")
    A <- rbind(c(O11[i],O12[i]), c(O21[i],O22[i]))
    result[i] <- chisq.test(A)$statistic
  }

# Fisher's test is even slower:
for (i in 1:nrow(Brown)) {
    if ((i %% 100) == 0) cat(i, " bigrams done\n")
    A <- rbind(c(O11[i],O12[i]), c(O21[i],O22[i]))
    result[i] <- fisher.test(A)$p.value
  }

detach()  # -- end of temporary attachment

# first, keep R from performing integer arithmetic (always keep this in mind!!)
Brown <- transform(Brown,
  O11=as.numeric(O11), O12=as.numeric(O12),
  O21=as.numeric(O21), O22=as.numeric(O22))

Brown <- transform(Brown,
  R1=O11+O12, R2=O21+O22,
  C1=O11+O21, C2=O12+O22,
  N=O11+O12+O21+O22)

# we could also have calculated them laboriously one by one:
Brown$R1 <- Brown$O11 + Brown$O12 # etc.

# also calculate expected frequencies
Brown <- transform(Brown,
  E11=(R1*C1)/N, E12=(R1*C2)/N,
  E21=(R2*C1)/N, E22=(R2*C2)/N)

## now check that E11, ..., E22 always add up to N
summary(Brown$E11 + Brown$E12 + Brown$E21 + Brown$E22 - Brown$N)
all.equal(Brown$E11 + Brown$E12 + Brown$E21 + Brown$E22, Brown$N) # tests for "near equality" (e.g. rounding errors)


# chi-squared statistic with Yates' correction
Brown <- transform(Brown,
    chisq = N * (abs(O11*O22 - O12*O21) - N/2)^2 / (R1 * R2 * C1 * C2))

## Compare this to the output of chisq.test() for some bigrams.
## What happens if you do not apply Yates' correction?

Brown <- transform(Brown,
  logl = 2 * (
    O11*log(O11/E11) + O12*log(O12/E12) +
    O21*log(O21/E21) + O22*log(O22/E22)
  ))

summary(Brown$logl)  # do you notice anything strange?

Brown <- transform(Brown,
  logl = 2 * (
    ifelse(O11>0, O11*log(O11/E11), 0) + 
    ifelse(O12>0, O12*log(O12/E12), 0) + 
    ifelse(O21>0, O21*log(O21/E21), 0) + 
    ifelse(O22>0, O22*log(O22/E22), 0)
  ))
# ifelse() is a vectorised if-conditional

## Can you compute the association scores for MI, Dice and log-odds without peeking ahead?

Brown <- transform(Brown,
  MI = log2(O11/E11),
  Dice = 2 * O11 / (R1 + C1),
  log.odds = log( ((O11 + .5) * (O22 + .5)) / ((O12 + .5) * (O21 + .5)) )
  )

## Now check summary(Brown): are there any more NA's?

sum(Brown$chisq > qchisq(.999,df=1)) # p < .001
sum(Brown$logl > qchisq(.999,df=1))

Brown <- transform(Brown,
  r.logl = rank(-logl, ties="min"), # rank by _decreasing_ score
  r.MI   = rank(-MI, ties="min"),   # see "?rank" for "ties" parameter
  r.Dice = rank(-Dice, ties="min"))

subset(Brown, r.logl <= 20, # 20-best list for log-likelihood
  c(word1,word2,O11,logl,r.logl,r.MI,r.Dice))

## Now do the same for MI and Dice.  What are your observations?

## How many anti-collocations are there among the 100 most collocational
## bigrams according to log-likelihood?

## sorting and sort index
x <- 10 * sample(10) # 10, 20, ..., 100 in random order
sort(x)  # sorting a vector is easy (default: ascending)
sort(x, decreasing=TRUE)

# But for sorting a data frame, we need an index vector that tells us
# in what _order_ to rearrange the rows of the table.
sort.idx <- order(x)  # also has "decreasing" option
sort.idx
x[sort.idx]


## Try to sort the bigram data set accordin to log-likelihood now

sort.idx <- order(Brown$logl, decreasing=TRUE)
Brown.logl <- Brown[sort.idx, ]
Brown.logl[1:20, 1:6]

## 1) Now construct a simple character vector with the first 100 bigrams,
##    or show only relevant columns of the data frame for the first 100 rows.

## 2) Show the first 100 noun-noun bigrams (POS code "N") and
##    the first 100 adjective-noun bigrams (POS codes "J" and "N").

## 3) If you know some programming, can you write a function that
##    displays the first n bigrams for a selected association measure?

# example solution for question 1)
paste(Brown.logl$word1, Brown.logl$word2)[1:100]
paste(Brown$word1, Brown$word2)[sort.idx[1:100]]

# example solution for question 3)
# -- advanced code ahead: make your life easy with some R knowledge --
show.nbest <- function(myData, AM=c("chisq","logl","MI","Dice","O11"), n=20) {
    AM <- match.arg(AM) # allows unique abbreviations
    idx <- order(myData[[AM]], decreasing=TRUE)
    myData[idx[1:n], c("word1","word2","O11",AM)]
  }

show.nbest(Brown, "chi")

## Can you construct a table that compares the measures side-by-side?

## Evaluation of association measures for MWE extraction: Brigitte Krenn's PNV data.
# load and attach the data set
PPV <- read.delim("krenn_pp_verb.tbl")
colnames(PPV)

attach(PPV)

## You should now be able to sort the data set and calculate precision for
## some association measures and n-best lists.  (hint: sum() counts TRUE
## entries in Boolean vector)

idx.logl <- order(log.like, decreasing=TRUE)
sum(is.colloc[idx.logl[1:500]]) / 500   # \(n = 500\)
sum(is.colloc[idx.logl[1:1000]]) / 1000 # \(n = 1000\)

# use cumsum() to calculate precision for all n-best lists
prec <- cumsum(is.colloc[idx.logl]) /  (1:nrow(PPV))
prec[c(100,200,500,1000,1500,2000)]    

show.prec <- function(myData, AM, n) {
  stopifnot(AM %in% colnames(myData)) # safety first!
  sort.idx <- order(myData[[AM]], decreasing=TRUE)
  prec <- cumsum(myData$is.colloc[sort.idx]) / (1:nrow(myData))
  result <- data.frame(100 * prec[n]) # percentages
  rownames(result) <- n  # add nice row/column labels
  colnames(result) <- AM
  result  # return single-column data frame with precision values
  }

show.prec(PPV, "chisq", c(100,200,500,1000))


# data frames of same height can be combined in this way
n.list <- c(100,200,500,1000,1500,2000)
prec.table <- cbind(
    show.prec(PPV, "log.like", n.list),
    show.prec(PPV, "Fisher", n.list),
    show.prec(PPV, "chisq", n.list),
    show.prec(PPV, "chisq.corr", n.list),
    show.prec(PPV, "z.score", n.list),
    show.prec(PPV, "t.score", n.list),
    show.prec(PPV, "MI", n.list),
    show.prec(PPV, "Dice", n.list),
    show.prec(PPV, "freq", n.list)
  )

# remember the lapply / do.call trick?
prec.list <- lapply(c("log.like", "Fisher", "chisq", "chisq.corr", "z.score", "t.score", "MI", "Dice", "freq"), function (AM) show.prec(PPV, AM, n.list))
prec.table <- do.call(cbind, prec.list)


round(prec.table, 1) # rounded values are more readable

## Final challenge: plotting precision graphs

# first, generate sort index for each association measure
idx.ll <- order(log.like, decreasing=TRUE)
idx.chisq <- order(chisq, decreasing=TRUE)
idx.t <- order(t.score, decreasing=TRUE)
idx.MI <- order(MI, decreasing=TRUE)
idx.Dice <- order(Dice, decreasing=TRUE)
idx.f <- order(freq, decreasing=TRUE)

# second, calculate precision for all n-best lists
n.vals <- 1:nrow(PPV)

prec.ll <- cumsum(is.colloc[idx.ll]) * 100 /  n.vals
prec.chisq <- cumsum(is.colloc[idx.chisq]) * 100 / n.vals
prec.t <- cumsum(is.colloc[idx.t]) * 100 / n.vals
prec.MI <- cumsum(is.colloc[idx.MI]) * 100 / n.vals
prec.Dice <- cumsum(is.colloc[idx.Dice]) * 100 / n.vals
prec.f <- cumsum(is.colloc[idx.f]) * 100 / n.vals

# increase font size, set plot margins (measured in lines of text)
par(cex=1.2, mar=c(4,4,1,1)+.1)

# third: plot as line, then add lines for further measures
plot(n.vals, prec.ll, type="l", 
  ylim=c(0,42), xaxs="i", # fit x-axis range tightly
  lwd=2, col="black",     # line width and colour
  xlab="n-best list", ylab="precision (%)")
lines(n.vals, prec.chisq, lwd=2, col="blue")
lines(n.vals, prec.t, lwd=2, col="red")
lines(n.vals, prec.MI, lwd=2, col="black", lty="dashed") # line type: solid, dashed, dotted, ...
lines(n.vals, prec.Dice, lwd=2, col="blue", lty="dashed")
lines(n.vals, prec.f, lwd=2, col="red", lty="dashed")

# add horizontal line for baseline precision
abline(h = 100 * sum(is.colloc) / nrow(PPV))

# and legend with labels for the precision lines
legend("topright", inset=.05, # easy positioning of box
  bg="white", # fill legend box so it may cover other graphics
  lwd=2,      # short vectors are recycled as necessary
  col=c("black", "blue", "red"), 
  lty=c("solid","solid","solid", # no default values here!
        "dashed","dashed","dashed"),
  # either string vector, or "expression" for mathematical typesetting
  legend=expression(G^2, X^2, t, "MI", "Dice", f))
