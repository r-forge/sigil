##
## Our solution to Exercise #4 (surface collocations and written/spoken keywords)
##

## 1) How to install add-on packages from CRAN

# a) use package installer in Windows / Mac OS X GUI

# b) download appropriate file from CRAN or SIGIL page,
#    then install from command-line with
#      R CMD INSTALL <package_file>

# c) install using R function install.packages()
install.packages("SIGIL")
# see ?install.packages for more options

## 2) Load library "SIGIL" and read documentation about data sets
library(SIGIL)
?BNCInChargeOf
head(BNCInChargeOf, 13) # get an impression
?BNCcomparison
head(BNCcomparison, 13)

## 3) Contingency tables for surface cooccurrences.

# According to lecture slides, first row represents word tokens within collocations spans,
# second row word tokens outside collocational spans (for given node word, here the phrase
# "in charge of"). For each collocate W, the first column of the corresponding contingency
# table contains the occurrences of W, and the second column contains all other tokens.

# This leads to the following equations for the observed frequencies O11, O12, O21, O22
COLL <- transform(BNCInChargeOf,
                  O11 = as.numeric(f.in), O12 = as.numeric(N.in - f.in),
                  O21 = as.numeric(f.out), O22 = as.numeric(N.out - f.out))
# note that we also convert the new variables to floating-point format, to avoid integer overflow

## 4) Marginal frequencies, sample size, expected frequencies, association measures.

# calculate marginal frequencies and sample size as in lecture slides
COLL <- transform(COLL,
                  R1 = O11 + O12, R2 = O21 + O22,
                  C1 = O11 + O21, C2 = O12 + O22,
                  N = O11 + O12 + O21 + O22)

summary(COLL$R1)  # should always be same value (N.in)
summary(COLL$R2)  # should always be same value (N.out)
summary(COLL$N)   # should always be same value (corpus size)

# calculate expected frequencies
COLL <- transform(COLL,
                  E11 = R1 * C1 / N, E12 = R1 * C2 / N,
                  E21 = R2 * C1 / N, E22 = R2 * C2 / N)

all.equal(COLL$E11 + COLL$E12 + COLL$E21 + COLL$E22, COLL$N) # check consistency

# calculate association scores: here we use the Dice coefficient for illustration
COLL <- transform(COLL,
                  Dice = 2 * O11 / (R1 + C1))
# (equation from www.collocations.de/AM)

# now rank data set by Dice scores (could also annotate ranks, which helps to compare different measures)
idx.Dice <- order(COLL$Dice, decreasing=TRUE)
COLL.Dice <- COLL[idx.Dice, ]
head(COLL.Dice[, c("collocate", "f.in", "f.out", "Dice")], 20) # select only relevant columns for better readability

# TMTOWDTI -- lecture slides show other ways for sorting the data set, calculating rankings and extracting n-best lists

## 5) + 6) Keyword identification.

# written and spoken sample sizes (NB: "OTHER" entry for all other nouns allows us to do this)
N.written <- sum(BNCcomparison$written)
N.spoken <- sum(BNCcomparison$spoken)

# Contingency table for frequency comparison, as explained in lecture slides:
#  - first column contains data for spoken sample, second column for written sample
#  - first row contains frequency counts for given noun in spoken/written sample
#  - second row contains sample size - frequency count
# NB: we have swapped order of columns, because some association measures work better
#     if the smaller sample is in the first column (as is the case for cooccurrence data)
KEY <- transform(BNCcomparison,
                 O11 = as.numeric(spoken), O21 = as.numeric(N.spoken - spoken),
                 O12 = as.numeric(written), O22 = as.numeric(N.written - written))

## 7) Marginal frequencies, N, expected frequencies, association measures (as above)

# calculate marginal frequencies and "total sample size" N as above
KEY <- transform(KEY,
                 R1 = O11 + O12, R2 = O21 + O22,
                 C1 = O11 + O21, C2 = O12 + O22,
                 N = O11 + O12 + O21 + O22)

summary(KEY$C1)  # should always be same value (spoken sample size)
summary(KEY$C2)  # should always be same value (written sample size)
summary(KEY$N)   # should always be same value ("total sample size")

# calculate expected frequencies
KEY <- transform(KEY,
                 E11 = R1 * C1 / N, E12 = R1 * C2 / N,
                 E21 = R2 * C1 / N, E22 = R2 * C2 / N)

all.equal(KEY$E11 + KEY$E12 + KEY$E21 + KEY$E22, KEY$N) # check consistency

# we don't really want to calculate keyness for the "OTHER" entry, so delete it
KEY <- subset(KEY, noun != "OTHER")

# calculate association scores: here we use MI as some researchers did in terminology extraction
KEY <- transform(KEY,
                 MI = log2(O11 / E11))
# (equation from www.collocations.de/AM)

# now rank data set by MI scores and list entries with
#  - large positive scores (spoken keywords, i.e. higher relative frequency in first column)
#  - large negative scores (written keywords, i.e. higher relative frequency in second column)
idx.MI <- order(KEY$MI, decreasing=TRUE)
KEY.MI <- KEY[idx.MI, ]
head(KEY.MI[, c("noun", "written", "spoken", "O11", "E11", "MI")], 20) 
tail(KEY.MI[, c("noun", "written", "spoken", "O11", "E11", "MI")], 20) 
