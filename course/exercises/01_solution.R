## 
## Solution to Exercise #1 (analysing the BNC metadata)
## 

## TIMTOWTDI ;-)
## "there is more than one way to do it"

## load the "SIGIL" package and attach the BNC metadata table

library(SIGIL)
nrow(BNCmeta)   # check that we can access the data
attach(BNCmeta) # so we can directly refe to the columns by name


## -- Question 1 --

dim(BNCmeta)

?BNCmeta # to display the help page, with brief explanations of variables

## How many metadata? It depends on how you define metadata.
##  - If you consider each column, id's included, to be metadata, then the number of
##    metadata variables is the same as the number of columns, i.e., 31.
##  - A corpus linguist might prefer to say that there are 31 variables, comprising
##    1 unique ID, 5 text statistics (n_words etc.) and 25 metadata items


## -- Question 2 --

## how many genres are there?
length(levels(genre)) 
## you could also use(summary(genre)) or length(table(genre))

## the "smallest" genres (fewest texts) 
head(sort(table(genre)))

counts <- sort(table(genre)) # pick out the name of the smallest genre directly
names(counts)[1] # returns string "W:fict:drama"


## -- Question 3 --

## it is good practice to ensure that files are written in UTF-8 encoding
write.csv(BNCmeta, file="bncmeta.csv", row.names=FALSE, fileEncoding="UTF-8")

## On a German Windows PC, you will have to use write.csv2() because
## Excel expects a different "continental" CSV format
write.csv2(BNCmeta, "bncmeta.csv", row.names=FALSE)


## -- Question 4 --

## summarise the distribution of text lengths
summary(n_words)
summary(n_s)

boxplot(n_words) # visual summary
boxplot(n_words, ylim=c(0,100000), yaxs="i") # zoom in to box + whiskers

## outliers in written-to-be-spoken
res <- boxplot(n_words[text_type == "written-to-be-spoken"], ylim=c(0, 50000), yaxs="i")
res$out # numerical values of the outlier points

## list the titles of the outlier texts
##  - with threshold obtained visually from boxplot
title[ text_type == "written-to-be-spoken" & n_words <= 20000 ]

##  - select exactly the outlier values using set operator %in%
title[ text_type == "written-to-be-spoken" & n_words %in% res$out ]


## -- Question 5 --

## do text lengths differ between text types?
boxplot(n_words ~ text_type, ylim=c(0,100000), col="lightgray") # filled boxes are easier to compare

## do they differ between male and female authors?
boxplot(n_words ~ author_sex, notch=TRUE, col="lightgray", ylim=c(0, 75000)) # if notches overlap, difference is usually not significant

## alternative solution:
boxplot(n_words[author_sex=="male"], n_words[author_sex=="female"], names=c("male","female"))


## -- Question 6 --

## produce a subset of the metadata table containing only texts for which author sex is known
temp <- subset(BNCmeta, author_sex=="male" | author_sex=="female")

## alternative version for Python programmers: 
temp <- subset(BNCmeta, author_sex %in% c("male", "female"))

## omit the title and irrelevant metadata columns (esp. those which have only a single value in the subset)
## let's check which columns have this property:
summary(temp)

## they appear to be:
## context interaction_type respondent_age respondent_class respondent_sex region mode

## subsetting with syntactic sugar:
male.female <- subset(temp, select=-c(title, context, interaction_type, mode, region, respondent_age, respondent_class, respondent_sex))

## because context ... region are adjacent in this order, we can shorten to range context:region
male.female <- subset(temp, select=-c(title, context:region, mode))


## -- Question 7 --

## how good is the linear correlation between text length in words and sentences?
plot(n_words, n_s, xlim=c(0, 75000), ylim=c(0, 7500)) # suitable range found by trial & error

## we can also visualize overlaps with translucent points ("alpha" channel of colour)
plot(n_words, n_s, xlim=c(0, 75000), ylim=c(0, 7500), pch=20, col="#0000AA30")

## Pearson correaltion r and confidence interval
cor.test(n_words, n_s)
##   corr coeff:    0.8687635
##   conf interval: 0.8610003 .. 0.8761220 

## The relation does not appear to be linear: plot strongly suggests that different sets
## of documents follow different patterns, some linerar, some not.
