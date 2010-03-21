##    Author: Stefan Evert
##   Purpose: R code examples from introduction session
##   Created: Mon Jul 21 09:48:26 2008
##  Modified: Sun Mar 21 20:28:28 2010 (severt)   

1+1
a <- 2     # assignment does \emph{not} print anything by default
a * 2
log(a)     # natural, i.e.\ base-\(e\) logarithm
log(a,2)   # base-2 logarithm

setwd("path/to/data")  # or use GUI menus
ls()                   # probably empty for now
ls                     # notice difference with previous line
quit()                 # or use GUI menus
quit(save="yes")
quit(save="no")

# NB: at least some interfaces support history recall, tab completion

a <- c(1,2,3) # \texttt{\textbf{c}} (for \emph{combine}) creates vectors
a * 2   # operators are applied to each element of a vector
log(a)  # also works for most standard functions
sum(a)  # basic vector operations: sum, length, product, \ldots
length(a)
sum(a)/length(a)

a <- 1:100            # integer sequence
a

a <- 10^(1:100)
a <- seq(from=0, to=10, by=0.1) # general sequence
a <- rnorm(100)       # 100 random numbers
a <- runif(100, 0, 5) # what you're used to from Java etc.

length(a)

summary(a)  # statistical summary of numeric vector
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02717 0.51770 1.05200 1.74300 2.32600 9.11100  

mean(a)
median(a)
sd(a)       # standard deviation is not included in summary

quantile(a)
#     0%    25%    50%    75%   100% 
# 0.0272 0.5177 1.0518 2.3261 9.1107 
quantile(a,.75)

a<-2^(1:100)        # don't forget the parentheses!
plot(a)

x<-1:100            # most often: plot \(x\) against \(y\)
y<-sqrt(x)
plot(x,y)

plot(x,a)
plot(x,a,log="y")   # various logarithmic plots
plot(x,a,log="x")
plot(x,a,log="xy")
plot(log(x),log(a))

hist(rnorm(100))    # histogram and density estimation
hist(rnorm(1000))
plot(density(rnorm(100000)))

a <- rbinom(10000,100,.5)
hist(a)

hist(a, probability=TRUE)
lines(density(a))

hist(a, probability=TRUE)
lines(density(a), col="red", lwd=3)

hist(a, probability=TRUE, main="Some Distribution", xlab="value", ylab="probability")
lines(density(a), col="red", lwd=3)

help("hist")  # R has excellent online documentation
?hist         # short, convenient form of the help command} 
help.search("histogram")
?help.search

help.start()  # searchable HTML documentation
# or use GUI menus to access \& search documentation

source("my_script.R") # more about files later
source(file.choose()) # select with file dialog box

brown <- read.table("brown.stats.txt", header=TRUE)
# if file is not in working directory, you must specify the full path
# (or use \texttt{setwd()} function we introduced before)

# exact behaviour of \texttt{file.choose()} depends on operating system
brown <- read.table(file.choose(), header=TRUE)

# more robust if you are sure file is in tab-delimited format
brown <- read.delim("brown.stats.txt")

# R can also read and write files in CSV format
write.csv(brown, "brown.stats.csv", row.names=FALSE)

# consistency check
brown.csv <- read.csv("brown.stats.csv")
all.equal(brown.csv, brown)

summary(brown)
colnames(brown)
dim(brown)       # number of rows and columns
head(brown)

plot(brown)

brown$to
head(brown$to)

# TASK: compute summary statistics (length, mean, max, etc.)
# for vectors in the Brown data frame

# what does the following do?
summary(brown$ty / brown$to)

attach(brown)   # attach data frame for convenient access
summary(ty/to)
detach()  # better to detach before you attach another frame

brown$ty[1]    # vector indexing starts with 1
brown[1,2]     # row, column

brown$ty[1:10] # use arbitrary vectors as indices
brown[1:10,2]

brown[1,]
brown[,2]
brown[brown$to < 2200, ]  # index with Boolean vector
length(mydata$ty[mydata$to >= 2200])
sum(mydata$to >= 2200)    # standard way to count matches

subset(brown, to < 2200)  # no need to attach here
lessdata <- subset(brown, to < 2200)

a <- brown$ty[brown$to >= 2200]

# equality: == (also works for strings)
# inequality: !=
# complex constraints: and &, or |, not !
# NB: always use single characters, not && or ||

## read in LOB data
lob <- read.delim("lob.stats.txt")

boxplot(brown$to,lob$to)
boxplot(brown$to,lob$to,names=c("brown","lob"))
boxplot(brown$to,lob$to,names=c("brown","lob"), ylim=c(1500,3000))
?boxplot

t.test(brown$to,lob$to)
wilcox.test(brown$to,lob$to)

brown.to.center <- brown$to[brown$to > 2200 & brown$to < 2400]
lob.to.center <- lob$to[lob$to > 2200 & lob$to < 2400]

t.test(brown.to.center, lob.to.center)

# how about sentence length?

# token and type wl are almost identical:
plot(brown$towl, brown$tywl)
cor.test(brown$towl, brown$tywl)
cor.test(brown$towl, brown$tywl, method="spearman")

# correlation with token count
plot(brown$to, brown$towl)
cor.test(brown$to, brown$towl)

