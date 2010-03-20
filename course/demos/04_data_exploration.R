##    Author: Stefan Evert
##   Purpose: R code examples from data exploration session (Italian NN compounds)
##   Created: Mon Jul 21 10:59:42 2008
##  Modified: Mon Jul 21 11:11:16 2008 (severt)   

d <- read.table("comp.stats.txt", header=TRUE)
attach(d)

# Compute basic statistics

# Look at the distribution of each cue among compounds of type attributive ("at") vs. relational ("re")

# Find out for which cues the distinction between attributive and relational is significant
# (using a t-test or Mann-Whitney ranks test)

# Also, which cues are correlated? (use "cor()" on the subset of the data-frame that contains the cues)

# cues are in columns 4 to 9


km <- kmeans(d[,4:9], 2, nstart=10)
km

# problem: extreme DELLL values dominate the clustering
# (relevant small cluster might be cluster 2 in your solution)
DELLL[km$cluster==1]
head(sort(DELLL, decreasing=TRUE))

scaled <- scale(d[,4:9])
summary(d[4:9])  # distribution of original data
summary(scaled)  # after scaling

km <- kmeans(scaled, 2, nstart=10)
km

table(km$cluster, d$TYPE) # confusion matrix

temp <- subset(d, select=c(HNPROP, NMPROP, DELLL, HDELPROP, DELMPROP, COS))
pr <- prcomp(temp, scale=TRUE)
pr

plot(pr)

biplot(pr)
biplot(pr, xlabs=TYPE, xlim=c(-.25,.25), ylim=c(-.25,.25))

plot(pr$x[,1:2], type="n", xlim=c(min(pr$x[,1]),4), ylim=c(min(pr$x[,2]),4))   # only sets up plot region
points(subset(pr$x, TYPE=="re"), col="blue", pch=19, lwd=2) # blue points for type ``re''
points(subset(pr$x, TYPE=="at"), col="red", pch=19, lwd=2)  # red points for type ``at''
legend("topright", inset=.05, fill=c("red","blue"), cex=1.5, legend=c("ATT","REL"), bg="white")     # legend explains colors

text(pr$rotation[1,1]*4, pr$rotation[1,2]*4, label="H N", cex=1.7)
text(pr$rotation[2,1]*4, pr$rotation[2,2]*4, label="N M", cex=1.7)
text(pr$rotation[3,1]*4, pr$rotation[3,2]*4, label="H DEL M", cex=1.7)
text(pr$rotation[4,1]*4, pr$rotation[4,2]*4, label="H DEL", cex=1.7)
text(pr$rotation[5,1]*4, pr$rotation[5,2]*4, label="DEL M", cex=1.7)
text(pr$rotation[6,1]*4, pr$rotation[6,2]*4, label="COS", cex=1.7)

km <- kmeans(pr$x[,1:4], 2, nstart=10)
table(km$cluster, d$TYPE)    

# what happens with more/fewer dimensions?
plot(pr$x[,1:2], type="n", xlim=c(min(pr$x[,1]),4), ylim=c(min(pr$x[,2]),4))
text(pr$x[,1], pr$x[,2], col=km$cluster, labels=TYPE)
# now refine this plot as on previous slides

