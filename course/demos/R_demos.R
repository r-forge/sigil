##    Author: Stefan Evert
##   Purpose: Some good demos to show off R
##   Created: Sun Jul 20 19:17:45 2008
##  Modified: Sun Jul 20 20:04:04 2008 (severt)   

## basic graphics capabilities of R
demo(graphics)

## perspective plots
demo(persp)

## alternatively, show off various types of plots with lattice package
library(lattice)
demo(lattice)

## very nice interactive and animated plots with tcltk
options(demo.ask=FALSE) # make sure this doesn't mess up our demos
library(tcltk)
demo(tkcanvas)  # interactive plot - points can be dragged
demo(tkdensity) # interactive, animated plot of density estimation

## 3D graphics with RGL
library(rgl)
options(demo.ask=TRUE)
demo(rgl) # shows off most cool stuff


## demos of (generalised) linear models
library(stats)
devAskNewPage(TRUE) # demos don't seem to set this automatically
demo(lm.glm)
demo(nlm) # looks nice, but does a lot of things that are very difficult to explain
devAskNewPage(FALSE)

## V&R glm demos are nice, but without graphics and don't stop between commands
demo(glm.vr)


##
## Machine learning with Decision Trees
##

data(iris)
iris[sort(sample(150,20)), ]
summary(iris)
attach(iris)
table(Species)

plot(Petal.Length, Petal.Width, pch=20, col=as.integer(Species))

plot(jitter(Petal.Length, amount=.03), jitter(Petal.Width, amount=.03), pch=20, col=as.integer(Species))
legend(1.0, 2.5, pch=20, col=1:nlevels(Species), legend=levels(Species))

library(rpart) # built-in library for partition/classification trees

## how well does sepal length separate the different species of iris

DT <- rpart(Species ~ Sepal.Length, data=iris, method="class")
print(DT)   # summary() gives very detailed information 
printcp(DT) # complexity parameter for various pruning depths

plot(DT)
text(DT, col="blue")

plotcp(DT)

## better classification achieved with petal length and width

DT <- rpart(Species ~ Petal.Length + Petal.Width, data=iris)
print(DT)
printcp(DT)

plot(DT)
text(DT, col="blue")

plotcp(DT)

## Weka offers a variety of ML algorithms, including Decision Trees
library(RWeka)   

weka.DT <- J48(Species ~ Petal.Length + Petal.Width, iris)
print(weka.DT)

evaluate_Weka_classifier(weka.DT)

evaluate_Weka_classifier(weka.DT, numFolds=10)

plot(weka.DT) # doesn't work at the moment because required package 'coin' is missing
