##
## Some examples and illustrations of regression & linear models
## Stefan Evert, 5 Feb 2009
##

## open graphics device with suitably scaled fonts and tight margins for screen presentation
dev.new()
par(cex=1.3, mar=c(4.1, 4.1, 0.5, 0.5))

##
## we need to be able to generate synthetic random data for two variables Y, X
## with a specified functional relationship (Y = f(X)) plus Gaussian random errors
##

## f ... relationship between X and Y: Y = f(X) + E [default: Y = X + E]
## n ... number of random samples
## sd ... standard deviation of Gaussian random error E
## x.min, x.max ... range of values for X (random uniform distribution)
## sort=TRUE ... order data points by X (increasing)
random.data <- function (f=function (x) x, n=50, sd=1, x.min=-2, x.max=2, sort=FALSE) {
	x <- runif(n, x.min, x.max)
	if (sort) x <- sort(X)
	y <- f(x) + rnorm(n, 0, sd)
	data.frame(x=x, y=y)
}

##
## simple linear regression with lm()
##

S <- random.data() # generate random data set (50 points)
plot(S$x, S$y, type="p", pch=20, xlab="X", ylab="Y") # always a good idea to plot data first

## linear regression uses lm() function, e.g. Y against X with y ~ x
R <- lm(y ~ x, data=S)
print(R)  # prints regression coefficients
abline(R$coefficients, lwd=2)  # plot regression line

## correlation coefficient has to be calculated with summary() method
print(summary(R)$r.squared)  # correlation coefficient = % of variance "explained"

##
## illustrate properties of linear regression with some examples on random data
##

## function to plot data set D with regression line
##  - regression for Y ~ X; if both=TRUE, also show inverse regression for X ~ Y
do.plot <- function (D, xlim=c(-2,2), ylim=c(-2,2), both=FALSE) {
	plot(D$x, D$y, type="p", pch=20, xlab="X", ylab="Y", xlim=xlim, ylim=ylim)

	## perform standard regression, plot regression line, obtain correlation
	regression <- lm(y ~ x, data=D) 
	abline(regression$coefficients, col="blue", lwd=2)
	r2 <- summary(regression)$r.squared

	## display correlation coefficient as legend box in bottom right corner
	## NB: bquote allows you to interpolate variables and expressions into a mathematical formula
	## (this is fairly advanced R stuff, look at ?plotmath for examples)
	legend("bottomright", inset=.02, bg="white", legend = bquote(R^2 == .(round(r2, 4))) )

	if (both) {
		## perform inverse regression (X ~ Y) and plot as above
		regression <- lm(x ~ y, data=D)
		ab <- regression$coefficients
		## regression line is X = a + bY --> Y = -(a/b) + (1/b)X
		abline(-ab[1]/ab[2], 1/ab[2], col="red", lwd=2)
		legend("topright", inset=.02, bg="white",
			lwd=3, col=c("blue", "red"), legend=c("Y ~ X", "X ~ Y"))
	}
}

## regression for linear relationship with different slopes (but same random error)
## (note how we pass an anonymous function to random.data() to specify relationship Y = f(X) + E)
do.plot(random.data(function (x) x)) # same as above: y = f(x) = x
do.plot(random.data(function (x) .5*x)) # y = f(x) = .5 * x
do.plot(random.data(function (x) .2*x)) # y = f(x) = .2 * x
do.plot(random.data(function (x) 2*x)) # y = f(x) = 2 * x
## --> correlation depends on ratio between slope (of linear relation) and amount of error (variance)

## here is the same effect for different amounts of random error (s.d. of Gaussian error)
do.plot(random.data(function (x) .3*x, sd=.2))
do.plot(random.data(function (x) .3*x, sd=.5))
do.plot(random.data(function (x) .3*x, sd=1))
do.plot(random.data(function (x) .3*x, sd=2))
## rerun each command several times to see how reliable regression line is!

## strong correlation --> regression lines for Y ~ X and X ~ Y are similar
do.plot(random.data(function (x) .3*x, sd=.2), both=TRUE)
do.plot(random.data(function (x) .3*x, sd=.5), both=TRUE)
do.plot(random.data(function (x) .3*x, sd=1), both=TRUE)

##
## the behaviour of linear regression for non-linear relationships is unpredictable
##

## monotonic non-linear function --> high correlation, regression line captures trend
## (hyperbolic asymptote, arbitrarily chosen for illustration purposes)
do.plot(random.data(function (x) 2 - 4/(x+3), sd=.5), both=FALSE)
do.plot(random.data(function (x) 2 - 4/(x+3), sd=.15), both=FALSE)
do.plot(random.data(function (x) 2 - 4/(x+3), sd=.15), both=TRUE)

## symmetric non-linear function (parabola) --> no linear correlation detected
do.plot(random.data(function (x) x^2 - 2, sd=.4), both=FALSE)
do.plot(random.data(function (x) x^2 - 2, sd=.4), both=TRUE)

##
## polynomial regression as a special case of the general linear model
##

## obtain random points S with parabolic relationship between X and Y
S <- random.data(function (x) x^2 - 2, sd=.4)
plot(S$x, S$y, type="p", pch=20, xlab="X", ylab="Y") # redo our standard regression plot
abline(lm(y ~ x, S)$coefficients, lwd=2, col="blue")

## for polynomial regression, we need to include a term x^2 = x*x in the model;
## since "*" has a special meaning in model formulae, we have to "protect" it with I()
R <- lm(y ~ x + I(x*x), data=S) # quadratic regression = 2-nd degree polynomial
print(R)
print(summary(R)$r.squared)
## --> polynomial regression achieves very good fit with R^2 > 90%

## we cannot use abline to plot the non-linear regression line (because it is not straight);
## the easiest and most general solution is to obtain model predictions for different values of X
L <- data.frame(x=seq(-2, 2, .1)) # L is a data frame with x- and y-coordinates of points on the line
L$y <- predict(R, newdata=L) # y-coordinates are model predictions for the specified x-coordinates
lines(L$x, L$y, lwd=2, col="red") # draw polynomial regression line in red
legend("top", inset=.05, bg="white", lwd=2, col=c("blue","red"), legend=c("linear","polynomial"))

## Exercise:
## Modify the code above to perform cubic polynomial regression (with 3-rd degree polynomial).
## Which additional terms do you need to include in the model formula?  Use summary(R) and 
## anova(R) to find out which of the regression coefficients are significant.  What do you expect?
## Does the regression line look plausible?
## [see below for a short explanation of summary(), anova() and predict() functions]

##
## multiple regression and linear models (U.S. cars in 1973-74)
##

data(mtcars) # a popular data set for simple demos
?mtcars # always find out what exactly the data set contains!

## before you perform a linear regression analysis, it's always a good idea to plot
## the data and look for obvious correlations and nonlinearities
pairs(~ disp + hp, mtcars) # displacement and power should be correlated

## can we predict fuel consumption (mpg) from displacement, power, weight, ... ?
pairs(~ mpg + disp + hp + wt, mtcars)

## we can define our own panel function e.g. to add smoothed trend lines to panels below diagonal
my.panel <- function (x, y, ...) {
	points(x, y, pch=20, ...) # this is what the default panel function does (except for pch=)
	lines(lowess(x,y), lwd=2, col="blue") # trend line from LOWESS smoother
}
pairs(~ mpg + disp + hp + wt, mtcars, lower.panel=my.panel)

## we should be able to predict acceleration from power and weight
pairs(~ qsec + hp + wt, mtcars, lower.panel=my.panel)

## regression of fuel consumption (mpg = miles per gallon) on power (hp) and weight (wt)
##  - 1 stands for the implicit intercept term and is usually omitted from the specification
##  - note how 1 + hp + wt corresponds to the columns of the design matrix
##  - this model assumes independent, additive linear effects from hp and wt!
R <- lm(mpg ~ 1 + hp + wt, data=mtcars)
print(R) # information about linear regression
## --> increasing hp and/or wt reduces mpg, but the coefficient for wt is much larger
## --> NB: coefficients (= slope) depend on measurement units and DO NOT show correlation strength!

## the summary() of a linear regression provides information about the statistical linear model
summary(R)
## --> shows that mpg is "explained" quite well by hp and wt (R^2 = 82.68%)
## --> standard errors for regression coefficients (as estimators = random variables)
## --> t-tests show whether each regression coefficient has a significant effect (i.e. is non-zero)

## analysis of variance (ANOVA) describes significance of predictor variables in an intuitive way; 
## it shows what proportion of the total variance (sum of squares) of mpg each predictor "explains"
anova(R)
sum((mtcars$mpg - mean(mtcars$mpg))^2) # total sum of squares of the mpg variable for comparison
## --> hp "explains" 678 out of 1126, wt another 252; residual error (ESS) is 195.05 
## --> F-tests compare nested models: mpg ~ 1 vs. mpg ~ 1 + hp vs. mpg ~ 1 + hp + wt
## --> all predictors have a highly significant effect

## NB: ANOVA results depend on ordering of predictors
anova(lm(mpg ~ 1 + wt + hp, data=mtcars))
## --> now wt explains the largest part of the variance, though complete model is the same as above
## --> weight is a much better predictor of fuel consumption by itself than power

## check for interaction between weight and power (wt:hp stands for product weight * power)
R <- lm(mpg ~ 1 + hp + wt + hp:wt, data=mtcars)
R <- lm(mpg ~ 1 + hp * wt, data=mtcars) # hp * wt is short-hand for terms plus interaction
anova(R)
## --> highly significant interaction explains another 65 of sum of squared
summary(R)
## --> small positive coefficient of hp:wt interaction: what does this imply?

## confidence intervals for (individual) regression coefficients
confint(R)

## we can now predict fuel consumption of new models of cars
##  - for illustration, predict Mercedes models from model trained on remaining cars
##  - to do so, we split mtcars into test (Mercedes) and training (all other cars) data
mtcars.test <- mtcars[8:14, ]
mtcars.train <- mtcars[-(8:14), ] # everything except for rows 8 .. 14
R <- lm(mpg ~ hp * wt, data=mtcars.train)
summary(R)
anova(R)

## prediction uses predict() method of regression model; values of predictor variables
## for new data points are passed in newdata= argument, as data frame with suitable names
## (mtcars.test has just the right form; columns hp and wt would be sufficient)
predicted.mpg <- predict(R, newdata=mtcars.test, interval="prediction")
## NB: interval="prediction" calculates confidence intervals for the predicted values

## print predictions for Mercedes models together with the correct values
print(cbind(predicted.mpg, mtcars.test))
## --> not quite spot on, but true values are contained in the (large) confidence intervals

## Exercise:
## Try to model other linear regressions (e.g. examples from pairs plots above) in a similar way.

##
## illustrate regression on two predictor variables with 3D plots
## [fuel consumption (mpg) against power (hp) and weight (wt) of cars in mtcars data set]
##

library(rgl) # ------------------------------ if you can't install RGL, skip this part
open3d() # RGL doesn't draw in the standard plot window

## plot3d() works almost like the standard plot() command, just in 3 dimensions
##  - use type="s" (small shaded spheres) if there aren't very many points in the scatterplot
##  - grab box with mouth to rotate, ctrl+click to zoom (or mouse wheel up/down)
plot3d(mtcars$hp, mtcars$wt, mtcars$mpg, # plot mpg against hp and wt
	type="s", size=.3, col="red",        # points are little (size=.2) red spheres
	xlim=c(0,400), ylim=c(0,8), zlim=c(0,40), # specify axis ranges, labels, etc.
	xlab="power", ylab="weight", zlab="miles/gallon")

## to plot regression plane as 3D surface, build regular grid G over x and y coordinates
grid.X <- seq(0, 400, length=21)  # X coordinates of grid points (horsepower)
grid.Y <- seq(0, 8, length=21)  # Y coordinates of grid points (weight)
## --> two-dimensional grid with 20 x 20 facets
##  - one facet would be enough for regression plane, but we'll plot curved surfaces later
##  - if your computer / graphics card isn't fast enough, reduce number of facets

## regression on power of engine (one predictor variable)
R <- lm(mpg ~ hp, data=mtcars)
summary(R) # --> R^2 = 60%

## z-coordinates for the surface3d() function are a bit tricky: we need a matrix that
## specifies z = f(x,y) for every combination of x and y values (from grid.X and grid.Y);
## the most convenient method uses outer(), which generates such a matrix from the 
## x- and y-coordinates of grid points and the "height function" f (such that z=f(x,y))

## f() must accept vectors x and y, and return the corresponding vector of height values z
f <- function (x,y) predict(R, newdata=data.frame(hp=x, wt=y))

## Z is the required matrix of height values ()
Z <- outer(grid.X, grid.Y, f)

## print Z to get a better intuitive grasp of the layout
rownames(Z) <- grid.X
colnames(Z) <- grid.Y
print(round(Z, 2))

## now plot the regression plane with surface3d()
par3d(ignoreExtent=TRUE) # it's important to disable automatic resizing of the 3D scene!
surface3d(grid.X, grid.Y, Z, col="blue", alpha=.7) # alpha < 1 makes surface semi-transparent
par3d(ignoreExtent=FALSE) # but turn it back on afterwards, otherwise next plot3d() will fail 

## compare this with regression on weight of the car (as single predictor)
pop3d() # remove previous regression plane from scene
R <- lm(mpg ~ wt, data=mtcars)
summary(R) # --> R^2 = 75%
Z <- outer(grid.X, grid.Y, f) # you know the drill ... 
par3d(ignoreExtent=TRUE)
surface3d(grid.X, grid.Y, Z, col="blue", alpha=.7)
par3d(ignoreExtent=FALSE)

## if we want to try more regression models, it's a good idea to wrap these lines in a function
lm.3d <- function (formula, data=mtcars) {
	R <- lm(formula, data=data)
	print(summary(R))
	## NB: global function f() can't see local variable R, so we need our own defintion of f!
	f <- function (x,y) predict(R, newdata=data.frame(hp=x, wt=y))
	Z <- outer(grid.X, grid.Y, f)
	par3d(ignoreExtent=TRUE)
	surface3d(grid.X, grid.Y, Z, col="blue", alpha=.7)
	par3d(ignoreExtent=FALSE)
	invisible(R) # return regression object, but make it invisible so it doesn't print automatically
}

## regression on both variables gives a much better fit
pop3d() 
lm.3d(mpg ~ hp + wt) # --> R^2 = 82.7%

## there is a significant interaction between power and weight -> curved surface
pop3d() 
lm.3d(mpg ~ hp * wt) # --> R^2 = 88.5%
## --> very good fit, but counterintuitive predictions outside data range! (why?)

## two-dimensional polynomial regression (need to list monomials explicitly)
pop3d()
R <- lm.3d(mpg ~ hp + wt + I(hp^2) + I(wt^2) + hp:wt, data=mtcars) # --> R^2 = 89.2%
## --> many factors are no longer significant ("causal competition")
## --> polynomial regression improves fit only slightly, but more plausible extrapolation 

## with a little fiddling, we can also plot a confidence region for the regression surface
## (NB: we are using interval="confidence" here to obtain a confidence region for the "true"
## regression surface; interval="prediction" as used above returns a larger confidence interval
## into which unseen data points are expected to fall with 95% confidence)
f.lower <- function (x,y) predict(R, newdata=data.frame(hp=x, wt=y), interval="confidence")[,2]
f.upper <- function (x,y) predict(R, newdata=data.frame(hp=x, wt=y), interval="confidence")[,3]
Z.lower <- outer(grid.X, grid.Y, f.lower)
Z.upper <- outer(grid.X, grid.Y, f.upper)
par3d(ignoreExtent=TRUE)
surface3d(grid.X, grid.Y, Z.lower, col="yellow", alpha=.3) # very light yellow
surface3d(grid.X, grid.Y, Z.upper, col="yellow", alpha=.3)
par3d(ignoreExtent=FALSE)
## --> good confidence along main diagonal (because power and weight are collinear)
## --> unreliable in orthogonal direction where there are hardly any data points

## Exercise:
## There is still one very implausible aspect of the regression surface.  What is it?
## Is polynomial regression an appropriate model for these data?  Is there a better transformation?

rgl.close()  # ------------------------------ end of RGL part
