library(SIGIL)
library(corpora)

LanguageCourse <- simulated.language.course()
sample.df(LanguageCourse, size=20, sort=TRUE)
summary(LanguageCourse)


## 1) Differences between classes?
boxplot(pre ~ class, data=LanguageCourse)
res <- aov(pre ~ class, data=LanguageCourse)
anova(res)

choose(7, 2) # number of pair-wise t-tests
TukeyHSD(res)

## 2) Is the course effective?
##    a) across all classes
with(LanguageCourse, t.test(post, pre, paired=TRUE))

##    b) for individual classes
classA <- subset(LanguageCourse, class=="A")
with(classA, t.test(post, pre, paired=TRUE))

##    c) explore differences in learning outcome btw classes
boxplot((post - pre) ~ class, data=LanguageCourse)
abline(h=0, col="blue")

res <- aov((post - pre) ~ class, data=LanguageCourse)
anova(res)

classD <- subset(LanguageCourse, class=="D")
with(classD, t.test(post - pre, mu=0)) # same as paired t-test for single class (cf. above)


## 3) Test normality assumptions etc.
with(LanguageCourse, { qqnorm(pre); qqline(pre) })
with(LanguageCourse, { qqnorm(post-pre); qqline(post-pre) })

with(LanguageCourse, boxplot(pre, post))
with(LanguageCourse, plot(pre, post, col=class, cex=1.5))
abline(0, 1)
