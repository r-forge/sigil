library(SIGIL)
library(corpora)

## results of a study testing the effectiveness of a computer-based language course
LC <- simulated.language.course()
sample.df(LC, size=20, sort=TRUE)
summary(LC) # what was the setup of the study?

## inappropriate t-test for independent samples
t.test(LC$post, LC$pre) # n.s.
boxplot(LC$pre, LC$post, names=c("pre", "post"), ylim=c(0,100), ylab="test score") # what is actually compared

## pre/post scores aren't independet samples because both depend on individual skills of students
plot(post ~ pre, data=LC, xlim=c(0,100), ylim=c(0,100)) # results in strong correlation
abline(0, 1, col="red") # suggests that most students actually improve

## appropriate test: t-test for paired samples
t.test(LC$post, LC$pre, data=LC, paired=TRUE) # ***, effect size >= 1.75 pts
t.test(LC$post - LC$pre) # simply a one-sample test of post-pre differences

## are there differences between the classes (in terms of language skills / level)?

## how many pairwise t-tests are carried out?

## what's the chance of at least one type I error if all H0 are true?
## --> family-wise error rate (FWER)

## Sidak correction, Bonferroni correction
## is Sidak too conservative or too liberal in this case?

## ANOVA + boxplot
## (perhaps mention var.test)

## post-hoc test for pair-wise diffs: TukeyHSD

## how do you test whether the course was successful in all classes?








## (ANOVA or multiple tests with Sidak correction)

## did you remember to verify the normality assumption?





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
