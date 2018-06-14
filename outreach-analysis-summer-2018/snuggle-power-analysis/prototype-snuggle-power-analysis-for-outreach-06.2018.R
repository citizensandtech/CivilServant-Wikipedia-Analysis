## LOAD LIBRARIES
library(ggplot2)
library(corrplot)
library(gmodels)
library(texreg)
library(stargazer)
library(rms)
library(blockrand)
library(lubridate)
library(psych)
library(MASS)

## CLEAR PROJECT
rm(list=ls())
set.seed(1528816215)

## METHODS
## randomSample pulls a random number of rows
## from a dataframe up to the number of rows
randomSample = function(df,n) { 
  return (df[sample(nrow(df), n,replace=TRUE),])
}

#######################################################################
## POWER ANALYSIS FOR OUTREACH JUNE 2018
## SNUGGLE EVALUATION
##
## ABOUT THIS PROJECT OVERALL
## This script conducts a power analysis in support of
## research by CivilServant with Wikimedia communities
## in 2018-2019. The purpose of this power analysis is
## to estimate in a naive way which communities are least
## likely to be able to support our research with Wikipedians
## so we do not unduly burden them with requests to work
## with us in our first pilot studies with Wikipedians
## More information is available at:
## https://meta.wikimedia.org/wiki/CivilServant%27s_Wikimedia_studies
##
## SNUGGLE EVALUATION POWER ANALYSIS 
## This document includes power analyses for a two-sided 
## experiment involving newcomers to Wikipedia
## and experienced contributors to Wikipedia.
## More information about this study is at:
## https://meta.wikimedia.org/wiki/CivilServant%27s_Wikimedia_studies#Snuggle


##########################################
## CONFIGURATION
## For now, start out with one language
## and then refactor to process multiple languages

language_wikipedias = c("es", "fr")

language = language_wikipedias[1]

#these thresholds have not been empirically validated
goodfaith.threshold = 0.65
damaging.threshold = 0.65


##########################################
##########################################
##           NEWCOMER DATAFRAME         ##
##########################################
##########################################

## This dataframe includes accounts that registered
## for a Wikipedia account in November 2017 for a given language
## along with information 
## Note that this power analysis is based on first-time editors who also
## registered their accounts in November 2017, so we undercount the number
## of first-time editors in this dataset
## FIELDS:
##   wiki: the language wiki in question
##   user.id: the user id of the account
##   registration.date: parseable date that the account registered 
##   registration: the date that the account registered (in Wikipedia timestamps)
##   edits.6.months: number of edits over a 6 month period
##   edits.2.weeks:  number of edits over 7 * 2 days
##   edits.4.weeks:  number of edits over 7 * 4 days
##   edits.8.weeks:  number of edits over 7 * 8 days
##   edits.12.weeks: number of edits over 7 * 12 days
##   ORES.rated.revisions: number of revisions that were rated by ORES
##                       revisions are excluded if ORES returns an error/NA
##                       or if the edit occurred after the first N edits (usually 5)
## IF AT LEAST ONE EDIT:
##   first.damaging.score: ORES damaging score, sampled late May 2018
##   first.goodfaith.score: ORES goodfaith score, sampled late May 2018
##   first.n.damaging.mean: ORES mean damaging score among first N edits, sampled late May 2018
##   first.n.goodfaith.mean: ORES mean goodfaith score among first N edits, sampled late May 2018
##   first.n.damaging.median: ORES median damaging score among first N edits, sampled late May 2018
##   first.n.goodfaith.median: ORES median goodfaith score among first N edits, sampled late May 2018
##
## SIMULATED VARIABLES FOR POWER ANALYSIS (SURVIVAL)
## These variables count the number of edits made in a particular interval,
## with the interval starting 48 hours after an account's first edit.
## We do this to imagine that someone receives a "treatment" for the experiment
## within the 48 hour period after making their first edit
## These simulated variables will be converted into binary values
## that indicate whether the account made 1 or more edits in the period
## and whether the account made 5 or more edits in the period
## Definitions via https://osf.io/preprints/socarxiv/8qsv6/
## FIELDS:
##   edits.3.4.weeks: number of edits made in the 1 week period starting on the 7 * 3 day
##   edits.4.8.weeks: number of edits made in the 4 week period starting on the 7 * 4 day
##   edits.8.24.weeks: number of edits made in the 16 week period starting on the 8 * 4 day
## SURVIVAL FIELDS: True or false, based on whether the account had made at least one edit
## in the period between the beginning of that week period and the end of the observation period
## (in this case 12 weeks + 48 hours after the first edit)
##   survival.week.period.N from 1 to 12

newcomers <- read.csv(paste("data/", language, "_newcomers_with_ores_scores_11.2017.csv", sep=""))
newcomers$registration.date <- as.POSIXct(newcomers$registration.date)
newcomers$registration.day <- as.numeric(floor(difftime(newcomers$registration.date,  min(newcomers$registration.date), unit="days")))

newcomers$edits.3.4.weeks.oneplus <- newcomers$edits.3.4.weeks > 0 
newcomers$edits.3.4.weeks.fiveplus <- newcomers$edits.3.4.weeks > 5

newcomers$edits.4.8.weeks.oneplus <- newcomers$edits.4.8.weeks > 0 
newcomers$edits.4.8.weeks.fiveplus <- newcomers$edits.4.8.weeks > 5

newcomers$edits.8.24.weeks.oneplus <- newcomers$edits.8.24.weeks > 0 
newcomers$edits.8.24.weeks.fiveplus <- newcomers$edits.8.24.weeks > 5

newcomers$goodfaith <- newcomers$first.goodfaith.score > goodfaith.threshold

summary(newcomers$edit.count > 0)
newcomer.editors <- subset(newcomers, edit.count > 0)

##########################################
## UNIVARIATE SUMMARY STATISTICS

## Edits over time
hist(newcomer.editors$edits.6.months)
hist(log1p(newcomer.editors$edits.6.months))
summary(newcomer.editors$edits.6.months>1)
summary(newcomer.editors$edits.12.weeks>1)
summary(newcomer.editors$edits.8.weeks>1)
summary(newcomer.editors$edits.4.weeks>1)
summary(newcomer.editors$edits.2.weeks>1)

## Wikipedian Survival

summary(newcomer.editors$edits.3.4.weeks.oneplus)
summary(newcomer.editors$edits.3.4.weeks.fiveplus)

summary(newcomer.editors$edits.4.8.weeks.oneplus)
summary(newcomer.editors$edits.4.8.weeks.fiveplus)

summary(newcomer.editors$edits.8.24.weeks.oneplus)
summary(newcomer.editors$edits.8.24.weeks.fiveplus)

## Days since earlest registration in dataset
summary(newcomer.editors$registration.day)
hist(newcomer.editors$registration.day)

## Revisions rated by ORES
summary(factor(newcomer.editors$ORES.rated.revisions))
hist(newcomer.editors$ORES.rated.revisions)

## ORES first damaging score
summary(newcomer.editors$first.damaging.score)
hist(newcomer.editors$first.damaging.score)
summary(newcomer.editors$first.damaging.score>damaging.threshold)

## ORES first.n.damaging.mean
summary(newcomer.editors$first.n.damaging.mean)
hist(newcomer.editors$first.n.damaging.mean)
summary(newcomer.editors$first.n.damaging.mean>damaging.threshold)

## ORES first.n.damaging.median
summary(newcomer.editors$first.n.damaging.median)
hist(newcomer.editors$first.n.damaging.median)
summary(newcomer.editors$first.n.damaging.median>damaging.threshold)

## ORES first goodfaith score
summary(newcomer.editors$first.goodfaith.score)
hist(newcomer.editors$first.goodfaith.score)
summary(newcomer.editors$first.goodfaith.score>damaging.threshold)

## ORES first goodfaith mean
summary(newcomer.editors$first.n.goodfaith.mean)
hist(newcomer.editors$first.n.goodfaith.mean)
summary(newcomer.editors$first.n.goodfaith.mean>goodfaith.threshold)

## ORES first goodfaith median
summary(newcomer.editors$first.n.goodfaith.median)
hist(newcomer.editors$first.n.goodfaith.median)
summary(newcomer.editors$first.n.goodfaith.median>goodfaith.threshold)

### survival PROBABILITY OF DROPPING OUT
### LOAD NEWCOMER SURVIVAL WEEK DATASET
### This dataset is an expanded version of the newcomer.editors dataframe
### FIELDS
###   user.id: Wikipedia user id
###   first.damaging score: first edit damaging score from ORES sampled in May 2018
###   first.goodfaith.score: first goodfaith score from ORES sampled in May 2018
###   week: which week period starting 48 hours after the first edit
###   survived: made at least one edit between this point and a point 20 weeks later
newcomer.survival.weeks <- read.csv(paste("data/", language, "_newcomer_survival_week_periods.2017.csv", sep=""))

print(paste("Are the length of newcomer editors and survival week editors equal? ", 
            (length(unique(newcomer.survival.weeks$user.id)) == nrow(newcomer.editors)), sep="1"))

newcomer.survival.weeks$goodfaith <- newcomer.survival.weeks$first.goodfaith.score > goodfaith.threshold 
newcomer.survival.weeks$damaging.goodfaith <-
  (newcomer.survival.weeks$first.goodfaith.score > goodfaith.threshold) &
  (newcomer.survival.weeks$first.damaging.score > damaging.threshold)

## ESTIMATE THE CHANCE OF SURVIVAL TO A GIVEN WEEK AND PLOT
summary(nsw1 <- glm(survived ~ factor(week) + damaging.goodfaith, data=newcomer.survival.weeks, family="binomial"))
nsw1.estimates.gf <- data.frame(week=c(seq(1,20)), damaging.goodfaith=TRUE)
nsw1.fit.gf <- predict(nsw1, nsw1.estimates.gf, type="link", se.fit=TRUE)
nsw1.estimates.gf$estimated.prob <- family(nsw1)$linkinv(nsw1.fit.gf$fit)
nsw1.estimates.gf$estimated.upper <- family(nsw1)$linkinv( nsw1.fit.gf$fit + 1.96*nsw1.fit.gf$se.fit)
nsw1.estimates.gf$estimated.lower <- family(nsw1)$linkinv( nsw1.fit.gf$fit - 1.96*nsw1.fit.gf$se.fit)

nsw1.estimates.nogf <- data.frame(week=c(seq(1,20)), damaging.goodfaith=FALSE)
nsw1.fit.nogf <- predict(nsw1, nsw1.estimates.nogf, type="link", se.fit=TRUE)
nsw1.estimates.nogf$estimated.prob <- family(nsw1)$linkinv(nsw1.fit.nogf$fit)
nsw1.estimates.nogf$estimated.upper <- family(nsw1)$linkinv( nsw1.fit.nogf$fit + 1.96*nsw1.fit.nogf$se.fit)
nsw1.estimates.nogf$estimated.lower <- family(nsw1)$linkinv( nsw1.fit.nogf$fit - 1.96*nsw1.fit.nogf$se.fit)

nsw1.estimates <- rbind(nsw1.estimates.gf, nsw1.estimates.nogf)

ggplot(nsw1.estimates, aes(week, estimated.prob, color=damaging.goodfaith)) +
  geom_line() +
  geom_ribbon(data=nsw1.estimates, 
              mapping=aes(x=week, ymax=estimated.upper, ymin=estimated.lower, 
                          fill=damaging.goodfaith, alpha=0.2)) +
  theme_bw(base_size = 15, base_family = "Helvetica") +
  scale_x_continuous(breaks=seq(1,20,1)) +
  xlab("Week period starting 48 hours after 1st edit") +
  ylab("Estimated prob of continuing to edit Wikipedia") +
  ggtitle(paste("Estimated newcomer survival curve by week for [",language, "] Wikipedia", sep=""))

## OBSERVED SURVIVAL PROBABILITIES
survival.prob <- function(col){
  summary(col)[['True']]/length(col)
}

editor.subset.survival <- function(df){
  nes <- data.frame(
    week = c(1,2,3,4,5,6,7,8,9,10,11,12, 13, 14, 15, 16, 17, 18 , 19, 20),
    survival = c(survival.prob(df$survival.week.period.1),
                 survival.prob(df$survival.week.period.2),
                 survival.prob(df$survival.week.period.3),
                 survival.prob(df$survival.week.period.4),
                 survival.prob(df$survival.week.period.5),
                 survival.prob(df$survival.week.period.6),
                 survival.prob(df$survival.week.period.7),
                 survival.prob(df$survival.week.period.8),
                 survival.prob(df$survival.week.period.9),
                 survival.prob(df$survival.week.period.10),
                 survival.prob(df$survival.week.period.11),
                 survival.prob(df$survival.week.period.12),
                 survival.prob(df$survival.week.period.13),
                 survival.prob(df$survival.week.period.14),
                 survival.prob(df$survival.week.period.15),
                 survival.prob(df$survival.week.period.16),
                 survival.prob(df$survival.week.period.17),
                 survival.prob(df$survival.week.period.18),
                 survival.prob(df$survival.week.period.19),
                 survival.prob(df$survival.week.period.20)
    )
  )
}

newcomer.editors.survival <- editor.subset.survival(newcomer.editors)

summary(lm(survival ~ log(week), data=newcomer.editors.survival))

ggplot(newcomer.editors.survival, aes(as.numeric(week), survival)) +
  geom_smooth(method="lm", formula=y ~ log(x)) +
  geom_point() +
  ylim(0,0.2) +
  scale_x_continuous(breaks=seq(1,20,1)) +
  theme_bw(base_size = 15, base_family = "Helvetica") +
  xlab("Week period starting 48 hours after 1st edit") +
  ylab("Prob of continuing to edit Wikipedia") +
  ggtitle(paste("Observed survival curve by week for [",language, "] Wikipedia", sep=""))
  
##########################################
## BIVARIATE SUMMARY STATISTICS

## edit count and first damaging score
describeBy(subset(newcomer.editors, edit.count<=20)$first.damaging.score, 
           subset(newcomer.editors, edit.count<=20)$edit.count)
# ggplot(subset(newcomer.editors, edit.count < 100), aes(factor(edit.count), first.damaging.score)) +
#   geom_boxplot()

## edit count and mean damaging score out of n
describeBy(subset(newcomer.editors, edit.count<=20)$first.n.damaging.mean, 
           subset(newcomer.editors, edit.count<=20)$edit.count)
# ggplot(subset(newcomer.editors, edit.count < 100), aes(factor(edit.count), first.n.damaging.mean)) +
#   geom_boxplot()


## edit count and first edit goodfaith score
describeBy(subset(newcomer.editors, edit.count<=20)$first.goodfaith.score, 
           subset(newcomer.editors, edit.count<=20)$edit.count)
# ggplot(subset(newcomer.editors, edit.count < 100), aes(factor(edit.count), first.goodfaith.score)) +
#   geom_boxplot()

## edit count and mean goodfaith score out of n
describeBy(subset(newcomer.editors, edit.count<=20)$first.n.goodfaith.mean, 
           subset(newcomer.editors, edit.count<=20)$edit.count)
# ggplot(subset(newcomer.editors, edit.count < 100), aes(factor(edit.count), first.n.goodfaith.mean)) +
#   geom_boxplot()

CrossTable(newcomer.editors$first.goodfaith.score>goodfaith.threshold, 
           newcomer.editors$first.damaging.score>damaging.threshold, 
           digits=3, max.width = 5, expected=FALSE, prop.r=TRUE, prop.c=TRUE, prop.t=FALSE, prop.chisq=FALSE) 

## Survival measures and goodfaith score
### three to four weeks
describeBy(newcomer.editors$first.damaging.score,
           newcomer.editors$edits.3.4.weeks.oneplus)
# ggplot(newcomer.editors, aes(edits.3.4.weeks.oneplus, first.damaging.score)) +
#   geom_boxplot()
# ggplot(newcomer.editors, aes(edits.3.4.weeks.oneplus, first.n.damaging.mean)) +
#   geom_boxplot()

describeBy(newcomer.editors$first.damaging.score,
           newcomer.editors$edits.3.4.weeks.fiveplus)
# ggplot(newcomer.editors, aes(edits.3.4.weeks.fiveplus, first.damaging.score)) +
#   geom_boxplot()
# ggplot(newcomer.editors, aes(edits.3.4.weeks.fiveplus, first.n.damaging.mean)) +
#   geom_boxplot()

### four to eight weeks
describeBy(newcomer.editors$first.goodfaith.score,
           newcomer.editors$edits.4.8.weeks.oneplus)
# ggplot(newcomer.editors, aes(edits.4.8.weeks.oneplus, first.goodfaith.score)) +
#   geom_boxplot()
# ggplot(newcomer.editors, aes(edits.4.8.weeks.oneplus, first.n.goodfaith.mean)) +
#   geom_boxplot()

describeBy(newcomer.editors$first.goodfaith.score,
           newcomer.editors$edits.4.8.weeks.fiveplus)
# ggplot(newcomer.editors, aes(edits.4.8.weeks.fiveplus, first.goodfaith.score)) +
#   geom_boxplot()
# ggplot(newcomer.editors, aes(edits.4.8.weeks.fiveplus, first.n.goodfaith.mean)) +
#   geom_boxplot()

### eight to twenty-four weeks
describeBy(newcomer.editors$first.goodfaith.score,
           newcomer.editors$edits.8.24.weeks.oneplus)
# ggplot(newcomer.editors, aes(edits.8.24.weeks.oneplus, first.goodfaith.score)) +
#   geom_boxplot()
# ggplot(newcomer.editors, aes(edits.8.24.weeks.oneplus, first.n.goodfaith.mean)) +
#   geom_boxplot()

describeBy(newcomer.editors$first.goodfaith.score,
           newcomer.editors$edits.8.24.weeks.oneplus)
# ggplot(newcomer.editors, aes(edits.8.24.weeks.oneplus, first.goodfaith.score)) +
#   geom_boxplot()
# ggplot(newcomer.editors, aes(edits.8.24.weeks.oneplus, first.n.goodfaith.mean)) +
#   geom_boxplot()

## Survival measures and damaging score
### three to four weeks
describeBy(newcomer.editors$first.damaging.score,
           newcomer.editors$edits.3.4.weeks.oneplus)
# ggplot(newcomer.editors, aes(edits.3.4.weeks.oneplus, first.damaging.score)) +
#   geom_boxplot()
# ggplot(newcomer.editors, aes(edits.3.4.weeks.oneplus, first.n.damaging.mean)) +
#   geom_boxplot()

describeBy(newcomer.editors$first.damaging.score,
           newcomer.editors$edits.3.4.weeks.fiveplus)
# ggplot(newcomer.editors, aes(edits.3.4.weeks.fiveplus, first.damaging.score)) +
#   geom_boxplot()
# ggplot(newcomer.editors, aes(edits.3.4.weeks.fiveplus, first.n.damaging.mean)) +
#   geom_boxplot()

### four to eight weeks
describeBy(newcomer.editors$first.damaging.score,
           newcomer.editors$edits.4.8.weeks.oneplus)
# ggplot(newcomer.editors, aes(edits.4.8.weeks.oneplus, first.damaging.score)) +
#   geom_boxplot()
# ggplot(newcomer.editors, aes(edits.4.8.weeks.oneplus, first.n.damaging.mean)) +
#   geom_boxplot()

describeBy(newcomer.editors$first.damaging.score,
           newcomer.editors$edits.4.8.weeks.fiveplus)
# ggplot(newcomer.editors, aes(edits.4.8.weeks.fiveplus, first.damaging.score)) +
#   geom_boxplot()
# ggplot(newcomer.editors, aes(edits.4.8.weeks.fiveplus, first.n.damaging.mean)) +
#   geom_boxplot()

### eight to twenty-four weeks
describeBy(newcomer.editors$first.damaging.score,
           newcomer.editors$edits.8.24.weeks.oneplus)
# ggplot(newcomer.editors, aes(edits.8.24.weeks.oneplus, first.damaging.score)) +
#   geom_boxplot()
# ggplot(newcomer.editors, aes(edits.8.24.weeks.oneplus, first.n.damaging.mean)) +
#   geom_boxplot()

describeBy(newcomer.editors$first.damaging.score,
           newcomer.editors$edits.8.24.weeks.oneplus)
# ggplot(newcomer.editors, aes(edits.8.24.weeks.oneplus, first.damaging.score)) +
#   geom_boxplot()
# ggplot(newcomer.editors, aes(edits.8.24.weeks.oneplus, first.n.damaging.mean)) +
#   geom_boxplot()

##########################################
## POWER ANALYSIS FOR NEWCOMERS         ##
##########################################
## In the following power analysis, we consider the folliwing two experiments: 
## a) among a 20% sample of newcomers, offering the teahouse versus an evaluation of removal
## b) among newcomers with high damaging, high goodfaith scores who are not enrolled in (a)
##    offer support, versus evaluation of removal
## in both cases, we look at the effect on edit count over time.
## In this power analysis, we assume (not accurately) that this data represents
## a state (a) without mentorship and (b) without an evaluation of edits for removal

## set the minimum observable effect to 1 percentage point increase
## this is what Morgan & Halfaker found from the Teahouse
## and it's likely to be an over-estimate since only 2 out of
## six of their comparisons were statistically-significant
min.observed.effect.no.ml = 1.1
min.observed.effect.no.ml.count = 1.1
ml.recruitment.sample = 0.7
number.of.comparisons = 2 ## 2 comparisons per group (omitting Morgan & Halfaker)
alpha.level = 0.05 / number.of.comparisons
power.goal = 0.8
num.simulations = 50

# I hypothesize that targeted interventions
# toward good faith contributors will have a 
# greater magnitude effect than a random attempt
# which is what the Teahouse intervention did
# so I am setting the min observable effect at 5 pct points
min.observed.effect.dg.ml = 1.5
min.observed.effect.dg.ml.count = 1.2
min.observed.effect.g.ml = 1.5
min.observed.effect.g.ml.count = 1.2
num.newcomer.editors <- nrow(newcomer.editors)

###### POWER CALCULATION FOR BINARY OUTCOMES
## Power calculation method by Alexander Coppock, licensed under the GPL 3.0
## https://egap.shinyapps.io/Power_Calculator/

power.calculator.binary <- function(p1, p0, alpha=0.05, N){
  lowertail <- (abs(p1 - p0) * sqrt(N/2))/sqrt(p1*(1-p1) + p0*(1-p0))
  uppertail <- -1*lowertail
  beta <- pnorm(lowertail- qnorm(1-alpha/2), lower.tail=TRUE) + 
    1 - pnorm(uppertail- qnorm(1-alpha/2), lower.tail=FALSE)
  return(beta)
}

## Iterate through power calculation to find minimum sample
min.sample.for.power.binary <- function(prop.1, prop.2, min.sample.size, max.sample.size, iterate.by, alpha.level, power.goal){
  pa.results <- data.frame(
    prop.1 = prop.1,
    prop.2 = prop.2,
    sample.size = min.sample.size,
    alpha = alpha.level,
    power=power.calculator.binary(prop.1, prop.2, N=min.sample.size))

  for(i in seq(1,(max.sample.size - min.sample.size)/iterate.by)){
    current.sample.size <- min.sample.size + i * iterate.by
    pa.results <- rbind(pa.results,
                        data.frame(
                          prop.1 = prop.1,
                          prop.2 = prop.2,
                          sample.size = current.sample.size,
                          alpha = alpha.level,
                          power=power.calculator.binary(prop.1, prop.2, N=current.sample.size))
                        )
  }
  pa.results
  min(subset(pa.results, power >= power.goal)$sample.size)
}

#### SET UP SUBSET DATAFRAMES
damaging.goodfaith.editors  <- 
  subset(newcomer.editors, 
         (first.damaging.score > damaging.threshold) & 
         (first.goodfaith.score > goodfaith.threshold))

goodfaith.editors  <- 
  subset(newcomer.editors, 
           (first.goodfaith.score > goodfaith.threshold))

damaging.editors  <- 
  subset(newcomer.editors, 
         (first.damaging.score > damaging.threshold))

num.damaging.goodfaith.editors <- nrow(damaging.goodfaith.editors)
num.damaging.editors <- nrow(damaging.editors)
num.goodfaith.editors <- nrow(goodfaith.editors)

### SURVIVAL POWER ANALYSIS
## CREATE SIMULATED DATASET OF SURVIVAL ANALYSIS
damaging.goodfaith.survival.weeks <- subset(newcomer.survival.weeks, damaging.goodfaith==TRUE)
damaging.goodfaith.survival <- editor.subset.survival(damaging.goodfaith.editors)


## CALCULATE MINIMUM SAMPLES PER WEEK FOR SURVIVAL ANALYSIS
## INVOLVING A PER-WEEK LOGISTIC REGRESSION

survival.min.samples <- function(df, sample.size, treat.multiplier, alpha.level, power.goal){
  min.samples = c()
  for(i in seq(1,nrow(df))){
    current.survival <- df$survival[df$week==i]
    min.sample <- min.sample.for.power.binary(current.survival,
                                       current.survival * treat.multiplier,
                                       sample.size  / 4,
                                       sample.size * 40,
                                       sample.size / 4,
                                       alpha.level,
                                       power.goal) / sample.size
    min.samples <- c(min.samples, min.sample)
  }
  min.samples
}

with.ml.g.sample <- num.goodfaith.editors * ml.recruitment.sample
with.ml.dg.sample <- num.damaging.goodfaith.editors * ml.recruitment.sample
all.sample <- num.newcomer.editors * (1 - ml.recruitment.sample)
without.ml.sample <- all.sample

all.survival <- editor.subset.survival(newcomer.editors)
ge.survival <- editor.subset.survival(goodfaith.editors)
dge.survival <- editor.subset.survival(damaging.goodfaith.editors)

ge.survival$min.sample <- survival.min.samples(ge.survival, with.ml.g.sample, 
                                               (min.observed.effect.g.ml), 
                                               alpha.level, power.goal)
dge.survival$min.sample <- survival.min.samples(dge.survival, with.ml.dg.sample, 
                                               (min.observed.effect.dg.ml), 
                                               alpha.level, power.goal)
ge.survival.12.week.survival.min.sample <- subset(ge.survival, week==12)$min.sample
dge.survival.12.week.survival.min.sample <- subset(dge.survival, week==12)$min.sample


all.survival$min.sample <- survival.min.samples(all.survival, all.sample, 
                                               (min.observed.effect.no.ml), 
                                               alpha.level, power.goal)
all.survival.12.week.survival.min.sample <- subset(all.survival, week==12)$min.sample




#######################################
### POWER ANALYSIS FOR COUNT OUTCOMES #
#######################################

simulate.study.from.negbin <- function(i, sample.size, model.nb, pct.diff, alpha.level){
  cat(".")
  coef = log(pct.diff) ## for example for a 12% increase, log(1.12) = 0.1133287
  control.group <- data.frame(TREAT=0, dv = rnegbin(sample.size/2, mu = model.nb$coefficients[1][['(Intercept)']], theta = model.nb$theta))
  treat.group <- data.frame(TREAT=1, dv = rnegbin(sample.size/2, mu = model.nb$coefficients[1][['(Intercept)']] + coef, theta = model.nb$theta))
  sim.obs <- rbind(control.group, treat.group)
  
  m.nb <- glm.nb(dv ~ TREAT, data=sim.obs)
  m.nb.intercept <- m.nb$coefficients[['(Intercept)']]
  m.nb.treat.effect <- m.nb$coefficients[['TREAT']]
  m.nb.treat.coef = coef(summary(m.nb))[2,]
  m.nb.stderr <- m.nb.treat.coef['Std. Error']
  m.nb.pvalue <- m.nb.treat.coef['Pr(>|z|)']
  m.nb.significant <- as.double(m.nb.pvalue) < alpha.level 
  data.frame(i=i,
             power.sim.comments.intercept      = m.nb.intercept, 
             power.sim.comments.treat.effect   = m.nb.treat.effect,
             power.sim.comments.stderr         = m.nb.stderr,
             power.sim.comments.pvalue         = m.nb.pvalue,
             power.sim.comments.significant    = m.nb.significant)
}


power.analysis.nb <- function(i, sample.size, model.nb, pct.diff, num.models, alpha.level){
  pm <- simulate.study.from.negbin(1,sample.size,model.nb, pct.diff, alpha.level)
  for(i in seq(2,num.models)){
    pm <- rbind(pm, simulate.study.from.negbin(i,sample.size,model.nb, pct.diff, alpha.level))
  }
  
  pct.significant <- sum(pm$power.sim.comments.significant)/nrow(pm)*100
  mean.effect <- mean(pm$power.sim.comments.treat.effect)
  data.frame(i=i, sample.size=sample.size, 
             pct.significant = pct.significant, 
             mean.effect = mean.effect)
}


nb.simulate.power <- function(df, subset.prop, sim.effect, num.models, alpha.level){
  maximum.sample.size <- nrow(df) * subset.prop * 10 ## maximum 6 months
  minimum.sample.size <- nrow(df) * subset.prop / 4 
  iterate.by <- nrow(df) * subset.prop / 4 ## iterate by one week periods
  summary(base.num.edits.nb  <- glm.nb(edits.8.weeks ~ 1, data=df))

  cat(paste("\nsample.size",minimum.sample.size))
  models.decision.nb <- power.analysis.nb(1,minimum.sample.size, base.num.edits.nb, sim.effect, num.models, alpha.level)
  for(i in seq(1,maximum.sample.size/iterate.by)){
    cat(paste("\nsample.size",minimum.sample.size + iterate.by * i))
    models.decision.nb <- rbind(models.decision.nb, power.analysis.nb(i,
                                                                      minimum.sample.size + iterate.by*i, base.num.edits.nb,
                                                                      sim.effect, num.models, alpha.level))
  }
  models.decision.nb
}



## RUN THE POWER ANALYSES FOR newcomer, damaging+goodfaith, goodfaith samples

all.count.power.analys.nb <- nb.simulate.power(newcomer.editors, 1 - ml.recruitment.sample, 
                                               min.observed.effect.no.ml.count, num.simulations, alpha.level) 

dg.count.power.analys.nb <- nb.simulate.power(damaging.goodfaith.editors, ml.recruitment.sample, 
                                              min.observed.effect.dg.ml.count, num.simulations, alpha.level) 

g.count.power.analys.nb <- nb.simulate.power(goodfaith.editors, ml.recruitment.sample, 
                                             min.observed.effect.g.ml.count, num.simulations, alpha.level) 



ggplot(all.count.power.analys.nb, aes(sample.size/nrow(newcomer.editors), pct.significant)) +
  geom_smooth() +
  geom_point ()

ggplot(dg.count.power.analys.nb, aes(sample.size/nrow(damaging.goodfaith.editors), pct.significant)) +
  geom_smooth() +
  geom_point ()

ggplot(g.count.power.analys.nb, aes(sample.size/nrow(goodfaith.editors), pct.significant)) +
  geom_smooth() +
  geom_point ()

all.count.power.analysis.nb.min.sample <- min(subset(all.count.power.analys.nb, pct.significant>80)$sample.size)
dg.count.power.analys.nb.min.sample <- min(subset(dg.count.power.analys.nb, pct.significant>80)$sample.size)
g.count.power.analys.nb.min.sample <- min(subset(g.count.power.analys.nb, pct.significant>80)$sample.size)

all.count.power.analysis.nb.min.sample.months <- all.count.power.analysis.nb.min.sample / nrow(newcomer.editors)
dg.count.power.analys.nb.min.sample.months <- dg.count.power.analys.nb.min.sample / nrow(damaging.goodfaith.editors)
g.count.power.analys.nb.min.sample.months <- g.count.power.analys.nb.min.sample / nrow(goodfaith.editors)

##################################
### REPORT POWER ANALYSIS RESULTS
##################################
cat(paste("\n",
          "####################################################################\n",
          "##      SNUGGLE STUDY LANGUAGE WIKIPEDIA POWER ANALYSIS REPORT   ##\n",
          "####################################################################\n",
          "\n",
          "Date: ", Sys.Date(),"\n",
          "Language: ", language,  "\n",
          num.newcomer.editors, " newcomer editors registered per month", "\n",
          num.damaging.goodfaith.editors, " damaging + goodfaith first-time editors per month: ", "\n",
          ml.recruitment.sample, " of newcomer editors retained for ML test", "\n",
          min.observed.effect.no.ml, " simulated min observed effect on survival, no-ML test", "\n",
          min.observed.effect.g.ml, " simulated min observed effect on survival, ML test", "\n",
          min.observed.effect.no.ml.count, " simulated min observed effect on edit count, no-ML test", "\n",
          min.observed.effect.g.ml.count, " simulated min observed effect on edit count, ML test (goodfaith)", "\n",
          min.observed.effect.dg.ml.count, " simulated min observed effect on edit count, ML test (damaging goodfaith)", "\n",
          power.goal , " goal for chance of observing effect", "\n",
          number.of.comparisons , " total comparisons", "\n", ## should this be per-subgroup?
          alpha.level , " bonferroni adjusted alpha level (0.05 / num comparisons)", "\n",
          "\n",
          "####################################################################\n",
          "##      RESULTS                                                   ##\n",
          "####################################################################\n",
          "Non-ML Test of Teahouse vs Revision Review (3 month survival):\n",
          "  ", all.survival.12.week.survival.min.sample, " months, ", 
                all.survival.12.week.survival.min.sample*without.ml.sample,
          " participants (for test within across editors retained for no ML intervention)", "\n",
          "Non-ML Test of Teahouse vs Revision Review (edit count over 3 months):\n",
          "  ", all.count.power.analysis.nb.min.sample.months, " months, ", 
                all.count.power.analysis.nb.min.sample,
          " participants (for test within across editors retained for no ML intervention)", "\n\n",
          "Test of Snuggle vs Standard Revision Review (3 month survival):\n",
          "  ", ge.survival.12.week.survival.min.sample, " months, ", 
                ge.survival.12.week.survival.min.sample*without.ml.sample,
                " participants (for test within retained goodfaith)", "\n",
          "  ", dge.survival.12.week.survival.min.sample, " months, ", 
                dge.survival.12.week.survival.min.sample*without.ml.sample,
                " participants (for test within retained damaging/goodfaith)", "\n",
          "Test of Snuggle vs Standard Revision Review (edit count over 3 months):\n",
          "  ", g.count.power.analys.nb.min.sample.months, " months, ", 
                g.count.power.analys.nb.min.sample,
          " participants (for test within retained goodfaith)", "\n",
          "  ", dg.count.power.analys.nb.min.sample.months, " months, ", 
                dg.count.power.analys.nb.min.sample,
          " participants (for test within retained damaging/goodfaith)", "\n",
          #         , "", "\n",
          sep=""))

##########################################
##########################################
##     NEWCOMER REVISION DATAFRAME      ##
##########################################
##########################################
## This dataframe includes revisions over 6 months for Wikipedia accounts
## that were registered in November 2017, for a given language:
##   wiki: language wikipedia
##   user.id: the user id of the account
##   registration: the date that the account registered
##   edits.6.months: number of edits over a 6 month long period
##   damaging: ORES score for that revision (0 to 1)
##   goodfaith: ORES score for goodfaith (0 to 1)
##   revision.id: ID in Wikipedia for the revision
##   revision.time: string containing the revision timestamp YYYY-MM-D HH:MM:SS
##   registration: string containing the account registration timestamp YYYY-MM-D HH:MM:SS
#newcomer.revisions <- read.csv(paste("data/", language, "_revisions_with_user_11.2017.csv", sep=""))



####################################################################
### SCRATCH AREA                  ##################################
####################################################################


# 
# ##########################################
# ## ESTIMATING PARTICIPATION VOLUME OVER TIME
# 
# summary(m.4.weeks.d <- glm.nb(edits.4.weeks ~ first.damaging.score, data=newcomer.editors))
# summary(m.8.weeks.d <- glm.nb(edits.8.weeks ~ first.damaging.score, data=newcomer.editors))
# 
# summary(m.4.weeks.g <- glm.nb(edits.4.weeks ~ first.goodfaith.score, data=newcomer.editors))
# summary(m.8.weeks.g <- glm.nb(edits.8.weeks ~ first.goodfaith.score, data=newcomer.editors))
# 
# summary(m.4.weeks.g.d <- glm.nb(edits.4.weeks ~ first.goodfaith.score + first.damaging.score, data=newcomer.editors))
# summary(m.8.weeks.g.d <- glm.nb(edits.8.weeks ~ first.goodfaith.score + first.damaging.score, data=newcomer.editors))
# 
# stargazer(m.4.weeks.d, m.4.weeks.g, m.4.weeks.g.d, type="text")
# 
# stargazer(m.8.weeks.d, m.8.weeks.g, m.8.weeks.g.d, type="text")
# 
# ##########################################
# ## ESTIMATING "SURVIVAL" (MORGAN AND HALFAKER)
# ### First edit goodfaith score
# summary(s.3.4.weeks.one.g <- glm(edits.3.4.weeks.oneplus ~ first.goodfaith.score, data=newcomer.editors))
# summary(s.4.8.weeks.one.g <- glm(edits.4.8.weeks.oneplus ~ first.goodfaith.score, data=newcomer.editors))
# summary(s.8.24.weeks.one.g <- glm(edits.8.24.weeks.oneplus ~ first.goodfaith.score, data=newcomer.editors))
# 
# summary(s.3.4.weeks.five.g <- glm(edits.3.4.weeks.fiveplus ~ first.goodfaith.score, data=newcomer.editors))
# summary(s.4.8.weeks.five.g <- glm(edits.4.8.weeks.fiveplus ~ first.goodfaith.score, data=newcomer.editors))
# summary(s.8.24.weeks.five.g <- glm(edits.8.24.weeks.fiveplus ~ first.goodfaith.score, data=newcomer.editors))
# 
# ### First edit damaging score
# summary(s.3.4.weeks.one.d <- glm(edits.3.4.weeks.oneplus ~ first.damaging.score, data=newcomer.editors))
# summary(s.4.8.weeks.one.d <- glm(edits.4.8.weeks.oneplus ~ first.damaging.score, data=newcomer.editors))
# summary(s.8.24.weeks.one.d <- glm(edits.8.24.weeks.oneplus ~ first.damaging.score, data=newcomer.editors))
# 
# summary(s.3.4.weeks.five.d <- glm(edits.3.4.weeks.fiveplus ~ first.damaging.score, data=newcomer.editors))
# summary(s.4.8.weeks.five.d <- glm(edits.4.8.weeks.fiveplus ~ first.damaging.score, data=newcomer.editors))
# summary(s.8.24.weeks.five.d <- glm(edits.8.24.weeks.fiveplus ~ first.damaging.score, data=newcomer.editors))
# 
# ### First edit goodfaith score + damaging score
# summary(s.3.4.weeks.one.g.d <- glm(edits.3.4.weeks.oneplus ~ first.goodfaith.score + first.damaging.score, data=newcomer.editors))
# summary(s.4.8.weeks.one.g.d <- glm(edits.4.8.weeks.oneplus ~ first.goodfaith.score + first.damaging.score, data=newcomer.editors))
# summary(s.8.24.weeks.one.g.d <- glm(edits.8.24.weeks.oneplus ~ first.goodfaith.score + first.damaging.score, data=newcomer.editors))
# 
# summary(s.3.4.weeks.five.g.d <- glm(edits.3.4.weeks.fiveplus ~ first.goodfaith.score + first.damaging.score, data=newcomer.editors))
# summary(s.4.8.weeks.five.g.d <- glm(edits.4.8.weeks.fiveplus ~ first.goodfaith.score + first.damaging.score, data=newcomer.editors))
# summary(s.8.24.weeks.five.g.d <- glm(edits.8.24.weeks.fiveplus ~ first.goodfaith.score + first.damaging.score, data=newcomer.editors))
# 
# ## RESULTS TABLES (MORGAN & HALFAKER)
# ## (using stargzer, which outputs DV)
# stargazer(s.3.4.weeks.one.g,s.4.8.weeks.one.g, s.8.24.weeks.one.g, type="text")
# stargazer(s.3.4.weeks.five.g,s.4.8.weeks.five.g, s.8.24.weeks.five.g, type="text")
# 
# stargazer(s.3.4.weeks.one.d,s.4.8.weeks.one.d, s.8.24.weeks.one.d, type="text")
# stargazer(s.3.4.weeks.five.d,s.4.8.weeks.five.d, s.8.24.weeks.five.d, type="text")
# 
# stargazer(s.3.4.weeks.one.g.d,s.4.8.weeks.one.g.d, s.8.24.weeks.one.g.d, type="text")
# stargazer(s.3.4.weeks.five.g.d,s.4.8.weeks.five.g.d, s.8.24.weeks.five.g.d, type="text")
# 
# 
# 
# ############ MORGAN AND HALFAKER POWER ANALYSIS
# dge.3.4.weeks.oneplus.prop <- nrow(subset(damaging.goodfaith.editors, edits.3.4.weeks.oneplus==TRUE)) / num.damaging.goodfaith.editors
# dge.3.4.weeks.fiveplus.prop <- nrow(subset(damaging.goodfaith.editors, edits.3.4.weeks.fiveplus==TRUE)) / num.damaging.goodfaith.editors
# dge.4.8.weeks.oneplus.prop <- nrow(subset(damaging.goodfaith.editors, edits.4.8.weeks.oneplus==TRUE)) / num.damaging.goodfaith.editors
# dge.4.8.weeks.fiveplus.prop <- nrow(subset(damaging.goodfaith.editors, edits.4.8.weeks.fiveplus==TRUE)) / num.damaging.goodfaith.editors
# dge.8.24.weeks.oneplus.prop <- nrow(subset(damaging.goodfaith.editors, edits.8.24.weeks.oneplus==TRUE)) / num.damaging.goodfaith.editors
# dge.8.24.weeks.fiveplus.prop <- nrow(subset(damaging.goodfaith.editors, edits.8.24.weeks.fiveplus==TRUE)) / num.damaging.goodfaith.editors
# 
# d.3.4.weeks.oneplus.prop <- nrow(subset(damaging.editors, edits.3.4.weeks.oneplus==TRUE)) / num.damaging.editors
# d.8.24.weeks.oneplus.prop <- nrow(subset(damaging.editors, edits.8.24.weeks.oneplus==TRUE)) / num.damaging.editors
# 
# g.3.4.weeks.oneplus.prop <- nrow(subset(goodfaith.editors, edits.3.4.weeks.oneplus==TRUE)) / num.goodfaith.editors
# g.8.24.weeks.oneplus.prop <- nrow(subset(goodfaith.editors, edits.8.24.weeks.oneplus==TRUE)) / num.goodfaith.editors
# 
# all.3.4.weeks.oneplus.prop <- nrow(subset(newcomer.editors, edits.3.4.weeks.oneplus==TRUE)) / num.newcomer.editors
# all.4.8.weeks.oneplus.prop <- nrow(subset(newcomer.editors, edits.4.8.weeks.oneplus==TRUE)) / num.newcomer.editors
# all.8.24.weeks.oneplus.prop <- nrow(subset(newcomer.editors, edits.8.24.weeks.oneplus==TRUE)) / num.newcomer.editors
# 
# 
# ## WITHOUT-ML EXPERIMENT INVOLVING A TEAHOUSE INVITE
# without.ml.sample <- num.newcomer.editors * (1 - ml.recruitment.sample)
# 
# ## min.sample as a proportion of one month
# 
# all.3.4.weeks.oneplus.min.sample.without.ml <- min.sample.for.power.binary(all.3.4.weeks.oneplus.prop,
#                                                                     all.3.4.weeks.oneplus.prop + 
#                                                                       min.observed.effect.no.ml,
#                                                                     without.ml.sample  / 4,
#                                                                     without.ml.sample * 12,
#                                                                     without.ml.sample / 4,
#                                                                     alpha.level,
#                                                                     power.goal) / without.ml.sample 
# 
# all.4.8.weeks.oneplus.min.sample.without.ml <- min.sample.for.power.binary(all.4.8.weeks.oneplus.prop,
#                                                                     all.4.8.weeks.oneplus.prop + 
#                                                                       min.observed.effect.no.ml,
#                                                                     without.ml.sample  / 4,
#                                                                     without.ml.sample * 12,
#                                                                     without.ml.sample / 4,
#                                                                     alpha.level,
#                                                                     power.goal) / without.ml.sample 
# 
# 
# all.8.24.weeks.oneplus.min.sample.without.ml <- min.sample.for.power.binary(all.8.24.weeks.oneplus.prop,
#                                                                      all.8.24.weeks.oneplus.prop + 
#                                                                        min.observed.effect.no.ml,
#                                                                      without.ml.sample  / 4,
#                                                                      without.ml.sample * 12,
#                                                                      without.ml.sample / 4,
#                                                                      alpha.level,
#                                                                      power.goal) / without.ml.sample 
# 
# ## WITH-ML EXPERIMENT INVOLVING A TEAHOUSE INVITE
# 
# dge.3.4.weeks.oneplus.min.sample.dg <- min.sample.for.power.binary(dge.3.4.weeks.oneplus.prop,
#                                                             dge.3.4.weeks.oneplus.prop + 
#                                                               min.observed.effect.dg.ml,
#                                                             with.ml.dg.sample  / 4,
#                                                             with.ml.dg.sample * 12,
#                                                             with.ml.dg.sample / 4,
#                                                             alpha.level,
#                                                             power.goal) / with.ml.dg.sample 
# 
# dge.4.8.weeks.oneplus.min.sample.dg <- min.sample.for.power.binary(dge.4.8.weeks.oneplus.prop,
#                                                             dge.4.8.weeks.oneplus.prop + 
#                                                               min.observed.effect.dg.ml,
#                                                             with.ml.dg.sample  / 4,
#                                                             with.ml.dg.sample * 12,
#                                                             with.ml.dg.sample / 4,
#                                                             alpha.level,
#                                                             power.goal) / with.ml.dg.sample 
# 
# 
# dge.8.24.weeks.oneplus.min.sample.dg <- min.sample.for.power.binary(dge.8.24.weeks.oneplus.prop,
#                                                              dge.8.24.weeks.oneplus.prop + 
#                                                                min.observed.effect.dg.ml,
#                                                              with.ml.dg.sample  / 4,
#                                                              with.ml.dg.sample * 12,
#                                                              with.ml.dg.sample / 4,
#                                                              alpha.level,
#                                                              power.goal) / with.ml.dg.sample 
