#~/usr/bin/env Rscript
## COMMAND TO RUN THIS SCRIPT
## Rscript --vanilla thanks-recipient-experiment-power-analysis-terminal.R <LANG>
## where <LANG> is "de" "fa" "pl" or "ar"
args = commandArgs(trailingOnly=TRUE)

scriptlang = args[1]
if((scriptlang %in% c("de", "fa", "pl", "ar"))!=TRUE){
    stop("Language argument needs to be de, fa, pl, or ar")
}

options("scipen"=9, "digits"=4)
library(dplyr)
library(MASS)
library(ggplot2)
library(rlang)
library(gmodels)
library(tidyverse)
library(viridis)
library(fabricatr)
library(DeclareDesign)
library(blockTools)
## Installed DeclareDesign 0.13 using the following command:
# install.packages("DeclareDesign", dependencies = TRUE,
#                 repos = c("http://R.declaredesign.org", "https://cloud.r-project.org"))

library(survminer)
library(survival)
## ^^ documentation: https://cran.r-project.org/web/packages/survminer/vignettes/Informative_Survival_Plots.html

## DOCUMENTATION AT: https://cran.r-project.org/web/packages/DeclareDesign/DeclareDesign.pdf
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
options(repr.plot.width=7, repr.plot.height=3.5)
sessionInfo()


# Return the difference in mu associated with an incidence rate ratio
# from a negative binomial model. This difference can then be used to
# simulate a negative binomial distribution for the effect of a given irr
#                                                                       
#` @param mu The baseline mu in question                               
#` @param irr The incidence rate ratio in question
mu.diff.from.mu.irr <- function(mu, irr){
    mu*(irr-1)
}

# Return the total sum of betas for a
# logistic regression, given a probability
#
#` @param p the probability in question
betas.logit.from.prob <- function(p){
    log(p/(1-p))
}


# Return the total sum of betas for a
# logistic regression, given a probability
#
#` @param Y list of observed Ys
betas.logit.from.mean <- function(Y){
    p = mean(Y)
    log(p/(1-p))
}

# Return the minimum power reported in a diagnosis
# 
#` @param diagnosis
min.diagnosis.power <- function(diagnosis){
    min(diagnosis$diagnosands_df['power'])
}

# Conduct a binary search for a certain level of statistical power
# within the constraints of a configuration file
#
#` @param ref.df:             The dataframe to draw from
#` @param config.df:          The configuration file in question
#` @param survival.tables:    Baseline data on survival rates
#` @parra d.lang:             What language is being tested (must be a value of lang in survival tables) 
#` @diagnosis.method:         The method that conducts a single DeclareDesign diagnosis and returns the diagnosis
#` @target.power:             The statistical power that ideally should be the minimum across the study
#` @target.tolerance:         How close to the desired statistical power is close enough?
#` @min.sample.diff:          If the search is close enough that the change is less than min.sample.diff, end the search
#` @start.sample.size:        If specified, use the starting value as the initial sample size to use

search.for.power <- function(ref.df, config.df, survival.tables, d.lang,
                             diagnosis.method = diagnose.experiment, 
                             target.power = 0.85, target.tolerance = 0.01, 
                             min.sample.diff = 100,
                             start.sample.size = NA){
    max.sample.size = config.df$n.max
    min.sample.size = config.df$n.min
    if(is.na(start.sample.size)){
        current.sample.size = as.integer(max.sample.size / 2)
    }else{
        current.sample.size = start.sample.size
    }
    current.power = 0.0

    ## Initialize first iteration
    num.prev.experience.groups <- length(unique(ref.df$prev_experience))
    print(paste("experience groups:", num.prev.experience.groups, 
                ".min:", min.sample.size * num.prev.experience.groups, 
                ".max:", max.sample.size * num.prev.experience.groups,
                ".current n per experience group:",  current.sample.size, 
                ".current total n:", current.sample.size * num.prev.experience.groups))
    flush.console()

    ptm = proc.time() #record time the simulation started
    ddf <- diagnosis.method(current.sample.size, ref.df, pa.config = config.df,
                            survival.tables, d.lang)
    ddf$diagnosands$n <- current.sample.size * num.prev.experience.groups
    diagnoses.df = ddf$diagnosands
    current.power <- min.diagnosis.power(ddf)

    ## output elapsed time for this iteration
    time.elapsed <- proc.time() -  ptm
    print(paste("     seconds:", as.integer(time.elapsed['elapsed'])))

    
    if(current.power < target.power){
        min.sample.size = current.sample.size
        print(paste(current.power, "<", target.power))
    }else{
        max.sample.size = current.sample.size
        print(paste(current.power, ">", target.power))
    }

    current.sample.size = min.sample.size + as.integer((max.sample.size - min.sample.size)/2)
    while(all.equal(target.power, current.power, tolerance = target.tolerance)!=TRUE){
        print(paste("experience groups:", num.prev.experience.groups, 
                    ".min:", min.sample.size * num.prev.experience.groups, 
                    ".max:", max.sample.size * num.prev.experience.groups,
                    ".current n per experience group:",  current.sample.size, 
                    ".current total n:", current.sample.size * num.prev.experience.groups))
        flush.console()

        ptm = proc.time() #record time the simulation started

        ## conduct simulations
        ddf <- diagnosis.method(current.sample.size, ref.df, pa.config = config.df,
                                survival.tables, d.lang)
        ddf$diagnosands$n <- current.sample.size * num.prev.experience.groups
        ## append simulation results to dataframe
        diagnoses.df <- rbind(diagnoses.df, ddf$diagnosands)
        
        ## output elapsed time for this iteration
        time.elapsed <- proc.time() -  ptm
        print(paste("     seconds:", as.integer(time.elapsed['elapsed'])))


        ## check the current statistical power and
        ## carry out the binary search by first updating the boundaries
        current.power <- min.diagnosis.power(ddf)
        if(current.power < target.power){
            min.sample.size = current.sample.size
            print(paste(current.power, "<", target.power))
        }else{
            max.sample.size = current.sample.size
            print(paste(current.power, ">", target.power))
        }
        ## update the current pointer, or break if
        ## the sample size difference is less than or equal
        ## to ten
        sample.size.diff <- as.integer((max.sample.size - min.sample.size)/2)
        if(abs(sample.size.diff) <= min.sample.diff){
            print(paste("Sample size difference ", abs(sample.size.diff), 
                        " <= ", min.sample.diff, ". Ending cycle.", sep=""))
            break
        }
        current.sample.size = min.sample.size + sample.size.diff
    }
    diagnoses.df
}

# Iterate linearly for a certain level of statistical power
# within the constraints of a configuration file
# at a certain sample size increment. Useful for
# illustrating ideas, or for comparing estimators with
# very different statistical power, where the binary search
# will optimize for the worst estimator but not show useful
# indormation about more efficient estimators
#
#` @param ref.df:             The reference dataframe to use
#` @param config.df:          The configuration file in question
#` @param survival.tables:    Baseline data on survival rates
#` @parra d.lang:            What language is being tested (must be a value of lang in survival tables) 
#` @param diagnosis.method:   The method that conducts a single DeclareDesign diagnosis and returns the diagnosis
#` @param iteration.interval: When iterating, use this interval between sample sizes


iterate.for.power <- function(ref.df, config.df, survival.tables, d.lang,
                              diagnosis.method = diagnose.experiment, 
                              iteration.interval){
    
    max.sample.size = config.df$n.max
    min.sample.size = config.df$n.min
    current.sample.size = min.sample.size
    
    iteration.count = ceiling((max.sample.size - min.sample.size) / iteration.interval)

    ## Initialize first iteration
    num.prev.experience.groups <- length(unique(ref.df$prev_experience))
    
    print(paste("experience groups:", num.prev.experience.groups, 
                ".min:", min.sample.size * num.prev.experience.groups, 
                ".max:", max.sample.size * num.prev.experience.groups,
                ".current n per experience group:",  current.sample.size, 
                ".current total n:", current.sample.size * num.prev.experience.groups))
    flush.console()

    ptm = proc.time()
    ddf <- diagnosis.method(current.sample.size, ref.df, pa.config = config.df,
                            survival.tables, d.lang)
    ddf$diagnosands$n <- current.sample.size
    diagnoses.df = ddf$diagnosands
    current.power <- min.diagnosis.power(ddf)
    time.elapsed <- proc.time() -  ptm
    print(paste("     seconds:", as.integer(time.elapsed['elapsed'])))
    
    for(i in seq(1, iteration.count)){
        current.sample.size = current.sample.size + iteration.interval
        print(paste("experience groups:", num.prev.experience.groups, 
                    ".min:", min.sample.size * num.prev.experience.groups, 
                    ".max:", max.sample.size * num.prev.experience.groups,
                    ".current n per experience group:",  current.sample.size, 
                    ".current total n:", current.sample.size * num.prev.experience.groups))
        flush.console()
    
        ptm = proc.time()
        ## conduct simulations
        ddf <- diagnosis.method(current.sample.size, ref.df, pa.config = config.df,
                                survival.tables, d.lang)
        ddf$diagnosands$n <- current.sample.size * num.prev.experience.groups
        ## append simulation results to dataframe
        diagnoses.df <- rbind(diagnoses.df, ddf$diagnosands)
        time.elapsed <- proc.time() -  ptm
        print(paste("     seconds:", as.integer(time.elapsed['elapsed'])))
    }
    diagnoses.df
}

# Create a plot of a power search or iteration output
# Especially useful in cases with multiple DVs or estimators
#
#` @param diagnoses Dataframe of diagnosis info
#` @param config.df the power analysis config dataframe

plot.power.results <- function(diagnoses, config.df){
    for(estimator_label in unique(diagnoses$estimator_label)){
        estimator.diagnoses <- diagnoses[diagnoses$estimator_label==estimator_label,]
        estimator_min_sample = min(estimator.diagnoses$n[estimator.diagnoses$power>0.8])
        p <- ggplot(data=estimator.diagnoses, aes(n, power)) +
                geom_point(color="coral") +
                xlab("sample size (horizontal line represents 85% chance of p < 0.05)") +
                ylim(0,1) +
                geom_hline(aes(yintercept=0.85), linetype="dashed") +
                theme_light() +
                ggtitle(paste(config.df$pa.label, ": statistical Power for Estimator ", estimator_label, "\n",
                              "Minimum sample: ", estimator_min_sample, sep="")) #+
                ggsave(paste("figures/power.analysis.", make.names(estimator_label), ".", config.df$pa.label, ".png", sep=""))
        print(p)
    }
}

data.path <- "~/Tresors/CivilServant/projects/wikipedia-integration/gratitude-study/datasets/power_analysis"
de.power.df <- read.csv(file.path(data.path, "de_gratitude_power-analysis_dataset_sim_date_20180306_with_email.csv"))
fa.ar.pl.power.df <- read.csv(file.path(data.path, "gratitude_power-analysis_dataset_sim_date_20180306_with_email.csv"))
fa.power.df <- subset(fa.ar.pl.power.df, lang=="fa")
ar.power.df <- subset(fa.ar.pl.power.df, lang=="ar")
pl.power.df <- subset(subset(fa.ar.pl.power.df, lang=="fa"))

simulated.treatment.date <- as.Date("20180306", "%Y%M%D")

subset.and.review.variables <- function(df.to.review){
    print(      "========================")
    print(paste("Review Variables For:", unique(df.to.review$lang)))
    print(      "========================")
    cat("\n")
    
    print(paste("Total rows:", nrow(df.to.review)))
    cat("\n")

    # VARIABLE: experience_level_pre_treatment: Definition: 
    #  the elapsed number of days between registration 
    #  and last edit in the observation period up to the simulated treatment date

    df.to.review$prev_experience <- as.integer(gsub("bin_", "", df.to.review$experience_level_pre_treatment))
    df.to.review$prev_experience <- factor(df.to.review$prev_experience, 
                                          levels = sort(unique(df.to.review$prev_experience)))
    
    df.to.review$num_prev_thanks_pre_treatment <- df.to.review$num_prev_thanks_before_treatment
    
    ## SHOW NUM EDITS BY EXPERIENCE GROUP:
    print("Aggregate num_edits_90_pre_treatment")
    print(aggregate(df.to.review[c("num_edits_90_pre_treatment")]>0,
              FUN=mean, by = list(df.to.review$prev_experience)))
    cat("\n")

    ## SHOW LABOR HOURS BY EXPERIENCE GROUP:
    print("Aggregate labor_hours_90_pre_treatment")
    print(aggregate(df.to.review[c("labour_hours_90_pre_treatment")],
              FUN=mean, by = list(df.to.review$prev_experience)))
    cat("\n")

    ## BECAUSE THIS DATASET INCLUDES INACTIVE EDITORS
    ## WE REDUCE THE DATASET ONLY TO EDITORS ACTIVE IN THE LAST 90 DAYS
    print(paste("Number of rows before removing inactive users:", nrow(df.to.review)))
    df.to.review <- subset(df.to.review, num_edits_90_pre_treatment > 0 )
    print(paste("Number of rows after removing inactive users:", nrow(df.to.review)))
    cat("\n")
    
    
    print("prev_experience")
    print(summary(factor(df.to.review$prev_experience)))
    cat("\n")

    # VARIABLE: newcomer: accounts that were created within in the last 90 days, 
    #   and which have made at least 3 edits since starting their account
    df.to.review$newcomer <- df.to.review$experience_level_pre_treatment == "bin_0" & 
        df.to.review$num_edits_90_pre_treatment >= 3

    
    print("BLOCKING VARIABLES")
    print("--------------------")
    
    
    print("Number of Newcomers with 3 or more edits")
    print(summary(subset(df.to.review, experience_level_pre_treatment == "bin_0")$newcomer))
    cat("\n")

    print("NEWCOMERS AND EMAILS")
    print("--------------------")
    print(CrossTable(df.to.review$has_email, df.to.review$newcomer, 
           prop.r = FALSE, prop.c=TRUE, prop.t = FALSE, prop.chisq = FALSE))
    
    df.to.review$has_email <- df.to.review$has_email == "True"
    
    
    # VARIABLE: num_prev_thanks_pre_treatment
    print("num_prev_thanks_pre_treatment")
    print(summary(df.to.review$num_prev_thanks_pre_treatment))
    cat("\n")
    
    ## SHOW PREVIOUS THANKS BY EXPERIENCE GROUP:
    print("num_prev_thanks_pre_treatment by prev_experience")
    print(aggregate(df.to.review[c("num_prev_thanks_pre_treatment")],
          FUN=mean, by = list(df.to.review$prev_experience)))
    cat("\n")

    # new section: behavioral variables
    print("BEHAVIORAL VARIABLES")
    print("--------------------")
    cat("\n")

    # VARIABLE: num.edit.diff: the number of edits
    # in the 90 days after treatment minutes the number in the 90 before treatment
    df.to.review$num.edit.diff <- df.to.review$num_edits_90_post_treatment - 
        df.to.review$num_edits_90_pre_treatment

    # VARIABLE: num.edit.diff: the number of edits
    # in the 90 days after treatment minutes the number in the 90 before treatment
    df.to.review$labor.hour.diff <- df.to.review$labour_hours_90_post_treatment - 
        df.to.review$labour_hours_90_pre_treatment
    
    ## SHOW NUM EDIT DIFF BY EXPERIENCE GROUP:
    print("Show pre/post difference in 90 day num edits by experience group")
    print(aggregate(df.to.review[c("num.edit.diff")],
          FUN=mean, by = list(df.to.review$prev_experience)))
    cat("\n")

    
    ## SHOW LABOR HOUR DIFF BY EXPERIENCE GROUP:
    print("Show pre/post difference in 90 day labor hours by experience group")
    print(aggregate(df.to.review[c("labor.hour.diff")],
              FUN=mean, by = list(df.to.review$prev_experience)))
    cat("\n")

    print(ggplot(df.to.review, aes(prev_experience, labor.hour.diff, color=prev_experience)) +
        geom_boxplot() +
        theme_bw() +
        ggtitle(paste("Difference in 90 day Pre/PostLabor Hours by Experience:", unique(df.to.review$lang), "dataset"))
    )
    
    df.to.review    
}

de.power.sub.df <- subset.and.review.variables(de.power.df)
fa.power.sub.df <- subset.and.review.variables(fa.power.df)
ar.power.sub.df <- subset.and.review.variables(ar.power.df)
pl.power.sub.df <- subset.and.review.variables(pl.power.df)

survival.tables <- read.csv("fa.ar.pl.de.20180306.survival.tables.csv")

# Method: apply.experience.labor.hours
#
#`@ref.power.df THe reference dataframe
#`@config.df 
apply.experience.labor.hours <- function(ref.power.df, config.df){
    for(experience.level in sort(unique(ref.power.df$prev_experience))){
        config.df[paste("exp.",experience.level,".labor.hour.90.day.before.mean", sep="")] <- mean(
            subset(ref.power.df, prev_experience == experience.level)$labour_hours_90_pre_treatment)
        config.df[paste("exp.",experience.level,".labor.hour.90.day.before.sd", sep="")] <- sd(
            subset(ref.power.df, prev_experience == experience.level)$labour_hours_90_pre_treatment)
    }
    config.df
}

pa.config.newcomer <- data.frame(
    n.max    = 2000, # max number of observations
    n.min    = 150, # min number of observations
    
    survey.week.interval = 2, 
    # simulated survey completion rate of all those
    # who were active on Wikipedia at that time, among the enrolled participants
    survey.participation.rate = 0.2, # 20%, via Julia
        
    # QUESTION: DO WE EXPECT THE EFFECT TO BE MULTIPLICATIVE OF THE LABOR HOURS
    #           OR THE SAME AMOUNT, NO MATTER HOW MANY LABOR HOURS A PERSON TYPICALLY PUTS IN?
    # ANSWER:   LET'S LOOK FOR A MINIMUM ADDITIVE LABOR HOURS, SINCE THERE IS A MINIMUM OBSERVABLE
    #           AMOUNT BASED ON THE WAY THE MEASURE IS CONSTRUCTED

    # the minimum effect we want to be able to observe for any group:
    # a thirty minute increase in the number of labor hours over 90 days
    # assumption, the treatment group has the same standard deviation
    labor.hour.diff.90.day.treat.mean.effect = 0.25,
    labor.hour.diff.90.day.treat.sd.effect = 0.5,
    
    ## EFFECT ON NUMBER OF THANKS GIVEN
    ## we are going to estimate the average treatment effect
    ## since negbin and logistic will be difficult to estimate
    thanks.given.90.day.placebo.value = 0,
    thanks.given.90.day.treat.theta = 0.2,
    thanks.given.90.day.treat.mu = 0.1,
    thanks.given.90.day.estimand = mean(rnbinom(10000, 0.2, mu = 0.1)),
    
    # QUESTIONS FOR SURVIVAL
    # These assumptions are based on the analysis in survival-analysis-of-power-analysis-dataset-R.ipynb
    # ASSUMPTION: THE EFFECT IS NOT CUMULATIVE; NO EFFECT ON THE SLOPE
    # ASSUMPTION: THE EFFECT ON RETENTION IS GREATER FOR LESS EXPERIENCED WIKIPEDIANS
    # ASSUMPTION: THE RELATIONSHIP BETWEEN EFFECT AND RETENTION IS CURVED/ASYMPTOTIC TO A MINIMUM EFFECT
    
    # we want to be able to observe at least a 1 percentage point increase in survival
    # for treated participants (knowing that not everyone will get treated)
    # survival.effect = 1.01,
    
    ## SELF-EFFICACY (-3 to 3)
    # We don't expect an effect but we want to see an effect
    # of at least 0.5 if any exists (<- TODO: review assumption)
    survey.efficacy.placebo.mean = -2,
    survey.efficacy.placebo.sd   = 1,
    survey.efficacy.treat.effect = 0.5,

    ## RELATIONSHIP WITH WIKIPEDIA COMMUNITY (1 to 6)
    # W want to observe an effect
    # of at least 0.5 if any exists (<- TODO: review assumption)
    survey.closeness.placebo.mean = 2,
    survey.closeness.placebo.sd   = 1,
    survey.closeness.treat.effect = 0.5,

    ## INDEX OF SOCIAL VALUE (1 to 6)
    # Combination of:
    # "My contributions are valued by other Wikipedians."
    # "My contributions have made positive difference for other Wikipedians."
    # "How much would you say the community overall is friendly?"
    # "How much would you say the community overall is supportive?"
    # We want to see an effect
    # of at least 0.5 if any exists (<- TODO: review assumption)
    survey.socialvalue.placebo.mean = 2,
    survey.socialvalue.placebo.sd   = 1,
    survey.socialvalue.treat.effect = 0.5    
    
)

pa.config.experienced <- data.frame(
    n.max    = 6000, # max number of observations per group
    n.min    = 100, # min number of observations per group (multiple experience groups)
    # iterate by 40?
    
    survey.week.interval = 2, 
    survey.participation.rate = 0.2, ## per Julia's research
        
    # QUESTION: DO WE EXPECT THE EFFECT TO BE MULTIPLICATIVE OF THE LABOR HOURS
    #           OR THE SAME AMOUNT, NO MATTER HOW MANY LABOR HOURS A PERSON TYPICALLY PUTS IN?
    # ANSWER:   LET'S LOOK FOR A MINIMUM ADDITIVE LABOR HOURS, SINCE THERE IS A MINIMUM OBSERVABLE
    #           AMOUNT BASED ON THE WAY THE MEASURE IS CONSTRUCTED

    # the minimum effect we want to be able to observe for any group:
    # a thirty minute increase in the number of labor hours over 90 days
    # assumption, the treatment group has the same standard deviation
    labor.hour.diff.90.day.treat.mean.effect = 0.25,
    labor.hour.diff.90.day.treat.sd.effect = 0.5,
    
    ## EFFECT ON NUMBER OF THANKS GIVEN
    ## we are going to estimate the average treatment effect
    ## since negbin and logistic will be difficult to estimate
    thanks.given.90.day.placebo.value = 0,
    thanks.given.90.day.treat.theta = 0.2,
    thanks.given.90.day.treat.mu = 0.1,
    thanks.given.90.day.estimand = mean(rnbinom(10000, 0.2, mu = 0.1)),
    
    # QUESTIONS FOR SURVIVAL
    # 
    # ASSUMPTION: THE EFFECT IS NOT CUMULATIVE; NO EFFECT ON THE SLOPE
    # ASSUMPTION: THE EFFECT ON RETENTION IS GREATER FOR LESS EXPERIENCED WIKIPEDIANS
    # ASSUMPTION: THE RELATIONSHIP BETWEEN EFFECT AND RETENTION IS CURVED/ASYMPTOTIC TO A MINIMUM EFFECT
    
    # we want to be able to observe at least a 1 percentage point increase in survival
    # for treated participants (knowing that not everyone will get treated)
    # survival.effect = 1.01,
    
    ## SELF-EFFICACY (-3 to 3)
    # We don't expect an effect but we want to see an effect
    # of at least 0.5 if any exists (<- TODO: review assumption)
    survey.efficacy.placebo.mean = 0,
    survey.efficacy.placebo.sd   = 1,
    survey.efficacy.treat.effect = 0.5,

    ## RELATIONSHIP WITH WIKIPEDIA COMMUNITY (1 to 7)
    # W want to observe an effect
    # of at least 0.5 if any exists (<- TODO: review assumption)
    survey.closeness.placebo.mean = 3,
    survey.closeness.placebo.sd   = 1,
    survey.closeness.treat.effect = 0.5,

    ## INDEX OF SOCIAL VALUE (1 to 7)
    # Combination of:
    # "My contributions are valued by other Wikipedians."
    # "My contributions have made positive difference for other Wikipedians."
    # "How much would you say the community overall is friendly?"
    # "How much would you say the community overall is supportive?"
    # We want to see an effect
    # of at least 0.5 if any exists (<- TODO: review assumption)
    survey.socialvalue.placebo.mean = 3,
    survey.socialvalue.placebo.sd   = 1,
    survey.socialvalue.treat.effect = 0.5    
    
)

# generate.sample.from.reference
# Return a sample of N size based on sampling with (or without) replacement
# from an existing dataframe
#                                                                       
#` @param ref.power.df The dataframe to sample from with replacement                               
#` @param sample.per.experience.group The sample size per prev_experience group
#` @param included.columns The columns to include in the dataframe
#` @param replacement Whether to sample with replacement

generate.sample.from.reference <- function(ref.power.df, 
                                             sample.per.experience.group, 
                                             included.columns = c('labor_hours_90_pre_treatment', 
                                                                  'labor_hours_90_post_treatment',
                                                                  'num_prev_thanks_pre_treatment',
                                                                  'has_email',
                                                                  'newcomer'),
                                             replacement=TRUE){

    ## if there's a misspelling in the column names
    ## correct the misspelling
    if(('labor_hours_90_pre_treatment' %in% colnames(ref.power.df))!=TRUE){
        ref.power.df$labor_hours_90_pre_treatment <- ref.power.df$labour_hours_90_pre_treatment
        ref.power.df$labor_hours_90_post_treatment <- ref.power.df$labour_hours_90_post_treatment
    }
    
    sim.df <- expand.grid(id = seq(sample.per.experience.group), prev_experience = unique(ref.power.df$prev_experience))
    sim.df$id <- seq.int(nrow(sim.df))
    
    for(colname in included.columns){
        sim.df[colname] <- NA
    }

    ## sample with replacement from observational data
    for(pe in unique(ref.power.df$prev_experience)){
        num.rows <- nrow(subset(ref.power.df, prev_experience==pe))
        subset.df <- subset(ref.power.df, prev_experience==pe)[sample(num.rows, sample.per.experience.group, replace=replacement),]

        for(colname in included.columns){
            sim.df[sim.df$prev_experience == pe,colname] <-  subset.df[,colname]
        }
    }
    
    sim.df
}


# Diagnose Experiment: take in a reference dataframe and configuration 
# and diagnose the design
#
#` @param sample.size.per.group: The sample size per prev_experience group (equal)
#` @param ref.df: the dataframe to draw from
#` @param pa.config: the configuration dataframe to use
#` @sims.count: the number of simulations to conduct and aggregate
#` @bootstrap.sims.count: the number of bootstraps to perform for estimating 
#                         confidence intervals for the generated diagnoses

diagnose.experiment <- function(sample.size.per.group, ref.df, pa.config, 
                                survival.tables, d.lang, 
                                sims.count=500, bootstrap.sims.count=500){
        
    pa.config <- apply.experience.labor.hours(ref.df, pa.config)
    d.df <- generate.sample.from.reference(ref.df, sample.size.per.group, replacement=TRUE)
    sample.size <- nrow(d.df)
    
    experienced.survival <- subset(survival.tables, lang==d.lang & week==pa.config$survey.week.interval & experience=="experienced")$active
    newcomer.survival.diff <- experienced.survival - subset(survival.tables, lang==d.lang & week==pa.config$survey.week.interval & experience=="newcomer")$active

    
    ## SET UP BLOCKS
    # do not include groups if there's only one
    if(length(unique(d.df$prev_experience))>1){

        d.df$blocks <- createBlockIDs(obj = block(data=d.df,
                                                  n.tr = 2,
                                                  id.vars="id",
                                                  groups="prev_experience",
                                                  block.vars = c("labor_hours_90_pre_treatment", "num_prev_thanks_pre_treatment"),
                                                  distance ="mahalanobis"
                                                  ),
                                       data=d.df,
                                       id.var = "id")
    }else{ 
        d.df$blocks <- createBlockIDs(obj = block(data=d.df,
                                                  n.tr = 2,
                                                  id.vars="id",
                                                  block.vars = c("labor_hours_90_pre_treatment", "num_prev_thanks_pre_treatment"),
                                                  distance ="mahalanobis"
                                                  ),
                                       data=d.df,
                                       id.var = "id")        
    }
    
    ## DEFINE THE DESIGN
    
    
    design <- 
        declare_population(
            data = d.df,
        ## explaining survey.participation:
        # IF: the account has an email address, we assume the messasge can reach them
        #     and assign their participation probability to be the configured rate
        # ELSE: we consider the probability of participation to be the configured rate
        #       among accounts that are still editing a given Wikipedia at that point, as follows:
        # proportion survived = experienced.survival + newcomer*newcomer.survival.diff
        # participation rate  = participation rate * proportion survived (per language)
        survey.participation = rbinom(nrow(d.df), 1, ifelse(d.df$has_email,
                                                            pa.config$survey.participation.rate,
                                                            pa.config$survey.participation.rate*(experienced.survival + 
                                                                                          newcomer.survival.diff*as.integer(d.df$newcomer))
                                                            ))
        ) +
        declare_potential_outcomes(

            ## survey efficacy on a -3 to 3 scale
            SE_Z_0 = draw_ordered(x=rnorm(nrow(d.df), 
                                          pa.config$survey.efficacy.placebo.mean, 
                                          pa.config$survey.efficacy.placebo.sd), 
                                  breaks = c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5))-4,
            SE_Z_1 = draw_ordered(x=rnorm(nrow(d.df), 
                                          pa.config$survey.efficacy.placebo.mean +
                                          pa.config$survey.efficacy.treat.effect, 
                                          pa.config$survey.efficacy.placebo.sd), 
                                  breaks = c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5))-4,
            ## Survey closeness on a 1 to 7 scale
            SC_Z_0 = draw_ordered(x=rnorm(nrow(d.df), 
                                          pa.config$survey.closeness.placebo.mean, 
                                          pa.config$survey.closeness.placebo.sd), 
                                  breaks = c(2, 3, 4, 5, 6)),
            SC_Z_1 = draw_ordered(x=rnorm(nrow(d.df), 
                                          pa.config$survey.closeness.placebo.mean +
                                          pa.config$survey.closeness.treat.effect, 
                                          pa.config$survey.closeness.placebo.sd), 
                                  breaks = c(2, 3, 4, 5, 6)),

            ## Survey socialvalue will be a continuous index
            SSV_Z_0 = rnorm(nrow(d.df), 
                                          pa.config$survey.socialvalue.placebo.mean, 
                                          pa.config$survey.socialvalue.placebo.sd),
            SSV_Z_1 = rnorm(nrow(d.df), 
                                          pa.config$survey.socialvalue.placebo.mean +
                                          pa.config$survey.socialvalue.treat.effect, 
                                          pa.config$survey.socialvalue.placebo.sd),

            ## labor hour difference
            LHD_Z_0 = labor_hours_90_post_treatment - labor_hours_90_pre_treatment,
            LHD_Z_1 = labor_hours_90_post_treatment - labor_hours_90_pre_treatment + 
                      rnorm(nrow(d.df), pa.config$labor.hour.diff.90.day.treat.mean.effect,
                            pa.config$labor.hour.diff.90.day.treat.sd.effect),
            ## Thanks Given
            TG_Z_0 = pa.config$thanks.given.90.day.placebo.value,
            TG_Z_1 = rnbinom(nrow(d.df), pa.config$thanks.given.90.day.treat.theta, 
                             mu = pa.config$thanks.given.90.day.treat.mu)

        ) +
        declare_assignment(prob = .5, blocks = blocks) +
        declare_estimand(ate_SE_1_0 = pa.config$survey.efficacy.treat.effect,
                         ate_SC_1_0 = pa.config$survey.closeness.treat.effect,
                         ate_SSV_1_0 = pa.config$survey.socialvalue.treat.effect,
                         ate_LHD_1_0 = pa.config$labor.hour.diff.90.day.treat.mean.effect,
                         ate_TG_1_0  = pa.config$thanks.given.90.day.estimand) +
        declare_reveal(outcome_variables = c("SE", "SC", "SSV", "LHD", "TG"), assignment_variables=c("Z")) +

        # in survey estimators, we include all that participated in the survey 
        declare_estimator(formula = SE ~ Z,  
            model    = difference_in_means,
            subset   = survey.participation == 1, 
            estimand = "ate_SE_1_0", 
            label    = "estimate-S_Efficacy_1_0-participated") +

        declare_estimator(formula = SC ~ Z,  
            model    = difference_in_means,
            subset   = survey.participation == 1, 
            estimand = "ate_SC_1_0", 
            label    = "estimate-S_Closeness_1_0-participated") +
    
        declare_estimator(formula = SSV ~ Z,  
            model    = difference_in_means,
            subset   = survey.participation == 1, 
            estimand = "ate_SSV_1_0", 
            label    = "estimate-S_SocialValue_1_0-participated") +

        declare_estimator(formula = LHD ~ Z,  
            model    = difference_in_means,
            blocks   = blocks,
            estimand = "ate_LHD_1_0", 
            label    = "estimate-LaborHoursDiff_1_0-blocked")  +
    
    ## for the power analysis, we only look at the effect among participants
    ## in blocks where no one had previously received thanks
    ## (we could do a followup where we do a linear interaction effect test on number of thanks)
    ## (if there are enough accounts that have previously received thanks)
    declare_estimator(label = "estimate-ThanksGiven_1_0-blocked-nothanks", 
        estimand="ate_TG_1_0",
        handler=tidy_estimator(function(data){
            eligible.blocks <- subset(aggregate(data$num_prev_thanks_pre_treatment, FUN=mean, by=list(data$blocks)), x==0)$Group.1
            data$eligible.block <- data$blocks %in% eligible.blocks
            tidy(difference_in_means(formula = TG ~ Z, data, 
                                     blocks = data$blocks,
                                     subset = (eligible.block==TRUE)))
        })) 
   
    ## CONDUCT THE DIAGNOSIS
    diagnosis <- diagnose_design(design, sims = sims.count, 
                                     bootstrap_sims = bootstrap.sims.count)
    diagnosis
}

print(paste("Running newcomer Power Analysis for ", scriptlang, sep=""))
pa.config.newcomer$pa.label <- paste(scriptlang, "-newcomer", sep="")
newcomer.pa.results <- iterate.for.power(subset(de.power.sub.df, prev_experience==0),
                                         pa.config.newcomer, 
                                         survival.tables, scriptlang,
                                         diagnose.experiment, 50) 
write.csv(newcomer.pa.results, file=paste("data/thankee-power-analysis-newcomer-", pa.config.newcomer$pa.label, ".csv", sep=""))
plot.power.results(newcomer.pa.results, pa.config.newcomer)


print(paste("Running Experienced Power Analysis for ", scriptlang, sep=""))
pa.config.experienced$pa.label <- "DE-experienced"
experienced.pa.results <- iterate.for.power(subset(de.power.sub.df, prev_experience!=0),
                                             pa.config.experienced, 
                                             survival.tables, scriptlang,
                                             diagnose.experiment, 300)
write.csv(experienced.pa.results, file=paste("data/thankee-power-analysis-experienced-", pa.config.newcomer$pa.label, ".csv", sep=""))
plot.power.results(experienced.pa.results, pa.config.experienced)
