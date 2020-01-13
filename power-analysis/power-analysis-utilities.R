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
