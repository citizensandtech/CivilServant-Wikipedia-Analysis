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


## METHOD FOR ADDING EXPERIENCE LEVEL APPROPRIATE LABOR HOURS FROM THE DATA TO A GIVEN
## CONFIGURATION DF
apply.experience.labor.hours <- function(ref.power.df, config.df){
    for(experience.level in sort(unique(ref.power.df$prev_experience))){
        config.df[paste("exp.",experience.level,".labor.hour.90.day.before.mean", sep="")] <- mean(
            subset(ref.power.df, prev_experience == experience.level)$labour_hours_90_pre_treatment)
        config.df[paste("exp.",experience.level,".labor.hour.90.day.before.sd", sep="")] <- sd(
            subset(ref.power.df, prev_experience == experience.level)$labour_hours_90_pre_treatment)
    }
    config.df
}



# generate.sample.from.reference
# Return a sample of N size based on sampling with (or without) replacement
# from an existing dataframe
#` @param ref.power.df The dataframe to sample from with replacement
#` @param sample.size  The number of observations to sample
#` @param included.columns The columns to include in the dataframe
#` @param replacement Whether to sample with replacement (default TRUE)
generate.sample.from.reference <- function(ref.power.df, sample.size, included.columns,
                                                      replacement=TRUE){
    sim.df <- data.frame(id = seq(sample.size))
    for(colname in included.columns){
        sim.df[colname] <- NA
    }

    num.rows <- nrow(ref.power.df)
    subset.df <- ref.power.df[sample(num.rows, sample.size, replace=replacement),]
   
    for(colname in included.columns){
        sim.df[,colname] <-  subset.df[,colname]
    }
    sim.df
}

# generate.stratified.sample.from.reference
# Return a sample of N size based on sampling with (or without) replacement
# from an existing dataframe
#                                                                       
#` @param ref.power.df The dataframe to sample from with replacement                               
#` @param sample.per.experience.group The sample size per prev_experience group
#` @param included.columns The columns to include in the dataframe
#` @param replacement Whether to sample with replacement (detault TRUE)

generate.stratified.sample.from.reference <- function(ref.power.df, 
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
