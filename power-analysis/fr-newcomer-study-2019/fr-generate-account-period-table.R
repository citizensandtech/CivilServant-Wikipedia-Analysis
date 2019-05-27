options("scipen"=9, "digits"=4)
library(dplyr)
library(MASS)
library(ggplot2)
library(rlang)


library(survminer)
library(survival)
## ^^ documentation: https://cran.r-project.org/web/packages/survminer/vignettes/Informative_Survival_Plots.html

## DOCUMENTATION AT: https://cran.r-project.org/web/packages/DeclareDesign/DeclareDesign.pdf
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
options(repr.plot.width=7, repr.plot.height=3.5)
sessionInfo()


### LOAD POWER ANALYSIS DATAFRAMES
data.path <- "~/Tresors/CivilServant/projects/wikipedia-integration/fr-newcomer-study/datasets/"
fr.power.df <- read.csv(file.path(data.path, "french_power-analysis_dataset_sim_date_20180307_v1.csv"))
fr.power.df$user_registration <- as.Date(fr.power.df$user_registration)
fr.power.df$lang <- "fr"

simulated.treatment.date <- as.Date("20180306", "%Y%M%D")

#CONSTRUCT AND REVIEW AN ACCOUNT-PERIOD DATASET FOR SURVIVAL ANALYSIS
initialize_account.periods.df <- function(){
    account.periods.df <- data.frame(
        user_id                       = NA,
        has_email                     = NA,
        has_gender                    = NA,
        registration_date             = NA,
        lang                          = NA,
        week                          = NA,
        num_edits_nonmain             = NA,
        any_edits_nonmain             = NA,
        num_edits_talk                = NA,
        any_edits_talk                = NA,
        num_edits                     = NA,
        any_edits                     = NA,
        labor_hours                   = NA,
        inactive                      = NA
    )
    account.periods.df
}

# Create a person period dataset from a single row of the power analysis data
#
#` @param row The row to convert into a series of rows
#` @param the number of weeks represented in this dataset
create.person.period.rows <- function(row, num.weeks=12){

    ## we use the definition that someone is active
    ## if they edited at least once in a one-month period
    ## this is why we rate active-ness as "was active in
    ## the four week period after the beginning of this week"
    last.person.period.week = 8

    num.edit.nonmain.post.treat.keys = c()
    num.edit.talk.post.treat.keys = c()
    num.edit.post.treat.keys = c()
    labor.hours.post.treat.keys = c()

    for(i in seq(1,12)){
        num.edit.post.treat.keys <- append(num.edit.post.treat.keys,
                                           paste("num_edits_week_",i,"_post_treatment", sep=""))
        num.edit.nonmain.post.treat.keys <- append(num.edit.nonmain.post.treat.keys,
                                           paste("num_edits_nonmain_week_",i,"_post_treatment", sep=""))
        num.edit.talk.post.treat.keys <- append(num.edit.talk.post.treat.keys,
                                           paste("num_edits_talk_week_",i,"_post_treatment", sep=""))
        labor.hours.post.treat.keys <- append(labor.hours.post.treat.keys,
                                           paste("labor_hours_week_",i,"_post_treatment", sep=""))
    }
    #print(paste("CREATED A MULTIWEEK DATASET:", length(num.edit.post.treat.keys)))

    account.periods.df <- initialize_account.periods.df()

    ## iterate backward from the 12th week
    ## and set still.active to true from the
    ## last week that we observe activity in the account
    still.active <- FALSE
    for(i in seq(num.weeks,1)){ 
        num.edit.key = num.edit.post.treat.keys[i]
        num.edit.nonmain.key = num.edit.nonmain.post.treat.keys[i]
        num.edit.talk.key = num.edit.talk.post.treat.keys[i]
        labor.hours.key = labor.hours.post.treat.keys[i]

        if(still.active == FALSE){
            if(row[[num.edit.key]] > 0){
                still.active = TRUE
            }
        }

        inactive <- NA
        if(i<=last.person.period.week){
            inactive <- still.active != TRUE
        }
        
        #print(paste("week: ", i,", inactive: ", inactive, ", still.active: ", still.active, sep=""))
        #print(paste("num.edit.key:", num.edit.key))
        #print(paste("num.edit.talk.key:", num.edit.talk.key))
        #print(paste("num.edit.nonmain.key:", num.edit.nonmain.key))
        #print(paste("num.edit.nonmain.key:", labor.hours.key))
        
        row.df <- data.frame(
            user_id                       = row$user_id,
            has_email                     = row$has_email,
            has_gender                    = row$has_gender,
            registration_date             = row$user_registration,
            lang                          = toString(row$lang),
            week                          = i,
            num_edits                     = row[[num.edit.key]],
            any_edits                     = row[[num.edit.key]] > 0,
            num_edits_talk                = row[[num.edit.talk.key]],
            any_edits_talk                = row[[num.edit.talk.key]] > 0,
            num_edits_nonmain             = row[[num.edit.nonmain.key]],
            any_edits_nonmain             = row[[num.edit.nonmain.key]] > 0,
            labor_hours                   = row[[labor.hours.key]],
            inactive                      = as.integer(inactive)
        )
        account.periods.df <- rbind(account.periods.df, row.df)
    }

    account.periods.df <- subset(account.periods.df, is.na(user_id)!=TRUE)
    account.periods.df
}

# Make a test row's edits and labor hours over N weeks blank
#
#` @row The row to make blank
create.blank.test.row <- function(row){
    for(i in seq(1,12)){
        row[paste("num_edits_week_",i,"_post_treatment", sep="")] <- 0
        row[paste("any_edits_week_",i,"_post_treatment", sep="")] <- FALSE
        row[paste("labor_hours_week_",i,"_post_treatment", sep="")] <- 0
        row[paste("any_labor_hours_week_",i,"_post_treatment", sep="")] <- FALSE

    }
    
    for(i in seq(1,12)){
    }
    row
}

## GENERATE AN ACCOUNT PERIOD DATAFRAME
# generate.account.period.dataframe
# Create an account period dataframe from a power analysis dataframe
# ` 
# `@power.df power analysis dataframe
generate.account.period.dataframe <- function(power.df){
    print(paste("Generating account period dataframe for", unique(power.df$lang)))
    account.periods.df <- initialize_account.periods.df()
    for(i in seq(1,nrow(power.df))){
        if(i %% 100 == 0){
            cat(".")
            flush.console()
        }
        row = power.df[i,]
        account.period.df  <- create.person.period.rows(row)
        account.periods.df <- rbind(account.periods.df, account.period.df)
    }
    cat("\n")
    account.periods.df <- subset(account.periods.df, is.na(user_id)!=TRUE)
    print(paste("Number of participants: ", nrow(power.df), sep=""))
    print(paste("Number of account periods: ", nrow(account.periods.df), sep="")) 
    print(paste("Number periods per account: ", nrow(account.periods.df) / nrow(power.df)), sep="")
    account.periods.df
}

fr.account.periods.df <- generate.account.period.dataframe(fr.power.df)

write.csv(fr.account.periods.df, file="fr.20180306.account.period.dataframe.csv")
