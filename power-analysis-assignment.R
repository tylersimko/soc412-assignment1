###########################
## SOC412: Power Analysis
## Emily Hedlund
## Tyler Simko
## 3/2/18
###########################

# setwd("~/Dropbox/princeton/spring18/soc412/power-analysis/")

## LOAD LIBRARIES
library(ggplot2)
library(corrplot)
library(gmodels)
library(texreg)
library(rms)
library(blockrand)
library(lubridate)
library(psych)

## CLEAR PROJECT
rm(list=ls())

###########################
## SUBREDDIT DATA ANALYSIS
##
## This dataset is a full set of posts that appeared in a subreddit in 2017
## Authors, post IDs, and dates have been anonymized or fuzzed
##
## COLUMNS:
## anon.author               - anonymized integer indicating the author
## id                        - anonymized post ID. Do not assume this is sorted by time
## author.prev.participation - previous comments and posts in the last six months by the author
## author.prev.posts         - previous posts by the author in the last six months
## front.page                - number of minutes that the post spent on the front page of reddit
## is.selftext               - a post with text rather than a link to a third party website
## newcomer.comments         - number of comments made by someone who hadn't previously commented in the last six months
## newcomer.comments.removed - number of comments by newcomers that were removed
## num.comments              - total number of comments in discussion
## num.comments.removed      - total number of comments removed in the discussion
## visible                   - final state of whether the post was allowed to remain visible by moderators or removed
## weekday                   - day of the week of the post
## weekend                   - whether the post was made on a weekend or not

####################################
## LOAD DATAFRAMES                ##
####################################
postsData <- read.csv("subreddit_posts.csv")
postsData$is.selftext <- postsData$is.selftext == "True"
postsData$weekend <- (postsData$weekday == 5 | postsData$weekday ==6)

hist(log1p(postsData$num.comments))

## control mean
control.num.comments <- mean(postsData$num.comments)
control.newcomer.comments <- mean(postsData$newcomer.comments)

####################################
## UTILITY METHODS                ##
####################################

## randomSample pulls a random number of rows
## from a dataframe up to the number of rows
## in that dataframe
randomSample = function(df,n) { 
  return (df[sample(nrow(df), n),])
}

## LOAD COLORBLIND SAFE PALETTE
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## replicate simulate.study function
set.seed(424)

## function to simulate a study
sample.size <- 1000
simulate.study <- function(i, sample.size, control.mean, treat.mean){
  cat(".")
  start.date <- "2017-08-10"
  posts <- data.frame(
    data = seq(as.Date(start.date), as.Date(start.date) + days(sample.size - 1), by = "day")
  )
  posts$wday <- wday(as.Date(posts$data))
  posts$weekend <- (posts$wday == 7 | posts$wday==1)
  num.obvs = nrow(posts)
  posts$day.num <- seq(1, num.obvs)
  
  # randomizations <- blockrand(n=num.obvs, num.levels = 2, block.sizes = c(4,4), id.prefix='post', block.prefix='block',stratum='post')
  #  print(paste("Poem row count", nrow(poems)))
  #  print(paste("Randomizations row count", nrow(randomizations)))
  
  ## randomly sample conidition A group
  randomizations <- randomSample(posts, sample.size / 2)
  ## which posts were selected for treatment
  conditionA <- row.names(randomizations)
  conditionB <- row.names(posts[row.names(posts) %in% conditionA,])
  posts$condition <- ifelse(row.names(posts) %in% conditionA, "A", "B")
  posts$interactions <- NA
  posts$interactions[posts$condition=="A"] <- rnorm(length(conditionA), control.mean) 
  posts$interactions[posts$condition=="B"] <- rnorm(length(conditionB), treat.mean) 
  
  ## will log-transform for our homework
  m1 <- lm(log1p(interactions) ~ condition, data=posts)
  
  power.sim.intercept <- m1$coefficients['(Intercept)']
  power.sim.treat.effect <- m1$coefficients['conditionB']
  treat.coef = coef(summary(m1))[2,]
  power.sim.stderr <- treat.coef['Std. Error']
  power.sim.pvalue <- treat.coef['Pr(>|t|)']
  power.sim.significant <- as.double(power.sim.pvalue) < 0.05
  data.frame(i=i, power.sim.ctl.fit = power.sim.intercept,power.sim.treat.effect=power.sim.treat.effect,
             power.sim.stderr=power.sim.stderr, power.sim.pvalue=power.sim.pvalue,power.sim.significant=power.sim.significant,
             true.effect = treat.mean - control.mean)
}

# sim <- simulate.study(1, sample.size, control.num.comments, control.num.comments * 1.1)
# sim$power.sim.treat.effect

####################################
## SUMMARY DATA                   ##
####################################

### TIME ELAPSED
max(date(postsData$created.utc)) - min(date(postsData$created.utc))

## POSTS PER DAY
nrow(postsData) / as.numeric(max(date(postsData$created.utc)) - min(date(postsData$created.utc)))

## POSTS REMOVED PER DAY
nrow(subset(postsData, visible==0)) / as.numeric(max(date(postsData$created.utc)) - min(date(postsData$created.utc)))

####################################
## UNIVARIATE ANALYSIS            ##
####################################

## NUMBER OF COMMENTS
summary(posts$num.comments)
hist(posts$num.comments)
hist(log1p(posts$num.comments))

## TODO: NUMBER OF COMMENTS REMOVED

## TODO: NUMBER OF NEWCOMER COMMENTS

## TODO: NUMBER OF NEWCOMER COMMENTS REMOVED

####################################
## BIVARIATE ANALYSIS             ##
####################################
corrplot(cor(posts[c('author.prev.participation', 'author.prev.posts', 'front_page', 
                     'newcomer.comments', 'newcomer.comments.removed', 'num.comments',
                     'num.comments.removed', 'visible')]))


#########################################################
### SIMULATE AND PLOT EXPERIMENTS WITH A SMALL SAMPLE  ##
### WHERE THERE iS A DIFFERENCE BETWEEN MEANS          ##
#########################################################
num.models.small.effect <- 100
mean.ctl.small.effect <- 10
mean.treat.small.effect <- 10.2
num.days <- 100

poem.models.small <- simulate.study(1,num.days,mean.ctl.small.effect,mean.treat.small.effect)
for(i in seq(2,num.models.small.effect)){
  poem.models.small <- rbind(poem.models.small, simulate.study(i,num.days,mean.ctl.small.effect,mean.treat.small.effect))
}

## REPORT HOW MANY RESULTS WERE STATISTICALLY SIGNIFICANT
paste(sprintf("%.01f", as.numeric(summary(poem.models.small$power.sim.significant)[['TRUE']])/nrow(poem.models.small)*100), 
      "% results statistically significant", sep="")

## SHOW STUDIES WHERE THE RESULT IS STATISTICALLY SIGNIFICANT
ggplot(poem.models.small, aes(power.sim.significant, power.sim.treat.effect, color=power.sim.significant)) +
  geom_jitter() +
  theme_bw(base_size = 15, base_family = "Helvetica") +
  ggtitle(paste(sprintf("%.0f", as.numeric(summary(poem.models.small$power.sim.significant)[['TRUE']])/nrow(poem.models.small)*100), 
                "% results statistically significant\nin simulated experiments", sep=""))

##########################################################
### POWER ANALYSIS SHOWING HOW MANY OBSERVATIONS TO GET ##
### AN EIGHTY PERCENT CHANCE OF OBSERVING THE EFFECT    ##
##########################################################

### constants
power.num.models <- 50

power.mean.ctl.comments <- 10
power.mean.treat.comments <- 10.1

power.mean.ctl.newcomer <- 10
power.mean.treat.newcomer <- 10.4

power.max.num.days <- 600
power.starting.num.days <- 12
power.increase.num.days.by <- 12

### power analysis function
power.analysis <- function(i, sample.size, control.mean, treat.mean, num.models){
  print(paste("Sample Size:", sample.size, "Control.mean:", control.mean, "Treat.mean:", treat.mean))
  pm <- simulate.study(1,sample.size,control.mean, treat.mean)
  for(i in seq(2,num.models)){
    pm <- rbind(pm, simulate.study(i,sample.size, control.mean,treat.mean))
  }
  #  print(paste(nrow(pm), " Models calculated"))
  pct.significant <- as.numeric(summary(pm$power.sim.significant)[['TRUE']])/nrow(pm)*100
  mean.effect.significant <- mean(subset(pm, power.sim.significant ==TRUE)$power.sim.treat.effect)
  mean.effect <-  mean(pm$power.sim.treat.effect)
  
  data.frame(i=i, sample.size=sample.size, pct.significant = pct.significant, 
             mean.effect.significant = mean.effect.significant, 
             mean.effect = mean.effect)
}

### run power analysis for number of comments
p.analyses <- power.analysis(1,power.starting.num.days,power.mean.ctl.comments,power.mean.treat.comments,power.num.models)
for(i in seq(power.starting.num.days/power.increase.num.days.by, power.max.num.days/power.increase.num.days.by)){
  sample.size <- i*power.increase.num.days.by
  p.analyses <- rbind(p.analyses, power.analysis(i,i*power.increase.num.days.by,power.mean.ctl.comments,power.mean.treat.comments,power.num.models))
}

### run power analysis for number of newcomer comments
p.analyses <- power.analysis(1,power.starting.num.days,power.mean.ctl.newcomer,power.mean.treat.newcomer,power.num.models)
for(i in seq(power.starting.num.days/power.increase.num.days.by, power.max.num.days/power.increase.num.days.by)){
  sample.size <- i*power.increase.num.days.by
  p.analyses <- rbind(p.analyses, power.analysis(i,i*power.increase.num.days.by,power.mean.ctl.comments,power.mean.treat.comments,power.num.models))
}

### PLOT THE RELATIONSHIP BETWEEN THE SAMPLE SIZE AND STATISTICAL POWER
ggplot(p.analyses, aes(sample.size, pct.significant, color=pct.significant>=80)) +
  geom_point() +
  scale_y_continuous(breaks = round(seq(0,100, by = 10),1)) +
  theme_bw(base_size = 15, base_family = "Helvetica") +
  ggtitle("The larger the sample size, the greater the chance of observing the effect")



#######################################
## EXAMPLE OF TAKING A RANDOM SAMPLE ##
#######################################

## EXAMPLE FOR TAKING A SAMPLE SMALLER THAN THE SOURCE
posts.10k <-randomSample(posts, 10000)

## EXAMPLE FOR TAKING A SAMPLE LARGER THAN THE SOURCE
posts.80k <- randomSample(posts, 40000)
posts.80k <- rbind(posts.80k, randomSample(posts, 40000))

#######################################################
## EXAMPLE OF AN ITERATION IN A POWER ANALYSIS       ##
#######################################################

## SETTINGS
num.observations = 10000

sim.posts <-randomSample(posts, num.observations)

## GENERATE RANDOMIZATIONS
randomizations <- blockrand(n=nrow(sim.posts), num.levels = 2, block.sizes = c(12,12), id.prefix='post', block.prefix='block',stratum='post')
sim.posts$condition <- head(randomizations$treatment, nrow(sim.posts))

## GENERATE AVERAGE TREATMENT EFFECT FOR NUM.COMMENTS

## model log-transformed number of comments
effect.multiplier = 1.5

sim.posts$num.comments.effect <- 1
sim.posts$num.comments.effect[sim.posts$condition=="B"] <- abs(rnorm(nrow(subset(sim.posts, condition=="B")), 
                                                                     effect.multiplier))
sim.posts$num.comments.sim <- sim.posts$num.comments * sim.posts$num.comments.effect

## PLOT RELATIONSHIP BETWEEN SIMULATED NUMBER AND OBSERVED NUMBER
#ggplot(sim.posts, aes(num.comments, num.comments.sim, color=condition)) +
#  geom_jitter()

## PLOT SIMULATED AVERAGE TREATMENT EFFECT
#ggplot(sim.posts, aes(condition, log1p(num.comments.sim), color=condition)) +
#  geom_violin()

## ESTIMATE AVERAGE TREATMENT EFFECT ON LOG-TRANSFORMED VARIABLE
summary(lm(log1p(num.comments.sim) ~ condition, data=sim.posts))


