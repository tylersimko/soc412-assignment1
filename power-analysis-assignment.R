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

##########################################################
### POWER ANALYSIS SHOWING HOW MANY OBSERVATIONS TO GET ##
### AN EIGHTY PERCENT CHANCE OF OBSERVING THE EFFECT    ##
##########################################################

### constants
power.num.models <- 50

power.max.num.days <- 150
power.starting.num.days <- 3
power.increase.num.days.by <- 3

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

### number of comments
# power.mean.ctl.comments <- 10
# effect.multiplier.comments <- 1.1
# power.mean.treat.comments <- power.mean.ctl.comments * effect.multiplier.comments
# 
# p.analyses <- power.analysis(1,power.starting.num.days,power.mean.ctl.comments,power.mean.treat.comments,power.num.models)
# for(i in seq(power.starting.num.days/power.increase.num.days.by, power.max.num.days/power.increase.num.days.by)){
#   sample.size <- i*power.increase.num.days.by
#   p.analyses <- rbind(p.analyses, power.analysis(i,i*power.increase.num.days.by,power.mean.ctl.comments,power.mean.treat.comments,power.num.models))
# }
# 
# ggplot(p.analyses, aes(sample.size, pct.significant, color=pct.significant>=80)) +
#   geom_point() +
#   scale_y_continuous(breaks = round(seq(0,100, by = 10),1)) +
#   theme_bw(base_size = 15, base_family = "Helvetica") +
#   ggtitle("Total number of comments")

# ### number of newcomer comments
power.mean.ctl.newcomer <- 10
effect.multiplier.newcomer <- 1.38
power.mean.treat.newcomer <- power.mean.ctl.newcomer * effect.multiplier.newcomer

p.analyses <- power.analysis(1,power.starting.num.days,power.mean.ctl.newcomer,power.mean.treat.newcomer,power.num.models)
for(i in seq(power.starting.num.days/power.increase.num.days.by, power.max.num.days/power.increase.num.days.by)){
  sample.size <- i*power.increase.num.days.by
  p.analyses <- rbind(p.analyses, power.analysis(i,i*power.increase.num.days.by,power.mean.ctl.comments,power.mean.treat.comments,power.num.models))
}

ggplot(p.analyses, aes(sample.size, pct.significant, color=pct.significant>=80)) +
  geom_point() +
  scale_y_continuous(breaks = round(seq(0,100, by = 10),1)) +
  theme_bw(base_size = 15, base_family = "Helvetica") +
  ggtitle("Total number of newcomer comments")