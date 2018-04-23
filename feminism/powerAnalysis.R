## LOAD LIBRARIES
library(ggplot2)
library(corrplot)
library(gmodels)
library(texreg)
library(rms)
library(blockrand)
library(lubridate)
library(psych)
library(RcppZiggurat)
library(MASS)

## CLEAR PROJECT
rm(list=ls())

####################################
## UTILITY METHODS                ##
####################################

## randomSample pulls a random number of rows
## from a dataframe up to the number of rows
## in that dataframe
randomSample = function(df,n) { 
  return (df[sample(nrow(df), n,replace=TRUE),])
  #return(sample(df, n, replace=TRUE))
}

## LOAD COLORBLIND SAFE PALETTE
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

####################################
## LOAD DATAFRAMES                ##
####################################
setwd("~/Documents/soc412-assignment1/")
fem <- read.csv("feminism/feminism_comments_on_posts_06.30.2017-12.31.2107.csv", stringsAsFactors = FALSE)

firstTimeNames <- fem$author[fem$previous.comments == 0]
firstTimeData <- fem[fem$author %in% firstTimeNames,]

## how many times did a first-time commenter comment again in 2 weeks
nums <- NULL

## for each first-time commenter, how many post a second time in two weeks? / 3 months?
for (i in 1:length(firstTimeNames)) {
  if (i %% 100 == 0) print(i)
  
  currentAuthor <- firstTimeNames[i]
  authorData <- firstTimeData[firstTimeData$author == currentAuthor,]
  ## transform times into UTC
  times <- as.POSIXct(authorData$created_utc, origin = "1970-01-01")
  earliestPost <- authorData[which(times == min(times)),]
  otherPosts <- subset(authorData, id != earliestPost$id)

  ## find earliest date
  earliestTime <- as.POSIXct(earliestPost$created_utc, origin = "1970-01-01")

  ## get posts within two weeks of this date
  twoWeekPosts <- which(as.POSIXct(otherPosts$created_utc, origin = "1970-01-01") >= earliestTime & 
        as.POSIXct(otherPosts$created_utc, origin = "1970-01-01") <= (as.Date(earliestTime) + 14))
  twoWeekPosts <- otherPosts[twoWeekPosts,]

  ## merge data
  currentAuthorData <- c(firstTimeNames[i], nrow(twoWeekPosts))
  nums <- rbind(nums, currentAuthorData)
}

row.names(nums) <- NULL
nums <- as.data.frame(nums)
colnames(nums) <- c("author", "count")
nums$author <- as.character(nums$author)
nums$num.comments <- as.numeric(as.character(nums$count))

##################################################
## MODEL NEGATIVE BINOMIAL DISTRIBUTIONS OF DVs ##
## (models to be used for power analysis)       ##
##################################################

base.num.comments                <- glm.nb(num.comments ~ 1, data=nums)
base.num.comments.log.lm         <- lm(log1p(num.comments) ~ 1, data=nums)
base.num.comments.p              <- glm(num.comments ~ 1, data=nums, family="poisson")

#######################################################
## EXAMPLE OF AN ITERATION IN A POWER ANALYSIS       ##
#######################################################

simulate.study.from.negbin <- function(i, sample.size, model.nb, pct.diff){
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
  m.nb.significant <- as.double(m.nb.pvalue) < 0.05  
  
  
  #m1 <- lm(log1p(dv) ~ TREAT, data=sim.obs)
  # m1.intercept <- m1$coefficients[['(Intercept)']]
  # m1.treat.effect <- m1$coefficients[['TREAT']]
  # m1.treat.coef = coef(summary(m1))[2,]
  # m1.stderr <- m1.treat.coef['Std. Error']
  # m1.pvalue <- m1.treat.coef['Pr(>|t|)']
  # m1.significant <- as.double(m1.pvalue) < 0.05
  data.frame(i=i,
             power.sim.comments.intercept      = m.nb.intercept, 
             power.sim.comments.treat.effect   = m.nb.treat.effect,
             power.sim.comments.stderr         = m.nb.stderr,
             power.sim.comments.pvalue         = m.nb.pvalue,
             power.sim.comments.significant    = m.nb.significant)
}


simulate.study <- function(i, sample.size, effect.size){
  cat(".")
  #print(paste("sample", sample.size, ". Effect size", effect.size, "."))
  
  sim.posts <-randomSample(posts, sample.size)
  sim.posts$wday <- wday(sim.posts$created.utc)
  sim.posts$weekend <- (sim.posts$wday == 5 | sim.posts$wday==6)
  num.obvs = nrow(sim.posts)
  sim.posts$post.num <- seq(1, num.obvs)
  time.end <- Sys.time()
  
  sim.posts$condition <- as.numeric(zrnorm(num.obvs)>0)
  sim.posts$num.comments.effect[sim.posts$condition==0] <- 1
  
  sim.posts$num.comments.effect[sim.posts$condition==1] <- abs(rnorm(nrow(subset(sim.posts, condition==1)), 0,.2) +
                                                                 effect.size)
  sim.posts$num.comments.sim <- round(sim.posts$num.comments * sim.posts$num.comments.effect)
  #ggplot(sim.posts, aes(factor(condition), log1p(num.comments.sim))) + geom_violin()
  #ggplot(sim.posts, aes(factor(condition), log1p(num.comments))) + geom_violin()
  
  #ggplot(sim.posts, aes(log1p(num.comments), log1p(num.comments.sim), color=factor(condition))) + geom_point()
  
  #sim.posts$num.comments.effect <- comments.effect.size
  #sim.posts$num.comments.sim <- ifelse(sim.posts$condition !=0, sim.posts$num.comments * sim.posts$num.comments.effect, sim.posts$num.comments) 
  #ggplot(sim.posts, aes(condition, log1p(newcomer.comments.sim), color=factor(condition)))+ geom_violin()
  
  m1 <- lm(log1p(num.comments.sim) ~ condition, data=sim.posts)
  power.sim.m1.intercept <- m1$coefficients['(Intercept)']
  power.sim.m1.treat.effect <- m1$coefficients['condition']
  power.sim.m1.treat.coef = coef(summary(m1))[2,]
  power.sim.m1.stderr <- power.sim.m1.treat.coef['Std. Error']
  power.sim.m1.pvalue <- power.sim.m1.treat.coef['Pr(>|t|)']
  power.sim.m1.significant <- as.double(power.sim.m1.pvalue) < 0.05
  
  m.nb <- glm.nb(num.comments.sim ~ condition, data=sim.posts)
  m.nb.intercept <- m.nb$coefficients[['(Intercept)']]
  m.nb.treat.effect <- m.nb$coefficients[['condition']]
  m.nb.treat.coef = coef(summary(m.nb))[2,]
  m.nb.stderr <- m.nb.treat.coef['Std. Error']
  m.nb.pvalue <- m.nb.treat.coef['Pr(>|z|)']
  m.nb.significant <- as.double(m.nb.pvalue) < 0.05  
  
  
  #print(paste("num comments significant", power.sim.m1.significant,"."))
  
  data.frame(i=i,
             power.sim.comments.intercept      = power.sim.m1.intercept, 
             power.sim.comments.treat.effect   = power.sim.m1.treat.effect,
             power.sim.comments.stderr         = power.sim.m1.stderr,
             power.sim.comments.pvalue         = power.sim.m1.pvalue,
             power.sim.comments.significant    = power.sim.m1.significant,
             
             m.nb.intercept                    = m.nb.intercept,
             m.nb.treat.effect                 = m.nb.treat.effect,
             m.nb.stderr                       = m.nb.stderr,
             m.nb.pvalue                       = m.nb.pvalue,
             m.nb.significant                  = m.nb.significant,
             comments.true.effect = effect.size)
}



power.analysis <- function(i, sample.size,  comments.effect.size, num.models){
  pm <- simulate.study(1,sample.size,comments.effect.size)
  for(i in seq(2,num.models)){
    pm <- rbind(pm, simulate.study(i,sample.size,comments.effect.size))
  }
  
  pct.comments.significant <- sum(pm$power.sim.comments.significant)/nrow(pm)*100
  mean.comments.effect <- mean(pm$power.sim.comments.treat.effect)
  
  mean.comments.nb.effect <- mean(pm$m.nb.treat.effect)
  pct.comments.nb.significant <- sum(pm$m.nb.significant)/nrow(pm)*100
  
  data.frame(i=i, sample.size=sample.size, 
             pct.comments.significant = pct.comments.significant,
             mean.comments.effect = mean.comments.effect,
             pct.comments.nb.significant = pct.comments.nb.significant,
             mean.comments.nb.effect = mean.comments.nb.effect)
}


power.analysis.nb <- function(i, sample.size, model.nb, pct.diff, num.models){
  pm <- simulate.study.from.negbin(1,sample.size,model.nb, pct.diff)
  for(i in seq(2,num.models)){
    pm <- rbind(pm, simulate.study.from.negbin(i,sample.size,model.nb, pct.diff))
  }
  
  pct.significant <- sum(pm$power.sim.comments.significant)/nrow(pm)*100
  mean.effect <- mean(pm$power.sim.comments.treat.effect)
  data.frame(i=i, sample.size=sample.size, 
             pct.significant = pct.significant, 
             mean.effect = mean.effect)
}

## SIMULATE ONE STUDY USING THE NEGATIVE BINOMIAL MODELING METHOD
x <- power.analysis.nb(1,1000, base.num.comments, 1.38, 100)

##################################################################
###### POWER ANALYSIS USING THE NEGATIVE BINOMIAL MODELING METHOD
##################################################################
maximum.sample.size <- 2000
minimum.sample.size <- 100
iterate.by          <- 20
sim.effect          <- 1.10
num.models          <- 100

models.decision.nb <- power.analysis.nb(1,minimum.sample.size, base.num.comments, sim.effect, num.models)
for(i in seq(1,maximum.sample.size/iterate.by)){
  print(paste("sample.size",minimum.sample.size + iterate.by*i))
  models.decision.nb <- rbind(models.decision.nb, power.analysis.nb(i,
                                                                    minimum.sample.size + iterate.by*i, base.num.comments,
                                                                    sim.effect, num.models))
}

ggplot(models.decision.nb, aes(sample.size, pct.significant)) +
  geom_smooth() +
  geom_point ()

######################################
### LOOK FURTHER AT SAMPLES BELOW 5000
maximum.sample.size <- 5000
minimum.sample.size <- 500
iterate.by          <- 500
sim.effect          <- 1.38
num.models          <- 100

for(i in seq(1,maximum.sample.size/iterate.by)){
  print(paste("sample.size",minimum.sample.size + iterate.by*i))
  models.decision.nb <- rbind(models.decision.nb, power.analysis.nb(i,
                                                                    minimum.sample.size + iterate.by*i, base.num.comments.weekend,
                                                                    sim.effect, num.models))
}

for(i in seq(1,maximum.sample.size/iterate.by)){
  print(paste("sample.size",minimum.sample.size + iterate.by*i))
  models.decision.mult <- rbind(models.decision.mult, power.analysis(i,
                                                                     minimum.sample.size + iterate.by*i,
                                                                     sim.effect, num.models))
  
}

######################################
### PLOT ALL POWER ANALYSIS RESULTS

ggplot(models.decision.mult, aes(sample.size, pct.comments.significant)) +
  geom_point (aes(color="lm sampled")) +
  geom_smooth(aes(color="lm sampled")) +
  geom_point(aes(x=models.decision.nb$sample.size, y=models.decision.nb$pct.significant, color="nb sim")) +
  geom_smooth(aes(x=models.decision.nb$sample.size, y=models.decision.nb$pct.significant, color="nb sim")) +
  geom_point(aes(x=models.decision.mult$sample.size, y=models.decision.mult$pct.comments.nb.significant, color="nb sampled")) +
  #  geom_smooth(aes(x=models.decision.mult$sample.size, y=models.decision.mult$pct.comments.nb.significant, color="nb sampled")) +
  ylim(0,106) +
  scale_y_continuous(breaks=c(seq(0,100, by=10))) +
  scale_color_manual(values=cbbPalette) +
  scale_x_continuous(breaks=c(seq(0,50000, by=10000))) +
  theme_bw(base_size = 15, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Statistical power for different sampling & modeling approaches\n(true effect: 38% increase in number of comments)")