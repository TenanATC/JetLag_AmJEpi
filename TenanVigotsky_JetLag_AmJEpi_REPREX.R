#This is the reproducible code file for the manuscript entitled
#"A Causal and Replication Analysis of Claims that Jet Lag Affects Team Sport Performance"
#Which was published in the American Journal of Epidemiology
#Preprint Link:
#Journal Article Link:

#The goal is for this code to both fully replicate the analysis but also "show the work"
#on code that was used to produce some of the files contributing to the replication
#For example, creation of the Randomization Inference dataset has a stochastic element
#And also can take a while, so we've commented out the code and instead posted our
#REPREX to the github page for download, the same is true for the weights...
#Though it sounds like the 'independenceweights' package should now be updated
#to be fully reproducible via a set.set of their optimizer, this was not available
#during the analysis and creation of our work.

library(lme4)
library(independenceWeights)
library(dplyr)
library(cobalt)
library(WeightIt)
library(wesanderson)
library(ggplot2)
library(patchwork)
library(mgcv)
library(pbapply)
library(permuco)
library(ggpubr)
library(gratia)




SSAC_dat <- readRDS(gzcon(url("https://github.com/TenanATC/JetLag_AmJEpi/raw/main/Tenan_JetLag_PrimaryData.rds")))

cores <- 8 #if you have a Windows machine and have less than 8 cores, you can just change this


#zero-ing out weather variables when elements are '0' for precip and wind speed
#public reports state that indoor temp is typically set to 70 degrees, so we'll use that for temperature
# https://www.cbsnews.com/minnesota/news/how-do-they-heat-u-s-bank-stadium/
# for the relative humidity, we'll set it at a random sample 50-60. This seems like a practical choice but not one
# where I have been able to find a good evidence base
SSAC_dat$precip <- ifelse(SSAC_dat$elements == 1, SSAC_dat$precip, 0)
SSAC_dat$humidity_high <- ifelse(SSAC_dat$elements == 1, SSAC_dat$humidity_high, sample(x= seq(45,65, by=0.1), 1))
SSAC_dat$temp_high <- ifelse(SSAC_dat$elements == 1, SSAC_dat$temp_high, sample(x= seq(70, 75, by = 0.01), 1))
SSAC_dat$wind_speed <- ifelse(SSAC_dat$elements == 1, SSAC_dat$wind_speed, 0)

SSAC_dat <- data.frame(SSAC_dat)
SSAC_dat$away_team <- as.factor(SSAC_dat$away_team)




#Demonstrating the need for accounting of team-effects in propensity scores
#Below code shows the team-level correlation for the treatments.
#While they're relatively low (0.16 and 0.29) they're definitely high enough to bias an analysis,
#On the other hand, the actual analysis, has very low ICC (0.00024198) so we won't do an MLM for final analysis
#Logit based on ICC method from Hox JJ, Moerbeek M, van de Schoot R (2018). 
#Multilevel Analysis: Techniques and Applications. Taylor and Francis. ISBN 9781138121362.  p. 107

lmer_varcor_time <- VarCorr(lmer(time_hrs_round ~ 1 + (1|away_team), data = SSAC_dat))
time_lmer <- as.data.frame(print(lmer_varcor_time, comp=c('Variance')))
time_lmer$vcov[1] / (time_lmer$vcov[1] + time_lmer$vcov[2])

lmer_varcor_tz <- VarCorr(lmer(away_tz_diff ~ 1 + (1|away_team), data = SSAC_dat))
tz_lmer <- as.data.frame(print(lmer_varcor_tz, comp=c('Variance')))
tz_lmer$vcov[1] / (tz_lmer$vcov[1] + tz_lmer$vcov[2])

glmer_logit <- glmer(spread_beat_away ~ 1 + (1|away_team), data = SSAC_dat, family = binomial)
summ_glmer_logit <- summary(glmer_logit)
summ_glmer_logit$varcor$away_team[[1]] / (summ_glmer_logit$varcor$away_team[[1]] + (pi^2 / 3))


#non-vegas outcome ICCs for informational purposes
glmer_wins <- glmer(away_win ~ 1 + (1|away_team), data = SSAC_dat, family = binomial)
summ_glmer_wins <- summary(glmer_wins)
summ_glmer_wins$varcor$away_team[[1]] / (summ_glmer_wins$varcor$away_team[[1]] + (pi^2 / 3))
#Above ICC corresponds to the Roy & Forest analysis

lmer_varcor_score <- VarCorr(lmer(postgame_spread ~ 1 + (1|away_team), data = SSAC_dat))
score_lmer <- as.data.frame(print(lmer_varcor_score, comp=c('Variance')))
score_lmer$vcov[1] / (score_lmer$vcov[1] + score_lmer$vcov[2])

lmer_varcor_spread <- VarCorr(lmer(spread_away ~ 1 + (1|away_team), data = SSAC_dat))
spread_lmer <- as.data.frame(print(lmer_varcor_spread, comp=c('Variance')))
spread_lmer$vcov[1] / (spread_lmer$vcov[1] + spread_lmer$vcov[2])


# ###################Energy Balancing Weights#####################################
covariates <- model.matrix( ~ 0 + precip +
                              humidity_high +
                              temp_high +
                              wind_speed +
                              grass_type + away_team,
                            data = SSAC_dat)
treatment <- model.matrix( ~ 0 + away_tz_diff + time_hrs_round:away_tz_diff,
                           data = SSAC_dat)
# set.seed(112233, kind = "L'Ecuyer-CMRG")
# (dcows <- independenceWeights::independence_weights(treatment, covariates))
# dcow_weights <- dcows$weights
# ggpubr::gghistogram(dcows$weights, bins = 40)
# 
# saveRDS(dcows, "dcows.rds")

dcows <- readRDS(gzcon(url("https://github.com/TenanATC/JetLag_AmJEpi/raw/main/dcows.rds")))
dcows

dcow_weights <- dcows$weights

#manually assessing the weights with the data
SSAC_dat$dcow_weights <- dcow_weights
percentiles <- quantile(dcow_weights, probs = c(0.05, 0.95))
SSAC_dat_low <- SSAC_dat %>% filter(dcow_weights <= percentiles[1])
prop.table(table(SSAC_dat_low$away_tz))
prop.table(table(SSAC_dat$away_tz))
prop.table(table(SSAC_dat_low$home_tz))
prop.table(table(SSAC_dat$home_tz))
prop.table(table(SSAC_dat_low$away_team))

SSAC_dat_high <- SSAC_dat %>% filter(dcow_weights >= percentiles[2])
prop.table(table(SSAC_dat_high$away_tz))
prop.table(table(SSAC_dat$away_tz))
prop.table(table(SSAC_dat_high$home_tz))
prop.table(table(SSAC_dat$home_tz))

#I don't see any evidence that any teams are being completely over-weighted or under-weighted by the PSW


#assessing these weights with cobalt's bal.tab, the 'standard' way of assessing covariate balance
covariates_noteams <- as.data.frame(covariates) %>% select(!starts_with('away_team'))
dcows_weighit_tz <- as.weightit(dcow_weights, covs=as.data.frame(covariates), treat= treatment[,1])
dcows_weighit_tzhrs <- as.weightit(dcow_weights, covs=as.data.frame(covariates), treat= treatment[,2])
dcows_weighit_tz_noteams <- as.weightit(dcow_weights, covs=covariates_noteams, treat= treatment[,1])
dcows_weighit_tzhrs_noteams <- as.weightit(dcow_weights, covs=covariates_noteams, treat= treatment[,2])
summ_weightit_tz <- summary(dcows_weighit_tz)
summary(dcows_weighit_tzhrs)
plot(summ_weightit_tz)


#evaulating balance manually for reporting
tz_baltab <- as.data.frame(bal.tab(dcows_weighit_tz, stats= c('m', 'cor'), thresholds = c(m=0.1, cor=0.1), un=T)$Balance)
tzhrs_baltab <- as.data.frame(bal.tab(dcows_weighit_tzhrs, stats= c('m', 'cor'), thresholds = c(m=0.1, cor=0.1), un=T)$Balance)

#making love plots for publication
v <- data.frame(old = c("precip", "humidity_high", "temp_high", "wind_speed", 
                        "grass_typeAT", "grass_typeFT", "grass_typeReal"),
                new = c("Precipitation", "Humidity", "High Temperature", 
                        "Wind Speed", "Artificial Turf Stadium", "Field Turf Stadium", "Real Grass Stadium"))
`-.gg` <- function(plot, layer) {
  if (missing(layer)) {
    stop("Cannot use `-.gg()` with a single argument. Did you accidentally put - on a new line?")
  }
  if (!is.ggplot(plot)) {
    stop('Need a plot on the left side')
  }
  plot$layers = c(layer, plot$layers)
  plot
}
(loveplot_plot1 <- love.plot(dcows_weighit_tz_noteams, stats= c('cor'),
                             var.order = 'unadjusted', 
                             size = 2.5,
                             var.names = v, 
                             colors = wes_palette("FantasticFox1",5,"discrete")[3:4], 
                             title = 'Time Zone',
                             themes = theme_classic()) +
    xlim(-.15,.15) -
    geom_rect(xmin = -0.1,
              xmax = 0.1,
              ymin = -1,
              ymax = 8,
              fill = "lightgrey",
              alpha = 0.05) +
    theme(axis.title.x = element_blank()))

(loveplot_plot2 <- love.plot(dcows_weighit_tzhrs_noteams, 
                             stats= c('cor'),
                             var.order = 'unadjusted',
                             size = 2.5,
                             var.names = v, 
                             colors = wes_palette("FantasticFox1",5,"discrete")[3:4],
                             title = 'Time Zone × Kickoff Time',
                             themes = theme_classic()) + 
    xlim(-.15,.15) -
    geom_rect(xmin = -0.1,
              xmax = 0.1,
              ymin = -1,
              ymax = 8,
              fill = "lightgrey",
              alpha = 0.05))
(loveplot_comb <- loveplot_plot1 / loveplot_plot2 + plot_layout(guides = "collect", heights = c(1,1,0.005)) & 
    theme(legend.position = "bottom", 
          legend.title = element_blank()))

# ggsave("loveplot.png", loveplot_plot, width = 3.5, height = 4, units = 'in', dpi = 500, scale = 1.5)


#main causal inference model of interest
tzone_l3_causal <- gam(spread_beat_away ~ s(away_tz_diff, k=5) + ti(away_tz_diff, time_hrs_round, k=7),  
                       data = SSAC_dat, method = 'REML', family = quasibinomial(link = "logit"), weights = dcow_weights) #
anova(tzone_l3_causal) #not valid p-values!
gam.check(tzone_l3_causal) # verifying smooth terms reasonable

#null logit model
tzone_l3_null <- gam(spread_beat_away ~ 1 ,
                     data = SSAC_dat, method = 'REML', family = quasibinomial(link = "logit"), weights = dcow_weights) #, weights = dcow_weights
#get F-stat for incorrect inference, but used to contrast with RI for correct inferences
real_logit_ftest <- anova(tzone_l3_null, tzone_l3_causal, test='F')

#main causal inference to assess a relationship between the Vegas Spread and "jet lag'
#using bam() in this case to go from 60 sec per model to less than 1 sec
vegas_l3_causal <- bam(spread_away ~ s(away_tz_diff, k=5) +  s(away_team, bs='re'), 
                       data = SSAC_dat, method = 'fREML', discrete=FALSE, family= gaussian(), weights = dcow_weights) 
gam.check(vegas_l3_causal, old.style = T)
#null logit model
vegas_l3_null <- bam(spread_away ~ 1 + s(away_team, bs='re'),
                     data = SSAC_dat, method = 'fREML', discrete=FALSE, family= gaussian(), weights = dcow_weights) 
#get F-stat for incorrect inference, but used to contrast with RI for correct inferences
real_vegas_ftest <- anova(vegas_l3_null, vegas_l3_causal, test='F')



#get grid of prediction and then two-tailed pvalues for each estimate
#this is later used to examine the p-value for each individual estimate as
#opposed to the omnibus effect of the main estimand
newdat_grid <- data.frame(expand.grid(away_tz_diff = seq.int(-3, 3), time_hrs_round = seq.int(11,24)))
real_logit_pred <- predict(tzone_l3_causal, newdat_grid, se.fit = T)
real_logit_zscoring <- real_logit_pred$fit/real_logit_pred$se.fit



##Because we're using frequency weights, the Standard Errors are invalid (though the estimates are right)
##Therefore we can't get valid p-values. Bootstrap won't work with splines and
##sampling posterior won't work with weights
##The below code employing randomization inference is commented out and has already been done and saved to github
away_tz <- SSAC_dat$away_tz
rand_int <- SSAC_dat$away_team


#' #Need to take the away team's usual time zone and use that to determine what their
#' #'possible' time changes could have been for RI (positivity assumption)
#' 
#' get_all_tz <- unique(SSAC_dat$away_tz)
#' eastern <- c(3,2,1,0)
#' central <- c(2,1,0,-1)
#' mtn <- c(1,0,-1,-2)
#' pac <- c(0,-1,-2,-3)
#' 
#' team_tz_randomizations <- vector(mode = 'list', length = 5000)
#' team_time_randomizations <- vector(mode = 'list', length = 5000)
#' 
#' time_teams <- data.frame(table(SSAC_dat$time_hrs_round, SSAC_dat$away_team))
# time_teams2<- time_teams %>% group_by(Var2) %>% mutate(frac = Freq/sum(Freq)) %>% arrange(Var1, .by_group = T)
# 
# 
# #This one creates the proportionally random kickoffs
# #It could be written to run faster but we only need to do it once
# for (g in 1:5000){
#   print(g)
#   tz_vec <- vector(mode = 'integer', length = length(away_tz))
#   time_vec <- vector(mode = 'integer', length = length(rand_int))
#   for (t in 1:length(away_tz)) {
#     tz_vec[t] <- ifelse(away_tz[t] == "America/New_York", sample(x= eastern, size = 1, replace = T),
#                         ifelse(away_tz[t] == "America/Kentucky/Louisville", sample(x= eastern, size = 1, replace = T),
#                                ifelse(away_tz[t] == "America/Indiana/Indianapolis", sample(x= eastern, size = 1, replace = T),
#                                       ifelse(away_tz[t] == "America/Detroit", sample(x= eastern, size = 1, replace = T),
#                                              ifelse(away_tz[t] == "America/Chicago", sample(x= central, size = 1, replace = T),
#                                                     ifelse(away_tz[t] == "America/Denver", sample(x= mtn, size = 1, replace = T),
#                                                            ifelse(away_tz[t] == "America/Phoenix", sample(x= mtn, size = 1, replace = T),
#                                                                   ifelse(away_tz[t] == "America/Boise", sample(x= mtn, size = 1, replace = T),
#                                                                          ifelse(away_tz[t] == "America/Los_Angeles", sample(x= pac, size = 1, replace = T),
#                                                                                 99999)))))))))
#     time_vec[t] <- sample(x = seq(11, 24, by=1), size = 1, replace = T,
#                           prob =time_teams2[ which (time_teams2$Var2== rand_int[t]),]$frac) #this one is for keeping proportions
#   }
#   team_tz_randomizations[[g]] <- tz_vec
#   team_time_randomizations[[g]] <- time_vec
# }
# 
# #This one creates fully random kickoff times
# #Note that it will over-write the above dataset if you run them both... so don't or change the names!
# for (g in 1:5000){
#   print(g)
#   tz_vec <- vector(mode = 'integer', length = length(away_tz))
#   time_vec <- vector(mode = 'integer', length = length(rand_int))
#   for (t in 1:length(away_tz)) {
#     tz_vec[t] <- ifelse(away_tz[t] == "America/New_York", sample(x= eastern, size = 1, replace = T),
#                         ifelse(away_tz[t] == "America/Kentucky/Louisville", sample(x= eastern, size = 1, replace = T),
#                                ifelse(away_tz[t] == "America/Indiana/Indianapolis", sample(x= eastern, size = 1, replace = T),
#                                       ifelse(away_tz[t] == "America/Detroit", sample(x= eastern, size = 1, replace = T),
#                                              ifelse(away_tz[t] == "America/Chicago", sample(x= central, size = 1, replace = T),
#                                                     ifelse(away_tz[t] == "America/Denver", sample(x= mtn, size = 1, replace = T),
#                                                            ifelse(away_tz[t] == "America/Phoenix", sample(x= mtn, size = 1, replace = T),
#                                                                   ifelse(away_tz[t] == "America/Boise", sample(x= mtn, size = 1, replace = T),
#                                                                          ifelse(away_tz[t] == "America/Los_Angeles", sample(x= pac, size = 1, replace = T),
#                                                                                 99999)))))))))
#     time_vec[t] <- sample(x = seq(11, 24, by=1), size = 1, replace = T) #for completely random times
#   }
#   team_tz_randomizations[[g]] <- tz_vec
#   team_time_randomizations[[g]] <- time_vec
# }


team_tz_randomizations_clust <- readRDS(gzcon(url("https://github.com/TenanATC/JetLag_AmJEpi/raw/main/team_tz_randomizations_clust.rds")))
team_time_randomizations_clust <- readRDS(gzcon(url("https://github.com/TenanATC/JetLag_AmJEpi/raw/main/team_time_randomizations_clust.rds")))




##While many Randomization Inference studies use the model coefficient, that does not work for splines
##What makes the most sense in our case is to use the 'drop 1' F-Statistic for Randomization Inference
##In this case, while we're using this particular type of interaction spline with 'mgcv'
##and we're making the noise data with those interaction variables, the model-wide metric should be valid
##This slightly changes our interpretation, so we're asking if time zone change and game time impact performance,
##Not a question of those individual main effects (though we could run those analyses separately)

##Vegas Spread of Interest
outcome_logit <- SSAC_dat$spread_beat_away
outcome_vegas <- SSAC_dat$spread_away


#F Statistic for Logit
effect_logit <- real_logit_ftest$F[[2]]


# build null distribution over a bunch of possible randomizations with kickoff clustered data
cl <- parallel::makePSOCKcluster(cores)
parallel::clusterExport(cl, c('team_tz_randomizations_clust', 'team_time_randomizations_clust', 'outcome_logit')) 

nulls_clust <- pblapply(1:length(team_tz_randomizations_clust), function(i) {
  library(mgcv)
  data_model <- data.frame(temp1 = team_tz_randomizations_clust[[i]], temp2 = team_time_randomizations_clust[[i]], 
                           outcome_logit)
  try(model <- gam(outcome_logit ~ s(temp1, k=5) + ti(temp1, temp2, k=7), data = data_model, 
                   family = quasibinomial(link = "logit"), method = 'REML')) 
  try(null_m <- gam(outcome_logit ~ 1, data = data_model, 
                    family = quasibinomial(link = "logit"), method = 'REML'))
  try(anova(null_m, model, test='F')$F[[2]])
  
},cl = cl)
parallel::stopCluster(cl)

ri_logit_clust <- unlist(nulls_clust) 
# plot the null and observed effect
hist(ri_logit_clust, breaks = 'fd')
abline(v = effect_logit)
# p-value
sum(abs(effect_logit) < abs(ri_logit_clust))/length(ri_logit_clust) # randomization test


#use RI to look at the statistical "significance" of individual point estimates
cl <- parallel::makePSOCKcluster(cores)
parallel::clusterExport(cl, c('team_tz_randomizations_clust', 'team_time_randomizations_clust', 'outcome_logit')) 

nulls_preds <- pblapply(1:length(team_tz_randomizations_clust), function(i) {
  library(mgcv)
  data_model <- data.frame(temp1 = team_tz_randomizations_clust[[i]], temp2 = team_time_randomizations_clust[[i]], 
                           outcome_logit)
  try(model <- gam(outcome_logit ~ s(temp1, k=5) + ti(temp1, temp2, k=7), data = data_model, 
                   family = quasibinomial(link = "logit"), method = 'REML')) 
  newdat_grid <- data.frame(expand.grid(temp1 = seq.int(-3, 3), temp2 = seq.int(11,24)))
  try(predx <- predict(model, newdat_grid, se.fit=T))
  #pred_pval <- 2*pnorm(abs(predx$fit/predx$se.fit), lower.tail = F)
  pred_zscore <- predx$fit/predx$se.fit
  
},cl = cl)
parallel::stopCluster(cl)

nulls_preds_z <- unname(as.matrix(cbind.data.frame(nulls_preds)))


#now implementing maxT correction from permuco library
real_logit_zscoring_t <- t(real_logit_zscoring)
nulls_preds_z_t <- t(nulls_preds_z)
mat_maxt <- rbind(real_logit_zscoring_t, nulls_preds_z_t)
maxt_p <- compute_stepdownmaxT(mat_maxt, 'two.sided')
newdat_grid <- data.frame(expand.grid(temp1 = seq.int(-3, 3), temp2 = seq.int(11,24)))
point_estimates_pvalue <- cbind.data.frame(newdat_grid, maxt_p)



#F Statistic for Vegas spread relationship to 'jet lag'
effect_vegas <- real_vegas_ftest$F[[2]]

# build null distribution over a bunch of possible randomizations with kickof clustered data, non-trimmed dataset
cl <- parallel::makePSOCKcluster(cores)
parallel::clusterExport(cl, c('team_tz_randomizations_clust', 'team_time_randomizations_clust', 'outcome_vegas', 'rand_int')) 

nulls_clust_vegas <- pblapply(1:length(team_tz_randomizations_clust), function(i) {
  library(mgcv)
  data_model <- data.frame(temp1 = team_tz_randomizations_clust[[i]], temp2 = team_time_randomizations_clust[[i]], 
                           outcome_vegas, rand_int)
  try(model <- bam(outcome_vegas ~ s(temp1, k=5) + s(rand_int, bs= 're'), data = data_model, #+ ti(temp1, temp2, k=7) 
                   family = gaussian(), method = 'fREML', discrete=FALSE)) 
  try(null_m <- bam(outcome_vegas ~ 1 + s(rand_int, bs= 're'), data = data_model, 
                    family = gaussian(), method = 'fREML', discrete=FALSE))
  try(anova(null_m, model, test='F')$F[[2]])
  
},cl = cl)
parallel::stopCluster(cl)

ri_vegas_clust <- unlist(nulls_clust_vegas) 
# plot the null and observed effect
hist(ri_vegas_clust, breaks = 'fd')
abline(v = effect_vegas, lwd=2)
# p-value
sum(abs(effect_vegas) < abs(na.omit(ri_vegas_clust)))/length(na.omit(ri_vegas_clust)) # randomization test
#That's awesome. No evidence that there is a relation between Vegas spread and our 'jet lag' data


#############################################################################################
####Making Figures of above analyses
#############################################################################################
#Making pretty ggplot of Logistic RI histogram
ri_logit_clust_df <- data.frame('F_statistic' = ri_logit_clust)

(distr <- ggplot() + 
    geom_histogram(data = ri_logit_clust_df, 
                   mapping = aes(x=F_statistic),
                   bins = nclass.FD(ri_logit_clust_df$F_statistic), 
                   fill='darkgray') +
    geom_vline(xintercept = effect_logit,
               color = wes_palette("FantasticFox1",5,"discrete")[4], 
               linewidth = 1.5) +
    annotate("text", x=3.5, y= 230, 
             label= "Observed",
             color = wes_palette("FantasticFox1",5,"discrete")[4]) +
    coord_cartesian(expand = F) +
    xlab("F-statistic") +
    theme_classic() +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
)



#Now the surface plot of logistic regression point estimates 
ds1 <- gratia::data_slice(tzone_l3_causal, away_tz_diff = evenly(away_tz_diff, 50), time_hrs_round = evenly(time_hrs_round, 50))
fv1 <- gratia::fitted_values(tzone_l3_causal, ds1, scale = "response")
pal <- wes_palette("Zissou1", 100, "continuous")
(surface_plot <-
    ggplot(data = fv1,
           aes(x = away_tz_diff,
               y = time_hrs_round,
               z = .fitted,
               fill = .fitted)) +
    ggrastr::geom_tile_rast() +
    metR::geom_contour2(aes(label = stat(level))) +
    coord_cartesian(xlim=c(-3,3),
                    ylim=c(11,24),
                    expand = F) +
    scale_y_continuous(breaks = seq(12,24,by=4)) +
    scale_fill_gradientn(colors = pal, lim = c(0,1), guide = guide_colorbar(theme = theme(legend.key.width = unit(0.75,"lines")))) +
    labs(fill = "Probability",
         x = "Hours lost/gained by away team travel",
         y = "Kickoff time (Eastern Time)")
)

design <- "
AB
AC
"
ci_results <-
  wrap_elements(full = loveplot_comb) + surface_plot + distr +
  plot_layout(ncol = 2, widths = c(1.2,1), heights = c(1,0.1), 
              design = design) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = 'bold'))

#slightly different layout from final manuscript, but same figures/results
ci_results




##Creating Supplement 2 figure for the histogram of relationship between pre-game spread and jet lag
#removing the few 'NA's that may occur due to model non-convergences with what is,
#essentially, a garbage model
ri_vegas_clust_df <- na.omit(data.frame('F_statistic' = ri_vegas_clust)) 

ri_vegasclust_plot <- ggplot(ri_vegas_clust_df, aes(x=F_statistic)) + geom_histogram(bins = 50, fill='steelblue1', color='black') +
  geom_vline(xintercept = effect_vegas, linewidth=1.2) + 
  geom_curve(aes(x=effect_vegas, y= 310, xend= 4.5, yend= 250), curvature = -0.5, linewidth=1, 
             arrow = arrow(length = unit(0.1, 'inches'), ends = 'first', type = 'closed')) +
  annotate("text", x= 5, y= 230, label= "F-Statistic from Observed Data\np = 0.142") +
  ggtitle('Randomization Inference for Effect Mediation on\nBeating the Las Vegas Spread') + 
  ylab('Count') + xlab('F-Statistic Distribution for All Potential Outcomes') + 
  theme_pubclean() + theme(plot.title = element_text(hjust = 0.5, size = 20),
                           axis.text=element_text(size=14),
                           axis.title=element_text(size=14))

(distr_vegas <- ggplot() + 
    geom_histogram(data = ri_vegas_clust_df, 
                   mapping = aes(x=F_statistic),
                   bins = nclass.FD(ri_vegas_clust_df$F_statistic), 
                   fill='darkgray') +
    geom_vline(xintercept = effect_vegas,
               color = wes_palette("FantasticFox1",5,"discrete")[4], 
               linewidth = 1.5) +
    annotate("text", x=3.5, y= 230, 
             label= "Observed",
             color = wes_palette("FantasticFox1",5,"discrete")[4]) +
    coord_cartesian(expand = F) +
    xlab("F-statistic") +
    theme_classic() +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
)




####Trying to re-create the analysis from Smith et al 1999
###"A variable x was defined based on the score of the west coast team minus the score of the east coast team 
###plus the point spread. Thus, the value of x would be positive if the west coast team
###beat the point spread and negative if the west coast team did not beat the point spread. One sample t-tests were performed to
###assess whether the x values were significantly greater than zero."

#overall
smith_west <- SSAC_dat %>% filter(away_tz_diff == -3) %>% mutate(smith_x = (away_score - home_score) - spread_away) %>%
  dplyr::select(home_team.y, away_team.y, home_score, away_score, smith_x, time_hrs_round)
smith_east <- SSAC_dat %>% filter(away_tz_diff == 3) %>% mutate(smith_x = (home_score - away_score) + spread_away) %>%
  dplyr::select(home_team.y, away_team.y, home_score, away_score, smith_x, time_hrs_round)
smith_dat <- rbind.data.frame(smith_east, smith_west)
t.test(smith_dat$smith_x)
t.test(smith_dat$smith_x)$stderr

#just evening games
smith_west_evening <- SSAC_dat %>% filter(away_tz_diff == -3) %>% mutate(smith_x = (away_score - home_score) - spread_away) %>%
  filter(time_hrs_round >=20) %>%  dplyr::select(home_team.y, away_team.y, home_score, away_score, smith_x)
smith_east_evening <- SSAC_dat %>% filter(away_tz_diff == 3) %>% mutate(smith_x = (home_score - away_score) + spread_away) %>%
  filter(time_hrs_round >=20) %>% dplyr::select(home_team.y, away_team.y, home_score, away_score, smith_x)
smith_evening_dat <- rbind.data.frame(smith_east_evening, smith_west_evening)
t.test(smith_evening_dat$smith_x)
t.test(smith_evening_dat$smith_x)$stderr
#non-significant, so we definitely can't replicate their results at all
#figure not in manuscript but just depicting the results graphically
smith_hist <- ggplot(smith_dat, aes(x=smith_x)) + geom_histogram(bins = 20, fill='steelblue1', color='black') +
  geom_vline(xintercept = 0, linewidth=2) + 
  ggtitle("Inability to Replicate Smith et al.'s Findings") + annotate("text", x=19, y=9.5, label= "P = 0.952", size=5) + 
  ylab('Count') + xlab("Number of Points West Coast Team Did or\nDid Not Beat The Las Vegas Point Spread By") + 
  theme_pubclean() + theme(plot.title = element_text(hjust = 0.5, size = 20),
                           axis.text=element_text(size=14),
                           axis.title=element_text(size=14))

smith_hist


####Now trying to re-create the Roy and Forest NFL analysis
### "The outcomes for each game were: game time (afternoon/17:00 hours Eastern Standard Time or earlier;
### evening/17:30 hours Eastern Standard Time or later); the result of the game 
### (then the winning percentage was calculated for all games); the number of time zones travelled;
### and the direction of travelling (westward, same time zone or eastward)."  ...
### Our work focused on replicating their Figure 1 for their simple linear regression

royforest_dat <- SSAC_dat %>% select(season.x, away_team, away_tz, away_win, time_hrs, away_tz_diff)
royforest_dat_regression <- royforest_dat %>% filter(time_hrs >17.5) %>% 
  mutate(timezone_royforest= away_tz_diff*-1)

#simple regression analysis, hopefully we can all agree that "linear probability models" are dumb
royforest_reg <- lm(away_win ~ timezone_royforest, data = royforest_dat_regression)
summary(royforest_reg)
cor.test(royforest_dat_regression$timezone_royforest, royforest_dat_regression$away_win)



