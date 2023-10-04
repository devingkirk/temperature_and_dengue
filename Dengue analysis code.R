# Load packages that may be required #
require(dplyr)
require(lmerTest)
require(tidyverse)
require(effects)
require(multcomp)
require(psych)
require(gridExtra)
library(MuMIn)
require(mgcv)
require(scales)
require(corrplot)


# Update directory before bringing in data #

# Bring in dengue and temperature data #
master_data <- read.csv("~/Dengue_database.csv") %>%
  filter(coef_or_cor == "correlation")

# Bring in R0 thermal performance curve to use in analysis #
dengueR0 <- read.csv("~/dengueR0.csv") %>%
  mutate(dR0 = lead(R0) - lag(R0),
         dT = lead(temp) - lag(temp),
         dR0_dT = dR0/dT)

# Join it with master_data #
master_data <- master_data %>%
  unique %>% 
  mutate(temp = round(mean_temp_C, 1)) %>% 
  left_join(dengueR0 %>% dplyr::select(temp, dR0_dT)) %>% 
  dplyr::select(-temp)


### PART 1: Temperature-specific analyses ###

# Subset data to only include observations used for temperature analysis #
subsetted_data <- master_data %>%
  filter(temperature_specific_analysis == 1) %>%
  mutate(broad_temperature_metric = as.factor(broad_temperature_metric )) %>%
  mutate(effect_type = as.factor(effect_type)) %>%
  mutate(study_code = as.factor(study_code)) %>%
  mutate(sqrt_sample_size = sqrt(sample_size)) %>%
  mutate(factor_temporal_scale = as.factor(data_temporal_scale))



# Set REML to false (and therefore fit using ML) so that we can do model selection #

# What we will call the true null model with no fixed effects #
null_model <- lmer(effect_size ~ (1|study_code), data=subsetted_data, REML = FALSE)

# Run null model (no mean temperature effect model) but does have broad temperature metric and random effect of study_code #
basic_model <- lmer(effect_size ~ broad_temperature_metric + (1|study_code), data=subsetted_data, REML = FALSE)

# Linear effect of temp #
linear_model  <- lmer(effect_size ~ broad_temperature_metric + mean_temp_C + (1|study_code), data=subsetted_data,REML = FALSE)

# Run model with additional quadratic effect of mean temp #
quadratic_model <- lmer(effect_size ~ broad_temperature_metric + mean_temp_C +  I(mean_temp_C^2)+
                          (1|study_code), data=subsetted_data,REML = FALSE)

# Model with dR0/dT instead of mean_temp and mean_temp^2
dR0dT_model <- lmer(effect_size ~ broad_temperature_metric + dR0_dT + 
                      (1|study_code), data=subsetted_data,REML = FALSE)


# Model comparison #
null.AIC <- AIC(null_model)
basic.AIC <- AIC(basic_model)
linear.AIC <- AIC(linear_model)
quadratic.AIC <- AIC(quadratic_model)
dR0dT.AIC <- AIC(dR0dT_model)

# Look at delta AIC between alternative models and null #
null.AIC - quadratic.AIC
null.AIC - dR0dT.AIC
null.AIC - basic.AIC
null.AIC - linear.AIC

# Pseudo-r2 #
rsq.basic <- r.squaredLR(basic_model, null=null_model)[1]
rsq.linear <- r.squaredLR(linear_model, null=null_model)[1]
rsq.quadratic <- r.squaredLR(quadratic_model, null=null_model)[1]
rsq.dR0dt <- r.squaredLR(dR0dT_model, null=null_model)[1]
rsqs <- c(rsq.basic,rsq.linear,rsq.quadratic,rsq.dR0dt);rsqs

# AICc #
null.AICc <- AICc(null_model)
basic.AICc <- AICc(basic_model)
linear.AICc <- AICc(linear_model)
quadratic.AICc <- AICc(quadratic_model)
dR0dT.AICc <- AICc(dR0dT_model)

null.AICc - quadratic.AIC
null.AICc - dR0dT.AIC
null.AICc - basic.AIC
null.AICc - linear.AIC



# For supplement, repeat analyses while weighting by square root of sampling size #
null_model_weighted <- lmer(effect_size ~ (1|study_code), data=subsetted_data, 
                             weights=sqrt_sample_size,
                             REML = FALSE)


# Run basic model with only broad temperature metric and random effect of study_code #
basic_model_weighted <- lmer(effect_size ~ broad_temperature_metric + (1|study_code), data=subsetted_data, 
                            weights=sqrt_sample_size,
                            REML = FALSE)

# Linear effect of temp #
linear_model_weighted <- lmer(effect_size ~ broad_temperature_metric + mean_temp_C + (1|study_code), data=subsetted_data, 
                              weights=sqrt_sample_size,REML = FALSE)

# Run model with additional quadratic effect of mean temp #
quadratic_model_weighted <- lmer(effect_size ~ broad_temperature_metric + mean_temp_C +  I(mean_temp_C^2)+
                          (1|study_code), data=subsetted_data, weights=sqrt_sample_size,REML = FALSE)


# Model with dR0/dT instead of mean_temp and mean_temp^2
dR0dT_model_weighted <- lmer(effect_size ~ broad_temperature_metric + dR0_dT + 
                          (1|study_code), data=subsetted_data, weights=sqrt_sample_size,REML = FALSE)


# Model comparison #
AIC(null_model_weighted)
AIC(basic_model_weighted)
AIC(linear_model_weighted)
AIC(quadratic_model_weighted)
AIC(dR0dT_model_weighted)


# Check for significant differences in temp predictors #
summary(glht(quadratic_model, linfct = mcp(broad_temperature_metric = "Tukey")), test = adjusted("holm"))


# Check confidence interval around quadratic peak temp #
car::deltaMethod(object=quadratic_model,g. = "-b3/(2*b4)", vcov. = vcov(quadratic_model), parameterNames = c("b0", "b1", "b2", "b3", "b4"))


# Now plot with residuals #

# For confidence intervals, use the effects package --> effect function. term= the fixed effect you want to get data on, mod= name of your model.
effects_quadratic <- effects::effect(term= "mean_temp_C", mod= quadratic_model,xlevels=100)

# Save the effects values as a df:
x_quadratic <- as.data.frame(effects_quadratic)

# Get residuals #
test <- data.frame(resid = residuals(quadratic_model), 
                   est_temp_effect = (effect("mean_temp_C", mod = quadratic_model, 
                                             xlevels = list(mean_temp_C = subsetted_data$mean_temp_C)) %>% 
                                        magrittr::extract2("fit")),
                   temp = subsetted_data$mean_temp_C) %>% 
  mutate(partial = resid + est_temp_effect)

# Plot with partial residuals #
quad_plot <- ggplot() + theme_classic(base_size=18) + theme(legend.position="none",axis.text=element_text(size=20),
                                                            axis.title.x = element_text(margin=margin(t=20,r=0,b=10,l=0),size=24),
                                                            axis.title.y= element_text(size=24)) +
  geom_ribbon(data=x_quadratic, aes(x=mean_temp_C, ymin=lower, ymax=upper), alpha= 0.2, fill="black") +
  geom_line(data=x_quadratic, aes(x=mean_temp_C, y=fit), colour="black",size=2)  +
  geom_point(data=test, aes(temp, partial), size=4 ,shape=21,stroke=1.5,fill="white") +
  labs(x=expression(paste("Mean temperature of study (",degree,"C)")), y="Model partial residuals") + 
  scale_x_continuous(breaks=seq(20,28,1)) + scale_y_continuous(breaks=seq(-0.75,0.75,.25))


quad_plot




### PART 2: Climatic, non-climatic, and study factors ###

# Adjust data for running PCA #
pca_df <- master_data %>%
  dplyr::select(c(study_code, study_location,
                  # outcome 
                  effect_size,
                  # covars to go into PCs
                  log_per_capita_GDP_PPP_2015,
                  log_pop_per_m2_weighted,
                  inapparent_infection_incidence_2010,
                  total_precipitation_stdDev_weighted,
                  total_precipitation_mean_weighted,
                  mean_2m_air_temperature_stdDev_weighted,
                  dR0_dT, 
                  # study factors 
                  average_lag_in_months,
                  broad_temperature_metric, 
                  data_temporal_scale,
                  broad_disease_metric, 
                  effect_type))  %>%
  dplyr::rename(temp_sd = mean_2m_air_temperature_stdDev_weighted, 
                precip_mean = total_precipitation_mean_weighted, 
                precip_sd = total_precipitation_stdDev_weighted,
                inf_incidence = inapparent_infection_incidence_2010, 
                log_pop_density = log_pop_per_m2_weighted, 
                log_gdp = log_per_capita_GDP_PPP_2015) %>% 
  mutate(factor_temporal_scale = as.factor(data_temporal_scale))

# Set seed #
set.seed(1723)

# Run PCA with varimax rotation and 4 components using psych package #
unique_pca_obs <- pca_df %>% 
  dplyr::select(study_code, study_location, temp_sd, precip_mean, precip_sd, 
                inf_incidence, log_pop_density, log_gdp, dR0_dT) %>% 
  unique


psych_pca <- psych::principal(unique_pca_obs %>% 
                                dplyr::select(-contains("study")),
                              rotate="varimax", nfactors=4, scores=TRUE)

# Combine scores from PCA with dataset above #
pca_df <- left_join(pca_df, 
                    unique_pca_obs %>% dplyr::select(contains("study")) %>% 
                      cbind(psych_pca$scores))


### Use new dataset with PCA scores to run bootstrapped regressions #

# Run 10000 regressions, sub-sampling once per study for each run #
num.reg <- 10000

# Names of model terms to save (not extracting splines associated with lags) #
names <- c("Intercept","RC1","RC2","RC3","RC4","Mean temp.",
           "Min. temp.","Incidence", "Daily data",
           "Monthly data","Weekly data",
           "Pearson","Spearman's")

num.var <- length(names)

# Save results here, in matrix with rows = number of runs, columns = number of coefficients from model #
model.results <- matrix(nrow=num.reg,ncol=num.var)
colnames(model.results) <- names

# Create matrices for fitted lag data: 1 for x (lag), one for y (fitted correlation) #
lag.x.matrix <- matrix(nrow=100,ncol=num.reg)
lag.y.matrix <- matrix(nrow=100,ncol=num.reg)

# Each regression should have sample size = original sample size of correlation database (358 for correlations)
sample.size <- dim(pca_df)[[1]]

for(i in 1:num.reg){
  
  temp.data <- pca_df %>%
    mutate(study_code = as.factor(study_code)) %>% 
    mutate(broad_temperature_metric = as.factor(broad_temperature_metric)) %>% 
    mutate(broad_disease_metric = as.factor(broad_disease_metric))
  
  # Randomly sample from the study codes, to prevent studies with many obs from being over represented 
  temp.data <- data.frame(study_code = sample(unique(temp.data$study_code), sample.size, replace = T)) %>% 
    mutate(boot_id = 1:n()) %>% 
    left_join(temp.data, by = c("study_code")) %>% 
    # for each study code that we samples, randomly select one obs
    group_by(boot_id) %>% slice_sample(n = 1)
  
  model <- mgcv::gam(effect_size ~ RC1 + RC2 + RC3 +RC4 + 
                       s(average_lag_in_months,k=3) + 
                       broad_temperature_metric + broad_disease_metric +
                       factor_temporal_scale + effect_type, 
                     data = temp.data)
  
  # Extract regression coefficients #
  model.results[i,] <- (summary(model)$p.coeff)
  
  # Extract fitted lags #
  pd <- plot.gam(model)

  lag.x.matrix[,i] <- pd[[1]]$x 
  lag.y.matrix[,i] <- pd[[1]]$fit
  
}


# Summarize results: get mean and 2.5 and 97.5 quantiles #
model.summary <- matrix(ncol=4,nrow=c(num.var))
rownames(model.summary) <- c(names)
colnames(model.summary) <- c("variable","mean","lower.CI","upper.CI")


for(i in 1:num.var){
  model.summary[i,2] <- mean(model.results[,i])
  model.summary[i,3] <- quantile(model.results[,i],probs=0.025)
  model.summary[i,4] <- quantile(model.results[,i],probs=0.975)
}

model.summary <- as.data.frame(model.summary)
model.summary[,1] <- names

model.summary


# Separate data for RC plot and study factor plot #
model.RC <- model.summary %>% 
  filter(variable %in% c("RC1","RC2","RC3","RC4"))

model.SF <- model.summary %>% 
  filter(variable %in% c("Mean temp.",
                         "Min. temp.","Incidence", "Daily data",
                         "Monthly data","Weekly data",
                         "Pearson","Spearman's"))
         


### RC plot ###

# RC labels 
RC_labs <- psych_pca$loadings %>% 
  unclass %>% 
  as.data.frame %>% 
  rownames_to_column() %>% 
  pivot_longer(starts_with("RC")) %>% 
  filter(abs(value) > 0.6) %>% 
  mutate(lab = recode(rowname, 
                      "temp_sd" = "Temp. SD", 
                      "precip_mean" = "Precip. Mean", 
                      "precip_sd" = "Precip. SD", 
                      "inf_incidence" = "Infection Burden", 
                      "log_pop_density" = "Population Density", 
                      "log_gdp" = "GDP", 
                      "dR0_dT" = "Marginal Temp.\n Suitability"),
         lab = paste0(ifelse(value < 0, "- ", "+ "), 
                      lab, " (",  round(value, 2), ")")) %>% 
  group_by(name) %>% 
  arrange(name, desc(abs(value))) %>% 
  summarise(lab = paste(lab, collapse = "\n")) %>% 
  mutate(xval = as.numeric(gsub("RC", "", name)))

RC.plot.only <- ggplot(data=model.RC,aes(x=variable,y=mean)) + 
  theme_classic(base_size=18) + theme(legend.position="bottom",axis.text=element_text(size=20, face="bold"),
                                      axis.title.x = element_text(margin=margin(t=20,r=0,b=10,l=0),size=24),
                                      legend.title=element_text(colour="black",size=24),
                                      legend.text=element_text(colour="black",size=24),
                                      axis.title.y= element_text(size=24),
                                      title=element_text(size=24))  +
  geom_hline(yintercept = 0, color="black",linetype="dashed",size=1) +
  geom_errorbar(aes(ymin=lower.CI,ymax=upper.CI),width=.1,size=1.2) +
  geom_point(size=5, shape=21, stroke=1.5,colour="black",fill="white") +
  geom_label(data = RC_labs, 
             aes(x = xval, label = lab), y = 0.145, size = 4, fontface = 2) + 
  ylim(-0.08,0.15) +
  labs(x="Rotated Component",y="Regression Coefficient")

RC.plot.only


# Study factor plot #
SF.plot.only <- ggplot(data=model.SF,aes(x=fct_inorder(variable),y=mean)) +
  theme_bw(base_size=22) + theme(legend.position="bottom",axis.text=element_text(size=18, face="bold"),
                                 axis.title.x = element_text(margin=margin(t=20,r=0,b=10,l=0),size=24),
                                 legend.title=element_text(colour="black",size=24),
                                 legend.text=element_text(colour="black",size=24),
                                 axis.title.y= element_text(size=24),
                                 title=element_text(size=24)) +
  geom_hline(yintercept = 0, color="black",linetype="dashed",size=1) +
  geom_errorbar(aes(ymin=lower.CI,ymax=upper.CI),width=.3,size=1.5)  +
  geom_point(size=6, stroke=2,fill="white", shape=21) + 
  labs(x="Study Factor",y="Effect of Study Factor on Correlation")

SF.plot.only


# Supplemental figure #
# Lag plot (in base R) #
# Calculate mean value across x #

mean.cor.over.lags <- matrix(nrow=100,ncol=2)
mean.cor.over.lags[,1] <- lag.x.matrix[,1]

for(i in 1:100){
  mean.cor.over.lags[i,2] <- mean(lag.y.matrix[i,])
}
  
par(mar=c(5,5,1,1))
plot(lag.x.matrix[,1], lag.y.matrix[,1], type="n", ylim=c(-0.25,0.25),
     ylab="Fitted temperature-dengue correlation",xlab="Lag (months)",cex.lab=1.5, cex.axis=1.5)
for(i in 1:num.reg){
  lines(lag.x.matrix[,i], lag.y.matrix[,i], col=alpha("black", alpha=0.03))
}
lines(mean.cor.over.lags[,2] ~ mean.cor.over.lags[,1], col="blue",lwd=3)



# Supplemental figure: quadratic mean temp model fit with raw data instead of residuals #
quad_plot_raw <- ggplot() + theme_classic(base_size=18) + theme(legend.position="none",axis.text=element_text(size=20),
                                                                axis.title.x = element_text(margin=margin(t=20,r=0,b=10,l=0),size=24),
                                                                axis.title.y= element_text(size=24)) +
  geom_ribbon(data=x_quadratic, aes(x=mean_temp_C, ymin=lower, ymax=upper), alpha= 0.2, fill="black") +
  geom_line(data=x_quadratic, aes(x=mean_temp_C, y=fit), colour="black",size=2)  +
  geom_point(data=subsetted_data, aes(mean_temp_C, effect_size), size=4 ,shape=21,stroke=1.5,fill="white") +
  labs(x=expression(paste("Mean temperature of study (",degree,"C)")), y="Reported correlation") + 
  scale_x_continuous(breaks=seq(20,28,1)) + scale_y_continuous(breaks=seq(-0.75,0.75,.25))


quad_plot_raw


# Supplemental figure: 
# Correlation plot for PCA factors #
cor.factors <- unique_pca_obs %>%
  dplyr::select(-(c(study_code,study_location))) %>%
  rename(Temperature_SD = temp_sd,
         Precipitation_SD = precip_sd,
         Mean_precipitation = precip_mean,
         Infection_burden = inf_incidence,
         Log_population_density = log_pop_density,
         Log_GDP = log_gdp,
         dR0dT = dR0_dT)

M <- cor(cor.factors)

corrplot(M, method="number")






