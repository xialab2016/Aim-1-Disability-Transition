#######################################
## Post-hoc analysis of msms results 
#########################################

##set working dir
setwd("C:/Users/vfngwa/OneDrive - University of Tasmania/Documents/InForMS Project/Hu chen analysis of pdds")


#load utility functions
source("utility-functions.R")

# Prepare model output for promote.fit
load("promote.fit.rdata")
eff <- f(fit = promote.fit)
coefs <- eff$beta
names(coefs) <- eff$cov
covmat <- as.data.frame(promote.fit$covmat)
rownames(covmat) <- eff$cov
colnames(covmat) <- eff$cov

# Transition counts for Promote

Q1 <- rbind(
  c(0, 563, 145, 61, 52),
  c(534, 0, 214, 151, 64),
  c(113, 229, 0, 142, 73),
  c(41, 129, 122, 0, 100),
  c(0, 0, 0, 0, 0)
)
# Get all parameter-level results
promote.sum <- get_parameter_estimates(coefs, covmat, Q1, k = 8)
print(promote.sum)

# Compute DMT effects by age group for PROMOTE
promote.dmteff <- compute_dmt_effects_by_age(coefs, covmat, Q1, k=8)
promote.dmteff <- promote.dmteff[order(promote.dmteff$AgeGroup, promote.dmteff$Treatment), ]
print(promote.dmteff, digits = 4)


# Prepare model output for AMSLS
load("amsls.fit.rdata")
eff <- f(fit = amsls.fit)
coefs <- eff$beta
names(coefs) <- eff$cov
covmat <- as.data.frame(amsls.fit$covmat)
rownames(covmat) <- eff$cov
colnames(covmat) <- eff$cov

# Transition counts for AMSLS
Q2 <- rbind(
  c(0, 824, 255, 198,  70),
  c(830,  0, 213, 574, 179),
  c(229, 216,   0, 121,  25),
  c(187, 495, 121,   0, 247),
  c(0, 0, 0, 0, 0)
)
amsls.sum <- get_parameter_estimates(coefs, covmat, Q2, k = 8)
print(amsls.sum)

# Compute DMT effects by age group for AMSLS
amsls.dmteff <- compute_dmt_effects_by_age(coefs, covmat, Q2, k=8)
amsls.dmteff <- amsls.dmteff[order(amsls.dmteff$AgeGroup, amsls.dmteff$Treatment), ]
print(amsls.dmteff, digits = 4)


##Meta-anlysis of promote and amsls
est_prom<-promote.dmteff[,1:5]
est_prom$cov<-ifelse(est_prom$AgeGroup=="Middle"&est_prom$Treatment=="Higher-efficacy", paste("bHM", est_prom$Transition, sep=""),NA) 
est_prom$cov<-ifelse(est_prom$AgeGroup=="Older"&est_prom$Treatment=="Higher-efficacy", paste("bHO", est_prom$Transition, sep=""),est_prom$cov) 
est_prom$cov<-ifelse(est_prom$AgeGroup=="Young"&est_prom$Treatment=="Higher-efficacy", paste("bHY", est_prom$Transition, sep=""),est_prom$cov) 

est_prom$cov<-ifelse(est_prom$AgeGroup=="Middle"&est_prom$Treatment=="Standard-efficacy", paste("bLM", est_prom$Transition, sep=""),est_prom$cov) 
est_prom$cov<-ifelse(est_prom$AgeGroup=="Older"&est_prom$Treatment=="Standard-efficacy", paste("bLO", est_prom$Transition, sep=""),est_prom$cov) 
est_prom$cov<-ifelse(est_prom$AgeGroup=="Young"&est_prom$Treatment=="Standard-efficacy", paste("bLY", est_prom$Transition, sep=""),est_prom$cov) 
est_prom$Covariate <- ifelse(
  est_prom$AgeGroup == "Middle" & est_prom$Treatment == "Higher-efficacy",
  "Middle-aged adults on higher-efficacy DMT vs untreated",
  ifelse(
    est_prom$AgeGroup == "Older" & est_prom$Treatment == "Higher-efficacy",
    "Older adults on higher-efficacy DMT vs untreated",
    ifelse(
      est_prom$AgeGroup == "Middle" & est_prom$Treatment == "Standard-efficacy",
      "Middle-aged adults on standard-efficacy DMT vs untreated",
      ifelse(
        est_prom$AgeGroup == "Older" & est_prom$Treatment == "Standard-efficacy",
        "Older adults on standard-efficacy DMT vs untreated",
        ifelse(
          est_prom$AgeGroup == "Young" & est_prom$Treatment == "Higher-efficacy",
          "Young adults on higher-efficacy DMT vs untreated",
          "Young adults on standard-efficacy DMT vs untreated"
        )
      )
    )
  )
)

est_prom<-est_prom[,c(1,7,6,4,5)]
names(est_prom)[c(1,3:5)]<-c("trans", "cov", "beta", "se")


est_amsls<-amsls.dmteff[,1:5]
est_amsls$cov<-ifelse(est_amsls$AgeGroup=="Middle"&est_amsls$Treatment=="Higher-efficacy", paste("bHM", est_amsls$Transition, sep=""),NA) 
est_amsls$cov<-ifelse(est_amsls$AgeGroup=="Older"&est_amsls$Treatment=="Higher-efficacy", paste("bHO", est_amsls$Transition, sep=""),est_amsls$cov) 
est_amsls$cov<-ifelse(est_amsls$AgeGroup=="Young"&est_amsls$Treatment=="Higher-efficacy", paste("bHY", est_amsls$Transition, sep=""),est_amsls$cov) 

est_amsls$cov<-ifelse(est_amsls$AgeGroup=="Middle"&est_amsls$Treatment=="Standard-efficacy", paste("bLM", est_amsls$Transition, sep=""),est_amsls$cov) 
est_amsls$cov<-ifelse(est_amsls$AgeGroup=="Older"&est_amsls$Treatment=="Standard-efficacy", paste("bLO", est_amsls$Transition, sep=""),est_amsls$cov) 
est_amsls$cov<-ifelse(est_amsls$AgeGroup=="Young"&est_amsls$Treatment=="Standard-efficacy", paste("bLY", est_amsls$Transition, sep=""),est_amsls$cov) 
est_amsls$Covariate <- ifelse(
  est_amsls$AgeGroup == "Middle" & est_amsls$Treatment == "Higher-efficacy",
  "Middle-aged adults on higher-efficacy DMT vs untreated",
  ifelse(
    est_amsls$AgeGroup == "Older" & est_amsls$Treatment == "Higher-efficacy",
    "Older adults on higher-efficacy DMT vs untreated",
    ifelse(
      est_amsls$AgeGroup == "Middle" & est_amsls$Treatment == "Standard-efficacy",
      "Middle-aged adults on standard-efficacy DMT vs untreated",
      ifelse(
        est_amsls$AgeGroup == "Older" & est_amsls$Treatment == "Standard-efficacy",
        "Older adults on standard-efficacy DMT vs untreated",
        ifelse(
          est_amsls$AgeGroup == "Young" & est_amsls$Treatment == "Higher-efficacy",
          "Young adults on higher-efficacy DMT vs untreated",
          "Young adults on standard-efficacy DMT vs untreated"
        )
      )
    )
  )
)
est_amsls<-est_amsls[,c(1,7,6,4,5)]
names(est_amsls)[c(1,3:5)]<-c("trans", "cov", "beta", "se")

# Default fixed-effect meta-analysis
meta.fixed <- meta.F(est_prom, est_amsls, method = "fixed")
meta.fixed$P_Meta<-as.numeric(meta.fixed$P_Meta)
# Random-effects meta-analysis
meta.random <- meta.F(est_prom, est_amsls, method = "random")
meta.random$P_Meta<-as.numeric(meta.random$P_Meta)


##make sure to cre


##Make forest plots of results

# 1. Forest Plot
library(ggplot2)
library(dplyr)
library(metafor)

res<-meta.fixed
res$Covariate<-est_amsls$Covariate
res$Outcome<-as.factor(ifelse(res$Transition=="1-2", "Mild disability",
                              ifelse(res$Transition%in%c("1-3", "2-3"), "Moderate disability",
                                     ifelse(res$Transition%in%c("1-4", "2-4", "3-4"), "Gait disability",
                                            ifelse(res$Transition%in%c("1-5", "2-5", "3-5", "4-5"), "Early cane", NA)))))
nm<-names(res)[3:11]
res<-res %>% mutate(across(all_of(nm), as.numeric))
res$Covariate <- sub("_[0-9]+$", "", res$Covariate)
res<-data.frame(res)
res=res[res$Transition%in%c("1-4", "1-5", "2-4", "2-5", "3-4", "3-5", "4-5" ),]
res=res[order(res$Covariate,res$Outcome,res$Transition, decreasing = T),]

plot.new()
tiff(filename = "meta-analysis.tif",width = 14, height = 10,units = "in",pointsize = 10,res = 300, compression = "lzw")
atx<-round(seq(from=min(log(range(res$HazardRatio_Meta)))+0.1,to=max(log(range(res$HazardRatio_Meta)))+0.5,by=0.3),2)
custom_forest_plot(res,atx = atx, cex.hr=1.0, rectp=3.5,cex.xp = 0.8, xlim=c(-4.5, 2.2),covp=4.5, txtcovp_wd=35, stdp=c(3.4, 2.9), trp=2.4, hrp=1.6,  outp=1.7, alpha=0.23)



####################################################################################
## Overal treatment effect taking into account teatment delay, duration on treatment
####################################################################################

keep_trans <- c("1-2","1-3", "1-4", "1-5", "2-3", "2-4", "2-5", "3-4", "3-5", "4-5")
keep_cov_prefix <- c("bMA", "bOA", "bH", "bL", "bdH", "bdL", "bHM", "bHO", "bLM", "bLO", "bTD")

prdf<-promote.sum[, 1:3]
names(prdf)<-c("cov", "beta", "se")
prdf$trans <- paste0(
  substr(prdf$cov, nchar(prdf$cov) - 1, nchar(prdf$cov) - 1),
  "-",
  substr(prdf$cov, nchar(prdf$cov), nchar(prdf$cov))
)
prdf$prefix <- sub("([a-zA-Z]+).*", "\\1", prdf$cov)
prdf <- prdf[prdf$trans %in% keep_trans & prdf$prefix %in% keep_cov_prefix, ]
prdf$prefix<-NULL
##AMSLS selected transition
amsdf<-amsls.sum[, 1:3]
names(amsdf)<-c("cov", "beta", "se")
amsdf$trans <- paste0(
  substr(amsdf$cov, nchar(amsdf$cov) - 1, nchar(amsdf$cov) - 1),
  "-",
  substr(amsdf$cov, nchar(amsdf$cov), nchar(amsdf$cov))
)
amsdf$prefix <- sub("([a-zA-Z]+).*", "\\1", amsdf$cov)
amsdf <- amsdf[amsdf$trans %in% keep_trans & amsdf$prefix %in% keep_cov_prefix, ]
amsdf$prefix<-NULL
meta.res <- meta.F(prdf, amsdf, method = "fixed")
meta.res <-meta.res[, c(1,2,4)]

##Extract the coefficient lists
coef_list <- extract_coef_list(meta.res, transitions = keep_trans)

# Names for age groups and DMT types
age_groups <- c("Young", "Middle-aged", "Older")
dmt_types <- c("high", "low")
delays <- c(0, -log(1)-1)   # No delay and 1-year delay

#Treatment duration vector: 0 to 5 years in 3-month steps
treatment_duration <- log(seq(1, exp(6), by = exp(0.25)))  #use log(seq(1, 6, by = 0.25)) instead

# Create empty data frame to store results
res.trt <- data.frame()
for (trans in names(coef_list)) {
  coef <- coef_list[[trans]]
  for (age in age_groups) {
    for (dmt in dmt_types) {
      for (delay in delays) {
        for (dur in treatment_duration) {
          HR <- compute_HR(age, delay, dur, dmt, coef=coef)
          res.trt <- rbind(res.trt, data.frame(
            Trans=trans,
            AgeGroup = age,
            DMT_Type = ifelse(dmt == "high", "Higher-efficacy DMTs", "Standard-efficacy DMTs"),
            Delay = ifelse(delay == 0, "No Delay", "1-Year Delay"),
            Duration = dur,
            HR = HR
          ))
        }
      }
    }
  }
}

res.trt$Outcome<-as.factor(ifelse(res.trt$Trans=="1-2", "Mild disability",
                        ifelse(res.trt$Trans%in%c("1-3", "2-3"), "Moderate diability",
                               ifelse(res.trt$Trans%in%c("1-4", "2-4", "3-4"), "Gait disability", "Early cane"))))

res.trt$Duration <- 1 + (res.trt$Duration/ 6) * (5 - 1)
library(ggplot2)
library(mgcv)  # required for 's()' formula in GAM
library(gamlss)
library(splines)


p1<-ggplot(res.trt, aes(x = Duration, y = HR, color = AgeGroup, linetype = Delay)) +
  geom_smooth(method = "loess", formula = y ~ ns(x, 1), se = FALSE, size = 1.2) +
  facet_grid(~ DMT_Type) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_color_manual(
    values = c(
      "Older" = "green4",
      "Middle-aged" = "orange3",
      "Young" = "red4"
    )
  ) +
  labs(
    title = " ",
    x = "Time on disease-modifying therapies (years)",
    y = "Hazard ratio of disability progression",
    color = "Age Group",
    linetype = "Treatment Delay"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "grey80", size = 0.3),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12)
  ) +
  guides(
    linetype = guide_legend(override.aes = list(color = "black"))
  )

if(1){plot.new()
tiff(filename = "overall_dmt_effects.tiff",width = 12, height = 14,units = "in",pointsize = 14,res = 300, compression = "lzw")
p1}


#####################################################
## Bubble plot of DMT PMR ratio across treatment age
#####################################################
###Get dynamic transition probabilities


# Load libraries
library(msm)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("C:/Users/vfngwa/OneDrive - University of Tasmania/Documents/InForMS Project/Hu chen analysis of pdds")

#load data
load("msdat.rdata")
#load("amsls.fit.rdata")
load("promote.fit.rdata")

source("utility-function.R")


msdat$tagecat <- cut(exp(msdat$tage),
                   breaks = c(0, 30, 50, Inf),
                   labels = c("0", "1", "2"),
                   right = FALSE)
msdat$tagecat<-relevel(msdat$tagecat, ref = "0")

msdat$tagecat1<-ifelse(msdat$tagecat=="1",1,0)
msdat$tagecat2<-ifelse(msdat$tagecat=="2",1,0)

res.prob <- predict_trans_prob_from_data(
  model = promote.fit,
  data = msdat,
  return_type = "probability",
  covariate_names = c("dd", "tdelay", "tagecat1", "tagecat2", "on_high", "on_low", "ddmt_high", "ddmt_low", "sex", "tsymp", "progstat"),
  time_var = "obstime",
  state_var = "state",
  ci="delta"
)


# Make bubble plots for different metrics and save data
#######################################
bubbleplot(
  model = promote.fit, #change fit amsls.fit if validating on amsls data.
  data = msdat,
  tp_results=res.prob, #if you specify tp_results, then function will be faster, otherwise, tp_results will be computed before plots are made=slower
  return_type = "probability",
  covariate_names = c("dd", "tdelay", "tagecat1", "tagecat2", "on_high", "on_low", "ddmt_high", "ddmt_low", "sex", "tsymp", "progstat"),
  time_var = "obstime",
  state_var = "state",
  saveplotdata = TRUE,
  title = "Promote PMR and CondWorsen",
  #metrics = c("PMR_norm", "CondWorsen", "CondWorsen_percent_reduction", "PMR_percent_relative", "TotalWorsenRate", "TotalImproveRate")
)



