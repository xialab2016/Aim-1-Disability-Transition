#######################################
## Post-hoc analysis of msms results 
#########################################
library(dplyr)
library(tableone)
library(reshape2)
library(ggplot2)
library(msm)
setwd("C:\\Users\\pumpk\\OneDrive - University of Pittsburgh\\XiaLab\\MSM\\2025")

#### Table 1 participant inclusion/exclusion ####
demo=read.csv("../DATA/2025/CleanDemographics.csv")
demo=demo %>% mutate(sex=ifelse(subject_sex==1, 1,0 ),
                     progstat=ifelse(subtype_enroll %in% c("PPMS","SPMS"), 1,0 ),# 0: CIS/RIS/RRMS; 1: PMS
                     datediag= as.Date(date_msdx),
                     datefde= as.Date(date_firstsx),
                     datediag = ifelse(is.na(datediag) & !is.na(datefde), as.character(datefde), as.character(datediag)),
                     datefde = ifelse(is.na(datefde) & !is.na(datediag), as.character(datediag), as.character(datefde)),
                     dage = round(as.numeric(as.Date(datediag) - as.Date(dob)) / 365.25, 2),
                     oage = round(as.numeric(as.Date(datefde) - as.Date(dob)) / 365.25, 2),
                     tsymp = abs(round(as.numeric(as.Date(datediag) - as.Date(datefde)) / 365.25, 2))
) %>%
  select(id,dob, sex, progstat, datediag, datefde, dage, oage, tsymp )
dim(demo)  #n1 in figure 1
disab = read.csv("../DATA/2025//Clean_PDDS_20241001.csv")
dim(disab[disab$id %in% demo$id,]) #n2 in figure 1
data=left_join(demo, disab, by="id")
#exclude ppts with missingness
data= data %>% 
  filter(is.na(dob)==F & is.na(sex)==F & is.na(progstat)==F & is.na(tsymp)==F&is.na(score)==F)
dim(data); length(unique(data$id)) #n3 and n4 in figure 1
#exclude ppts had no MS diagnosis or onset age <10 
data=data %>% filter(oage>=10 )
dim(data); length(unique(data$id))#n5 and n6 in figure 1
#exclude pdds recorded <3 months aprt 
data = data %>%
  arrange(id, reviewdate) %>%
  group_by(id) %>%
  mutate(
    lagrev = lag(reviewdate),
    t_interval = as.numeric(as.Date(reviewdate) - as.Date(lagrev)) / 365.25,
    t_interval=ifelse(is.na(t_interval)==T, 0, t_interval),
    years=cumsum(t_interval),
    state = case_when(
      score == 0 ~ 1,
      score == 1 ~ 2,
      score == 2 ~ 3,
      score == 3 ~ 4,
      TRUE ~ 5
    )
  ) %>%
  filter( (t_interval == 0 | t_interval >= 0.25))   %>%
  group_by(id) %>%
  mutate(n = n())      
dim(data); length(unique(data$id))#n7 and n8 in figure 1
#exclude ppts had <2 PDDS 
data=data %>% mutate(include = ifelse(n <2, 0,1)) 
length(unique(data$id[data$include==1]))   #n9 in figure 1
length(data$score[data$include==1])#n10 in figure 1


dmt <- read.csv("../DATA/2025//CleanDMT.csv")
colnames(dmt)[1]="id"
dmt=dmt %>% group_by(id) %>% 
  summarise(ever_high=sum(high=ifelse(sum(Efficacy==2)>=1,1,0)),
            ever_low=ifelse(sum(Efficacy==1)>=1,1,0),
            dmtstart= min(as.Date(dmtstartdate[Efficacy !=0])))
data= left_join(data, dmt, by ="id")
data= data %>% 
  group_by(id)%>%
  mutate(
    fup = max(years),
    dd= as.numeric(min(as.Date(reviewdate))-as.Date(datefde))/365.25,
    trtdelay = as.numeric(as.Date(dmtstart)-as.Date(datediag))/365.25,
    trtdelay = ifelse(is.finite(trtdelay)==F, NA, trtdelay),
    ever_high=ifelse(is.na(ever_high)==T | ever_high==0,0,1),
    ever_low=ifelse(is.na(ever_low)==T | ever_low==0,0,1),
    never=ifelse(ever_high==0 & ever_low==0,1, 0))

data_bl= data %>% distinct(id, .keep_all = T)
data_bl_include= data_bl %>% filter(include==1)
tab1=CreateTableOne(data=data_bl_include,
                    vars = c("sex","oage","tsymp","progstat","fup","dd","score","ever_high","ever_low","never","trtdelay"),
                    factorVars = c("sex","progstat","ever_high","ever_low","never"))
print(tab1,nonnormal ="score" ) # Table 1

stab1=CreateTableOne(data=data_bl,
                     strata = "include",
                     vars = c("sex","oage","tsymp","progstat","fup","dd","score","ever_high","ever_low","never","trtdelay"),
                     factorVars = c("sex","progstat","ever_high","ever_low","never"))
print(stab1,nonnormal ="score" ) #Table 2
statetable.msm(state = state, subject = id, data = data[which(data$include==1),]) #transition matrix




# Step 1: Create the raw count matrix
transition_counts <- matrix(c(
  2331, 600, 160, 68, 89,
  573, 891,242, 176, 151,
  122, 264, 394, 178, 156,
  49, 152, 167, 437, 169,
  7, 10, 10, 15, 2827
), nrow = 5, byrow = TRUE)

sum(transition_counts)+1897


rownames(transition_counts) <- c("No","Mild","Moderate","Gait","Early cane or worse")
colnames(transition_counts) <- c("No","Mild","Moderate","Gait","Early cane or worse")
transition_props <- prop.table(transition_counts, margin = 1)


df_counts <- melt(transition_counts)
df_props <- melt(transition_props)


df_long <- merge(df_counts, df_props, by = c("Var1", "Var2"))
colnames(df_long) <- c("From", "To", "Count", "Proportion")
df_long= df_long %>% mutate(To=factor(To, levels=c("No","Mild","Moderate","Gait","Early cane or worse")),
                            From=factor(From, levels=c("Early cane or worse","Gait","Moderate","Mild","No")))


# Step 4: Plot
p1=ggplot(df_long, aes(x = To, y = From, fill = Proportion)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(Count, "\n(", round(Proportion*100, 1), "%)")), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Proportion") +
  labs(title = "A. Number of transitions in PROMOTE",
       x = "To State", y = "From State") +
  scale_x_discrete(position = "top")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0),
        legend.position = "none",
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())


transition_counts <- matrix(c(
  5008, 999, 294, 213, 61,
  985, 2656,256, 896, 181,
  261, 281, 786, 201, 115,
  176, 708, 220, 2929, 705,
  8, 21, 21, 29, 7060
), nrow = 5, byrow = TRUE)
sum(transition_counts)+2286

rownames(transition_counts) <- c("No","Mild","Moderate","Gait","Early cane or worse")
colnames(transition_counts) <- c("No","Mild","Moderate","Gait","Early cane or worse")
transition_props <- prop.table(transition_counts, margin = 1)


df_counts <- melt(transition_counts)
df_props <- melt(transition_props)


df_long <- merge(df_counts, df_props, by = c("Var1", "Var2"))
colnames(df_long) <- c("From", "To", "Count", "Proportion")
df_long= df_long %>% mutate(To=factor(To, levels=c("No","Mild","Moderate","Gait","Early cane or worse")),
                            From=factor(From, levels=c("Early cane or worse","Gait","Moderate","Mild","No")))


# Step 4: Plot
p2=ggplot(df_long, aes(x = To, y = From, fill = Proportion)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(Count, "\n(", round(Proportion*100, 1), "%)")), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Proportion") +
  labs(title = "B. Number of transitions in AMSLS",
       x = "To State", y = "From State") +
  scale_x_discrete(position = "top")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0),,
        legend.position = "none",axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
p1/
  p2
ggsave("Results/Figure3_Transition matrix.png",width = 10, height = 10, units = "in")

data3= data2 %>% group_by(id) %>% summarise(fup=max(years))
summary(data3$fup); sd(data3$fup)
data3 = data2 %>% distinct(id, sex,tsymp,progstat)
data3= data2 %>% group_by(id) %>% summarise(dd= min(as.Date(reviewdate)-as.Date(datefde)),
                                            dd=as.numeric(dd/365.25))

summary(data3$dd)
sd(data3$dd)

dmt <- read.csv("../DATA/2025//CleanDMT.csv")
dmt=dmt %>%filter(id_participant %in% data2$id)
length(unique(data2[which(data2$id %in% dmt$id_participant==F),"id"]$id)) #128 no treted
dmt2=dmt %>% group_by(id_participant) %>% summarise(high=sum(high=ifelse(sum(Efficacy==2)>=1,1,0)),
                                                   std=ifelse(sum(Efficacy==1)>=1,1,0))
table(dmt2$high)
##set working dd##set working dir
#setwd("~/OneDrive - University of Pittsburgh/XiaLab/MSM/2025/")


#load utility functions
source("utility-functions.R")

# Prepare model output for promote.fit
load("Results/promote.fit.rdata")
eff <- f(fit = promote.fit)
coefs <- eff$beta
names(coefs) <- eff$cov
covmat <- as.data.frame(promote.fit$covmat)
rownames(covmat) <- eff$cov
colnames(covmat) <- eff$cov

# Transition counts for Promote

# Q1 <- rbind(
#   c(0, 563, 145, 61, 52),
#   c(534, 0, 214, 151, 64),
#   c(113, 229, 0, 142, 73),
#   c(41, 129, 122, 0, 100),
#   c(0, 0, 0, 0, 0)
# )


 Q1 <- rbind(
   c(0, 600, 160, 68, 89),
   c(573, 0, 242, 176, 151),
   c(122, 264, 0, 178, 156),
   c(49, 152, 167, 0, 169),
   c(0, 0, 0, 0, 0)
 )
# Get all parameter-level results
promote.sum <- get_parameter_estimates(coefs, covmat, Q1, k = 8)
print(promote.sum)
promote.sum = promote.sum %>%
  mutate(include = ifelse(grepl("12",Parameter)==T |grepl("13",Parameter)==T |grepl("14",Parameter)==T |grepl("15",Parameter)==T |grepl("23",Parameter)==T |grepl("24",Parameter)==T |grepl("25",Parameter)==T |grepl("34",Parameter)==T |grepl("35",Parameter)==T |grepl("45",Parameter)==T  ,1,0))%>%
  filter(include==1) %>% 
  select(Parameter, Estimate, StdError)
#write.csv(promote.sum.out, "Results/promote.sum.csv")
# Prepare model output for AMSLS
load("Results/amsls.fit.250430.rdata")
eff <- f(fit = amsls.fit)
coefs <- eff$beta
names(coefs) <- eff$cov
covmat <- as.data.frame(amsls.fit$covmat)
rownames(covmat) <- eff$cov
colnames(covmat) <- eff$cov

# Transition counts for AMSLS
# Q2 <- rbind(
#   c(0, 824, 255, 198,  70),
#   c(830,  0, 213, 574, 179),
#   c(229, 216,   0, 121,  25),
#   c(187, 495, 121,   0, 247),
#   c(0, 0, 0, 0, 0)
# )

Q2 <- rbind(
  c(0, 999, 294, 213,  61),
  c(985,  0, 256, 896, 181),
  c(261, 281,   0, 201,  115),
  c(176, 708, 220,   0, 705),
  c(0, 0, 0, 0, 0)
)
amsls.sum <- get_parameter_estimates(coefs, covmat, Q2, k = 8)

amsls.sum= amsls.sum %>%
  mutate(include = ifelse(grepl("12",Parameter)==T |grepl("13",Parameter)==T |grepl("14",Parameter)==T |grepl("15",Parameter)==T |grepl("23",Parameter)==T |grepl("24",Parameter)==T |grepl("25",Parameter)==T |grepl("34",Parameter)==T |grepl("35",Parameter)==T |grepl("45",Parameter)==T  ,1,0))%>%
  filter(include==1)%>% 
  select(Parameter, Estimate, StdError)
#write.csv(amsls.sum.out, "Results/amsls.sum.csv")
merged_data <- merge(
  promote.sum, amsls.sum,
  by = "Parameter",
  suffixes = c("_promote", "_amsls")
)
write.csv(merged_data, "Results/Table2_merged.csv")
merged_data = read.csv("Results/Table2_merged.csv")
##Table 2 ##
results <- data.frame(
    Covariate = merged_data$Parameter,
    Beta_Promote = merged_data$Estimate_promote,
    SE_Promote = merged_data$StdError_promote,
    Beta_AMSLS = merged_data$Estimate_amsls,
    SE_AMSLS = merged_data$StdError_amsls,
    HazardRatio_Meta = NA,
    LogHR_Meta = NA,
    LogHR_Meta_SE = NA,
    LogHR_LowerCI = NA,
    LogHR_UpperCI = NA,
    P_Meta = NA,
    Q_Statistic = NA,
    I2_Percent = NA,
    Tau2 = NA,
    PROMOTE = NA,
    AMSLS = NA
  )
  
  # Perform meta-analysis row by row
  for (i in seq_len(nrow(merged_data))) {
    # Extract log HRs and SEs
    b.est <- c(merged_data$Estimate_promote [i], merged_data$Estimate_amsls [i])
    se <- c(merged_data$StdError_promote [i], merged_data$StdError_amsls[i])
    
    # Fixed-effect weights
    w.fixed <- 1 / se^2
    b.fixed <- sum(w.fixed * b.est) / sum(w.fixed)
    se.fixed <- sqrt(1 / sum(w.fixed))
    
    # Q-statistic and heterogeneity
    q.stat <- sum(w.fixed * (b.est - b.fixed)^2)
    df.q <- length(b.est) - 1
    i2 <- max(0, (q.stat - 10) / q.stat) * 100
    
    # DerSimonian-Laird tau² for random effects
    tau2 <- max(0, (q.stat - df.q) / (sum(w.fixed) - sum(w.fixed^2) / sum(w.fixed)))
    w.random <- 1 / (se^2 + tau2)
    
    # Random effects pooled estimate
    b.random <- sum(w.random * b.est) / sum(w.random)
    se.random <- sqrt(1 / sum(w.random))
    
   
      b.meta <- b.fixed
      se.meta <- se.fixed
    
    # CI and p-value
    ci.lower <- b.meta - 1.96 * se.meta
    ci.upper <- b.meta + 1.96 * se.meta
    p.meta <- pchisq((b.meta / se.meta)^2, df = 1, lower.tail = FALSE)
    
    # Hazard ratios
    hr.meta <- exp(b.meta)
    hr.lower <- exp(ci.lower)
    hr.upper <- exp(ci.upper)
    hr_promote <- exp(b.est[1])
    hr_amsls <- exp(b.est[2])
    
    # Store results
    results[i, ] <- c(
      merged_data$Parameter [i],
      merged_data$Estimate_promote[i],
      merged_data$StdError_promote[i],
      merged_data$Estimate_amsls [i],
      merged_data$StdError_amsls [i],
      round(hr.meta, 2),
      round(b.meta, 4),
      round(se.meta, 4),
      round(hr.lower, 2),
      round(hr.upper, 2),
      formatC(p.meta, format = "e", digits = 2),
      round(q.stat, 2),
      round(i2, 2),
      round(tau2, 4),
      paste0(round(hr_promote, 2), " (", round(merged_data$StdError_promote [i], 2), ")"),
      paste0(round(hr_amsls, 2), " (", round(merged_data$StdError_amsls[i], 2), ")")
    )
  }
  
 

meta.fixed <- results
meta.fixed  = meta.fixed %>% 
  mutate(Beta_Promote=as.numeric(Beta_Promote),
         SE_Promote=as.numeric(SE_Promote),
         Beta_AMSLS=as.numeric(Beta_AMSLS),
         SE_AMSLS=as.numeric(SE_AMSLS),
         HR_PROMOTE= paste0(round(exp(Beta_Promote),2), " [",round(exp(Beta_Promote-qnorm(0.975)*SE_Promote),2),", ",round(exp(Beta_Promote+qnorm(0.975)*SE_Promote),2), "]"),
         HR_AMSLS= paste0(round(exp(Beta_AMSLS),2), " [",round(exp(Beta_AMSLS-qnorm(0.975)*SE_AMSLS),2),", ",round(exp(Beta_AMSLS+qnorm(0.975)*SE_AMSLS),2), "]"),
         HR_meta= paste0(HazardRatio_Meta," [", round(as.numeric(LogHR_LowerCI),2),", ",round(as.numeric(LogHR_UpperCI),2),"]"))%>%
  select(Covariate, Beta_Promote, SE_Promote, HR_PROMOTE,
         Beta_AMSLS, SE_AMSLS, HR_AMSLS,
         HazardRatio_Meta,LogHR_LowerCI,LogHR_UpperCI, HR_meta, P_Meta, Q_Statistic, I2_Percent, Tau2)%>% 
  filter(Covariate %in% promote.sum$Parameter) %>% 
  mutate(Covariate=factor(Covariate, levels=unique(promote.sum$Parameter)))%>%
  arrange(Covariate)

write.csv(meta.fixed, "Results/Table2.csv")



### Figure 4 ###

# Compute DMT effects by age group for PROMOTE
load("Results/promote.fit.rdata")
eff <- f(fit = promote.fit)
coefs <- eff$beta
names(coefs) <- eff$cov
covmat <- as.data.frame(promote.fit$covmat)
rownames(covmat) <- eff$cov
colnames(covmat) <- eff$cov
promote.dmteff <- compute_dmt_effects_by_age(coefs, covmat, Q1, k=8)
promote.dmteff <- promote.dmteff[order(promote.dmteff$AgeGroup, promote.dmteff$Treatment), ]
print(promote.dmteff, digits = 4)

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
          "Younger adults on higher-efficacy DMT vs untreated",
          "Younger adults on standard-efficacy DMT vs untreated"
        )
      )
    )
  )
)

est_prom<-est_prom[,c(1,7,6,4,5)]
names(est_prom)[c(1,3:5)]<-c("trans", "cov", "beta", "se")
est_prom



# Compute DMT effects by age group for AMSLS
load("Results/amsls.fit.250430.rdata")
eff <- f(fit = amsls.fit)
coefs <- eff$beta
names(coefs) <- eff$cov
covmat <- as.data.frame(amsls.fit$covmat)
rownames(covmat) <- eff$cov
colnames(covmat) <- eff$cov
amsls.dmteff <- compute_dmt_effects_by_age(coefs, covmat, Q2, k=8)
amsls.dmteff <- amsls.dmteff[order(amsls.dmteff$AgeGroup, amsls.dmteff$Treatment), ]
print(amsls.dmteff, digits = 4)
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
          "Younger adults on higher-efficacy DMT vs untreated",
          "Younger adults on standard-efficacy DMT vs untreated"
        )
      )
    )
  )
)
est_amsls<-est_amsls[,c(1,7,6,4,5)]
names(est_amsls)[c(1,3:5)]<-c("trans", "cov", "beta", "se")

meta.fixed <- meta.F(est_prom, est_amsls, method = "fixed")

merged = merge(est_prom, est_amsls, by=c("trans","Covariate","cov"), suffixes = c("_prom","_amsls"))
write.csv(merged, "Results/DMT_by_Age.csv", row.names = F)


# Default fixed-effect meta-analysis for DMT effect
merged=read.csv("Results/DMT_by_Age.csv")
est_prom = merged %>% select(trans, Covariate, cov, beta_prom, se_prom)
colnames(est_prom)[4:5]=c("beta","se")
est_amsls = merged %>% select(trans, Covariate, cov, beta_amsls, se_amsls)
colnames(est_amsls)[4:5]=c("beta","se")
est_prom=est_prom %>% mutate(HR = round(exp(beta),2),
                             Lower= round(exp(beta-qnorm(0.975)*se),2),
                             Upper= round(exp(beta+qnorm(0.975)*se),2),
                             est=paste0(HR, " [",Lower, ", ", Upper,"]"))
est_amsls=est_amsls %>% mutate(HR = round(exp(beta),2),
                               Lower= round(exp(beta-qnorm(0.975)*se),2),
                               Upper= round(exp(beta+qnorm(0.975)*se),2),
                               est=paste0(HR, " [",Lower, ", ", Upper,"]"))
meta.fixed <- meta.F(est_prom, est_amsls, method = "fixed")
colnames(meta.fixed)[1]="cov"
meta.fixed = left_join(meta.fixed, est_prom[,c("Covariate","cov","est")], by="cov")
meta.fixed$PROMOTE=meta.fixed$est
meta.fixed = left_join(meta.fixed, est_amsls[,c("cov","est")], by="cov")
meta.fixed$AMSLS=meta.fixed$est.y




##Make forest plots of results

# 1. Forest Plot
library(ggplot2)
library(dplyr)
library(metafor)

res<-meta.fixed %>% select(-c(est.x, est.y))
res$Outcome<-as.factor(ifelse(res$Transition=="1-2", "Mild disability",
                              ifelse(res$Transition%in%c("1-3", "2-3"), "Moderate disability",
                                     ifelse(res$Transition%in%c("1-4", "2-4", "3-4"), "Gait disability",
                                            ifelse(res$Transition%in%c("1-5", "2-5", "3-5", "4-5"), "Early cane", NA)))))
nm<-names(res)[3:11]
res<-res %>% mutate(across(all_of(nm), as.numeric))
res$Covariate <- sub("_[0-9]+$", "", res$Covariate)
res<-data.frame(res)
res=res %>% mutate(META=paste0(HazardRatio_Meta ," [",LogHR_LowerCI ,", ",LogHR_UpperCI ,"]") )
write.csv(res,"Results/meta_results.csv")
res=res[res$Transition%in%c("1-4", "1-5", "2-4", "2-5", "3-4", "3-5", "4-5" ),]
res=res[order(res$Covariate,res$Outcome,res$Transition, decreasing = T),]
res$Covariate=factor(res$Covariate, levels =  c("Younger adults on higher-efficacy DMT vs untreated",
                                                "Younger adults on standard-efficacy DMT vs untreated",
                                                "Middle-aged adults on higher-efficacy DMT vs untreated",
                                                "Middle-aged adults on standard-efficacy DMT vs untreated",
                                                "Older adults on higher-efficacy DMT vs untreated" ,
                                                "Older adults on standard-efficacy DMT vs untreated"     ))
res=res %>% arrange(desc(Covariate))


custom_forest_plot <- function(data, cex.hr, cex.xp, atx, xlim=c(-4.0, 3.2), rectp=3, covp=10, txtcovp_wd=10, stdp=c(1,1), trp=1, hrp=1,  outp=1,alpha=0.2) {
  tryCatch({
    # Ensure necessary columns exist
    required_columns <- c("Covariate", "Transition","PROMOTE", "AMSLS",   "HazardRatio_Meta", "LogHR_Meta", "LogHR_LowerCI", "LogHR_UpperCI", "Outcome", "LogHR_Meta_SE")
    if (!all(required_columns %in% colnames(data))) {
      stop("Data frame must contain required columns: ", paste(required_columns, collapse = ", "))
    }
    
    # Identify groups for each covariate
    covariate_groups <- rle(as.character(data$Covariate))
    group_starts <- cumsum(c(1, head(covariate_groups$lengths, -1)))
    group_lengths <- covariate_groups$lengths
    
    # Prepare display columns
    data$HR_CI_Display <- paste0(
      round(data$HazardRatio_Meta, 2), " [",
      round(data$LogHR_LowerCI, 2), ", ",
      round(data$LogHR_UpperCI, 2), "]"
    )
    data$I2_Outcome <-data$Outcome
    
    # Get row positions
    row_positions <- seq_len(nrow(data))
    
    # Generate the forest plot without the "Study" header
    metafor::forest(
      annotate = FALSE,
      x = data$LogHR_Meta,
      sei = data$LogHR_Meta_SE,
      slab = NA,               # No slab to prevent the "Study" header
      rows = row_positions,
      xlab = "Hazard Ratio",
      refline = 0,
      atransf = exp,
      at = atx,
      digits = c(2, 2),
      cex = cex.hr,
      pch = 15,
      xlim = xlim,
      lwd = 1.2,
      col = "darkgreen",
      header = FALSE           # Suppress the "Study" label
    )
    
    # Add column headers
    text(-covp, max(row_positions) + 2, "Covariate", pos = 4, font = 2, cex = cex.xp)
    text(-trp-0.1, max(row_positions) + 2, "Transition", pos = 4, font = 2, cex=cex.xp)

    text(hrp, max(row_positions) + 2, "HR [95% CI]", pos = 2, font = 2, cex=cex.xp)
    
    
    # Add separators and display values aligned to row positions
    for (i in row_positions) {
      if (i %in% group_starts) {
        lines(x = xlim, y = rep(i - 0.5, 2), lty = "solid", col = "black")
      }
      
      # Display covariate name only for the first row of each group
      #if (i %in% group_starts) {
      # text(-3.0, i, data$Covariate[i], pos = 4, cex = cex.xp, font = 2)
      #}
      # Display covariate name in the middle of the group
      for (g in seq_along(group_starts)) {
        group_center <- group_starts[g] + floor(group_lengths[g] / 2) - 1
        if (i == group_center) {
          wrapped_label <- paste(strwrap(data$Covariate[i], width = txtcovp_wd), collapse = "\n")
          text(-covp, i, wrapped_label, pos = 4, cex = cex.xp, font = 2)
        }
      }
      # Display other columns
      text(-trp, i, data$Transition[i], pos = 4, cex = cex.xp)
    
      text(hrp, i, data$HR_CI_Display[i], pos = 2, cex = cex.xp)
     
    }
    
    # Add background shading for rows with transparency
    for (i in row_positions) {
      if (i %% 2 == 0) {
        rect(
          -rectp, i - 0.5, 3.2, i + 0.5,
          col = adjustcolor("#297B52", alpha.f = alpha),  # 0.5 is 50% transparency
          border = NA
        )
      }
    }
  }, error = function(e) {
    message("Error in forest plot generation: ", e$message)
  })
}


plot.new()
tiff(filename = "Results/DMT by age2.tiff",width = 14, height = 8,units = "in",pointsize = 10,res = 300, compression = "lzw")
atx<-round(seq(from=min(log(range(res$HazardRatio_Meta)))+0.1,to=max(log(range(res$HazardRatio_Meta)))+0.5,by=0.3),2)
custom_forest_plot(res,atx = atx, cex.hr=1.0, rectp=3.5,cex.xp = 0.8, xlim=c(-4.5, 2.2),covp=4.5, txtcovp_wd=35, stdp=c(3.4, 2.9), trp=2.4, hrp=1.6,  outp=1.7, alpha=0.23)

plot.new()
tiff(filename = "Results/DMT by age.tiff",width = 10, height = 8,units = "in",pointsize = 10,res = 300, compression = "lzw")
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

res.trt =res.trt %>% mutate(AgeGroup2=ifelse(AgeGroup=="Middle-aged","Young",ifelse(AgeGroup=="Young","Middle-aged","Older")),
                            AgeGroup2=factor(AgeGroup2, levels=c("Young","Middle-aged","Older")) )
p1<-ggplot(res.trt, aes(x = Duration, y = HR, color = AgeGroup2, linetype = Delay)) +
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
ggsave("Results/overall_dmt_effects.png",width = 12, height = 14,units = "in")
if(1){plot.new()
  tiff(filename = "Results/overall_dmt_effects.tiff",width = 12, height = 14,units = "in",pointsize = 14,res = 300, compression = "lzw")
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


#load data
load("Results/amsls.fit.250430.rdata")
#load("amsls.fit.rdata")
load("Results/promote.fit.rdata")

source("utility-functions.R")

msdat=read.csv("../DATA/2025/MSDAT_20250507.csv") #1070 ID and 3803 PDDS
summary(msdat$obstime)
msdat = msdat %>% group_by(id) %>% 
  mutate(visit = 1:n(),
         dd = as.numeric(as.Date(reviewdate[visit==1])-as.Date(datefde)))
msdat= msdat %>% mutate(dd=ifelse(dd<0, 0, dd), dd=log(dd/365.25+1))

dmt=read.csv("../DATA/2025/CleanDMT.csv")
colnames(dmt)[1]="id"
dmt=dmt %>%
  group_by(id) %>% 
  mutate(trt_start = min(as.Date(dmtstartdate))) %>%
  distinct(id, trt_start)
msdat2=left_join(msdat, dmt, by="id")
msdat2=msdat2 %>% mutate(tage= round(as.numeric((as.Date(trt_start)-as.Date(dob))/365.25),6),
                         rage=round(as.numeric((as.Date(reviewdate)-as.Date(dob))/365.25),6),#review age
                         oage= round(as.numeric((as.Date(datefde)-as.Date(dob))/365.25),6), #oage
                         tage= ifelse(is.na(tage)==T, NA, ifelse(tage<oage, oage, tage)),
                         tage=ifelse(is.na(tage)==T, rage, tage),
                         tagecat= ifelse(tage >50, "2", ifelse(tage>30,"1", "0")),
                         tagecat = factor(tagecat, levels = c("0","1","2")))

msdat2=msdat2 %>% mutate(tdelay=round(tage-oage, 6),
                         tdelay= ifelse(tdelay<=0, 0, tdelay),
                         tdelay=log(tdelay+1))


msdat=msdat2
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



