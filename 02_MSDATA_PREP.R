library(data.table)
library(dplyr)
library(msm)
library(reshape)
library(readxl)
library(lubridate)
library(purrr)
library(tidyr)

#setwd("C:\\Users\\pumpk\\OneDrive - University of Pittsburgh\\XiaLab\\MSM\\2025\\")
setwd("~/OneDrive - University of Pittsburgh/XiaLab/MSM/2025/")
##### Prepare demographic data#####
demo=read.csv("../DATA/2025/CleanDemographics.csv")
demo2=demo %>% mutate(sex=ifelse(subject_sex==1, 1,0 ),
                      race_eth=ifelse(race_cmb=="White" & eth_cmb=="Non-Hispanic or Latino", 0,1) ,# 0- non-hispanic white, 1-other
                      progstat=ifelse(subtype_enroll %in% c("PPMS","SPMS"), 1,0 ),# 0: CIS/RIS/RRMS; 1: PMS
                      datediag= as.Date(date_msdx),
                      datefde= as.Date(date_firstsx),
                      datediag = ifelse(is.na(datediag) & !is.na(datefde), as.character(datefde), as.character(datediag)),
                      datefde = ifelse(is.na(datefde) & !is.na(datediag), as.character(datediag), as.character(datefde)),
                      dage = round(as.numeric(as.Date(datediag) - as.Date(dob)) / 365.25, 2),
                      oage = round(as.numeric(as.Date(datefde) - as.Date(dob)) / 365.25, 2),
                      tsymp = abs(round(as.numeric(as.Date(datediag) - as.Date(datefde)) / 365.25, 2))
) %>%
  select(id,dob, sex, race_eth, progstat, datediag, datefde, dage, oage, tsymp )
demo2= demo2 %>% filter(is.na(tsymp)==F & oage>=18)
dim(demo2) #2391 ppts have demographic data


##### prepare pdds data #####
disab = read.csv("../DATA/2025//Clean_PDDS_20241001.csv")
colnames(disab)[3]="pdds"
dim(disab); length(unique(disab$id)) #2847 ppts with 23592 pdds
disab = disab %>%
  arrange(id, reviewdate) %>%
  group_by(id) %>%
  mutate(
    lagrev = lag(reviewdate),
    t_interval = as.numeric(as.Date(reviewdate) - as.Date(lagrev)) / 365.25,
    t_interval=ifelse(is.na(t_interval)==T, 0, t_interval),
    years=cumsum(t_interval),
    state = case_when(
      pdds == 0 ~ 1,
      pdds == 1 ~ 2,
      pdds == 2 ~ 3,
      pdds == 3 ~ 4,
      TRUE ~ 5
    )
  ) 
summary(disab$t_interval)  # median is 0.24 years

disab= disab %>%
  filter(!is.na(pdds) & (t_interval == 0 | t_interval >= 0.25))   # drop pdds that are <3 month apart 
dim(disab); length(unique(disab$id)) #2847 ppts with 13982 pdds


disab2 <- disab %>%
  group_by(id) %>%
  mutate(
    maxstate = max(state, na.rm = TRUE), 
    absorb_year = min(years[state == 5], na.rm = TRUE)
  )

# Remove observations occured after absorb_time
disab2 <- disab2 %>%
  filter(years<=absorb_year)
dim(disab2); length(unique(disab2$id)) #2847 ppts with 10280 pdds

# Arrange and remove consecutive duplicate states
disab2 <- disab2 %>%
  arrange(id, years) %>%
  group_by(id) %>%
  mutate(lagstate = lag(state)) %>%
  filter(state != lagstate | is.na(lagstate)) 
dim(disab2); length(unique(disab2$id)) #2847 ppts with 5956 pdds

disab2 = disab2 %>%
  group_by(id) %>%
  mutate(n = n()) %>%
  filter(n > 1)                                  
dim(disab2); length(unique(disab2$id)) #1247 ppts with 4356 pdds
disab2=disab2 %>% filter(id %in% demo2$id)
dim(disab2); length(unique(disab2$id)) #1070 ppts with 3803 pdds

disab2 = disab2 %>%
  arrange(id, reviewdate) %>%
  group_by(id) %>%
  mutate(
    lagrev = lag(reviewdate),
    t_interval = as.numeric(as.Date(reviewdate) - as.Date(lagrev)) / 365.25,
    t_interval=ifelse(is.na(t_interval)==T, 0, t_interval),
    years=cumsum(t_interval),
    obstime = cumsum(log(t_interval + 1)),
    n=n(), 
    state = case_when(
      pdds == 0 ~ 1,
      pdds == 1 ~ 2,
      pdds == 2 ~ 3,
      pdds == 3 ~ 4,
      TRUE ~ 5
    )
  ) %>%
  select(id, reviewdate, pdds, state, lagrev, lagstate, t_interval, years, n, maxstate, obstime)
summary(disab2$obstime) # min =0, max=5.8, median = 0.9
disab=disab2
dim(disab); length(unique(disab$id)) #1070 ppts with 3803 pdds

Q <- statetable.msm(state = state, subject = id, data = disab); Q
Q<-round(Q/rowSums(Q),6)
Q=rbind(Q,"5"=rep(0,5)); Q


######################################
## Step 3: Read and Process DMT Data ##
######################################

# Read and clean DMT data
dmt <- read.csv("../DATA/2025//CleanDMT.csv")
colnames(dmt)[2]="dmt"
colnames(dmt)[4]="dmtenddate"
dmt2 <- dmt %>%
  mutate(across(c(dmtstartdate, dmtenddate), as.Date))%>%
  arrange(id_participant, dmtstartdate, dmtenddate ) 
##replace DMT names with common names
# Define a lookup table for common names
dmt_common_names <- c(
  "azathioprine"="Azathioprine",
  "cladribine" = "Cladribine",
  "cyclophosphamide"="Cyclophosphamide",
  "daclizumab"="Daclizumab",
  "dimethyl fumarate" = "Dimethyl Fumarate",
  "dioroximel fumarate"="Dioroximel Fumarate",
  "fingolimod" = "Fingolimod",
  "glatiramer acetate" = "Glatiramer Acetate",
  "interferon beta-1a"="Interferon beta-1a",
  "peginterferon beta-1a"="Interferon beta-1a",
  "interferon beta-1b" = "Interferon beta-1b",
  "mitoxantrone" = "Mitoxantrone",
  "methotrexate"="Methotrexate",
  "monomethyl fumarat"="Monomethyl Fumarat",
  "natalizumab" = "Natalizumab",
  "ocrelizumab" = "Ocrelizumab",
  "ofatumumab" = "Ofatumumab",
  "ozanimod"="Ozanimod",
  "rituximab"="Rituximab",
  "siponimod" = "Siponimod",
  "teriflunomide" = "Teriflunomide"
)
dmt2 <- dmt2 %>%
    mutate(dmt_cname = ifelse(dmt %in% names(dmt_common_names), dmt_common_names[dmt], "other"))
# Define high/low efficacy class
high_efficacy_dmts <- c("Cladribine","Cyclophosphamide","Daclizumab","Mitoxantrone",
                        "Natalizumab","Ocrelizumab","Ofatumumab","Rituximab","Siponimod")
low_efficacy_dmts <- c("Dimethyl Fumarate","Dioroximel Fumarate","Fingolimod","Glatiramer Acetate",
                       "Interferon beta-1b", "Interferon beta-1a","Methotrexate","Ozanimod","Teriflunomide","Azathioprine")
                       


dmt_daily <- dmt2 %>%
  mutate(

    # Classify DMT
    dmt_class = case_when(
      dmt_cname %in% high_efficacy_dmts ~ "high",
      dmt_cname %in% low_efficacy_dmts ~ "low",
      TRUE ~ NA_character_
    )
  ) %>%
  # Remove problematic or unclassified entries
  filter(!is.na(dmt_class), !is.na(dmtstartdate), !is.na(dmtenddate)) %>%
  # Remove rows where end < start (just in case)
  filter(dmtenddate >= dmtstartdate) %>%
  rowwise() %>%
  do(data.frame(
    id = .$id,
    date = seq(.$dmtstartdate, .$dmtenddate, by = "day"),
    dmt_class = .$dmt_class
  )) %>%
  ungroup()
# Compute cumulative exposure up to lagrev (start of interval)
cum_exposure_start <- disab2 %>%
  filter(!is.na(lagrev)) %>%
  select(id, lagrev) %>%
  left_join(dmt_daily, by = "id") %>%
  filter(date < lagrev) %>%
  group_by(id, lagrev, dmt_class) %>%
  summarise(days = n(), .groups = "drop") %>%
  pivot_wider(names_from = dmt_class, values_from = days, names_prefix = "cum_") %>%
  mutate(
    cum_high = round(cum_high / 365.25, 2),
    cum_low = round(cum_low / 365.25, 2)
  ) %>%
  rename_with(~ paste0("ddmt_", gsub("cum_", "", .x)), starts_with("cum_"))

# Merge into disab
disab2 <- disab2 %>%
  left_join(cum_exposure_start, by = c("id", "lagrev")) %>%
  mutate(
    ddmt_high=replace_na(ddmt_high,0),
    ddmt_low=replace_na(ddmt_low,0),
    ddmt_high=log(ddmt_high+1),
    ddmt_low=log(ddmt_low+1)
  )

# Add binary treatment at start of interval
dmt_start <- dmt2 %>%
  mutate(
    dmt_class = case_when(
      dmt_cname %in% high_efficacy_dmts ~ "high",
      dmt_cname %in% low_efficacy_dmts ~ "low",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(dmt_class))



dmt_active_name <- dmt2 %>%
  filter(!is.na(dmt_cname))
colnames(dmt_active_name)[1]="id"
# Merge into disab
disab2 <- disab2 %>%
  rowwise() %>%
  mutate(
    dmt_name = dmt_active_name %>%
      filter(id == cur_data()$id,
             dmtstartdate <= lagrev,
             dmtenddate >= lagrev) %>%
      pull(dmt_cname) %>%
      dplyr::first(default = NA_character_)
  ) %>%
  ungroup()

dim(disab2)

# Check if on treatment at lagrev
disab2 <- disab2 %>%
  rowwise() %>%
  mutate(
    on_high = as.integer(any(dmt_start$id == id & dmt_start$dmt_class == "high" &
                               dmt_start$dmtstartdate <= lagrev & dmt_start$dmtenddate >= lagrev)),
    on_low = as.integer(any(dmt_start$id == id & dmt_start$dmt_class == "low" &
                              dmt_start$dmtstartdate <= lagrev & dmt_start$dmtenddate >= lagrev)),
    on_high=replace_na(on_high,0),
    on_low=replace_na(on_low,0)
  ) %>%
  ungroup()

msdat = disab2 %>%
    left_join(demo2, by = "id") 
write.csv(msdat, "../DATA/2025/MSDAT_20250507.csv", row.names = F)

# 
# #### Merge pdds and DMT #####
# dmt=read.csv("../DATA/2025//CleanDMT.csv")
# colnames(dmt)[1:3]=c("id","dmtstartdate" ,"dmtenddate") 
# disab2 <- disab %>%
#   left_join(dmt, by = "id", relationship = "many-to-many") %>%
#   arrange (id, reviewdate) %>%
#   group_by(id) %>%
#   mutate(first_review=min(reviewdate),
#          last_review=max(reviewdate),
#          first_dmt_start=min(dmtstartdate),
#          last_dmt_end=max(dmtenddate) )
# 
# disab2_a=disab2 %>% 
#   filter((as.Date(first_review)<as.Date(first_dmt_start) &
#             as.Date(last_review)<as.Date(first_dmt_start))|
#            (is.na(first_dmt_start)==T) ) 
# id_a= unique(disab2_a$id); length(id_a)  ##83
# dmt2_a=disab2_a %>% select(id, first_review, last_review) %>% distinct(id, .keep_all = T)
# dmt2_a=dmt2_a %>% mutate(dmtstartdate = first_review, dmtenddate = last_review,
#                          Efficacy=0, 
#                          dur=as.numeric(as.Date(dmtenddate)-as.Date(dmtstartdate))) %>%
#   select(id, dmtstartdate, dmtenddate, Efficacy, dur)
# dmt2_a=dmt2_a %>% mutate(dmtstartdate=as.character(dmtstartdate),
#                          dmtenddate=as.character(dmtenddate))
# 
# disab2_b=disab2 %>% 
#   filter( as.Date(first_review)<as.Date(first_dmt_start) &
#             as.Date(last_review)>=as.Date(first_dmt_start) &
#             as.Date(last_review)<=as.Date(last_dmt_end) &
#             is.na(first_dmt_start)==F ) 
# id_b= unique(disab2_b$id); length(id_b)  #46
# dmt2_b=disab2_b %>% select(id, first_review, first_dmt_start) %>% distinct(id, .keep_all = T)
# dmt2_b=dmt2_b %>% mutate(dmtstartdate = first_review, dmtenddate = as.Date(first_dmt_start)-1,
#                          Efficacy=0, 
#                          dur=as.numeric(as.Date(dmtenddate)-as.Date(dmtstartdate))) %>%
#   select(id, dmtstartdate, dmtenddate, Efficacy, dur)
# dmt2_b=dmt2_b %>% mutate(dmtstartdate=as.character(dmtstartdate),
#                          dmtenddate=as.character(dmtenddate))
# 
# 
# disab2_c=disab2 %>% 
#   filter( as.Date(first_review)<as.Date(first_dmt_start) &
#             as.Date(last_review)>as.Date(last_dmt_end) &
#             is.na(first_dmt_start)==F ) 
# id_c= unique(disab2_c$id); length(id_c)  #6
# dmt2_c1=disab2_c %>% select(id, first_review, first_dmt_start) %>% distinct(id, .keep_all = T)
# dmt2_c1=dmt2_c1 %>% mutate(dmtstartdate = first_review, dmtenddate = as.Date(first_dmt_start)-1,
#                            Efficacy=0, 
#                            dur=as.numeric(as.Date(dmtenddate)-as.Date(dmtstartdate))) %>%
#   select(id, dmtstartdate, dmtenddate, Efficacy, dur)
# dmt2_c1=dmt2_c1 %>% mutate(dmtstartdate=as.character(dmtstartdate),
#                            dmtenddate=as.character(dmtenddate))
# dmt2_c2=disab2_c %>% select(id, last_dmt_end, last_review) %>% distinct(id, .keep_all = T)
# dmt2_c2=dmt2_c2 %>% mutate(dmtstartdate = as.Date(last_dmt_end)+1, dmtenddate = last_review,
#                            Efficacy=0, 
#                            dur=as.numeric(as.Date(dmtenddate)-as.Date(dmtstartdate))) %>%
#   select(id, dmtstartdate, dmtenddate, Efficacy, dur)
# dmt2_c2=dmt2_c2 %>% mutate(dmtstartdate=as.character(dmtstartdate),
#                            dmtenddate=as.character(dmtenddate))
# dmt2_c=rbind(dmt2_c1, dmt2_c2)
# 
# disab2_d=disab2 %>% 
#   filter( as.Date(first_review)>=as.Date(first_dmt_start) &
#             as.Date(first_review)<=as.Date(last_dmt_end) &
#             as.Date(last_review)>=as.Date(first_dmt_start) &
#             as.Date(last_review)<=as.Date(last_dmt_end) &
#             is.na(first_dmt_start)==F ) 
# id_d= unique(disab2_d$id); length(id_d)  #782
# 
# disab2_e=disab2 %>% 
#   filter( as.Date(first_review)>=as.Date(first_dmt_start) &
#             as.Date(first_review)<=as.Date(last_dmt_end) &
#             as.Date(last_review)>as.Date(last_dmt_end)  &
#             is.na(first_dmt_start)==F ) 
# id_e= unique(disab2_e$id); length(id_e)  #169
# dmt2_e=disab2_e %>% select(id, last_dmt_end, last_review) %>% distinct(id, .keep_all = T)
# dmt2_e=dmt2_e %>% mutate(dmtstartdate = as.Date(last_dmt_end)+1, dmtenddate = last_review,
#                          Efficacy=0, 
#                          dur=as.numeric(as.Date(dmtenddate)-as.Date(dmtstartdate))) %>%
#   select(id, dmtstartdate, dmtenddate, Efficacy, dur)
# dmt2_e=dmt2_e %>% mutate(dmtstartdate=as.character(dmtstartdate),
#                          dmtenddate=as.character(dmtenddate))
# 
# 
# disab2_f=disab2 %>% 
#   filter( as.Date(first_review)>as.Date(last_dmt_end) &
#             as.Date(last_review)>as.Date(last_dmt_end)  &
#             is.na(first_dmt_start)==F ) 
# id_f= unique(disab2_f$id); length(id_f)  #43
# dmt2_f=disab2_f %>% select(id, first_review, last_review) %>% distinct(id, .keep_all = T)
# dmt2_f=dmt2_f %>% mutate(dmtstartdate = first_review, dmtenddate = last_review,
#                          Efficacy=0, 
#                          dur=as.numeric(as.Date(dmtenddate)-as.Date(dmtstartdate))) %>%
#   select(id, dmtstartdate, dmtenddate, Efficacy, dur)
# dmt2_f=dmt2_f %>% mutate(dmtstartdate=as.character(dmtstartdate),
#                          dmtenddate=as.character(dmtenddate))
# 
# dmt2=rbind(dmt, dmt2_a, dmt2_b, dmt2_c, dmt2_e, dmt2_f)
# dmt2= dmt2 %>% arrange(id, dmtstartdate)
# disab2 <- disab %>%
#   left_join(dmt2, by = "id", relationship = "many-to-many") %>%
#   arrange (id, reviewdate) %>%
#   mutate(flag= ifelse(as.Date(reviewdate)>=as.Date(dmtstartdate) &
#                         as.Date(reviewdate)<= as.Date(dmtenddate)  ,1, 0))
# 
# disab2_a= disab2 %>% filter(id %in% id_a) %>% 
#   mutate(total_dmt_years=0 )%>% filter(flag==1)
# dim(disab2_a); dim(disab %>% filter(id %in% id_a))
# 
# 
# disab2_b_tmp= disab2 %>% filter(id %in% id_b ) %>%
#   group_by(id) %>%
#   mutate(tmpv=ifelse(reviewdate==max(reviewdate),1,0))%>%
#   filter(tmpv==1) %>% 
#   mutate(dur= as.numeric(pmin(as.Date(dmtenddate), as.Date(reviewdate))-as.Date(dmtstartdate)))%>%
#   filter(dur>=0)%>% 
#   summarise(total_dmt_years=  sum(dur*(Efficacy!=0)))
# 
# disab2_b= disab2 %>% filter(id %in% id_b) %>% 
#   left_join(disab2_b_tmp, by="id")%>% filter(flag==1)
# dim(disab2_b); dim(disab %>% filter(id %in% id_b))
# 
# disab2_c= disab2 %>% filter(id %in% id_c) %>% 
#   group_by(id) %>%
#   mutate(tmpv=ifelse(reviewdate==max(reviewdate),1,0))%>%
#   mutate(total_dmt_years=  sum(dur[tmpv==1]*((Efficacy!=0)[tmpv==1])))%>%
#   filter(flag==1)%>% select(-tmpv)
# dim(disab2_c); dim(disab %>% filter(id %in% id_c))
# 
# disab2_d_tmp1= disab2 %>% filter(id %in% id_d) %>% 
#   group_by(id) %>%
#   mutate(tmpv=ifelse(reviewdate==max(reviewdate),1,0))%>%
#   filter(tmpv==1) %>% 
#   mutate(dur= as.numeric(pmin(as.Date(dmtenddate), as.Date(reviewdate))-as.Date(dmtstartdate)))%>%
#   filter(dur>=0)%>% 
#   summarise(total_dmt_years1=  sum(dur*(Efficacy!=0)))
# disab2_d_tmp2= disab2 %>% filter(id %in% id_d) %>% 
#   group_by(id) %>%
#   mutate(tmpv=ifelse(reviewdate==min(reviewdate),1,0))%>%
#   filter(tmpv==1) %>% 
#   mutate(dur= as.numeric(pmin(as.Date(dmtenddate), as.Date(reviewdate))-as.Date(dmtstartdate)))%>%
#   filter(dur>=0)%>% 
#   summarise(total_dmt_years2=  sum(dur*(Efficacy!=0)))
# disab2_d_tmp=merge(disab2_d_tmp1, disab2_d_tmp2, by="id")
# disab2_d_tmp= disab2_d_tmp %>% mutate(total_dmt_years= total_dmt_years1-total_dmt_years2) %>% select(id, total_dmt_years)
# disab2_d= disab2 %>% filter(id %in% id_d) %>% 
#   left_join(disab2_d_tmp, by="id")%>% filter(flag==1)
# dim(disab2_d); dim(disab %>% filter(id %in% id_d))
# 
# 
# 
# disab2_e_tmp1= disab2 %>% filter(id %in% id_e) %>% 
#   group_by(id) %>%
#   mutate(tmpv=ifelse(reviewdate==max(reviewdate),1,0))%>%
#   filter(tmpv==1) %>% 
#   mutate(dur= as.numeric(pmin(as.Date(dmtenddate), as.Date(reviewdate))-as.Date(dmtstartdate)))%>%
#   filter(dur>=0)%>% 
#   summarise(total_dmt_years1=  sum(dur*(Efficacy!=0)))
# disab2_e_tmp2= disab2 %>% filter(id %in% id_e) %>% 
#   group_by(id) %>%
#   mutate(tmpv=ifelse(reviewdate==min(reviewdate),1,0))%>%
#   filter(tmpv==1) %>% 
#   mutate(dur= as.numeric(pmin(as.Date(dmtenddate), as.Date(reviewdate))-as.Date(dmtstartdate)))%>%
#   filter(dur>=0)%>% 
#   summarise(total_dmt_years2=  sum(dur*(Efficacy!=0)))
# disab2_e_tmp=merge(disab2_e_tmp1, disab2_e_tmp2, by="id")
# disab2_e_tmp= disab2_e_tmp %>% mutate(total_dmt_years= total_dmt_years1-total_dmt_years2) %>% select(id, total_dmt_years)
# disab2_e= disab2 %>% filter(id %in% id_e) %>% 
#   left_join(disab2_e_tmp, by="id")%>% filter(flag==1)
# dim(disab2_e); dim(disab %>% filter(id %in% id_e))
# 
# 
# disab2_f= disab2 %>% filter(id %in% id_f) %>% 
#   mutate(total_dmt_years=0 )%>% filter(flag==1)
# dim(disab2_f); dim(disab %>% filter(id %in% id_f))
# 
# 
# disab2=rbind(disab2_a,disab2_b, disab2_c, disab2_d, disab2_e, disab2_f )
# dim(disab2) #4805, same as disab

#####DMT-related variables #####
# disab2= disab2 %>% arrange(id, reviewdate) %>%
#   group_by(id) %>%
#   mutate(total_years_up_to_review = round(as.numeric(as.Date(reviewdate)-min(as.Date(reviewdate)))/365.25,2),
#          total_dmt_years=round(total_dmt_years/365.25,2),
#          dmt_high=ifelse(Efficacy==2, 1, 0),
#          dmt_low=ifelse(Efficacy==1, 1, 0),
#          ddmt_frac=total_years_up_to_review / total_dmt_years,
#          ddmt_frac= ifelse(total_dmt_years==0, 0, ddmt_frac),
#          ddmt = cumsum(log(ddmt_frac + 1)),
#          ddmt_high=dmt_high*ddmt,
#          ddmt_low=dmt_low*ddmt)
# tmp= disab2 %>% group_by(id) %>%
#   summarise(total_dmt_years= max(total_dmt_years),
#             total_review_years= max(total_years_up_to_review)) # check - total_dmt_years <= total_review_years
# disab3= disab2 %>% select(id, reviewdate, pdds, state, lagrev, lagstate, t_interval, n, maxstate, obstime, 
#                           dmtstartdate, dmtenddate, dmt_high, dmt_low, total_years_up_to_review, total_dmt_years, 
#                           ddmt_frac, ddmt, ddmt_high, ddmt_low)
# summary(disab3$ddmt_high)

# 
# #### Merge demo data with (PDDS and DMT) #####
# msdat = disab3 %>%
#   left_join(demo2, by = "id") 
# write.csv(msdat, "../DATA/2025/MSDAT.csv", row.names = F)
