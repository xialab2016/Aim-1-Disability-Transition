library(dplyr)
library(stringr)
library(tidyr)

setwd("~/OneDrive - University of Pittsburgh/XiaLab/RedCap/20241001/")
######1. PDDS from PROMOTE ######
database <- read.csv("PROMOTEDatabase_DATA_2024-10-01_1545.csv") #Export of the whole database project
database <- database %>% filter(str_detect(id_participant, "PRT"))
database <- database %>% filter(!str_detect(id_participant, "_x"))
colnames(database)[which(grepl("pdds",colnames(database), ignore.case = T)==T)] #names of PDDS-related
######1.1 PDDS history######
pdds_hx <-  database[,c("id_participant", "pdds_number", paste0("pdds_date_",1:25),paste0("pdds_value_",1:25))] #maximum 25 records
pdds_hx <- pdds_hx %>% filter(!is.na(pdds_number))
table(pdds_hx$pdds_number, useNA = "ifany")
colnames(pdds_hx)
pdds_hxScore <- pdds_hx %>% gather(scoreNumber, score,28:52 )
pdds_hxScore <- pdds_hxScore[,c("id_participant", "scoreNumber", "score")]
pdds_hxDate <- pdds_hx %>% gather(dateNumber, date, 3:27)
pdds_hxDate <- pdds_hxDate[,c("id_participant", "dateNumber", "date")]
pdds_hx <- cbind(pdds_hxDate, pdds_hxScore)
pdds_hx <- pdds_hx[,c("id_participant", "date", "score")]
pdds_hx <- pdds_hx %>% filter(!is.na(score))

table(pdds_hx$score, useNA = "ifany")
summary(as.Date(pdds_hx$date)) #check date
pdds_hx=pdds_hx %>% filter(is.na(as.Date(date))==F)
which(as.Date(pdds_hx$date)>Sys.Date() )
pdds_hx$date[which(as.Date(pdds_hx$date)>Sys.Date() )]
pdds_hx$date[3599]="2018-05-28"
pdds_hx$date[6611]="2019-06-13"
pdds_hx$date[6650]="2023-01-31"
dim(pdds_hx); length(unique(pdds_hx$id_participant)) #13,583 PDDS; 2080 ID
#####1.2 PDDS Blood sample ######
pdds_blood <- database[,c("id_participant", "blood1_date", "blood1_pddsv2")]
pdds_blood <- pdds_blood %>% filter(!is.na(blood1_pddsv2))
names(pdds_blood)[2] <- 'date'
names(pdds_blood)[3] <- 'score'
table(pdds_blood$score, useNA = "ifany")
summary(as.Date(pdds_blood$date)) #check date
###### 1.3 PDDS from Stool Sample######
pdds_stool <- database[,c("id_participant", "stool0a_collection_date", "stool0b_collection_date", "collection1a_date", "collection1b_date", "stool_pddsv2_date", "stool_pddsv2")]

paste_noNA <- function(x,sep=", ") {
  gsub(", " ,sep, toString(x[!is.na(x) & x!="" & x!="NA"] ) ) }
sep=", "
pdds_stool$date <- apply( pdds_stool[ , c(2:6) ] , 1 , paste_noNA , sep=sep)

pdds_stool <- pdds_stool[,c("id_participant", "date", "stool_pddsv2")]
pdds_stool <- pdds_stool %>% filter(!is.na(stool_pddsv2))
pdds_stool$date <- sub("\\,.*", "", pdds_stool$date)
pdds_stool <- pdds_stool %>% filter(str_detect(date, "20"))
names(pdds_stool)[3] <- 'score'
table(pdds_stool$score, useNA = "ifany")
summary(as.Date(pdds_stool$date)) #check date
###### 1.4 PDDS of CSF sample ######
pdds_csf <- database[,c("id_participant", "csf_date", "csf_pddsv2")]
pdds_csf <- pdds_csf %>% filter(!is.na(csf_pddsv2))
names(pdds_csf)[2] <- 'date'
names(pdds_csf)[3] <- 'score'
table(pdds_csf$score, useNA = "ifany")
summary(as.Date(pdds_csf$date)) #check date

######2. PDDS from Legacy ######
legacy <- read.csv("PROMOTELegacyPrimary_DATA_2024-10-01_1551.csv") #Export of the whole PQ/SRO Legacy project
colnames(legacy)[which(grepl("pdds",colnames(legacy), ignore.case = T)==T)] #names of PDDS-related
pddsLeg <- legacy[,c("redcap_survey_identifier", "promote_self_reported_outcomes_timestamp", "pdds")]
names(pddsLeg)[2] <- 'date'
pddsLeg <- pddsLeg %>% filter(str_detect(date, "20"))
pddsLeg <- pddsLeg %>% filter((!str_detect(pdds, "U")))
pddsLeg=pddsLeg %>% filter(redcap_survey_identifier!="")
names(pddsLeg)[1] <- 'id_participant'
pddsLeg$date=as.Date(pddsLeg$date)
names(pddsLeg)[3] <- 'score'
table(pddsLeg$score, useNA = "ifany")
summary(as.Date(pddsLeg$date)) 

######3. PDDS from SRO & Update- Long1&2 ######
longitudinal1<- read.csv("PROMOTELongitudinalQ_DATA_2024-10-01_1549.csv") 
longitudinal2<- read.csv("PROMOTELongitudinalQ_DATA_2024-10-01_1648.csv") #Combined longitudinal projects
longitudinal=rbind(longitudinal1, longitudinal2)
longitudinal <- longitudinal %>% filter(!str_detect(id_participant_l, "_x"))
colnames(longitudinal)[which(grepl("pdds",colnames(longitudinal), ignore.case = T)==T)] #names of PDDS-related
######3.1 from longitudinal ######
pddsSROUpd <- longitudinal[,c("id_participant_l", "promote_self_reported_outcomes_timestamp", "pddsv2", "update_questionnaire_timestamp", "pddsv2_4d3bcb")]
pddsSROUpdDate <- pddsSROUpd %>% gather(dateNumber, date, 2,4)
pddsSROUpdScore <- pddsSROUpd %>% gather(scoreNumber, score, 3,5)
pddsSROUpd <- cbind(pddsSROUpdDate, pddsSROUpdScore)
names(pddsSROUpd)[1] <- 'id_participant'

pddsSROUpd <- pddsSROUpd[,c("id_participant", "date", "score")]
pddsSROUpd <- pddsSROUpd %>% filter(!is.na(score))
pddsSROUpd <- pddsSROUpd %>% filter(str_detect(date, "20"))
pddsSROUpd$date=as.Date(pddsSROUpd$date)
table(pddsSROUpd$score, useNA = "ifany")
summary(as.Date(pddsSROUpd$date)) 
######3.2 from SRO ######
pddsOut <- longitudinal[,c("id_participant_l", "outcome_measures_timestamp", "pddsv2_outcomes")]
names(pddsOut)[2] <- 'date'
pddsOut <- pddsOut %>% filter(str_detect(date, "20"))
names(pddsOut)[1] <- 'id_participant'
names(pddsOut)[3] <- 'score'
pddsOut$date=as.Date(pddsOut$date)
table(pddsOut$score, useNA = "ifany")
summary(as.Date(pddsOut$date)) 

###### 4 PDDS from social network legacy #####
snq=read.csv("PROMOTELegacySocialN_DATA_2024-10-01_1643.csv")
pddsSNQ <- snq[,c("redcap_survey_identifier", "social_network_questionnaire_timestamp","social_pdds")]
pddsSNQ <- pddsSNQ %>% filter(str_detect(social_network_questionnaire_timestamp, "20"))
pddsSNQ <- pddsSNQ %>% filter(str_detect(redcap_survey_identifier, "PRT"))
names(pddsSNQ)[1] <- 'id_participant'
names(pddsSNQ)[2] <- 'date'
names(pddsSNQ)[3] <- 'score'
pddsSNQ$date=as.Date(pddsSNQ$date)
table(pddsSNQ$score, useNA = "ifany")
pddsSNQ=pddsSNQ %>% filter(is.na(score)==F)
summary(as.Date(pddsSNQ$date)) 

###### 5 PDDS from sensor#####
sensor=read.csv("PROMOTESensor_DATA_2024-10-01_1644.csv")
pddsSensor <- sensor[,c("record_id", "monthlyq_timestamp", "pddsv2")]
pddsSensor <- pddsSensor %>% filter(!str_detect(record_id, "_0"))
pddsSensor <- pddsSensor %>% filter(str_detect(monthlyq_timestamp, "20"))
pddsSensor <- pddsSensor %>% filter(!is.na(pddsv2))
names(pddsSensor)[1] <- 'id_participant'
names(pddsSensor)[2] <- 'date'
names(pddsSensor)[3] <- 'score'
pddsSensor$date=as.Date(pddsSensor$date)
table(pddsSensor$score, useNA = "ifany")
summary(as.Date(pddsSensor$date)) 
###### 5 PDDS from covid long#####
covid=read.csv("COVID19LongitudinalQ_DATA_2024-10-01_1647.csv")
pddsCovid <- covid[,c("id_participant_l", "covid19_baseline_timestamp", "pddsv2_baseline", "covid19_quarterly_timestamp", "pddsv2_3m")]
pdds_CovidDate <- pddsCovid %>% gather(dateNumber, date, 2,4)
pdds_CovidScore <- pddsCovid %>% gather(scoreNumber, score, 3,5)
pddsCovid <- cbind(pdds_CovidDate, pdds_CovidScore)
names(pddsCovid)[1] <- 'id_participant'
pddsCovid <- pddsCovid[,c("id_participant", "date", "score")]
pddsCovid <- pddsCovid %>% filter(!is.na(score))
pddsCovid$date=as.Date(pddsCovid$date)
table(pddsCovid$score, useNA = "ifany")
summary(as.Date(pddsCovid$date)) 
pddsCovid=pddsCovid %>% filter(is.na(date)==F)
###### 6 PDDS from vaccine#####
vaccine=read.csv("COVID19VaccineLongit_DATA_2024-10-01_1648.csv")
pddsVaccine <- vaccine[,c("id_participant_l", "covid19_vaccine_questions_timestamp", "pddsv2_outcomes")]
pddsVaccine <- pddsVaccine %>% filter(str_detect(covid19_vaccine_questions_timestamp, "20"))
pddsVaccine <- pddsVaccine %>% filter(!is.na(pddsv2_outcomes))
names(pddsVaccine)[1] <- 'id_participant'
names(pddsVaccine)[2] <- 'date'
names(pddsVaccine)[3] <- 'score'
pddsVaccine$date=as.Date(pddsVaccine$date)
table(pddsVaccine$score, useNA = "ifany")
summary(as.Date(pddsVaccine$date)) 
###### 7. PDDS from prebiotic ######
prepro=read.csv("PreProbiotic2_DATA_2024-10-01_1645.csv")
pddsPrePro <- prepro[,c("record_id", "promote_id", "v1_survey_timestamp", "v2_survey_timestamp", "v3_survey_timestamp", "v4_survey_timestamp", "v5_survey_timestamp", "msrsr_and_pdds_questionnaire_timestamp", "pdds")]
pddsPrePro <- pddsPrePro %>%
  group_by(record_id) %>%
  mutate(promote_id = first(promote_id))




paste_noNA <- function(x,sep=", ") {
  gsub(", " ,sep, toString(x[!is.na(x) & x!="" & x!="NA"] ) ) }
sep=", "
pddsPrePro$date <- apply( pddsPrePro[ , c(3:8) ] , 1 , paste_noNA , sep=sep)

pddsPrePro <- pddsPrePro %>% filter(!is.na(pdds))
names(pddsPrePro)[2] <- 'id_participant'
names(pddsPrePro)[9] <- 'score'
pddsPrePro <- pddsPrePro[,c("id_participant", "date", "score")]
pddsPrePro$date=as.Date(pddsPrePro$date)
table(pddsPrePro$score, useNA = "ifany")
summary(as.Date(pddsPrePro$date)) 
pddsPrePro=pddsPrePro %>%filter(is.na(date)==F)
###### 8. PDDS from pasc ######
pasc=read.csv("PASCSurvey_DATA_2024-10-01_1645.csv")
pasc <- pasc %>% filter(str_detect(id_participant, "PRT"))

pddsPASC <- pasc[,c("id_participant", "covid_date1", "pdds_precovid", "postacute_sequelae_of_sarscov2_timestamp", "pdds_current")]
pddsPASC_date <- pddsPASC %>% gather(dateType, date, 2,4)
pddsPASC_date <- pddsPASC_date[,c(1,4,5)]

pddsPASC_score <- pddsPASC %>% gather(scoreType, score, 3,5)
pddsPASC_score <- pddsPASC_score[,c(1,4,5)]

pddsPASC2 <- cbind(pddsPASC_date, pddsPASC_score)
pddsPASC2 <- pddsPASC2[,c(1,3,6)]
pddsPASC2$date=as.Date(pddsPASC2$date)
table(pddsPASC2$score, useNA = "ifany")
summary(as.Date(pddsPASC2$date)) 
pddsPASC2=pddsPASC2 %>% filter(is.na(date)==F & is.na(score)==F)

###### 9. commbine all scores #######
pdds_all <- rbind(pdds_hx, pdds_blood, pdds_stool,  pdds_csf,pddsLeg, pddsSROUpd, pddsOut, pddsSNQ, pddsSensor, pddsCovid, pddsVaccine, pddsPASC2, pddsPrePro)

pdds_all <- pdds_all %>% mutate(flag=grepl("-",date))
tmp1=as.Date(pdds_all$date[which(pdds_all$flag==T)])
tmp2=lubridate::as_date(as.numeric(pdds_all$date[which(pdds_all$flag==F)]))
pdds_all$date2=as.Date("2024-12-31")
pdds_all$date2[which(pdds_all$flag==T)]=tmp1
pdds_all$date2[which(pdds_all$flag==F)]=tmp2
summary(pdds_all$date2)
pdds_all=pdds_all%>% select(-c(date, flag))
colnames(pdds_all)[3]="date"

table(pdds_all$score) # 9 is missing
pdds_all <- pdds_all %>% filter(score!=9)
pdds_all <- unique(pdds_all)

pdds_all$score <- as.numeric(pdds_all$score)
pdds_all$date <- as.Date(pdds_all$date)
pdds_all <- pdds_all %>% filter(str_detect(id_participant, "PRT"))
pdds_all=pdds_all %>% arrange(id_participant, date)
pdds_all=pdds_all %>% group_by(id_participant, date) %>% summarise(score=round(mean(score),0))

dim(pdds_all) #23,592 PDDS
length(unique(pdds_all$id_participant)) ##2847 ids

write.csv(pdds_all, "~/OneDrive - University of Pittsburgh/XiaLab/MSM/DATA/20241001 update/Clean_PDDS_20241001.csv", row.names = F)



######## 10 query enrollment date ########
enrollment <- longitudinal %>% filter(str_detect(redcap_event_name, "year_01"))
enrollment <- enrollment[,c("id_participant_l", "doe")]
names(enrollment)[1] <- 'id_participant'
enrollment <- enrollment %>% filter(str_detect(id_participant, "PRT")) %>% 
  filter(is.na(as.Date(doe))==F)
enrollment=unique(enrollment)
summary(as.Date(enrollment$doe))
dim(enrollment) #2971 ppts enrollment

######## 11 query demographics ########
database <- read.csv("PROMOTEDatabase_DATA_2024-10-01_1545.csv") #Export of the whole database project
database <- database %>% filter(str_detect(id_participant, "PRT"))
database <- database %>% filter(!str_detect(id_participant, "_x"))
demographics <- database %>% filter(str_detect(redcap_event_name, "consent"))
demographics <- demographics[,c("id_participant","dob", "subject_group","subtype_enroll","subject_sex", "race", "ethnicity", "date_firstsx", "date_msdx","height_feet","height_inches","weight")] #3901 ppts

######## 12 merge enrollment date and demographics ########
data=merge(enrollment, demographics, by="id_participant", all=T) #3904 ppts
table(data$subject_group, useNA = "ifany") # most NAs are enrolled in 2024
data=data %>% filter(subject_group ==1)  ##fiter out non-MS ppts;1	MS
table(data$subtype_enroll, useNA = "ifany")
data=data %>% mutate(subtype_enroll=ifelse(is.na(subtype_enroll)==T, "Unknown",
                                           ifelse(subtype_enroll%in%c(1, 2, 99),"CIS/RIS",
                                                  ifelse(subtype_enroll==3, "RRMS",
                                                         ifelse(subtype_enroll==4, "SPMS","PPMS" )))))
table(data$subtype_enroll, useNA = "ifany") #26 unknown MS subtype
summary(as.Date(data$dob))
summary(data)
data$height_feet=as.numeric(data$height_feet)
which(data$height_feet>10)
data$height_feet[1634]=4
data$height_inches[1634]=5
data$weight[1634]=NA
data$height_feet[2157]=5
data$height_inches[2157]=4
data$weight[2157]=NA
which(data$height_inches>10)
data=data %>% mutate(height_inches=ifelse(height_inches<10, height_inches, height_inches-10),
                     height_feet=ifelse(height_inches<10, height_feet, height_feet+1))
which(data$weight<50 | data$weight>500)
data$weight[1869]=NA
data=data %>% mutate(cm=(height_feet+0.1*height_inches)*30.48,
                     kg=weight/2.205,
                     bmi=kg/(cm/100)^2)
summary(data$bmi)
data=data%>% mutate(onset_to_confirm=as.numeric(as.Date(date_msdx)-as.Date(date_firstsx))/365.25)
ind=which(data$onset_to_confirm<0)
data[ind, c("date_firstsx", "date_msdx")] = data[ind, c("date_msdx", "date_firstsx")] 
data=data%>% mutate(onset_to_confirm=as.numeric(as.Date(date_msdx)-as.Date(date_firstsx))/365.25)
which(data$onset_to_confirm<0) #msdx after first_sx
######### 13 query in EHR to make up for some missingness #######
linked=read.csv("~/OneDrive - UPMC/EHR/UPMC_20240827/ms_promote_linked_20240522.csv")
colnames(linked)[2]="id_participant"
data2=data
#race #
table(data$race, useNA="ifany") #8 not sure, 18 N/A
id=data2$id_participant[which(data2$race==99 | is.na(data2$race)==T)]
tmp1=data %>% filter(id_participant %in% id) %>% select(id_participant, race)
tmp2=linked %>% filter(id_participant %in% id) %>% select(id_participant,RACE_TITLE )
tmp=left_join(tmp1, tmp2, by="id_participant")
table(tmp$RACE_TITLE, useNA="ifany")
data2$race_cmb=as.character(data2$race)
data2$race_cmb[data2$id_participant %in% tmp$id_participant]=tmp$RACE_TITLE
table(data2$race_cmb, useNA = "ifany")
data2=data2 %>% mutate(race_cmb=ifelse(is.na(race_cmb)==T, "Unknown",
                                              ifelse(race_cmb %in% c("1","BLACK/AFRICAN AMERICAN"),"Black/African American",
                                                     ifelse(race_cmb=="2", "American Indian or Alaska Native",
                                                            ifelse(race_cmb =="3", "Asian",
                                                                   ifelse(race_cmb%in% c("4","WHITE" ),"White",
                                                                          ifelse(race_cmb =="5", "Multiracial",
                                                                                 ifelse(race_cmb =="6", "Native Hawaiian or Other Pacific Islander",
                                                                                        ifelse(race_cmb %in% c("90","OTHER"),"Other","Unknown")))))))))
table(data2$race_cmb) #16 Unknown

#ethnicity #
table(data$ethnicity, useNA="ifany") #30 not sure, 20 N/A
id=data2$id_participant[which(data2$ethnicity==99 | is.na(data2$ethnicity)==T)]
tmp1=data %>% filter(id_participant %in% id) %>% select(id_participant, ethnicity)
tmp2=linked %>% filter(id_participant %in% id) %>% select(id_participant,ETHNIC_TITLE )
tmp=left_join(tmp1, tmp2, by="id_participant")
table(tmp$ETHNIC_TITLE, useNA="ifany")
data2$eth_cmb=as.character(data2$ethnicity)
data2$eth_cmb[data2$id_participant %in% tmp$id_participant]=tmp$ETHNIC_TITLE
table(data2$eth_cmb, useNA = "ifany")
data2=data2 %>% mutate(eth_cmb=ifelse(is.na(eth_cmb)==T, "Unknown",
                                       ifelse(eth_cmb =="1","Hispanic or Latino",
                                              ifelse(eth_cmb %in% c("2","NON-HISPANIC OR LATINO/A"), "Non-Hispanic or Latino",  "Unknown"))))
table(data2$eth_cmb) #42 Unknown

#first sx and msdx
summary(as.Date(data2$date_firstsx))
summary(as.Date(data2$date_msdx))


tmp=data %>% select(id_participant,date_firstsx, date_msdx ) %>%
  mutate(sx_missing=is.na(as.Date(date_firstsx)),
         dx_missing=is.na(as.Date(date_msdx)))
tmp=left_join(tmp, linked[,c("id_participant","PATIENT_STUDY_ID")], by="id_participant")
which(is.na(tmp$PATIENT_STUDY_ID)==T) ##163 ppts not in EHR
for_ehr_querry=tmp %>% filter((sx_missing==T | dx_missing==T )& is.na(PATIENT_STUDY_ID)==F)
dim(for_ehr_querry) #160 ppts missing msdx and/or sx date
write.csv(for_ehr_querry, "~/OneDrive - University of Pittsburgh/XiaLab/MSM/DATA/20241001 update/msdx_sx_to_querry_in_EHR.csv", row.names = F)


#load("ehr.Rdta")
# ehr2=ehr[which(ehr$patient_num %in%  toq$PATIENT_STUDY_ID),]
# ehr3=ehr2[which(ehr2$feature_id=="PheCode:335"),]
# ehr4=ehr3 %>% mutate(start_date=as.Date(start_date))
# ehr4=ehr4 %>% arrange(patient_num, start_date)
# ehr5=ehr4[!duplicated(ehr4[,c(1,3)]),]
# ehr5=ehr5 %>% group_by(patient_num) %>% mutate(n=1:n())
# colnames(toq)[6]
# dd=merge(ehr5, toq, by="patient_num")
# write.csv(dd,"msdx_sx_querried_in_EHR.csv", row.names=F)

queried=read.csv("~/OneDrive - University of Pittsburgh/XiaLab/MSM/DATA/20241001 update/msdx_sx_querried_in_EHR.csv")
queried=queried %>% group_by(patient_num) %>% mutate(n=1:n())
queried2=queried %>% group_by( patient_num)%>%
  mutate(date_msdx_impute=start_date[n==1],
         date_msdx=ifelse(dx_missing==T,date_msdx_impute, date_msdx ))
summary(as.Date(queried2$date_msdx)) #no missing msdx

queried2=queried2%>% mutate(onset_to_confirm=as.numeric(as.Date(date_msdx)-as.Date(date_firstsx))/365.25)
tmp_id=unique(queried2$id_participant[which(queried2$onset_to_confirm<0)]) #check msdx  before sx
tmp_id
tmp=queried2 %>% filter(id_participant %in% tmp_id) %>% 
  mutate(gap=as.numeric(as.Date(date_firstsx)-as.Date(start_date)))%>%
  group_by(id_participant) %>%
  mutate(date_msdx_impute=start_date[gap==max(gap[gap<=0])],
    date_msdx=date_msdx_impute)

queried2$date_msdx[which(queried2$onset_to_confirm<0)]=tmp$date_msdx
queried2$date_msdx_impute[which(queried2$onset_to_confirm<0)]=tmp$date_msdx_impute
queried2=queried2%>% mutate(onset_to_confirm=as.numeric(as.Date(date_msdx)-as.Date(date_firstsx))/365.25)
unique(queried2$id_participant[which(queried2$onset_to_confirm<0)]) #check msdx  before sx
queried2=queried2 %>% mutate(date_firstsx=ifelse(sx_missing==F,date_firstsx, date_msdx ))
queried2=queried2%>% mutate(onset_to_confirm=as.numeric(as.Date(date_msdx)-as.Date(date_firstsx))/365.25)
unique(queried2$id_participant[which(queried2$onset_to_confirm<0)]) #check msdx  before sx
queried_final=queried2 %>% select(id_participant, date_firstsx, date_msdx)
queried_final=unique(queried_final)
for (i in 1:nrow(queried_final)){
  ind=which(data2$id_participant == queried_final$id_participant[i])
  data2$date_firstsx[ind]=queried_final$date_firstsx[i]
  data2$date_msdx[ind]=queried_final$date_msdx[i]
}

data2=data2%>% mutate(onset_to_confirm=as.numeric(as.Date(date_msdx)-as.Date(date_firstsx))/365.25)
which(data2$onset_to_confirm<0) #msdx after first_sx
summary(as.Date(data2$date_firstsx)) #20 missing sx date
summary(as.Date(data2$date_msdx)) ##26 missing dx date
data2=data2 %>% mutate(date_firstsx=ifelse(is.na(as.Date(date_firstsx))==T, date_msdx, date_firstsx))
data2=data2%>% mutate(onset_to_confirm=as.numeric(as.Date(date_msdx)-as.Date(date_firstsx))/365.25)
which(data2$onset_to_confirm<0) #msdx after first_sx
summary(as.Date(data2$date_firstsx)) #9 missing sx date
summary(as.Date(data2$date_msdx)) ##26 missing dx date (need to drop later)

linked2=linked %>% filter(id_participant %in% data2$id_participant) %>% select(id_participant, PATIENT_STUDY_ID)
data2=data2 %>% left_join(linked2, by="id_participant")
data2=data2%>% mutate(inEHR=ifelse(is.na(as.numeric(PATIENT_STUDY_ID))==T, 0, 1))
data2=data2 %>% select(id_participant, PATIENT_STUDY_ID,inEHR,
                       doe, dob, date_firstsx, date_msdx, onset_to_confirm, 
                       subtype_enroll, subject_sex, race_cmb, eth_cmb, bmi )

write.csv(data2, "~/OneDrive - University of Pittsburgh/XiaLab/MSM/DATA/20241001 update/CleanDemographics.csv", row.names = F)

