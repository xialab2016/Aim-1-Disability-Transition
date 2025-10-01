library(dplyr)
library(reshape)
library(readxl)
library(lubridate)
#trt=read.csv("C:\\Users\\pumpk\\OneDrive - University of Pittsburgh\\XiaLab\\Aim2/Clean data/UPMC/treatment.csv")
trt=read.csv("~/OneDrive - University of Pittsburgh/XiaLab/Aim2/Clean data/UPMC/treatment.csv")

trt=trt %>% filter(is.na(start)== F & is.na(Efficacy)==F & Treatment_GenericName !="None" )
trt$end=NA
trt=trt %>% group_by(id_participant ) %>%
  mutate(end=as.Date(lead(start))-1,
         dur=as.numeric(end-as.Date(start)))%>%
  filter(dur!=0 ) %>%
  mutate(end=as.Date(lead(start))-1,
         dur=as.numeric(end-as.Date(start))) %>%
  filter(is.na(end)==F)

tmp=trt %>%
  group_by(id_participant)  %>% summarise( rle_result = list(rle(Efficacy))  # Apply rle() within each group
  ) %>%
  mutate(
    lengths = map(rle_result, ~ .x$lengths),  # Extract lengths
    values = map(rle_result, ~ .x$values)     # Extract values
  ) %>%
  select(-rle_result) %>% 
  left_join(trt, tmp,by="id_participant")
tmp$group_inid=NA
for (i in tmp$id_participant){
  tmp$group_inid[which(tmp$id_participant==i)]= rep(1:length(tmp$values[which(tmp$id_participant==i)][[1]]),as.numeric(tmp$lengths[which(tmp$id_participant==i)][[1]]))
  
}


tmp2=tmp %>% group_by(id_participant, group_inid) %>% 
  mutate(dmtstartdate=min(start), 
         dmtendated=max(end),
         dur=sum(dur))%>%
  distinct(id_participant, group_inid,.keep_all = T)%>%
  ungroup()%>%
  select(id_participant, Treatment_GenericName,dmtstartdate, dmtendated, Efficacy, dur)

write.csv(tmp2, "../DATA/2025/CleanDMT.csv", row.names = F)



