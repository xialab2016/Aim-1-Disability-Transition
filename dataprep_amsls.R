##################################
##  Prepared by Valery Fuh-Ngwa ##
##################################

# Set working directory
setwd("C:/Users/vfngwa/OneDrive - University of Tasmania/Documents/InForMS Project/InForMS Data")

# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyverse)
library(msm)
library(expss)
library(lubridate)
library(caret)
library(mstate)
library(zoo)

########################################
## Step 1: Read and Process Demographic Data ##
########################################

demo <- read_excel("inforMS_PortalDemographics.xlsx")
demo <- demo[,-4]
names(demo) <- c("id", "dob", "sex", "yearfde", "yeardiag", "mstype")

# Recode and clean demographic data
demo <- demo %>%
  mutate(
    sex = 2 - sex,  # 1 = Male, 0 = Female
    across(4:6, ~ replace(., . %in% -9999, NA)),
    datediag = as.Date(paste(yeardiag, "/07/01", sep = "")),
    datefde = as.Date(paste(yearfde, "/07/01", sep = "")),
    mstype = case_when(
      mstype == 1 ~ "ROMS",
      mstype == 2 ~ "POMS",
      mstype %in% c(6, 7) | is.na(mstype) ~ "CIS|UNKNOWN",
      TRUE ~ as.character(mstype)
    ),
    progstat = ifelse(mstype == "ROMS", 0, 1),
    datediag = ifelse(is.na(datediag) & !is.na(datefde), as.character(datefde), as.character(datediag)),
    datefde = ifelse(is.na(datefde) & !is.na(datediag), as.character(datediag), as.character(datefde)),
    dage = round(as.numeric(as.Date(datediag) - as.Date(dob)) / 365.25, 2),
    oage = round(as.numeric(as.Date(datefde) - as.Date(dob)) / 365.25, 2),
    tsymp = abs(round(as.numeric(as.Date(datediag) - as.Date(datefde)) / 365.25, 2))
  ) %>%
  select(-yearfde, -yeardiag)

######################################
## Step 2: Read and Process Disability Data ##
######################################

disab <- read_excel("inforMS_PortalDisabilityLevel.xlsx")
names(disab) <- c("id", "reviewdate", "pdds")

# Merge disability data with demographics and clean it
disab <- disab %>%
  left_join(demo, by = "id") %>%
  filter(!is.na(id)) %>%
  arrange(id, reviewdate) %>%
  group_by(id) %>%
  mutate(
    lagrev = lag(reviewdate),
    t_interval = as.numeric(as.Date(reviewdate) - as.Date(lagrev)) / 365.25,
    state = case_when(
      pdds == 0 ~ 1,
      pdds == 1 ~ 2,
      pdds == 2 ~ 3,
      pdds == 3 ~ 4,
      TRUE ~ 5
    )
  ) %>%
  #filter(!is.na(pdds) & (t_interval == 0 | t_interval >= 0.25)) %>%
  ungroup() %>%
  filter(id %in% id[duplicated(id)])

disab$reviewdate[disab$reviewdate=="9999-09-09"]=as.Date.character("1999-09-09")
# Add time-related features and filter for 20 years
disab <- disab %>%
  group_by(id) %>%
  mutate(
    visit = row_number(),
    t_interval=ifelse(visit==1, 0, t_interval),
    obstime = cumsum(log(t_interval + 1)),
    years = cumsum(t_interval)
  ) %>%
  distinct()


######################################
## Step 3: Read and Process DMT Data ##
######################################

# Read and clean DMT data
dmt <- read_excel("DMTHistoryFinal.xlsx") %>%
  rename(id = 1) %>%
  arrange(id, DMTStartDate) %>%
  group_by(id) %>%
  mutate(
    DMT = zoo::na.locf(DMT, na.rm = FALSE),
    DMT = zoo::na.locf(DMT, fromLast = TRUE, na.rm = FALSE)
  ) %>%
  select(id, DMTStartDate, DMTEndDate, DMT) %>%
  rename(dmtstartdate = DMTStartDate, dmtenddate = DMTEndDate, dmt = DMT) %>%
  mutate(across(c(dmtstartdate, dmtenddate), as.Date))%>%
  arrange(id, dmtstartdate, dmtenddate ) %>%
  ungroup()


##replace DMT names with common names
# Define a lookup table for common names
dmt_common_names <- c(
  "Betaferon" = "Interferon beta-1b",
  "Aubagio" = "Teriflunomide",
  "Gilenya" = "Fingolimod",
  "Tecfidera" = "Dimethyl Fumarate",
  "Avonex" = "Interferon beta-1a",
  "Rebif" = "Interferon beta-1a",
  "Copaxone" = "Glatiramer Acetate",
  "Mavenclad" = "Cladribine",
  "Ocrevus" = "Ocrelizumab",
  "Plegridy" = "Interferon beta-1a",
  "Tysabri" = "Natalizumab",
  "Clinical Trial" = "Other",
  "Lemtrada" = "Alemtuzumab",
  "Mayzent" = "Siponimod",
  "Azasan" = "Azathioprine",
  "Trexall" = "Methotrexate",
  "Mitoxantrone" = "Mitoxantrone",
  "Kesimpta" = "Ofatumumab",
  "UsedZep" = "Other",
  "Zinbryta" = "Daclizumab",
  "Tysabir" = "Natalizumab",
  "Stem Cell Transplant" = "Other"
)

# Harmonised DMT Processing Block (ddmt_high/low)

# Replace DMT names with common names
if (!"dmt_cname" %in% names(dmt)) {
  dmt <- dmt %>%
    mutate(dmt_cname = ifelse(dmt %in% names(dmt_common_names), dmt_common_names[dmt], "other"))
}

# Define high/low efficacy class
high_efficacy_dmts <- c("Natalizumab", "Ocrelizumab", "Alemtuzumab", "Cladribine", "Ofatumumab", "Mitoxantrone", "Rituximab", "Daclizumab", "Siponimod", "Autologous Hematopoietic Stem Cell Transplantation")
low_efficacy_dmts <- c("Interferon beta-1b", "Interferon beta-1a", "Glatiramer Acetate", "Fingolimod", "Teriflunomide", "Azathioprine", "Dimethyl Fumarate", "Methotrexate")

# Replace missing dmtenddate with dmtstartdate + 1 day
dmt_daily <- dmt %>%
  mutate(
    # Replace missing end dates with start date + 1 day
    dmtenddate = if_else(is.na(dmtenddate), dmtstartdate + 1, dmtenddate),
    
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
cum_exposure_start <- disab %>%
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
disab <- disab %>%
  left_join(cum_exposure_start, by = c("id", "lagrev")) %>%
  mutate(
    ddmt_high=replace_na(ddmt_high,0),
    ddmt_low=replace_na(ddmt_low,0),
    ddmt_high=log(ddmt_high+1),
    ddmt_low=log(ddmt_low+1)
  )

# Add binary treatment at start of interval
dmt_start <- dmt %>%
  mutate(
    dmtenddate = if_else(is.na(dmtenddate), dmtstartdate + 1, dmtenddate),
    dmt_class = case_when(
      dmt_cname %in% high_efficacy_dmts ~ "high",
      dmt_cname %in% low_efficacy_dmts ~ "low",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(dmt_class))

# Check if on treatment at lagrev
disab <- disab %>%
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


disab <- disab %>%
  group_by(id) %>%
  mutate(
    maxstate = max(state, na.rm = TRUE), 
    absorb_time = suppressWarnings(min(obstime[state == 5], na.rm = TRUE)),
  )

# Remove observations where maxstate = 5, state < 5, and obstime > absorb_time
disab <- disab %>%
  filter(!(maxstate == 5 & state < 5 & obstime > absorb_time))

# Arrange and remove consecutive duplicate states
disab <- disab %>%
  arrange(id, obstime) %>%
  group_by(id) %>%
  mutate(
    lagstate = lag(state),
    lagst = ifelse(is.na(lagstate), state, lagstate),
    trans = paste(lagst, "-", state, sep = "")
  ) %>%
  filter(state != lagstate | is.na(lagstate)) %>%
  distinct(id, obstime, .keep_all = TRUE)

disab = disab %>%
  group_by(id) %>%
  mutate(
    n = n(),
    on_high=ifelse(ddmt_high>0&ddmt_low>0, 1, on_high),
    on_low=ifelse(ddmt_low>0&ddmt_low>0, 1, on_low),
    visit = row_number(),
    logt=log(years+1),
    oage=log(oage),
    tsymp=log(tsymp+1),
    maxstate = max(state, na.rm = TRUE), 
    absorb_time = min(obstime[state == 5], na.rm = TRUE),
    lagstate = lag(state),
  ) %>%
  filter(n > 1&years<=25) %>%                                       # Keep only ids with multiple observations
  filter(!(maxstate == 5 & state < 5 & obstime > absorb_time)) %>%  # Remove invalid states after absorption
  filter(!(state == 5 & lagstate == 5)) %>%                         # Remove repeated state 5 transitions
  dplyr::select(-n, -maxstate, -absorb_time, -lagstate, -lagst)     # Remove temporary columns
disab <- disab[duplicated(disab$id) | duplicated(disab$id, fromLast = TRUE), ]
msdat<-data.frame(disab)
setwd("C:/Users/vfngwa/OneDrive - University of Tasmania/Documents/InForMS Project/Hu chen analysis of pdds")
save(msdat, file="msdat.rdata")
