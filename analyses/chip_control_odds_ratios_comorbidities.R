#script to get list of cases and controls from wave 8 and 9 for calculation of odds ratios for common comorbidities
library(readxl)
library(dplyr)

#import wave 8/9 comorbidities
setwd("N:/My Documents/ELSA_data/upgrade")
w_8_9_comorbs <- read.csv("ELSA_wave_8_9_comorbidities_all.csv", header = TRUE)

#import wave 8/9 drugs
w8_9_drug_classes <- read.csv("wave_8_9_drug_classes.csv", header = TRUE)

#import list of cases and select only those in wave 8/9
setwd("N:/My Documents/ELSA_data")
chip_cases <- read_xlsx("Waves 2,4,6,8,9 Lab IDs and Idauniq -Chip cases.xlsx")
chip_cases_8_9 <- chip_cases %>% dplyr::filter(elsa_wave == 8 | elsa_wave == 9)

#import list of controls and select only those in wave 8/9
controls <- read_xlsx("Waves 2,4,6,8,9 Lab IDs and Idauniq -Control cases.xlsx")
controls_8_9 <- controls %>% dplyr::filter(elsa_wave == 8 | elsa_wave == 9)

#combine w_8_9_comorbs and drug classes
w_8_9_metadata <- merge(w_8_9_comorbs, w8_9_drug_classes, by = "idauniq", all.x = TRUE)

#wave 8 9 cases metadata
w_8_9_metadata_chip_cases <- w_8_9_metadata %>% dplyr::filter(idauniq %in% chip_cases_8_9$idauniq)

#wave 8 9 controls metadata
w_8_9_metadata_control_cases <- w_8_9_metadata %>% dplyr::filter(idauniq %in% controls_8_9$idauniq)

table(w_8_9_metadata_control_cases$idauniq %in% w_8_9_metadata_chip_cases$idauniq)

#overlap cases - from pilot with varying ocncentrations of dna and different platforms - exclude from cases and controls
overlaps <- w_8_9_metadata_control_cases %>% dplyr::filter(w_8_9_metadata_control_cases$idauniq %in% w_8_9_metadata_chip_cases$idauniq)

w_8_9_metadata_chip_cases_post_excl <- w_8_9_metadata_chip_cases %>% dplyr::filter(! idauniq %in% overlaps$idauniq)

w_8_9_metadata_control_cases_post_excl <- w_8_9_metadata_control_cases %>% dplyr::filter(! idauniq %in% overlaps$idauniq)

#add chip or control variable
w_8_9_metadata_chip_cases_post_excl$status <- "chip"
w_8_9_metadata_control_cases_post_excl$status <- "control"

#rbind together
late_wave_cases_and_controls <- rbind(w_8_9_metadata_chip_cases_post_excl, w_8_9_metadata_control_cases_post_excl)

#remove those with blood disorder
late_wave_cases_and_controls <- late_wave_cases_and_controls %>% dplyr::filter(blood_dis == 0)

#replace age -7 (indicating age > 90) with 90
late_wave_cases_and_controls$indager <- replace(late_wave_cases_and_controls$indager, which(late_wave_cases_and_controls$indager == -7), 90)


#####now look at odds ratios
library(epitools)
#match for age, sex and smoking status
library(MatchIt)

m_out1 <- matchit(status_num ~ indager + indsex + current_smoker, data = late_wave_cases_and_controls_dementia, method = "nearest", distance = "glm")
m_data <- match.data(m_out1)

#ihd 
table(m_data$status, m_data$ihd)
oddsratio(as.matrix(table(m_data$ihd, m_data$status)), rev = "c")

#stroke 
oddsratio(as.matrix(table(m_data$stroke, m_data$status)), rev = "c")

#htn 
table(m_data$status, m_data$htn)
oddsratio(as.matrix(table(m_data$htn, m_data$status)), rev = "c")

#asthma 
oddsratio(as.matrix(table(m_data$asthma, m_data$status)), rev = "c")

#lung_dis 
oddsratio(as.matrix(table(m_data$lung_disease, m_data$status)), rev = "c")

#inflam_dis 
oddsratio(as.matrix(table(m_data$inflam_dis, m_data$status)), rev = "c")

#other heart disease 
oddsratio(as.matrix(table(m_data$heart_dis, m_data$status)), rev = "c")


#osteoporosis 
oddsratio(as.matrix(table(m_data$osteoporosis, m_data$status)), rev = "c")

#diabetes
oddsratio(as.matrix(table(m_data$diabetes, m_data$status)), rev = "c")




