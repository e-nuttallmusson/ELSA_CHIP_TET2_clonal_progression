#script to summarise comorbidities in ELSA per wave

library(dplyr)
library(tidyr)

#setwd
setwd("~/Documents/UCL/PhD/ELSA_data/clean_data")

#read in survey data
wave_2 <- read.csv("wave_2_clean.txt", header = TRUE, sep = "\t")
wave_4 <- read.csv("wave_4_clean.txt", header = TRUE, sep = "\t")
wave_6 <- read.csv("wave_6_clean.txt", header = TRUE, sep = "\t")
wave_8 <- read.csv("wave_8_clean_new.txt", header = TRUE, sep = "\t")
wave_9 <- read.csv("wave_9_clean.txt", header = TRUE, sep = "\t")

#wave 8 and 9 data
#check for degree of overlap
table(wave_8$idauniq %in% wave_9$idauniq)

#FALSE  TRUE 
#1299  7146 

#therefore 1299 individuals who aren't in wave 9

wave_8$w8_only <- ifelse(wave_8$idauniq %in% wave_9$idauniq, FALSE, TRUE)

wave_8_not_9 <- wave_8 %>% dplyr::filter(w8_only == TRUE)
wave_8_not_9 <- subset(wave_8_not_9, select = -c(w8_only))

#make wave 9 df with composite comorbidities coded
wave_9_comorbidities <- wave_9

#consider wave 9 data - 8736 individuals
#diagnosis of malignant blood disorder - hedbwbl

wave_9_blood <- wave_9 %>% dplyr::filter(hedbwbl == 10 | hedibbl == 1 | hedbdbl == 1 | hedbsbl == 1)
nrow(wave_9_blood)
table(wave_9_blood$hedbwbl)
table(wave_9_blood$hedbdbl)
table(wave_9_blood$hedbsbl)
table(wave_9_blood$hedibbl)

table(wave_9$hedbwbl)
table(wave_9$hedibbl)
table(duplicated(wave_9_blood))

wave_9_comorbidities$blood_dis <- ifelse(wave_9_comorbidities$hedbwbl == 10 | wave_9_comorbidities$hedibbl == 1 | wave_9_comorbidities$hedbdbl == 1 | wave_9_comorbidities$hedbsbl == 1, 1, 0)



#inflammatory conditions
#ra - heartra. other (non oa) arthritis - heartot
#ms/mnd - hedbsms
#no variable for ibd/psoriasis/autoimmune liver disease
#therefore just inflammatory arthritis
wave_9_arthritis <- wave_9 %>% dplyr::filter(heartra ==1 | heartot == 1)
nrow(wave_9_arthritis)
table(wave_9$heartra)
table(wave_9$heartoa, wave_9$heartot)

wave_9_comorbidities$inflam_dis <- ifelse(wave_9_comorbidities$heartra ==1 | wave_9_comorbidities$heartot == 1, 1, 0)
table(wave_9_comorbidities$inflam_dis)

#whether has cancer 
#whether still has cancer - hedbsca
#new diagnosis of cancer - hedibca
#received cancer treatment in last 2 years - hecanb
wave_9_cancer <- wave_9 %>% dplyr::filter(hedbsca == 1 | hedibca ==1 | hecanb ==1)
nrow(wave_9_cancer)

wave_9_comorbidities$cancer <- ifelse(wave_9_comorbidities$hedbsca == 1 | wave_9_comorbidities$hedibca ==1 | wave_9_comorbidities$hecanb ==1, 1 , 0)
table(wave_9_comorbidities$cancer)

#simplify variables for key main comorbidities: diabetes, heart disease, stroke, osteoporosis, asthma, copd, smoking, hypertension, alcohol consumption
#diabetes:
#confirms diagnosis of diabetes/high bm: hedacdi
#told has diabetes by doctor: heacd
#on insulin: heins
#on medication for diabetes: hemdb
#hediadi == 1, hedimdi == 1 - these refer to bp and diabetes and shouldn't be included
wave_9_diabetes <- wave_9 %>% dplyr::filter(hedacdi == 1 | heacd == 1 | heins == 1 | hemdb == 1)
nrow(wave_9_diabetes)

wave_9_comorbidities$diabetes <- ifelse(wave_9_comorbidities$hedacdi == 1 | wave_9_comorbidities$heacd == 1 | wave_9_comorbidities$heins == 1 | wave_9_comorbidities$hemdb == 1, 1, 0)

#ischaemic heart disease
#diagnosis of heart attack fed forward: hedawmi ==3 or 2
#new diagnosis of heart attack: hediami == 1, hedimmi = 1
#angina fed forward: hedawan == 2 or 3
#new diagnosis of angina: hediman == 1, hediaan == 1

wave_9_ihd <- wave_9 %>% dplyr::filter(hedawmi == 3 | hedawmi == 2 | hedimmi == 1 | hediami == 1 | hedawan == 2 | hedawan == 3 | hediman == 1 | hediaan == 1)
nrow(wave_9_ihd)
wave_9_comorbidities$ihd <- ifelse(wave_9_comorbidities$hedawmi == 3 | wave_9_comorbidities$hedawmi == 2 | wave_9_comorbidities$hedimmi == 1 | wave_9_comorbidities$hediami == 1 | wave_9_comorbidities$hedawan == 2 | wave_9_comorbidities$hedawan == 3 | wave_9_comorbidities$hediman == 1 | wave_9_comorbidities$hediaan == 1, 1, 0)

#heart disease other: arrhythmias, heart failure, valvular heart disease
#confirms chf: hedachf == 1
#confirms murmur: hedachm == 1
#confirms arrhythmia: hedacar == 1
#confirms other heart disease: hedac95 == 1
#new diagnosis of heart failure: hediahf == 1, hedimhf ==1
#new diagnosis of heart murmur: hediahm == 1, hedimhm == 1
#new diagnosis of arrhythmia: hediaar == 1, hedimar ==1
#hedim85 == 1
wave_9_heart_disease_other <- wave_9 %>% dplyr::filter(hedachf == 1 | hedachm == 1 | hedacar == 1 | hedac95 == 1 | hediahf == 1 | hedimhf == 1 | hediahm == 1 | hedimhm == 1 | hediaar == 1 | hedimar == 1)
nrow(wave_9_heart_disease_other)

wave_9_comorbidities$heart_dis <- ifelse(wave_9_comorbidities$hedachf == 1 | wave_9_comorbidities$hedachm == 1 | wave_9_comorbidities$hedacar == 1 | wave_9_comorbidities$hedac95 == 1 | wave_9_comorbidities$hediahf == 1 | wave_9_comorbidities$hedimhf == 1 | wave_9_comorbidities$hediahm == 1 | wave_9_comorbidities$hedimhm == 1 | wave_9_comorbidities$hediaar == 1 | wave_9_comorbidities$hedimar == 1, 1, 0)

#stroke
#hedacst == 1 | hediast == 1 | hedimst == 1
wave_9_stroke <- wave_9 %>% dplyr::filter(hedacst == 1 | hediast == 1 | hedimst == 1)
nrow(wave_9_stroke) #304 individuals

wave_9_comorbidities$stroke <- ifelse(wave_9_comorbidities$hedacst == 1 | wave_9_comorbidities$hediast == 1 | wave_9_comorbidities$hedimst == 1, 1, 0)

#osteoporosis
#hedbdos == 1 | hedibos == 1

wave_9_osteoporosis <- wave_9 %>% dplyr::filter(hedbdos == 1| hedibos == 1)
nrow(wave_9_osteoporosis) #712 individuals

wave_9_comorbidities$osteoporosis <- ifelse(wave_9_comorbidities$hedbdos == 1| wave_9_comorbidities$hedibos == 1, 1, 0)

#asthma
#hedbdas == 1 | hedibas == 1 | heama ==1 | heamb == 1
wave_9_asthma <- wave_9 %>% dplyr::filter(hedbdas == 1 | hedibas == 1 | heama == 1 | heamb == 1)
nrow(wave_9_asthma) #983 individuals

wave_9_comorbidities$asthma <- ifelse(wave_9_comorbidities$hedbdas == 1 | wave_9_comorbidities$hedibas == 1 | wave_9_comorbidities$heama == 1 | wave_9_comorbidities$heamb == 1, 1, 0)

#lung disease
#hedbdlu == 1 | hedblu == 1 | hediblu == 1
wave_9_lung_disease <- wave_9 %>% dplyr::filter(hedbdlu ==1 | hedblu ==1 | hediblu ==1)
nrow(wave_9_lung_disease) #458 individuals

wave_9_comorbidities$lung_disease <- ifelse(wave_9_comorbidities$hedbdlu ==1 | wave_9_comorbidities$hedblu ==1 | wave_9_comorbidities$hediblu ==1, 1, 0)

#smoking
#current smoker (1): heska == 1
#ex smoker (2): heskf == 1
wave_9_current_smoker <- wave_9 %>% dplyr::filter(heska == 1 | heskf == 1)
nrow(wave_9_current_smoker) #815 individuals

wave_9_comorbidities$current_smoker <- ifelse(wave_9_comorbidities$heska == 1 | wave_9_comorbidities$heskf == 1, 1, 0)

#hypertension
# hedacbp == 1 | hediabp == 1 | hedimbp == 1 | hehbpb == 1 | hemda == 1 | hemdab == 1
wave_9_htn <- wave_9 %>% dplyr::filter(hedacbp == 1 | hediabp == 1 | hedimbp == 1 | hehbpb == 1| hemdab == 1)
nrow(wave_9_htn) #3241 individuals

wave_9_comorbidities$htn <- ifelse(wave_9_comorbidities$hedacbp == 1 | wave_9_comorbidities$hediabp == 1 | wave_9_comorbidities$hedimbp == 1 | wave_9_comorbidities$hehbpb == 1| wave_9_comorbidities$hemdab == 1, 1, 0)

#alcohol consumption
#heavy = 1, moderate = 2, occasional = 3, none = 4
#heavy = scalcm == 1 | scalcm == 2
#mod: scalcm == 3 | scalcm == 4 
#occ: scalcm == 5 | scalcm == 6 | scalcm == 7
#none: scalcm == 8 | scalcm == -1
#only missing values are those who didn't answer - 14 individuals
wave_9_heavy_alcohol <- wave_9 %>% dplyr::filter(scalcm == 1 | scalcm == 2)
nrow(wave_9_heavy_alcohol) #1381 individuals


wave_9_mod_alcohol <- wave_9 %>% dplyr::filter(scalcm == 3 | scalcm == 4)
nrow(wave_9_mod_alcohol) #2871

wave_9_occ_alcohol <- wave_9 %>% dplyr::filter(scalcm == 5 | scalcm == 6 | scalcm == 7)
nrow(wave_9_occ_alcohol)  #2136                                   
                                     
wave_9_no_alcohol <- wave_9 %>% dplyr::filter(scalcm == 8 | scalcm == -1)
nrow(wave_9_no_alcohol) #2277

wave_9_comorbidities$alcohol <- ifelse(wave_9_comorbidities$scalcm == 1 | wave_9_comorbidities$scalcm == 2, "heavy", NA)
wave_9_comorbidities$alcohol <- ifelse(wave_9_comorbidities$scalcm == 3 | wave_9_comorbidities$scalcm == 4, "moderate", wave_9_comorbidities$alcohol)
wave_9_comorbidities$alcohol <- ifelse(wave_9_comorbidities$scalcm == 5 | wave_9_comorbidities$scalcm == 6 | wave_9_comorbidities$scalcm == 7, "occasional", wave_9_comorbidities$alcohol)
wave_9_comorbidities$alcohol <- ifelse(wave_9_comorbidities$scalcm == 8 | wave_9_comorbidities$scalcm == -1, "none", wave_9_comorbidities$alcohol)

table(wave_9_comorbidities$alcohol)


###################################wave 8 data
#diagnosis of malignant blood disorder - hedbwbl
wave_8_not_9_comorbs <- wave_8_not_9

wave_8_not_9_comorbs$blood_dis <- ifelse(wave_8_not_9_comorbs$hedbwbl == 10 | wave_8_not_9_comorbs$hedibbl == 1 | wave_8_not_9_comorbs$hedbdbl == 1 | wave_8_not_9_comorbs$hedbsbl == 1, 1, 0)

wave_8_blood <- wave_8 %>% dplyr::filter(grepl(10, hedbwbl) | grepl(1, hedibbl))
nrow(wave_8_blood)
table(wave_9_blood$idauniq %in% wave_8_blood$idauniq)


#inflammatory conditions
#ra - heartra. other (non oa) arthritis - heartot
#ms/mnd - hedbsms
#no variable for ibd/psoriasis/autoimmune liver disease
wave_8_not_9_comorbs$inflam_dis <- ifelse(wave_8_not_9_comorbs$heartra ==1 | wave_8_not_9_comorbs$heartot == 1, 1, 0)
table(wave_8_not_9_comorbs$inflam_dis)

#whether has cancer 
#whether still has cancer - hedbsca
#new diagnosis of cancer - hedibca
#received cancer treatment in last 2 years - hecanb
wave_8_not_9_comorbs$cancer <- ifelse(wave_8_not_9_comorbs$hedbsca == 1 | wave_8_not_9_comorbs$hedibca ==1 | wave_8_not_9_comorbs$hecanb ==1, 1 , 0)
table(wave_8_not_9_comorbs$cancer)


#simplify variables for key main comorbidities: diabetes, heart disease, stroke, osteoporosis, asthma, copd, smoking, hypertension, alcohol consumption
#diabetes:
#confirms diagnosis of diabetes/high bm: hedacdi
#told has diabetes by doctor: heacd
#on insulin: heins
#on medication for diabetes: hemdb
wave_8_not_9_comorbs$diabetes <- ifelse(wave_8_not_9_comorbs$hedacdi == 1 | wave_8_not_9_comorbs$heacd == 1 | wave_8_not_9_comorbs$heins == 1 | wave_8_not_9_comorbs$hemdb == 1, 1, 0)


#ischaemic heart disease
#diagnosis of heart attack fed forward: hedawmi ==3
#new diagnosis of heart attack: hediami == 1, hedimmi = 1
#angina fed forward: hedawan == 2
#new diagnosis of angina: hediman == 1, hediaan == 1
wave_8_not_9_comorbs$ihd <- ifelse(wave_8_not_9_comorbs$hedawmi == 3 | wave_8_not_9_comorbs$hedawmi == 2 | wave_8_not_9_comorbs$hedimmi == 1 | wave_8_not_9_comorbs$hediami == 1 | wave_8_not_9_comorbs$hedawan == 2 | wave_8_not_9_comorbs$hedawan == 3 | wave_8_not_9_comorbs$hediman == 1 | wave_8_not_9_comorbs$hediaan == 1, 1, 0)

#heart disease other: arrhythmias, heart failure, valvular heart disease
#confirms chf: hedachf == 1
#confirms murmur: hedachm == 1
#confirms arrhythmia: hedacar == 1
#confirms other heart disease: hedac95 == 1
#new diagnosis of heart failure: hediahf == 1, hedimhf ==1
#new diagnosis of heart murmur: hediahm == 1, hedimhm == 1
#new diagnosis of arrhythmia: hediaar == 1, hedimar ==1
#hedim85 == 1
wave_8_not_9_comorbs$heart_dis <- ifelse(wave_8_not_9_comorbs$hedachf == 1 | wave_8_not_9_comorbs$hedachm == 1 | wave_8_not_9_comorbs$hedacar == 1 | wave_8_not_9_comorbs$hedac95 == 1 | wave_8_not_9_comorbs$hediahf == 1 | wave_8_not_9_comorbs$hedimhf == 1 | wave_8_not_9_comorbs$hediahm == 1 | wave_8_not_9_comorbs$hedimhm == 1 | wave_8_not_9_comorbs$hediaar == 1 | wave_8_not_9_comorbs$hedimar == 1 | wave_8_not_9_comorbs$hedim85 == 1, 1, 0)

#stroke
#hedacst == 1 | hediast == 1 | hedimst == 1

wave_8_not_9_comorbs$stroke <- ifelse(wave_8_not_9_comorbs$hedacst == 1 | wave_8_not_9_comorbs$hediast == 1 | wave_8_not_9_comorbs$hedimst == 1, 1, 0)

#osteoporosis
#hedbdos == 1 | hedibos == 1

wave_8_not_9_comorbs$osteoporosis <- ifelse(wave_8_not_9_comorbs$hedbdos == 1| wave_8_not_9_comorbs$hedibos == 1, 1, 0)


#asthma
#hedbdas == 1 | hedibas == 1 | heama ==1 | heamb == 1
wave_8_not_9_comorbs$asthma <- ifelse(wave_8_not_9_comorbs$hedbdas == 1 | wave_8_not_9_comorbs$hedibas == 1 | wave_8_not_9_comorbs$heama == 1 | wave_8_not_9_comorbs$heamb == 1, 1, 0)


#lung_disease
#hedbdlu == 1 | hedblu == 1 | hediblu == 1
wave_8_not_9_comorbs$lung_disease <- ifelse(wave_8_not_9_comorbs$hedbdlu ==1 | wave_8_not_9_comorbs$hedblu ==1 | wave_8_not_9_comorbs$hediblu ==1, 1, 0)


#smoking
#current smoker (1): heska == 1
#ex smoker (2): heskf == 1
wave_8_not_9_comorbs$current_smoker <- ifelse(wave_8_not_9_comorbs$heska == 1 | wave_8_not_9_comorbs$heskf == 1, 1, 0)

#hypertension
# hedacbp == 1 | hediabp == 1 | hedimbp == 1 | hehbpb == 1 | hemda == 1 | hemdab == 1
wave_8_not_9_comorbs$htn <- ifelse(wave_8_not_9_comorbs$hedacbp == 1 | wave_8_not_9_comorbs$hediabp == 1 | wave_8_not_9_comorbs$hedimbp == 1 | wave_8_not_9_comorbs$hehbpb == 1| wave_8_not_9_comorbs$hemdab == 1 , 1, 0)


#alcohol consumption
#heavy = 1, moderate = 2, occasional = 3, none = 4
#heavy = scalcm == 1 | scalcm == 2
#mod: scalcm == 3 | scalcm == 4 
#occ: scalcm == 5 | scalcm == 6 | scalcm == 7
#none: scalcm == 8 | scalcm == -1
#only missing values are those who didn't answer - 14 individuals
wave_8_not_9_comorbs$alcohol <- ifelse(wave_8_not_9_comorbs$scako == 1 | wave_8_not_9_comorbs$scako == 2, "heavy", NA)
wave_8_not_9_comorbs$alcohol <- ifelse(wave_8_not_9_comorbs$scako == 3 | wave_8_not_9_comorbs$scako == 4, "moderate", wave_8_not_9_comorbs$alcohol)
wave_8_not_9_comorbs$alcohol <- ifelse(wave_8_not_9_comorbs$scako == 5 | wave_8_not_9_comorbs$scako == 6 | wave_8_not_9_comorbs$scako == 7, "occasional", wave_8_not_9_comorbs$alcohol)
wave_8_not_9_comorbs$alcohol <- ifelse(wave_8_not_9_comorbs$scako == 8 | wave_8_not_9_comorbs$scako == -1, "none", wave_8_not_9_comorbs$alcohol)



#now merge key comorbidities into late wave comorb df

wave_9_comorbidities2 <- subset(wave_9_comorbidities, select = c(idauniq, indager, indsex, blood_dis, inflam_dis, cancer, diabetes, ihd, heart_dis, stroke, osteoporosis, asthma, lung_disease, current_smoker, htn, alcohol))
wave_8_not_9_comorbidities2 <- subset(wave_8_not_9_comorbs, select = c(idauniq, indager, indsex, blood_dis, inflam_dis, cancer, diabetes, ihd, heart_dis, stroke, osteoporosis, asthma, lung_disease, current_smoker, htn, alcohol))

#bind together

late_wave_comorbidities <- rbind(wave_9_comorbidities2, wave_8_not_9_comorbidities2)


setwd("~/Documents/UCL/PhD/Supervision/upgrade/upgrade_data_files/ELSA")

write.csv(late_wave_comorbidities, "ELSA_wave_8_9_comorbidities_all.csv", row.names = FALSE)
