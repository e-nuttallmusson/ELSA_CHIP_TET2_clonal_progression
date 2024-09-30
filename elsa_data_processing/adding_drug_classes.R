#script to clean and categorise elsa drug data from nurse visit wave 6

setwd("~/Documents/UCL/PhD/ELSA_data/drug_information")

#import BNF codes
bnf_data <- read.csv("20230901_1693566430911_BNF_Code_Information.csv")

#import elsa wave 6 nurse data

setwd("~/Documents/UCL/PhD/ELSA_data/wave_data")

wave_6_nurse <- read.table("wave_6_elsa_nurse_data_v2.tab", header = TRUE, sep = "\t")

vars <- c("idauniq")

for (i in (1:27)) {
  a <- paste0("DrC",i)
  vars <- append(vars, a)
  rm(a)
}

#simplify to just drugs
wave_6_drugs <- subset(wave_6_nurse, select = vars)

#creat key value pairs of bnf drug codes
drug_codes <- append(bnf_data$BNF.Paragraph.Code, "-1")

drug_classes <- append(bnf_data$BNF.Paragraph, NA)

#some common drugs missing from classification - manually reviewed and coded

# 20551: Angiotensin-converting enzyme (ACE) inhibitors
# 20552: Angiotensin II receptor antagonists
# 21201: Statins
# 21202: Other lipid-lowering drugs
# 30100: Corticosteroids - asthma
# 50200: Antifungal
# 60105: Analgesia for diabetic neuropathy
# 60121: Sulphonylureas
# 60122: Biguanides (e.g. Metformin)
# 60123: Other antidiabetic medications
# 999996: Unable to code
# 50300: Antivirals
# 40803: Anticonvulsants (febrile)
# 131300: Hirudoid
# 20553: Renin inhibitors


drug_codes <- append(drug_codes, c("20551", "20552",  "21201",  "21202",  "30100",  "50200",  "60105",  "60121",  "60122", "60123", "999996"))
drug_codes <- append(drug_codes, "50300")
drug_codes <- append(drug_codes, "131300")
drug_codes <- append(drug_codes, "40803")
drug_codes <- append(drug_codes, "20553")

drug_classes <- append(drug_classes, c("Angiotensin-converting enzyme (ACE) inhibitors", "Angiotensin II receptor antagonists", "Statins", "Other lipid-lowering drugs", "Corticosteroids - asthma", "Antifungal", "Analgesia for diabetic neuropathy", "Sulphonylureas", "Biguanides (e.g. Metformin)", "Other antidiabetic medications", "Unable to code"))                    
drug_classes <- append(drug_classes, "Antivirals")
drug_classes <- append(drug_classes, "Hirudoid")
drug_classes <- append(drug_classes, "Anticonvulsants (febrile)")
drug_classes <- append(drug_classes, "Renin inhibitors")

drug_list <- as.data.frame(cbind(drug_codes, drug_classes))

drug_list <- distinct(drug_list)

#for each drug column, add column of drug class
#drug 1
drug_1 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC1[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_1 <- append(drug_1, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs, drug_1)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_1, .after = DrC1)

#drug 2

drug_2 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC2[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_2 <- append(drug_2, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_2)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_2, .after = DrC2)

#drug 3

drug_3 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC3[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_3 <- append(drug_3, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_3)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_3, .after = DrC3)

#drug 4

drug_4 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC4[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_4 <- append(drug_4, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_4)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_4, .after = DrC4)

#60112 missing - ?antidiabetic but doesn't seem to correspond to identifiable drug

#drug 5

drug_5 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC5[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_5 <- append(drug_5, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_5)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_5, .after = DrC5)

#drug 6
drug_6 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC6[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_6 <- append(drug_6, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_6)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_6, .after = DrC6)

#40742 missing - not clear what this is

#drug 7
drug_7 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC7[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_7 <- append(drug_7, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_7)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_7, .after = DrC7)


#drug 8
drug_8 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC8[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_8 <- append(drug_8, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_8)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_8, .after = DrC8)

#drug 9
drug_9 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC9[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_9 <- append(drug_9, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_9)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_9, .after = DrC9)

#drug 10
drug_10 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC10[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_10 <- append(drug_10, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_10)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_10, .after = DrC10)

#drug 11
drug_11 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC11[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_11 <- append(drug_11, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_11)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_11, .after = DrC11)

#drug 12
drug_12 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC12[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_12 <- append(drug_12, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_12)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_12, .after = DrC12)

#drug 13
drug_13 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC13[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_13 <- append(drug_13, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_13)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_13, .after = DrC13)

#drug 14
drug_14 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC14[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_14 <- append(drug_14, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_14)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_14, .after = DrC14)

#30105 missing - probably some form of bronchodilator but no seciton in bnf

#drug 15
drug_15 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC15[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_15 <- append(drug_15, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_15)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_15, .after = DrC15)

#drug 16
drug_16 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC16[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_16 <- append(drug_16, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_16)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_16, .after = DrC16)

#drug 17
drug_17 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC17[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_17 <- append(drug_17, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_17)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_17, .after = DrC17)

#drug 18
drug_18 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC18[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_18 <- append(drug_18, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_18)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_18, .after = DrC18)

#drug 19
drug_19 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC19[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_19 <- append(drug_19, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_19)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_19, .after = DrC19)

#drug 20
drug_20 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC20[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_20 <- append(drug_20, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_20)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_20, .after = DrC20)

#drug 21
drug_21 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC21[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_21 <- append(drug_21, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_21)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_21, .after = DrC21)

#drug 22
drug_22 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC22[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_22 <- append(drug_22, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_22)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_22, .after = DrC22)

# 150200 missing - no bnf code

#drug 23
drug_23 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC23[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_23 <- append(drug_23, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_23)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_23, .after = DrC23)


#drug 24
drug_24 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC24[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_24 <- append(drug_24, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_24)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_24, .after = DrC24)

#drug 25
drug_25 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC25[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_25 <- append(drug_25, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_25)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_25, .after = DrC25)

#drug 26
drug_26 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC26[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_26 <- append(drug_26, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_26)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_26, .after = DrC26)

#drug 27
drug_27 <- c()

for (i in (1:nrow(wave_6_drugs))) {
  a <- wave_6_drugs$DrC27[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_27 <- append(drug_27, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_6_drugs_1 <- cbind(wave_6_drugs_1, drug_27)

wave_6_drugs_1 <- wave_6_drugs_1 %>% relocate(drug_27, .after = DrC27)


#add in data about whether taken in the last week
vars <- c("idauniq", "MedBIA")

for (i in (2:27)) {
  a <- paste0("MedBIA",i)
  vars <- append(vars, a)
  rm(a)
}

#create df of just whether taken in the last week

w6_week <- subset(wave_6_nurse, select = vars)

#combine with wave_6_drugs_1

wave_6_drugs_take <- merge(wave_6_drugs_1, w6_week, by="idauniq")

#rename MedBIA to MedBIA1

names(wave_6_drugs_take)[names(wave_6_drugs_take) == 'MedBIA'] <- 'MedBIA1'

#reorder columns so can see if taking

for (i in (1:27)) {
  wave_6_drugs_take <- wave_6_drugs_take %>% relocate(paste0("MedBIA",i), .after = paste0("drug_",i))
}

#check which classes of medications represented

medication_classes <- table(wave_6_drugs_take$drug_1)

#then need to simplify into key categories as to whether taking or not
#nsaids - Non-steroidal anti-inflammatory drugs, Other anti-inflammatory preparations
#antiplatelets - Antiplatelet drugs
#immunosuppressants - Corticosteroids, Corticosteroids - asthma, Corticosteroids and other immunosuppressants, Drugs affecting the immune response, Glucocorticoid therapy, Other drugs for rheumatic diseases, Other immunomodulating drugs, Rheumatic disease suppressant drugs, Aminosalicylates, Antiproliferative immunosuppressants
#metformin - Biguanides (e.g. Metformin)

nsaids <- c(rep(FALSE, nrow(wave_6_drugs_take)))

wave_6_drugs_take_class <- cbind(wave_6_drugs_take, nsaids)

wave_6_drugs_take_class <- wave_6_drugs_take_class %>% 
  mutate(nsaids = case_when(
    drug_1 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% c("Non-steroidal anti-inflammatory drugs", "Other anti-inflammatory preparations") & 1 %in% MedBIA27 ~ TRUE,
    TRUE ~ as.logical(nsaids)))


#antiplatelets

antiplatelets <- c(rep(FALSE, nrow(wave_6_drugs_take)))

wave_6_drugs_take_class <- cbind(wave_6_drugs_take_class, antiplatelets)

wave_6_drugs_take_class <- wave_6_drugs_take_class %>% 
  mutate(antiplatelets = case_when(
    drug_1 %in% "Antiplatelet drugs" & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% "Antiplatelet drugs" & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% "Antiplatelet drugs" & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% "Antiplatelet drugs" & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% "Antiplatelet drugs" & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% "Antiplatelet drugs" & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% "Antiplatelet drugs" & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% "Antiplatelet drugs" & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% "Antiplatelet drugs" & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% "Antiplatelet drugs" & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% "Antiplatelet drugs" & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% "Antiplatelet drugs" & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% "Antiplatelet drugs" & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% "Antiplatelet drugs" & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% "Antiplatelet drugs" & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% "Antiplatelet drugs" & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% "Antiplatelet drugs" & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% "Antiplatelet drugs" & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% "Antiplatelet drugs" & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% "Antiplatelet drugs" & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% "Antiplatelet drugs" & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% "Antiplatelet drugs" & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% "Antiplatelet drugs" & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% "Antiplatelet drugs" & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% "Antiplatelet drugs" & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% "Antiplatelet drugs" & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% "Antiplatelet drugs" & 1 %in% MedBIA27 ~ TRUE,
    TRUE ~ as.logical(antiplatelets)))


#immunosuppressants

immunosuppressants <- c(rep(FALSE, nrow(wave_6_drugs_take)))

wave_6_drugs_take_class <- cbind(wave_6_drugs_take_class, immunosuppressants)

wave_6_drugs_take_class <- wave_6_drugs_take_class %>% 
  mutate(immunosuppressants = case_when(
    drug_1 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% c("Corticosteroids", "Corticosteroids - asthma", "Corticosteroids and other immunosuppressants", "Drugs affecting the immune response", "Glucocorticoid therapy", "Other drugs for rheumatic diseases", "Other immunomodulating drugs", "Rheumatic disease suppressant drugs", "Aminosalicylates", "Antiproliferative immunosuppressants") & 1 %in% MedBIA27 ~ TRUE,
    TRUE ~ as.logical(immunosuppressants)))



#metformin

metformin <- c(rep(FALSE, nrow(wave_6_drugs_take)))

wave_6_drugs_take_class <- cbind(wave_6_drugs_take_class, metformin)

wave_6_drugs_take_class <- wave_6_drugs_take_class %>% 
  mutate(metformin = case_when(
    drug_1 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% "Biguanides (e.g. Metformin)" & 1 %in% MedBIA27 ~ TRUE,
    TRUE ~ as.logical(metformin)))


#RAAS system

RAAS <- c(rep(FALSE, nrow(wave_6_drugs_take)))

wave_6_drugs_take_class <- cbind(wave_6_drugs_take_class, RAAS)


wave_6_drugs_take_class <- wave_6_drugs_take_class %>% 
  mutate(RAAS = case_when(
    drug_1 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% MedBIA27 ~ TRUE,
    drug_1 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% "Angiotensin II receptor antagonists" & 1 %in% MedBIA27 ~ TRUE,
    TRUE ~ as.logical(RAAS)))

#calcium channel blockers

ccb <- c(rep(FALSE, nrow(wave_6_drugs_take)))

wave_6_drugs_take_class <- cbind(wave_6_drugs_take_class, ccb)


wave_6_drugs_take_class <- wave_6_drugs_take_class %>% 
  mutate(ccb = case_when(
    drug_1 %in% "Calcium-channel blockers" & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% "Calcium-channel blockers" & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% "Calcium-channel blockers" & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% "Calcium-channel blockers" & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% "Calcium-channel blockers" & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% "Calcium-channel blockers" & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% "Calcium-channel blockers" & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% "Calcium-channel blockers" & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% "Calcium-channel blockers" & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% "Calcium-channel blockers" & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% "Calcium-channel blockers" & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% "Calcium-channel blockers" & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% "Calcium-channel blockers" & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% "Calcium-channel blockers" & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% "Calcium-channel blockers" & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% "Calcium-channel blockers" & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% "Calcium-channel blockers" & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% "Calcium-channel blockers" & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% "Calcium-channel blockers" & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% "Calcium-channel blockers" & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% "Calcium-channel blockers" & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% "Calcium-channel blockers" & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% "Calcium-channel blockers" & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% "Calcium-channel blockers" & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% "Calcium-channel blockers" & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% "Calcium-channel blockers" & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% "Calcium-channel blockers" & 1 %in% MedBIA27 ~ TRUE,
    TRUE ~ as.logical(ccb)))

#beta blockers

bb <- c(rep(FALSE, nrow(wave_6_drugs_take)))

wave_6_drugs_take_class <- cbind(wave_6_drugs_take_class, bb)


wave_6_drugs_take_class <- wave_6_drugs_take_class %>% 
  mutate(bb = case_when(
    drug_1 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% MedBIA27 ~ TRUE,
    TRUE ~ as.logical(bb)))


#thiazides

thiazides <- c(rep(FALSE, nrow(wave_6_drugs_take)))

wave_6_drugs_take_class <- cbind(wave_6_drugs_take_class, thiazides)


wave_6_drugs_take_class <- wave_6_drugs_take_class %>% 
  mutate(thiazides = case_when(
    drug_1 %in% "Thiazides and related diuretics" & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% "Thiazides and related diuretics" & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% "Thiazides and related diuretics" & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% "Thiazides and related diuretics" & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% "Thiazides and related diuretics" & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% "Thiazides and related diuretics" & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% "Thiazides and related diuretics" & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% "Thiazides and related diuretics" & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% "Thiazides and related diuretics" & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% "Thiazides and related diuretics" & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% "Thiazides and related diuretics" & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% "Thiazides and related diuretics" & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% "Thiazides and related diuretics" & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% "Thiazides and related diuretics" & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% "Thiazides and related diuretics" & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% "Thiazides and related diuretics" & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% "Thiazides and related diuretics" & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% "Thiazides and related diuretics" & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% "Thiazides and related diuretics" & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% "Thiazides and related diuretics" & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% "Thiazides and related diuretics" & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% "Thiazides and related diuretics" & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% "Thiazides and related diuretics" & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% "Thiazides and related diuretics" & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% "Thiazides and related diuretics" & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% "Thiazides and related diuretics" & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% "Thiazides and related diuretics" & 1 %in% MedBIA27 ~ TRUE,
    TRUE ~ as.logical(thiazides)))

#sulphonylureas

sulphonylureas <- c(rep(FALSE, nrow(wave_6_drugs_take)))

wave_6_drugs_take_class <- cbind(wave_6_drugs_take_class, sulphonylureas)


wave_6_drugs_take_class <- wave_6_drugs_take_class %>% 
  mutate(sulphonylureas = case_when(
    drug_1 %in% "Sulphonylureas" & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% "Sulphonylureas" & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% "Sulphonylureas" & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% "Sulphonylureas" & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% "Sulphonylureas" & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% "Sulphonylureas" & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% "Sulphonylureas" & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% "Sulphonylureas" & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% "Sulphonylureas" & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% "Sulphonylureas" & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% "Sulphonylureas" & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% "Sulphonylureas" & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% "Sulphonylureas" & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% "Sulphonylureas" & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% "Sulphonylureas" & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% "Sulphonylureas" & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% "Sulphonylureas" & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% "Sulphonylureas" & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% "Sulphonylureas" & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% "Sulphonylureas" & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% "Sulphonylureas" & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% "Sulphonylureas" & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% "Sulphonylureas" & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% "Sulphonylureas" & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% "Sulphonylureas" & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% "Sulphonylureas" & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% "Sulphonylureas" & 1 %in% MedBIA27 ~ TRUE,
    TRUE ~ as.logical(sulphonylureas)))

#statins and other lipid lowering drugs

lld <- c(rep(FALSE, nrow(wave_6_drugs_take)))

wave_6_drugs_take_class <- cbind(wave_6_drugs_take_class, lld)


wave_6_drugs_take_class <- wave_6_drugs_take_class %>% 
  mutate(lld = case_when(
    drug_1 %in% "Statins" & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% "Statins" & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% "Statins" & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% "Statins" & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% "Statins" & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% "Statins" & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% "Statins" & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% "Statins" & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% "Statins" & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% "Statins" & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% "Statins" & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% "Statins" & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% "Statins" & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% "Statins" & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% "Statins" & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% "Statins" & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% "Statins" & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% "Statins" & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% "Statins" & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% "Statins" & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% "Statins" & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% "Statins" & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% "Statins" & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% "Statins" & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% "Statins" & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% "Statins" & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% "Statins" & 1 %in% MedBIA27 ~ TRUE,
    drug_1 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% "Other lipid-lowering drugs" & 1 %in% MedBIA27 ~ TRUE,
    TRUE ~ as.logical(lld)))


#ssris

ssri <- c(rep(FALSE, nrow(wave_6_drugs_take)))

wave_6_drugs_take_class <- cbind(wave_6_drugs_take_class, ssri)


wave_6_drugs_take_class <- wave_6_drugs_take_class %>% 
  mutate(ssri = case_when(
    drug_1 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% MedBIA27 ~ TRUE,
    TRUE ~ as.logical(ssri)))


#loop diuretics

loop <- c(rep(FALSE, nrow(wave_6_drugs_take)))

wave_6_drugs_take_class <- cbind(wave_6_drugs_take_class, loop)

wave_6_drugs_take_class <- wave_6_drugs_take_class %>% 
  mutate(loop = case_when(
    drug_1 %in% "Loop diuretics" & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% "Loop diuretics" & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% "Loop diuretics" & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% "Loop diuretics" & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% "Loop diuretics" & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% "Loop diuretics" & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% "Loop diuretics" & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% "Loop diuretics" & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% "Loop diuretics" & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% "Loop diuretics" & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% "Loop diuretics" & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% "Loop diuretics" & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% "Loop diuretics" & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% "Loop diuretics" & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% "Loop diuretics" & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% "Loop diuretics" & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% "Loop diuretics" & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% "Loop diuretics" & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% "Loop diuretics" & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% "Loop diuretics" & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% "Loop diuretics" & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% "Loop diuretics" & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% "Loop diuretics" & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% "Loop diuretics" & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% "Loop diuretics" & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% "Loop diuretics" & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% "Loop diuretics" & 1 %in% MedBIA27 ~ TRUE,
    TRUE ~ as.logical(loop)))


#ppis

ppi <- c(rep(FALSE, nrow(wave_6_drugs_take)))

wave_6_drugs_take_class <- cbind(wave_6_drugs_take_class, ppi)


wave_6_drugs_take_class <- wave_6_drugs_take_class %>% 
  mutate(ppi = case_when(
    drug_1 %in% "Proton pump inhibitors" & 1 %in% MedBIA1 ~ TRUE,
    drug_2 %in% "Proton pump inhibitors" & 1 %in% MedBIA2 ~ TRUE, 
    drug_3 %in% "Proton pump inhibitors" & 1 %in% MedBIA3 ~ TRUE,
    drug_4 %in% "Proton pump inhibitors" & 1 %in% MedBIA4 ~ TRUE,
    drug_5 %in% "Proton pump inhibitors" & 1 %in% MedBIA5 ~ TRUE,
    drug_6 %in% "Proton pump inhibitors" & 1 %in% MedBIA6 ~ TRUE,
    drug_7 %in% "Proton pump inhibitors" & 1 %in% MedBIA7 ~ TRUE,
    drug_8 %in% "Proton pump inhibitors" & 1 %in% MedBIA8 ~ TRUE,
    drug_9 %in% "Proton pump inhibitors" & 1 %in% MedBIA9 ~ TRUE,
    drug_10 %in% "Proton pump inhibitors" & 1 %in% MedBIA10 ~ TRUE,
    drug_11 %in% "Proton pump inhibitors" & 1 %in% MedBIA11 ~ TRUE,
    drug_12 %in% "Proton pump inhibitors" & 1 %in% MedBIA12 ~ TRUE,
    drug_13 %in% "Proton pump inhibitors" & 1 %in% MedBIA13 ~ TRUE,
    drug_14 %in% "Proton pump inhibitors" & 1 %in% MedBIA14 ~ TRUE,
    drug_15 %in% "Proton pump inhibitors" & 1 %in% MedBIA15 ~ TRUE,
    drug_16 %in% "Proton pump inhibitors" & 1 %in% MedBIA16 ~ TRUE,
    drug_17 %in% "Proton pump inhibitors" & 1 %in% MedBIA17 ~ TRUE,
    drug_18 %in% "Proton pump inhibitors" & 1 %in% MedBIA18 ~ TRUE,
    drug_19 %in% "Proton pump inhibitors" & 1 %in% MedBIA19 ~ TRUE,
    drug_20 %in% "Proton pump inhibitors" & 1 %in% MedBIA20 ~ TRUE,
    drug_21 %in% "Proton pump inhibitors" & 1 %in% MedBIA21 ~ TRUE,
    drug_22 %in% "Proton pump inhibitors" & 1 %in% MedBIA22 ~ TRUE,
    drug_23 %in% "Proton pump inhibitors" & 1 %in% MedBIA23 ~ TRUE,
    drug_24 %in% "Proton pump inhibitors" & 1 %in% MedBIA24 ~ TRUE,
    drug_25 %in% "Proton pump inhibitors" & 1 %in% MedBIA25 ~ TRUE,
    drug_26 %in% "Proton pump inhibitors" & 1 %in% MedBIA26 ~ TRUE,
    drug_27 %in% "Proton pump inhibitors" & 1 %in% MedBIA27 ~ TRUE,
    TRUE ~ as.logical(ppi)))




#drop unnecessary variables

wave_6_drugs_simp <- subset(wave_6_drugs_take_class, select = c(idauniq, nsaids, antiplatelets, immunosuppressants, metformin, RAAS, ccb, bb, thiazides, sulphonylureas, lld, ssri, loop, ppi))

#output file to combine with mutational data in DSH

write.table(wave_6_drugs_simp, "wave_6_drug_categories.xlsx", row.names = FALSE, sep = "\t")



##############################wave 8/9####################################

#import combined w8 and w9 nurse data

setwd("//idhs.ucl.ac.uk/user/User_Data/regmenu/My Documents/ELSA_data/w8_w9_modelling")

wave_8_9_nurse_data <- read.table("elsa_nurse_w8w9_data_sl_protect.tab", header = TRUE, sep = "\t")


#import BNF codes
bnf_data <- read.csv("20230901_1693566430911_BNF_Code_Information.csv")

#get just drug data
vars <- c("idauniq")

for (i in (1:22)) {
  a <- paste0("drc",i)
  vars <- append(vars, a)
  rm(a)
}

#simplify to just drugs
wave_8_9_drugs <- subset(wave_8_9_nurse_data, select = vars)

#creat key value pairs of bnf drug codes
drug_codes <- append(bnf_data$BNF.Paragraph.Code, "-1")

drug_classes <- append(bnf_data$BNF.Paragraph, NA)

drug_codes <- append(drug_codes, c("20551", "20552",  "21201",  "21202",  "30100",  "50200",  "60105",  "60121",  "60122", "60123", "999996"))
drug_codes <- append(drug_codes, "50300")
drug_codes <- append(drug_codes, "131300")
drug_codes <- append(drug_codes, "40803")
drug_codes <- append(drug_codes, "20553")

drug_classes <- append(drug_classes, c("Angiotensin-converting enzyme (ACE) inhibitors", "Angiotensin II receptor antagonists", "Statins", "Other lipid-lowering drugs", "Corticosteroids - asthma", "Antifungal", "Analgesia for diabetic neuropathy", "Sulphonylureas", "Biguanides (e.g. Metformin)", "Other antidiabetic medications", "Unable to code"))                    
drug_classes <- append(drug_classes, "Antivirals")
drug_classes <- append(drug_classes, "Hirudoid")
drug_classes <- append(drug_classes, "Anticonvulsants (febrile)")
drug_classes <- append(drug_classes, "Renin inhibitors")

drug_list <- as.data.frame(cbind(drug_codes, drug_classes))

drug_list <- distinct(drug_list)


#for each drug column, add column of drug class
#drug 1
drug_1 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc1[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_1 <- append(drug_1, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs, drug_1)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_1, .after = drc1)

#drug 2

drug_2 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc2[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_2 <- append(drug_2, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_2)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_2, .after = drc2)

#drug 3

drug_3 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc3[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_3 <- append(drug_3, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_3)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_3, .after = drc3)

#drug 4

drug_4 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc4[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_4 <- append(drug_4, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_4)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_4, .after = drc4)

#60112 missing - ?antidiabetic but doesn't seem to correspond to identifiable drug

#drug 5

drug_5 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc5[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_5 <- append(drug_5, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_5)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_5, .after = drc5)

#drug 6
drug_6 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc6[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_6 <- append(drug_6, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_6)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_6, .after = drc6)

#40742 missing - not clear what this is

#drug 7
drug_7 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc7[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_7 <- append(drug_7, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_7)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_7, .after = drc7)


#drug 8
drug_8 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc8[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_8 <- append(drug_8, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_8)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_8, .after = drc8)

#drug 9
drug_9 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc9[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_9 <- append(drug_9, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_9)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_9, .after = drc9)

#drug 10
drug_10 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc10[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_10 <- append(drug_10, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_10)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_10, .after = drc10)

#drug 11
drug_11 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc11[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_11 <- append(drug_11, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_11)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_11, .after = drc11)

#drug 12
drug_12 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc12[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_12 <- append(drug_12, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_12)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_12, .after = drc12)

#drug 13
drug_13 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc13[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_13 <- append(drug_13, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_13)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_13, .after = drc13)

#drug 14
drug_14 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc14[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_14 <- append(drug_14, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_14)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_14, .after = drc14)

#30105 missing - probably some form of bronchodilator but no seciton in bnf

#drug 15
drug_15 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc15[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_15 <- append(drug_15, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_15)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_15, .after = drc15)

#drug 16
drug_16 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc16[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_16 <- append(drug_16, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_16)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_16, .after = drc16)

#drug 17
drug_17 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc17[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_17 <- append(drug_17, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_17)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_17, .after = drc17)

#drug 18
drug_18 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc18[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_18 <- append(drug_18, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_18)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_18, .after = drc18)

#drug 19
drug_19 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc19[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_19 <- append(drug_19, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_19)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_19, .after = drc19)

#drug 20
drug_20 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc20[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_20 <- append(drug_20, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_20)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_20, .after = drc20)

#drug 21
drug_21 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc21[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_21 <- append(drug_21, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_21)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_21, .after = drc21)

#drug 22
drug_22 <- c()

for (i in (1:nrow(wave_8_9_drugs))) {
  a <- wave_8_9_drugs$drc22[i]
  b <- drug_list %>% dplyr::filter(grepl(a, drug_codes))
  c <- b$drug_classes[1]
  drug_22 <- append(drug_22, c)
  rm(a)
  rm(b)
  rm(c)
}

wave_8_9_drugs_1 <- cbind(wave_8_9_drugs_1, drug_22)

wave_8_9_drugs_1 <- wave_8_9_drugs_1 %>% relocate(drug_22, .after = drc22)

#add in data about whether taken in the last week
vars <- c("idauniq", "medbia")

for (i in (2:22)) {
  a <- paste0("medbia",i)
  vars <- append(vars, a)
  rm(a)
}

#create df of just whether taken in the last week

w_8_9_week <- subset(wave_8_9_nurse_data, select = vars)

#combine with wave_6_drugs_1

wave_8_9_drugs_take <- merge(wave_8_9_drugs_1, w_8_9_week, by="idauniq")

#rename medbia to medbia1

names(wave_8_9_drugs_take)[names(wave_8_9_drugs_take) == 'medbia'] <- 'medbia1'

#reorder columns so can see if taking

for (i in (1:22)) {
  wave_8_9_drugs_take <- wave_8_9_drugs_take %>% relocate(paste0("medbia",i), .after = paste0("drug_",i))
}


#statins

statins <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take, statins)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(statins = case_when(
    drug_1 %in% "Statins" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Statins" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Statins" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Statins" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Statins" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Statins" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Statins" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Statins" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Statins" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Statins" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Statins" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Statins" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Statins" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Statins" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Statins" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Statins" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Statins" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Statins" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Statins" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Statins" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Statins" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Statins" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(statins)))


#calcium channel blockers

ccb <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, ccb)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(ccb = case_when(
    drug_1 %in% "Calcium-channel blockers" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Calcium-channel blockers" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Calcium-channel blockers" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Calcium-channel blockers" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Calcium-channel blockers" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Calcium-channel blockers" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Calcium-channel blockers" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Calcium-channel blockers" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Calcium-channel blockers" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Calcium-channel blockers" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Calcium-channel blockers" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Calcium-channel blockers" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Calcium-channel blockers" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Calcium-channel blockers" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Calcium-channel blockers" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Calcium-channel blockers" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Calcium-channel blockers" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Calcium-channel blockers" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Calcium-channel blockers" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Calcium-channel blockers" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Calcium-channel blockers" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Calcium-channel blockers" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(ccb)))

#ppi

ppi <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, ppi)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(ppi = case_when(
    drug_1 %in% "Proton pump inhibitors" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Proton pump inhibitors" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Proton pump inhibitors" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Proton pump inhibitors" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Proton pump inhibitors" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Proton pump inhibitors" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Proton pump inhibitors" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Proton pump inhibitors" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Proton pump inhibitors" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Proton pump inhibitors" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Proton pump inhibitors" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Proton pump inhibitors" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Proton pump inhibitors" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Proton pump inhibitors" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Proton pump inhibitors" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Proton pump inhibitors" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Proton pump inhibitors" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Proton pump inhibitors" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Proton pump inhibitors" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Proton pump inhibitors" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Proton pump inhibitors" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Proton pump inhibitors" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(ppi)))

#Angiotensin-converting enzyme (ACE) inhibitors

acei <- c(rep(FALSE, nrow(wave_8_9_drugs_take_class)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, acei)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(acei = case_when(
    drug_1 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Angiotensin-converting enzyme (ACE) inhibitors" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(acei)))

#Antiplatelet drugs

antiplatelets <- c(rep(FALSE, nrow(wave_8_9_drugs_take_class)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, antiplatelets)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(antiplatelets = case_when(
    drug_1 %in% "Antiplatelet drugs" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Antiplatelet drugs" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Antiplatelet drugs" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Antiplatelet drugs" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Antiplatelet drugs" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Antiplatelet drugs" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Antiplatelet drugs" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Antiplatelet drugs" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Antiplatelet drugs" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Antiplatelet drugs" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Antiplatelet drugs" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Antiplatelet drugs" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Antiplatelet drugs" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Antiplatelet drugs" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Antiplatelet drugs" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Antiplatelet drugs" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Antiplatelet drugs" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Antiplatelet drugs" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Antiplatelet drugs" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Antiplatelet drugs" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Antiplatelet drugs" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Antiplatelet drugs" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(antiplatelets)))

#thyroid hormones- nb this could include meds for hypo or hyper so will ignore

#Beta-adrenoceptor blocking drugs

bb <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, bb)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(bb = case_when(
    drug_1 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Beta-adrenoceptor blocking drugs" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(bb)))

#arb

arb <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, arb)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(arb = case_when(
    drug_1 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Angiotensin II receptor antagonists" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(arb)))


#Oral anticoagulants

anticoag <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, anticoag)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(anticoag = case_when(
    drug_1 %in% "Oral anticoagulants" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Oral anticoagulants" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Oral anticoagulants" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Oral anticoagulants" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Oral anticoagulants" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Oral anticoagulants" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Oral anticoagulants" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Oral anticoagulants" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Oral anticoagulants" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Oral anticoagulants" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Oral anticoagulants" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Oral anticoagulants" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Oral anticoagulants" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Oral anticoagulants" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Oral anticoagulants" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Oral anticoagulants" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Oral anticoagulants" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Oral anticoagulants" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Oral anticoagulants" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Oral anticoagulants" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Oral anticoagulants" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Oral anticoagulants" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(anticoag)))

#Selective serotonin re-uptake inhibitors

ssri <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, ssri)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(ssri = case_when(
    drug_1 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Selective serotonin re-uptake inhibitors" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(ssri)))

#Vitamin D

vit_d_drug <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, vit_d_drug)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(vit_d_drug = case_when(
    drug_1 %in% "Vitamin D" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Vitamin D" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Vitamin D" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Vitamin D" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Vitamin D" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Vitamin D" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Vitamin D" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Vitamin D" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Vitamin D" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Vitamin D" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Vitamin D" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Vitamin D" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Vitamin D" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Vitamin D" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Vitamin D" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Vitamin D" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Vitamin D" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Vitamin D" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Vitamin D" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Vitamin D" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Vitamin D" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Vitamin D" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(vit_d_drug)))

#Biguanides (e.g. Metformin)

metformin <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, metformin)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(metformin = case_when(
    drug_1 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Biguanides (e.g. Metformin)" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(metformin)))

#Thiazides and related diuretics

thiazides <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, thiazides)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(thiazides = case_when(
    drug_1 %in% "Thiazides and related diuretics" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Thiazides and related diuretics" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Thiazides and related diuretics" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Thiazides and related diuretics" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Thiazides and related diuretics" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Thiazides and related diuretics" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Thiazides and related diuretics" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Thiazides and related diuretics" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Thiazides and related diuretics" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Thiazides and related diuretics" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Thiazides and related diuretics" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Thiazides and related diuretics" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Thiazides and related diuretics" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Thiazides and related diuretics" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Thiazides and related diuretics" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Thiazides and related diuretics" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Thiazides and related diuretics" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Thiazides and related diuretics" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Thiazides and related diuretics" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Thiazides and related diuretics" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Thiazides and related diuretics" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Thiazides and related diuretics" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(thiazides)))

#Adrenoceptor agonists

adr_ag <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, adr_ag)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(adr_ag = case_when(
    drug_1 %in% "Adrenoceptor agonists" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Adrenoceptor agonists" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Adrenoceptor agonists" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Adrenoceptor agonists" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Adrenoceptor agonists" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Adrenoceptor agonists" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Adrenoceptor agonists" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Adrenoceptor agonists" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Adrenoceptor agonists" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Adrenoceptor agonists" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Adrenoceptor agonists" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Adrenoceptor agonists" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Adrenoceptor agonists" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Adrenoceptor agonists" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Adrenoceptor agonists" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Adrenoceptor agonists" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Adrenoceptor agonists" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Adrenoceptor agonists" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Adrenoceptor agonists" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Adrenoceptor agonists" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Adrenoceptor agonists" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Adrenoceptor agonists" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(adr_ag)))

#Tricyclic and related antidepressant drugs

tca <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, tca)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(tca = case_when(
    drug_1 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Tricyclic and related antidepressant drugs" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(tca)))

#Non-opioid analgesics and compound preparations

non_op_an <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, non_op_an)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(non_op_an = case_when(
    drug_1 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Non-opioid analgesics and compound preparations" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(non_op_an)))

#Corticosteroids (respiratory)
inh_cort <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, inh_cort)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(inh_cort = case_when(
    drug_1 %in% "Corticosteroids (respiratory)" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Corticosteroids (respiratory)" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Corticosteroids (respiratory)" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Corticosteroids (respiratory)" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Corticosteroids (respiratory)" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Corticosteroids (respiratory)" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Corticosteroids (respiratory)" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Corticosteroids (respiratory)" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Corticosteroids (respiratory)" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Corticosteroids (respiratory)" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Corticosteroids (respiratory)" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Corticosteroids (respiratory)" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Corticosteroids (respiratory)" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Corticosteroids (respiratory)" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Corticosteroids (respiratory)" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Corticosteroids (respiratory)" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Corticosteroids (respiratory)" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Corticosteroids (respiratory)" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Corticosteroids (respiratory)" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Corticosteroids (respiratory)" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Corticosteroids (respiratory)" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Corticosteroids (respiratory)" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(inh_cort)))

#Non-steroidal anti-inflammatory drugs

nsaids <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, nsaids)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(nsaids = case_when(
    drug_1 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Non-steroidal anti-inflammatory drugs" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(nsaids)))

#Bisphosphonates and other drugs
bisphos <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, bisphos)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(bisphos = case_when(
    drug_1 %in% "Bisphosphonates and other drugs" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Bisphosphonates and other drugs" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Bisphosphonates and other drugs" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Bisphosphonates and other drugs" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Bisphosphonates and other drugs" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Bisphosphonates and other drugs" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Bisphosphonates and other drugs" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Bisphosphonates and other drugs" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Bisphosphonates and other drugs" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Bisphosphonates and other drugs" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Bisphosphonates and other drugs" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Bisphosphonates and other drugs" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Bisphosphonates and other drugs" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Bisphosphonates and other drugs" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Bisphosphonates and other drugs" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Bisphosphonates and other drugs" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Bisphosphonates and other drugs" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Bisphosphonates and other drugs" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Bisphosphonates and other drugs" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Bisphosphonates and other drugs" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Bisphosphonates and other drugs" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Bisphosphonates and other drugs" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(bisphos)))

#Gout treatment

gout_rx <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, gout_rx)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(gout_rx = case_when(
    drug_1 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Gout and cytotoxic induced hyperuicaemia" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(gout_rx)))

#Drugs for urinary retention

ur_ret <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, ur_ret)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(ur_ret = case_when(
    drug_1 %in% "Drugs for urinary retention" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Drugs for urinary retention" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Drugs for urinary retention" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Drugs for urinary retention" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Drugs for urinary retention" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Drugs for urinary retention" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Drugs for urinary retention" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Drugs for urinary retention" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Drugs for urinary retention" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Drugs for urinary retention" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Drugs for urinary retention" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Drugs for urinary retention" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Drugs for urinary retention" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Drugs for urinary retention" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Drugs for urinary retention" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Drugs for urinary retention" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Drugs for urinary retention" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Drugs for urinary retention" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Drugs for urinary retention" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Drugs for urinary retention" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Drugs for urinary retention" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Drugs for urinary retention" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(ur_ret)))

#Treatment of glaucoma

glauc <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, glauc)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(glauc = case_when(
    drug_1 %in% "Treatment of glaucoma" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Treatment of glaucoma" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Treatment of glaucoma" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Treatment of glaucoma" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Treatment of glaucoma" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Treatment of glaucoma" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Treatment of glaucoma" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Treatment of glaucoma" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Treatment of glaucoma" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Treatment of glaucoma" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Treatment of glaucoma" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Treatment of glaucoma" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Treatment of glaucoma" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Treatment of glaucoma" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Treatment of glaucoma" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Treatment of glaucoma" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Treatment of glaucoma" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Treatment of glaucoma" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Treatment of glaucoma" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Treatment of glaucoma" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Treatment of glaucoma" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Treatment of glaucoma" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(glauc)))

#Opioid analgesics

opioids <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, opioids)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(opioids = case_when(
    drug_1 %in% "Opioid analgesics" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Opioid analgesics" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Opioid analgesics" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Opioid analgesics" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Opioid analgesics" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Opioid analgesics" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Opioid analgesics" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Opioid analgesics" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Opioid analgesics" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Opioid analgesics" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Opioid analgesics" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Opioid analgesics" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Opioid analgesics" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Opioid analgesics" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Opioid analgesics" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Opioid analgesics" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Opioid analgesics" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Opioid analgesics" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Opioid analgesics" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Opioid analgesics" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Opioid analgesics" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Opioid analgesics" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(opioids)))

#Control of epilepsy
epilep <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, epilep)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(epilep = case_when(
    drug_1 %in% "Control of epilepsy" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Control of epilepsy" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Control of epilepsy" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Control of epilepsy" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Control of epilepsy" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Control of epilepsy" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Control of epilepsy" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Control of epilepsy" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Control of epilepsy" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Control of epilepsy" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Control of epilepsy" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Control of epilepsy" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Control of epilepsy" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Control of epilepsy" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Control of epilepsy" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Control of epilepsy" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Control of epilepsy" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Control of epilepsy" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Control of epilepsy" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Control of epilepsy" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Control of epilepsy" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Control of epilepsy" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(epilep)))

#Female sex hormones and their modulators

oestr <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, oestr)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(oestr = case_when(
    drug_1 %in% "Female sex hormones and their modulators" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Female sex hormones and their modulators" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Female sex hormones and their modulators" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Female sex hormones and their modulators" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Female sex hormones and their modulators" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Female sex hormones and their modulators" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Female sex hormones and their modulators" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Female sex hormones and their modulators" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Female sex hormones and their modulators" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Female sex hormones and their modulators" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Female sex hormones and their modulators" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Female sex hormones and their modulators" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Female sex hormones and their modulators" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Female sex hormones and their modulators" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Female sex hormones and their modulators" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Female sex hormones and their modulators" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Female sex hormones and their modulators" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Female sex hormones and their modulators" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Female sex hormones and their modulators" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Female sex hormones and their modulators" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Female sex hormones and their modulators" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Female sex hormones and their modulators" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(oestr)))

#Loop diuretics
loop <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, loop)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(loop = case_when(
    drug_1 %in% "Loop diuretics" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Loop diuretics" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Loop diuretics" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Loop diuretics" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Loop diuretics" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Loop diuretics" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Loop diuretics" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Loop diuretics" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Loop diuretics" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Loop diuretics" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Loop diuretics" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Loop diuretics" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Loop diuretics" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Loop diuretics" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Loop diuretics" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Loop diuretics" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Loop diuretics" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Loop diuretics" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Loop diuretics" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Loop diuretics" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Loop diuretics" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Loop diuretics" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(loop)))


#sulphonylureas
sulph <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, sulph)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(sulph = case_when(
    drug_1 %in% "Sulphonylureas" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Sulphonylureas" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Sulphonylureas" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Sulphonylureas" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Sulphonylureas" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Sulphonylureas" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Sulphonylureas" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Sulphonylureas" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Sulphonylureas" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Sulphonylureas" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Sulphonylureas" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Sulphonylureas" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Sulphonylureas" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Sulphonylureas" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Sulphonylureas" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Sulphonylureas" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Sulphonylureas" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Sulphonylureas" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Sulphonylureas" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Sulphonylureas" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Sulphonylureas" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Sulphonylureas" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(sulph)))


#Antihistamines
antihist <- c(rep(FALSE, nrow(wave_8_9_drugs_take)))

wave_8_9_drugs_take_class <- cbind(wave_8_9_drugs_take_class, antihist)

wave_8_9_drugs_take_class <- wave_8_9_drugs_take_class %>% 
  mutate(antihist = case_when(
    drug_1 %in% "Antihistamines" & 1 %in% medbia1 ~ TRUE,
    drug_2 %in% "Antihistamines" & 1 %in% medbia2 ~ TRUE, 
    drug_3 %in% "Antihistamines" & 1 %in% medbia3 ~ TRUE,
    drug_4 %in% "Antihistamines" & 1 %in% medbia4 ~ TRUE,
    drug_5 %in% "Antihistamines" & 1 %in% medbia5 ~ TRUE,
    drug_6 %in% "Antihistamines" & 1 %in% medbia6 ~ TRUE,
    drug_7 %in% "Antihistamines" & 1 %in% medbia7 ~ TRUE,
    drug_8 %in% "Antihistamines" & 1 %in% medbia8 ~ TRUE,
    drug_9 %in% "Antihistamines" & 1 %in% medbia9 ~ TRUE,
    drug_10 %in% "Antihistamines" & 1 %in% medbia10 ~ TRUE,
    drug_11 %in% "Antihistamines" & 1 %in% medbia11 ~ TRUE,
    drug_12 %in% "Antihistamines" & 1 %in% medbia12 ~ TRUE,
    drug_13 %in% "Antihistamines" & 1 %in% medbia13 ~ TRUE,
    drug_14 %in% "Antihistamines" & 1 %in% medbia14 ~ TRUE,
    drug_15 %in% "Antihistamines" & 1 %in% medbia15 ~ TRUE,
    drug_16 %in% "Antihistamines" & 1 %in% medbia16 ~ TRUE,
    drug_17 %in% "Antihistamines" & 1 %in% medbia17 ~ TRUE,
    drug_18 %in% "Antihistamines" & 1 %in% medbia18 ~ TRUE,
    drug_19 %in% "Antihistamines" & 1 %in% medbia19 ~ TRUE,
    drug_20 %in% "Antihistamines" & 1 %in% medbia20 ~ TRUE,
    drug_21 %in% "Antihistamines" & 1 %in% medbia21 ~ TRUE,
    drug_22 %in% "Antihistamines" & 1 %in% medbia22 ~ TRUE,
    TRUE ~ as.logical(antihist)))


drugs <- c(wave_8_9_drugs_take_class$drug_1, wave_8_9_drugs_take_class$drug_2, wave_8_9_drugs_take_class$drug_3, wave_8_9_drugs_take_class$drug_4, wave_8_9_drugs_take_class$drug_5, wave_8_9_drugs_take_class$drug_6, wave_8_9_drugs_take_class$drug_7, wave_8_9_drugs_take_class$drug_8, wave_8_9_drugs_take_class$drug_9, wave_8_9_drugs_take_class$drug_10, wave_8_9_drugs_take_class$drug_11, wave_8_9_drugs_take_class$drug_12, wave_8_9_drugs_take_class$drug_13, wave_8_9_drugs_take_class$drug_14, wave_8_9_drugs_take_class$drug_15, wave_8_9_drugs_take_class$drug_16, wave_8_9_drugs_take_class$drug_17, wave_8_9_drugs_take_class$drug_18, wave_8_9_drugs_take_class$drug_19, wave_8_9_drugs_take_class$drug_20, wave_8_9_drugs_take_class$drug_21, 
           wave_8_9_drugs_take_class$drug_22)

#get df with drug classes
w8_9_drug_classes <- subset(wave_8_9_drugs_take_class, select = c("idauniq", "statins", "ccb", "ppi",         
                                                                  "acei", "antiplatelets", "bb", "arb", "anticoag",   
                                                                  "ssri", "vit_d_drug", "metformin", "adr_ag", "tca",         
                                                                  "non_op_an", "inh_cort", "nsaids", "bisphos", "gout_rx",   
                                                                  "ur_ret", "glauc", "opioids", "epilep", "oestr",        
                                                                  "loop"))

#export wave_8_9_drug_classes
setwd("N:/My Documents/ELSA_data/upgrade")

write.csv(w8_9_drug_classes, "wave_8_9_drug_classes.csv", row.names = FALSE)



