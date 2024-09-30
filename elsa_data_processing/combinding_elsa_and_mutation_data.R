#script to combine wave 8 and 9 comorbidity data with mutation data
library(dplyr)
library(tidyr)
library(readxl)
library(writexl)


setwd("N:/My Documents/ELSA_data/upgrade/")

late_wave_comorbs <- read.csv("ELSA_wave_8_9_comorbidities_all.csv", header = TRUE)

#import chip case mutation data and chip sample data
setwd("N:/My Documents/ELSA_data/")

chip_variants <- readxl::read_xlsx("chip_cases_whitelist_elsa.xlsx")

cases <- readxl::read_xlsx("Waves 2,4,6,8,9 Lab IDs and Idauniq -Chip cases.xlsx")

#import control data
controls <- readxl::read_xlsx("Waves 2,4,6,8,9 Lab IDs and Idauniq -Control cases.xlsx")

#table chip cases per wave
table(cases$elsa_wave)

cases_simp <- subset(cases, select=-c(bldrec, consn, CONBST))

#consider multiallelic calls separately

multiallelic <- chip_variants %>% dplyr::filter(grepl("multiallelic", Otherinfo10))

chip_variants_mono <- chip_variants %>% dplyr::filter(!grepl("multiallelic", Otherinfo10))

#adjust multiallelic to take maximum AF
multiallelic <- tidyr::separate(multiallelic, Otherinfo13_2, into = c("allele", "depth", "vaf", "extra1", "extra2", "extra3", "extra4"), sep = ":", remove = FALSE, fill = "left")
multiallelic <- subset(multiallelic, select = -c(allele, depth, extra1, extra2, extra3, extra4))

#separate depth colunm
multiallelic <- tidyr::separate(multiallelic, vaf, into = c("vaf_1", "vaf_2", "vaf_3", "vaf_4"), sep = ",", remove = FALSE, fill = "right")

#fill in AF and derived depth variables
multiallelic$AF <- as.numeric(pmax(as.numeric(multiallelic$vaf_1), as.numeric(multiallelic$vaf_2), as.numeric(multiallelic$vaf_3), as.numeric(multiallelic$vaf_4), na.rm = TRUE))

multiallelic$AD1 <- as.integer(multiallelic$DP*(1-multiallelic$AF))
multiallelic$AD2 <- as.integer(multiallelic$DP*multiallelic$AF)

#drop extra variables created
multiallelic <- subset(multiallelic, select = -c(vaf, vaf_1, vaf_2, vaf_3, vaf_4))

#recombine with mono calls

chip_variants <- rbind(chip_variants_mono, multiallelic)

#export these data
setwd("N:/My Documents/ELSA_data/upgrade")
write.table(chip_variants, "anon_chip_variant_data_clean.txt", row.names = FALSE, sep = "\t")

#select subset of chip variant variables

chip_simp_prog <- subset(chip_variants, select = c("Sample", "Gene.refGene", "NonsynOI", "Otherinfo10", "AF", "AD1", "AD2"))

#merge with elsa IDs

chip_simp_merged_prog <- merge(chip_simp_prog, cases_simp, by="Sample")
chip_simp_merged_prog2 <- subset(chip_simp_merged_prog, select = c("idauniq", "NonsynOI", "Gene.refGene", "Sample", "elsa_wave", "AF", "AD1", "AD2"))

#drop those with NA for elsa id
chip_simp_merged_prog2 <- chip_simp_merged_prog2 %>% dplyr::filter(!is.na(idauniq))


#reshape long to wide
chip_merged_prog_wide <- pivot_wider(chip_simp_merged_prog2, id_cols = c("idauniq", "NonsynOI", "Gene.refGene"), names_from = "elsa_wave", values_from = c("AF", "Sample", "AD1", "AD2"))

#convert nulls to NAs

chip_merged_prog_wide <- chip_merged_prog_wide %>% mutate(across(AF_9:AD2_6, ~replace(., lengths(.)== 0, NA)))

#move variables into correct order

chip_merged_prog_wide <- chip_merged_prog_wide %>% relocate(AF_8, .after = AF_6)
chip_merged_prog_wide <- chip_merged_prog_wide %>% relocate(AF_9, .after = AF_8)

chip_merged_prog_wide <- chip_merged_prog_wide %>% relocate(Sample_8, .after = Sample_6)
chip_merged_prog_wide <- chip_merged_prog_wide %>% relocate(Sample_9, .after = Sample_8)

chip_merged_prog_wide <- chip_merged_prog_wide %>% relocate(AD1_8, .after = AD1_6)
chip_merged_prog_wide <- chip_merged_prog_wide %>% relocate(AD1_9, .after = AD1_8)

chip_merged_prog_wide <- chip_merged_prog_wide %>% relocate(AD2_8, .after = AD2_6)
chip_merged_prog_wide <- chip_merged_prog_wide %>% relocate(AD2_9, .after = AD2_8)

#remove cases of chip identified in pilot using differing dna conc and different ifc
cases_to_check <- subset(cases, select = c(idauniq, elsa_wave))
colnames(cases_to_check)[2] <- "elsa_wave_case"
cases_to_check <- cases_to_check %>% dplyr::filter(!is.na(idauniq))

ctls_to_check <- subset(controls, select = c(idauniq, elsa_wave))
colnames(ctls_to_check)[2] <- "elsa_wave_control"
ctls_to_check <- ctls_to_check %>% dplyr::filter(!is.na(idauniq))

case_ctl_overlap <- merge(cases_to_check, ctls_to_check, by = "idauniq")

#same wave sample identified as case and control due to pilot optimisation - can't be sure of status
case_ctl_overlap$overlap <- ifelse(case_ctl_overlap$elsa_wave_case == case_ctl_overlap$elsa_wave_control, TRUE, FALSE)

case_ctl_overlap <- case_ctl_overlap %>% dplyr::filter(overlap == TRUE)

case_ctl_overlap <- unique(case_ctl_overlap)


chip_merged_post_excl <- chip_merged_prog_wide %>% dplyr::filter(!idauniq %in% case_ctl_overlap$idauniq)

#add in control data - individuals identified as controls in a particular wave, designate vaf (conservatively) as 0.01
#get overlaps
chip_merged_post_excl$overlap <- ifelse(chip_merged_post_excl$idauniq %in% controls$idauniq, TRUE, FALSE)

overlaps <- chip_merged_post_excl %>% dplyr::filter(overlap == TRUE)

non_overlaps_cases <- chip_merged_post_excl %>% dplyr::filter(overlap == FALSE)
#seperate controls into waves
controls_w2 <- controls %>% dplyr::filter(elsa_wave == 2)
controls_w4 <- controls %>% dplyr::filter(elsa_wave == 4)
controls_w6 <- controls %>% dplyr::filter(elsa_wave == 6)
controls_w8 <- controls %>% dplyr::filter(elsa_wave == 8)
controls_w9 <- controls %>% dplyr::filter(elsa_wave == 9)

overlaps$ctl_w2 <- ifelse(overlaps$idauniq %in% controls_w2$idauniq, TRUE, FALSE)
overlaps$ctl_w4 <- ifelse(overlaps$idauniq %in% controls_w4$idauniq, TRUE, FALSE)
overlaps$ctl_w6 <- ifelse(overlaps$idauniq %in% controls_w6$idauniq, TRUE, FALSE)
overlaps$ctl_w8 <- ifelse(overlaps$idauniq %in% controls_w8$idauniq, TRUE, FALSE)
overlaps$ctl_w9 <- ifelse(overlaps$idauniq %in% controls_w9$idauniq, TRUE, FALSE)

#replace vaf with 0.01 for relevant control waves
overlaps$AF_2 <- ifelse(overlaps$ctl_w2 == TRUE, 0.01, overlaps$AF_2)
overlaps$AF_4 <- ifelse(overlaps$ctl_w4 == TRUE, 0.01, overlaps$AF_4)
overlaps$AF_6 <- ifelse(overlaps$ctl_w6 == TRUE, 0.01, overlaps$AF_6)
overlaps$AF_8 <- ifelse(overlaps$ctl_w8 == TRUE, 0.01, overlaps$AF_8)
overlaps$AF_9 <- ifelse(overlaps$ctl_w9 == TRUE, 0.01, overlaps$AF_9)

#drop extra vars so can merge
overlaps_clean <- subset(overlaps, select = colnames(non_overlaps_cases))

#merge
chip_cases_complete <- rbind(non_overlaps_cases, overlaps_clean)

#now need to deal with those with nested character for vaf due to multiple measurements
library(stringr)
chip_cases_list <- chip_cases_complete %>% dplyr::filter(if_any(.cols = c(AF_2, AF_4, AF_6, AF_8, AF_9), ~ grepl("c", .)))

chip_cases_no_list <- dplyr::anti_join(chip_cases_complete, chip_cases_list)

#unlist afs in chip_cases_list
#wave 2
new_AF_2 <- c()
for (i in (1:nrow(chip_cases_list))) {
  a <- mean(unlist(chip_cases_list$AF_2[i]))
  new_AF_2 <- append(new_AF_2, a)
  rm(a)
}

chip_cases_list$AF_2 <- new_AF_2 

#wave 4
new_AF_4 <- c()
for (i in (1:nrow(chip_cases_list))) {
  a <- mean(unlist(chip_cases_list$AF_4[i]))
  new_AF_4 <- append(new_AF_4, a)
  rm(a)
}

chip_cases_list$AF_4 <- new_AF_4 

#wave 6

new_AF_6 <- c()
for (i in (1:nrow(chip_cases_list))) {
  a <- mean(unlist(chip_cases_list$AF_6[i]))
  new_AF_6 <- append(new_AF_6, a)
  rm(a)
}

chip_cases_list$AF_6 <- new_AF_6 

#wave 8

new_AF_8 <- c()
for (i in (1:nrow(chip_cases_list))) {
  a <- mean(unlist(chip_cases_list$AF_8[i]))
  new_AF_8 <- append(new_AF_8, a)
  rm(a)
}

chip_cases_list$AF_8 <- new_AF_8 

#wave 9

new_AF_9 <- c()
for (i in (1:nrow(chip_cases_list))) {
  a <- mean(unlist(chip_cases_list$AF_9[i]))
  new_AF_9 <- append(new_AF_9, a)
  rm(a)
}

chip_cases_list$AF_9 <- new_AF_9 

#now re combine with other non listed data

chip_cases_clean <- rbind(chip_cases_no_list, chip_cases_list)
chip_cases_clean_simp <- subset(chip_cases_clean, select = c(idauniq, Gene.refGene, NonsynOI, AF_2, AF_4, AF_6, AF_8, AF_9))

#get individuals with repeat variants 
#count the NAs

chip_cases_clean_simp$missing <- rowSums(is.na(chip_cases_clean_simp))

#combine with late wave comorbidities

chip_cases_comorb <- merge(chip_cases_clean_simp, late_wave_comorbs, by = "idauniq", all.x = TRUE)

#combine with wave 6 medication data
setwd("N:/My Documents/ELSA_data/drug_regression_analysis")
wave_6_drug_info <- read.table("wave_6_drug_categories.txt", sep = "\t", header = TRUE)

chip_cases_comorb_drugs <- merge(chip_cases_comorb, wave_6_drug_info, by = "idauniq", all.x = TRUE)

#add variant type
chip_cases_comorb_drugs$variant_type[grep("fs*", chip_cases_comorb_drugs$NonsynOI)] <- "nonsense"
chip_cases_comorb_drugs$variant_type[grep("X", chip_cases_comorb_drugs$NonsynOI)] <- "nonsense"
chip_cases_comorb_drugs$variant_type[grep("nan", chip_cases_comorb_drugs$NonsynOI)] <- "splice_variant"
chip_cases_comorb_drugs$variant_type[grep("del", chip_cases_comorb_drugs$NonsynOI)] <- "indel"
chip_cases_comorb_drugs$variant_type[grep("\\*$", chip_cases_comorb_drugs$NonsynOI)] <- "nonsense"
chip_cases_comorb_drugs$variant_type[is.na(chip_cases_comorb_drugs$variant_type)] <- "snv"

#add variable of most recent vaf
chip_cases_comorb_drugs$recent_vaf <- chip_cases_comorb_drugs$AF_9
chip_cases_comorb_drugs$recent_vaf <- ifelse(is.na(chip_cases_comorb_drugs$recent_vaf), chip_cases_comorb_drugs$AF_8, chip_cases_comorb_drugs$recent_vaf)
chip_cases_comorb_drugs$recent_vaf <- ifelse(is.na(chip_cases_comorb_drugs$recent_vaf), chip_cases_comorb_drugs$AF_6, chip_cases_comorb_drugs$recent_vaf)
chip_cases_comorb_drugs$recent_vaf <- ifelse(is.na(chip_cases_comorb_drugs$recent_vaf), chip_cases_comorb_drugs$AF_4, chip_cases_comorb_drugs$recent_vaf)
chip_cases_comorb_drugs$recent_vaf <- ifelse(is.na(chip_cases_comorb_drugs$recent_vaf), chip_cases_comorb_drugs$AF_2, chip_cases_comorb_drugs$recent_vaf)


#remove duplicate ASXL1 calls at 646 
chip_cases_comorb_drugs_nonasxl1 <- chip_cases_comorb_drugs %>% dplyr::filter(!grepl("ASXL1", Gene.refGene))
chip_cases_comorb_drugs_asxl1 <- chip_cases_comorb_drugs %>% dplyr::filter(grepl("ASXL1", Gene.refGene))

chip_cases_comorb_drugs_asxl1 <- chip_cases_comorb_drugs_asxl1 %>% distinct(idauniq, Gene.refGene, recent_vaf, .keep_all = TRUE)

chip_cases_comorb_drugs_asxl1fixed <- rbind(chip_cases_comorb_drugs_nonasxl1, chip_cases_comorb_drugs_asxl1)


#now just get info on those with repeat variants
chip_repeats <- dplyr::filter(chip_cases_comorb_drugs_asxl1fixed, missing <= 3)

table(chip_repeats$Gene.refGene)

#add in rate of clonal progression for repeat variants 
#delta vaf per year for different combinations - elsa wave 9 2019, elsa wave 8 2017, elsa wave 6 2012, elsa wave 4 2008, elsa wave 2 2004)
#chip_repeats$n1 <- 15
#chip_repeats$n2 <- 13
#chip_repeats$n3 <- 11
#chip_repeats$n4 <- 9
#chip_repeats$n5 <- 5
chip_repeats$AF_2 <- as.numeric(chip_repeats$AF_2)
chip_repeats$AF_4 <- as.numeric(chip_repeats$AF_4)
chip_repeats$AF_6 <- as.numeric(chip_repeats$AF_6)
chip_repeats$AF_8 <- as.numeric(chip_repeats$AF_8)
chip_repeats$AF_9 <- as.numeric(chip_repeats$AF_9)

#calculate i, percentage growth rate, for each combination

chip_repeats$i_1 <- ((chip_repeats$AF_9 / chip_repeats$AF_2)^{1/15} - 1)*100
chip_repeats$i_2 <- ((chip_repeats$AF_8 / chip_repeats$AF_2)^{1/13} - 1)*100
chip_repeats$i_3 <- ((chip_repeats$AF_9 / chip_repeats$AF_4)^{1/11} - 1)*100
chip_repeats$i_4 <- ((chip_repeats$AF_8 / chip_repeats$AF_4)^{1/9} - 1)*100
chip_repeats$i_5 <- ((chip_repeats$AF_8 / chip_repeats$AF_6)^{1/5} - 1)*100
chip_repeats$i_6 <- ((chip_repeats$AF_6 / chip_repeats$AF_2)^{1/8} - 1)*100 


#coalesce delta vaf per year into one column
chip_repeats$growth_rate <- chip_repeats$i_1
chip_repeats$growth_rate <- ifelse(is.na(chip_repeats$growth_rate), chip_repeats$i_2, chip_repeats$growth_rate)
chip_repeats$growth_rate <- ifelse(is.na(chip_repeats$growth_rate), chip_repeats$i_3, chip_repeats$growth_rate)
chip_repeats$growth_rate <- ifelse(is.na(chip_repeats$growth_rate), chip_repeats$i_4, chip_repeats$growth_rate)
chip_repeats$growth_rate <- ifelse(is.na(chip_repeats$growth_rate), chip_repeats$i_5, chip_repeats$growth_rate)
chip_repeats$growth_rate <- ifelse(is.na(chip_repeats$growth_rate), chip_repeats$i_6, chip_repeats$growth_rate)

#drop unecessary variables
chip_repeats <- subset(chip_repeats, select = -c(i_1, i_2, i_3, i_4, i_5, i_6))

#save this data frame - all chip repeats including haem malignancy and wave 6 drug data
setwd("N:/My Documents/ELSA_data/upgrade/")


write_xlsx(chip_repeats, "chip_repeats_all_comorbs_wave6_drug.xlsx")

#now combine with wave 9 metadata

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

all_drugs <- as.data.frame(table(drugs))

#merge drug classes with the chip case data

#get df with drug classes
w8_9_drug_classes <- subset(wave_8_9_drugs_take_class, select = c("idauniq", "statins", "ccb", "ppi",         
                                                                  "acei", "antiplatelets", "bb", "arb", "anticoag",   
                                                                  "ssri", "vit_d_drug", "metformin", "adr_ag", "tca",         
                                                                  "non_op_an", "inh_cort", "nsaids", "bisphos", "gout_rx",   
                                                                  "ur_ret", "glauc", "opioids", "epilep", "oestr",        
                                                                  "loop"))


#merge with mutation and comorbidity data for those with samples in wave 8 or 9
chip_cases_comorb_8 <- chip_cases_comorb %>% dplyr::filter(!is.na(AF_8))
chip_cases_comorb_9 <- chip_cases_comorb %>% dplyr::filter(!is.na(AF_9))

chip_cases_comorb_8_and_9 <- unique(rbind(chip_cases_comorb_8, chip_cases_comorb_9))


chip_cases_comorb_drugs_w8_9 <- merge(chip_cases_comorb_8_and_9, w8_9_drug_classes, by = "idauniq", all.x = TRUE)


#now merge with blood parameters
wave_8_9_bloods <- subset(wave_8_9_nurse_data, select = c("idauniq", "sysval", "diaval", "pulval", "mapval", "chol", "hdl", "trig", "ldl", "rtin", "hscrp", "vitd", "igf1", "hba1c", "cfib", "hgb", "mch", "wbc"))

chip_cases_comorb_drugs_w8_9 <- merge(chip_cases_comorb_drugs_w8_9, wave_8_9_bloods, by = "idauniq", all.x = TRUE)

#replace negative numbers (missing data) with NA for blood and bp parameters

chip_cases_comorb_drugs_w8_9$sysval[chip_cases_comorb_drugs_w8_9$sysval <0] <- NA
chip_cases_comorb_drugs_w8_9$diaval[chip_cases_comorb_drugs_w8_9$diaval <0] <- NA
chip_cases_comorb_drugs_w8_9$pulval[chip_cases_comorb_drugs_w8_9$pulval <0] <- NA
chip_cases_comorb_drugs_w8_9$chol[chip_cases_comorb_drugs_w8_9$chol <0] <- NA
chip_cases_comorb_drugs_w8_9$hdl[chip_cases_comorb_drugs_w8_9$hdl <0] <- NA
chip_cases_comorb_drugs_w8_9$trig[chip_cases_comorb_drugs_w8_9$trig <0] <- NA
chip_cases_comorb_drugs_w8_9$ldl[chip_cases_comorb_drugs_w8_9$ldl <0] <- NA
chip_cases_comorb_drugs_w8_9$rtin[chip_cases_comorb_drugs_w8_9$rtin <0] <- NA
chip_cases_comorb_drugs_w8_9$hscrp[chip_cases_comorb_drugs_w8_9$hscrp <0] <- NA
chip_cases_comorb_drugs_w8_9$vitd[chip_cases_comorb_drugs_w8_9$vitd <0] <- NA
chip_cases_comorb_drugs_w8_9$igf1[chip_cases_comorb_drugs_w8_9$igf1 <0] <- NA
chip_cases_comorb_drugs_w8_9$hba1c[chip_cases_comorb_drugs_w8_9$hba1c <0] <- NA
chip_cases_comorb_drugs_w8_9$cfib[chip_cases_comorb_drugs_w8_9$cfib <0] <- NA
chip_cases_comorb_drugs_w8_9$hgb[chip_cases_comorb_drugs_w8_9$hgb <0] <- NA
chip_cases_comorb_drugs_w8_9$mch[chip_cases_comorb_drugs_w8_9$mch <0] <- NA
chip_cases_comorb_drugs_w8_9$wbc[chip_cases_comorb_drugs_w8_9$wbc <0] <- NA

#now save this df - all wave 8 and 9 cases with most recent comorbidities and wave 8 and 9 drugs

#need to convert AF vars to numeric
chip_cases_comorb_drugs_w8_9$AF_2 <- as.numeric(chip_cases_comorb_drugs_w8_9$AF_2)
chip_cases_comorb_drugs_w8_9$AF_4 <- as.numeric(chip_cases_comorb_drugs_w8_9$AF_4)
chip_cases_comorb_drugs_w8_9$AF_6 <- as.numeric(chip_cases_comorb_drugs_w8_9$AF_6)
chip_cases_comorb_drugs_w8_9$AF_8 <- as.numeric(chip_cases_comorb_drugs_w8_9$AF_8)
chip_cases_comorb_drugs_w8_9$AF_9 <- as.numeric(chip_cases_comorb_drugs_w8_9$AF_9)

setwd("N:/My Documents/ELSA_data/upgrade/")
write_xlsx(chip_cases_comorb_drugs_w8_9, "wave_8_and_9_chip_cases_comorbs_drugs_bloods.xlsx")

