#script to combine wave 8 and 9 comorbidity data with mutation data
library(dplyr)
library(tidyr)

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


#select subset of chip variant variables

chip_simp_prog <- subset(chip_variants, select = c("Sample", "Gene.refGene", "NonsynOI", "Otherinfo10", "AF", "AD1", "AD2"))

#merge with elsa IDs

chip_simp_merged_prog <- merge(chip_simp_prog, cases_simp, by="Sample")
chip_simp_merged_prog2 <- subset(chip_simp_merged_prog, select = c("idauniq", "NonsynOI", "Gene.refGene", "Sample", "elsa_wave", "AF", "AD1", "AD2"))

#drop those without elsa id
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

#assign vaf with 0.01 for relevant control waves
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

#remove duplicates/splice variants in different loci without protein change identifier - these have been manually reviewed from raw mutation data

ids_to_remove <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11) #censored for GDPR

chip_cases_list2 <- chip_cases_list %>% dplyr::filter(!idauniq %in% ids_to_remove)
         
#separate those with multiple variants as observations not independent
mult <- c(12, 13) #censored for GDPR
chip_cases_list3 <- chip_cases_list2 %>% dplyr::filter(idauniq %in% mult)

chip_cases_list3$num <- NA

for (i in (1:nrow(chip_cases_list3))) {
  chip_cases_list3$num[i] <- length(unlist(chip_cases_list3$AF_8[i]))
}

#keep records from samples with multiple variants with most replicates

chip_cases_mult_max <- chip_cases_list3 %>% group_by(idauniq)%>% slice(which.max(num))

chip_cases_sing <- chip_cases_list2 %>% dplyr::filter(!idauniq %in% mult)

chip_cases_mult_max <- subset(chip_cases_mult_max, select = -c(num))

#recombine
chip_cases_ind <- rbind(chip_cases_sing, chip_cases_mult_max)

#now calculated pooled sd
reps1 <- unlist(chip_cases_ind$AF_8[1])
sd1 <- sd(reps1)
n1 <- length(reps1)

reps2 <- unlist(chip_cases_ind$AF_8[2])
sd2 <- sd(reps2)
n2 <- length(reps2)

reps3 <- unlist(chip_cases_ind$AF_8[3])
sd3 <- sd(reps3)
n3 <- length(reps3)

reps4 <- unlist(chip_cases_ind$AF_6[4])
sd4 <- sd(reps4)
n4 <- length(reps4)

reps5 <- unlist(chip_cases_ind$AF_6[5])
sd5 <- sd(reps5)
n5 <- length(reps5)

reps6 <- unlist(chip_cases_ind$AF_6[6])
sd6 <- sd(reps6)
n6 <- length(reps6)

reps7 <- unlist(chip_cases_ind$AF_6[7])
sd7 <- sd(reps7)
n7 <- length(reps7)

reps8 <- unlist(chip_cases_ind$AF_8[8])
sd8 <- sd(reps8)
n8 <- length(reps8)

reps9 <- unlist(chip_cases_ind$AF_8[9])
sd9 <- sd(reps9)
n9 <- length(reps9)

#calculate pooled sd

pooled_sd <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2 + (n3-1)*sd3^2 + (n4-1)*sd4^2 + (n5-1)*sd5^2 + (n6-1)*sd6^2 + (n7-1)*sd7^2 + (n8-1)*sd8^2 + (n9-1)*sd9^2)/(n1 + n2 + n3 + n4 + n5 + n6 + n7 + n8 + n9 - 9))



