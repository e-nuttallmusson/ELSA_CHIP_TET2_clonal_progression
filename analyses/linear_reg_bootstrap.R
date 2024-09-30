#script to assess TET2 clonal growth rate in linear regression 
library(dplyr)
library(readxl)
library(olsrr)
library(jtools)
library(broom.mixed)
library(lmtest)
library(Metrics)
library(caret)
library(ggplot2)
library(glmtoolbox)
library(truncnorm)
library(broom)
library(data.table)

setwd("N:/My Documents/ELSA_data/upgrade/tet2_vaf_progression_linear_regression/")

#read in repeat chip mutation data with wave 8/9 comorbidities and wave 6 medications

chip_repeats_comorbs_drugs <- read_xlsx("chip_repeats_all_comorbs_wave6_drug.xlsx")

#add variable of number of mutations
chip_repeats_comorbs_drugs <- chip_repeats_comorbs_drugs %>% group_by(idauniq) %>% mutate(count = n())

#select tet2 repeats
tet2_chip_repeats_comorbs_drugs <- chip_repeats_comorbs_drugs %>% dplyr::filter(grepl("TET2", Gene.refGene))

#remove haem disorder
tet2_chip_repeats_comorbs_drugs <- tet2_chip_repeats_comorbs_drugs %>% dplyr::filter(blood_dis == 0)

#remove those with missing drug data - this gives 112 observations
tet2_chip_repeats_comorbs_drugs <- tet2_chip_repeats_comorbs_drugs %>% dplyr::filter(!is.na(nsaids))

tet2_chip_repeats_comorbs_drugs <- as.data.frame(tet2_chip_repeats_comorbs_drugs)


#calculate growth rate
tet2_chip_repeats_comorbs_drugs$i_1 <- ((tet2_chip_repeats_comorbs_drugs$AF_9 / tet2_chip_repeats_comorbs_drugs$AF_2)^{1/15} - 1)*100
tet2_chip_repeats_comorbs_drugs$i_2 <- ((tet2_chip_repeats_comorbs_drugs$AF_8 / tet2_chip_repeats_comorbs_drugs$AF_2)^{1/13} - 1)*100
tet2_chip_repeats_comorbs_drugs$i_3 <- ((tet2_chip_repeats_comorbs_drugs$AF_9 / tet2_chip_repeats_comorbs_drugs$AF_4)^{1/11} - 1)*100
tet2_chip_repeats_comorbs_drugs$i_4 <- ((tet2_chip_repeats_comorbs_drugs$AF_8 / tet2_chip_repeats_comorbs_drugs$AF_4)^{1/9} - 1)*100
tet2_chip_repeats_comorbs_drugs$i_5 <- ((tet2_chip_repeats_comorbs_drugs$AF_8 / tet2_chip_repeats_comorbs_drugs$AF_6)^{1/5} - 1)*100
tet2_chip_repeats_comorbs_drugs$i_6 <- ((tet2_chip_repeats_comorbs_drugs$AF_6 / tet2_chip_repeats_comorbs_drugs$AF_2)^{1/8} - 1)*100 


#coalesce delta vaf per year into one column
tet2_chip_repeats_comorbs_drugs$growth_rate <- tet2_chip_repeats_comorbs_drugs$i_1
tet2_chip_repeats_comorbs_drugs$growth_rate <- ifelse(is.na(tet2_chip_repeats_comorbs_drugs$growth_rate), tet2_chip_repeats_comorbs_drugs$i_2, tet2_chip_repeats_comorbs_drugs$growth_rate)
tet2_chip_repeats_comorbs_drugs$growth_rate <- ifelse(is.na(tet2_chip_repeats_comorbs_drugs$growth_rate), tet2_chip_repeats_comorbs_drugs$i_3, tet2_chip_repeats_comorbs_drugs$growth_rate)
tet2_chip_repeats_comorbs_drugs$growth_rate <- ifelse(is.na(tet2_chip_repeats_comorbs_drugs$growth_rate), tet2_chip_repeats_comorbs_drugs$i_4, tet2_chip_repeats_comorbs_drugs$growth_rate)
tet2_chip_repeats_comorbs_drugs$growth_rate <- ifelse(is.na(tet2_chip_repeats_comorbs_drugs$growth_rate), tet2_chip_repeats_comorbs_drugs$i_5, tet2_chip_repeats_comorbs_drugs$growth_rate)
tet2_chip_repeats_comorbs_drugs$growth_rate <- ifelse(is.na(tet2_chip_repeats_comorbs_drugs$growth_rate), tet2_chip_repeats_comorbs_drugs$i_6, tet2_chip_repeats_comorbs_drugs$growth_rate)




#for independent observations select single variant - keep largest growth rate
tet2_chip_repeats_ordered <- tet2_chip_repeats_comorbs_drugs[order(tet2_chip_repeats_comorbs_drugs$idauniq, abs(tet2_chip_repeats_comorbs_drugs$growth_rate)),]
tet2_chip_repeats_single <- tet2_chip_repeats_ordered[ !duplicated(tet2_chip_repeats_ordered$idauniq), ]


#create variable indicating which elsa_wave control sample was
wave_2_ctls <- tet2_chip_repeats_single %>% dplyr::filter(tet2_chip_repeats_single$AF_2 == 0.01)
wave_4_ctls <- tet2_chip_repeats_single %>% dplyr::filter(tet2_chip_repeats_single$AF_4 == 0.01)
wave_6_ctls <- tet2_chip_repeats_single %>% dplyr::filter(tet2_chip_repeats_single$AF_6 == 0.01)
wave_8_ctls <- tet2_chip_repeats_single %>% dplyr::filter(tet2_chip_repeats_single$AF_8 == 0.01)
wave_9_ctls <- tet2_chip_repeats_single %>% dplyr::filter(tet2_chip_repeats_single$AF_9 == 0.01)

wave_2_ctls$elsa_wave <- 2
wave_4_ctls$elsa_wave <- 4
wave_6_ctls$elsa_wave <- 6
wave_9_ctls$elsa_wave <- 9

#recombine
tet2_chip_repeats_single_ctlwave <- do.call("rbind", list(wave_2_ctls, wave_4_ctls, wave_6_ctls, wave_9_ctls))

#import control details with sample number of control
setwd("N://My Documents/ELSA_data")
control_samples <- readxl::read_xlsx("Waves 2,4,6,8,9 Lab IDs and Idauniq -Control cases.xlsx")
case_samples <- read_xlsx("Waves 2,4,6,8,9 Lab IDs and Idauniq -Chip cases.xlsx")
#simplify
control_samples <- subset(control_samples, select = c(Sample, idauniq, elsa_wave))

case_samples <- subset(case_samples, select = c(idauniq, Sample))
#merge
tet2_chip_repeats_single_ctlwave <- merge(tet2_chip_repeats_single_ctlwave, control_samples, by = c("idauniq", "elsa_wave"), all.x = TRUE)


#import mutation data to get mutation coordinates
setwd("N://My Documents/ELSA_data/elsa_variant_calling/")

variant_data <- read_xlsx("chip_cases_whitelist_elsa.xlsx")
#merge with sample info
variant_data <- merge(variant_data, case_samples, by = "Sample", all.x = TRUE)
#simplify
vars_simp <- subset(variant_data, select = c(idauniq, Gene.refGene, NonsynOI, Chr, Start, End))

#merge mutation data
tet2_chip_repeats_mutation_coords <- merge(tet2_chip_repeats_single_ctlwave, vars_simp, by = c("idauniq", "Gene.refGene", "NonsynOI"), all.x = TRUE)

#now drop idauniq to get anonymous data
tet2_chip_repeats_mutation_coords_anon <- subset(tet2_chip_repeats_mutation_coords, select = -c(idauniq))

#export this data
setwd("N://My Documents/ELSA_data/upgrade/tet2_vaf_progression_linear_regression")
write_xlsx(tet2_chip_repeats_mutation_coords_anon, "tet2_chip_repeats_mutation_anon.xlsx")

#import df with locus specific depths manually added
tet2_chip_repeats_depth <- read_xlsx("tet2_chip_repeats_mutation_anon_depth.xlsx")
#simplify
tet2_chip_repeats_depth <- subset(tet2_chip_repeats_depth, select = c(Gene.refGene, NonsynOI, Sample, growth_rate, Depth))

#add idauniq to depth df
tet2_chip_repeats_mutation_coords_simp <- subset(tet2_chip_repeats_mutation_coords, select = c(idauniq, Gene.refGene, NonsynOI, Sample, growth_rate))
tet2_chip_repeats_depth_idauniq <- merge(tet2_chip_repeats_depth, tet2_chip_repeats_mutation_coords_simp, by = "Sample")
tet2_chip_repeats_depth_idauniq <- unique(tet2_chip_repeats_depth_idauniq)
tet2_chip_repeats_depth_idauniq <- tet2_chip_repeats_depth_idauniq %>% dplyr::filter(!is.na(Depth))
tet2_chip_repeats_depth_idauniq <- subset(tet2_chip_repeats_depth_idauniq, select = c(idauniq, Depth))

tet2_chip_repeats_mutations_coords_depth <- merge(tet2_chip_repeats_mutation_coords, tet2_chip_repeats_depth_idauniq, by = "idauniq")

#one individual has technical repeats all confirming wave 6 sample is control - remove duplicates
tet2_chip_repeats_mutations_coords_depth1 <- tet2_chip_repeats_mutations_coords_depth[ !duplicated(tet2_chip_repeats_mutations_coords_depth$idauniq), ]

#make new variable which is maximum allele frequency that could be present at limit of detection. Filtering strategy: 20 alt reads needed
tet2_chip_repeats_mutations_coords_depth1$lod <- 20/tet2_chip_repeats_mutations_coords_depth1$Depth
#simplify
lod <- subset(tet2_chip_repeats_mutations_coords_depth1, select = c(idauniq, lod))

#make new df with lod 
tet2_chip_repeats_single_lod <- merge(tet2_chip_repeats_single, lod, by = "idauniq", all.x = TRUE)

#drop those obs with a lod > 0.02 (depth <1000 when designated control sample)

tet2_chip_repeats_single_lod <- tet2_chip_repeats_single_lod %>% dplyr::filter(is.na(lod) | lod < 0.02)

#108 individuals############



#add in wave 6 nurse variables - as mid point, consistent with medication data
setwd("N:/My Documents/ELSA_data/elsa_survey_data")
w6_nurse_data <- read.table("wave_6_nurse_clean.txt", sep = "\t", header = TRUE)
w6_nurse_data_simp <- subset(w6_nurse_data, select = c(idauniq, hgb, wbc, mch, hscrp, rtin, BMIVAL, chol, cfib))
tet2_chip_repeats_single_lod <- merge(tet2_chip_repeats_single_lod, w6_nurse_data_simp, by = "idauniq", all.x = TRUE)

#replace neg numbers (missing data) with median
tet2_chip_repeats_single_lod$hgb <- ifelse(tet2_chip_repeats_single_lod$hgb < 0, median(tet2_chip_repeats_single_lod$hgb, na.rm = TRUE), tet2_chip_repeats_single_lod$hgb)
tet2_chip_repeats_single_lod$wbc <- ifelse(tet2_chip_repeats_single_lod$wbc < 0, median(tet2_chip_repeats_single_lod$wbc, na.rm = TRUE), tet2_chip_repeats_single_lod$wbc)
tet2_chip_repeats_single_lod$mch <- ifelse(tet2_chip_repeats_single_lod$mch < 0, median(tet2_chip_repeats_single_lod$mch, na.rm = TRUE), tet2_chip_repeats_single_lod$mch)
tet2_chip_repeats_single_lod$hscrp <- ifelse(tet2_chip_repeats_single_lod$hscrp < 0, median(tet2_chip_repeats_single_lod$hscrp, na.rm = TRUE), tet2_chip_repeats_single_lod$hscrp)
tet2_chip_repeats_single_lod$rtin <- ifelse(tet2_chip_repeats_single_lod$rtin < 0, median(tet2_chip_repeats_single_lod$rtin, na.rm = TRUE), tet2_chip_repeats_single_lod$rtin)
tet2_chip_repeats_single_lod$chol <- ifelse(tet2_chip_repeats_single_lod$chol < 0, median(tet2_chip_repeats_single_lod$chol, na.rm = TRUE), tet2_chip_repeats_single_lod$chol)
tet2_chip_repeats_single_lod$BMIVAL <- ifelse(tet2_chip_repeats_single_lod$BMIVAL < 0, median(tet2_chip_repeats_single_lod$BMIVAL, na.rm = TRUE), tet2_chip_repeats_single_lod$BMIVAL)
tet2_chip_repeats_single_lod$cfib <- ifelse(tet2_chip_repeats_single_lod$cfib < 0, median(tet2_chip_repeats_single_lod$cfib, na.rm = TRUE), tet2_chip_repeats_single_lod$cfib)


#now do regression
#need to recheck frequency of independent variables

#check which vars sufficiently powered
table(tet2_chip_repeats_single_lod$indsex) #ok
table(tet2_chip_repeats_single_lod$inflam_dis) #too few
table(tet2_chip_repeats_single_lod$cancer) #too few
table(tet2_chip_repeats_single_lod$diabetes) #ok
table(tet2_chip_repeats_single_lod$ihd) #too few
table(tet2_chip_repeats_single_lod$heart_dis) #ok
table(tet2_chip_repeats_single_lod$stroke) #too few
table(tet2_chip_repeats_single_lod$osteoporosis) #ok
table(tet2_chip_repeats_single_lod$asthma) #ok
table(tet2_chip_repeats_single_lod$lung_disease) #too few
table(tet2_chip_repeats_single_lod$current_smoker) #too few
table(tet2_chip_repeats_single_lod$htn) #ok
table(tet2_chip_repeats_single_lod$dementia) #too few
table(tet2_chip_repeats_single_lod$nsaids) #too few
table(tet2_chip_repeats_single_lod$antiplatelets) #ok
table(tet2_chip_repeats_single_lod$immunosuppressants) #too few
table(tet2_chip_repeats_single_lod$metformin) #ok
table(tet2_chip_repeats_single_lod$RAAS) #ok
table(tet2_chip_repeats_single_lod$ccb) #ok
table(tet2_chip_repeats_single_lod$bb) #ok
table(tet2_chip_repeats_single_lod$thiazides) #ok
table(tet2_chip_repeats_single_lod$sulphonylureas) #too few
table(tet2_chip_repeats_single_lod$lld) #ok
table(tet2_chip_repeats_single_lod$ssri) #too few
table(tet2_chip_repeats_single_lod$loop) #too few
table(tet2_chip_repeats_single_lod$ppi) #ok

#replace age -7 (indicating >= 90) with 90
tet2_chip_repeats_single_lod$indager <- replace(tet2_chip_repeats_single_lod$indager, tet2_chip_repeats_single_lod$indager == -7, 90)

#add m or f sex
tet2_chip_repeats_single_lod$sex <- ifelse(tet2_chip_repeats_single_lod$indsex == 1, "male", "female")
##############################################
#use olsrr to step forward p
#var select with cholesterol
var_select <- glm(growth_rate ~ indager + indsex + diabetes + heart_dis + osteoporosis + asthma + htn + antiplatelets + metformin + RAAS + ccb + bb + thiazides + lld + ppi + hgb + wbc + mch + hscrp + rtin + BMIVAL + chol + count + cfib, data = tet2_chip_repeats_single_lod)

k <- ols_step_forward_p(var_select, p_val = 0.35)
plot(k)
k$metrics
k$model

#include chol in base model
tet2_gr_model <- lm(growth_rate ~ indager + hscrp + lld + bb + hgb + antiplatelets + thiazides + sex + htn + chol, data = tet2_chip_repeats_single_lod)
summary(tet2_gr_model)
confint(tet2_gr_model, level = 0.95)


#look for autocorrelation
dw_result <- dwtest(tet2_gr_model)
print(dw_result)
#data:  tet2_gr_model
#DW = 1.8791, p-value = 0.2407
#alternative hypothesis: true autocorrelation is greater than 0
#no evidence of autocorrelation in residuals

#loocv analysis
#extract RMSE from model
sqrt(mean(tet2_gr_model$residuals^2))
# 5.854274

#mae
mae(tet2_chip_repeats_single_lod$growth_rate, predict(tet2_gr_model))
#4.525417

#do LOOCV
train_control <- caret::trainControl(method = "LOOCV")

model_loocv <- caret::train(growth_rate ~ indager + hscrp + lld + bb + hgb + antiplatelets + thiazides + indsex + htn + chol, data = tet2_chip_repeats_single_lod, method = 'lm', trControl = train_control)

print(model_loocv)

#plot model
tet2_chip_repeats_single_lod$p <- predict(tet2_gr_model)
ggplot(tet2_chip_repeats_single_lod) +
  geom_point(aes(x=p, y=growth_rate)) +
  geom_abline(intercept = 0, slope = 1, color = "green")

#check vifs
gvif(tet2_gr_model)

#regression diagnostic plots
plot(tet2_gr_model)

#statin only model 
statin_only_model <- lm(growth_rate ~ lld, data = tet2_chip_repeats_single_lod)
summary(statin_only_model)
confint(statin_only_model)

#thiazide only model 
thiazide_only_model <- lm(growth_rate ~ thiazides, data = tet2_chip_repeats_single_lod)
summary(thiazide_only_model)
confint(thiazide_only_model)

#antiplatelet only model 
antiplatelet_only_model <- lm(growth_rate ~ antiplatelets, data = tet2_chip_repeats_single_lod)
summary(antiplatelet_only_model)
confint(antiplatelet_only_model)

#age 
age_only_model <- lm(growth_rate ~ indager, data = tet2_chip_repeats_single_lod)
summary(age_only_model)
confint(age_only_model)

#crp 
crp_only_model <- lm(growth_rate ~ hscrp, data = tet2_chip_repeats_single_lod)
summary(crp_only_model)
confint(crp_only_model)

#bb 
bb_only_model <- lm(growth_rate ~ bb, data = tet2_chip_repeats_single_lod)
summary(bb_only_model)
confint(bb_only_model)

#hgb 
hgb_only_model <- lm(growth_rate ~ hgb, data = tet2_chip_repeats_single_lod)
summary(hgb_only_model)
confint(hgb_only_model)

#sex
sex_only_model <- lm(growth_rate ~ sex, data = tet2_chip_repeats_single_lod)
summary(sex_only_model)
confint(sex_only_model)

#htn
htn_only_model <- lm(growth_rate ~ htn, data = tet2_chip_repeats_single_lod)
summary(htn_only_model)
confint(htn_only_model)

chol_only_model <- lm(growth_rate ~ chol, data = tet2_chip_repeats_single_lod)
summary(chol_only_model)
confint(chol_only_model)

#now do boot strap analysis##############################################
#sd of vaf measurements calculated previously at 2.7%

#AF_2
AF_2 <- subset(tet2_chip_repeats_single_lod, select = c(idauniq, AF_2))

for (i in (1:1000)) {
  a <- rep(NA, nrow(AF_2))
  AF_2 <- cbind(AF_2, a)
  names(AF_2)[i+2] <- paste0("af_2_", i)
  rm(a)
}

set.seed(214)

for (j in (1:nrow(AF_2))){
  a <- rtruncnorm(1000, a=0.001, b=55, mean = AF_2$AF_2[j], sd = 0.027)
  AF_2[j,3:1002] <- a
  rm(a)
}


#AF_4
AF_4 <- subset(tet2_chip_repeats_single_lod, select = c(idauniq, AF_4))

for (i in (1:1000)) {
  a <- rep(NA, nrow(AF_4))
  AF_4 <- cbind(AF_4, a)
  names(AF_4)[i+2] <- paste0("af_4_", i)
  rm(a)
}


for (j in (1:nrow(AF_4))){
  a <- rtruncnorm(1000, a=0.001, b=55, mean = AF_4$AF_4[j], sd = 0.027)
  AF_4[j,3:1002] <- a
  rm(a)
}

#AF_6
AF_6 <- subset(tet2_chip_repeats_single_lod, select = c(idauniq, AF_6))

for (i in (1:1000)) {
  a <- rep(NA, nrow(AF_6))
  AF_6 <- cbind(AF_6, a)
  names(AF_6)[i+2] <- paste0("af_6_", i)
  rm(a)
}


for (j in (1:nrow(AF_6))){
  a <- rtruncnorm(1000, a=0.001, b=55, mean=AF_6$AF_6[j], sd = 0.027)
  AF_6[j,3:1002] <- a
  rm(a)
}


#AF_8
AF_8 <- subset(tet2_chip_repeats_single_lod, select = c(idauniq, AF_8))

for (i in (1:1000)) {
  a <- rep(NA, nrow(AF_8))
  AF_8 <- cbind(AF_8, a)
  names(AF_8)[i+2] <- paste0("af_8_", i)
  rm(a)
}


for (j in (1:nrow(AF_8))){
  a <- rtruncnorm(1000, a=0.001, b=55, mean = AF_8$AF_8[j], sd = 0.027)
  AF_8[j,3:1002] <- a
  rm(a)
}


#AF_9
AF_9 <- subset(tet2_chip_repeats_single_lod, select = c(idauniq, AF_9))

for (i in (1:1000)) {
  a <- rep(NA, nrow(AF_9))
  AF_9 <- cbind(AF_9, a)
  names(AF_9)[i+2] <- paste0("af_9_", i)
  rm(a)
}


for (j in (1:nrow(AF_9))){
  a <- rtruncnorm(1000, a=0.001, b=55, mean = AF_9$AF_9[j], sd = 0.027)
  AF_9[j,3:1002] <- a
  rm(a)
}


#calculate i_1 etc, each in different df
#i_1
tet2_roc_boot_i1 <- subset(tet2_chip_repeats_single_lod, select = c(idauniq))

for (i in (1:1000)) {
  a <- rep(NA, nrow(tet2_roc_boot_i1))
  tet2_roc_boot_i1 <- cbind(tet2_roc_boot_i1, a)
  names(tet2_roc_boot_i1)[i+1] <- paste0("roc_", i)
  rm(a)
}

for (j in (1:nrow(tet2_roc_boot_i1))){
  for (i in (1:1000)){
    tet2_roc_boot_i1[j, (i+1)] <- ((AF_9[j, (i+2)] / AF_2[j, (i+2)])^{1/15} - 1)*100
  }
}

#i_2
tet2_roc_boot_i2 <- subset(tet2_chip_repeats_single_lod, select = c(idauniq))

for (i in (1:1000)) {
  a <- rep(NA, nrow(tet2_roc_boot_i2))
  tet2_roc_boot_i2 <- cbind(tet2_roc_boot_i2, a)
  names(tet2_roc_boot_i2)[i+1] <- paste0("roc_", i)
  rm(a)
}

for (j in (1:nrow(tet2_roc_boot_i2))){
  for (i in (1:1000)){
    tet2_roc_boot_i2[j, (i+1)] <- ((AF_8[j, (i+2)] / AF_2[j, (i+2)])^{1/13} - 1)*100
  }
}


#i_3
tet2_roc_boot_i3 <- subset(tet2_chip_repeats_single_lod, select = c(idauniq))

for (i in (1:1000)) {
  a <- rep(NA, nrow(tet2_roc_boot_i3))
  tet2_roc_boot_i3 <- cbind(tet2_roc_boot_i3, a)
  names(tet2_roc_boot_i3)[i+1] <- paste0("roc_", i)
  rm(a)
}

for (j in (1:nrow(tet2_roc_boot_i3))){
  for (i in (1:1000)){
    tet2_roc_boot_i3[j, (i+1)] <- ((AF_9[j, (i+2)] / AF_4[j, (i+2)])^{1/11} - 1)*100
  }
}


#i_4
tet2_roc_boot_i4 <- subset(tet2_chip_repeats_single_lod, select = c(idauniq))

for (i in (1:1000)) {
  a <- rep(NA, nrow(tet2_roc_boot_i4))
  tet2_roc_boot_i4 <- cbind(tet2_roc_boot_i4, a)
  names(tet2_roc_boot_i4)[i+1] <- paste0("roc_", i)
  rm(a)
}

for (j in (1:nrow(tet2_roc_boot_i4))){
  for (i in (1:1000)){
    tet2_roc_boot_i4[j, (i+1)] <- ((AF_8[j, (i+2)] / AF_4[j, (i+2)])^{1/9} - 1)*100
  }
}


#i_5
tet2_roc_boot_i5 <- subset(tet2_chip_repeats_single_lod, select = c(idauniq))

for (i in (1:1000)) {
  a <- rep(NA, nrow(tet2_roc_boot_i5))
  tet2_roc_boot_i5 <- cbind(tet2_roc_boot_i5, a)
  names(tet2_roc_boot_i5)[i+1] <- paste0("roc_", i)
  rm(a)
}

for (j in (1:nrow(tet2_roc_boot_i5))){
  for (i in (1:1000)){
    tet2_roc_boot_i5[j, (i+1)] <- ((AF_8[j, (i+2)] / AF_6[j, (i+2)])^{1/5} - 1)*100
  }
}


#i_6
tet2_roc_boot_i6 <- subset(tet2_chip_repeats_single_lod, select = c(idauniq))

for (i in (1:1000)) {
  a <- rep(NA, nrow(tet2_roc_boot_i6))
  tet2_roc_boot_i6 <- cbind(tet2_roc_boot_i6, a)
  names(tet2_roc_boot_i6)[i+1] <- paste0("roc_", i)
  rm(a)
}

for (j in (1:nrow(tet2_roc_boot_i6))){
  for (i in (1:1000)){
    tet2_roc_boot_i6[j, (i+1)] <- ((AF_6[j, (i+2)] / AF_2[j, (i+2)])^{1/8} - 1)*100
  }
}


##########################

#merge the tet2_roc_boot dfs together

tet2_roc_boot <- rows_patch(tet2_roc_boot_i1, tet2_roc_boot_i2, by="idauniq")
tet2_roc_boot <- rows_patch(tet2_roc_boot, tet2_roc_boot_i3, by="idauniq")
tet2_roc_boot <- rows_patch(tet2_roc_boot, tet2_roc_boot_i4, by="idauniq")
tet2_roc_boot <- rows_patch(tet2_roc_boot, tet2_roc_boot_i5, by="idauniq")
tet2_roc_boot <- rows_patch(tet2_roc_boot, tet2_roc_boot_i6, by="idauniq")


#now merge with tet2_chip_repeats_single_lod

tet2_single_roc_boot <- merge(tet2_chip_repeats_single_lod, tet2_roc_boot, by = "idauniq")


#separating out coefficients
simulations <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(simulations) <- c("(intercept)", "indager", "hscrp", "lldTRUE", "bbTRUE", "hgb", "antiplateletsTRUE", "thiazidesTRUE", "sexmale", "htnTRUE", "chol")


for (i in (1:1000)) {
  temp <- subset(tet2_single_roc_boot, select = c("indager", "hscrp", "lld", "bb", "hgb", "antiplatelets", "thiazides", "sex", "htn", "chol", paste0("roc_",i)))
  names(temp)[11] <- "run"
  model <- lm(formula = run ~ indager + hscrp + lld + bb + hgb + antiplatelets + thiazides + sex + htn + chol, data = temp)
  simulations[i,] <- coefficients(model) 
  rm(model)
  rm(temp)
}


#generate p value from bootstrap
#this tests the hypothesis that the coefficients are equal to 0
simulations2 <- transpose(simulations)
rownames(simulations2) <- c("(intercept)", "indager", "hscrp", "lldTRUE", "bbTRUE", "hgb", "antiplateletsTRUE", "thiazidesTRUE", "sexmale", "htnTRUE", "chol")

simulations2m <- as.matrix(simulations2)

pvals <- sapply(1:nrow(simulations2m), function(x) {
  distribution <- ecdf(simulations2m[x,])
  qt0 <- distribution(0)
  if(qt0 < 0.5){
    return(2*qt0)
  } else {
    return(2*(1-qt0))
  }
})

boot_strap_p_values <- data.frame(rownames(simulations2), pvals)



#for lld confidence interval
lld_ordered <- sort(simulations$lldTRUE)
lld_ordered[25] 
lld_ordered[975] 
lld_ordered[500]

#indager
age_ordered <- sort(simulations$indager)
age_ordered[25] 
age_ordered[975] 
age_ordered[500] 

#hscrp
crp_ordered <- sort(simulations$hscrp)
crp_ordered[25] 
crp_ordered[975] 
crp_ordered[500] 

#bbTRUE
bb_ordered <- sort(simulations$bbTRUE)
bb_ordered[25] 
bb_ordered[975] 
bb_ordered[500] 

#hgb
hgb_ordered <- sort(simulations$hgb)
hgb_ordered[25] 
hgb_ordered[975] 
hgb_ordered[500] 

#antiplateletsTRUE
antiplatelets_ordered <- sort(simulations$antiplateletsTRUE)
antiplatelets_ordered[25] 
antiplatelets_ordered[975] 
antiplatelets_ordered[500] 

#thiazidesTRUE
thi_ordered <- sort(simulations$thiazidesTRUE)
thi_ordered[25] 
thi_ordered[975] 
thi_ordered[500] 

#indsex
sex_ordered <- sort(simulations$sexmale)
sex_ordered[25] 
sex_ordered[975] 
sex_ordered[500] 

#htnTRUE
htn_ordered <- sort(simulations$htnTRUE)
htn_ordered[25] 
htn_ordered[975] 
htn_ordered[500] 

#chol
chol_ordered <- sort(simulations$chol)
chol_ordered[25] 
chol_ordered[975] 
chol_ordered[500] 

