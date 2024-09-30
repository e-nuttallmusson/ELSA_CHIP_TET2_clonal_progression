#script to look at tet2 high vaf vs low vaf logistic regression

#import wave_8_9 chip mutation data, wave_8_9_comorbidities, drugs and blood parameters

setwd("N:/My Documents/ELSA_data/upgrade/")
chip_cases_w8_9_metadata <- read_xlsx("wave_8_and_9_chip_cases_comorbs_drugs_bloods.xlsx")

#add recent_vaf variable
chip_cases_w8_9_metadata$recent_vaf <- chip_cases_w8_9_metadata$AF_9
chip_cases_w8_9_metadata$recent_vaf <- ifelse(is.na(chip_cases_w8_9_metadata$recent_vaf), chip_cases_w8_9_metadata$AF_8, chip_cases_w8_9_metadata$recent_vaf)

#add mutation count variable
chip_cases_w8_9_metadata <- chip_cases_w8_9_metadata %>% group_by(idauniq) %>% mutate(count = n())

#just consider TET2 chip - consider largest vaf if more than one tet2 mutation present
chip_tet2_mut <- chip_cases_w8_9_metadata %>% dplyr::filter(grepl("TET2", Gene.refGene))

chip_tet2_mut_dups <- chip_tet2_mut %>% group_by(idauniq) %>% filter(n() > 1) %>% arrange(idauniq, desc(recent_vaf))

chip_tet2_mut_dups2 <- chip_tet2_mut_dups[!duplicated(chip_tet2_mut_dups$idauniq),] 

#get non duplicated data
chip_tet2_mut_no_dups <- chip_tet2_mut %>% group_by(idauniq) %>% filter(n() == 1) %>% ungroup()

#recombine

chip_tet2_mut_no_dups <- rbind(chip_tet2_mut_no_dups, chip_tet2_mut_dups2)

#add variable for variant type - nonsense, fs etc
chip_tet2_mut_no_dups$variant_type <- NA
chip_tet2_mut_no_dups$variant_type[grep("fs*", chip_tet2_mut_no_dups$NonsynOI)] <- "nonsense"
chip_tet2_mut_no_dups$variant_type[grep("X", chip_tet2_mut_no_dups$NonsynOI)] <- "nonsense"
chip_tet2_mut_no_dups$variant_type[grep("nan", chip_tet2_mut_no_dups$NonsynOI)] <- "splice_variant"
chip_tet2_mut_no_dups$variant_type[grep("del", chip_tet2_mut_no_dups$NonsynOI)] <- "indel"
chip_tet2_mut_no_dups$variant_type[grep("\\*$", chip_tet2_mut_no_dups$NonsynOI)] <- "nonsense"
chip_tet2_mut_no_dups$variant_type[is.na(chip_tet2_mut_no_dups$variant_type)] <- "snv"

#replace -7 (indicating age > 90) in indager
chip_tet2_mut_no_dups$indager <- replace(chip_tet2_mut_no_dups$indager, which(chip_tet2_mut_no_dups$indager == -7), 90)
chip_tet2_mut_no_dups$mapval <- replace(chip_tet2_mut_no_dups$mapval, which(chip_tet2_mut_no_dups$mapval < 0), NA)

#convert boolean to 1/0
cols <- sapply(chip_tet2_mut_no_dups, is.logical)
chip_tet2_mut_no_dups[,cols]<- lapply(chip_tet2_mut_no_dups[,cols], as.numeric)

#drop wave AF vals as working from recent_vaf
chip_tet2_mut_no_dups <- subset(chip_tet2_mut_no_dups, select = -c(AF_2, AF_4, AF_6, AF_8, AF_9))

#replace nas with median
#for alcohol - none is 0, occasional is 1, moderate is 2, heavy is 3
chip_tet2_mut_no_dups$alcohol <- ifelse(chip_tet2_mut_no_dups$alcohol == "none", 0, chip_tet2_mut_no_dups$alcohol)
chip_tet2_mut_no_dups$alcohol <- ifelse(chip_tet2_mut_no_dups$alcohol == "occasional", 1, chip_tet2_mut_no_dups$alcohol)
chip_tet2_mut_no_dups$alcohol <- ifelse(chip_tet2_mut_no_dups$alcohol == "moderate", 2, chip_tet2_mut_no_dups$alcohol)
chip_tet2_mut_no_dups$alcohol <- ifelse(chip_tet2_mut_no_dups$alcohol == "heavy", 3, chip_tet2_mut_no_dups$alcohol)

chip_tet2_mut_no_dups$alcohol <- as.numeric(chip_tet2_mut_no_dups$alcohol)

imps <- lapply(chip_tet2_mut_no_dups, median, na.rm = TRUE)

for (i in colnames(chip_tet2_mut_no_dups)) {
  chip_tet2_mut_no_dups[,i][is.na(chip_tet2_mut_no_dups[,i])] <- as.numeric(imps[i])
}

#drop those with blood disorder and cancer/cancer treatment
chip_tet2_mut_no_dups <- chip_tet2_mut_no_dups %>% dplyr::filter(blood_dis == 0)

#there are 3 samples with designated proxy VAF 0.01 who were cases samples that became controls in w9 - remove 
summary(chip_tet2_mut_no_dups$recent_vaf)
chip_tet2_mut_no_dups <- chip_tet2_mut_no_dups %>% dplyr::filter(recent_vaf != 0.01)

##################logistic regression
#binary outcome - low vaf high vaf
chip_tet2_mut_no_dups$vaf_cat <- ifelse(chip_tet2_mut_no_dups$recent_vaf <0.1, 0, 1)

#split
set.seed(245)
split_tet2_mut_cat <- sample(c(rep(0, 0.7 *nrow(chip_tet2_mut_no_dups)), rep(1, 0.3*nrow(chip_tet2_mut_no_dups))))
split_tet2_mut_cat <- append(split_tet2_mut_cat, 1)

train_tet2_mut_cat <- chip_tet2_mut_no_dups[split_tet2_mut_cat == 0, ]
test_tet2_mut_cat <- chip_tet2_mut_no_dups[split_tet2_mut_cat == 1, ]

#look at numbers in training cohort
table(train_tet2_mut_cat$vaf_cat)

#0   1 
#328  60 

table(test_tet2_mut_cat$vaf_cat)
# 0   1 
#131  36  

#look for categorical covariates with sufficient numbers to include in analysis

table(train_tet2_mut_cat$indsex, train_tet2_mut_cat$vaf_cat)
table(train_tet2_mut_cat$diabetes, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$ihd, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$heart_dis, train_tet2_mut_cat$vaf_cat) #ok
table(train_tet2_mut_cat$stroke, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$osteoporosis, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$asthma, train_tet2_mut_cat$vaf_cat) #ok
table(train_tet2_mut_cat$lung_disease, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$current_smoker, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$htn, train_tet2_mut_cat$vaf_cat) #ok
table(train_tet2_mut_cat$alcohol, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$dementia, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$inflam_dis, train_tet2_mut_cat$vaf_cat) #too few?
table(train_tet2_mut_cat$statins, train_tet2_mut_cat$vaf_cat) #ok
table(train_tet2_mut_cat$ccb, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$acei, train_tet2_mut_cat$vaf_cat) #ok
table(train_tet2_mut_cat$ppi, train_tet2_mut_cat$vaf_cat) #ok
table(train_tet2_mut_cat$antiplatelets, train_tet2_mut_cat$vaf_cat) #ok
table(train_tet2_mut_cat$bb, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$arb, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$anticoag, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$ssri, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$vit_d_drug, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$metformin, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$adr_ag, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$tca, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$non_op_an, train_tet2_mut_cat$vaf_cat) #ok
table(train_tet2_mut_cat$inh_cort, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$nsaids, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$bisphos, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$gout_rx, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$ur_ret, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$glauc, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$opioids, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$epilep, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$oestr, train_tet2_mut_cat$vaf_cat) #too few
table(train_tet2_mut_cat$loop, train_tet2_mut_cat$vaf_cat) #too few


#step wise variable selection
library(olsrr)
forward_selection <- glm(vaf_cat ~ heart_dis + asthma + htn + statins + ccb + acei + ppi + antiplatelets + wbc + hgb + hscrp + indager + indsex + chol + ldl + trig + rtin + cfib + sysval + diaval + pulval + mch + count + variant_type, data = train_tet2_mut_cat, family = binomial(link = 'logit'))
k <- ols_step_forward_p(forward_selection, p_val = 0.35)
plot(k)
k$metrics
k$model

#take top 5 predictors given data structure - don't want to overfit
#indager, asthma, statins, hgb, chol
#test for linearity between continuous vars and logits - none significant so assumption not violated
library(car)
boxTidwell(formula = vaf_cat ~ indager, other.x = ~ wbc + antiplatelets + statins + hscrp, data = train_tet2_mut_cat)
boxTidwell(formula = vaf_cat ~ hgb, other.x = ~ indager + antiplatelets + statins + hscrp, data = train_tet2_mut_cat)
boxTidwell(formula = vaf_cat ~ chol, other.x = ~ wbc + antiplatelets + statins + indager + hscrp, data = train_tet2_mut_cat)


forward_selection2 <- glm(vaf_cat ~ indager + asthma + statins + hgb + chol, family = binomial(link = "logit"), data = train_tet2_mut_cat)
summary(forward_selection2)
pR2(forward_selection2)

exp(cbind(OR=coef(forward_selection2), confint(forward_selection2)))

#check for multicolinearity
library(glmtoolbox)
gvif(forward_selection2)

#check for linearity in logits
#select just numeric vars
model_num <- subset(train_tet2_mut_cat, select = c(indager, hgb, chol))
predictors <- colnames(model_num)
model_num$probabilities <- forward_selection2$fitted.values
library(tidyr)
model_num <- model_num %>% mutate(logit = log(probabilities/(1-probabilities))) %>% select(-probabilities) %>% gather(key = "predictors", value = "predictor.value", -logit)

library(ggplot2)

ggplot(model_num, aes(y = logit, x = predictor.value))+
  geom_point(size = 0.5, alpha = 0.5)+
  geom_smooth(method = "loess") +
  theme_bw() +
  facet_wrap(~predictors, scale = "free_x")

#linearity assumptions OK

#assess accuracy of model

fitted.results <- predict(forward_selection2, newdata = test_tet2_mut_cat, type = 'response')
fitted.results <- as.numeric(ifelse(fitted.results> 0.5, 1, 0))

misClasError <- mean(fitted.results!= test_tet2_mut_cat$vaf_cat)
print(paste('Accuraacy is', 1-misClasError))

#Accuracy is 0.790419161676647


#rocr

library(ROCR)

p <- predict(forward_selection2, newdata = test_tet2_mut_cat, type = 'response')
pr <- prediction(p, test_tet2_mut_cat$vaf_cat)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc
#0.7232824

#check age adjusted univariate logistic regressions

library(oddsratio)
#asthma 
asthma_age <- glm(vaf_cat ~ indager + asthma, family = binomial(link = "logit"), data = train_tet2_mut_cat)
summary(asthma_age)
#age adjusted odds ratio
exp(cbind(OR=coef(asthma_age), confint(asthma_age)))

#hgb
hgb_age <- glm(vaf_cat ~ indager + hgb, family = binomial(link = "logit"), data = train_tet2_mut_cat)
summary(hgb_age)
#age adjusted odds ratio
exp(cbind(OR=coef(hgb_age), confint(hgb_age)))


#statins 
statins_age <- glm(vaf_cat ~ indager + statins, family = binomial(link = "logit"), data = train_tet2_mut_cat)
summary(statins_age)
#age adjusted odds ratio
exp(cbind(OR=coef(statins_age), confint(statins_age)))

#chol 
chol_age <- glm(vaf_cat ~ indager + chol, family = binomial(link = "logit"), data = train_tet2_mut_cat)
summary(chol_age)
#age adjusted odds ratio
exp(cbind(OR=coef(chol_age), confint(chol_age)))

