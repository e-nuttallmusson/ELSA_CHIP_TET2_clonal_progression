library(dplyr)
library(tidyr)
library(splitstackshape)

#import combined IFC calls post whitelist filter
setwd("~/Documents/UCL/PhD/ELSA_variant_calling/")
ifc_no <- 67
          
###################
setwd(paste0('~/Documents/UCL/PhD/ELSA_variant_calling/mutect2/IFC_',ifc_no))

ifc_original <- read.csv(paste0("IFC_", ifc_no, "_whitelist_combined.csv"))

#get data in workable format

ifc_pgtpid<-dplyr::filter(ifc_original, grepl("GT:AD:AF:DP:F1R2:F2R1:PGT:PID:PS:SB",Otherinfo12))
ifc_original<-dplyr::filter(ifc_original, !grepl("GT:AD:AF:DP:F1R2:F2R1:PGT:PID:PS:SB",Otherinfo12))

#make copy of Otherinfo13 to not lose info
ifc_pgtpid$Otherinfo13_2 <- ifc_pgtpid$Otherinfo13
ifc_original$Otherinfo13_2 <- ifc_original$Otherinfo13

#separate out columns for pgt pid variants
ifc_pgtpid <- separate(ifc_pgtpid, Otherinfo13, str_split_1(ifc_original$Otherinfo12[1], ":"), sep = ":", convert = TRUE, extra = "merge")
ifc_pgtpid <- cSplit(ifc_pgtpid, "SB", sep = ":")
ifc_pgtpid <- ifc_pgtpid %>% rename("PGT" = "SB_1", "PID" = "SB_2", "PS" = "SB_3", "SB" = "SB_4")
ifc_pgtpid <- subset(ifc_pgtpid, select = -c(PGT, PID, PS))
ifc_pgtpid <- cSplit(ifc_pgtpid, "SB", sep = ",")

#separate out columns for non pgt variants
ifc_original <- separate(ifc_original, Otherinfo13, str_split_1(ifc_original$Otherinfo12[1], ":"), sep = ":", convert = TRUE, extra = "merge")
ifc_original <- cSplit(ifc_original, "SB", sep = ",")

#recombine
ifc_original <- rbind(ifc_original, ifc_pgtpid)


#move multiallelic to manual review
ifc_multi<-dplyr::filter(ifc_original, grepl("multiallelic",Otherinfo10))
ifc<-dplyr::filter(ifc_original, !grepl("multiallelic",Otherinfo10))


#flag for removal if 0 alt/ref reads on one strand
ifc$no_strand <- c(rep(FALSE, nrow(ifc)))
for (i in (1:nrow(ifc))) {
  if(ifc$SB_1[i] == 0 || ifc$SB_2[i] == 0 || ifc$SB_3[i] == 0 || ifc$SB_4[i] == 0) {
    ifc$no_strand[i] <-TRUE
  }
}


#remove if alternate allele depth < 20
ifc$AD_copy <- ifc$AD
ifc <- separate(ifc, AD, c("AD1", "AD2"), sep = ",", convert = TRUE)
low_alt_depth <- c(rep("FALSE", length(ifc$AD2)))
ifc$low_alt_depth <- low_alt_depth
ifc$low_alt_depth[ifc$AD2 < 20] <- "TRUE"

#section to remove if VAF < 0.02 and move to another dataframe of filtered calls
ifc$AF <- as.numeric(ifc$AF)
ifc$AF_copy <- ifc$AF
ifc$lowvaf <- c(rep("FALSE", nrow(ifc)))
ifc$lowvaf[ifc$AF <0.02] <- "TRUE"


#remove asxl1 non multiallelic calls for specific frameshifts if "G645Vfs*58" or and "G646Wfs*12" | "G646Vfs*58" with VAF<0.1
asxl1_fail <- c(rep(FALSE, nrow(ifc)))
for (row in 1:nrow(ifc)) {
  if (ifc$NonsynOI[row] == "G645Vfs*58") {
    asxl1_fail[row] <- TRUE
  }
  if (ifc$NonsynOI[row] == "G646Wfs*12" && ifc$AF[row] <0.1) {
    asxl1_fail[row] <- TRUE
  }
  if (ifc$NonsynOI[row] == "G646Vfs*58" && ifc$AF[row] <0.1) {
    asxl1_fail[row] <- TRUE
  }
}

ifc$asxl1_fail <- asxl1_fail

#binomial test for VAF
var_binomial <- c()
for (row in 1:nrow(ifc)) {
  a <- binom.test(c(ifc$AD2[row], ifc$AD1[row]), p = 0.5, alternative = "two.sided")$p.value
  var_binomial <- append(var_binomial, a)
}
ifc$var_binomial <- var_binomial 

ifc$fail_binomial <- ifelse(ifc$var_binomial > 0.01, "TRUE", "FALSE")

#TET2 exceptions to binomial
for (row in 1:nrow(ifc)) {
  ifc$fail_binomial[row] <- ifelse(ifc$Gene.refGene[row] == "TET2" && ifc$NonsynOI[row] == "H1904R", "FALSE", "FALSE")
  ifc$fail_binomial[row] <- ifelse(ifc$Gene.refGene[row] == "TET2" && ifc$NonsynOI[row] == "I1873T", "FALSE", "FALSE")
  ifc$fail_binomial[row] <- ifelse(ifc$Gene.refGene[row] == "TET2" && ifc$NonsynOI[row] == "T1884A", "FALSE", "FALSE")
}




#section to remove variants flagged as artifacts in cv's paper
#import dataframes
setwd("~/Documents/UCL/PhD/ELSA_variant_calling/cv_chip_pipeline")
cv_artifacts <- read.xlsx("1-s2.0-S000649712300143X-BLOOD_BLD-2022-018825-mmc1.xlsx", sheet = "S2. Variants removed – artifact")
cv_fail_binomial <- read.xlsx("1-s2.0-S000649712300143X-BLOOD_BLD-2022-018825-mmc1.xlsx", sheet = "S4. Variants failing binomial")

#separate variants by ";" delimiter
#if present then remove call and add to another dataframe of filtered calls

#need to do the same for the manual filter output, and combine with filtered calls
artifact <- c()
vec <- c()
for (row in 1:nrow(ifc)) {
  x <- ifc$transcriptOI[row]
  for (i in 1:nrow(cv_artifacts)) {
    b <- identical(x, cv_artifacts$variant[i])
    vec <- append(vec, b)
  }
  c <- TRUE%in%vec 
  artifact <- append(artifact, c)
  vec <- c()
}
ifc$cv_artifact <- artifact

#binomial filter
binomial <- c()
vec2 <- c()
for (row in 1:nrow(ifc)) {
  x <- ifc$transcriptOI[row]
  for (i in 1:nrow(cv_fail_binomial)) {
    d <- identical(x, cv_fail_binomial$Varianta[i])
    vec2 <- append(vec, d)
  }
  e <- TRUE%in%vec 
  binomial <- append(binomial, e)
  vec2 <- c()
}
ifc$cv_binomial_fail <- binomial



#make df for bl manual review and for filtered out calls
blacklist_manual_review <- data.frame(matrix(ncol = 55, nrow = 0))
colnames(blacklist_manual_review) <- colnames(ifc)

blacklist_filtered_out <- data.frame(matrix(ncol = 55, nrow = 0))
colnames(blacklist_filtered_out) <- colnames(ifc)


#for remaining calls if fails because of no_read_mate, VAF, cv_artifact, failed binomial then move calls to separate filtered dataframe
#no strand
blacklist_filtered_out <- rbind(blacklist_filtered_out, ifc[ifc$no_strand == TRUE,]) 
ifc<-ifc[!(ifc$no_strand == TRUE),]


#VAF
blacklist_filtered_out <- rbind(blacklist_filtered_out, ifc[ifc$lowvaf == TRUE,]) 
ifc<-ifc[!(ifc$lowvaf == TRUE),]

#alternate read depth
blacklist_filtered_out <- rbind(blacklist_filtered_out, ifc[ifc$low_alt_depth == TRUE,]) 
ifc<-ifc[!(ifc$low_alt_depth == TRUE),]

#cv artifact
blacklist_filtered_out <- rbind(blacklist_filtered_out, ifc[ifc$cv_artifact == TRUE,]) 
ifc<-ifc[!(ifc$cv_artifact == TRUE),]

#non multiallelic asxl1 fails
blacklist_filtered_out <- rbind(blacklist_filtered_out, ifc[ifc$asxl1_fail == TRUE,]) 
ifc<-ifc[!(ifc$asxl1_fail == TRUE),]

#remove calls failing base quality filter
ifc<-dplyr::filter(ifc, !grepl("base_qual",Otherinfo10))

#binomial - add to manual review, but remove those specified in the cv paper
blacklist_manual_review <-rbind(blacklist_manual_review, ifc[ifc$fail_binomial == TRUE,])
ifc<-ifc[!(ifc$fail_binomial == TRUE),]

blacklist_filtered_out <- rbind(blacklist_filtered_out, blacklist_manual_review[blacklist_manual_review$cv_binomial_fail == TRUE,])
blacklist_manual_review <- blacklist_manual_review[!(blacklist_manual_review$cv_binomial_fail == TRUE),]


#consider multiallelic for manual review################

#split AF column by ","
ifc_multi <- ifc_multi %>% separate(AF, paste0("AF_", c(1:3)), sep = ",", convert = TRUE, extra = "merge", fill = "right")

#artifacts - flag
artifact <- c()
vec <- c()
for (row in 1:nrow(ifc_multi)) {
  x <- ifc_multi$transcriptOI[row]
  for (i in 1:nrow(cv_artifacts)) {
    b <- identical(x, cv_artifacts$variant[i])
    vec <- append(vec, b)
  }
  c <- TRUE%in%vec 
  artifact <- append(artifact, c)
  vec <- c()
}
ifc_multi$cv_artifact <- artifact


#filter reamining multiallelic calls for only those with VAF <0.1
ifc_multi$AF_1 <- as.numeric(ifc_multi$AF_1)
ifc_multi$AF_2 <- as.numeric(ifc_multi$AF_2)
ifc_multi$AF_3 <- as.numeric(ifc_multi$AF_3)

#creat vaf flag
multi_lowvaf <- c()
for (row in 1:nrow(ifc_multi)) {
  rm(x)
  if (is.na(ifc_multi$AF_3[row]) == TRUE) {
    x <- ifc_multi$AF_1[row] < 0.1 && ifc_multi$AF_2[row] < 0.1
  }
  else {
    x <- ifc_multi$AF_1[row] < 0.1 && ifc_multi$AF_2[row] < 0.1 && ifc_multi$AF_3[row] < 0.1
  }
  multi_lowvaf <- append(multi_lowvaf, x)
}
ifc_multi$multi_lowvaf <- multi_lowvaf



#create dataframe for filtered multiallelic calls

multiallelic_filtered_out <- data.frame(matrix(ncol = 48, nrow = 0))
colnames(multiallelic_filtered_out) <- colnames(ifc_multi)

#remove cv artifacts
#cv artifact
multiallelic_filtered_out <- rbind(multiallelic_filtered_out, ifc_multi[ifc_multi$cv_artifact == TRUE,]) 
ifc_multi<-ifc_multi[!(ifc_multi$cv_artifact == TRUE),]

#remove variants if all vafs <0.02
#creat vaf flag - general vaf
gen_lowvaf <- c()
for (row in 1:nrow(ifc_multi)) {
  rm(x)
  if (is.na(ifc_multi$AF_3[row]) == TRUE) {
    x <- ifc_multi$AF_1[row] < 0.02 && ifc_multi$AF_2[row] < 0.02
  }
  else {
    x <- ifc_multi$AF_1[row] < 0.02 && ifc_multi$AF_2[row] < 0.02 && ifc_multi$AF_3[row] < 0.02
  }
  gen_lowvaf <- append(gen_lowvaf, x)
}
ifc_multi$gen_lowvaf <- gen_lowvaf


#remove low vaf asxl1 variants
asxl1_fs<-subset(ifc_multi, NonsynOI == "G646Wfs*12" | NonsynOI == "G646Vfs*58")
asxl1_fs_lowvaf <- subset(asxl1_fs, multi_lowvaf == TRUE)
asxl1_fs_highvaf <- subset(asxl1_fs, multi_lowvaf == FALSE)
ifc_multi <- subset(ifc_multi, NonsynOI != "G646Wfs*12" & NonsynOI != "G646Vfs*58")

#add high vaf variants back in
ifc_multi <- rbind(ifc_multi, asxl1_fs_highvaf)

#remove low gen vaf

multiallelic_filtered_out <- rbind(multiallelic_filtered_out, ifc_multi[ifc_multi$gen_lowvaf == TRUE,], fill = TRUE) 
ifc_multi<-ifc_multi[!(ifc_multi$gen_lowvaf == TRUE),]



######################################review variants flagged for manual review
print("number of rows in blacklist manual reivew:") 
nrow(blacklist_manual_review)



ifc <- rbind(ifc, blacklist_manual_review)

######################################review variants flagged for multiallelic manual review in ifc_multi

#allvariants <- dplyr::bind_rows(ifc, ifc_multi)
#allvariants <- allvariants[order(allvariants$Sample),]


###########consider calls for manual filtration outputted by whitelist----------------------------------------------------------------

setwd(paste0('~/Documents/UCL/PhD/ELSA_variant_calling/mutect2/IFC_',ifc_no))
ifc_man_original <- read.csv(paste0("IFC_", ifc_no, "_manualreview_combined.csv"))


ifc_man_pgtpid<-dplyr::filter(ifc_man_original, grepl("GT:AD:AF:DP:F1R2:F2R1:PGT:PID:PS:SB",Otherinfo12))
ifc_man_original<-dplyr::filter(ifc_man_original, !grepl("GT:AD:AF:DP:F1R2:F2R1:PGT:PID:PS:SB",Otherinfo12))

#make copy of Otherinfo13 to not lose info
ifc_man_pgtpid$Otherinfo13_2 <- ifc_man_pgtpid$Otherinfo13
ifc_man_original$Otherinfo13_2 <- ifc_man_original$Otherinfo13

#separate out columns for pgt pid variants
ifc_man_pgtpid <- separate(ifc_man_pgtpid, Otherinfo13, str_split_1(ifc_man_original$Otherinfo12[1], ":"), sep = ":", convert = TRUE, extra = "merge")
ifc_man_pgtpid <- cSplit(ifc_man_pgtpid, "SB", sep = ":")
ifc_man_pgtpid <- ifc_man_pgtpid %>% rename("PGT" = "SB_1", "PID" = "SB_2", "PS" = "SB_3", "SB" = "SB_4")
ifc_man_pgtpid <- subset(ifc_man_pgtpid, select = -c(PGT, PID, PS))
ifc_man_pgtpid <- cSplit(ifc_man_pgtpid, "SB", sep = ",")

#separate out columns for non pgt variants
ifc_man_original <- separate(ifc_man_original, Otherinfo13, str_split_1(ifc_man_original$Otherinfo12[1], ":"), sep = ":", convert = TRUE, extra = "merge")
ifc_man_original <- cSplit(ifc_man_original, "SB", sep = ",")

#recombine
ifc_man_original <- rbind(ifc_man_original, ifc_man_pgtpid)


#move multiallelic to manual review
ifc_man_multi<-dplyr::filter(ifc_man_original, grepl("multiallelic",Otherinfo10))
ifc_man<-dplyr::filter(ifc_man_original, !grepl("multiallelic",Otherinfo10))


#flag for removal if 0 alt/ref reads on one strand
ifc_man$no_strand <- c(rep(FALSE, nrow(ifc_man)))
for (i in (1:nrow(ifc_man))) {
  if(ifc_man$SB_1[i] == 0 || ifc_man$SB_2[i] == 0 || ifc_man$SB_3[i] == 0 || ifc_man$SB_4[i] == 0) {
    ifc_man$no_strand[i] <-TRUE
  }
}


#remove if alternate allele depth < 20
ifc_man$AD_copy <- ifc_man$AD
ifc_man <- separate(ifc_man, AD, c("AD1", "AD2"), sep = ",", convert = TRUE)
low_alt_depth <- c(rep("FALSE", length(ifc_man$AD2)))
ifc_man$low_alt_depth <- low_alt_depth
ifc_man$low_alt_depth[ifc_man$AD2 < 20] <- "TRUE"

#section to remove if VAF < 0.02 and move to another dataframe of filtered calls
ifc_man$AF <- as.numeric(ifc_man$AF)
ifc_man$AF_copy <- ifc_man$AF
ifc_man$lowvaf <- c(rep("FALSE", nrow(ifc_man)))
ifc_man$lowvaf[ifc_man$AF <0.02] <- "TRUE"


#remove asxl1 non multiallelic calls for specific frameshifts if "G645Vfs*58" or and "G646Wfs*12" | "G646Vfs*58" with VAF<0.1
asxl1_fail <- c(rep(FALSE, nrow(ifc_man)))
for (row in 1:nrow(ifc_man)) {
  if (ifc_man$NonsynOI[row] == "G645Vfs*58") {
    asxl1_fail[row] <- TRUE
  }
  if (ifc_man$NonsynOI[row] == "G646Wfs*12" && ifc_man$AF[row] <0.1) {
    asxl1_fail[row] <- TRUE
  }
  if (ifc_man$NonsynOI[row] == "G646Vfs*58" && ifc_man$AF[row] <0.1) {
    asxl1_fail[row] <- TRUE
  }
}

ifc_man$asxl1_fail <- asxl1_fail

#binomial test for VAF
var_binomial <- c()
for (row in 1:nrow(ifc_man)) {
  a <- binom.test(c(ifc_man$AD2[row], ifc_man$AD1[row]), p = 0.5, alternative = "two.sided")$p.value
  var_binomial <- append(var_binomial, a)
}
ifc_man$var_binomial <- var_binomial 

ifc_man$fail_binomial <- ifelse(ifc_man$var_binomial > 0.01, "TRUE", "FALSE")

#TET2 exceptions to binomial
for (row in 1:nrow(ifc_man)) {
  ifc_man$fail_binomial[row] <- ifelse(ifc_man$Gene.refGene[row] == "TET2" && ifc_man$NonsynOI[row] == "H1904R", "FALSE", "FALSE")
  ifc_man$fail_binomial[row] <- ifelse(ifc_man$Gene.refGene[row] == "TET2" && ifc_man$NonsynOI[row] == "I1873T", "FALSE", "FALSE")
  ifc_man$fail_binomial[row] <- ifelse(ifc_man$Gene.refGene[row] == "TET2" && ifc_man$NonsynOI[row] == "T1884A", "FALSE", "FALSE")
}




#section to remove variants flagged as artifacts in cv's paper
#import dataframes
#setwd("~/Documents/UCL/PhD/ELSA_variant_calling/cv_chip_pipeline")
#cv_artifacts <- read.xlsx("1-s2.0-S000649712300143X-BLOOD_BLD-2022-018825-mmc1.xlsx", sheet = "S2. Variants removed – artifact")
#cv_fail_binomial <- read.xlsx("1-s2.0-S000649712300143X-BLOOD_BLD-2022-018825-mmc1.xlsx", sheet = "S4. Variants failing binomial")

#separate variants by ";" delimiter
#if present then remove call and add to another dataframe of filtered calls

#need to do the same for the manual filter output, and combine with filtered calls
artifact <- c()
vec <- c()
for (row in 1:nrow(ifc_man)) {
  x <- ifc_man$transcriptOI[row]
  for (i in 1:nrow(cv_artifacts)) {
    b <- identical(x, cv_artifacts$variant[i])
    vec <- append(vec, b)
  }
  c <- TRUE%in%vec 
  artifact <- append(artifact, c)
  vec <- c()
}
ifc_man$cv_artifact <- artifact

#binomial filter
binomial <- c()
vec2 <- c()
for (row in 1:nrow(ifc_man)) {
  x <- ifc_man$transcriptOI[row]
  for (i in 1:nrow(cv_fail_binomial)) {
    d <- identical(x, cv_fail_binomial$Varianta[i])
    vec2 <- append(vec, d)
  }
  e <- TRUE%in%vec 
  binomial <- append(binomial, e)
  vec2 <- c()
}
ifc_man$cv_binomial_fail <- binomial



#make df for bl manual review and for filtered out calls
bl_man_manual_review <- data.frame(matrix(ncol = 55, nrow = 0))
colnames(bl_man_manual_review) <- colnames(ifc_man)

bl_man_filtered_out <- data.frame(matrix(ncol = 55, nrow = 0))
colnames(bl_man_filtered_out) <- colnames(ifc_man)


#for remaining calls if fails because of no_read_mate, VAF, cv_artifact, failed binomial then move calls to separate filtered dataframe
#no strand
bl_man_filtered_out <- rbind(bl_man_filtered_out, ifc_man[ifc_man$no_strand == TRUE,]) 
ifc_man<-ifc_man[!(ifc_man$no_strand == TRUE),]

#bl_man_filtered_out <- rbind(bl_man_filtered_out, ifc_man[ifc_man$orientation_bias == TRUE,]) 
#ifc_man<-ifc_man[!(ifc_man$orientation_bias == TRUE),]

#VAF
bl_man_filtered_out <- rbind(bl_man_filtered_out, ifc_man[ifc_man$lowvaf == TRUE,]) 
ifc_man<-ifc_man[!(ifc_man$lowvaf == TRUE),]

#alternate read depth
bl_man_filtered_out <- rbind(bl_man_filtered_out, ifc_man[ifc_man$low_alt_depth == TRUE,]) 
ifc_man<-ifc_man[!(ifc_man$low_alt_depth == TRUE),]

#cv artifact
bl_man_filtered_out <- rbind(bl_man_filtered_out, ifc_man[ifc_man$cv_artifact == TRUE,]) 
ifc_man<-ifc_man[!(ifc_man$cv_artifact == TRUE),]

#non multiallelic asxl1 fails
bl_man_filtered_out <- rbind(bl_man_filtered_out, ifc_man[ifc_man$asxl1_fail == TRUE,]) 
ifc_man<-ifc_man[!(ifc_man$asxl1_fail == TRUE),]

#remove calls failing base quality filter
ifc_man<-dplyr::filter(ifc_man, !grepl("base_qual",Otherinfo10))

#binomial - add to manual review, but remove those specified in the cv paper
bl_man_manual_review <-rbind(bl_man_manual_review, ifc_man[ifc_man$fail_binomial == TRUE,])
ifc_man<-ifc_man[!(ifc_man$fail_binomial == TRUE),]

bl_man_filtered_out <- rbind(bl_man_filtered_out, bl_man_manual_review[bl_man_manual_review$cv_binomial_fail == TRUE,])
bl_man_manual_review <- bl_man_manual_review[!(bl_man_manual_review$cv_binomial_fail == TRUE),]


#consider multiallelic for manual review################

#split AF column by ","
ifc_man_multi <- ifc_man_multi %>% separate(AF, paste0("AF_", c(1:3)), sep = ",", convert = TRUE, extra = "merge", fill = "right")

#artifacts - flag
artifact <- c()
vec <- c()
for (row in 1:nrow(ifc_man_multi)) {
  x <- ifc_man_multi$transcriptOI[row]
  for (i in 1:nrow(cv_artifacts)) {
    b <- identical(x, cv_artifacts$variant[i])
    vec <- append(vec, b)
  }
  c <- TRUE%in%vec 
  artifact <- append(artifact, c)
  vec <- c()
}
ifc_man_multi$cv_artifact <- artifact


#filter reamining multiallelic calls for only those with VAF <0.1
ifc_man_multi$AF_1 <- as.numeric(ifc_man_multi$AF_1)
ifc_man_multi$AF_2 <- as.numeric(ifc_man_multi$AF_2)
ifc_man_multi$AF_3 <- as.numeric(ifc_man_multi$AF_3)

#creat vaf flag
multi_lowvaf <- c()
for (row in 1:nrow(ifc_man_multi)) {
  rm(x)
  if (is.na(ifc_man_multi$AF_3[row]) == TRUE) {
    x <- ifc_man_multi$AF_1[row] < 0.1 && ifc_man_multi$AF_2[row] < 0.1
  }
  else {
    x <- ifc_man_multi$AF_1[row] < 0.1 && ifc_man_multi$AF_2[row] < 0.1 && ifc_man_multi$AF_3[row] < 0.1
  }
  multi_lowvaf <- append(multi_lowvaf, x)
}
ifc_man_multi$multi_lowvaf <- multi_lowvaf



#create dataframe for filtered multiallelic calls

multiallelic_man_filtered_out <- data.frame(matrix(ncol = 48, nrow = 0))
colnames(multiallelic_man_filtered_out) <- colnames(ifc_man_multi)

#remove cv artifacts
#cv artifact
multiallelic_man_filtered_out <- rbind(multiallelic_man_filtered_out, ifc_man_multi[ifc_man_multi$cv_artifact == TRUE,]) 
ifc_man_multi<-ifc_man_multi[!(ifc_man_multi$cv_artifact == TRUE),]

#remove variants if all vafs <0.02
#creat vaf flag - general vaf
gen_lowvaf <- c()
for (row in 1:nrow(ifc_man_multi)) {
  rm(x)
  if (is.na(ifc_man_multi$AF_3[row]) == TRUE) {
    x <- ifc_man_multi$AF_1[row] < 0.02 && ifc_man_multi$AF_2[row] < 0.02
  }
  else {
    x <- ifc_man_multi$AF_1[row] < 0.02 && ifc_man_multi$AF_2[row] < 0.02 && ifc_man_multi$AF_3[row] < 0.02
  }
  gen_lowvaf <- append(gen_lowvaf, x)
}
ifc_man_multi$gen_lowvaf <- gen_lowvaf


#remove low vaf asxl1 variants
asxl1_man_fs<-subset(ifc_man_multi, NonsynOI == "G646Wfs*12" | NonsynOI == "G646Vfs*58")
asxl1_man_fs_lowvaf <- subset(asxl1_man_fs, multi_lowvaf == TRUE)
asxl1_man_fs_highvaf <- subset(asxl1_man_fs, multi_lowvaf == FALSE)
ifc_man_multi <- subset(ifc_man_multi, NonsynOI != "G646Wfs*12" & NonsynOI != "G646Vfs*58")

#add high vaf variants back in
ifc_man_multi <- rbind(ifc_man_multi, asxl1_man_fs_highvaf)

#remove low gen vaf

multiallelic_man_filtered_out <- rbind(multiallelic_man_filtered_out, ifc_man_multi[ifc_man_multi$gen_lowvaf == TRUE,], fill = TRUE) 
ifc_man_multi<-ifc_man_multi[!(ifc_man_multi$gen_lowvaf == TRUE),]

######################################review manual review variants flagged for further manual review
print("number of rows in bl_man manual reivew:") 
nrow(bl_man_manual_review)


#REVIEW POINT - if appropriate do this step
ifc_man <- rbind(ifc_man, bl_man_manual_review)

#REVIEW POINT - combine unfiltered monoalleleic and multiallelic calls after manual filtering
man_review_variants <- dplyr::bind_rows(ifc_man, ifc_man_multi)




######################################review variants flagged for multiallelic manual review in ifc_multi
#REVIEW POINT - if appropriate add the ifc_multi calls
allvariants <- dplyr::bind_rows(ifc, ifc_multi)

#REVIEW POINT - if appropriate add the man_review_variants calls - add if they are in cosmic?
#allvariants <- dplyr::bind_rows(allvariants, man_review_variants)

#############################################
#check non whitelist genes - STAT3, STAT5, MYD88, TERT
setwd(paste0('~/Documents/UCL/PhD/ELSA_variant_calling/mutect2/IFC_',ifc_no))

all_calls <- read.csv(paste0("IFC_", ifc_no, "_allvariants_combined.csv"))

stat3 <- length(grep("STAT3", all_calls$AAChange.refGene))
print(stat3)
stat5 <- length(grep("STAT5", all_calls$AAChange.refGene)) 
print(stat5)
myd88 <- length(grep("MYD88", all_calls$AAChange.refGene)) 
print(myd88)
tert <- length(grep("TERT", all_calls$AAChange.refGene)) 
print(tert)
zbtb33 <- length(grep("ZBTB33", all_calls$AAChange.refGene))
print(zbtb33)

#############################################
#order variant calls
allvariants <- allvariants[order(allvariants$Sample),]


###################
#annotate with filters applied
setwd(paste0('~/Documents/UCL/PhD/ELSA_variant_calling/mutect2/IFC_',ifc_no))

#state filters applied
filters <- paste0("DP=50, Alt_DP=20, min 1 read on all strands, base_qual applied, manual review variants filtered and added if in cosmic, u2af1 fixed genome, stat3 variants:", stat3, ", stat5 variants:", stat5, ", tert variants:", tert, ", zbtb33 variants", zbtb33)
filters <- append(filters, paste0("samples with variants requiring manual review passing filters but not in cosmic:", man_review_variants$Sample))

list_of_datasets <- list("variants" = allvariants, "filters" = filters, "manual_review" = man_review_variants)
write.xlsx(list_of_datasets, file = paste0("IFC_",ifc_no, "_filtered_post_wl.xlsx"), rowNames = FALSE)

