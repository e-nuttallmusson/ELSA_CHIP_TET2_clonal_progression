library(openxlsx)
library(plyr)
library(dplyr)
library(stringr)

#import and combine mean coverage per sample for each IFC

IFC_list <- c(2:38)
IFC_list <- append(IFC_list, c(42:77))
IFC_list <- append(IFC_list, 100)

setwd("~/Documents/UCL/PhD/ELSA_variant_calling/mutect2/coverage_summaries")
coverage_all <- read.csv("IFC_1_coverage_summary.csv")

for (i in IFC_list) {
  a <- read.csv(paste0("IFC_", i, "_coverage_summary.csv"))
  coverage_all <- rbind(coverage_all, a)
  rm(a)
}


#import and combine all calls post whitelist filtering
setwd("~/Documents/UCL/PhD/ELSA_variant_calling/mutect2")
post_wl_all <- read.xlsx("./IFC_1/IFC_1_filtered_post_wl.xlsx")

IFC_list <- c(2:38)
IFC_list <- append(IFC_list, c(42:60))
IFC_list <- append(IFC_list, c(62:77))
IFC_list <- append(IFC_list, 100)

for (i in IFC_list) {
  a <- read.xlsx(paste0("./IFC_", i,"/IFC_", i, "_filtered_post_wl.xlsx"))
  post_wl_all <- rbind(post_wl_all, a)
  rm(a)
}

write.xlsx(post_wl_all, "post_wl_filter_all_IFCs.xlsx")

#make variable indicating if sample has variants called by whitelist to add to coverage_all
whitelist_variant <- c()

for (row in 1:nrow(coverage_all)) {
  x <- coverage_all$Sample[row] %in% post_wl_all$Sample
  whitelist_variant <- append(whitelist_variant, x)
}

whitelist_variant <- as.character(whitelist_variant)

coverage_all$whitelist_variant <- whitelist_variant


#import and combine all manual review calls post whitelist filtering

setwd("~/Documents/UCL/PhD/ELSA_variant_calling/mutect2")
post_wl_man_rev_all <- read.xlsx("./IFC_4/IFC_4_filtered_post_wl.xlsx", sheet=3)

IFC_list <- c(8:38)
IFC_list <- append(IFC_list, c(42:60))
IFC_list <- append(IFC_list, c(62:77))
IFC_list <- append(IFC_list, 100)

for (i in IFC_list) {
  a <- read.xlsx(paste0("./IFC_", i,"/IFC_", i, "_filtered_post_wl.xlsx"), sheet=3)
  post_wl_man_rev_all <- rbind(post_wl_man_rev_all, a)
  rm(a)
}

manrev_variant <- c()

for (row in 1:nrow(coverage_all)) {
  x <- coverage_all$Sample[row] %in% post_wl_man_rev_all$Sample
  manrev_variant <- append(manrev_variant, x)
}

manrev_variant <- as.character(manrev_variant)

coverage_all$manrev_variant <- manrev_variant


#filter all_calls by filtering strategy in blacklist filter without filtering by whitelist. Set asxl1 filter for 646fs to <5%, VAF limit 0.01


setwd("~/Documents/UCL/PhD/ELSA_variant_calling/mutect2/")
IFC_list <- c(2:38)
IFC_list <- append(IFC_list, c(42:60))
IFC_list <- append(IFC_list, c(62:77))
IFC_list <- append(IFC_list, 100)

all_calls <- read.csv("./mono_all_calls/IFC_1_all_mono_variants.csv")

for (i in IFC_list) {
  a <- read.csv(paste0("./mono_all_calls/IFC_", i, "_all_mono_variants.csv"))
  all_calls <- rbind.fill(all_calls, a)
  rm(a)
  b <- read.csv(paste0("./multi_all_calls/IFC_", i, "_all_multi_variants.csv"))
  all_calls <- rbind.fill(all_calls, b)
  rm(b)
}

#does sample have variant called in all_calls - add to coverage dataframe
any_variant <- c()

for (row in 1:nrow(coverage_all)) {
  x <- coverage_all$Sample[row] %in% all_calls$Sample
  any_variant <- append(any_variant, x)
}

any_variant <- as.character(any_variant)

coverage_all$any_variant <- any_variant



#remove germline variants
all_calls<-dplyr::filter(all_calls, !grepl("germline",Otherinfo10))

#a lot of commonly occuring variants e.g. in CBL which fail slippage or contamination
#filters based on mutect2 output
##########
all_calls_filtered2 <- dplyr::filter(all_calls, !grepl("contamination",Otherinfo10))
all_calls_filtered2 <- dplyr::filter(all_calls_filtered2, !grepl("^synonymous SNV",ExonicFunc.refGene))
all_calls_filtered2 <- dplyr::filter(all_calls_filtered2, !grepl("base_qual",Otherinfo11))
all_calls_filtered2 <- dplyr::filter(all_calls_filtered2, !grepl("weak_evidence",Otherinfo11))
all_calls_filtered2 <- dplyr::filter(all_calls_filtered2, !grepl("contamination",Otherinfo11))
########

#only consider variants supported on forward and reverse
###########
all_calls_filtered2 <- dplyr::filter(all_calls_filtered2, SB_1>0)
all_calls_filtered2 <- dplyr::filter(all_calls_filtered2, SB_2>0)
all_calls_filtered2 <- dplyr::filter(all_calls_filtered2, SB_3>0)
all_calls_filtered2 <- dplyr::filter(all_calls_filtered2, SB_4>0)


#Filter out slippage variants that are present at a higher frequency than ASXL1:G646Wfs*12 and VAF <0.05 

CBL_D460del <- dplyr::filter(slippage, grepl("CBL:D460del", gene_NonsynOI))
EZH2_D189del <- dplyr::filter(slippage, grepl("EZH2:D189del", gene_NonsynOI))
TET2_P401del <- dplyr::filter(slippage, grepl("TET2:P401del", gene_NonsynOI))
TET2_Q1542del <- dplyr::filter(slippage, grepl("TET2:Q1542del", gene_NonsynOI))
TET2_Q1548del <- dplyr::filter(slippage, grepl("TET2:Q1548del", gene_NonsynOI))

CBL_D460del_hivaf <- dplyr::filter(CBL_D460del, AF > 0.05)
EZH2_D189del_hivaf <- dplyr::filter(EZH2_D189del, AF > 0.05)
TET2_P401del_hivaf <- dplyr::filter(TET2_P401del, AF > 0.05)
TET2_Q1542del_hivaf <- dplyr::filter(TET2_Q1542del, AF > 0.05)
TET2_Q1548del_hivaf <- dplyr::filter(TET2_Q1548del, AF > 0.05)

#remove all slippage calls and then add those with vaf >0.05 back in
all_calls_filtered2_noslip <- dplyr::filter(all_calls_filtered2, !grepl("slippage",Otherinfo10))
all_calls_filtered2_noslip$gene_NonsynOI <- gsub(" ", "", paste(all_calls_filtered2_noslip$Gene.refGene, ":", all_calls_filtered2_noslip$NonsynOI))

all_calls_hivaf_slip <- rbind(all_calls_filtered2_noslip, CBL_D460del_hivaf, EZH2_D189del_hivaf, TET2_P401del_hivaf, TET2_Q1542del_hivaf, TET2_Q1548del_hivaf)

###########
#create separate dfs for whitelist vs non whitelist

whitelist <- dplyr::filter(all_calls_hivaf_slip, grepl("TRUE",whitelist))

non_whitelist <- dplyr::filter(all_calls_hivaf_slip, grepl("FALSE",whitelist))


non_whitelist_cosmic <- dplyr::filter(non_whitelist, grepl("ID",cosmic70))


length(unique(all_calls_hivaf_slip$NonsynOI))

#number of repeated variants observed

var_freq <- dplyr::count(non_whitelist, Gene.refGene, ExonicFunc.refGene, NonsynOI, sort = TRUE)

var_freq_snv <- dplyr::filter(var_freq, grepl("nonsynonymous SNV", ExonicFunc.refGene))

#Use polyphen for non whitelist snvs post filtering to identify likely benign variants which are unlikely to have a phenotype
pos <- c()
AA1 <- c()
AA2 <- c()

for (i in (1:nrow(var_freq_snv))) {
  var_freq_snv$pos[i] <- substr(var_freq_snv$NonsynOI[i], 2, nchar(var_freq_snv$NonsynOI[i])-1)
  var_freq_snv$AA1[i] <- substr(var_freq_snv$NonsynOI[i], 1, 1)
  var_freq_snv$AA2[i] <- substr(var_freq_snv$NonsynOI[i], nchar(var_freq_snv$NonsynOI[i]), nchar(var_freq_snv$NonsynOI[i]))
}

for (i in nrow(var_freq_snv)) {
  a <- substr(var_freq_snv$NonsynOI[i], 2, nchar(var_freq_snv$NonsynOI[i])-1)
  pos <- append(pos, a)
}

#select necessary parameters for polyphen

polyphen_query <- var_freq_snv[,c("Gene.refGene", "pos", "AA1", "AA2")]
polyphen_query <- polyphen_query %>% dplyr::mutate(Gene.refGene = str_replace(Gene.refGene, "U2AF1;U2AF1L5", "U2AF1"))

write.table(polyphen_query, "all_calls_snv_polyphen.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


#import polyphen output
setwd("~/Documents/UCL/PhD/ELSA_variant_calling/mutect2/Polyphen_files")
polyphen_output <- read.table("all_IFC_polyphen_snv_output.txt", header = FALSE, sep = "\t")
polyphen_output<- polyphen_output %>% dplyr::rename(Gene.refGene = V1, pos = V2, AA1 = V3, AA2 = V4)

#merge polyphen output and var_freq
polyphen_output$pos <- as.character(polyphen_output$pos)
polyphen_output$NonsynOI <- gsub(" ", "", paste(polyphen_output$AA1, polyphen_output$pos, polyphen_output$AA2))
var_freq_snv_polyphen <- polyphen_query %>% right_join(polyphen_output, by = c("Gene.refGene", "pos", "AA1", "AA2"))


#annotate non whitelist variants with polyphen 2 output

#make new variable concatenanting gene and NonsynOI
var_freq_snv_polyphen$gene_NonsynOI <- gsub(" ", "", paste(var_freq_snv_polyphen$Gene.refGene, ":", var_freq_snv_polyphen$NonsynOI))
non_whitelist$gene_NonsynOI <- gsub(" ", "", paste(non_whitelist$Gene.refGene, ":", non_whitelist$NonsynOI))

var_polyphen <- var_freq_snv_polyphen[, c("gene_NonsynOI", "V10")]

#change nonwhitelist U2AF1 notation prior to merging
non_whitelist <- non_whitelist %>% dplyr::mutate(Gene.refGene = str_replace(Gene.refGene, "U2AF1;U2AF1L5", "U2AF1"))

non_whitelist <- merge(non_whitelist, var_polyphen, by = "gene_NonsynOI", all = TRUE)

non_whitelist_damaging <- dplyr::filter(non_whitelist, !grepl("benign",V10))

#add variable to coverage df detailing if sample has call in filtered call (filtered for hard characteristics and mutect2 filters)
filtered_variant <- c()

for (row in 1:nrow(coverage_all)) {
  x <- coverage_all$Sample[row] %in% all_calls_hivaf_slip$Sample
  filtered_variant <- append(filtered_variant, x)
}

filtered_variant <- as.character(filtered_variant)

coverage_all$filtered_variant <- filtered_variant


#add variable to coverage df detailing if sample has likely pathogenic non whitelist variant

damaging_nonwhitelist <- c()

for (row in 1:nrow(coverage_all)) {
  x <- coverage_all$Sample[row] %in% non_whitelist_damaging$Sample
  damaging_nonwhitelist <- append(damaging_nonwhitelist, x)
}

damaging_nonwhitelist <- as.character(damaging_nonwhitelist)

coverage_all$damaging_nonwhitelist <- damaging_nonwhitelist

#add variable to coverage df detailing if non whitelist variants have a cosmicID

cosmic_nonwhitelist <- c()

for (row in 1:nrow(coverage_all)) {
  x <- coverage_all$Sample[row] %in% non_whitelist_cosmic$Sample
  cosmic_nonwhitelist <- append(cosmic_nonwhitelist, x)
}

cosmic_nonwhitelist <- as.character(cosmic_nonwhitelist)

coverage_all$cosmic_nonwhitelist <- cosmic_nonwhitelist


#####section of code to import variants published in IPSS-M MDS paper
setwd("~/Documents/UCL/PhD/cosmic_mds/")
ipssm_all <- read.xlsx("ipssm_cosmic.xlsx")

ipssm_all$protein_annovar <- ipssm_all$PROTEIN_CHANGE

for (i in (1:nrow(ipssm_all))) {
  if (ipssm_all$VT[i] == "Sub") {
    ipssm_all$protein_annovar[i] <- gsub("\\*", "X", ipssm_all$protein_annovar[i])
  }
}

ipssm_all$gene_protein <- gsub(" ", "", paste(ipssm_all$GENE, ":", ipssm_all$protein_annovar))
ipssm_all$gene_protein <- gsub("p.", "", ipssm_all$gene_protein)

#create variable detailing if variant is in ipssm cohort

ipssm_variants <- ipssm_all$gene_protein

all_calls_hivaf_slip$in_ipssm <- c(rep(NA, nrow(all_calls_hivaf_slip)))

for (row in (1:nrow(all_calls_hivaf_slip))) {
  all_calls_hivaf_slip$in_ipssm[row] <- all_calls_hivaf_slip$gene_NonsynOI[row] %in% ipssm_variants
}


#create vector detailing if sample has variant in ipssm cohort

sample_ipssm_list <- c()

for (row in (1:nrow(coverage_all))) {
  trial <- all_calls_hivaf_slip[all_calls_hivaf_slip$Sample == coverage_all$Sample[row],]
  a <- TRUE %in% trial$in_ipssm
  sample_ipssm_list <- append(sample_ipssm_list, a)
  rm(trial)
  rm(a)
}

coverage_all$in_ipssm <- sample_ipssm_list

#find samples without whitelist variant, manual review variant or damaging non whitelist variant:

no_variants <- dplyr::filter(coverage_all, grepl("FALSE", damaging_nonwhitelist))
no_variants <- dplyr::filter(no_variants, grepl("FALSE", whitelist_variant))
no_variants <- dplyr::filter(no_variants, grepl("FALSE", manrev_variant))
no_variants <- dplyr::filter(no_variants, grepl("FALSE", cosmic_nonwhitelist))
no_variants <- dplyr::filter(no_variants, grepl("FALSE", in_ipssm))

no_variants_deep <- dplyr::filter(no_variants, Mean>800)

#also create a list of samples without any variants called, regardless of polyphen output - ideal candidate controls

no_variants_ideal <- dplyr::filter(coverage_all, grepl("FALSE", filtered_variant))
no_variants_ideal <- dplyr::filter(no_variants_ideal, grepl("FALSE", whitelist_variant))
no_variants_ideal <- dplyr::filter(no_variants_ideal, grepl("FALSE", manrev_variant))
no_variants_ideal <- dplyr::filter(no_variants_ideal, grepl("FALSE", cosmic_nonwhitelist))
no_variants_ideal <- dplyr::filter(no_variants_ideal, grepl("FALSE", in_ipssm))

no_variants_ideal_deep <- dplyr::filter(no_variants_ideal, Mean>800)

#output no_variants_deep as possible control candidates

write.xlsx(no_variants_deep, "control_candidates_post_filter_polyphen_all_calls.xlsx")

write.xlsx(no_variants_ideal_deep, "control_candidates_post_filter_ideal.xlsx")

#check
cosmic <- c()

for (row in 1:nrow(no_variants_ideal)) {
  x <- no_variants_ideal$Sample[row] %in% non_whitelist_cosmic$Sample
  cosmic <- append(cosmic, x)
}

table(cosmic)

###exploring stopgain calls erroneously filtered out by whitelist

stopgains <- dplyr::filter(all_calls_hivaf_slip, grepl("stopgain", ExonicFunc.refGene))

#remove stopgains that were identified by whitelist filter and ones that don't appear on the whitelist (have X but not on whitelist)
stopgains_missed <- dplyr::filter(stopgains, grepl("FALSE", whitelist))
stopgains_missed <- dplyr::filter(stopgains_missed, !grepl("X", NonsynOI))
stopgains_missed <- dplyr::filter(stopgains_missed, !grepl("PPM1D", Gene.refGene))
stopgains_missed <- dplyr::filter(stopgains_missed, !grepl("CBL", Gene.refGene))

setwd("~/Documents/UCL/PhD/ELSA_variant_calling/mutect2")
write.xlsx(stopgains_missed, "./missed_nonsense/whitelist_missed_nonsense.xlsx")


#####

#add in stopgains missed by whitelist to list of post_wl_all
stopgains_missed_merge <- subset(stopgains_missed, select=-c(gene_NonsynOI, in_ipssm))
post_wl_all_w_stopgains <- rbind(post_wl_all, stopgains_missed_merge)

post_wl_all_w_stopgains <- post_wl_all_w_stopgains[order(post_wl_all_w_stopgains$Sample),]

write.xlsx(post_wl_all_w_stopgains, "all_wl_all_IFC_with_missed_stopgains.xlsx")


###import sample id with elsa lab number and wave number

setwd("~/Documents/UCL/PhD/ELSA_variant_calling/mutect2/")

sample_info <- read.xlsx("ELSA_sample_info.xlsx")

#find missing value in all_calls - S165-100 (not demultiplexed by sequncing company, )
missing <- sample_info$`sample-id`[!(sample_info$`sample-id` %in% coverage_all$Sample)]
s165 <- c(0, 0, 0, 0, 0, 0, "S165-100", rep(FALSE, 7))
coverage_all <- rbind(coverage_all, s165)

elsa_wave <- c(rep(NA, nrow(sample_info)))
sample_info <- cbind(sample_info, elsa_wave)

#some samples have ugly strings and not classifiable

setwd("~/Documents/UCL/PhD/ELSA_sample_qc/wave_samples")
wave2 <- read.xlsx("ELSA Wave 2 DNA Samples Source Bioscience Plates 1-72 .xlsx")
wave4 <- read.xlsx("ELSA Wave 4 DNA Samples Source Bioscience Plates 73-93.xlsx")
wave6 <- read.xlsx("ELSA Wave 6 EDTA Blood Samples.xlsx")
wave8 <- read.xlsx("elsa_waves_8and9.xlsx", sheet = "Wave 8")
wave9 <- read.xlsx("elsa_waves_8and9.xlsx", sheet = "Wave 9")


#annotate for wave 2
for (i in (1:nrow(sample_info))) {
  if (sample_info$Sample[i] %in% wave2$Sample.ID) {
    sample_info$elsa_wave[i] <- 2
  }
}

#annotate for wave 4
for (i in (1:nrow(sample_info))) {
  if (sample_info$Sample[i] %in% wave4$Sample.ID) {
    sample_info$elsa_wave[i] <- 4
  }
}

#annotate for wave 6
for (i in (1:nrow(sample_info))) {
  if (sample_info$Sample[i] %in% wave6$Name) {
    sample_info$elsa_wave[i] <- 6
  }
}

#annoate for wave 8
for (i in (1:nrow(sample_info))) {
  if (sample_info$Sample[i] %in% wave8$Name) {
    sample_info$elsa_wave[i] <- 8
  }
}

#annotate for wave 9
for (i in (1:nrow(sample_info))) {
  if (sample_info$Sample[i] %in% wave9$Name) {
    sample_info$elsa_wave[i] <- 9
  }
}


sample_info <- sample_info %>% select(Sample, `sample-id`, elsa_wave, duplicate)
names(sample_info)[2] <- "Sample"
names(sample_info)[1] <- "lab_no"

#merge coverage_all and sample info
coverage_all2 <- merge(coverage_all, sample_info, by="Sample")

#remove non ELSA samples from cases and controls
#coverage
coverage_elsa <- coverage_all2[!(is.na(coverage_all2$elsa_wave)), ]

#cases
post_wl_all_w_stopgains_samples <- merge(post_wl_all_w_stopgains, sample_info, by="Sample")
cases_elsa <- post_wl_all_w_stopgains_samples[!(is.na(post_wl_all_w_stopgains_samples$elsa_wave)), ]


#controls
#candidate controls - no damaging polyphen variants
no_variants_deep_samples <- merge(no_variants_deep, sample_info, by="Sample")
controls_elsa_deep <- no_variants_deep_samples[!(is.na(no_variants_deep_samples$elsa_wave)), ]

#ideal candidate controls
no_variants_ideal_deep_samples <- merge(no_variants_ideal_deep, sample_info, by="Sample")
controls_elsa_ideal_deep <- no_variants_ideal_deep_samples[!(is.na(no_variants_ideal_deep_samples$elsa_wave)), ]


#Ensure no variants with vaf < 0.02

for (i in (1:nrow(cases_elsa))) {
  if (!is.na(cases_elsa$AF[i])) {
    if (cases_elsa$AF[i] < 0.02) {
      cases_elsa$lowvaf[i] <- TRUE
    }
  }
}

cases_elsa <- dplyr::filter(cases_elsa, !grepl(TRUE, lowvaf))

#write out case and control files



setwd("~/Documents/UCL/PhD/ELSA_variant_calling/mutect2/summary_files/elsa_samples")
write.xlsx(coverage_elsa, "coverage_elsa.xlsx")
write.xlsx(cases_elsa, "chip_cases_whitelist_elsa.xlsx")
write.xlsx(controls_elsa_deep, "controls_elsa_deep.xlsx")
write.xlsx(controls_elsa_ideal_deep, "ideal_controls_elsa_deep.xlsx")

#to create files to export to ELSA
elsa_chip_samples <- subset(cases_elsa, select = c("Sample", "lab_no", "elsa_wave"))
elsa_chip_samples <- dplyr::distinct(elsa_chip_samples)

elsa_control_samples <- subset(elsa_controls, select = c("Sample", "lab_no", "elsa_wave"))


list_of_datasets <- list("chip_cases" = elsa_chip_samples, "controls" = elsa_control_samples)
write.xlsx(list_of_datasets, file = "elsa_desired_samples.xlsx")
