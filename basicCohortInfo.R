##TCGA
pheno = data.table::fread("Data/tcga_clinical_data.tsv")

TCGAStudies = levels(as.factor(pheno$`TCGA PanCanAtlas Cancer Type Acronym`))

phenoCT = split(pheno, pheno$`TCGA PanCanAtlas Cancer Type Acronym`)

num = unlist(lapply(phenoCT, nrow))

nFem = unlist(lapply(phenoCT, function(x){length(which(x$Sex == "Female"))}))

medAge = unlist(lapply(phenoCT, function(x){median(x$`Diagnosis Age`, na.rm = T)}))

bInfo = rbind(num, nFem, medAge)

write.csv(bInfo, "Data/BasicCTInfoTCGA.csv")


##GTEx
sampleAnn = data.table::fread("Data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
sampleAnnSub = sampleAnn[, c(1,6,7)]

summary(as.factor(sampleAnn$SMAFRZE))
sampleAnnSub = sampleAnnSub[which(sampleAnn$SMAFRZE == "RNASEQ"),]
nSampTissue = summary(as.factor(sampleAnnSub$SMTS))

write.csv(nSampTissue, "Data/GTExSamplesPerTissue.csv")

##GENIE
GENIEdata = data.table::fread("Data/genie_public_clinical_data.tsv")

#remove NA mutation count  samples
sum(is.na(GENIEdata$`Mutation Count`))

GENIEdata2 = GENIEdata[-which(is.na(GENIEdata$`Mutation Count`)),]

#filter to cancer types in which ICB therapies are commonly used
summary(as.factor(GENIEdata2$`Cancer Type`))

GENIEdata3 = GENIEdata2[which(GENIEdata2$`Cancer Type` %in% c("Bladder Cancer", "Breast Cancer", "Colorectal Cancer", "Esophagogastric Cancer", "Head and Neck Cancer", "Melanoma", "Non-Small Cell Lung Cancer", "Renal Cell Carcinoma")),]

#use only one sample per patient
GENIEdata4 = GENIEdata3[which(GENIEdata3$`Number of Samples Per Patient` == 1),]

#change column names
colnames(GENIEdata4)[c(4,5, 10)] = c("Age", "CancerType", "MutationCount")

#remove samples with incompletely specified or unspecified ages
summary(as.factor(GENIEdata4$Age))
GENIEdata5 = GENIEdata4[-which(GENIEdata4$Age %in% c("<18", ">89", "Unknown")),]

#make age numeric
GENIEdata5$Age = as.numeric(GENIEdata5$Age)

#split by cancer type
GENIECT = split(GENIEdata5, GENIEdata5$CancerType)

nPatients = unlist(lapply(GENIECT, nrow))
nFem = unlist(lapply(GENIECT, function(x){length(which(x$Sex == "Female"))}))
medAge = unlist(lapply(GENIECT, function(x){median(x$Age)}))

bInfo = rbind(nPatients, nFem, medAge)

write.csv(bInfo, "Data/BasicCTInfoGENIE.csv")


#METABRIC
pheno = data.table::fread("D:/FertigLabLargeData/ImmuneAging/data_clinical_patient.txt")

pheno = pheno[-c(1:4),]

length(which(pheno$`Age at Diagnosis` < 30))
length(which(pheno$`Age at Diagnosis` >= 30 & pheno$`Age at Diagnosis` <= 50))
length(which(pheno$`Age at Diagnosis` > 50 & pheno$`Age at Diagnosis` <= 65))
length(which(pheno$`Age at Diagnosis` > 65 & pheno$`Age at Diagnosis` <= 80))
length(which(pheno$`Age at Diagnosis` > 80))

summary(as.factor(pheno$`3-Gene classifier subtype`))
