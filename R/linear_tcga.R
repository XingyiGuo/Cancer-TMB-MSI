######## MSI-TMB project #########
######## Data preparation ########
library(data.table)
library(tidyverse)
library(dplyr)

# 3. TCGA
# We do note have Sample_type in TCGA data.
# cliS$`American Joint Committee on Cancer Metastasis Stage Code` Excluding all Mx data resulted in a loss of approximately 50% of the data. Therefore, for the TCGA dataset, we chose to include all M0, M1, and Mx data.

a<-fread("/Users/kumc/Dropbox/MSI-TMB/Data/MSI_core.csv",head=TRUE,stringsAsFactors=F)
a$Tumor_Sample_Barcode <- a$Barcode
a1<-a[a$MSI=="msi-h",]
a1$MSI_status<-"MSI-H"
a2<-a[a$MSI!="msi-h",]
a2$MSI_status<-"Non-MSI-H"
a<-rbind(a1, a2)

b<-fread("/Users/kumc/Dropbox/MSI-TMB/Data/mutation-load_updated.txt",head=TRUE,stringsAsFactors=F)
b$Tumor_Sample_Barcode<-b$Tumor_Sample_ID
b$Tumor_Sample_Barcode<-substr(b$Tumor_Sample_Barcode,1,12)
b$tmb <- b$`Non-silent per Mb`
b1<-b[(b$`Non-silent per Mb`)>=10,]
b1$TMB_status <- "TMB-high"
b2<-b[(b$`Non-silent per Mb`)<10,]
b2$TMB_status <- "TMB-low"
b<-rbind(b1,b2)

cliM<-fread("/Users/kumc/Dropbox/MSI-TMB/Data/TCGA_raw_genotype/STAD_TCGA.csv",head=TRUE,stringsAsFactors=F) #COREAD, SKCM, KIRC, KIRP ####
cliM$Tumor_Sample_Barcode<-substr(cliM$Tumor_Sample_Barcode,1,12) 

cliM<-merge(cliM,a,by="Tumor_Sample_Barcode")
cliM<-merge(cliM,b,by="Tumor_Sample_Barcode")

cliS<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/stad_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE) ####COADREAD..csv
cliS<-cliS[!duplicated(cliS$`Patient ID`),] ###########
cliS$AGE <- cliS$`Diagnosis Age`
table(cliS$AGE)
cliS<-cliS[!(cliS$AGE %in% NA),]
cliS$AGE <- as.numeric(as.character(cliS$AGE))

cliS$PRIMARY_RACE <- cliS$`Race Category`
cliS$PRIMARY_RACE <- ifelse(cliS$PRIMARY_RACE %in% "Asian","Asian",ifelse(cliS$PRIMARY_RACE %in% "Black or African American","BLack",ifelse(cliS$PRIMARY_RACE %in% "White","White","Others")))

cliS<-cliS[!(cliS$PRIMARY_RACE %in% 'Others'),]

table(cliS$PRIMARY_RACE)

# Genotype data mining
mutation<-c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_DEL","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation") ###GENIE

mut.dat<-cliM
sub.mut<- mut.dat[mut.dat$Tumor_Sample_Barcode %in% cliS$`Patient ID`,]
sub.mut <- as.data.frame(sub.mut)
table(sub.mut$Variant_Classification)
sub.mut$SAMPLE_ID<-sub.mut$Tumor_Sample_Barcode

sub.mut<-sub.mut[sub.mut$Variant_Classification %in% mutation,]
sub.mut$SAMPLE_ID
cliS$SAMPLE_ID<-cliS$`Patient ID`
A.matrix.1 <- merge(sub.mut, cliS, by = "SAMPLE_ID")

#### add samples without mutations detected ####
x2 <- cliS[!(cliS$SAMPLE_ID %in% sub.mut$SAMPLE_ID),] ## to add samples without mutations 
sub.mut.1 <- sub.mut[1:dim(x2)[1],] ### total samples 
sub.mut.1$SAMPLE_ID  <- x2$SAMPLE_ID  
sub.mut.1$Hugo_Symbol<-"ADD"
A.matrix.2 <- merge(sub.mut.1,x2,by = "SAMPLE_ID")    
#### end #######################################

#### merge A.matrix + A.matrix.1 ###############
A.matrix <- rbind(A.matrix.2,A.matrix.1)
A.matrix <- A.matrix[(A.matrix$Tumor_Sample_Barcode %in% a$Tumor_Sample_Barcode),]######

#### main analysis I: Pan-cancers / baselines ####
A.base <- A.matrix
A.base1 <- A.base
#### end #######################################

A.base2 <- A.base1[!duplicated(A.base1$SAMPLE_ID),]
A.base2 <- A.base1[!duplicated(A.base1$Tumor_Sample_Barcode),]
A.base1 <- A.base1[!duplicated(A.base1$Tumor_Sample_Barcode),]
#### grouping ##################################
#TMB=0 and MSI=0 then group=0
#TMB=0 and MSI=1 then group=1
#TMB=1 and MSI=0 then group=2
#TMB=1 and MSI=1 then group=3

a1<-A.base1[A.base1$TMB_status=="TMB-low" & A.base1$MSI_status=="Non-MSI-H",]
a1$group<-0
b<-A.base1[A.base1$TMB_status=="TMB-low" & A.base1$MSI_status=="MSI-H",]
b$group<-1
c<-A.base1[A.base1$TMB_status=="TMB-high" & A.base1$MSI_status=="Non-MSI-H",]
c$group<-2
d<-A.base1[A.base1$TMB_status=="TMB-high" & A.base1$MSI_status=="MSI-H",]
d$group<-3
total<-rbind(a1,b,c,d)

t<-total[!duplicated(total$Tumor_Sample_Barcode),]

#### Linear association ########################
test<-t[(t$`Cancer Type`=="Esophagogastric Cancer"&t$PRIMARY_RACE=="Asian"),] #Asian BLack White 

Aout <- summary(lm(tmb ~ relevel(as.factor(MSI_status),ref="Non-MSI-H") + `Diagnosis Age` + `Person Gender`+`Cancer Type Detailed`+`Center of sequencing`, data = test))

Aout
Aout$coefficients[2, 4] #P-vale
