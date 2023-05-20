######## MSI-TMB project #########
######## Data preparation ########
library(data.table)
library(tidyverse)
library(dplyr)

# 2. K-MASTER
a<-fread("/Users/kumc/Dropbox/MSI-TMB/Data/a_km_22pancancer.csv",head=TRUE,stringsAsFactors=F) #MSI,TMB: tmb cut-off 10

#a<-subset(a,select=-c(V1,TMB_status)) #n=5418
#if TMB-H cut-off set as >10
#a1<-a[a$tmb>=10,]
#a1$TMB_status<-"TMB-high"
#a2<-a[a$tmb<=10,]
#a2$TMB_status<-"TMB-low"
#a<-rbind(a1,a2)
#min(a$tmb[a$TMB_status=="TMB-high"])

a<-subset(a,select=c(SAMPLE_ID,MSI_status,tmb,TMB_status,CANCER_TYPE,CANCER_TYPE_DETAILED))
cliM<-fread("/Users/kumc/Dropbox/MSI-TMB/Data/KM_genotype2.csv",head=TRUE,stringsAsFactors=F) #Genotype
cliS<-fread("/Users/kumc/Dropbox/MSI-TMB/Data/KM_clinical_simple_updated.csv",head=TRUE,stringsAsFactors=F) #Clinical data

cliM<-subset(cliM,select=-c(V1))
cliM<-subset(cliM,select=-c(V1))#2
cliS<-subset(cliS,select=-c(V1))

# Clinical data mining
cliS<-subset(cliS,select=-c(CANCER_TYPE,CANCER_TYPE_DETAILED))
cliS<-cliS[!duplicated(cliS$PATIENT_ID),]

cliS$AGE <- cliS$AGE_AT_SEQ_REPORT
cliS<-cliS[!(cliS$AGE %in% 'Unknown'),]
cliS$AGE <- as.numeric(as.character(cliS$AGE))

cliS$PRIMARY_RACE <- ifelse(cliS$PRIMARY_RACE %in% "Asian","Asian",ifelse(cliS$PRIMARY_RACE %in% "Black","BLack",ifelse(cliS$PRIMARY_RACE %in% "White","White","Others")))
cliS<-cliS[!(cliS$PRIMARY_RACE %in% 'Others'),]

cliS$SAMPLE_TYPE <- ifelse(cliS$SAMPLE_TYPE %in% "Primary","Primary",ifelse(cliS$SAMPLE_TYPE %in% "Metastasis","Metastasis","Others"))
cliS<-cliS[!(cliS$SAMPLE_TYPE %in% 'Others'),]

# Genotype data mining
mutation<-c("frameshift_deletion","frameshift_insertion","nonframeshift_deletion","nonframeshift_insertion","nonsynonymous_SNV","stopgain","stoploss") ###KM

mut.dat<-cliM
sub.mut<- mut.dat[mut.dat$Tumor_Sample_Barcode %in% cliS$SAMPLE_ID,]
sub.mut <- as.data.frame(sub.mut)
table(sub.mut$Variant_Classification)
sub.mut$SAMPLE_ID<-sub.mut$Tumor_Sample_Barcode

sub.mut<-sub.mut[sub.mut$Variant_Classification %in% mutation,]
A.matrix.1 <- merge(sub.mut, cliS, by = "SAMPLE_ID")

#### add samples without mutations detected ####
x2 <- cliS[!(cliS$SAMPLE_ID %in% sub.mut$SAMPLE_ID),] ## to add samples without mutations 
sub.mut.1 <- sub.mut[1:dim(x2)[1],] ### total of samples: 385 -338 
sub.mut.1$SAMPLE_ID  <- x2$SAMPLE_ID  
sub.mut.1$Hugo_Symbol<-"ADD"
A.matrix.2 <- merge(sub.mut.1,x2,by = "SAMPLE_ID")    
#### end #######################################

#### merge A.matrix + A.matrix.1 ###############
A.matrix <- rbind(A.matrix.2,A.matrix.1)
A.matrix <- A.matrix[(A.matrix$SAMPLE_ID %in% a$SAMPLE_ID),]######

#### main analysis I: Pan-cancers / baselines ####
A.base <- A.matrix
A.base1 <- A.base
#### end #######################################

A.base2 <- A.base1[!duplicated(A.base1$SAMPLE_ID),]
A.base2 <- A.base1[!duplicated(A.base1$Tumor_Sample_Barcode),]

a$Tumor_Sample_Barcode<-a$SAMPLE_ID
#### grouping ##################################
#TMB=0 and MSI=0 then group=0
#TMB=0 and MSI=1 then group=1
#TMB=1 and MSI=0 then group=2
#TMB=1 and MSI=1 then group=3

A.base1<-merge(A.base1,a,by="Tumor_Sample_Barcode") #merge with tmb+msi

a1<-A.base1[A.base1$TMB_status=="TMB-low" & A.base1$MSI_status=="Non-MSI-H",]
a1$group<-0
b<-A.base1[A.base1$TMB_status=="TMB-low" & A.base1$MSI_status=="MSI-H",]
b$group<-1
c<-A.base1[A.base1$TMB_status=="TMB-high" & A.base1$MSI_status=="Non-MSI-H",]
c$group<-2
d<-A.base1[A.base1$TMB_status=="TMB-high" & A.base1$MSI_status=="MSI-H",]
d$group<-3
total<-rbind(a1,b,c,d)

t<-total[!duplicated(total$SubjectID),] 
t<-total[!duplicated(total$Tumor_Sample_Barcode),] #n=5356

#### Linear association ########################
table(t$CANCER_TYPE)

test<-t[t$CANCER_TYPE=="Stomach cancer",]

Aout <- summary(lm(tmb ~ relevel(as.factor(MSI_status),ref="Non-MSI-H") + AGE + SEX + CANCER_TYPE_DETAILED + SAMPLE_TYPE + SEQ_ASSAY_ID, data = test)) 

Aout
Aout$coefficients[2, 4] #P-vale
