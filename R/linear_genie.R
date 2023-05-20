######## MSI-TMB project #########
######## Data preparation ########
library(data.table)
library(tidyverse)
library(dplyr)

# 1. GENIE
a<-fread("/Users/kumc/Dropbox/MSI-TMB/Data/a_genie_22pancancer.csv",head=TRUE,stringsAsFactors=F) #MSI,TMB: tmb cut-off 10
cliM<-fread("/Users/kumc/Dropbox/MSI-TMB/Data/GENIE_genotype2.csv",head=TRUE,stringsAsFactors=F) #Genotype
cliS<-fread("/Users/kumc/Dropbox/MSI-TMB/Data/GENIE_clinical_simple.csv",head=TRUE,stringsAsFactors=F) #Clinical data

a<-subset(a,select=-c(V1))
cliM<-subset(cliM,select=-c(V1))
cliM<-subset(cliM,select=-c(V1))#2
cliS<-subset(cliS,select=-c(V1))

# Clinical data mining
cliS<-subset(cliS,select=-c(CANCER_TYPE,CANCER_TYPE_DETAILED))
cliS<-cliS[!duplicated(cliS$PATIENT_ID),]

cliS$AGE_AT_SEQ_REPORT <- gsub(">89","Unknown",cliS$AGE_AT_SEQ_REPORT,perl =TRUE)
cliS$AGE_AT_SEQ_REPORT <- gsub("<18","Unknown",cliS$AGE_AT_SEQ_REPORT,perl =TRUE)
cliS$AGE <- cliS$AGE_AT_SEQ_REPORT
cliS<-cliS[!(cliS$AGE %in% 'Unknown'),]
cliS$AGE <- as.numeric(as.character(cliS$AGE))

cliS$PRIMARY_RACE <- ifelse(cliS$PRIMARY_RACE %in% "Asian","Asian",ifelse(cliS$PRIMARY_RACE %in% "Black","BLack",ifelse(cliS$PRIMARY_RACE %in% "White","White","Others")))
cliS<-cliS[!(cliS$PRIMARY_RACE %in% 'Others'),]

cliS$SAMPLE_TYPE <- ifelse(cliS$SAMPLE_TYPE %in% "Primary","Primary",ifelse(cliS$SAMPLE_TYPE %in% "Metastasis","Metastasis","Others"))
cliS<-cliS[!(cliS$SAMPLE_TYPE %in% 'Others'),]

# Genotype data mining
mutation<-c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_DEL","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation") ###GENIE

mut.dat<-cliM
sub.mut<- mut.dat[mut.dat$Tumor_Sample_Barcode %in% cliS$SAMPLE_ID,]
sub.mut <- as.data.frame(sub.mut)
table(sub.mut$Variant_Classification)
sub.mut$SAMPLE_ID<-sub.mut$Tumor_Sample_Barcode

sub.mut<-sub.mut[sub.mut$Variant_Classification %in% mutation,]
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

# deleted unmatched samples
total <- total[!((total$CANCER_TYPE=="Ovarian Cancer"&total$SEX=="Male")|(total$CANCER_TYPE=="Prostate Cancer"&total$SEX=="Female")),] #unmatched=5
t<-total[!duplicated(total$SubjectID),] #n=46320
t<-total[!duplicated(total$Tumor_Sample_Barcode),] #n=46320

#### Linear association ########################
table(t$CANCER_TYPE)
table(test$CANCER_TYPE_DETAILED)
table(test$SEQ_ASSAY_ID)
table(test$SAMPLE_TYPE)
table(test$AGE)
table(test$TMB_status)
table(test$MSI_status)

test<-t[t$CANCER_TYPE=="Cervical Cancer",]

Aout <- summary(lm(tmb ~ relevel(as.factor(MSI_status),ref="Non-MSI-H") + AGE + SEX + PRIMARY_RACE + CANCER_TYPE_DETAILED + SAMPLE_TYPE + SEQ_ASSAY_ID, data = test)) 

Aout <- summary(lm(tmb ~ relevel(as.factor(MSI_status),ref="Non-MSI-H") + AGE + PRIMARY_RACE + CANCER_TYPE_DETAILED+ SAMPLE_TYPE + SEQ_ASSAY_ID, data = test)) 

Aout
Aout$coefficients[2, 4] #P-vale  

#### Binary association ########################
table(t$CANCER_TYPE)

test<-t[t$CANCER_TYPE=="Soft Tissue Sarcoma",]
#test$TMB_status<-relevel(as.factor(test$TMB_status),ref="TMB-low")

Aout <- summary(glm(relevel(as.factor(test$TMB_status),ref="TMB-low") ~ relevel(as.factor(MSI_status),ref="Non-MSI-H") + AGE + SEX + PRIMARY_RACE + CANCER_TYPE_DETAILED + SAMPLE_TYPE + SEQ_ASSAY_ID, data = test, family="binomial")) 

Aout <- summary(glm(relevel(as.factor(test$TMB_status),ref="TMB-low") ~ relevel(as.factor(MSI_status),ref="Non-MSI-H") + AGE + PRIMARY_RACE + CANCER_TYPE_DETAILED + SAMPLE_TYPE + SEQ_ASSAY_ID, data = test, family="binomial")) 

Aout
Aout$coefficients[2, 4] #P-vale  

#### Stratification by race: Linear association ####################################### 
table(t$CANCER_TYPE)

out.m<-matrix("na",25,13)
test1<-t[(t$CANCER_TYPE=="Endometrial Cancer"&t$PRIMARY_RACE=="White"),] #Asian BLack White
test2<-t[(t$CANCER_TYPE=="Endometrial Cancer"&t$PRIMARY_RACE=="BLack"),]
test3<-t[(t$CANCER_TYPE=="Endometrial Cancer"&t$PRIMARY_RACE=="Asian"),]

cancer<-"Endometrial"
i<-1
j<-0

table(test1$CANCER_TYPE_DETAILED)
table(test1$SEQ_ASSAY_ID)
table(test2$SAMPLE_TYPE)
table(test1$AGE)
table(test1$TMB_status)
table(test1$MSI_status)
table(test2$SEX)

Aout <- summary(lm(tmb ~ relevel(as.factor(MSI_status),ref="Non-MSI-H") + AGE + SEX + CANCER_TYPE_DETAILED + SAMPLE_TYPE + SEQ_ASSAY_ID, data = test1))

for(i in 1)
{
  Aout <- summary(lm(tmb ~ relevel(as.factor(MSI_status),ref="Non-MSI-H") + AGE + CANCER_TYPE_DETAILED +  SAMPLE_TYPE + SEQ_ASSAY_ID, data = test1))
  Aout2 <- summary(lm(tmb ~ relevel(as.factor(MSI_status),ref="Non-MSI-H") + AGE + CANCER_TYPE_DETAILED +  SAMPLE_TYPE + SEQ_ASSAY_ID, data = test2))
  Aout3 <- summary(lm(tmb ~ relevel(as.factor(MSI_status),ref="Non-MSI-H") + AGE + CANCER_TYPE_DETAILED +  SAMPLE_TYPE + SEQ_ASSAY_ID, data = test3))

#Aout <- summary(lm(tmb ~ relevel(as.factor(MSI_status),ref="Non-MSI-H") + AGE + SEX + CANCER_TYPE_DETAILED + SAMPLE_TYPE + SEQ_ASSAY_ID, data = test))

frq<- dim(test1)[1] 
frq2<- dim(test2)[1]  
frq3<- dim(test3)[1]  

sum.coef<- Aout$coef
sum.coef2<- Aout2$coef
sum.coef3<- Aout3$coef

j<-j+1

out.m[j,]<- c(cancer,frq,Aout$coeff[2,c(1,2,4)],frq2,Aout2$coeff[2,c(1,2,4)],frq3,Aout3$coeff[2,c(1,2,4)])
}	

out.m1<-out.m[1:25,]
out.m1<-as.data.frame(out.m1)

names(out.m1) <- c("cancer type","N_white","Beta","SE","P","N_black","Beta","SE","P","N_asian","Beta","SE","P")
out.m1$P<- as.numeric(as.character(out.m1$P))

write.csv(out.m1, "/Users/kumc/Desktop/genie.csv")
