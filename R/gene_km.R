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

#### grouping ##################################
#TMB=0 and MSI=0 then group=0
#TMB=0 and MSI=1 then group=1
#TMB=1 and MSI=0 then group=2
#TMB=1 and MSI=1 then group=3

A.base1<-merge(A.base1,a,by="SAMPLE_ID") #merge with tmb+msi

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
t<-total[!duplicated(total$SAMPLE_ID),] #n=5356


#**********************************************#
######## Gene analysis #########################
#**********************************************#
# Repeat this process for each cancer type #

table(total$CANCER_TYPE)
t<-total[total$CANCER_TYPE=="Stomach cancer",]
test<-t[!duplicated(t$Tumor_Sample_Barcode),]

j <- 0
A.matrix<-t

### gene analysis ### 
gene<-fread("/Users/kumc/Dropbox/MSI-TMB/Data/gene_list/gene_km_stad.csv",head=TRUE,stringsAsFactors=F)

#KM-PANEL
PANEL<-fread("/Users/kumc/Dropbox/MSI-TMB/Data/PANEL_KMupdate.csv")

out.m<-matrix("na",40,20)

for(g1 in gene$Gene) 
{ 
  panel <- as.character(PANEL[PANEL$Hugo_Symbol %in% g1,]$SEQ_ASSAY_ID)
  
  Asub <- A.matrix[(A.matrix$Center %in% panel),] #KM Center=seq assay
  
  Asub$y<-ifelse((Asub$Hugo_Symbol %in% g1) & (Asub$Variant_Classification %in% mutation),1,0)
  Asub.1 <- Asub[Asub$y == 1,]
  Asub.2 <- Asub[!(Asub$SAMPLE_ID %in% Asub.1$SAMPLE_ID),]
  Asub.1$SAMPLE_ID.y
  Asub.1<-Asub.1[!duplicated(Asub.1$SAMPLE_ID),]
  Asub.2<-Asub.2[!duplicated(Asub.2$SAMPLE_ID),]
  Bsub<-rbind(Asub.1,Asub.2)
  
  frq<- dim(Asub.1)[1]/dim(Bsub)[1]
  frq1<- dim(Asub.1[Asub.1$group ==0,])[1]/dim(Bsub[Bsub$group==0,])[1]
  frq2<- dim(Asub.1[Asub.1$group ==1,])[1]/dim(Bsub[Bsub$group==1,])[1]
  frq3<- dim(Asub.1[Asub.1$group ==2,])[1]/dim(Bsub[Bsub$group==2,])[1]
  frq4<- dim(Asub.1[Asub.1$group ==3,])[1]/dim(Bsub[Bsub$group==3,])[1]
  My1<-dim(Asub.1[Asub.1$group == 0,])[1]
  My0<-dim(Bsub[Bsub$group == 0,])[1]
  Mo1<-dim(Asub.1[Asub.1$group == 1,])[1]
  Mo0<-dim(Bsub[Bsub$group == 1,])[1]
  Mt1<-dim(Asub.1[Asub.1$group == 2,])[1]
  Mt0<-dim(Bsub[Bsub$group == 2,])[1]
  Ms1<-dim(Asub.1[Asub.1$group == 3,])[1]
  Ms0<-dim(Bsub[Bsub$group == 3,])[1]
  
  Aout <- summary(glm(relevel(as.factor(MSI_status),ref="Non-MSI-H") ~ relevel(as.factor(y),ref="0") + AGE + SEX + CANCER_TYPE_DETAILED + SAMPLE_TYPE + SEQ_ASSAY_ID + tmb, data = Bsub, family="binomial")) #CANCER_TYPE_DETAILED, SEX
  
  sum.coef<- Aout$coef
  j<-j+1
  
  est<-exp(sum.coef[2,1])
  upper.ci<-exp(sum.coef[2,1]+1.96*sum.coef[2,2])
  lower.ci<-exp(sum.coef[2,1]-1.96*sum.coef[2,2])
  
  out.m[j,]<- c(g1,Aout$coeff[2,c(1,2,4)],est,lower.ci,upper.ci,frq,frq1,frq2,frq3,frq4,My1,My0,Mo1,Mo0,Mt1,Mt0,Ms1,Ms0)
}	

j
table(Bsub$y)
out.m1<-out.m[1:40,]
out.m1<-as.data.frame(out.m1)

names(out.m1) <- c("Gene","Beta","SE","P","OR","95CI1","95CI2","Freq","Freq_0","Freq_1","Freq_2","Freq_3","#_Carrie_0","#_0All","#_Carrie_1","#_1All","#_Carrie_2","#_2All","#_Carrie_3","#_3All")

out.m1$P<- as.numeric(as.character(out.m1$P))
out.m1$P_BF<-p.adjust(out.m1$P,method="bonferroni")

write.csv(out.m1, "/Users/kumc/Desktop/gene_raw_genie/stad_km.csv")
