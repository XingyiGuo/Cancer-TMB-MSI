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
A.base1 <- A.base #n=444684
#### end #######################################

#A.base1$SAMPLE_ID[A.base1$CENTER=="DFCI"]
#A.base1$SubjectID[A.base1$CENTER=="DFCI"]
#A.base1$Tumor_Sample_Barcode[A.base1$CENTER=="DFCI"]
#a$Tumor_Sample_Barcode

A.base2 <- A.base1[!duplicated(A.base1$SAMPLE_ID),] #n=47292
A.base2 <- A.base1[!duplicated(A.base1$Tumor_Sample_Barcode),] #n=46324

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
total<-rbind(a1,b,c,d) #n=444684

# deleted unmatched samples
total <- total[!((total$CANCER_TYPE=="Ovarian Cancer"&total$SEX=="Male")|(total$CANCER_TYPE=="Prostate Cancer"&total$SEX=="Female")),] #unmatched=5 #n=444667
t<-total[!duplicated(total$SubjectID),] #n=46320
t<-total[!duplicated(total$Tumor_Sample_Barcode),] #n=46320


#**********************************************#
######## Gene analysis #########################
#**********************************************#

table(total$CANCER_TYPE)
t<-total[total$CANCER_TYPE=="Colorectal Cancer",]
j <- 0
A.matrix<-t

### gene analysis ### 
gene <- as.data.frame(table(A.matrix$Hugo_Symbol))

gene <- gene[((gene$Freq > 350) | 
                gene$Var1=="MLH1" | gene$Var1=="MLH3" | gene$Var1=="MSH2"| gene$Var1=="MSH3"| gene$Var1=="MSH6"
              | gene$Var1=="PMS1" | gene$Var1=="PMS2" | gene$Var1=="POLD1"| gene$Var1=="POLE"| gene$Var1=="PRKDC"
              | gene$Var1=="APC"| gene$Var1=="BRAF"),] #200
gene <- gene[!(as.character(gene$Var1) %in% 'ADD'),]

#GENIE-PANEL
cliS$Tumor_Sample_Barcode<-cliS$SAMPLE_ID
p<-merge(cliM,cliS,by="Tumor_Sample_Barcode")
PANEL<-p

out.m<-matrix("na",60,20)

for(g1 in gene$Var1) 
{ 
panel <- as.character(PANEL[PANEL$Hugo_Symbol %in% g1,]$SEQ_ASSAY_ID)
#Asub <- A.matrix[(A.matrix$Center %in% panel),] #KM

Asub <- A.matrix[(A.matrix$SEQ_ASSAY_ID %in% panel),] #GENIE

Asub$y<-ifelse((Asub$Hugo_Symbol %in% g1) & (Asub$Variant_Classification %in% mutation),1,0)
Asub.1 <- Asub[Asub$y == 1,]
Asub.2 <- Asub[!(Asub$SAMPLE_ID %in% Asub.1$SAMPLE_ID),]
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

Aout <- summary(glm(relevel(as.factor(MSI_status),ref="Non-MSI-H") ~ relevel(as.factor(y),ref="0") + AGE + SEX + PRIMARY_RACE + CANCER_TYPE_DETAILED + SAMPLE_TYPE + SEQ_ASSAY_ID + tmb, data = Bsub, family="binomial")) #tmb

sum.coef<- Aout$coef
j<-j+1

est<-exp(sum.coef[2,1])
upper.ci<-exp(sum.coef[2,1]+1.96*sum.coef[2,2])
lower.ci<-exp(sum.coef[2,1]-1.96*sum.coef[2,2])

out.m[j,]<- c(g1,Aout$coeff[2,c(1,2,4)],est,lower.ci,upper.ci,frq,frq1,frq2,frq3,frq4,My1,My0,Mo1,Mo0,Mt1,Mt0,Ms1,Ms0)
}	

j
out.m1<-out.m[1:58,]
out.m1<-as.data.frame(out.m1)

names(out.m1) <- c("Gene","Beta","SE","P","OR","95CI1","95CI2","Freq","Freq_0","Freq_1","Freq_2","Freq_3","#_Carrie_0","#_0All","#_Carrie_1","#_1All","#_Carrie_2","#_2All","#_Carrie_3","#_3All")

out.m1$P<- as.numeric(as.character(out.m1$P))
out.m1$P_BF<-p.adjust(out.m1$P,method="bonferroni")

write.csv(out.m1, "/Users/kumc/Desktop/gene_raw_genie/coread_genie.csv")
