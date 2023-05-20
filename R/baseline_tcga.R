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

cliM<-fread("/Users/kumc/Dropbox/MSI-TMB/Data/TCGA_raw_genotype/PAAD_TCGA.csv",head=TRUE,stringsAsFactors=F) #COREAD, SKCM, KIRC, KIRP ####
cliM$Tumor_Sample_Barcode<-substr(cliM$Tumor_Sample_Barcode,1,12) 

cliM<-merge(cliM,a,by="Tumor_Sample_Barcode")
cliM<-merge(cliM,b,by="Tumor_Sample_Barcode")

cliS<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/paad_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE) ####COADREAD
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

tt<-table(total$Cancer_type,total$group)
tt
table(total$`Person Gender`)

#### Demographics #############################
table(t$CANCER_TYPE)
mean(t$AGE)
sd(t$AGE)
table(t$PRIMARY_RACE)

#### Total AGE mean and SD ####################

cliM1<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/BLAD_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM2<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/BRCA_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM3<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/CESC_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM4<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/COREAD_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM4<-subset(cliM4,select=-c(V1))
cliM5<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/UCEC_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM6<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/ESCA_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM7<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/GBM_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM8<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/HNSC_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM9<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/LIHC_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM10<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/LUAD_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM11<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/LUSC_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM12<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/SKCM_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM13<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/OV_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM14<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/PAAD_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM15<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/PRAD_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM16<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/KIRC_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM16.2<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/KIRP_TCGA.csv",head=TRUE,stringsAsFactors=F)
cliM17<-fread("/Users/kumc/Desktop/MSI-TMB/Data/TCGA_raw_genotype/STAD_TCGA.csv",head=TRUE,stringsAsFactors=F)

cliM<-rbind(cliM1,cliM2,cliM3,cliM4,cliM5,cliM6,cliM7,cliM8,cliM9,cliM10,
            cliM11,cliM12,cliM13,cliM14,cliM15,cliM16,cliM16.2,cliM17)
cliM$Tumor_Sample_Barcode<-substr(cliM$Tumor_Sample_Barcode,1,12) 

cliM<-merge(cliM,a,by="Tumor_Sample_Barcode")
cliM<-merge(cliM,b,by="Tumor_Sample_Barcode")

cliS1<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/blca_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F) 
cliS2<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/brca_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F) 
cliS3<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/cesc_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F) 
cliS4<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/coad_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE,stringsAsFactors=F)
cliS4.2<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/read_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F)
cliS5<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/ucec_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F) 
cliS6<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/esca_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F) 
cliS7<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/gbm_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F) 
cliS8<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/hnsc_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F) 
cliS9<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/lihc_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F) 
cliS10<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/luad_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F) 
cliS11<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/lusc_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F)
cliS12<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/skcm_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE,stringsAsFactors=F) 
cliS13<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/ov_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F) 
cliS14<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/paad_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F) 
cliS15<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/prad_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F) 
cliS16<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/kirc_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F) 
cliS16.2<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/kirp_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F)
cliS17<-fread("/Users/kumc/Dropbox/Immune response/TCGAcli/stad_tcga_pan_can_atlas_2018_clinical_data.tsv",head=TRUE,stringsAsFactors=F)

cliS<-rbind(cliS1,cliS2,cliS3,cliS4,cliS4.2,cliS5,cliS6,cliS7,cliS8,cliS9,cliS10,
            cliS11,cliS12,cliS13,cliS14,cliS15,cliS16,cliS16.2,cliS17)
cliS<-cliS[!duplicated(cliS$`Patient ID`),] ###########
cliS$AGE <- cliS$`Diagnosis Age`
table(cliS$AGE)
cliS<-cliS[!(cliS$AGE %in% NA),]
cliS$AGE <- as.numeric(as.character(cliS$AGE))

cliS$PRIMARY_RACE <- cliS$`Race Category`
cliS$PRIMARY_RACE <- ifelse(cliS$PRIMARY_RACE %in% "Asian","Asian",ifelse(cliS$PRIMARY_RACE %in% "Black or African American","BLack",ifelse(cliS$PRIMARY_RACE %in% "White","White","Others")))
cliS<-cliS[!(cliS$PRIMARY_RACE %in% 'Others'),]

table(cliS$`Tissue Source Site`)
table(cliS$`American Joint Committee on Cancer Metastasis Stage Code`)
table(cliS$`American Joint Committee on Cancer Tumor Stage Code`)
table(cliS$`Tumor Disease Anatomic Site`)
table(cliS$`Tumor Type`)

table(cliS$SAMPLE_TYPE)
cliS1<-cliS
cliS1a<-cliS1[(cliS1$`American Joint Committee on Cancer Metastasis Stage Code`=="CM0 (I+)"|cliS1$`American Joint Committee on Cancer Metastasis Stage Code`=="M0"),]
cliS1a$SAMPLE_TYPE<-"Primary"
cliS1b<-cliS1[(cliS1$`American Joint Committee on Cancer Metastasis Stage Code`=="M1"|cliS1$`American Joint Committee on Cancer Metastasis Stage Code`=="M1A"|cliS1$`American Joint Committee on Cancer Metastasis Stage Code`=="M1B"|cliS1$`American Joint Committee on Cancer Metastasis Stage Code`=="M1C"),]
cliS1b$SAMPLE_TYPE<-"Metastasis"
cliS1<-rbind(cliS1a,cliS1b)
cliS$SAMPLE_TYPE <- ifelse(cliS$SAMPLE_TYPE %in% "Primary","Primary",ifelse(cliS$SAMPLE_TYPE %in% "Metastasis","Metastasis","Others"))
cliS<-cliS[!(cliS$SAMPLE_TYPE %in% 'Others'),]

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
A.base1 <- A.base1[!duplicated(A.base1$Tumor_Sample_Barcode),] #n=5118

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

total <- A.base1
t<-total[!duplicated(total$Tumor_Sample_Barcode),]
t<-total[!duplicated(total$Tumor_Sample_Barcode),]

tt<-table(total$Cancer_type,total$group)
tt
table(total$`Person Gender`)

mean(total$AGE)
sd(total$AGE)
