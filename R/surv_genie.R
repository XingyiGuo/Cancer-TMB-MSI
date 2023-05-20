######## MSI-TMB project #########
######## Data preparation ########

library(data.table)
library(tidyverse)
library(dplyr)
library("survival")
library("survminer")

#address: /Users/kumc/Dropbox/MSI-TMB/Data/crc-surv.csv
#Colorectal
crc_group <- fread("/Users/kumc/Desktop/t_crc.csv",head=TRUE,stringsAsFactors=F) 
crc_surv <- fread("/Users/kumc/Desktop/crc-surv.csv",head=TRUE,stringsAsFactors=F) 

#Lung	
lung_group <- fread("/Users/kumc/Dropbox/MSI-TMB/Data/t_lung.csv",head=TRUE,stringsAsFactors=F) 
lung_surv <- fread("/Users/kumc/Dropbox/MSI-TMB/Data/lung-surv.csv",head=TRUE,stringsAsFactors=F)

#Data mining
lung_group <- subset(lung_group,select=-c(V1))
lung_surv <- subset(lung_surv,select=-c(V1))
lung_group <- lung_group[!duplicated(lung_group$Tumor_Sample_Barcode),]
lung_surv$Tumor_Sample_Barcode <- lung_surv$`Sample ID`
lung_surv <- lung_surv[!duplicated(lung_surv$Tumor_Sample_Barcode),]

table(lung_group$CANCER_TYPE,lung_group$CANCER_TYPE_DETAILED)	
tt <- lung_surv
tt$Tumor_Sample_Barcode <- tt$`Sample ID`
	
  nivolumab<-tt[(tt$`OS Status from start of Nivolumab`==0 | tt$`OS Status from start of Nivolumab`==1),]
	nivolumab$`OS Status from start of IO` <- nivolumab$`OS Status from start of Nivolumab`
	nivolumab$`OSMonths from start of IO`<-nivolumab$`OSMonths from start of Nivolumab`
	nivolumab$`PFS-I Status from start of IO` <- nivolumab$`PFS-I Status from start of Nivolumab`
	nivolumab$`PFS-IMonths from start of IO`<-nivolumab$`PFS-IMonths from start of Nivolumab`
	nivolumab$`PFS-M Status from start of IO` <- nivolumab$`PFS-M Status from start of Nivolumab`
	nivolumab$`PFS-MMonths from start of IO`<-nivolumab$`PFS-MMonths from start of Nivolumab`
	
	pembrolizumab<-tt[(tt$`OS Status from start of Pembrolizumab`==0 | tt$`OS Status from start of Pembrolizumab`==1),]
	pembrolizumab$`OS Status from start of IO` <- pembrolizumab$`OS Status from start of Pembrolizumab`
	pembrolizumab$`OSMonths from start of IO`<-pembrolizumab$`OSMonths from start of Pembrolizumab`
	pembrolizumab$`PFS-I Status from start of IO` <- pembrolizumab$`PFS-I Status from start of Pembrolizumab`
	pembrolizumab$`PFS-IMonths from start of IO`<-pembrolizumab$`PFS-IMonths from start of Pembrolizumab`
	pembrolizumab$`PFS-M Status from start of IO` <- pembrolizumab$`PFS-M Status from start of Pembrolizumab`
	pembrolizumab$`PFS-MMonths from start of IO`<-pembrolizumab$`PFS-MMonths from start of Pembrolizumab`
	
	atezolizumab<-tt[(tt$`OS Status from start of Atezolizumab`==0 | tt$`OS Status from start of Atezolizumab`==1),]
	atezolizumab$`OS Status from start of IO` <- atezolizumab$`OS Status from start of Atezolizumab`
	atezolizumab$`OSMonths from start of IO`<-atezolizumab$`OSMonths from start of Atezolizumab`
	atezolizumab$`PFS-I Status from start of IO` <- atezolizumab$`PFS-I Status from start of Atezolizumab`
	atezolizumab$`PFS-IMonths from start of IO`<-atezolizumab$`PFS-IMonths from start of Atezolizumab`
	atezolizumab$`PFS-M Status from start of IO` <- atezolizumab$`PFS-M Status from start of Atezolizumab`
	atezolizumab$`PFS-MMonths from start of IO`<-atezolizumab$`PFS-MMonths from start of Atezolizumab`
	
	IO<-rbind(nivolumab,pembrolizumab,atezolizumab)
	
	lung_merge<-merge(lung_group, IO, by="Tumor_Sample_Barcode")

#ID matched test		
	table(lung_merge$CANCER_TYPE)
	
	lung_merge<-merge(lung_group, tt, by="Tumor_Sample_Barcode")
	
	table(lung_merge$Center.x)
  table(tt$Center)	
  table(lung_group$Center)
  table(lung_surv$Center)
  tt<-lung_surv
  t<-lung_group
  tt_test<-tt[tt$Center=="VICC",]
  t_test<-t[t$Center=="VICC",]
  tt_test$`Sample ID`
  t_test$SAMPLE_ID
######## Survival analysis OS and PFS ########

#Lung	
lung<-lung_merge
	
	lung$smoke<-lung$`Non-Small Cell Lung Cancer: Cigarette Use at Time of Diagnosis`
	lung$stage <- lung$`Stage at Diagnosis`
	lung$PDL1 <- lung$`Positive PD-L1 Result at the Time of Sample Acquisition`
	lung$PDL1test <- lung$`PD-L1 Testing at Time of Sample Acquisition`
	#table(lung$PDL1test)
	#table(lung$PDL1test, lung$PDL1)
	
	lung$time<-lung$`OSMonths from start of IO`
	lung$status<-lung$`OS Status from start of IO`
	
	lung$time<-lung$`PFS-IMonths from start of IO`
	lung$status<-lung$`PFS-I Status from start of IO`

lung_yes<-lung[lung$PDL1=="Yes",]	
lung_no<-lung[lung$PDL1=="No",]	
lung_notest<-lung[lung$PDL1test=="No"]

lung<-rbind(lung_no,lung_yes)

#COXPH	
	res.cox <- coxph(Surv(time, status) ~ group + AGE + SEX + relevel(as.factor(PRIMARY_RACE),ref="White") + CANCER_TYPE + SAMPLE_TYPE + SEQ_ASSAY_ID + smoke, data = lung) #group MSI_status TMB_status

	res.cox
	summary(res.cox)

	#relevel(as.factor(TMB_status),ref="TMB-low")
	Aout <- summary(glm(as.factor(PDL1) ~ group + AGE + SEX + PRIMARY_RACE + CANCER_TYPE_DETAILED + SAMPLE_TYPE + SEQ_ASSAY_ID + smoke, data = lung, family="binomial")) 
	Aout
	
	#http://www.sthda.com/english/wiki/cox-proportional-hazards-model
	
	group_df <- with(lung,data.frame(group = c(0,1,2,3),
	                                 AGE=mean(lung$AGE),
	                                 SEX = 'Male', 
	                                 SAMPLE_TYPE="Primary",
	                                 SEQ_ASSAY_ID="MSK-IMPACT410",
	                                 CANCER_TYPE="Lung Adenocarcinoma",
	                                 PRIMARY_RACE="White",
	                                 smoke="Former user (quit >1 year)"))
	# group = c(0,1,2,3),
	group_df
	fit <- survfit(res.cox, newdata=group_df)
	
	ggsurvplot(fit, data=lung, conf.int = TRUE, legend.labs=c("group=0", "group=1", "group=2", "group=3"),
	           ggtheme = theme_minimal()) + ylab("") + xlab("") 
		
	 
  fit  
	
	# https://rdrr.io/cran/survminer/man/surv_median.html
	fit2 <- surv_fit(Surv(time, status) ~ group, data = lung)
	surv_median(fit2)
	# "group=0", "group=1", "group=2", "group=3"
	# "MSI_status=MSI-H", "MSI_status=Non-MSI-H"
	#	"TMB_status=TMB-high", "TMB_status=TMB-low"
	
	#calculate number at risk
	lung$time<-lung$`PFS-MMonths from start of IO`
	lung$status<-lung$`PFS-M Status from start of IO`
	
	km_trt_fit <- survfit(Surv(time, status) ~ TMB_status, data=lung)
	summary(km_trt_fit, times = c(0,10,20,30,40,50,60))
	#************************************************************#
	#************************************************************#
	