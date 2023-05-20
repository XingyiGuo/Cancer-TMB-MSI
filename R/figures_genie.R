######## MSI-TMB project #########
######## Data preparation ########
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)

##################################

# Figures 
# Figure 2A
library(lattice)
# Create data
gfg <- data.frame(x = c(31.2991388992887,	28.5660801198053, 18.324055965097,	17.6019256807582, 17.1195652173913,	11.9565217391304, 13.7959866220736,	37.1655518394649, 13.7333333333333,	11.4666666666667, 13.4594397962896,	4.25609312477264, 12.9672897196262,	31.0747663551402, 11.9642191576593,	3.65262765560939, 11.5107913669065,	23.021582733813, 10.6440071556351, 2.14669051878354, 10.2591792656587,	4.21166306695464, 9.8622763970239,	21.5450372012031, 9.39007092198582,	6.45390070921986, 8.8339222614841,	20.8480565371025, 8.72093023255814,	2.90697674418605, 8.22746521476104,	4.65819721718088, 8.09859154929578,	3.16901408450704, 7.97024442082891,	3.18809776833156, 7.53768844221105,	2.77859887673662, 6.96055684454756,	2.32018561484919, 6.38297872340426,	46.4095744680851, 5.51508844953174,	3.53798126951093, 12.5302245250432,	14.7754749568221), grp = rep(c('g01','g02','g03','g04','g05','g06','g07','g08','g09','g10','g11','g12','g13','g14','g15','g16','g17','g18','g19','g20','g21','g22','g23'), each = 2), subgroup = LETTERS[1:2])
gfg
# 'UCEC', 'COREAD', 'STAD', 'BLAD', 'ESCA', 'OV', 'LUSC', 'PRAD', 'CESC', 'RCC', 'CC', 'LUAD', 'BRCA', 'HNSC', 'GINE', 'SARC', 'LIHC', 'GBM', 'PAAD', 'MESO', 'MELA', 'THYR', 'Total'

# Create grouped barplot using lattice
barchart(x ~ grp, data = gfg, groups = subgroup)

#####################################################

# Figure 2B
# dot plot http://rstudio-pubs-static.s3.amazonaws.com/1810_f2741ac792984c7eb5803eff9fcfebf3.html
test1<-t[t$MSI_status=="MSI-H",]

x1  = factor(test1$CANCER_TYPE, levels=c('Endometrial Cancer', 'Colorectal Cancer', 'Stomach Cancer', 'Bladder Cancer', 'Esophageal Cancer', 'Ovarian Cancer', 'Lung Squamous Cell Carcinoma', 'Prostate Cancer', 'Cervical Cancer', 'Renal Cell Carcinoma', 'Cholangiocarcinoma', 'Lung Adenocarcinoma', 'Breast Cancer', 'Head and Neck Cancer', 'Gastrointestinal Neuroendocrine Tumor', 'Soft Tissue Sarcoma', 'Hepatocellular Carcinoma', 'Glioma', 'Pancreatic Cancer', 'Mesothelioma', 'Melanoma', 'Thyroid Cancer'))

p <- ggplot(test1, aes(x=x1, y=tmb, fill=MSI_status))
p + geom_point(shape=20, color="black", fill="black", position=position_jitter(width=.15, height=0)) + ylim(0,100) + geom_hline(yintercept=c(10), linetype="dashed", color = "red") + theme_classic() + theme(legend.position="none") + ylab("") + xlab("")

##################################################

# Figure 2C
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)

table(t$CANCER_TYPE)
t_sample<-t[t$CANCER_TYPE=="Prostate Cancer",]

# Plot https://r-graph-gallery.com/135-stacked-density-graph.html
t_sample %>%
  filter(MSI_status %in% c("MSI-H", "Non-MSI-H")) %>%
  ggplot( aes(x=tmb, color=MSI_status, fill=MSI_status)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values=c("red", "skyblue")) +
  ylab("") + xlab("") +
   xlim(0,100) + #ylim(0,0.2, by=0.05) +
  theme(panel.background = element_rect(fill='white', color='black', linetype='solid')) +
  theme(legend.position="none") + geom_vline(xintercept = 10) +
  scale_y_continuous(limits = c(0, 0.28), breaks = seq(0, 0.28, 0.05))

############################################

#Stacked + Percentage
# create a dataset
specie <- c(rep("G", 4) , rep("K", 4) , rep("T", 4))
condition <- rep(c("Non-MSI-H and TMB-low", "MSI-H and TMB-low", "Non-MSI-H and TMB-high", "MSI-H and TMB-high"), 3)
value <- c(1162,47,950,97,85,0,1,0,56,0,24,0)
data <- data.frame(specie,condition,value)

# Stacked + percent
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=c("firebrick1","gold","olivedrab3","lightskyblue1"))+
  theme(panel.background = element_rect(fill='white', color='black', linetype='solid'))+
  theme(legend.position="left") +
  ylab("") +
  xlab("") + theme(legend.position="none")

#EC: c(1640, 268, 195, 568, 57,	1, 7, 5, 153, 24, 18, 43)
#CRC: c(4956, 521, 473, 697, 1323, 6, 102, 56, 222, 1, 3, 38)
#STAD: c(586, 62, 24, 64, 539, 5, 31, 28, 192, 1, 4, 53)
#BLAD: c(1359, 144, 703, 186, 191,	2, 44, 4, 286, 0,	63,	0)
#ESCA: c(1199,129,95,77,95,0,4,0,153,0,0,3)

#OV: c(2324,308,55,62,279,0,14,3,317,0,0,0)
#LUSC: c(536,54,209,57,55,0,10,0,233,0,51,0)
#PRAD: c(2344,241,18,80,101,0,9,3,153,0,0,0)
#LUAD: c(4584,372,1110,251,171,0,35,1,302,0,111,0)
#BRCA: c(6017,578,371,84,678,0,45,4,754,0,18,0)

#SARC: c(1470,106,47,30,286,0,11,1,0,0,0,0)
#GBM: c(847,64,19,11,14,0,1,2,201,0,4,0)
#PAAD: c(3071,218,57,37,198,0,5,1,157,0,1,0)
#MELA: c(1162,47,950,97,85,0,1,0,56,0,24,0)

############################################

# Supplementary Figure 2
library(gcookbook)
library(dplyr)
data<-fread("/Users/kumc/Desktop/repel.csv",head=TRUE)

corr_sp <- ggplot(data, aes(x = Beta, y = log10P),pch = 21,fill = data$Cancer) +
  geom_point() +
  ylab("") +
  xlab("")+
  theme(panel.background = element_rect(fill='white', color='black', linetype='solid'))

library(ggrepel)
corr_sp + geom_text_repel(aes(label = Cancer), size = 4) + geom_hline(yintercept = 2.64, linetype="dashed", color = "red")
#corr_sp + geom_text(aes(label = Cancer), size = 3)
#corr_sp + geom_label_repel(aes(label = Cancer), size = 3)

est<-exp(sum.coef[2,1])
upper.ci<-exp(sum.coef[2,1]+1.96*sum.coef[2,2])
lower.ci<-exp(sum.coef[2,1]-1.96*sum.coef[2,2])

#https://stackoverflow.com/questions/58657802/forest-plot-with-subgroups-in-ggplot2
# Forest plot
df<-fread("/Users/kumc/Desktop/other_s3.csv",head=TRUE) 

Outcome_order <- c('KM_Asian_STAD','GENIE_Asian_STAD',
                   'KM_Asian_COREAD','GENIE_Asian_COREAD',
                   'KM_Asian_UCEc','GENIE_Asian_UCEC')
df$Outcome = factor (df$name, level=Outcome_order)

#define colours for dots and bars
dotCOLS = c("darkgray","blue","green")
barCOLS = c("black","skyblue","lightgreen")

p <- ggplot(df, aes(x=name2, y=`log(OR)`, ymin=`log(Lower)`, ymax=`log(Upper)`, col="black", fill="black")) + 
  #specify position here
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  #specify position here too
  geom_point(size=2, shape=21, colour="black", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS)+
  #scale_x_discrete(name="(Post)operative outcomes") +
  scale_y_continuous(name="Log OR (95% CI)", limits = c(0, 80)) +
  coord_flip() +
  theme_minimal()
p+
  ylab("") +
  xlab("")+
  theme(panel.background = element_rect(fill='white', color='black', linetype='solid'))+
  theme(legend.position="none")

############################################

# Figure 4A
# create a dataset
x1  = factor(test1$CANCER_TYPE, levels=c('Endometrial Cancer', 'Colorectal Cancer', 'Stomach Cancer', 'Bladder Cancer', 'Esophageal Cancer', 'Ovarian Cancer', 'Lung Squamous Cell Carcinoma', 'Prostate Cancer', 'Lung Adenocarcinoma', 'Breast Cancer', 'Soft Tissue Sarcoma', 'Glioma', 'Pancreatic Cancer', 'Melanoma'))

specie <- c(rep("g01",2), rep("g02",2), rep("g03",2), rep("g04",2),rep("g05",2), rep("g06",2), rep("g07",2), rep("g08",2),
            rep("g09",2), rep("g10",2), rep("g11",2), rep("g12",2),rep("g13",2), rep("g14",2))
condition <- rep(c("norminal P" , "Bonferroni-corrected P"), 14)
value <- c(29,41,20,38,1,8,12,20,1,5,1,6,0,6,
           2,10,5,12,3,15,0,9,1,3,2,8,2,10)
data <- data.frame(specie,condition,value)

# Stacked
p<-ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="stack", stat="identity") +
  ylab("") +
  xlab("")+
  theme(panel.background = element_rect(fill='white', color='black', linetype='solid'))+
  theme(legend.position="none")
  p+ scale_fill_manual(values =c("red", "lightskyblue1"))
#####################################################

# create a dataset
specie <- c(rep("01UCEC",2), rep("02COREAD",2), rep("03STAD",2), rep("04BLAD",2), rep("05ESCA",2), rep("06OV",2), rep("07LUSC",2),
            rep("08PRAD",2), rep("09LUAD",2), rep("10BRCA",2), rep("11SARC",2), rep("12GBM",2), rep("13PAAD",2), rep("14MELA",2))
condition <- rep(c("BF P" , "Nominal P"), 14)
value <- c(29,41,20,38,1,8,12,20,1,5,1,6,0,6,
           2,10,5,12,3,15,0,9,1,3,2,8,2,10)
data <- data.frame(specie,condition,value)
  
# Grouped
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="dodge", stat="identity") +
  ylab("") +
  xlab("")+
  theme(panel.background = element_rect(fill='white', color='black', linetype='solid'))+
  theme(legend.position="none") +
  scale_fill_manual(values =c("red", "lightskyblue1"))
#####################################################

# Figure 4B
#https://stackoverflow.com/questions/58657802/forest-plot-with-subgroups-in-ggplot2
# Forest plot
df<-fread("/Users/kumc/Desktop/other_s4.csv",head=TRUE) 
df<-df[!(df$Cancer=="Stomach"|df$Cancer=="Endometrium"|df$Cancer=="Colorectum"),]
Outcome_order <- c('MELA_IL7R',
                   'MELA_ATRX',
                   'PAAD_TSC2',
                   'PAAD_ARID1B',
                   'GBM_ARID1B',
                   'BRCA_PIK3CA',
                   'BRCA_RB1',
                   'BRCA_ARID1B',
                   'LUAD_SMARCA4',
                   'LUAD_ARID1A',
                   'LUAD_KEAP1',
                   'LUAD_ARID1B',
                   'LUAD_STK11',
                   'PRAD_FOXA1',
                   'PRAD_ARID1B',
                   'OV_ARID1B',
                   'ESCA_ARID1A',
                   'BLAD_ARID1A',
                   'BLAD_KMT2D',
                   'BLAD_ERCC2',
                   'BLAD_NOTCH1',
                   'BLAD_NSD1',
                   'BLAD_FBXW7',
                   'BLAD_ASXL1',
                   'BLAD_NOTCH3',
                   'BLAD_CREBBP',
                   'BLAD_ARID1B',
                   'BLAD_ZFHX3',
                   'BLAD_MSH2')
df$Outcome = factor (df$name2, level=Outcome_order) #Gene, name

#define colours for dots and bars
dotCOLS = c("gray","blue","green")
barCOLS = c("black","skyblue","lightgreen")

p <- ggplot(df, aes(x=Outcome, y=OR, ymin=Lower, ymax=Upper, col=group, fill=group)) + 
  #specify position here
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  #specify position here too
  geom_point(size=2, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS)+
  #scale_x_discrete(name="(Post)operative outcomes") +
  scale_y_continuous(name="OR (95% CI)", limits = c(0, 20)) +
  coord_flip() +
  theme_minimal()
p+
  ylab("") +
  xlab("")+
  theme(panel.background = element_rect(fill='white', color='black', linetype='solid'))+
  theme(legend.position="none")

#####################################################

# Bar graph
ucec<-fread("/Users/kumc/Desktop/stad.csv",head=TRUE) 
  ucec<-data.frame(ucec)
  ucec$log10P
library(ggplot2)
ggplot(ucec,aes(log10P,reorder(GENIE,+log10P),fill=Beta))+geom_bar(stat="identity")+
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0)+
  ylab("")+
  xlab("")+
  theme(panel.background = element_rect(fill='white', color='black', linetype='solid'))+
  theme(legend.position="none")+
  geom_col(size = 0.01, color = "black") 
#####################################################

# Figure 5
# Ben diagram
t<-fread("/Users/kumc/Desktop/sig_genes.csv",head=TRUE) 
ucec<-t[(t$Cancer=="Endometrium" & t$`BF P`<0.05),]
coread<-t[(t$Cancer=="Colorectum" & t$`BF P`<0.05),]
stad<-t[(t$Cancer=="Stomach"&t $ `BF P`<0.05),]

#1)
ucec$ucec_p<-ucec$`BF P`
coread$coread_p<-coread$`BF P`
stad$stad_p<-stad$`BF P`

ucec<-subset(ucec,select=c(Gene, ucec_p))
coread<-subset(coread,select=c(Gene, coread_p)) 
stad<-subset(stad,select=c(Gene, stad_p)) 

a<-merge(coread,stad,by="Gene")
b<-merge(ucec,stad,by="Gene")
c<-merge(ucec,coread,by="Gene")

b<-ucec[!(ucec$Gene %in% c$Gene),]
b<-coread[!(coread$Gene %in% c$Gene),]

########################################## The End May 10, 2023 #####################
#2)
ucect1<-ucec[!(ucec$Gene %in% coread$Gene),]
ucect2<-ucect1[!(ucect1$Gene %in% stad$Gene),]

coreadt1<-coread[!(coread$Gene %in% stad$Gene),]
coreadt2<-coreadt1[!(coreadt1$Gene %in% ucec$Gene),]

stadt1<-stad[!(stad$Gene %in% ucec$Gene),]
stadt2<-stadt1[!(stadt1$Gene %in% coread$Gene),]

#3)
ucec<-subset(ucec,select=c(Gene, Beta, P, cancer, OR, Lower, Upper))
ucec<-ucec[ucec$Gene %in% total$Gene,]
coread<-subset(coread,select=c(Gene, Beta, P, cancer, OR, Lower, Upper)) 
coread<-coread[coread$Gene %in% total$Gene,]
stad<-subset(stad,select=c(Gene, Beta, P, cancer, OR, Lower, Upper)) 
stad<-stad[stad$Gene %in% total$Gene,]

total<-rbind(ucec,coread,stad)
write.csv(total, "/Users/choij/Desktop/meta_sample.csv")

#4)
Cohort<-fread("/Users/choij/Desktop/bubble.csv",head=TRUE)

Cohort$stars <- cut(Cohort$P, breaks=c(-Inf, 0.0005, 0.05, Inf), label=c("", "*", "**"))  # Create column of significance labels

Outcome_order <- c('CREBBP', 'MSH6', 'KMT2C', 'MLH1','SMARCA4','NOTCH3','POLD1','ARID1B','ARID1A','MSH3','RNF43')
Cohort$Outcome = factor (Cohort$Gene, level=Outcome_order)

ggplot(Cohort, aes(x=cancer, y=Outcome, fill=Beta_modify)) + 
  geom_tile()+
  #scale_fill_gradient2(low="blue", mid="white", high="red",midpoint=0) +
  scale_fill_gradient2(low="blue", mid="white", high="red",
                       breaks=c(-3,0,3),labels=c("Minimum",0,"Maximum"),
                       limits=c(-3,3)) +
  geom_text(aes(label=stars), color="gray33", size=3)
#+
  ylab("") +
  xlab("")+
  theme(panel.background = element_rect(fill='white', color='black', linetype='solid'))+
  theme(legend.position="none")

#####################################################
ucec1<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_110622/dataset/ucec.csv",head=TRUE) 
  coread1<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_110622/dataset/coread.csv",head=TRUE) 
  stad1<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_110622/dataset/stad.csv",head=TRUE) 
  ov<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_110622/dataset/ov.csv",head=TRUE) 
  blad<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_110622/dataset/blad.csv",head=TRUE) 
  lusc<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_110622/dataset/lusc.csv",head=TRUE) 
  esca<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_110622/dataset/esca.csv",head=TRUE) 
  prad<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_110622/dataset/prad.csv",head=TRUE) 
  luad<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_110622/dataset/luad.csv",head=TRUE) 
  gbm<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_110622/dataset/gbm.csv",head=TRUE) 
  brca<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_110622/dataset/brca.csv",head=TRUE) 
  sarc<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_110622/dataset/sarc.csv",head=TRUE) 
  paad<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_110622/dataset/paad.csv",head=TRUE) 
  mela<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_110622/dataset/mela.csv",head=TRUE) 
  
  ucec<-subset(ucec,select=-c(V1))  

### p-value #############################  
ucec$ucec_p<-ucec$P
coread$coread_p<-coread$P
stad$stad_p<-stad$P
ov$ov_p<-ov$P
blad$blad_p<-blad$P

lusc$lusc_p<-lusc$P
esca$esca_p<-esca$P
prad$prad_p<-prad$P
luad$luad_p<-luad$P
gbm$gbm_p<-gbm$P

brca$brca_p<- brca$P
sarc$sarc_p<-sarc$P
paad$paad_p<-paad$P
mela$mela_p<-mela$P

ucec<-subset(ucec,select=c(Gene, ucec_p))
coread<-subset(coread,select=c(Gene, coread_p)) 
stad<-subset(stad,select=c(Gene, stad_p)) 
ov<-subset(ov,select=c(Gene, ov_p)) 
blad<-subset(blad,select=c(Gene, blad_p)) 

lusc<-subset(lusc,select=c(Gene, lusc_p))
esca<-subset(esca,select=c(Gene, esca_p)) 
prad<-subset(prad,select=c(Gene, prad_p)) 
luad<-subset(luad,select=c(Gene, luad_p)) 
gbm<-subset(gbm,select=c(Gene, gbm_p)) 

brca<-subset(brca,select=c(Gene, brca_p)) 
sarc<-subset(sarc,select=c(Gene, sarc_p)) 
paad<-subset(paad,select=c(Gene, paad_p)) 
mela<-subset(mela,select=c(Gene, mela_p)) 

c<-merge(c,mela,by="Gene")
write.csv(c, "/Users/choij/Desktop/c.csv")


### beta #############################
ucec$ucec_b<-ucec$Beta
coread$coread_b<-coread$Beta
stad$stad_b<-stad$Beta
ov$ov_b<-ov$Beta
blad$blad_b<-blad$Beta

lusc$lusc_b<-lusc$Beta
esca$esca_b<-esca$Beta
prad$prad_b<-prad$Beta
luad$luad_b<-luad$Beta
gbm$gbm_b<-gbm$Beta

brca$brca_b<- brca$Beta
sarc$sarc_b<-sarc$Beta
paad$paad_b<-paad$Beta
mela$mela_b<-mela$Beta

ucec<-subset(ucec,select=c(Gene, ucec_b))
coread<-subset(coread,select=c(Gene, coread_b)) 
stad<-subset(stad,select=c(Gene, stad_b)) 
ov<-subset(ov,select=c(Gene, ov_b)) 
blad<-subset(blad,select=c(Gene, blad_b)) 

lusc<-subset(lusc,select=c(Gene, lusc_b))
esca<-subset(esca,select=c(Gene, esca_b)) 
prad<-subset(prad,select=c(Gene, prad_b)) 
luad<-subset(luad,select=c(Gene, luad_b)) 
gbm<-subset(gbm,select=c(Gene, gbm_b)) 

brca<-subset(brca,select=c(Gene, brca_b)) 
sarc<-subset(sarc,select=c(Gene, sarc_b)) 
paad<-subset(paad,select=c(Gene, paad_b)) 
mela<-subset(mela,select=c(Gene, mela_b)) 

b<-merge(b,mela,by="Gene")
write.csv(b, "/Users/choij/Desktop/b.csv")
#####################################################

install.packages("gcookbook")
library(gcookbook)
library(dplyr)
data<-fread("/Users/choij/Dropbox/VUMC project/P11_KMVU/result_121822/bu.csv",head=TRUE)
data$`-LOG10(P)`<-data$`#NAME?`
corr_sp <- ggplot(data, aes(x = Beta, y = `-LOG10(P)`),pch = 21) +
  geom_point()+
  #ylab("") +
  #xlab("")+
  theme(panel.background = element_rect(fill='white', color='black', linetype='solid'))

library(ggrepel)
corr_sp +
  geom_text_repel(aes(label = Cancer), size = 3)
corr_sp +
  geom_text(aes(label = Cancer), size = 3)+
  geom_hline(yintercept=2.64, lty=2) 
corr_sp +
  geom_label_repel(aes(label = Cancer), size = 3)+
  geom_hline(yintercept=2.64, lty=2) 
###
corr_sp +
  geom_label_repel(aes(label = Name), size = 4)+
  scale_fill_manual(values = setNames(c("lightblue", "darkgreen", "yellow", "red"), levels(data$race)))

#https://r-graphics.org/recipe-scatter-labels



#####################################################
library(dendextend)
library(factoextra)
library(cluster)
require(data.table)
library(aod)
library('tidyverse')
library(ISLR)
library("survival")
library("survminer")
library(reshape)
library(dplyr)
library(broom)
library(ggplot2)
install.packages("hrbrthemes")
library(hrbrthemes)
library(viridis)
library(dplyr)
#Bubble plot
data<-fread("C:/Users/choij/Desktop/bubble.csv",head=TRUE)
data$P<- -log10(data$P)
data$P<- data$`""-LOG10(P)""`
ggplot(data, aes(x=cancer, y=Gene, size=P, color=beta)) +
  geom_point(alpha=0.8)+
  scale_size(range = c(5, 12), name="-Log 10 P")+
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0)+
  #scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A") +
  theme_ipsum() +
  theme(legend.position="bottom") 
  #ylab("-Log10 P") +
  #xlab("Effect size") 


library(ggcorrplot)
library("ggplot2")

Cohort <- read.csv("/Users/kumc/Dropbox/PC/Desktop/corr.csv", head=T) 

# Compute a correlation matrix
corr <- round(cor(Cohort), 10)
corr
head(corr[, 1:17])

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(Cohort)
p.mat
head(p.mat[, 1:17])

# Tile graph
head(Cohort)
Cohort<-data
Cohort$stars <- cut(Cohort$P, breaks=c(-Inf, 0.05, 0.0005, Inf), label=c("*", "**", ""))  # Create column of significance labels

ggplot(Cohort, aes(x=cancer, y=Outcome, fill=Beta_modify)) + 
  geom_tile()+
  #scale_fill_gradient2(low="blue", mid="white", high="red",midpoint=0) +
  scale_fill_gradient2(low="blue", mid="white", high="red",
                       breaks=c(-3,0,3),labels=c("Minimum",0,"Maximum"),
                       limits=c(-3,3)) +
  geom_text(aes(label=stars), color="gray33", size=3)

# geom_text(aes(label = round(coefficient, 3))) +

# Cohort$variable <- factor(Cohort$Breast.cancer, levels=levels(with(Cohort[Cohort$Breast.cancer=="Breast cancer",], reorder(Cohort$Breast.cancer, -Cohort$Breast.cancer))))  # Sort coefficients by value


Outcome_order <- c("MSH6", "MLH1", "POLD1", "MSH3", "CREBBP", "KMT2C","SMARCA4", "NOTCH3","ARID1B","ARID1A", "RNF43")
Cohort$Outcome = factor(Cohort$Gene, level=Outcome_order)

#######################
