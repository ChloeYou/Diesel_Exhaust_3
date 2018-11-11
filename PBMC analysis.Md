---
title: "DE3-PBMC Analysis"
author: "cyou"
date: "6/27/2018"
output: html_document
---
DE3: PBMC Analysis 
=====================
# some functions used 
```{r}
Mval<-function(beta) log2(beta/(1-beta))
```

PBMC have 108 samples(excludes the missing samples) collected from 9 individuals.
each exposure corresponds to three time point: before exposure, 24 hours after expsoure and 48 hours after expsoure 

**to potentionally simplify things, we will compare the methylation profile of the 24 hr and 48 hr samples first**
# compare methylation profile of the 24 hr and 48 hr samples using the invariable probes filtered betas
```{r}
load("DE3meta.RData") # filtermeta
load("de3_bmiq_PBMC_cell_corrected.RData") #de3_bmiq_PBMC
load("de3_PBMC_filtered_invariable_probes.RData") #de3_PBMC_betas

## extract PBMC from filtermeta 
meta_PBMC<-filtermeta[filtermeta$SampleType=="PBMC",]

## create M value file from betas 
mval_PBMC<-Mval(de3_PBMC_betas)

## check data and save
dim(mval_PBMC) # 689498    108
stopifnot(!is.na(de3_PBMC_betas)) # no NAs
stopifnot(!is.na(mval_PBMC))
save(mval_PBMC,file="mval_PBMC.RData")

## set filtered air-salin as control group
meta_PBMC$Exposure<-relevel(meta_PBMC$Exposure,ref = "FA-S")

## check sample order
identical(as.character(colnames(de3_PBMC_betas)),as.character(meta_PBMC$SampleID)) # TRUE 

## extract samples collected at 24 and 48 hrs
meta_PBMC_subset<-meta_PBMC[meta_PBMC$SampleTimeHours%in%c("24","48"),]  ## 71 samples

## one pair of replicates in 24 and 48 hrs, check correlation 310-17  
### get PBMC betas for the replicates
checkRep<- betas(de3_bmiq_PBMC)[,grep("rep",colnames(betas(de3_bmiq_PBMC)))]
cor <- cor.test(checkRep[,"310-17_rep1"], checkRep[,"310-17_rep2"])  #0.9957836 
cor <- cor.test(checkRep[,"317-18_rep1"], checkRep[,"317-18_rep2"])  #0.9982665
### conclusion: highly correlated, seems fine

## extract Mval for 24 hrs and 48 hrs for lmer later 
mval_PBMC_subset<-mval_PBMC[,colnames(mval_PBMC)%in%meta_PBMC_subset$SampleID]
dim(mval_PBMC_subset) # 689498     71

## since 310-17 replicate seems fine, we will just randomly pick one and remove the other from meta and mval
meta_PBMC_subset<-meta_PBMC_subset[meta_PBMC_subset$SampleID!="310-17_rep2",] # 70 31
mval_PBMC_subset<-mval_PBMC_subset[,colnames(mval_PBMC_subset)!="310-17_rep2"] # 689498     70
identical(as.character(colnames(mval_PBMC_subset)),as.character(meta_PBMC_subset$SampleID)) # True
### change 310-17 sample name to 310-17
meta_PBMC_subset$SampleID<-as.character(meta_PBMC_subset$SampleID)
meta_PBMC_subset[2,1]<-"310-17"
colnames(mval_PBMC_subset)[2]<-"310-17"

## remove 310-17_rep2 from beta file as well
## DE3-BRUSH analysis.Rmd line 56
de3_PBMC_betas<-as.data.frame(de3_PBMC_betas)
colnames(de3_PBMC_betas)[which(names(de3_PBMC_betas) == "310-17_rep1")] <- "310-17"

## pair data 
### 347-14 and 355-14 were not paired. 347-14 was identified as outlier and 355-13 was not available according to gDNA chart
### we will remove 347-14 and 355-14 from meta and mval 
meta_PBMC_subset<-meta_PBMC_subset[!meta_PBMC_subset$SampleID%in%c("347-14","355-14"),] # 68 31
mval_PBMC_subset<-mval_PBMC_subset[,!colnames(mval_PBMC_subset)%in%c("347-14","355-14")] # 689498     68

## overview of the data
as.factor(meta_PBMC_subset$VolunteerNumber) #9 individuals:  310 315 317 318 347 355 361 365 368
### there should be 9*4(exposures)*2(two timepoints)=72 samples (matches the two samples we removed earlier because of missing paired data. Perfect!)

save(meta_PBMC_subset,mval_PBMC_subset,file = "PBMC_compare_timepoints.RData")

# null hypothesis: DNAm profile at 24hrs and 48 hrs after expsoure are the same within the same individual  

## idea: create a new variable that assigns 34 pairs 
## add new column to the meta 
library(xlsx)
write.csv(meta_PBMC_subset, "meta_PBMC_subset.csv")

## did this in excel forgive me god of RStudio
meta_PBMC_subset<-read.csv("meta_PBMC_subset.csv",header =TRUE,sep = ",")
meta_PBMC_subset<-meta_PBMC_subset[,colnames(meta_PBMC_subset)!="X"] #68 18
### removed some column

## reorder meta SampleID
meta_PBMC_subset<-meta_PBMC_subset[order(match(meta_PBMC_subset$SampleID,colnames(mval_PBMC_subset))),]

save(meta_PBMC_subset,mval_PBMC_subset,file = "PBMC_compare_timepoints_2.RData")

## LINEAR MODEL
PBMC_tp_lmer<-pbsapply(1:nrow(mval_PBMC_subset),function(CpG) {
  metaex<-meta_PBMC_subset
  metaex$Mval<-mval_PBMC_subset[CpG,]
  mod<-lmer(Mval ~ as.factor(SampleTimeHours)+Exposure+(1|VolunteerNumber),data=metaex)
  mod_2<-lmer(Mval ~ Exposure+(1|VolunteerNumber),data=metaex)
  anova(mod,mod_2)[2,8]
    }
  )

PBMC_tp_lmer<-as.data.frame(PBMC_tp_lmer)
PBMC_tp_lmer$CpG<-rownames(mval_PBMC_subset)
save(PBMC_tp_lmer,file = "PBMC_two_timepoints_compare.RData")

## plot p-value distribution 
ggplot(PBMC_tp_lmer, aes(PBMC_tp_lmer))+geom_histogram(fill="grey90", color="black")+xlab("Nominal P Value")+xlim(min(PBMC_tp_lmer$PBMC_tp_lmer),max(PBMC_tp_lmer$PBMC_tp_lmer))+theme_minimal()


### try this one just for fun!
### running this in terminal rn
 PBMC_tp_lmer_2<-pbsapply(1:nrow(mval_PBMC_subset),function(CpG) {
  metaex<-meta_PBMC_subset
  metaex$Mval<-mval_PBMC_subset[CpG,]
  mod<-lmer(Mval ~ as.factor(SampleTimeHours)+(1|Exposure)+(1|VolunteerNumber),data=metaex)
  mod_2<-lmer(Mval ~ (1|Exposure)+(1|VolunteerNumber),data=metaex)
  anova(mod,mod_2)[2,8]
    }
  )
PBMC_tp_lmer_2<-as.data.frame(PBMC_tp_lmer_2)
PBMC_tp_lmer_2$CpG<-rownames(mval_PBMC_subset)
save(PBMC_tp_lmer_2,file = "PBMC_two_timepoints_compare_2.RData")
## plot p-value distribution 
ggplot(PBMC_tp_lmer_2, aes(PBMC_tp_lmer_2))+geom_histogram(fill="grey90", color="black")+xlab("Nominal P Value")+xlim(min(PBMC_tp_lmer_2$PBMC_tp_lmer_2),max(PBMC_tp_lmer_2$PBMC_tp_lmer_2))+theme_minimal()
```

## use correlation and root mean square error to examine the two timepoints
### 1. correlation 
```{r}
## construct two matrices for beta values, one for 24 hrs and one for 48 hrs, for each matrix, the rows are CpG sites and the columns are SampleID
## construct the 24 hrs matrix
str(meta_PBMC_subset$SampleTimeHours)
meta_PBMC_24<-meta_PBMC_subset[meta_PBMC_subset$SampleTimeHours==24,]
meta_PBMC_48<-meta_PBMC_subset[meta_PBMC_subset$SampleTimeHours==48,]
## should match the sample ID in the meta data 
meta_PBMC_24<-meta_PBMC_24[order(meta_PBMC_24$SampleID),]
meta_PBMC_48<-meta_PBMC_48[order(meta_PBMC_48$SampleID),]


head(colnames(de3_PBMC_betas))
beta_pbmc_24<-de3_PBMC_betas[,colnames(de3_PBMC_betas)%in%meta_PBMC_24$SampleID]
beta_pbmc_24<-beta_pbmc_24[,order(match(colnames(beta_pbmc_24),meta_PBMC_24$SampleID))]
beta_pbmc_48<-de3_PBMC_betas[,colnames(de3_PBMC_betas)%in%meta_PBMC_48$SampleID]
beta_pbmc_48<-beta_pbmc_48[,order(match(colnames(beta_pbmc_48),meta_PBMC_48$SampleID))]
## test correlation 
mapply(cor,as.data.frame(beta_pbmc_48),as.data.frame(beta_pbmc_24)) ## 34 correlation coefficients, each representing a pair of 24/48 hr column(each column has 689498 rows/entriese)

```
### 2. root mean square error
```{r}
library(Metrics)
mapply(rmse,as.data.frame(beta_pbmc_48),as.data.frame(beta_pbmc_24)) ## within individual 
```
### use correlation and rmse on 0 hr and 24 hr
```{r}
meta_PBMC_0<-meta_PBMC[meta_PBMC$SampleTimeHours==-4,]
meta_PBMC_0<-meta_PBMC_0[order(meta_PBMC_0$SampleLabel),]

## check replication 317-18_rep1 & 317-18_rep2
checkRep<- betas(de3_bmiq_PBMC)[,grep("rep",colnames(betas(de3_bmiq_PBMC)))]
cor <- cor.test(checkRep[,"317-18_rep1"], checkRep[,"317-18_rep2"])  #0.9982665

## remove rep from meta, beta, mval
### we will remove 317-18_rep2
meta_PBMC_0<-meta_PBMC_0[!meta_PBMC_0$SampleID=="317-18_rep2",]

### and then change the name of 317-18_rep1 into 317-18
meta_PBMC_0$SampleID<-as.character(meta_PBMC_0$SampleID)
meta_PBMC_0$SampleID[12]<-"317-18"

### after investigation, 347-12 and 355-12 are the corresponding samples don't have a match in the 24 hr and 48  hr file 
meta_PBMC_0<-meta_PBMC_0[!meta_PBMC_0$SampleID%in%c("347-12","355-12"),]

## beta 
colnames(de3_PBMC_betas)[50]<-"317-18"
beta_pbmc_0<-de3_PBMC_betas[,colnames(de3_PBMC_betas)%in%meta_PBMC_0$SampleID]
beta_pbmc_0<-beta_pbmc_0[,order(match(colnames(beta_pbmc_0),meta_PBMC_0$SampleID))]

## check match
colnames(beta_pbmc_0)
colnames(beta_pbmc_24)

mapply(cor,as.data.frame(beta_pbmc_0),as.data.frame(beta_pbmc_24))
mapply(rmse,as.data.frame(beta_pbmc_0),as.data.frame(beta_pbmc_24))
```

## use correlation and rmse on 0 hr and 48 hrs
```{r}
mapply(cor,as.data.frame(beta_pbmc_0),as.data.frame(beta_pbmc_48))
mapply(rmse,as.data.frame(beta_pbmc_0),as.data.frame(beta_pbmc_48))
```
## make a summary for all these values 
```{r}
summary(mapply(cor,as.data.frame(beta_pbmc_0),as.data.frame(beta_pbmc_48)))
summary(mapply(cor,as.data.frame(beta_pbmc_0),as.data.frame(beta_pbmc_24)))
summary(mapply(cor,as.data.frame(beta_pbmc_48),as.data.frame(beta_pbmc_24)))

summary(mapply(rmse,as.data.frame(beta_pbmc_0),as.data.frame(beta_pbmc_48)))
summary(mapply(rmse,as.data.frame(beta_pbmc_0),as.data.frame(beta_pbmc_24)))
summary(mapply(rmse,as.data.frame(beta_pbmc_48),as.data.frame(beta_pbmc_24)))
```

# Split the timepoint data and lmer
## 1. PBMC timepoint: 0hr and 24 hrs
```{r}
## SampleTimeHour will be tested first 
## merged meta file of 0 hour and 24 hours 
meta_pbmc_0_24<-rbind(meta_PBMC_24,meta_PBMC_0)

## merged beta file and then convert to Mval 
beta_pbmc_0_24<-cbind(beta_pbmc_0,beta_pbmc_24)
mval_pbmc_0_24<-Mval(beta_pbmc_0_24)
save(mval_pbmc_0_24,meta_pbmc_0_24,file="pbmc_first_lm.RData")

mval_pbmc_0_24<-as.matrix(mval_pbmc_0_24)

## already running in terminal: 2018-07-12
lm_0_24<-sapply(1:nrow(mval_pbmc_0_24),function(CpG) {
  metaex<-meta_pbmc_0_24
  metaex$Mval<-mval_pbmc_0_24[CpG,]
  mod_pbmc<-lmer(Mval ~ Exposure+(1|VolunteerNumber)+as.factor(SampleTimeHours),data=metaex)
  mod_pbmc_2<-lmer(Mval ~ Exposure+(1|VolunteerNumber),data=metaex)
  anova(mod_pbmc,mod_pbmc_2)[2,8]
    }
  )
## saved as "pbmc_first_lme.RData"

## create volcano plot thingzzzz
## create effect size file (beta)
## use beta_pbmc_0_24
dim(beta_pbmc_0_24) #689498     68
## creating a paired dataset to minus the values within pairs
colnames(beta_pbmc_0_24)
meta_pbmc_0_24$SampleID
## yees it matches 

## create two matrix, one for each time point 
## matrix for 0 hr, rownames are SampleID/VolunteerNumber, colnames are CpG sites(will be reused later)

colnames(beta_pbmc_0)
colnames(beta_pbmc_24)

deltabeta_0_24<- as.matrix(beta_pbmc_24)-as.matrix(beta_pbmc_0)
dim(deltabeta_0_24) #689498     34

delbet<-rowMeans(deltabeta_0_24)
delbet<-as.data.frame(delbet)
#### VOLCANO PLOT
head(rownames(delbet))
CpGsites<-rownames(delbet)
CpGsites<-as.data.frame(CpGsites)
lm_0_24<-as.data.frame(lm_0_24)
rownames(lm_0_24)<-CpGsites$CpGsites

makeVolcano(lm_0_24$lm_0_24,delbet$delbet,CpGsites$CpGsites,0.01,0.0001,"PBMC 24hr-0hr Volcano Plot")

#multiple test correction
Multi_test_corr_relaxed<-p.adjust(lm_0_24$lm_0_24, method = "fdr", n = length(lm_0_24$lm_0_24))

summary(Multi_test_corr_relaxed)
# min is 0.9979, nothing passed the FDR threshold :'(
```

## interaction plot: exposure, SampleTime DNAm
  
```{r}
## pick random probes to test interaction
metaex<-meta_pbmc_0_24
metaex$Mval<-mval_pbmc_0_24[239056,]
interaction.plot(metaex$Exposure, metaex$SampleTimeHours,metaex$Mval)
  
```

## 2. PBMC timepoint: 0hr and 48 hrs
```{r}
# run the lmm
## SampleTimeHour will be tested first 
## merged meta file of 0 hour and 48 hours 
meta_pbmc_0_48<-rbind(meta_PBMC_48,meta_PBMC_0)
## merged beta file and then convert to Mval 
beta_pbmc_0_48<-cbind(beta_pbmc_0,beta_pbmc_48)
mval_pbmc_0_48<-Mval(beta_pbmc_0_48)
save(mval_pbmc_0_48,meta_pbmc_0_48,file="pbmc_second_lm.RData")

mval_pbmc_0_48<-as.matrix(mval_pbmc_0_48)
## already running in terminal: 2018-07-12
lm_0_48<-sapply(1:nrow(mval_pbmc_0_48),function(CpG) {
  metaex<-meta_pbmc_0_48
  metaex$Mval<-mval_pbmc_0_48[CpG,]
  mod_pbmc<-lmer(Mval ~ Exposure+(1|VolunteerNumber)+as.factor(SampleTimeHours),data=metaex)
  mod_pbmc_2<-lmer(Mval ~ Exposure+(1|VolunteerNumber),data=metaex)
  anova(mod_pbmc,mod_pbmc_2)[2,8]
    }
  )
save(lm_0_48,file="pbmc_second_lme.RData")
lm_0_48<-as.data.frame(lm_0_48)
lm_0_48$CpG<-rownames(beta_pbmc_48)
ggplot(lm_0_48, aes(lm_0_48))+geom_histogram(fill="grey90", color="black")+xlab("Nominal P Value")+xlim(min(lm_0_48$lm_0_48),max(lm_0_48$lm_0_48))+theme_minimal() 
## saved as pval_0_48_original
## saved as "pbmc_second_lme.RData" 

## method two
lm_0_48<-sapply(1:nrow(mval_pbmc_0_48),function(CpG) { 
  metaex<-meta_pbmc_0_48
  metaex$Mval<-mval_pbmc_0_48[CpG,]
  Model.REML = lmer(Mval ~ Exposure+as.factor(SampleTimeHours)+(1|VolunteerNumber), REML = TRUE, data = metaex)
  anova(Model.REML)[2,6] # pval for SampleTimeHours
    }
  )
lm_0_48_lmerTest<-lm_0_48
save(lm_0_48_lmerTest,file="pbmc_0_48_lmerTest_version.RData")
## terminal

## compare lm_0_48 and lm_0_48_lmerTest!!!
## plot p-value distribution 

lm_0_48_lmerTest<-as.data.frame(lm_0_48_lmerTest)
lm_0_48_lmerTest$CpG<-rownames(beta_pbmc_48)
ggplot(lm_0_48_lmerTest, aes(lm_0_48_lmerTest))+geom_histogram(fill="grey90", color="black")+xlab("Nominal P Value")+xlim(min(lm_0_48_lmerTest$lm_0_48_lmerTest),max(lm_0_48_lmerTest$lm_0_48_lmerTest))+theme_minimal() 
## saved as pval_0_48_lmerTest 

## create volcano plot thingzzzz
## create effect size file (beta)
## use beta_pbmc_0_48
dim(beta_pbmc_0_48) #689498     68
## creating a paired dataset to minus the values within pairs
colnames(beta_pbmc_0_48)
meta_pbmc_0_48$SampleID
## yees it matches 

## create two matrix, one for each time point 
## matrix for 0 hr, rownames are SampleID/VolunteerNumber, colnames are CpG sites(will be reused later)

colnames(beta_pbmc_0)
colnames(beta_pbmc_48)

deltabeta_0_48<- as.matrix(beta_pbmc_48)-as.matrix(beta_pbmc_0)
dim(deltabeta_0_48) #689498     34

delbet<-rowMeans(deltabeta_0_48)
delbet<-as.data.frame(delbet)
#### VOLCANO PLOT
head(rownames(delbet))
CpGsites<-rownames(delbet)
CpGsites<-as.data.frame(CpGsites)
lm_0_48<-as.data.frame(lm_0_48)
rownames(lm_0_48)<-CpGsites$CpGsites

makeVolcano(lm_0_48$lm_0_48,delbet$delbet,CpGsites$CpGsites,0.01,0.00001,"PBMC 48hr-0hr Volcano Plot")

#multiple test correction
Multi_test_corr_relaxed<-p.adjust(lm_0_48$lm_0_48, method = "fdr", n = length(lm_0_48$lm_0_48))

summary(Multi_test_corr_relaxed)
# min is 0.002749
```

## PBMC next step:
### trying linear mixed effect alternatives
 cross check the hits from two groups   
```{r}
str(meta_PBMC$SampleTimeHours)# num
meta_PBMC$SampleTimeHours<-as.factor(meta_PBMC$SampleTimeHours)
str(meta_PBMC$SampleTimeHours) # factors w/ 3 levels
## run in terminal 
save(meta_pbmc_0_24,mval_pbmc_0_48,mval_pbmc_0_24,meta_pbmc_0_48,file="terminal_pbmc_exposure.RData")
library(lmerTest)
lm_0_24<-sapply(1:nrow(mval_pbmc_0_24),function(CpG) { 
  metaex<-meta_pbmc_0_24
  metaex$Mval<-mval_pbmc_0_24[CpG,]
  Model.REML = lmer(Mval ~ Exposure*SampleTimeHours+(1|VolunteerNumber), REML = TRUE, data = metaex)
  anova(Model.REML)[1,6] # pval for Exposure
    }
  )
## saved as pbmc_0_24_lmer_ver0821.RData
lm_0_24<-as.data.frame(lm_0_24)
lm_0_24$CpG<-rownames(beta_pbmc_24)

lm_0_48<-sapply(1:nrow(mval_pbmc_0_48),function(CpG) { 
  metaex<-meta_pbmc_0_48
  metaex$Mval<-mval_pbmc_0_48[CpG,]
  Model.REML = lmer(Mval ~ Exposure*SampleTimeHours+(1|VolunteerNumber), REML = TRUE, data = metaex)
  anova(Model.REML)[1,6] # pval for Exposure 
    }
  )
## saved as pbmc_0_48_lmer_ver0821.RData
lm_0_48<-as.data.frame(lm_0_48)
lm_0_48$CpG<-rownames(beta_pbmc_48)

#  anova(Model.REML)[1,6] # pval for Exposure 
#  anova(Model.REML)[2,6] # pval for SampleTimeHours 
#  anova(Model.REML)[3,6] # pval for interaction 
```
## next step: check if the small pvals have overlaps 
```{r}
summary(lm_0_24$lm_0_24)
summary(lm_0_48$lm_0_48)
## summary of nominal pval doesn't seem too bad, maybe FDR?
#multiple test correction
Multi_test_corr_relaxed_0_24<-p.adjust(lm_0_24$lm_0_24, method = "fdr", n = length(lm_0_24$lm_0_24))
summary(Multi_test_corr_relaxed_0_24) ## report this to mike
stat_hits<-as.data.frame(beta_pbmc_0_24)[which(Multi_test_corr_relaxed_0_24<=0.35),] # 1 at fdr 0.05; 1 at fdr 0.1; 2 at fdr 0.3, 3 at fdr 0.35

Multi_test_corr_relaxed_0_48<-p.adjust(lm_0_48$lm_0_48, method = "fdr", n = length(lm_0_48$lm_0_48))
summary(Multi_test_corr_relaxed_0_48) ## report this to mike
stat_hits<-as.data.frame(beta_pbmc_0_48)[which(Multi_test_corr_relaxed_0_48<=0.35),] # NOTHING SURVIVED FDR

## hypothesis: you can see some chanes after 24 hrs of exposure, but not much after 48 hrs
```
## next step: create effect size and check if that has overlaps 
```{r}
## deltabeta_0_24 and deltabeta_0_48 have already been created 
deltabeta_0_24<-as.data.frame(deltabeta_0_24)
head(rownames(deltabeta_0_24)) # it's a go
delbet_0_24<-rowMeans(deltabeta_0_24)
delbet_0_24<-as.data.frame(delbet_0_24)

stat_hits<-as.data.frame(beta_pbmc_0_24)[which(abs(delbet_0_24$delbet_0_24)>=0.05),] # this gets the list of CpG names, 18 
stat_hits_2<-as.data.frame(delbet_0_24)[which(abs(delbet_0_24$delbet_0_24)>=0.05),] # this gets the effect size corresponding to each CpG hit



deltabeta_0_48<-as.data.frame(deltabeta_0_48)
delbet_0_48<-rowMeans(deltabeta_0_48)
delbet_0_48<-as.data.frame(delbet_0_48)

stat_hits_3<-as.data.frame(beta_pbmc_0_48)[which(abs(delbet_0_48$delbet_0_48)>=0.05),] # this gets the list of CpG names, 6
stat_hits_4<-as.data.frame(delbet_0_48)[which(abs(delbet_0_48$delbet_0_48)>=0.05),] # this gets the effect size corresponding to each CpG hit

length(intersect(rownames(stat_hits_3),rownames(stat_hits))) # 0

```

## check for gene name overlaps in effect size
```{r}
stat_hits_0_24_CpG<-as.data.frame(rownames(stat_hits))
stat_hits_0_48_CpG<-as.data.frame(rownames(stat_hits_3))

fd_pbmc<-fData(de3_bmiq_PBMC)
fd_pbmc_0_24<-fd_pbmc[rownames(fd_pbmc)%in%stat_hits_0_24_CpG$`rownames(stat_hits)`,]
fd_pbmc_0_48<-fd_pbmc[rownames(fd_pbmc)%in%stat_hits_0_48_CpG$`rownames(stat_hits_3)`,]
###
fd_de3<-fData(de3_bmiq)
fd_de3_subset<-fd_de3[rownames(fd_de3)=="cg13808185",]
fd_de3_subset_2<-fd_de3[rownames(fd_de3)=="cg05266460",]
```
## HLA-DQB1 seems like a pretty interesting gene
### extract all CpG sites associated with HLA-DQB1 and HLA-DRB1 and make boxplot grouped by exposure


## linear model might be missing asthma status
```{r}
meta_pbmc_0_24$AsthmaStatus <- ifelse(meta_pbmc_0_24$AsthmaMeth8==1, "Asthma", "Healthy")
meta_pbmc_0_48$AsthmaStatus <- ifelse(meta_pbmc_0_48$AsthmaMeth8==1, "Asthma", "Healthy")
library(lme4)
library(lmerTest)
lm_0_24<-sapply(1:nrow(mval_pbmc_0_24),function(CpG) { 
  metaex<-meta_pbmc_0_24
  metaex$Mval<-mval_pbmc_0_24[CpG,]
  Model.REML = lmer(Mval ~ Exposure*SampleTimeHours+AsthmaStatus+(1|VolunteerNumber), REML = TRUE, data = metaex)
  anova(Model.REML)[1,6] # pval for Exposure
    }
  )
## 2018-08-29 check on the result above 
## rename the object nd turn into data frame
lm_0_24_asthstat_added<-lm_0_24
lm_0_24_asthstat_added<-as.data.frame(lm_0_24_asthstat_added)
lm_0_24_asthstat_added$CpG<-rownames(beta_pbmc_24)

summary(lm_0_24_asthstat_added$lm_0_24_asthstat_added) ## doesn't look bad
## check FDR
Multi_test_corr_relaxed<-p.adjust(lm_0_24_asthstat_added$lm_0_24_asthstat_added, method = "fdr", n = length(lm_0_24_asthstat_added$lm_0_24_asthstat_added))

summary(Multi_test_corr_relaxed)
# min is 0.044, not bad
length(Multi_test_corr_relaxed[which(Multi_test_corr_relaxed<=0.1)]) # 1
length(Multi_test_corr_relaxed[which(Multi_test_corr_relaxed<=0.2)]) # 1
length(Multi_test_corr_relaxed[which(Multi_test_corr_relaxed<=0.3)]) # 1
## LMAO 

## take a look at the pval distribution and see if the weird dip at the front gets better
ggplot(lm_0_24_asthstat_added,aes(lm_0_24_asthstat_added))+geom_histogram(fill="grey90", color="black")+xlab("Nominal P Value")+xlim(min(lm_0_24_asthstat_added$lm_0_24_asthstat_added),max(lm_0_24_asthstat_added$lm_0_24_asthstat_added))+theme_minimal()


### OMG YES YES LOOK AT THIS PRETTY PLOT!!!! THIS IS THE MODEL WE WILL USE FROM NOW ON OKAY!!!!!

```

Now that we have both the effect size and pval figured out, maybe redo volcano plot to find both significant sites and hopefully something interesting pops up!!
```{r}

```






