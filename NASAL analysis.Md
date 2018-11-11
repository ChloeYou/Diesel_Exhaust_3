---
title: "DE3-NASAL"
author: "cyou"
date: "8/7/2018"
output: html_document
---

## try SVA spike in of blood data
### subset meta data 
```{r}
## subset meta data for both tissues and save file 
NASAL_meta<-filtermeta[filtermeta$SampleType=="NASAL",]
de3_bmiq_NASAL<-de3_bmiq[,sampleNames(de3_bmiq)%in%NASAL_meta$SampleID]
dim(de3_bmiq_NASAL)
```

### seperate meta of PBMC and NASAL 
```{r}
PBMC_meta<-filtermeta[filtermeta$SampleType=="PBMC",]
## combine PBMC AND NASAL meta data
NP_meta<-rbind(PBMC_meta,NASAL_meta) # 145 
## extract PBMC and NASAL de3_bmiq
de3_bmiq_NP<-de3_bmiq[,sampleNames(de3_bmiq)%in%NP_meta$SampleID]

## reorder NP_meta$SampleID
 for (k in 1:length(sampleNames(de3_bmiq_NP))) {
  for (i in 1:nrow(NP_meta)) {
    if(NP_meta$SampleID[i]== sampleNames(de3_bmiq_NP)[k]){
      temp<-NP_meta[k,]
      NP_meta[k,]<-NP_meta[i,]
      NP_meta[i,]<-temp
    }
  }
}

#### 2018-08-07
## sup or unsup(?) SVA
## SVA, courtesy of Rachel's code on GitHub
beta<-as.data.frame(betas(de3_bmiq_NP))

## need Mval for SVA
mval_variable = apply(as.data.frame(beta), 1, Mval) # need mvalues for combat
mval_variable = as.data.frame(mval_variable)
mval_variable = t(mval_variable)

#Surrogate variable analysis: https://www.bioconductor.org/packages/release/bioc/html/sva.html
library(sva)

#Null model matrix must be nested in the full model matrix
mod = model.matrix(~NP_meta$VolunteerNumber)
mod0 = model.matrix(~1, data.frame(NP_meta$VolunteerNumber))

save(mval_variable,mod,mod0,file = "NASAL_sva.RData")
## surrogates
svobj = sva(as.matrix(mval_variable),mod,mod0)
sv_unsup_ind_NP<-svobj$sv

#save(sv_unsup_ind_NP,file="NASAL and PBMC SVA.RData")
dim(sv_unsup_ind_NP) # 145(PBMC+NASAL) 7

## Checking for SV's that are related to tissue types
sv_unsup_ind_NP<-as.data.frame(sv_unsup_ind_NP)
rownames(sv_unsup_ind_NP)<-colnames(betas(de3_bmiq_NP))

############ CHECK MORE SV OTHER THAN SV1 

## plot SV1 against SV2 and colour by Tissue Type
ggplot(sv_unsup_ind_NP, aes(V1, V2, color=as.factor(NP_meta$SampleType)))+geom_point(shape=19)+theme_bw()+xlab("SV1")+ylab("SV2")

ggplot(sv_unsup_ind_NP, aes(V1, V7, color=as.factor(NP_meta$SampleType)))+geom_point(shape=19)+theme_bw()+xlab("SV1")+ylab("SV3")
## conclusion: tissue difference is only reflected in SV1

sv_unsup_ind_NP$SampleID<-rownames(sv_unsup_ind_NP)
SV1_NP<-sv_unsup_ind_NP[,c("V1","SampleID")]
SV1_NP<-SV1_NP[SV1_NP$SampleID%in%NASAL_meta$SampleID,]

## this step chunk can change depending on how many/which SV will be used for the regression model 
## correct usng SV1
diff <- as.data.frame(SV1_NP)
identical(colnames(betas(de3_bmiq_NASAL)), rownames(diff)) # TRUE
betas <- betas(de3_bmiq_NASAL)

save(betas,diff,file="nasal_terminal.RData")
## load"nasal_terminal.RData" and run in terminal
avebeta.lm<-apply(betas, 1, function(x){
  nasal <- diff[colnames(betas),] ## extract NASAL SV from the combined data
  lm(x~V1, data=nasal) ## run lm only on NASAL data 
})

## extract matrix of residuals from resulting linear models
residuals <- t(sapply(avebeta.lm, function(x)residuals(summary(x)))) ## the residual here is biological differences? 
colnames(residuals) <- colnames(betas) # re-name residuals columns with sample names
save(residuals,file="residual_nasal.RData")

## add the residuals of each regression model to the mean methylation value of each probe (mean across all samples) to obtain the “adjusted” methylation data.
adj.residuals_nasal <- residuals+matrix(apply(betas, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals)) ## this step is key
head(adj.residuals_nasal)
save(adj.residuals_nasal,file ="adjusted_nasal_data.RData") # will save later too but just to be safe yknow

### IS THIS STEP NECESSARY FOR THIS DATA????????????
# To make sure we do not induce any NAs into the dataset when we convert the beta values back M-values (by log2 transformation), we need to ensure we do not have any corrected beta values that are greater or equal to zero or any beta values that are greater than 1.

adj.residuals_nasal[adj.residuals_nasal<=0] <- 0.001 # convert any values that are less than or equal to zero to 0.001
adj.residuals_nasal[adj.residuals_nasal>1] <- 0.999 # convert any values that are greater than 1 to 0.999

## compare the before and after adjustment beta value distribution
### before
## BetaValue distribtuions combined into a format ggplot can understand 
Beta_sample_melted<- melt(betas(de3_bmiq_NASAL))
#add meta
Beta_Plot<-merge(Beta_sample_melted,NASAL_meta,by.x="X2" ,by.y="SampleID")

ggplot(Beta_Plot, aes(value,colour=as.factor(substr(Beta_Plot_adj$X2,1,3)))) + geom_density(size=1)+theme(legend.position='none')

ggplot(Beta_Plot, aes(value, group=X2,color=VolunteerNumber))+ geom_density(size=1)+theme_bw()+facet_wrap(~VolunteerNumber)

### after
## BetaValue distribtuions combined into a format ggplot can understand 
Beta_sample_melted<- melt(adj.residuals_nasal)
dim(adj.residuals_nasal) # 37
#add meta
Beta_Plot_adj<-merge(Beta_sample_melted,NASAL_meta,by.x="X2" ,by.y="SampleID")

ggplot(Beta_Plot_adj, aes(value,colour=as.factor(substr(Beta_Plot_adj$X2,1,3)))) + geom_density(size=1)+facet_wrap(~as.factor(substr(Beta_Plot_adj$X2,1,3)),nrow=3)+theme(legend.position='none')

ggplot(Beta_Plot_adj, aes(value, group=X2,color=VolunteerNumber))+ geom_density(size=1)+theme_bw()+facet_wrap(~VolunteerNumber)

save(adj.residuals_nasal,file ="adjusted_nasal_data.RData")
```

### replot the four tissue against PC1 stuff and see where NASAL falls into 
 
```{r}
## redoing this PCA to pull out the PC1 for NASAL
## PC1 against all four tissue types 
subset<-na.omit(betas(de3_bmiq))
PCA_full<- prcomp(subset)
Loadings<-as.data.frame(PCA_full$rotation) 
Loadings$TissueType<- filtermeta$SampleType

Loadings_nasal<-Loadings[Loadings$TissueType=="NASAL",] # 37 226

ggplot(Loadings, aes(PC1, PC2, color=as.factor(TissueType)))+geom_point(shape=19)+theme_bw()+ylab("PC2")+xlab("PC1 for all four tissue types") ## blood is already corrected 

ggplot(Loadings, aes(PC1, PC3, color=as.factor(TissueType)))+geom_point(shape=19)+theme_bw()+ylab("PC3")+xlab("PC1 for all four tissue types") ## blood is already corrected 
## PC3 also shows clusters of BAL,PBMC,BRUSH
## NASAL is still scattered

# boxplot of NASAL PC1 grouped by asthma status
NASAL_meta<-filtermeta[filtermeta$SampleType=="NASAL",]
identical(as.character(NASAL_meta$SampleID),as.character(rownames(Loadings_nasal))) #True
NASAL_meta$AsthmaStatus <- ifelse(NASAL_meta$AsthmaMeth8==1, "Asthma", "Healthy")
Loadings_nasal$AsthmaStatus<-NASAL_meta$AsthmaStatus

ggplot(Loadings_nasal, aes(as.factor(AsthmaStatus),PC1, color=as.factor(NASAL_meta$AsthmaStatus)))+theme_minimal()+ylab("NASAL PC1")+xlab("Asthma Status")+geom_boxplot()+geom_point(shape=19)+theme(legend.position='none')

# boxplot of NASAL PC2 grouped by asthma status
ggplot(Loadings_nasal, aes(as.factor(AsthmaStatus),PC2, color=as.factor(NASAL_meta$AsthmaStatus)))+theme_minimal()+ylab("NASAL PC2")+xlab("Asthma Status")+geom_boxplot()+geom_point(shape=19)+theme(legend.position='none')

## perform t-test on the two gorups?
asthma<-Loadings_nasal$PC1[which(Loadings_nasal$AsthmaStatus=="Asthma")]
healthy<-Loadings_nasal$PC1[which(Loadings_nasal$AsthmaStatus=="Healthy")]
t.test(asthma,healthy) # p-value = 0.0004768
wilcox.test(asthma,healthy) #p-value = 0.0001225

asthma<-Loadings_nasal$PC2[which(Loadings_nasal$AsthmaStatus=="Asthma")]
healthy<-Loadings_nasal$PC2[which(Loadings_nasal$AsthmaStatus=="Healthy")]
t.test(asthma,healthy) # p-value =  0.005436
wilcox.test(asthma,healthy) #p-value = 0.005938

## PC1 to PC5, asthma and healthy are *statistically* significantly different
```
![NASAL SV1 against SV2 plotted seperated by tissue](https://github.com/kobor-lab/Co-op_Projects/blob/master/ChloeYou/DE3/NASAL%20SV1-SV2.png)

![NASAL beta value distribution before and after correction](https://github.com/kobor-lab/Co-op_Projects/blob/master/ChloeYou/DE3/NASAL_correction_individual.png)

![NASAL_PC1_Asthma_Status_boxplot](https://github.com/kobor-lab/Co-op_Projects/blob/master/ChloeYou/DE3/nasal_pc1_asthma%20status.png)
![NASAL_PC2_Asthma_Status_boxplot](https://github.com/kobor-lab/Co-op_Projects/blob/master/ChloeYou/DE3/nasal_pc2_asthma%20status.png)
- During lab meeting, Mike Suggested that the nasal profile might be confounded with asthma status of the individual, so we pulled out the PC1 and PC2 score and made boxplot and performed statistical tests to see if the average PC score of healthy and asthmatic individuals are significantly different.  (Yes they were!)
