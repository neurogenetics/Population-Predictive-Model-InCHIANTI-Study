#ordinal model for medication risk Final analysis 2/27/17
getwd()
setwd("/Users/sargentlj/Google Drive/NIH/Databases/DataBases/InCHIANTI/InC-Run/Baseline_1998-2000/Medication-base")

#MMSE Cog1 Phenotype1
phenotype1 <- read.csv("/Users/sargentlj/Google Drive/NIH/Databases/DataBases/InCHIANTI/InC-Run/Baseline_1998-2000/Medication-base/Ordinal_Medication.Phenotype_Cog1.csv")
head(phenotype1)
dim(phenotype1)

#load in medications
Medications <- read.csv("/Users/sargentlj/Google Drive/NIH/Databases/DataBases/InCHIANTI/InC-Run/Baseline_1998-2000/Medication-base/Code_Medications.csv")
head(Medications)
dim(Medications)

#merge risk scores and phenotype1
med_ordinal_complete_data1 <- merge(Medications,phenotype1,by ="IID", all=TRUE)
head(med_ordinal_complete_data1)

med_ordinal_complete_data <- med_ordinal_complete_data1[,c(1:3,5,7,10:17,22:24)]
head(med_ordinal_complete_data)

summary(med_ordinal_complete_data)
write.csv(med_ordinal_complete_data,file="ordinal_age_meds.csv")

#65 and older
med_ordinal_complete_data <- read.csv("/Users/sargentlj/Google Drive/NIH/Databases/DataBases/InCHIANTI/InC-Run/Baseline_1998-2000/Medication-base/ordinal_age_meds.csv", header=TRUE)
head(med_ordinal_complete_data)
nrow(med_ordinal_complete_data)

##look at structure of data
str(med_ordinal_complete_data)
#Convert data from int to ordered data
med_ordinal_complete_data$COG1 <- as.ordered(med_ordinal_complete_data$COG1)
med_ordinal_complete_data$FRAILCODE<-as.ordered(med_ordinal_complete_data$FRAILCODE)
med_ordinal_complete_data$COGFRAIL<-as.ordered(med_ordinal_complete_data$COGFRAIL)


library(MASS)
library(ggplot2)

#Total participant
nrow(med_ordinal_complete_data)

# Demographics

#All Participants

#mean age 
mean(med_ordinal_complete_data[,"IXAGE"])
sd(med_ordinal_complete_data[,"IXAGE"])

mode(med_ordinal_complete_data[,"SEX"])
S <- med_ordinal_complete_data[,c("SEX")]
S_sum <- table(S)

#Covariates
+ SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease

##Model 1 - MMSE

#ordinal test ACB Medications 
#Hess=TRUE to have the model return the observed information matrix from optimization (called the Hessian) which is used to get standard errors
ordinal_Meds_model1 <- polr(COG1 ~ Total_ACB + SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data = med_ordinal_complete_data, Hess=TRUE)
summary(ordinal_Meds_model1)

#chi-square test  - goodness of fit Ho: current model is good enough 
1-pchisq(deviance(ordinal_Meds_model1),df.residual(ordinal_Meds_model1))

#store table 
ctable <- data.frame(coef(summary(ordinal_Meds_model1)))
ctable$pval = round((pnorm(abs(ctable$t.value), lower.tail = FALSE) * 2),2)
ctable

#cofidence intervals
ci <- confint(ordinal_Meds_model1)
confint.default(ordinal_Meds_model1)
exp(coef(ordinal_Meds_model1))
exp(cbind(OR = coef(ordinal_Meds_model1), ci))

#probabilities 
m1.pred <- predict(ordinal_Meds_model1, type="probs")
summary(m1.pred)

#number of cognitive decline particpants 
mode(med_ordinal_complete_data[,"COG1"])
Cog1 <- med_ordinal_complete_data[,c("COG1")]
Cog1_sum <- table(Cog1)


#ordinal test Benzo Medications 
ordinal_MedsB_model1 <- polr(COG1 ~ Benzodiazepine + SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data = med_ordinal_complete_data, Hess=TRUE)
summary(ordinal_MedsB_model1)

#chi-square test  - goodness of fit Ho: current model is good enough 
1-pchisq(deviance(ordinal_MedsB_model1),df.residual(ordinal_MedsB_model1))

#store table 
ctable <- data.frame(coef(summary(ordinal_MedsB_model1)))
ctable$pval = round((pnorm(abs(ctable$t.value), lower.tail = FALSE) * 2),2)
ctable

#cofidence intervals
ci <- confint(ordinal_MedsB_model1)
confint.default(ordinal_MedsB_model1)
exp(coef(ordinal_MedsB_model1))
exp(cbind(OR = coef(ordinal_MedsB_model1), ci))

#probabilities 
m1.pred <- predict(ordinal_MedsB_model1, type="probs")
summary(m1.pred)

#ACB-Frail

#Ordinal Regression 
Med_frailmodel1 <- polr(FRAILCODE ~ Total_ACB + SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data = med_ordinal_complete_data, Hess=TRUE)
summary(Med_frailmodel1)

#chi-square test  - goodness of fit Ho: current model is good enougth 
1-pchisq(deviance(Med_frailmodel1),df.residual(Med_frailmodel1))
#Full table with p-values 
ctable <- data.frame(coef(summary(Med_frailmodel1)))
ctable$pval = round((pnorm(abs(ctable$t.value), lower.tail = FALSE) * 2),2)
ctable

#cofidence intervals
ci <- confint(Med_frailmodel1)
confint.default(Med_frailmodel1)
exp(coef(Med_frailmodel1))
exp(cbind(OR = coef(Med_frailmodel1), ci))

#probabilities 
m2.pred <- predict(Med_frailmodel1, type="probs")
summary(m2.pred)

mode(med_ordinal_complete_data[,"FRAILCODE"])
frail <- med_ordinal_complete_data[,c("FRAILCODE")]
frail_sum <- table(frail)

#Frail Benzo 
Med_Benzo_frailmodel1 <- polr(FRAILCODE ~ Benzodiazepine + SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data = med_ordinal_complete_data, Hess=TRUE)
summary(Med_Benzo_frailmodel1)

#chi-square test  - goodness of fit Ho: current model is good enougth 
1-pchisq(deviance(Med_Benzo_frailmodel1),df.residual(Med_Benzo_frailmodel1))
#Full table with p-values 
ctable <- data.frame(coef(summary(Med_Benzo_frailmodel1)))
ctable$pval = round((pnorm(abs(ctable$t.value), lower.tail = FALSE) * 2),2)
ctable

#cofidence intervals
ci <- confint(Med_Benzo_frailmodel1)
confint.default(Med_Benzo_frailmodel1)
exp(coef(Med_Benzo_frailmodel1))
exp(cbind(OR = coef(Med_Benzo_frailmodel1), ci))

#probabilities 
m2.pred <- predict(Med_Benzo_frailmodel1, type="probs")
summary(m2.pred)

mode(Med_Benzo_frailmodel1[,"FRAILCODE"])
frail <- Med_Benzo_frailmodel1[,c("FRAILCODE")]
frail_sum <- table(frail)

#Cogfrail -ACB 

Med_cogfrailmodel1 <- polr(COGFRAIL ~ Total_ACB + SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data = med_ordinal_complete_data, Hess=TRUE)
summary(Med_cogfrailmodel1)

#chi-square test  - goodness of fit Ho: current model is good enougth 
1-pchisq(deviance(Med_cogfrailmodel1),df.residual(Med_cogfrailmodel1))

#Full table with p-values 
ctable <- data.frame(coef(summary(Med_cogfrailmodel1)))
ctable$pval = round((pnorm(abs(ctable$t.value), lower.tail = FALSE) * 2),2)
ctable

#cofidence intervals
ci <- confint(Med_cogfrailmodel1)
confint.default(Med_cogfrailmodel1)
exp(coef(Med_cogfrailmodel1))
exp(cbind(OR = coef(Med_cogfrailmodel1), ci))

#probabilities 
m3.pred <- predict(Med_cogfrailmodel1, type="probs")
summary(m3.pred)

mode(med_ordinal_complete_data[,"COGFRAIL"])
cogfrail <- med_ordinal_complete_data[,c("COGFRAIL")]
cogfrail_sum <- table(cogfrail)

#Cogfrail - Benzo 
Med_Benzo_cogfrailmodel1 <- polr(COGFRAIL ~ Benzodiazepine + SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data = med_ordinal_complete_data, Hess=TRUE)
summary(Med_Benzo_cogfrailmodel1)

#chi-square test  - goodness of fit Ho: current model is good enougth 
1-pchisq(deviance(Med_Benzo_cogfrailmodel1),df.residual(Med_Benzo_cogfrailmodel1))

#Full table with p-values 
ctable <- data.frame(coef(summary(Med_Benzo_cogfrailmodel1)))
ctable$pval = round((pnorm(abs(ctable$t.value), lower.tail = FALSE) * 2),2)
ctable

#cofidence intervals
ci <- confint(Med_Benzo_cogfrailmodel1)
confint.default(Med_Benzo_cogfrailmodel1)
exp(coef(Med_Benzo_cogfrailmodel1))
exp(cbind(OR = coef(Med_Benzo_cogfrailmodel1), ci))

#probabilities 
m3.pred <- predict(Med_Benzo_cogfrailmodel1, type="probs")
summary(m3.pred)

#######Model 2 ordinal analysis######
##Important note - this data.frame would need to be recoded if ordinal model is to be used. All > 300 scores for Trail A & B are noted in VXTRAIAS (n=6) and VXTRAIBS (n=69)
Model2.ordinal <- read.csv("/Users/sargentlj/Documents/InChiantiSNP/AllMergedandRisk/Risk-Score-results/Recode_Cog2_cli_rawe.csv")
head(Model2.ordinal)

#merge risk scores and phenotype2
med_ordinal_complete_data2 <- merge(Medications,Model2.ordinal,by ="IID", all=TRUE)
head(med_ordinal_complete_data2)

med_ordinal_complete_data2$NTrailB.Code <- as.ordered(med_ordinal_complete_data2$NTrailB.Code)

#Cog2 ACB Trail B ordinal analysis 
ordinal_Meds_model2 <- polr(NTrailB.Code ~ Total_ACB, data = med_ordinal_complete_data2, Hess=TRUE)
summary(ordinal_Meds_model2)

#chi-square test  - goodness of fit Ho: current model is good enough 
1-pchisq(deviance(ordinal_Meds_model2),df.residual(ordinal_Meds_model2))

#store table 
ctable <- data.frame(coef(summary(ordinal_Meds_model2)))
ctable$pval = round((pnorm(abs(ctable$t.value), lower.tail = FALSE) * 2),2)
ctable

#cofidence intervals
ci <- confint(ordinal_Meds_model2)
confint.default(ordinal_Meds_model2)
exp(coef(ordinal_Meds_model2))
exp(cbind(OR = coef(ordinal_Meds_model2), ci))

#probabilities 
m2.pred <- predict(ordinal_Meds_model2, type="probs")
summary(m2.pred)

#Cog2 ACB Trail A ordinal analysis 

med_ordinal_complete_data2$NTrailA.Code <- as.ordered(med_ordinal_complete_data2$NTrailA.Code)

#Cog2 ACB Trail A ordinal analysis 
TrailA.ordinal_Meds_model2 <- polr(NTrailA.Code ~ Total_ACB, data = med_ordinal_complete_data2, Hess=TRUE)
summary(TrailA.ordinal_Meds_model2)

#chi-square test  - goodness of fit Ho: current model is good enough 
1-pchisq(deviance(TrailA.ordinal_Meds_model2),df.residual(TrailA.ordinal_Meds_model2))

#store table 
ctable <- data.frame(coef(summary(TrailA.ordinal_Meds_model2)))
ctable$pval = round((pnorm(abs(ctable$t.value), lower.tail = FALSE) * 2),2)
ctable

#cofidence intervals
ci <- confint(TrailA.ordinal_Meds_model2)
confint.default(TrailA.ordinal_Meds_model2)
exp(coef(TrailA.ordinal_Meds_model2))
exp(cbind(OR = coef(TrailA.ordinal_Meds_model2), ci))

#probabilities 
m2.pred <- predict(TrailA.ordinal_Meds_model2, type="probs")
summary(m2.pred)

##CogFrail ACB Trail B ordinal analysis 

med_ordinal_complete_data2$NCOGFRAIL_TrailB <- as.ordered(med_ordinal_complete_data2$NCOGFRAIL_TrailB)

CogFrail_B.ordinal_Meds_model2 <- polr(NCOGFRAIL_TrailB ~ Total_ACB, data = med_ordinal_complete_data2, Hess=TRUE)
summary(CogFrail_B.ordinal_Meds_model2)

#chi-square test  - goodness of fit Ho: current model is good enough 
1-pchisq(deviance(CogFrail_B.ordinal_Meds_model2),df.residual(CogFrail_B.ordinal_Meds_model2))

#store table 
ctable <- data.frame(coef(summary(CogFrail_B.ordinal_Meds_model2)))
ctable$pval = round((pnorm(abs(ctable$t.value), lower.tail = FALSE) * 2),2)
ctable

#cofidence intervals
ci <- confint(CogFrail_B.ordinal_Meds_model2)
confint.default(CogFrail_B.ordinal_Meds_model2)
exp(coef(CogFrail_B.ordinal_Meds_model2))
exp(cbind(OR = coef(CogFrail_B.ordinal_Meds_model2), ci))

#probabilities 
m2.pred <- predict(CogFrail_B.ordinal_Meds_model2, type="probs")
summary(m2.pred)