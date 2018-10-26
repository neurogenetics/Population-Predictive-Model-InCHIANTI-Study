#ordinal model for medication risk Final analysis and updates to tables 3-4 in manuscript on 10-3-18 adding pseudo r-squared McFadden Reported Odds conversion Log_eff
getwd()
setwd("load_file/create container")

#MMSE Cog1 Phenotype1
#load phenotype1 
head(phenotype1)
dim(phenotype1)

#load in medications
Medications <- read.csv("/Volumes/GoogleDrive/My Drive/NIH/Databases/DataBases/InCHIANTI/InC-Run/Baseline_1998-2000/Medication-base/Code_Medications.csv")
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
med_ordinal_complete_data <- read.csv("load ordinal dataframe", header=TRUE)
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
# evaluate model fit 
ordinal_Meds_model1 <- polr(COG1 ~ Total_ACB + SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data = med_ordinal_complete_data, Hess=TRUE)
summary(ordinal_Meds_model1)
library(pscl)
pR2(ordinal_Meds_model1)

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


#ACB-Frail

#Ordinal Regression 
#Hess=TRUE to have the model return the observed information matrix from optimization (called the Hessian) which is used to get standard errors
# evaluate model fit 
Med_frailmodel1 <- polr(FRAILCODE ~ Total_ACB + SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data = med_ordinal_complete_data, Hess=TRUE)
summary(Med_frailmodel1)
pR2(Med_frailmodel1)
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


#Cogfrail -ACB 
#Hess=TRUE to have the model return the observed information matrix from optimization (called the Hessian) which is used to get standard errors
# evaluate model fit 
Med_cogfrailmodel1 <- polr(COGFRAIL ~ Total_ACB + SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data = med_ordinal_complete_data, Hess=TRUE)
summary(Med_cogfrailmodel1)
pR2(Med_cogfrailmodel1)

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

