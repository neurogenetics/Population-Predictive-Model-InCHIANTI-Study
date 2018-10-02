##glm model for medication risk Final analysis 5/13/2017 second analysis done with publication 
setwd("/Users/sargentlj/Google Drive/NIH/Databases/DataBases/InCHIANTI/InC-Run/Baseline_1998-2000/Medication-base")
getwd()

#Fulldata.frame used for boosting analyses
fulldata.frame <- read.csv("/Users/sargentlj/Google Drive/NIH/Databases/DataBases/InCHIANTI/InC-Run/Baseline_1998-2000/Boosting_Analysis/Recode.Final.Data.Frame.csv",header=TRUE)
head(fulldata.frame)
nrow(fulldata.frame)

#mean age 
mean(fulldata.frame[,"IXAGE"])
sd(fulldata.frame[,"IXAGE"])

mode(fulldata.frame[,"SEX"])
S <- fulldata.frame[,c("SEX")]
S_sum <- table(S)

#Total ACB burden
mode(fulldata.frame[,"Total_ACB"])
S_ACB <- fulldata.frame[,c("Total_ACB")]
S_sum_ACB <- table(S_ACB)
sum(is.na(fulldata.frame$Total_ACB))

#table ACB burden by phenotype 
#Cognition
Cog1 <- fulldata.frame[,c("COG1","Total_ACB")]
S_sum_Cog1<- table(Cog1)

Frail <- fulldata.frame[,c("FRAILCODE","Total_ACB")]
S_sum_Frail<- table(Frail)

CogFrail <- fulldata.frame[,c("COGFRAIL","Total_ACB")]
S_sum_CogFrail<- table(CogFrail)

cogB <- fulldata.frame[,c("TrailB.Code","Total_ACB")]
S_sum_cogB <- table(cogB)

cogA <- fulldata.frame[,c("TrailA.Code","Total_ACB")]
S_sum_cogA <- table(cogA)

CogFrailB <- fulldata.frame[,c("binomial_COGFRAIL_B","Total_ACB")]
S_sum_CogFrailB <- table(CogFrailB)

CogFrailA <- fulldata.frame[,c("COGFRAIL_TrailA","Total_ACB")]
S_sum_CogFrailA <- table(CogFrailA)



library(MASS)
library(aod)

str(fulldata.frame)
#Convert data from int to factor
fulldata.frame$COG1 <- factor(fulldata.frame$COG1)
fulldata.frame$COGFRAIL <- factor(fulldata.frame$COGFRAIL)
fulldata.frame$FRAILCODE <- factor(fulldata.frame$FRAILCODE)
fulldata.frame$TrailB.Code <-factor(fulldata.frame$TrailB.Code)
fulldata.frame$TrailA.Code <- factor(fulldata.frame$TrailA.Code)
fulldata.frame$binomial_COGFRAIL_B  <-factor(fulldata.frame$binomial_COGFRAIL_B )
fulldata.frame$COGFRAIL_TrailA <- factor(fulldata.frame$COGFRAIL_TrailA)


#Covariates
+ SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease

#Analysis#

#Medication Analysis Cognition1- ACB
binomial_cog1_Med <- glm(COG1 ~ Total_ACB + SEX + IXAGE + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data = fulldata.frame, family = binomial)
summary(binomial_cog1_Med)
confint(binomial_cog1_Med)
wald.test(b = coef(binomial_cog1_Med), Sigma = vcov(binomial_cog1_Med), Terms = 4:6)

#Check number of Cog 
mode(fulldata.frame[,"COG1"])
Cog1 <- fulldata.frame[,c("COG1")]
Cog1_sum <- table(Cog1)

#Medication Analysis Frail1
binomial_Med_Frail1 <- glm(FRAILCODE ~ Total_ACB + SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data = fulldata.frame, family = binomial)
summary(binomial_Med_Frail1)
confint(binomial_Med_Frail1)
wald.test(b = coef(binomial_Med_Frail1), Sigma = vcov(binomial_Med_Frail1), Terms = 4:6)

#Check number of FRAILCODE 
mode(fulldata.frame[,"FRAILCODE"])
FRAILCODE1 <- fulldata.frame[,c("FRAILCODE")]
FRAILCODE1_sum <- table(FRAILCODE1)

#Medication Analysis Cognitive-Frailty1
binomial_Med_cogfrail1 <- glm(COGFRAIL ~ Total_ACB + SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data = fulldata.frame, family=binomial)
summary(binomial_Med_cogfrail1)
confint(binomial_Med_cogfrail1)
wald.test(b = coef(binomial_Med_cogfrail1), Sigma = vcov(binomial_Med_cogfrail1), Terms = 4:6)

#Check number of COGFRAIL 
mode(fulldata.frame[,"COGFRAIL"])
COGFRAIL1 <- fulldata.frame[,c("COGFRAIL")]
COGFRAIL1_sum <- table(COGFRAIL1)

##Model 2

#Trail B Meds
#no covariates if you add the covariates the model become insignificant 
##too many variables  in the model causing over-fitting can not run CI or wald test for TMT A & B only 
trailB_Meds <- glm(TrailB.Code ~ Total_ACB, data=fulldata.frame, family=binomial)
summary(trailB_Meds)
confint(trailB_Meds)
wald.test(b = coef(trailB_Meds_C), Sigma = vcov(trailB_Meds_C), Terms = 4:6)

trailB_Meds_C <- glm(TrailB.Code ~ Total_ACB + SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data=fulldata.frame, family=binomial)
summary(trailB_Meds_C)

#Check number of TrailB.Code 
mode(fulldata.frame[,"TrailB.Code"])
TrailB.Code <- fulldata.frame[,c("TrailB.Code")]
TrailB.Code_sum <- table(TrailB.Code)

#no covariates if you add the covariates the model become insignificant 
#Trail A Meds
trailA_Meds <- glm(TrailA.Code ~ Total_ACB, data=fulldata.frame, family=binomial)
summary(trailA_Meds)
confint(trailA_Meds)
wald.test(b = coef(trailA_Meds), Sigma = vcov(trailA_Meds), Terms = 4:6)

trailA_Meds_C <- glm(TrailA.Code ~ Total_ACB + SEX + IXAGE + AXDEMENT + AXVASDEM_VAS + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data=fulldata.frame, family=binomial)
summary(trailA_Meds_C)

#Check number of TrailA.Code 
mode(fulldata.frame[,"TrailA.Code"])
TrailA.Code <- fulldata.frame[,c("TrailA.Code")]
TrailA.Code_sum <- table(TrailA.Code)

#Cogfrail B Meds
CogFrail_trail_Meds <- glm(binomial_COGFRAIL_B ~ Total_ACB + SEX + IXAGE + AXDEMENT + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data=fulldata.frame, family=binomial)
summary(CogFrail_trail_Meds)
confint(CogFrail_trail_Meds)
wald.test(b = coef(CogFrail_trail_Meds), Sigma = vcov(CogFrail_trail_Meds), Terms = 4:6)

#Check number of binomial_COGFRAIL_B 
mode(fulldata.frame[,"binomial_COGFRAIL_B"])
binomial_COGFRAIL_B <- fulldata.frame[,c("binomial_COGFRAIL_B")]
binomial_COGFRAIL_B_sum <- table(binomial_COGFRAIL_B)


#Cogfrail A Meds

CogFrail_trailA_Meds <- glm(COGFRAIL_TrailA ~ Total_ACB + SEX + IXAGE + AXDEMENT + IX1_V26_level_education + IXCESD_T_Depression + VX7_V12_Parkinson.s_Disease, data=fulldata.frame, family=binomial)
summary(CogFrail_trailA_Meds)
confint(CogFrail_trailA_Meds)
wald.test(b = coef(CogFrail_trailA_Meds), Sigma = vcov(CogFrail_trailA_Meds), Terms = 4:6)

#Check number of COGFRAIL_TrailA 
mode(fulldata.frame[,"COGFRAIL_TrailA"])
COGFRAIL_TrailA <- fulldata.frame[,c("COGFRAIL_TrailA")]
COGFRAIL_TrailA_sum <- table(COGFRAIL_TrailA)






