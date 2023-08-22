library(caret)
library(doParallel)
library(ebmc)

library(pROC)
library(ResourceSelection)
# library(randomForest)
library(e1071)
library(generalhoslem)
library(ggplot2)
library(reshape2)
library(plyr)
library(Hmisc)


#loading MODEL

load("MODEL_TAREA.RData")

#TAREA_BOOST is the GLMBOOST MODEL
#PLATT_REEC is the  PLATT recalibration

#loading patieent data
#current thershold a=0.4
a<-0.45

P<-import("PATIENT_NIBP.csv")

names(P)[names(P) == "V1"] <- "DAY"
names(P)[names(P) == "V2"] <- "Age"
names(P)[names(P) == "V3"] <- "GCS"
names(P)[names(P) == "V4"] <- "HR_NPEAKS"
names(P)[names(P) == "V5"] <- "HR_HFD"
names(P)[names(P) == "V6"] <- "HR-HF"
names(P)[names(P) == "V7"] <- "SBP_70_90"
names(P)[names(P) =="V8"] <- "DBP_0_50"
names(P)[names(P) =="V9"] <- "MBP_0_65"
names(P)[names(P) =="V10"] <- "SAPS2"
names(P)[names(P) =="V11"] <- "OUTCOME"

P_in<-P[1,c(2:9)]


pred11<- predict(TAREA_BOOST,newdata=P_in,type = "prob")[,"Class2"]
br<-data.frame(forecast=pred11)

pred12<-predict(PLATT_REC,newdata=br)
pred12<-1-1/(1+exp(pred12))





p2<-ifelse(pred12 >a,1,0)# 





igs<-P$SAPS2[1]
igs<--7.7631+0.0737*igs+0.9971*log(igs+1) 
igs<-exp(igs)/(1+exp(igs))

igs_DEC<-ifelse(igs <0.5,0,1)# 



message("the BOOST PREDICTION AFTER RECALIBRATION FOR NON-SURVIVING IS ",pred12," the BOOST PREDICTION IS THEREFORE ",p2);
message("the SAPS2 PREDICTION  FOR NON-SURVIVING IS ",igs, "the SAPS2 PREDICTION IS THEREFORE ",igs_DEC);
message("The true outcome is",P$OUTCOME)
message("O means SURVIVOR and 1 NON-SURVIVOR")