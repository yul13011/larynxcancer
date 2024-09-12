
library(haven)
LarynxSCC_DSS <- read_sav("C:/Users/yul13011/Dropbox/2017-2018/Linda's paper/Revisions/SEER2000on_LarynxSCC_2019Feb22.sav")

LarynxSCC_DSS=LarynxSCC_DSS[LarynxSCC_DSS$`filter_$`==1, ]

LarynxSCC_DSS=LarynxSCC_DSS[, c("DSS_Larynx","SurvivalMonths","Age","Sex","Race_NEW","PrimarySiteLabeled","PrimarySite","LNDx", "AJCC6thStage_NEW","TreatmentGrp","Marital_NEW")]

LarynxSCC_DSS=LarynxSCC_DSS[LarynxSCC_DSS$LNDx!="99",]
LarynxSCC_DSS=LarynxSCC_DSS[LarynxSCC_DSS$Marital_NEW!="6",]
LarynxSCC_DSS=LarynxSCC_DSS[LarynxSCC_DSS$AJCC6thStage_NEW!="8",]


LarynxSCC_DSS$PrimarySite=factor(LarynxSCC_DSS$PrimarySite, 
                                 levels=c("320","321","322"),
                                 labels=c("C32.0-Glottis","C32.1-Supraglottis","C32.2-Subglottis"))
LarynxSCC_DSS$Race_NEW=factor(LarynxSCC_DSS$Race_NEW,
                              levels=c("5","1","2","3","4"),
                              labels=c("White",
                                       "American Indian/Alaska Native",
                                       "Asian or Pacific Islander",
                                       "Black",
                                       "Unknown"))
                              
LarynxSCC_DSS$Marital_NEW=factor(LarynxSCC_DSS$Marital_NEW,
                                 levels=c("4","1","2","3","5","7"),
                                 labels=c("Single","Divorced","Married","Separated","Unknown","widowed"))

LarynxSCC_DSS$AJCC6thStage_NEW=factor(LarynxSCC_DSS$AJCC6thStage_NEW,
                                      levels = c("1","2","3","4"))
LarynxSCC_DSS$LNDx=factor(LarynxSCC_DSS$LNDx,
                          levels = c("0","1"),
                          labels = c("not examined","examined"))

LarynxSCC_DSS$TreatmentGrp=factor(LarynxSCC_DSS$TreatmentGrp,
                                  levels = c("0","1","2","3"),
                                  labels = c("Radiation only","Local Surgery",
                                             "Local surgery + radiation",
                                             "Open surgery"))
LarynxSCC_DSS$female=ifelse(LarynxSCC_DSS$Sex=="Female",1,0) # create a variable to match on #
LarynxSCC_DSS$DSS_Larynx=ifelse(LarynxSCC_DSS$Sex=="Female",1,0) # create a variable to match on #

LarynxSCC_DSS=data.frame(LarynxSCC_DSS)



## Run Cox regression on the unmatched data and check proportional hazards assumption 

library(survival)
library(survminer)
library(dplyr)
library(timereg)

surv=Surv(time=LarynxSCC_DSS$SurvivalMonths,event= LarynxSCC_DSS$DSS_Larynx)


# K-M curve for males vs. females
fit1=survfit(surv~female,data=LarynxSCC_DSS)
ggsurvplot(fit1,data=LarynxSCC_DSS,pval=TRUE,
           legend="bottom",legend.title="Sex",legend.labs=c("Female","Male"))


# K-M curve for males vs. females with glottic tumors
LarynxSCC_DSS1=LarynxSCC_DSS[LarynxSCC_DSS$PrimarySite=="C32.0-Glottis",]
surv_object1=Surv(time=LarynxSCC_DSS1$SurvivalMonths,event=LarynxSCC_DSS1$DSS_Larynx)
fit11=survfit(surv_object1~female,data=LarynxSCC_DSS1)
ggsurvplot(fit11, pval=TRUE,title = "Glottis", legend.title="Sex",legend.labs=c("Female","Male"),
           risk.table = TRUE, tables.theme = theme_cleantable(),xlab="Time (Month)")


# K-M curve for males vs. females with supraglottis tumors
LarynxSCC_brief2=LarynxSCC_DSS[LarynxSCC_DSS$PrimarySite=="C32.1-Supraglottis",]
surv_object2=Surv(time=LarynxSCC_brief2$SurvivalMonths,event=LarynxSCC_brief2$DSS_Larynx)
fit12=survfit(surv_object2~female,data=LarynxSCC_brief2)
ggsurvplot(fit12, pval=TRUE,title = "Supraglottis", legend.title="Sex",legend.labs=c("Female","Male"),
           risk.table = TRUE, tables.theme = theme_cleantable())

# K-M curve for males vs. females with subglottic tumors
LarynxSCC_brief3=LarynxSCC_DSS[LarynxSCC_DSS$PrimarySite=="C32.2-Subglottis",]
surv_object3=Surv(time=LarynxSCC_brief3$SurvivalMonths,event=LarynxSCC_brief3$DSS_Larynx)
fit13=survfit(surv_object3~female,data=LarynxSCC_brief3)
ggsurvplot(fit13, pval=TRUE,title = "Subglottis", legend.title="Sex",legend.labs=c("Female","Male"),
           risk.table = TRUE, tables.theme = theme_cleantable())


# Fit cox regression 
fit1=coxph(surv~Age+female+Race_NEW+Marital_NEW+
             PrimarySite+AJCC6thStage_NEW+TreatmentGrp+LNDx, data=LarynxSCC_DSS)
summary(fit1)

# Evaluate model assumption 
temp1=cox.zph(fit1)
print(temp1)



### We find out that the data does not satisfy a proportional hazards assumption
### Try a different model 
### Fit nonparametric additive hazards model which allows every independent variable to have time-varying effects


LarynxSCC_DSS$age.c = LarynxSCC_DSS$Age-mean(LarynxSCC_DSS$Age)


add1=aalen(surv~I(exp(age.c/10))+female+Race_NEW+Marital_NEW
           +PrimarySite+AJCC6thStage_NEW+TreatmentGrp+LNDx,
           data=LarynxSCC_DSS,residuals = 1)
summary(add1)
plot(add1,specific.comps = 3)



### AJCC6thStage_NEW failed the tests for constant effect
### Fit semi-parametric additive hazards model 


LarynxSCC_DSS$age.c=LarynxSCC_DSS$Age-mean(LarynxSCC_DSS$Age)

add2=aalen(surv~const(I(exp(age.c/10)))+const(female)+const(Race_NEW)+const(Marital_NEW)
           +const(PrimarySite)+AJCC6thStage_NEW+const(TreatmentGrp)+const(LNDx),
           data=LarynxSCC_DSS)
summary(add2)

### PS Matching and check balance ###


library(MatchIt)
library("optmatch")
library(cobalt)

# Now try full matching # 
m.out=matchit(female~Age+TreatmentGrp+Race_NEW+Marital_NEW+LNDx+AJCC6thStage_NEW+PrimarySite, data=LarynxSCC_DSS, 
              method="full")
summary(m.out,standardize = T)

m.data_full=match.data(m.out,group="all",distance = "distance",weights="weights",subclass = "subclass")

# Check balance by graphical evidence
bal.plot(m.out,var.name = "distance",which="both",type="histogram")+
  scale_fill_discrete(name="Sex",labels=c("Male","Female"))+
  ggtitle("")+xlab("Propensity Score")

# Balance is good # 

### Run additive hazards models on the post-matching data 
surv2=Surv(time=m.data_full$SurvivalMonths,event= m.data_full$DSS_Larynx)

fit2=coxph(surv2~Age+female+Race_NEW+Marital_NEW+
             PrimarySite+AJCC6thStage_NEW+TreatmentGrp+LNDx, data=m.data_full,weights=weights)
summary(fit2)
ggforest(fit1,data=LarynxSCC_DSS)

temp1=cox.zph(fit2)
print(temp1)

m.data_full$age.c=m.data_full$Age-mean(m.data_full$Age)

add3=aalen(surv2~I(exp(age.c/10))+female+Race_NEW+Marital_NEW+
             PrimarySite+AJCC6thStage_NEW+TreatmentGrp+LNDx,
           data=m.data_full, m.data_full$weights,residuals = 1)
summary(add3)


### Stage and LDX failed test for constant effect 
### Fit semi-parametric additive hazards model 


add4=aalen(surv2~const(I(exp(age.c/10)))+const(female)+const(Race_NEW)+const(Marital_NEW)
           +const(PrimarySite)+AJCC6thStage_NEW+const(TreatmentGrp)+const(LNDx),
           data=m.data_full,weights=m.data_full$weights)
summary(add4)

