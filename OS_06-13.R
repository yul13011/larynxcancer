
################################### Load data ###########################################

library(readxl)
SEER2000on_LarynxSCC_012117_2004onwards = read_excel("C:/Users/yul13011/Dropbox/2017-2018/Linda's paper/SEER2000on_LarynxSCC_012117_2004onwards.xlsx")


# Subset data to exclude some cases 
LarynxSCC=SEER2000on_LarynxSCC_012117_2004onwards[SEER2000on_LarynxSCC_012117_2004onwards$`filter_$`==1, ]

# Subset data to included variables on age, gender, race, diseases site, LNB, NM, T-stage, Treatment, Marital status
LarynxSCC_brief=LarynxSCC[, c("VitalStatus_NEW","SurvivalMonths","Age","Sex_NEW","Race_NEW","PrimarySiteLabeled","PrimarySite","LNDx", "AJCC6thStage_NEW","TreatmentGrp","Marital_NEW")]
dim(LarynxSCC_brief)
names(LarynxSCC_brief)

LarynxSCC_brief=LarynxSCC_brief[LarynxSCC_brief$LNDx!="99",]
LarynxSCC_brief=LarynxSCC_brief[LarynxSCC_brief$Marital_NEW!="6",]
LarynxSCC_brief=LarynxSCC_brief[LarynxSCC_brief$AJCC6thStage_NEW!="8",]

# Convert variables to categorical variables 
LarynxSCC_brief$Sex_NEW=factor(LarynxSCC_brief$Sex_NEW,
                               levels=c("2","1"),
                               labels=c("male","female"))
LarynxSCC_brief$PrimarySite=factor(LarynxSCC_brief$PrimarySite, 
                                   levels=c("320","321","322"),
                                   labels=c("C32.0-Glottis","C32.1-Supraglottis","C32.2-Subglottis"))
LarynxSCC_brief$Race_NEW=factor(LarynxSCC_brief$Race_NEW,
                                levels=c("5","1","2","3","4"),
                                labels=c("White",
                                         "American Indian/Alaska Native",
                                         "Asian or Pacific Islander",
                                         "Black",
                                         "Unknown"
                                ))
LarynxSCC_brief$Marital_NEW=factor(LarynxSCC_brief$Marital_NEW,
                                   levels=c("4","1","2","3","5","7"),
                                   labels=c("Single","Divorced","Married","Separated","Unknown","widowed"))

LarynxSCC_brief$AJCC6thStage_NEW=factor(LarynxSCC_brief$AJCC6thStage_NEW,
                                        levels = c("1","2","3","4"))
LarynxSCC_brief$LNDx=factor(LarynxSCC_brief$LNDx,
                            levels = c("0","1"),
                            labels = c("not examined","examined"))

LarynxSCC_brief$TreatmentGrp=factor(LarynxSCC_brief$TreatmentGrp,
                                    levels = c("0","1","2","3"),
                                    labels = c("Radiation only","Local Surgery",
                                               "Local surgery + radiation",
                                               "Open surgery"))
LarynxSCC_brief$female=ifelse(LarynxSCC_brief$Sex_NEW=="female",1,0) # create a variable to match on #
LarynxSCC_brief=data.frame(LarynxSCC_brief)

########################### KM Curve ###############################


library(survival)
library(survminer)
library(dplyr)
library(timereg)


surv_object=Surv(time=LarynxSCC_brief$SurvivalMonths,event=LarynxSCC_brief$VitalStatus_NEW)

# K-M curve for males vs. females
fit1=survfit(surv_object~Sex_NEW,data=LarynxSCC_brief)
ggsurvplot(fit1,data=LarynxSCC_brief,pval=TRUE,
           legend="bottom",legend.title="Sex",legend.labs=c("Male","Female"))

fit2=survfit(surv_object~Marital_NEW,data=LarynxSCC_brief)
ggsurvplot(fit2,data=LarynxSCC_brief,pval=TRUE,xlim=c(0,60))

# K-M curve for males vs. females with glottic tumors
LarynxSCC_brief1=LarynxSCC_brief[LarynxSCC_brief$PrimarySite=="C32.0-Glottis",]
surv_object1=Surv(time=LarynxSCC_brief1$SurvivalMonths,event=LarynxSCC_brief1$VitalStatus_NEW)
fit11=survfit(surv_object1~Sex_NEW,data=LarynxSCC_brief1)
ggsurvplot(fit11, pval=TRUE,title = "Glottis", legend.title="Sex",legend.labs=c("Male","Female"),
           risk.table = TRUE, tables.theme = theme_cleantable(),xlab="Time (Month)")

ggsurvplot(fit11, pval=TRUE,title = "Glottis", legend.title="Sex",legend.labs=c("Male","Female"),
           xlab="Time (Month)",xlim=c(0,60))

# K-M curve for males vs. females with supraglottis tumors
LarynxSCC_brief2=LarynxSCC_brief[LarynxSCC_brief$PrimarySite=="C32.1-Supraglottis",]
surv_object2=Surv(time=LarynxSCC_brief2$SurvivalMonths,event=LarynxSCC_brief2$VitalStatus_NEW)
fit12=survfit(surv_object2~Sex_NEW,data=LarynxSCC_brief2)
ggsurvplot(fit12, pval=TRUE,title = "Supraglottis", legend.title="Sex",legend.labs=c("Male","Female"),
           risk.table = TRUE, tables.theme = theme_cleantable())

# K-M curve for males vs. females with subglottic tumors
LarynxSCC_brief3=LarynxSCC_brief[LarynxSCC_brief$PrimarySite=="C32.2-Subglottis",]
surv_object3=Surv(time=LarynxSCC_brief3$SurvivalMonths,event=LarynxSCC_brief3$VitalStatus_NEW)
fit13=survfit(surv_object3~Sex_NEW,data=LarynxSCC_brief3)
ggsurvplot(fit13, pval=TRUE,title = "Supglottis", legend.title="Sex",legend.labs=c("Male","Female"),
           risk.table = TRUE, tables.theme = theme_cleantable())



############################### Mutiplicative and Additive Hazards Regression ######################



# Fit a cox-regression 
surv_object=Surv(time=LarynxSCC_brief$SurvivalMonths,event=LarynxSCC_brief$VitalStatus_NEW)

fit_cox=coxph(Surv(SurvivalMonths,VitalStatus_NEW==1)~Age+female+Race_NEW+Marital_NEW+
                PrimarySite+AJCC6thStage_NEW+TreatmentGrp+LNDx,
              data=LarynxSCC_brief)
summary(fit_cox)
ggforest(fit_cox,data=LarynxSCC_brief)

temp1=cox.zph(fit_cox,global=T)
print(temp1)

LarynxSCC_brief$resid_mart=residuals(fit_cox,type="martingale")
LarynxSCC_brief$resid_coxsnell=LarynxSCC_brief$resid_mart-LarynxSCC_brief$VitalStatus_NEW
fitcoxsnell=coxph(Surv(resid_coxsnell,VitalStatus_NEW)~1,data=LarynxSCC_brief)
df_base_haz <- basehaz(fitcoxsnell, centered = FALSE)



# Fit nonparametric additive hazards model 
LarynxSCC_brief$age.c=LarynxSCC_brief$Age-mean(LarynxSCC_brief$Age)


add1=aalen(surv_object~I(exp(age.c/10))+Sex_NEW+Race_NEW+Marital_NEW
           +PrimarySite+AJCC6thStage_NEW+TreatmentGrp+LNDx,
           data=LarynxSCC_brief,residuals = 1)
summary(add1)

# Fit semi-parametric additive hazards model 

add2=aalen(surv_object~const(I(exp(age.c/10)))+const(female)+const(Race_NEW)+const(Marital_NEW)
           +const(PrimarySite)+AJCC6thStage_NEW+const(TreatmentGrp)+const(LNDx),
           data=LarynxSCC_brief)
summary(add2)


###################################### PS Matching #################################### 

library(MatchIt)
library("optmatch")
library(cobalt)


# full matching # 
m.out5=matchit(female~Age+TreatmentGrp+Race_NEW+Marital_NEW+LNDx+PrimarySite+AJCC6thStage_NEW, 
               data=LarynxSCC_brief, method="full")
summary(m.out5,standardize = T)


m.data_full=match.data(m.out5,group="all",distance = "distance",weights="weights",subclass = "subclass")

# Check balance by graphical evidence
bal.plot(m.out5,var.name = "distance",which="both",type="histogram")+
  scale_fill_discrete(name="Sex",labels=c("Male","Female"))+
  ggtitle("")+xlab("Propensity Score")


####################### Semiparametric Additive Hazards Regression on Matched Data #######################


m.data_full$female=ifelse(m.data_full$Sex_NEW=="female",1,0)

m.data_full$age.c=m.data_full$Age-mean(m.data_full$Age)

surv_full=Surv(time=m.data_full$SurvivalMonths,event=m.data_full$VitalStatus_NEW)

fit_cox=coxph(Surv(SurvivalMonths,VitalStatus_NEW==1)~Age+female+Race_NEW+Marital_NEW+
                PrimarySite+AJCC6thStage_NEW+TreatmentGrp+LNDx,
              data=m.data_full,weights=m.data_full$weights)
summary(fit_cox)


temp1=cox.zph(fit_cox,global=T)
print(temp1)



add3=aalen(surv_full~I(exp(age.c/10))+female+Race_NEW+Marital_NEW+
             TreatmentGrp+LNDx+AJCC6thStage_NEW+PrimarySite,
           data=m.data_full, weights=m.data_full$weights,residuals = 1)

summary(add3)

plot(add3,specific.comps = 3)

res1=cum.residuals(add3,data=m.data_full,
                   modelmatrix=model.matrix(~-1+factor(cut(exp(age.c/10),4)),m.data_full))
plot(res1)
summary(res1)


resids=cum.residuals(add3,m.data_full,cum.resid = 1,n.sim=50) 
summary(resids)



add4=aalen(surv_full~const(I(exp(age.c/10)))+const(female)+const(Race_NEW)+const(Marital_NEW)+
             const(TreatmentGrp)+const(LNDx)+AJCC6thStage_NEW+const(PrimarySite),
           data=m.data_full, weights=m.data_full$weights,residuals = 1)

summary(add4)



