
# 03/21/2021
# JL
# Survival
# 10/14/2022
 

rm(list = ls(all = TRUE))

library(ggplot2)
library(dplyr)
library(survival)
library(surviminer)
library(readxl)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(Rcpp)
library(stringr)
library(reshape)
library(data.table)
library(ggsurvplot)

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rlang")
install.packages("rlang")

#if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/survminer")
install.packages("survminer")
library("survminer")


library(survival)
library(rpart)
library(sas7bdat)
library(Hmisc)
library(plotrix)

options(stringsAsFactors = FALSE)
setwd("C:/")

surv = read_excel("TCGA Colon and Rectal Cancer_survival.xlsx")
rnaseq = read.csv('tcga_sox9_phen_surv_RNAseq.csv')
surv_phen = read.csv("tcga_sur_phen_sox9all.csv")


set.seed(4)



########################################################### 



library(rpart)
install.packages("OIsurv")
library(OIsurv)

#removed<-c(2,20,11,32,45,69,95,97,100,125,130,132,135,160,174,181,183,184,186,190,197,198)
#dat2<-dat1[!(dat1$No %in% removed), ]

removed <- c(0)
surv_phen <- surv_phen[!(surv_phen$OS.time %>% removed),]
surv_phen = filter(surv_phen, surv_phen$OS.time != 0) 

tfit = rpart(formula = Surv(as.numeric(surv_phen$OS.time), as.numeric(surv_phen$OS))~SOX9, 
             data = surv_phen, cp = 0.001)
print(tfit)


plot(tfit, uniform = TRUE, branch = 1, compress = TRUE)
text(tfit, use.n = TRUE)



############################################################# 



install.packages("devtools") # (if not already installed)
library(devtools)


res.cut <- surv_cutpoint(surv_phen, time = "OS.time", event = "OS", variables = c("SOX9"))
                         
res.cut <- surv_cutpoint(surv_phen, time = "DSS.time", event = "DSS", variables = c("SOX9"))
                        
res.cut <- surv_cutpoint(surv_phen, time = "DFI.time", event = "DFI", variables = c("SOX9"))
                        

summary(res.cut)



# 2. Plot cutpoint for sox9
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(res.cut, "SOX9", palette = "npg")
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# 4. Fit survival curves and visualize
library("survival")

fit <- survfit(Surv(as.numeric(surv_phen$OS.time/30), as.numeric(surv_phen$OS)) ~ SOX9, data = res.cat)
fit <- survfit(Surv(as.numeric(surv_phen$DSS.time/30), as.numeric(surv_phen$DSS)) ~ SOX9, data = res.cat)
fit <- survfit(Surv(as.numeric(surv_phen$DFI.time/30), as.numeric(surv_phen$DFI)) ~ SOX9, data = res.cat)

ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE)

os_p <- ggsurvplot(
  fit,                  # survfit object with calculated statistics.
  size = 1.25,
  linetype = 1,
  data = res.cat,       # data used to fit survival curves.
  risk.table = FALSE,         # show risk table.
  pval = TRUE,                # show p-value of log-rank test.
  conf.int = FALSE,           # show confidence intervals for 
  pval.coord = c(0.15,0.15),
  pval.size = 9,
  # pval.text = element_text(color="black", size=24,  angle = 60),
  # palette = c("#E7B800", "#2E9FDF"),
  palette = c("red", "blue"),
  xlim = c(0,150),            # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time (Month)",      # customize X axis label.
  break.time.by = 30,         # break X axis in time intervals by 500.
  # ggtheme = theme_light(),  # customize plot and risk table with a theme.
  risk.table.pos = "out",
  risk.table.col = "black",
  risk.table.y.text = FALSE,
  # risk.table.y.text.col=TRUE,
  tables.theme = theme_cleantable(),
  font.tickslab = c(16),
  font.labs = "black",
  ggtheme = theme_classic2(base_size = 16, base_family = "Arial"),
  font.family = "Arial",
  # ggtheme = theme_classic(),
  # risk.table.y.text.col = T, # colour risk table text annotations.
  # risk.table.height = 0.25,  # the height of the risk table
  # risk.table.y.text = FALSE, # show bars instead of names in text annotations
  # in legend of risk table.
  # ncensor.plot = TRUE,       # plot the number of censored subjects at time t
  # ncensor.plot.height = 0.25,
  conf.int.style = "step",     # customize style of confidence intervals
  # surv.median.line = "hv",   # add the median survival pointer.
  legend.title = "",
  legend.labs = c("High SOX9", "Low SOX9"), # change legend labels.
  legend.coord = c(3,3),
  axis.ticks = element_line(size = 2),
  axis.line = element_line(size = 2),
  axis.labs = "black",
  axis.text.labs =  "black",
  cumevents.text.col = "black" 
) 

ggpar(os_p, 
      font.main = c(26, "bold","black"),
      font.x = c(28, "black"),
      font.y = c(28, "black"),
      font.caption = c(26, "bold","black"), 
      font.legend = c(20, "black"), 
      font.tickslab = c(24, "black"))
# jpeg('km_aneu.jpeg',width = 10, height = 7.5, units = 'in',res=300)

os_p
dev.off()

surv_phen$sox9Ep = ""

for(i in 1:length(surv_phen$sampleID))
{
  if (surv_phen$SOX9[i] > 11.66)
  {
    surv_phen$sox9Ep[i] = "high"
  }
  else 
  {
    surv_phen$sox9Ep[i] = "low"
  }
}


surv_phen9 <- surv_phen[order(as.numeric(surv_phen$SOX9)),]

# subset of first quartile and the bottom quartile of patients
surv_phena = surv_phen9 [57:376,] 
surv_phena = mutate(surv_phena, pecentage = 1)
surv_phenb = surv_phen9[1:56,] 
surv_phenb = mutate(surv_phenb, pecentage = 0)

surv_phen9 = rbind(surv_phena, surv_phenb)

os <- Surv(time = as.numeric(surv_phen9$DSS.time/30), event = as.numeric(surv_phen9$DSS))
fit_km <- survfit(os ~ pecentage, data = surv_phen9)

# summary(fit_km)

plot(fit_km)

os_p <- ggsurvplot(
  fit_km,                     # survfit object with calculated statistics.
  data = surv_phen9,          # survfit object with calculated statistics.
  size = 1.25,
  linetype = 1,
  risk.table = FALSE,         # show risk table.
  pval = TRUE,                # show p-value of log-rank test.
  conf.int = FALSE,           # show confidence intervals for 
  pval.coord = c(0.15,0.15),
  pval.size = 9,
  # pval.text = element_text(color="black", size=24,  angle = 60),
  # palette = c("#E7B800", "#2E9FDF"),
  palette = c("red", "blue"),
  xlim = c(0,150),            # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time (Month)",      # customize X axis label.
  break.time.by = 30,         # break X axis in time intervals by 500.
  # ggtheme = theme_light(),  # customize plot and risk table with a theme.
  risk.table.pos = "out",
  risk.table.col = "black",
  risk.table.y.text = FALSE,
  # risk.table.y.text.col=TRUE,
  tables.theme = theme_cleantable(),
  font.tickslab = c(16),
  font.labs = "black",
  ggtheme = theme_classic2(base_size = 16, base_family = "Arial"),
  font.family = "Arial",
  # ggtheme = theme_classic(),
  # risk.table.y.text.col = T, # colour risk table text annotations.
  # risk.table.height = 0.25,  # the height of the risk table
  # risk.table.y.text = FALSE, # show bars instead of names in text annotations
  # in legend of risk table.
  # ncensor.plot = TRUE,       # plot the number of censored subjects at time t
  # ncensor.plot.height = 0.25,
  conf.int.style = "step",     # customize style of confidence intervals
  # surv.median.line = "hv",   # add the median survival pointer.
  legend.title = "",
  legend.labs = c("Low SOX9", "High SOX9"), # change legend labels.
  legend.coord = c(3,3),
  axis.ticks = element_line(size = 2),
  axis.line = element_line(size = 2),
  axis.labs = "black",
  axis.text.labs =  "black",
  cumevents.text.col = "black" 
) 

ggpar(os_p, 
      font.main = c(26, "bold","black"),
      font.x = c(28, "black"),
      font.y = c(28, "black"),
      font.caption = c(26, "bold","black"), 
      font.legend = c(20, "black"), 
      font.tickslab = c(24, "black"))
# jpeg('km_aneu.jpeg',width = 10, height = 7.5, units = 'in',res=300)

os_p
dev.off()


###################################################################

# sox9 spline curve

org_par <- par(xaxs = "i", ask = TRUE)
library(survival)
library(rms)
library(dplyr)
library(Gmisc)
library( htmlTable)

 
install.packages('ggforestplot')
library(ggforestplot)
library(rms) 

surv_phen = read.csv("surv_phen.csv")
surv_phen9 <- surv_phen[order(as.numeric(surv_phen$SOX9)),]

# subset of first quartile and the bottom quartile of patients

surv_phena = surv_phen9 [57:376,] 
surv_phena = mutate(surv_phena, SOX9 = 'High')
surv_phenb = surv_phen9[1:56,] 
surv_phenb = mutate(surv_phenb, SOX9 = 'Low')

surv_phen9 = rbind(surv_phena, surv_phenb)

surv_phen9$Lymphatic_invasion <- factor(surv_phen9$Lymphatic_invasion, levels = c("Yes", "No"))
surv_phen9$Lymphatic_invasion = relevel(surv_phen9$Lymphatic_invasion, ref = "NO")

os <- Surv(time = surv_phen9$OS.time/30,event = surv_phen9$OS)
cox <- coxph(os ~ `Age`  + `Sex` + `SOX9`, data = surv_phen9)
cox <- coxph(os ~ `SOX9` + `Lymphatic_invasion`, data = surv_phen9)

summary(cox)
cox.zph(cox)

# plot(cox.zph(cox))
# ggforestplot::forestplot(cox,data=tcga_cox,fontsize = 1.8, noDigits = 4)

cox_plot <- ggforest(cox, data = surv_phen9, fontsize = 2, noDigits = 4)
forest_model(cox)

# forestplot(cox)
# ('cox_aneu_hyper_rrmm_bm2.#',width = 10, height = 7.5)

cox_plot


######################################################### 


dd <- datadist(surv_phen9$SOX9)
options(datadist='dd')

f <- cph(os ~ `SOX9`, data = surv_phen9,x = TRUE, y = TRUE)
cox.zph(f,"rank") # tests of PH
anova(f)
plot(Predict(f,`SOX9`))

########################################################### 

temp <- cox.zph(cox) 
plot(temp, var=4)
plot(temp[4])
abline(0, 0, lty=4)
abline(lm(temp$y[,4] ~ temp$x)$coefficients, lty=4, col=3)  
title(main="SOX9 spline") 

########################################################### 

library(splines)
library(Greg)
fit.cph <- cph(os ~ SOX9 , data = surv_phen9, x = TRUE, y = TRUE)

termplot(cox)
os <- Surv(time = surv_phen9$OS.time/30,event = surv_phen9$OS)
cox <- coxph(os ~ `Age`  + `Gender` + `Lymphatic_invasion` + bs(`SOX9`,4), data = surv_phen9)
cox <- coxph(os ~ `Age`  + `Gender` + bs(`SOX9`,3), data = surv_phen9)
cox <- coxph(os ~  bs(`SOX9`,4), data = surv_phen9)

plotHR(cox, term = "Age", plot.bty = "o", xlim = c(30, 70), xlab = "Age")
plotHR(cox, term = "SOX9", plot.bty = "o", xlim = c(5, 20), xlab = "SOX9")

########################################################## 

library(pspline)
library(survival)
library(Greg)

## create the survival object
os <- Surv(time = surv_phen$OS.time/30,event = surv_phen$OS )
ds <- Surv(time = surv_phen$DSS.time/30,event = surv_phen$DSS )
os <- Surv(time = surv_phen$DFI.time/30,event = surv_phen$DFI )


cox <- coxph(os ~ pspline(`SOX9`,df = 3), data = surv_phen)
cox <- coxph(os ~ pspline(`SOX9`,df = 3) +`Pathologic_stage` , data = surv_phen)
cox <- coxph(os ~ pspline(`SOX9`,df = 3) + `Age`, data = surv_phen)
cox <- coxph(os ~ pspline(`SOX9`,df = 3) + `Sex` + `Age`, data = surv_phen)

cox <- coxph(os ~ pspline(`SOX9`, df = 3) + `Age` + `Pathologic_stage` , data = surv_phen)
cox <- coxph(ds ~ pspline(`SOX9`,df = 3) + `Age`+ `Sex`+ `Pathologic_stage`, data = surv_phen)

## fit survival model
 

## get predicted values for fitted spline
predicted <- predict(cox , type = "terms" , se.fit = TRUE , terms = 4)

## plot lines
plotHR(cox, term = "SOX9", plot.bty = "o", xlim = c(5, 20), ylim = 2^c(-3,4),rug = "" , xlab = "Log2 SOX9 gene expression", ylab = "Hazard ratio")
plotHR(cox, term = "SOX9", plot.bty = "o", xlim = c(5, 20), xlab = "SOX9",  ylog=FALSE)
plotHR(cox, terms = "SOX9" , se = TRUE , rug = "" , x.log = FALSE ,   ylab = "Hazard Ratio" ,
       main = NULL ,  xlim = c(5, 20), ylim = 2^c(-3,3),  col.term = "#08519C",  col.se = "#DEEBF7", 
       cex = 1 , bty = "n" , axes = TRUE, xlab = "Log2 SOX9 gene expression" )
abline( h =0 , lty = 2 )
#abline(h = 0)

termplot(cox, term= 1, se=TRUE, col.term=2, col.se=1, xlab="Log2 SOX9 gene expression", ylab="Log relative death rate")
termplot(cox, term= 1, se=TRUE, col.term=2, col.se=1)

ptemp = termplot(cox, se=TRUE, plot=FALSE)

ageterm <- ptemp$Age # this will be a data frame
center <- with(ageterm, y[x == 50])
ytemp <- ageterm$y + outer(ageterm$se, c(0, -1.96, 1.96), '*')
matplot(ageterm$x, exp(ytemp - center), log='y',
        type='l', lty=c(1,2,2), col=1,
        xlab="Age at diagnosis", ylab ="Relative death rate")


sox9term <- ptemp$SOX9 # this will be a data frame
center <- with(sox9term, y[x == 12])
ytemp <- sox9term$y + outer(sox9term$se, c(0, -1.96, 1.96), '*')
matplot(sox9term$x, exp(ytemp - center), log='y',
        type='l', lty=c(1,2,2), col.term=2 ,
        xlab="Log2 SOX9 gene expression", ylab="Hazard ratio")

############################################################## 

plotHR(cox, term = "Age", plot.bty = "o", xlim = c(20, 100), ylim = 2^c(-3,3), xlab = "Age")

############################################################### 

##Scaled Schoenfeld Residuals Chart

##Schoenfeld residuals are used to test the assumption of proportional hazards. Schoenfeld residuals 
##"can essentially be thought of as the observed minus the expected values of the covariates at each failure time"
##(Steffensmeier & Jones, 2004: 121). There is a Schoenfeld residual for each subject for each covariate. The
##image below shows the SSR for the RACECAUCASIAN coefficient. The plot of Schoenfeld residuals against time 
##for any covariate should not show a pattern of changing residuals for that covariate. If there is a pattern, 
##that covariate is time-dependent. As a rule of thumb, a non-zero slope is an indication of a violation of the
##proportional hazard assumption. The dotted lines outline the 95% confidence interval.

################################################################ 

#sas<- read.sas7bdat("Y:/ORAL/EPOC_IP/data/final_short06feb2019more.sas7bdat")
#names(sas)
#sas_tp <- subset(sas, time_CFS > 0)

#=== Survival Fit:
#tfitc <- rpart(Surv(time_CFS, event_CFS) ~ HISTBASELINER3_2 + PDL1_pct_mean1, cp=0.001, data=sas_tp)
#plot(tfitc, uniform = TRUE, branch = 1, compress = TRUE)
#text(tfitc, use.n = TRUE)

#tfit2c <- prune(tfitc, cp = 0.001)
#par(mfrow=c(1,2), mar = rep(3,4), xpd=NA)
#plot(tfit2c, uniform = TRUE, branch = 0.7, compress = TRUE)
#text(tfit2c, use.n = TRUE)

############################################################## 

n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('Male','Female'), n, 
                     rep=TRUE, prob=c(.6, .4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex == 'Female'))
dt <- -log(runif(n))/h
label(dt) <- 'Follow-up Time'
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
units(dt) <- "Year"
dd <- datadist(age, sex)
options(datadist='dd')
S <- Surv(dt,e)

f <- cph(S ~ rcs(age,4) + sex, x=TRUE, y=TRUE)
cox.zph(f, "rank") # tests of PH
anova(f)
plot(Predict(f, age, sex)) # plot age effect, 2 curves for 2 sexes

################################################# 



#Compare the normal survival with SOX9 survival

surv_phen <- surv_phen[order(as.numeric(surv_phen$SOX9)),]

# subset of first quartile and the bottom quartile of patients
surv_phena = surv_phen [57:376,] 
surv_phena = mutate(surv_phena, pecentage = 1)
surv_phenb = surv_phen[1:56,] 
surv_phenb = mutate(surv_phenb, pecentage = 0)

surv_phen9 = rbind(surv_phena, surv_phenb)

################################################

normal = read.csv('tcga_sox9_RNAseq_tumor_normal.csv')

normal = filter(normal, sample_type == "Normal Tissue")
colnames(surv)[1] = "sampleID"

surv_phenN = merge(normal, surv, by = "sampleID")

surv_phen1 = surv_phenN[, -(2:4)]
surv_phenN1 = mutate(surv_phen1, pecentage = 0)


surv_phena1 =surv_phena[, -(10:21)]


surv_phen9 = rbind(surv_phena1, surv_phenN1)

 

##############################################





