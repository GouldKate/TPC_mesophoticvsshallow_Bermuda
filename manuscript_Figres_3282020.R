library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggpubr)

theme_set(theme_sjplot())

All.CIT<-read.csv("C:/Github/TPC_BDA_github/R_output/Final/October_RmarkdownFinal/All.CI.csv")
All.fittedT<-read.csv("C:/Github/TPC_BDA_github/R_output/Final/October_RmarkdownFinal/All.fitted.csv")

All.fittedT%<>%
  mutate(group=factor(paste(treatment,rate.type)))

All.CIT%<>%
  mutate(group=factor(paste(treatment,rate.type)))
### TPC
STPCall<-ggplot() +
  geom_ribbon(data=subset(All.CIT, lwr_CI>0),
              aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedT, aes(x=(K - 273.15), y=log.rate, shape=rate.type)) +
  geom_smooth(data=All.fittedT, aes(x=(K - 273.15), y=fitted, colour=group,se=F),size=1.25) +
  ylab('Rate()') +
  xlab('Temperature (ºC)') +
  theme_classic()+
  scale_color_manual(breaks = c("deepGP","deepR", "shallowGP", "shallowR"), values= c("green4","#003366", "green1","#00CCFF"))+
  facet_wrap(~species)

STPCall

ggexport(STPCall, filename = "TPC_final.pdf")

#####
##### Significant parameter estimates
res.aov3 <- aov(lnc~ species, data = testGPlnc)

summary(res.aov3) 
TukeyHSD(res.aov3)

model.tables(res.aov3, type="means", se= TRUE)

sumlncGP = summarySE(testGPlnc,
                  measurevar="lnc",            groupvars=c("species"))

sumlncGP

testGPPo<- na.omit(testGPP) #remove NaNs for means


res.aov3 <- aov(Pmax~ species, data = testGPPo)

summary(res.aov3) 
TukeyHSD(res.aov3)

model.tables(res.aov3, type="means", se= TRUE)


sumPmaxGP = summarySE(testGPPo,
                  measurevar="Pmax",            groupvars=c("species"))

sumPmaxGP


#MCAV-DLAB -0.32391599 -0.45858743 -0.18924455 0.0000127
#OFRA-DLAB -0.06566590 -0.20033735  0.06900554 0.5281746
#PAST-DLAB -0.26008747 -0.40074720 -0.11942773 0.0003054
#OFRA-MCAV  0.25825009  0.12984590  0.38665427 0.0001172
#PAST-MCAV  0.06382853 -0.07084292  0.19849997 0.5510111
#PAST-OFRA -0.19442156 -0.32909301 -0.05975012 0.0035731


### Respiration Pmax
testRo<- na.omit(testR) #remove NaNs for means


res.aov3 <- aov(Pmax~ species, data = testRo)

summary(res.aov3) 
TukeyHSD(res.aov3)

model.tables(res.aov3, type="means", se= TRUE)


sumPmaxGP = summarySE(testRo,
                      measurevar="Pmax",            groupvars=c("species"))

sumPmaxGP



Fitwinlm<-lm(lnc~species, testGPlnc)
Fitwinlm
summary(Fitwinlm)
anova(Fitwinlm)


################
logGP<- ggplot(mydataP, aes(x=temp.Cat, log.rate, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, shape=species, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  ggtitle("")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4","green1"))+
  facet_wrap("species")

logGP

#respiration
logR<- ggplot(mydataR, aes(x=temp.Cat, log.rate, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, shape=species, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  ggtitle("")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))+
  facet_wrap("species")

logR


ggexport(logGP, filename = "fig2GP.pdf")
ggexport(logR, filename = "fig2r.pdf")
#
library(ggpubr)
#coefficient plots of diff between depths at specific tempreatures
#Gross photo diff- DLAB, MCAv, OFRA 19,19,34

GPdiffDMO<-ggpubr::ggarrange(gGP19DP, G19MPP, G34OPP + rremove("x.text"),
                           ncol = 2, nrow = 3)

GPdiffDMO
ggexport(GPdiffDMO, filename = "GPdiffDMO.pdf")


#MCav resp- PAST respi 19 32


RMPdepDiff<-ggpubr::ggarrange(R23DRP,R19MRP,R32PRP + rremove("x.text"),
                           ncol = 2, nrow = 3)

RMPdepDiff
ggexport(RMPdepDiff, filename = "RMPdepDiff.pdf")




###################### Figure 3. Depth Coefficient curves for E, Eh, Th, Topt, lnc, Pmax
###################### GP E
E.dep <- lm(E~treatment,testGPE)
summary(E.dep)
Anova(E.dep)
summ(E.dep, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coefEd <- as.data.frame(summary(E.dep)$coefficients)
coefEd
names(coefEd)[2] = "se" 
coral <- c(treatmentshallow="")
coefEd$vars = coral
coefEd<-coefEd[-1,]

GpEd<-ggplot(coefEd, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003300", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#33FF33", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab("Energy (eV)")+
  xlab("")+
  ggtitle('E')
GpEd


#### Eh
Eh.dep <- lm(Eh~treatment,GPEh)
summary(Eh.dep)
Anova(Eh.dep)
summ(Eh.dep, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coefEhd <- as.data.frame(summary(Eh.dep)$coefficients)
coefEhd
names(coefEhd)[2] = "se" 
coral <- c(treatmentshallow="")
coefEhd$vars = coral
coefEhd<-coefEhd[-1,]

GpEhd<-ggplot(coefEhd, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003300", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#33FF33", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab("Energy (eV)")+
  xlab("")+
  ggtitle('Eh')
GpEhd

#DLAB#### Eh
DEh.dep <- lm(Eh~treatment,DtestGP)
summary(DEh.dep)
Anova(DEh.dep)
summ(DEh.dep, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
DcoefEhd <- as.data.frame(summary(DEh.dep)$coefficients)
DcoefEhd
names(DcoefEhd)[2] = "se" 
coral <- c(treatmentshallow="")
DcoefEhd$vars = coral
DcoefEhd<-DcoefEhd[-1,]

DGpEhd<-ggplot(DcoefEhd, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003300", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#33FF33", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab("Energy (eV)")+
  xlab("")+
  ggtitle('DLAB Eh')
DGpEhd
ggexport(DGpEhd, filename = "DGpEhd.pdf")

#####Th
Th.dep <- lm(Th_C~treatment,testGPTh_C)
summary(Th.dep)
Anova(Th.dep)
summ(Th.dep, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coefThd <- as.data.frame(summary(Th.dep)$coefficients)
coefThd
names(coefThd)[2] = "se" 
coral <- c(treatmentshallow="")
coefThd$vars = coral
coefThd<-coefThd[-1,]

GpThd<-ggplot(coefThd, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003300", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#33FF33", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab('Temperature (ºC)')+
  xlab("")+
  ggtitle('Th')
GpThd


### Topt

To.dep <- lm(topt_C~treatment,testGPt)
summary(To.dep)
Anova(To.dep)
summ(To.dep, confint = TRUE, digits = 1)

To.dep <- lm(topt_C~treatment,testGPt)

#Code for Coefficents without densities treat
coefTod <- as.data.frame(summary(To.dep)$coefficients)
coefTod
names(coefTod)[2] = "se" 
coral <- c(treatmentshallow="")
coefTod$vars = coral
coefTod<-coefTod[-1,]

GpTod<-ggplot(coefTod, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003300", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#33FF33", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab('Temperature (ºC)')+
  xlab("")+
  ggtitle('Topt')
GpTod


#gpubr::ggarrange(GpEd, GpEhd,GpThd,GpTod + rremove("x.text"),
  #                ncol = 2, nrow = 4)

#####Th
Pmax.dep <- lm(Pmax~treatment,testGPP)
summary(Pmax.dep)
Anova(Pmax.dep)
summ(Pmax.dep, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coefpmaxd <- as.data.frame(summary(Pmax.dep)$coefficients)
coefpmaxd
names(coefpmaxd)[2] = "se" 
coral <- c(treatmentshallow="")
coefpmaxd$vars = coral
coefpmaxd<-coefpmaxd[-1,]

Gppmaxd<-ggplot(coefpmaxd, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003300", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#33FF33", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  xlab("")+
  ggtitle('Pmax')
Gppmaxd

##### Summs Fpr P:R ratios
#####
sumpmaxP = summarySE(testGPP, na.rm=TRUE,measurevar="Pmax",
      groupvars=c("treatment"))

sumpmaxR = summarySE(testRP, na.rm=TRUE,measurevar="Pmax",
                    groupvars=c("treatment"))
sumpmaxP
sumpmaxR

sumpmaxPs = summarySE(testGP, na.rm=TRUE,measurevar="Pmax",
                     groupvars=c("treatment", "species"))

sumpmaxRs = summarySE(testR, na.rm=TRUE,measurevar="Pmax",
                     groupvars=c("treatment", "species"))

sumpmaxPs
sumpmaxRs

sumlncPs = summarySE(testGP, na.rm=TRUE,measurevar="lnc",
                      groupvars=c("treatment", "species"))

sumlncRs = summarySE(testR, na.rm=TRUE,measurevar="lnc",
                      groupvars=c("treatment", "species"))

sumlncPs
sumlncRs

### lnc
lnc.dep <- lm(lnc~treatment,testGPlnc
)
summary(lnc.dep)
Anova(lnc.dep)
summ(lnc.dep, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coeflncd <- as.data.frame(summary(lnc.dep)$coefficients)
coeflncd
names(coeflncd)[2] = "se" 
coral <- c(treatmentshallow="")
coeflncd$vars = coral
coeflncd<-coeflncd[-1,]

Gplncd<-ggplot(coeflncd, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003300", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#33FF33", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  xlab("")+
  ggtitle('lnc')
Gplncd

ggpubr::ggarrange(Gplncd,Gppmaxd + rremove("x.text"),
                  ncol = 2, nrow = 2)

Figure3<-ggpubr::ggarrange(GpEd, GpEhd,GpThd,GpTod, Gplncd,Gppmaxd + rremove("x.text"),
                  ncol = 2, nrow = 3)

Figure3
ggexport(Figure3, filename = "figure3.pdf")

####### Respiration 
#deprth Coefficient curves for E, Eh, Th, Topt, lnc, Pmax
###################### R E
E.depr <- lm(E~treatment,testR)
summary(E.depr)
Anova(E.depr)
summ(E.depr, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
rcoefEd <- as.data.frame(summary(E.depr)$coefficients)
rcoefEd
names(rcoefEd)[2] = "se" 
coral <- c(treatmentshallow="")
rcoefEd$vars = coral
rcoefEd<-rcoefEd[-1,]

REd<-ggplot(rcoefEd, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003366", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#00CCFF", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab("Energy (eV)")+
  xlab("")+
  ggtitle('E')
REd
ggexport(REd, filename = "Erespircoeff.pdf")
##### Box plot

ResE<-ggplot(E.depr, aes(x=treatment, E, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab('Energy (eV)')+
  ggtitle("E")+
  stat_compare_means(label.x= "deep", label.y=4.5, method = "anova")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF")) 
ResE

ggexport(ResE, filename = "resp_logEtreat.pdf")
## Spp PAST  has no different E across depth

#PastlogE<-subset(testR, species=="PAST")
#PLET<-lm(E~treatment,PastlogE)
#Anova(PLET)
#
#RlogEtrP<-ggplot(PastlogE, aes(x=treatment, E, fill= treatment))+
#  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
#  geom_boxplot(position = position_dodge(1))+
#  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
#  guides(fill= guide_legend(override.aes = list(shape= NA)))+
#  theme_classic()+
#  xlab(NULL)+
#  ylab('Energy (eV)')+
#  ggtitle("PAST E")+
#  stat_compare_means(label.x= "deep", label.y=7.5, method = "anova")+
#  scale_fill_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF")) 
#RlogEtrP
#
#ggexport(RlogEtrP, filename = "RlogEtrP.pdf")
#### Eh
Eh.depr <- lm(Eh~treatment,testREh)
summary(Eh.depr)
Anova(Eh.depr)
summ(Eh.depr, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
rcoefEhd <- as.data.frame(summary(Eh.depr)$coefficients)
rcoefEhd
names(rcoefEhd)[2] = "se" 
coral <- c(treatmentshallow="")
rcoefEhd$vars = coral
rcoefEhd<-rcoefEhd[-1,]

REhd<-ggplot(rcoefEhd, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003366", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#00CCFF", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab("Energy (eV)")+
  xlab("")+
  ggtitle('Eh')
REhd
#####Th
Th.depr <- lm(log.Th_C~treatment,testRTh_C)
summary(Th.depr)
Anova(Th.depr)
summ(Th.depr, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
rcoefThd <- as.data.frame(summary(Th.depr)$coefficients)
rcoefThd
names(rcoefThd)[2] = "se" 
coral <- c(treatmentshallow="")
rcoefThd$vars = coral
rcoefThd<-rcoefThd[-1,]

RThd<-ggplot(rcoefThd, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003366", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#00CCFF", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab('Temperature (ºC)')+
  xlab("")+
  ggtitle('log Th')
RThd

ggexport(RThd, filename = "RThd.pdf")

### Topt
To.depr <- lm(topt_C~treatment,testRt)
summary(To.depr)
Anova(To.depr)
summ(To.depr, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
rcoefTod <- as.data.frame(summary(To.depr)$coefficients)
rcoefTod
names(rcoefTod)[2] = "se" 
coral <- c(treatmentshallow="")
rcoefTod$vars = coral
rcoefTod<-rcoefTod[-1,]

RTod<-ggplot(rcoefTod, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003366", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#00CCFF", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab('Temperature (ºC)')+
  xlab("")+
  ggtitle('Topt')
RTod

#ggpubr::ggarrange(REd, REhd,RThd,RTod + rremove("x.text"),
#                  ncol = 2, nrow = 4)

#####Th
Pmax.depr <- lm(Pmax~treatment,testR)
summary(Pmax.depr)
Anova(Pmax.depr)
summ(Pmax.depr, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
rcoefpmaxd <- as.data.frame(summary(Pmax.depr)$coefficients)
rcoefpmaxd
names(rcoefpmaxd)[2] = "se" 
coral <- c(treatmentshallow="")
rcoefpmaxd$vars = coral
rcoefpmaxd<-rcoefpmaxd[-1,]

Rpmaxd<-ggplot(rcoefpmaxd, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003366", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#00CCFF", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  xlab("")+
  ggtitle('Pmax')
Rpmaxd


### lnc
lnc.depr <- lm(lnc~treatment, testRlnc)
summary(lnc.depr)
Anova(lnc.depr)

summ(lnc.depr, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
rcoeflncd <- as.data.frame(summary(lnc.depr)$coefficients)
rcoeflncd
names(rcoeflncd)[2] = "se" 
coral <- c(treatmentshallow="")
rcoeflncd$vars = coral
rcoeflncd<-rcoeflncd[-1,]

Rlncd<-ggplot(rcoeflncd, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003366", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#00CCFF", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  xlab("")+
  ggtitle('lnc')
Rlncd

ggexport(Rlncd, filename = "Rlncd.pdf")


Figure4<-ggpubr::ggarrange(REd, REhd,RThd,RTod, Rlncd,Rpmaxd + rremove("x.text"),
                           ncol = 2, nrow = 3)

Figure4
ggexport(Figure4, filename = "figure4.pdf")

###################### GP E
###################### 
######################  Gross photo coefficient plots for species 
E.depS <- lm(E~species,testGPEs)
summary(E.depS)
Anova(E.depS)
summ(E.depS, confint = TRUE, digits = 3)

coefEGPS <- as.data.frame(summary(E.depS)$coefficients[-1,1:2])
names(coefEGPS)[2] = "se" 
coral <- c(speciesMCAV="MCAV", speciesOFRA="OFRA",speciesPAST="PAST")
coefEGPS$vars = coral
#

EgpS<-ggplot(coefEGPS, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003300", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#33FF33", width=0) +
  geom_point(size=3, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab("Energy (eV)")+
  xlab("")+
  ggtitle('E')

EgpS

#### Eh
Eh.deps <- lm(Eh~species,testGPEh)
summary(Eh.deps)
Anova(Eh.deps)
summ(Eh.deps, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coefEhGPS <- as.data.frame(summary(Eh.deps)$coefficients[-1,1:2])
names(coefEhGPS)[2] = "se" 
coral <- c(speciesMCAV="MCAV", speciesOFRA="OFRA",speciesPAST="PAST")
coefEhGPS$vars = coral
#

EhgpS<-ggplot(coefEhGPS, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003300", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#33FF33", width=0) +
  geom_point(size=3, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab("Energy (eV)")+
  xlab("")+
  ggtitle('Eh')

EhgpS

#####Th
Th.deps <- lm(Th_C~species,testGPTh_C)
summary(Th.deps)
Anova(Th.deps)
summ(Th.deps, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coefThGPS <- as.data.frame(summary(Th.deps)$coefficients[-1,1:2])
names(coefThGPS)[2] = "se" 
coral <- c(speciesMCAV="MCAV", speciesOFRA="OFRA",speciesPAST="PAST")
coefThGPS$vars = coral
#

ThgpS<-ggplot(coefThGPS, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003300", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#33FF33", width=0) +
  geom_point(size=3, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab("Temperature (ºC)")+
  xlab("")+
  ggtitle('Th')

ThgpS


### Topt
To.deps <- lm(topt_C~species,testGPt)
summary(To.deps)
Anova(To.deps)
summ(To.deps, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coefToGPS <- as.data.frame(summary(To.deps)$coefficients[-1,1:2])
names(coefToGPS)[2] = "se" 
coral <- c(speciesMCAV="MCAV", speciesOFRA="OFRA",speciesPAST="PAST")
coefToGPS$vars = coral
#

TogpS<-ggplot(coefToGPS, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003300", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#33FF33", width=0) +
  geom_point(size=3, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab("Temperature (ºC)")+
  xlab("")+
  ggtitle('Topt')

TogpS


#####Th
Pm.deps <- lm(Pmax~species,testGPP)
summary(Pm.deps)
Anova(Pm.deps)
summ(Pm.deps, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coefPGPS <- as.data.frame(summary(Pm.deps)$coefficients[-1,1:2])
names(coefPGPS)[2] = "se" 
coral <- c(speciesMCAV="MCAV", speciesOFRA="OFRA",speciesPAST="PAST")
coefPGPS$vars = coral
#

PmgpS<-ggplot(coefPGPS, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003300", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#33FF33", width=0) +
  geom_point(size=3, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  xlab("")+
  ggtitle('Pmax')

PmgpS

### lnc
lnc.deps <- lm(lnc~species,testGPlnc)
summary(lnc.deps)
Anova(lnc.deps)


summ(lnc.deps, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coeflncGPS <- as.data.frame(summary(lnc.deps)$coefficients[-1,1:2])
names(coeflncGPS)[2] = "se" 
coral <- c(speciesMCAV="MCAV", speciesOFRA="OFRA",speciesPAST="PAST")
coeflncGPS$vars = coral
#

lncgpS<-ggplot(coeflncGPS, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003300", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#33FF33", width=0) +
  geom_point(size=3, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  xlab("")+
  ggtitle('lnc')

lncgpS


Figure5<-ggpubr::ggarrange(EgpS,EhgpS,ThgpS,TogpS,lncgpS,PmgpS + rremove("x.text"),
                           ncol = 2, nrow = 3)

Figure5
ggexport(Figure5, filename = "figure5.pdf")

##### coefficients for sppecies Respiration
##### 
###################### R E
###################### 
######################  Gross photo coefficient plots for species 
E.depSr <- lm(E~species,testREs)
summary(E.depSr)
Anova(E.depSr)
summ(E.depSr, confint = TRUE, digits = 3)

coefERS <- as.data.frame(summary(E.depSr)$coefficients[-1,1:2])
names(coefERS)[2] = "se" 
coral <- c(speciesMCAV="MCAV", speciesOFRA="OFRA",speciesPAST="PAST")
coefERS$vars = coral
s#


ERS<-ggplot(coefERS, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003366", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#00CCFF", width=0) +
  geom_point(size=3, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab("Energy (eV)")+
  xlab("")+
  ggtitle('E')

ERS

#### Eh
Eh.depsr <- lm(Eh~species,testREhs)
summary(Eh.depsr)
Anova(Eh.depsr)
summ(Eh.depsr, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coefEhRS <- as.data.frame(summary(Eh.depsr)$coefficients[-1,1:2])
names(coefEhRS)[2] = "se" 
coral <- c(speciesMCAV="MCAV", speciesOFRA="OFRA",speciesPAST="PAST")
coefEhRS$vars = coral
#

EhRS<-ggplot(coefEhRS, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003366", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#00CCFF", width=0) +
  geom_point(size=3, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab("Energy (eV)")+
  xlab("")+
  ggtitle('Eh')

EhRS

#####Th
Th.depsr <- lm(log.Th_C~species,testRTh_C)
summary(Th.depsr)
Anova(Th.depsr)
summ(Th.depsr, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coefThRS <- as.data.frame(summary(Th.depsr)$coefficients[-1,1:2])
names(coefThRS)[2] = "se" 
coral <- c(speciesMCAV="MCAV", speciesOFRA="OFRA",speciesPAST="PAST")
coefThRS$vars = coral
#

ThRS<-ggplot(coefThRS, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003366", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#00CCFF", width=0) +
  geom_point(size=3, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab("Temperature (ºC)")+
  xlab("")+
  ggtitle('Th')

ThRS


### Topt
To.depsr <- lm(topt_C~species,testRt)
summary(To.depsr)
Anova(To.depsr)
summ(To.depsr, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coefToRS <- as.data.frame(summary(To.depsr)$coefficients[-1,1:2])
names(coefToRS)[2] = "se" 
coral <- c(speciesMCAV="MCAV", speciesOFRA="OFRA",speciesPAST="PAST")
coefToRS$vars = coral
#

ToRS<-ggplot(coefToRS, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003366", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#00CCFF", width=0) +
  geom_point(size=3, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab("Temperature (ºC)")+
  xlab("")+
  ggtitle('Topt')

ToRS


#####Pmax
Pm.depsr <- lm(Pmax~species,testRP)
summary(Pm.depsr)
Anova(Pm.depsr)
summ(Pm.depsr, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coefPRS <- as.data.frame(summary(Pm.depsr)$coefficients[-1,1:2])
names(coefPRS)[2] = "se" 
coral <- c(speciesMCAV="MCAV", speciesOFRA="OFRA",speciesPAST="PAST")
coefPRS$vars = coral
#

PmRS<-ggplot(coefPRS, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003366", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#00CCFF", width=0) +
  geom_point(size=3, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  xlab("")+
  ggtitle('Pmax')

PmRS

### lnc
lnc.depsr <- lm(lnc~species,testRlncs)
summary(lnc.depsr)
Anova(lnc.depsr)
summ(lnc.depsr, confint = TRUE, digits = 1)

#Code for Coefficents without densities treat
coeflncRS <- as.data.frame(summary(lnc.depsr)$coefficients[-1,1:2])
names(coeflncRS)[2] = "se" 
coral <- c(speciesMCAV="MCAV", speciesOFRA="OFRA",speciesPAST="PAST")
coeflncRS$vars = coral
#

lncRS<-ggplot(coeflncRS, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003366", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#00CCFF", width=0) +
  geom_point(size=3, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  xlab("")+
  ggtitle('lnc')

lncRS


Figure6<-ggpubr::ggarrange(ERS,EhRS,ThRS,ToRS,lncRS,PmRS + rremove("x.text"),
                           ncol = 2, nrow = 3)

Figure6
ggexport(Figure6, filename = "figure6.pdf")



#photosysnthesis parameters in Fig. 5 and respiration parameters in Fig. 6. 

# GP E, Eh,lnc, Th, Topt, Pmax
# E
E<-lm(E~species,testGPE)
summary(E)
Anova(E)
ET<-lm(E~treatment,testGPE)
Anova(ET)

ETs<-aov(E~species*treatment,testGPE)
summary(ETs)
TukeyHSD(ETs)

GPEspp<-ggplot(testGPE, aes(x=species, E, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab('Energy (eV)')+
  ggtitle("E")+
  stat_compare_means(label.x= "MCAV", label.y=1, method = "anova")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
GPEspp

  

GPEtr<-ggplot(testGPE, aes(x=species, E, fill= species))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(shape= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab('Energy (eV)')+
  ggtitle("E Photosynthesis")+
  stat_compare_means(label.x= "MCAV", label.y=1, method = "anova")+
  scale_fill_manual(breaks = c("DLAB", "MCAV", "OFRA", "PAST"), values= c( "#3399FF", "#FF9900", "#339900","#CC0000"))
GPEtr

ggexport(GPEspp, filename = "GPEspp.pdf")
ggexport(GPEtr, filename = "GPEtr.pdf")

## NExt parameters below
## 
## GP Eh
#DLAB is SIG BETWEEN DEPTHS:
# Eh
Eh<-lm(Eh~species,testGPEh)
summary(Eh)
Anova(Eh)
EhT<-lm(Eh~treatment,testGPEh)
Anova(EhT)

EhTs<-aov(Eh~species*treatment,testGPEh)
summary(EhTs)
TukeyHSD(EhTs)

GPEhtr<-ggplot(testGPEh, aes(x=species, Eh, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab('Energy (eV)')+
  ggtitle("Eh")+
  stat_compare_means(label.x= "MCAV", label.y=7, method = "t.test")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
GPEhtr

GPEhspp<-ggplot(testGPEh, aes(x=species, Eh, fill= species))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(shape= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab('Energy (eV)')+
  ggtitle("Eh Photosynthesis")+
  stat_compare_means(label.x= "MCAV", label.y=7, method = "anova")+
  scale_fill_manual(breaks = c("DLAB", "MCAV", "OFRA", "PAST"), values= c( "#3399FF", "#FF9900", "#339900","#CC0000"))
GPEhspp


ggexport(GPEhtr, filename = "GPEhspp.pdf")
ggexport(GPEhspp, filename = "GPEhtr.pdf")

## Spp DLAB has different Eh across depth

Dlabeh<-subset(testGPEh, species=="DLAB")
DEhT<-lm(Eh~treatment,Dlabeh)
Anova(DEhT)

GPEhtrD<-ggplot(Dlabeh, aes(x=treatment, Eh, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab('Energy (eV)')+
  ggtitle("DLAB Eh")+
  stat_compare_means(label.x= "deep", label.y=7.5, method = "anova")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
GPEhtrD

ggexport(GPEhtrD, filename = "DGPEtr.pdf")

## GP
## lnc
# lnc
lncTs<-aov(lnc~species*treatment,testGPlnc)
summary(lncTs)
TukeyHSD(lncTs)
lncspp<-aov(lnc~species,testGPlnc)
TukeyHSD(lncspp)

lnc<-lm(lnc~species,testGPlnc)
summary(lnc)
Anova(lnc)
lncT<-lm(lnc~treatment,testGPlnc)
Anova(lncT)


GPlnctr<-ggplot(testGPlnc, aes(x=species, lnc, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  ggtitle("lnc")+
  stat_compare_means(label.x= "MCAV", label.y=2, method = "anova")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
GPlnctr


GPlncspp<-ggplot(testGPlnc, aes(x=species, lnc, fill= species))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(shape= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  ggtitle("lnc Photosynthesis")+
  stat_compare_means(label.x= "MCAV", label.y=2, method = "anova")+
  scale_fill_manual(breaks = c("DLAB", "MCAV", "OFRA", "PAST"), values= c( "#3399FF", "#FF9900", "#339900","#CC0000"))
GPlncspp


ggexport(GPlnctr, filename = "GPlncspp.pdf")
ggexport(GPlncspp, filename = "GPlnctr.pdf")


### Checkign DLAB... 
#Dlablnc<-subset(testGP, species=="DLAB")
#DlncT<-lm(lnc~treatment,Dlablnc)
#Anova(DlncT)
#leveneTest(lnc ~ treatment, data = Dlablnc)
#res.Dlnc <- aov(lnc ~ treatment, data = Dlablnc)
#aov_residualDLNC <- residuals(object = res.Dlnc)
## Run Shapiro-Wilk test
#shapiro.test(x = aov_residualDLNC )
#
#GPlnctGPD<-ggplot(Dlablnc, aes(x=treatment, lnc, fill= treatment))+
#  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
#  geom_boxplot(position = position_dodge(1))+
#  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
#  guides(fill= guide_legend(override.aes = list(shape= NA)))+
#  theme_classic()+
#  xlab(NULL)+
#  ylab('log rate')+
#  ggtitle("DLAB lnc")+
#  stat_compare_means(label.x= "deep", label.y=1.8, method = "anova")+
#  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
#GPlnctGPD

## Th
# Th_C
Th_CTs<-aov(Th_C~species*treatment,testGPTh_C)

summary(Th_CTs)

Th_C<-lm(Th_C~species,testGPTh_C)
summary(Th_C)

Th_CT<-lm(Th_C~treatment,testGPTh_C)
Anova(Th_CT)

GPTh_Ctr<-ggplot(testGPTh_C, aes(x=species, Th_C, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab("Temperature (ºC)")+
  ggtitle("Th_C")+
  stat_compare_means(label.x= "MCAV", label.y=37, method = "anova")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
GPTh_Ctr

GPTh_Cspp<-ggplot(testGPTh_C, aes(x=species, Th_C))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(shape= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab("Temperature (ºC)")+
  ggtitle("Th")+
  stat_compare_means(label.x= "MCAV", label.y=37.5, method = "anova")
GPTh_Cspp


ggexport(GPTh_Ctr, filename = "GPTh_Cspp.pdf")
ggexport(GPTh_Cspp, filename = "GPTh_Ctr.pdf")

Thbwgpspp<-ggplot(testGPTh_C, aes(x=species, Th_C))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(shape= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab("Temperature (ºC)")+
  ggtitle("Th")+
  stat_compare_means(label.x= "MCAV", label.y=37.5, method = "anova")

ggexport(Thbwgpspp, filename = "Thbwgpsppr.pdf")

### Checkign DLAB... 
#DlabTh<-subset(testGP, species=="DLAB")
#DlncTHC<-lm(Th_C~treatment,DlabTh)
#Anova(DlncTHC)
#leveneTest(Th_C ~ treatment, data = DlabTh)
#res.Dlnc <- aov(Th_C ~ treatment, data = DlabTh)
#aov_residualDLNC <- residuals(object = res.Dlnc)
## Run Shapiro-Wilk test
#shapiro.test(x = aov_residualDLNC )
#
#GPThtGPD<-ggplot(DlabTh, aes(x=treatment, Th_C, fill= treatment))+
#  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
#  geom_boxplot(position = position_dodge(1))+
#  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
#  guides(fill= guide_legend(override.aes = list(shape= NA)))+
#  theme_classic()+
#  xlab(NULL)+
#  ylab('log rate')+
#  ggtitle("DLAB Th_C")+
#  stat_compare_means(label.x= "deep", label.y=38, method = "anova")+
#  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
#GPThtGPD
##### Topt
# topt_C
topt_CTs<-aov(topt_C~species*treatment,testGPt)

summary(topt_CTs)

topt_C<-lm(topt_C~species,testGPt)
summary(topt_C)

topt_CT<-lm(topt_C~treatment,testGPt)
Anova(topt_CT)

GPtopt_Ctr<-ggplot(testGPt, aes(x=species, topt_C, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab("Temperature (ºC)")+
  ggtitle("Topt")+
  stat_compare_means(label.x= "MCAV", label.y=37, method = "anova")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
GPtopt_Ctr

GPtopt_Cspp<-ggplot(testGPt, aes(x=species, topt_C, fill= species))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(shape= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab("Temperature (ºC)")+
  ggtitle("Topt Photosynthesis")+
  stat_compare_means(label.x= "MCAV", label.y=37, method = "anova")+
  scale_fill_manual(breaks = c("DLAB", "MCAV", "OFRA", "PAST"), values= c( "#3399FF", "#FF9900", "#339900","#CC0000"))
GPtopt_Cspp

ggexport(GPtopt_Ctr, filename = "GPtopt_Cspp.pdf")
ggexport(GPtopt_Cspp, filename = "GPtopt_Ctr.pdf")
### Checkign DLAB... 
#Dlabtop<-subset(testGPt, species=="DLAB")
#DlabPto<-lm(topt_C~treatment,Dlabtop)
#Anova(DlabPto)
#leveneTest(topt_C ~ treatment, data = Dlabtop)
#res.Dto <- aov(topt_C ~ treatment, data = Dlabtop)
#aov_residualto <- residuals(object = res.Dto)
## Run Shapiro-Wilk test
#shapiro.test(x = aov_residualto )
#
#GPThtGPD<-ggplot(Dlabtop, aes(x=treatment, topt_C, fill= treatment))+
#  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
#  geom_boxplot(position = position_dodge(1))+
#  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
#  guides(fill= guide_legend(override.aes = list(shape= NA)))+
#  theme_classic()+
#  xlab(NULL)+
#  ylab('log rate')+
#  ggtitle("DLAB topt_C")+
#  stat_compare_means(label.x= "deep", label.y=33, method = "anova")+
#  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
#GPThtGPD

##### Pmax 
# Pmax
PmaxTs<-aov(Pmax~species*treatment,testGPP)

summary(PmaxTs)

Pmax<-lm(Pmax~species,testGPP)
summary(Pmax)

PmaxT<-lm(Pmax~treatment,testGPP)
Anova(PmaxT)

GPPmaxtr<-ggplot(testGPP, aes(x=species, Pmax, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  ggtitle("Pmax")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
GPPmaxtr


GPPmaxspp<-ggplot(testGPP, aes(x=species, Pmax, fill= species))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(shape= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  ggtitle("Pmax Photosynthesis")+
  stat_compare_means(label.x= "MCAV", label.y=1.5, method = "anova")+
  scale_fill_manual(breaks = c("DLAB", "MCAV", "OFRA", "PAST"), values= c( "#3399FF", "#FF9900", "#339900","#CC0000"))
GPPmaxspp


GPPmaxspp<-ggplot(testGPP, aes(x=species, Pmax))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  ggtitle("Pmax")+
  stat_compare_means(label.x= "MCAV", label.y=1.3, method = "anova")
GPPmaxspp

ggexport(GPPmaxtr, filename = "GPPmaxtr.pdf")
ggexport(GPPmaxspp, filename = "GPPmaxspp.pdf")

### Checkign DLAB... 
#DlabPmax<-subset(testGPP, species=="DLAB")
#DlabPma<-lm(Pmax~treatment,DlabPmax)
#Anova(DlabPma)
#leveneTest(Pmax ~ treatment, data = DlabPmax)
#res.Dpmac <- aov(Pmax ~ treatment, data = DlabPmax)
#aov_residualpma <- residuals(object = res.Dpmac)
## Run Shapiro-Wilk test
#shapiro.test(x = aov_residualpma )
#
#GPThtGPD<-ggplot(DlabPma, aes(x=treatment, Pmax, fill= treatment))+
#  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
#  geom_boxplot(position = position_dodge(1))+
#  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
#  guides(fill= guide_legend(override.aes = list(shape= NA)))+
#  theme_classic()+
#  xlab(NULL)+
#  ylab('log rate')+
#  ggtitle("DLAB Pmax")+
#  stat_compare_means(label.x= "deep", label.y=1.3, method = "anova")+
#  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
#GPThtGPD

#### Add all  GP graphs together E, E, lnc, Th, Topt, Pmax


GPtreatall<-ggpubr::ggarrange(GPEtr,GPEhspp,GPTh_Cspp,GPtopt_Cspp,GPlncspp,GPPmaxspp
 + rremove("x.text"),
                           ncol = 2, nrow = 3)
GPtreatall
ggexport(GPtreatall, filename = "GPtreatall.pdf")

GPsppall<-ggpubr::ggarrange(GPEspp,GPEhtr,GPTh_Ctr,GPtopt_Ctr,GPlnctr,GPPmaxtr
 + rremove("x.text"),
                               ncol = 2, nrow = 3)
GPsppall
ggexport(GPsppall, filename = "GPsppall.pdf")
## no Anovas
## Graphs
GPEsppA<-ggplot(testGPEs, aes(x=species, E, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab('Energy (eV)')+
  ggtitle("E")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
GPEsppA

GPEhtrA<-ggplot(testGPEh, aes(x=species, Eh, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab('Energy (eV)')+
  ggtitle("Eh")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
GPEhtrA

GPTh_CtrA<-ggplot(testGPTh_Ct, aes(x=species, Th_C, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab("Temperature (ºC)")+
  ggtitle("Th")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
GPTh_CtrA

GPtopt_CtrA<-ggplot(testGPt, aes(x=species, topt_C, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab("Temperature (ºC)")+
  ggtitle("Topt")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
GPtopt_CtrA

GPlnctrA<-ggplot(testGPlnct, aes(x=species, lnc, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  ggtitle("lnc")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
GPlnctrA

GPPmaxtrA<-ggplot(testGPP, aes(x=species, Pmax, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  ggtitle("Pmax")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("green4", "green1"))
GPPmaxtrA
## save graphs together for manuscript
finalfigure3<-ggpubr::ggarrange(GPEsppA,GPEhtrA,GPTh_CtrA,GPtopt_CtrA,GPlnctrA,GPPmaxtrA
                            + rremove("x.text"),
                            ncol = 2, nrow = 3)
finalfigure3
ggexport(finalfigure3, filename = "finalfigure3.pdf")
finalfigure3
ggexport(GPEhtrD, filename = "GPEhtrD.pdf")
##### FIRST ALL GP&R coefficient plots

##### 
##### - then do GP & Respiration boxplots don't forget log.E for Respiration
##

#MAKE below ALL for respiration then copy into manuscript figure 3282020
#then ALL the box plot figures will be done... Then to coefindMethodSignatures(
#  for all species interactions for R AND R
  
  
  # R E, Eh,lnc, Th, Topt, Pmax
  # E
  Er<-lm(E~species,testR)
  summary(Er)
  Anova(Er)
  
  ErT<-lm(E~treatment,testR)
  Anova(ErT)
  
  ErTs<-lm(E~species*treatment,testR)
  Anova(ErTs)
  
 # logET<-aov(log.E~treatment,testRE)
 # summary(logET)
  
  REspp<-ggplot(testR, aes(x=species, E, fill= treatment))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab('Energy (eV)')+
    ggtitle("E")+
    scale_fill_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  REspp
  
  
  
  REtr<-ggplot(testR, aes(x=species, E, fill= species))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(shape= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab('Energy (eV)')+
    ggtitle("E")+
    stat_compare_means(label.x= "MCAV", label.y=40, method = "anova")+
    scale_fill_manual(breaks = c("DLAB", "MCAV", "OFRA", "PAST"), values= c( "#3399FF", "#FF9900", "#339900","#CC0000"))
  REtr
  
  ggexport(REspp, filename = "REspp.pdf")
  ggexport(REtr, filename = "REtr.pdf")
  
  
   ## NExt parameters below
  ## 
  ## R Eh
  #DLAB is SIG BETWEEN DEPTHS:
  # Eh
  EhR<-lm(Eh~species,testREh)
  summary(EhR)
  Anova(EhR)
  
  EhTr<-lm(Eh~treatment,testREh)
  Anova(EhTr)
  
  EhTsr<-lm(Eh~species*treatment,testREh)
  Anova(EhTsr)
  
  
  REhtr<-ggplot(testREh, aes(x=species, Eh, fill= treatment))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab('Energy (eV)')+
    ggtitle("Eh")+
    stat_compare_means(label.x= "MCAV", label.y=45, method = "anova")+
    scale_fill_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  REhtr
  
  REhspp<-ggplot(testREhs, aes(x=species, Eh, fill= species))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(shape= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab('Energy (eV)')+
    ggtitle("Eh Respiration")+
    stat_compare_means(label.x= "MCAV", label.y=45, method = "anova")+
    scale_fill_manual(breaks = c("DLAB", "MCAV", "OFRA", "PAST"), values= c( "#3399FF", "#FF9900", "#339900","#CC0000"))
  REhspp
  
  ggexport(REhtr, filename = "REhspp.pdf")
  ggexport(REhspp, filename = "REhtr.pdf")
  
  ## R
  ## lnc
  # lnc
  lncTsr<-aov(lnc~species*treatment,testRlnc)
  Anova(lncTsr)
  
  lnc<-lm(lnc~species,testRlnc)
  Anova(lnc)
  
  lncT<-lm(lnc~treatment,testRlnc)
  Anova(lncT)
  
  
  Rlnctr<-ggplot(testRlnc, aes(x=species, lnc, fill= treatment))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
    ggtitle("lnc")+
    stat_compare_means(label.x= "MCAV", label.y=42, method = "anova")+
    scale_fill_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  Rlnctr
  
  Rlncspp<-ggplot(testRlnc, aes(x=species, lnc, fill= species))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(shape= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
    ggtitle("lnc Respiration")+
    stat_compare_means(label.x= "MCAV", label.y=46, method = "anova")+
    scale_fill_manual(breaks = c("DLAB", "MCAV", "OFRA", "PAST"), values= c( "#3399FF", "#FF9900", "#339900","#CC0000"))
  Rlncspp
  
  ggexport(Rlnctr, filename = "Rlncspp.pdf")
  ggexport(Rlncspp, filename = "Rlnctr.pdf")
  
  ## Th
  # Th_C
  Th_CTs<-aov(log.Th_C~species*treatment,testRTh_C)
  summary(Th_CTs)
  
  Th_C<-lm(log.Th_C~species,testRTh_C)
  summary(Th_C)
  
  Th_CT<-lm(log.Th_C~treatment,testRTh_C)
  Anova(Th_CT)
  
  RTh_Ctr<-ggplot(testRTh_C, aes(x=species, log.Th_C, fill= treatment))+
    #stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab("Temperature (ºC)")+
    ggtitle("log Th")+
    scale_fill_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  RTh_Ctr
  
  RTh_Cspp<-ggplot(testRTh_C, aes(x=species, log.Th_C))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(shape= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab("Temperature (ºC)")+
    ggtitle("log Th")+
    stat_compare_means(label.x= "MCAV", label.y=5, method = "anova")
  RTh_Cspp
  
  ggexport(RTh_Ctr, filename = "RTh_Ctr.pdf")
  ggexport(RTh_Cspp, filename = "RTh_spp.pdf")
  ##### Topt
  # topt_C
  topt_CTs<-aov(topt_C~species*treatment,testRt)
  summary(topt_CTs)
  
  topt_C<-lm(topt_C~species,testRt)
  summary(topt_C)
  
  topt_CT<-lm(topt_C~treatment,testRt)
  Anova(topt_CT)
  
  Rtopt_Ctr<-ggplot(testRt, aes(x=species, topt_C, fill= treatment))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab("Temperature (ºC)")+
    ggtitle("Topt")+
    stat_compare_means(label.x= "MCAV", label.y=40, method = "anova")+
    scale_fill_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  Rtopt_Ctr
  
  Rtopt_Cspp<-ggplot(testRt, aes(x=species, topt_C, fill= species))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(shape= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab("Temperature (ºC)")+
    ggtitle("Topt Respiration")+
    stat_compare_means(label.x= "MCAV", label.y=38, method = "anova")+
    scale_fill_manual(breaks = c("DLAB", "MCAV", "OFRA", "PAST"), values= c( "#3399FF", "#FF9900", "#339900","#CC0000"))
  Rtopt_Cspp
  
  ggexport(Rtopt_Ctr, filename = "Rtopt_Cspp.pdf")
  ggexport(Rtopt_Cspp, filename = "Rtopt_Ctr.pdf")
  ##### Pmax 
  # Pmax
  PmaxTs<-aov(Pmax~species*treatment,testR)
  summary(PmaxTs)
  
  Pmax<-lm(Pmax~species,testR)
  summary(Pmax)
  
  PmaxT<-lm(Pmax~treatment,testR)
  Anova(PmaxT)
  
  RPmaxtr<-ggplot(testR, aes(x=species, Pmax, fill= treatment))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
    ggtitle("Pmax")+
    scale_fill_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  RPmaxtr
  
  RPmaxspp<-ggplot(testR, aes(x=species, Pmax))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(shape= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
    ggtitle("Pmax")+
    stat_compare_means(label.x= "MCAV", label.y=.8, method = "anova")
  RPmaxspp
  
  ggexport(RPmaxtr, filename = "RPmaxtr.pdf")
  ggexport(RPmaxspp, filename = "RPmaxspp.pdf")
  
  
  #### Add all  R graphs together E, E, lnc, Th, Topt, Pmax
  
  
  Rtreatall<-ggpubr::ggarrange(REtr,REhspp,RTh_Cspp,Rtopt_Cspp,Rlncspp,RPmaxspp
                               + rremove("x.text"),
                               ncol = 2, nrow = 3)
  Rtreatall
  ggexport(Rtreatall, filename = "Rtreatall.pdf")
  
  Rsppall<-ggpubr::ggarrange(REspp,REhtr,RTh_Ctr,Rtopt_Ctr,Rlnctr,RPmaxtr
                             + rremove("x.text"),
                             ncol = 2, nrow = 3)
  Rsppall
  ggexport(Rsppall, filename = "Rsppall.pdf")
  ## no Anovas
  ## Graphs
  REsppA<-ggplot(testRE, aes(x=species, log.E, fill= treatment))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab('Energy (eV)')+
    ggtitle("log E")+
    scale_fill_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  REsppA
  
  REhtrA<-ggplot(testREht, aes(x=species, Eh, fill= treatment))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab('Energy (eV)')+
    ggtitle("Eh")+
    scale_fill_manual(breaks = c("deep", "shallow"),values= c("#003366","#00CCFF"))
  REhtrA
  RTh_CtrA<-ggplot(testRTh_Cs, aes(x=species, Th_C, fill= treatment))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab("Temperature (ºC)")+
    ggtitle("Th")+
    scale_fill_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  RTh_CtrA
  Rtopt_CtrA<-ggplot(testRt, aes(x=species, topt_C, fill= treatment))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab("Temperature (ºC)")+
    ggtitle("Topt")+
    scale_fill_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  Rtopt_CtrA
  
  RlnctrA<-ggplot(testRlnc, aes(x=species, lnc, fill= treatment))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
    ggtitle("lnc")+
    scale_fill_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  RlnctrA
  RPmaxtrA<-ggplot(testRP, aes(x=species, Pmax, fill= treatment))+
    stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
    geom_boxplot(position = position_dodge(1))+
    #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
    guides(fill= guide_legend(override.aes = list(shape= NA)))+
    theme_classic()+
    xlab(NULL)+
    ylab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
    ggtitle("Pmax")+
    scale_fill_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  RPmaxtrA
  ## save graphs together for manuscript
  finalfigure4<-ggpubr::ggarrange(REsppA,REhtrA,RTh_CtrA,Rtopt_CtrA,RlnctrA,RPmaxtrA
                              + rremove("x.text"),
                              ncol = 2, nrow = 3)
  
  
  finalfigure4
  ggexport(finalfigure4, filename = "finalfigure4.pdf")
  
  
### Coefficient plot for the interaction spp*depth
### 
  library(sjPlot)
  library(sjmisc)
  library(jtools)
  library(ggplot2)
#  fitGPE <- lm(E ~ species*treatment, data = testGPE)
#  
#summary(fitGPE)
#  anova(fitGPE)
#  GPEplot<-plot_model(fitGPE, type = "int", mdrt.values = "meansd", dodge= 0.5)+
#    scale_color_manual(breaks = c("deep", "shallow"), values= c("green4","green1"))
#  
 
  fitRE <- lm(E ~ species*treatment, data = testR)
  summary(fitRE)
  anova(fitRE)
  REplot<- plot_model(fitRE, type = "int", mdrt.values = "meansd", dodge= 0.5 )+
    scale_color_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  fitRE
   #### 
  
  ###### EH
  fitGPEh <- lm(Eh ~ species*treatment, data = GPEh)
  summary(fitGPEh)
  anova(fitGPEh)
  GPEhplot<-plot_model(fitGPEh, type = "int", mdrt.values = "meansd", dodge= 0.5 )+
    scale_color_manual(breaks = c("deep", "shallow"), values= c("green4","green1"))
 
  fitREh <- lm(Eh ~ species*treatment, data = testREh)
  summary(fitREh)
  anova(fitREh)
  REhplot<-plot_model(fitREh, type = "int", mdrt.values = "meansd",  dodge= 0.5 )+
    scale_color_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  REhplot
  
  ###### Th
  fitGPTh <- lm(Th_C ~ species*treatment, data = testGP)
  
  summary(fitGPTh)
  anova(fitGPTh)
  GPthplot<-plot_model(fitGPTh, type = "int", mdrt.values = "meansd", dodge= 0.5)+
    scale_color_manual(breaks = c("deep", "shallow"), values= c("green4","green1"))
  
  
  fitRTh <- lm(log.Th_C ~ species*treatment, data = testRTh_C)
  summary(fitRTh)
  anova(fitRTh)
  RThplot<- plot_model(fitRTh, type = "int", mdrt.values = "meansd", dodge= 0.5)+
    scale_color_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  

  
  ###### Topt
  fitGPtopt <- lm(topt_C ~ species*treatment, data = testGPt)
  summary(fitGPtopt)
  anova(fitGPtopt)
  GPtoptplot<-plot_model(fitGPtopt, type = "int", mdrt.values = "meansd", dodge= 0.5 )+
    scale_color_manual(breaks = c("deep", "shallow"), values= c("green4","green1"))
  
  fitRtopt <- lm(topt_C ~ species*treatment, data = testRt)
  summary(fitRtopt)
  anova(fitRtopt)
  Rtoptplot<-plot_model(fitRtopt, type = "int", mdrt.values = "meansd",  dodge= 0.5 )+
    scale_color_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  
  
  ###### lnc

  fitGPlnc <- lm(lnc ~ species*treatment, data = testGPlnc)
  
  summary(fitGPlnc)
  anova(fitGPlnc)
  GPlncplot<-plot_model(fitGPlnc, type = "int", mdrt.values = "meansd",  dodge= 0.5)+
    scale_color_manual(breaks = c("deep", "shallow"), values= c("green4","green1"))
  
  
  fitRlnc <- lm(lnc ~ species*treatment, data = testRlnc)
  summary(fitRlnc)
  anova(fitRlnc)
  Rlncplot<- plot_model(fitRlnc, type = "int", mdrt.values = "meansd", dodge= 0.5 )+
    scale_color_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  
  
  
  ###### Pmax
  fitGPPmax <- lm(Pmax ~ species*treatment, data = testGPP)
  summary(fitGPPmax)
  anova(fitGPPmax)
  GPpmaxplot<-plot_model(fitGPPmax, type = "int", mdrt.values = "meansd", dodge= 0.5 )+
    scale_color_manual(breaks = c("deep", "shallow"), values= c("green4","green1"))
  
  fitRPmax <- lm(Pmax ~ species*treatment, data = testR)
  summary(fitRPmax)
  anova(fitRPmax)
  Rpmaxplot<-plot_model(fitRPmax, type = "int", mdrt.values = "meansd", dodge= 0.5 )+
    scale_color_manual(breaks = c("deep", "shallow"), values= c("#003366","#00CCFF"))
  
  
  
  
  
  
  GP_spptreat_coef<-ggpubr::ggarrange(  GPEplot,  GPEhplot,  GPthplot,  GPtoptplot,  GPlncplot,  GPpmaxplot + rremove("x.text"),
                                ncol = 2, nrow = 3)
  
 GP_spptreat_coef
 ggexport(GP_spptreat_coef, filename = "GP_spptreat_coef.pdf")
  
  

  
  
  R_spptreat_coef<-ggpubr::ggarrange(REplot,  REhplot,  RThplot,  Rtoptplot,  Rlncplot,  Rpmaxplot + rremove("x.text"),
                                ncol = 2, nrow = 3)
  
  R_spptreat_coef
  ggexport(R_spptreat_coef, filename = "R_spptreat_coef.pdf")
  
  
  ##### individual Topts across depth for each species
  ##### 
  
Pasttopttreat<-subset(testGPt, species=="PAST")
PTopttr<-lm(topt_C~treatment,Pasttopttreat)
Anova(PTopttr)
To.dep <- lm(topt_C~treatment,Pasttopttreat)
  
#Code for Coefficents without densities treat
coefTod <- as.data.frame(summary(To.dep)$coefficients)
coefTod
names(coefTod)[2] = "se" 
coral <- c(treatmentshallow="")
coefTod$vars = coral
coefTod<-coefTod[-1,]

GpToP<-ggplot(coefTod, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#003300", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#33FF33", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab('Temperature (ºC)')+
  xlab("")+
    ggtitle('Topt PAST')
GpToP
  
#respiration Topt
PasttopttreatR<-subset(testR, species=="PAST")
PTopttrr<-lm(topt_C~treatment,PasttopttreatR)
Anova(PTopttrr)
To.dep <- lm(topt_C~treatment,PasttopttreatR)

#Code for Coefficents without densities treat
coefTod <- as.data.frame(summary(To.dep)$coefficients)
coefTod
names(coefTod)[2] = "se" 
coral <- c(treatmentshallow="")
coefTod$vars = coral
coefTod<-coefTod[-1,]

RGpToP<-ggplot(coefTod, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se),
                lwd=1, colour="#00CCFF", width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour="#003366", width=0) +
  geom_point(size=5, pch=21, fill="yellow") +
  theme_bw()+
  coord_flip()+
  ylab('Temperature (ºC)')+
  xlab("")+
  ggtitle('Topt PAST')
RGpToP
