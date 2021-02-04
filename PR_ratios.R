### P/R ration
library(readr)
raw.data <- read_csv("R_output/Final/final_NP_R_GP_edited.csv")
#raw.data<-read.csv("C:/Github/TPC_BDA_github/R_output/Final/final_NP_R_NP_edited.csv")
#View(raw.data)


#filterng out NP, removes it from the list 
mydata<-raw.data
mydata<-filter(mydata, rate.type !="NP") 

mydatanP<-subset(mydata,rate.type=="GP")
mydataR<-subset(mydata,rate.type=="R")

#Take absolute values of Respiration to log it, then combine dataframes and log both R and GP
mydataR$umol.cm2.hr<-abs(mydataR$umol.cm2.hr)
mydatanP$umol.cm2.hr<-abs(mydatanP$umol.cm2.hr)
Mydata<-rbind(mydatanP,mydataR)

mydata<-Mydata

mydata$log.rate <- log(mydata$umol.cm2.hr + 1) 

mydata%<>%
  mutate(K=mydata$Temp.C + 273.15)

mydata$temp.Cat<-factor(mydata$temp.Cat)

mydataPR<-subset(mydata,rate.type=="GP")%>%
  dplyr::select(temp.Cat,treatment, genotype,fragment.ID, species,log.rate)
mydataRP<-subset(mydata,rate.type=="R")%>%
  dplyr::select(temp.Cat,treatment, genotype,fragment.ID, species,log.rate)


mydataPR$P<-mydataPR$log.rate
mydataPR$log.rate=NULL

mydataRP$R<-mydataRP$log.rate
mydataRP$log.rate=NULL


all<-left_join(mydataPR,mydataRP)

All <- all %>%
  group_by(genotype) %>%
  mutate(PR = P/R)

#write.csv(All, "C:/Github/TPC_BDA_github/AllPR.csv")
#AllD<-subset(All, treatment=="deep")
#AllS<-subset(All, treatment=="shallow")
#######
A1<- aov(PR~species, All)
A2<- aov(PR~treatment, All)
A3<- aov(PR~species*treatment, All)
A4<- aov(PR~species*treatment*temp.Cat, All)
A5<- aov(PR~species*treatment*genotype, All)
A6<- aov(PR~genotype, All)

summary(A1)
#0.0631
summary(A2)
summary(A3)
summary(A4)
summary(A5)
summary(A6)

Anova(A1)
Anova(A2)
Anova(A3)
Anova(A4)
Anova(A5)
Anova(A6)
#######
#######
mPGPsp<-lm(PR~species*treatment,All)
summary(mPGPsp)
mPGPsp2<-aov(PR~species*treatment,All)
summary(mPGPsp2)
anova(mPGPsp)
anova(mPGPsp2)

fiti<-lm(PR~species*treatment,All)
summary(fiti)
Anova(fiti)

fit<-lm(PR~species*treatment*temp.Cat,All)
summary(fit)
Anova(fit)

fitg<-lm(PR~species*treatment*genotype,All)
summary(fitg)
Anova(fitg)


AIC(fit,fiti,fitil)
set_theme(
  base = theme_classic())

PRplot2<-plot_model(fit, type = "pred", dodge=1,show.p, colors = c("magenta", "maroon"), terms = c("treatment", "species"))
PRplot3<-plot_model(fit, type = "pred", dodge=1,show.p, colors = c("maroon", "magenta"), terms = c("temp.Cat","treatment", "species"))
PRplot4<-plot_model(fit, type = "pred", dodge=1,show.p, colors = c("yellow","orange", "green", "blue"), terms = c("temp.Cat","species"))
PRplot5<-plot_model(fit, type = "pred", dodge=1,show.p, colors = c("magenta", "maroon"), terms = c("temp.Cat","treatment"))
PRplot6<-plot_model(fit, type = "pred", dodge=1,show.p, terms = "species")
PRplot7<-plot_model(fit, type = "pred", dodge=1, terms = "treatment")

par(mfrow=c(4,2))
PRplot1<-plot_model(fit, type = "int", dodge=1,show.p, colors = c("maroon", "magenta"), terms = c("treatment", "species"))

ggexport(PRplot1, filename = " PRplot1.pdf")
ggexport(PRplot2, filename = " PRplot2.pdf")
ggexport(PRplot3, filename = " PRplot3.pdf")
ggexport(PRplot4, filename = " PRplot4.pdf")
ggexport(PRplot5, filename = " PRplot5.pdf")
ggexport(PRplot6, filename = " PRplot6.pdf")
ggexport(PRplot7, filename = " PRplot7.pdf")

library(ggeffects)
#ggeffect(fit, terms, ci.lvl = 0.95)
predfit<-ggpredict(fit, terms = c("temp.Cat","species","treatment"),ci.lvl = 0.95)

ggplot(predfit, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)

All$logPR<-log(All$PR +1)

fitil<-lm(logPR~species*treatment+temp.Cat,All)
summary(fitil)
Anova(fitil)


plot_model(fitil, type = "int", dodge=1, show.p, colors = c("magenta", "maroon"), terms = c("treatment", "species"))
plot_model(fitil, type = "eff", dodge=1,show.p,colors = c("magenta", "maroon"), terms = c("temp.Cat","treatment", "species"))
plot_model(fitil, type = "pred", dodge=1,colors = c("magenta", "maroon"), terms = c("temp.Cat","species"))
plot_model(fitil, type = "pred", dodge=1,colors = c("magenta", "maroon"), terms = c("temp.Cat","treatment"))
plot_model(fitil, type = "pred", dodge=1,show.p, terms = "species")
plot_model(fitil, type = "pred", dodge=1,show.p, terms = "treatment")

#
#ggeffect(model, terms, ci.lvl = 0.95,
predfitil<-ggpredict(fitil, terms = ~temp.Cat*species*treatment,ci.lvl = 0.95)

get_model_data(fitil, type = "pred")



fitili<-lm(logPR~species*treatment*temp.Cat,All)
summary(fitili)

plot_model(fitili, type = "int", dodge=1, terms = c("temp.Cat","treatment", "species"))


fit1<-lm(PR~treatment, All)
summary(fit1)
TukeyHSD(aov(PR~treatment, All))

plot_model(fit1, type = "est", dodge=1)
plot_model(fit, type = "int", dodge=1, terms = c("treatment", "species"))

#ttest of slopes across treatments and species

#fitCI<-confint(fit)
#fitCI2<-predict(fit, interval = "prediction")
plot_model(fit, type = "pred", terms = c("temp.Cat", "treatment", "species"))

ggplot(All, aes(x=temp.Cat, y=PR, color=treatment)) +
  geom_point() + 
  geom_smooth(method=lm)+
  scale_color_manual(breaks = c("deep", "shallow"), values= c("maroon","magenta"))+
  facet_wrap(c("species","treatment"))

##################
library(lsmeans)
m.interaction <- lm(PR ~ species*treatment, data = All)
anova(m.interaction)
m.interaction <- glm(PR ~ species*treatment, data = All)
anova(m.interaction)##################

ggplot(AllD, aes(x = temp.Cat, y = PR, color = species) ) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap("species")

ggplot(AllS, aes(x = temp.Cat, y = PR, color = species) ) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap("species")



ggplot(AllD, aes(x = temp.Cat, y = PR, color = species))+
  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x)+
  facet_wrap("species")

########
m.species <- lm(PR ~ species, data = All)
Anova(m.species)
m.species <- aov(PR ~ species, data = All)
summary(m.species)
TukeyHSD(m.species)
m.treat <- lm(PR ~ treatment, data = All)
Anova(m.treat)
m.treat <- aov(PR ~ treatment, data = All)
summary(m.treat)
TukeyHSD(m.treat)

res.aov3 <- aov(PR~ species, data = All)

summary(res.aov3) 
TukeyHSD(res.aov3)

model.tables(res.aov3, type="means", se= TRUE)

sumlncGP = summarySE(testGPlnc,
                     measurevar="lnc",            groupvars=c("treatment","species"))

sumlncGP

######
######
sumpmaxPR = summarySE(All, na.rm=TRUE,measurevar="PR",
                     groupvars=c("treatment", "species", "temp.Cat"))
sumpmaxPR

write.csv(sumpmaxPR,"C:/Github/TPC_BDA_github/PRsums.csv")

sumpmaxPRt = summarySE(All, na.rm=TRUE,measurevar="PR",
                      groupvars=c("treatment", "temp.Cat"))
sumpmaxPRt

write.csv(sumpmaxPRt,"C:/Github/TPC_BDA_github/PRsumst.csv")


AllPR<-ggplot(All, aes(x=temp.Cat, PR, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab('P:R')+
  ggtitle("")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("maroon", "magenta"))
AllPR

ggexport(AllPR, filename = "AllPR.pdf")


AllPRs<-ggplot(All, aes(x=temp.Cat, PR, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab('P:R')+
  ggtitle("")+
  facet_wrap("species")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("maroon", "magenta"))
AllPRs

ggexport(AllPRs, filename = "AllPRs.pdf")
###############################
###############################
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df3 <- data_summary(All, varname="PR",
                    groupnames=c( "treatment","species", "temp.Cat"))
# Convert dose to a factor variable
head(df3)
#method = lm, se = FALSE
p<- ggplot(df3, aes(x=temp.Cat, y=PR, group=treatment, color=treatment)) + 
  geom_smooth(se=TRUE, span= 1.5) +
  geom_point(size=2,position=position_dodge(0.5))+
  #geom_errorbar(aes(ymin=PR-sd, ymax=PR+sd), width=.2,
               # position=position_dodge(0.5))+
  labs(y="P:R", x = "Temperature (ºC)")+
  theme_classic() +
  scale_color_manual(values=c('maroon','magenta'))+
  facet_wrap("species")
p
p2<-p+ geom_hline(yintercept=1, linetype="dashed", color = "red")
p2
ggexport(p, filename = "PR_spp.pdf")
ggexport(p2, filename = "PR_sppwithcut.pdf")


###
df2t <- data_summary(All, varname="PR", 
                    groupnames=c("treatment", "temp.Cat"))
# Convert dose to a factor variable
head(df2t)

pt<- ggplot(df2t, aes(x=temp.Cat, y=PR, group=treatment, color=treatment)) + 
  geom_smooth(se=TRUE, span= 1.5) +
  geom_point(size=2,position=position_dodge(0.5))+
  #geom_errorbar(aes(ymin=PR-sd, ymax=PR+sd), width=.2,
   #             position=position_dodge(0.5))+
  labs(y="P:R", x = "Temperature (ºC)")+
  theme_classic() +
  scale_color_manual(values=c('maroon','magenta'))
pt

pt2<-pt+ geom_hline(yintercept=1, linetype="dashed", color = "red")


ggexport(pt, filename = "PR_treat.pdf")
ggexport(pt2, filename = "PR_treatcut.pdf")


###### Just up to 29C to see slope
###### 
All$temp.Cat<-as.numeric(as.character(All$temp.Cat))
Allslope<-subset(All, temp.Cat<= 29)
df2tslope <- data_summary(All, varname="PR", 
                     groupnames=c("species","treatment", "temp.Cat"))
# Convert dose to a factor variable
head(df2tslope)

pslop<- ggplot(df2tslope, aes(x=temp.Cat, y=PR, group=treatment, color=treatment)) + 
  geom_smooth(se=TRUE, span= 1.5) +
  geom_point(size=2,position=position_dodge(0.5))+
  #geom_errorbar(aes(ymin=PR-sd, ymax=PR+sd), width=.2,
              # position=position_dodge(0.5))+
  labs(y="P:R", x = "Temperature (ºC)")+
  theme_classic() +
  scale_color_manual(values=c('maroon','magenta'))+
  facet_wrap(~species)
pslop
pslop2spp<-pslop+ geom_hline(yintercept=2, linetype="dashed", color = "red")
pslop2spp

## Diifferent color for manuscript
pslop<- ggplot(df2tslope, aes(x=temp.Cat, y=PR, group=treatment, color=treatment)) + 
  geom_smooth(se=TRUE, span= 1.5) +
  geom_point(size=2,position=position_dodge(0.5))+
  #geom_errorbar(aes(ymin=PR-sd, ymax=PR+sd), width=.2,
  # position=position_dodge(0.5))+
  labs(y="P:R", x = "Temperature (ºC)")+
  theme_classic() +
  scale_color_manual(values=c('blue','red'))+
  facet_wrap(~species)
pslop
pslop2spp<-pslop+ geom_hline(yintercept=2, linetype="dashed", color = "black")
pslop2spp
ggexport(pslop2spp, filename = "Fig7.pdf")


slofit<-lm(PR~treatment,Allslope)
summary(slofit)
Anova(slofit)
ggexport(pslop2spp, filename = " pslop2spp2.pdf")


# Finished line plot
### GAMS of data
All$treatment<-factor(All$treatment)
All$genotype<-factor(All$genotype)
All$species<-factor(All$species)
All$fragment.ID<-factor(All$fragment.ID)

m<-nls(y~a*x/(1+x))
m<-nls(PR~species*treatment/(b+treatment), All)
m
library(mgcv)
gamst<-gam(PR~s(temp.Cat), na.rm=TRUE,data=All, method="REML")
gm<-gam(PR~ s(temp.Cat), All)


############################
AllPR<-ggplot(sumpmaxPRt, aes(x=temp.Cat, PR, fill= treatment))+
  stat_boxplot(geom= "errorbar", width= 0.5, position = position_dodge(1))+
  geom_boxplot(position = position_dodge(1))+
  #geom_point(aes(fill= treatment, group= fragment.ID2), position = position_jitterdodge(jitter.width = 0, jitter.height = 0))+
  guides(fill= guide_legend(override.aes = list(shape= NA)))+
  theme_classic()+
  xlab(NULL)+
  ylab('P:R')+
  ggtitle("")+
  scale_fill_manual(breaks = c("deep", "shallow"), values= c("maroon", "magenta"))
AllPR


### subset by species
Dall<-subset(All, species=="DLAB")
DDall<-subset(Dall,treatment=="deep")
SDall<-subset(Dall,treatment=="shallow")
Mall<-subset(All, species=="MCAV")
DMall<-subset(Mall,treatment=="deep")
SMall<-subset(Mall,treatment=="shallow")
Oall<-subset(All, species=="OFRA")
DOall<-subset(Oall,treatment=="deep")
SOall<-subset(Oall,treatment=="shallow")
Pall<-subset(All, species=="PAST")
DPall<-subset(Pall,treatment=="deep")
SPall<-subset(Pall,treatment=="shallow")


Dfit<-lm(PR~treatment,Dall)
summary(Dall)
Anova(Dfit)

#0.3523
##df= 1
###F= 0.8843
ggplot(Dall, aes(x=temp.Cat, y=PR, color=treatment)) +
  geom_point() + 
  geom_smooth(method=lm)+
  scale_color_manual(breaks = c("deep", "shallow"), values= c("maroon","magenta"))+
  facet_wrap("treatment")

Dall$LoCI <- predict(Dfit, newdata = Dall, 
                      interval = "confidence", 
                      level = 0.9)[, 2]
Dall$HiCI <- predict(Dfit, newdata = Dall, 
                      interval = "confidence", 
                      level = 0.9)[, 3]

ggplot(Dall,
       aes(x = temp.Cat,
           y = PR)) +
  # Add a ribbon with the confidence band
  geom_smooth(
    aes(
      # lower and upper bound of the ribbon
      ymin = LoCI, ymax = HiCI,
      # Different colour for men/women
      fill = treatment, colour = treatment
    ),
    stat = "identity") +
  ylab("P:R") +
  xlab("Temperature (ºC)")+
  facet_wrap("treatment")

######
Ddfit<-lm(PR~temp.Cat,DDall)
summary(Ddfit)
Anova(Ddfit)

DDall$LoCI <- predict(Ddfit, newdata = DDall, 
                        interval = "confidence", 
                        level = 0.9)[, 2]
DDall$HiCI <- predict(Ddfit, newdata = DDall, 
                        interval = "confidence", 
                        level = 0.9)[, 3]


ggplot(DDall,
       aes(x = temp.Cat,
           y = PR)) +
  # Add a ribbon with the confidence band
  geom_smooth(
    aes(
      # lower and upper bound of the ribbon
      ymin = LoCI, ymax = HiCI,
      # Different colour for men/women
      fill = treatment, colour = treatment
    ),
    stat = "identity") +
  ylab("P:R") +
  xlab("Temperature (ºC)")

###### 2-way ANOVAS for depth for each species
Dfit<-lm(PR~treatment,Dall)
summary(Dall)
Anova(Dfit)
#0.3523
Mfit<-lm(PR~treatment,Mall)
summary(Mall)
Anova(Mfit)
#0.5543
Ofit<-lm(PR~treatment,Oall)
summary(Oall)
Anova(Ofit)
#0.9979
Pfit<-lm(PR~treatment,Pall)
summary(Pall)
Anova(Pfit)
#0.7298

#####
mydata$fragment.ID2<-mydata$fragment.ID
#
mydata$fragment.ID2 <- as.character(mydata$fragment.ID2)
#
mydata$fragment.ID2 = substr(mydata$fragment.ID2,1,nchar(mydata$fragment.ID2)-2)

mydata$fragment.ID2<-factor(mydata$fragment.ID2)
##########Photosynthesis to respiration (P/R) ratios estimate the
#degree to which production (photosynthesis rate) by the
#zooxanthellae exceeds maintenance (respiration rate)
#requirements of both the zooxanthellae and the coral
#host (Coles and Jokiel 1977).
#In this study, P/R was
#calculated from the ratio of gross production (GP) to
#264
#respiration (R). Essentially, this measured what proportion
#of the respiratory needs of the coral-algal association
#was met by algal production.



