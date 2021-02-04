##Photosynthesis and Respiration code

rm(list=ls())



library(magrittr) # magrittr::
library(ggplot2)
library(plyr)
library(dplyr) #mutate
library(boot)
#-----Save plot in PP-------##
library(esquisse)
library(rvg)
library(mime)
library(Hmisc)
library(checkmate)
library(car)
library(curl)
##Install packages
# load packages
library(nls.multstart)
library(broom) #augment
library(boot)
library(purrr) #map
library(tidyverse)
library(nlstools)
library(proto)
library(nls2)
library(here)
library(lme4)
library(Matrix)
library(lmerTest)
library(tidyr) #unnest

#set wd
# get the file path

setwd("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final")
getwd()




# had to adjust respiration- to reflect Pgross and respiration (saved here on desktop as FINAL_resp_data)
#read in your rawdata
raw_data <- read.csv("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/final_NP_R_GP_edited.csv")
View(raw_data)

#filterng out NP, removes it from the list 
mydata<-raw_data
mydata<-filter(mydata, rate.type !="NP") 

mydataP<-subset(mydata,rate.type=="GP")
mydataR<-subset(mydata,rate.type=="R")

#### Take absolute value of Respiration then log bot
#### 
#### h respiration and Photosynthesis for Sharpesschcoolfield equation
mydataR$umol.cm2.hr<-abs(mydataR$umol.cm2.hr)

Mydata<-rbind(mydataP,mydataR)

mydata<-Mydata

mydata$log.rate <- log(mydata$umol.cm2.hr + 1) 

#mydata$log.rate <- log(mydata$umol.cm2.hr + 1)  #logging and adding 0.3(-2 was smallest in data set) because a log of zero does not exist

#convert your C temperature to K
mydata%<>%
  mutate(K=mydata$Temp.C + 273.15)



View(mydata)

write.csv(mydata,"C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/mydata.csv")



##############################################-----



#Define the Schoolfield equation:
schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
  Tc <- 273.15 + Tc
  k <- 8.62e-5
  boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp)))
  inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
  return(boltzmann.term + inactivation.term)
}



##renameing fragment ID's to represent ALL temperatures
mydata$fragment.ID2<-mydata$fragment.ID
#
mydata$fragment.ID2 <- as.character(mydata$fragment.ID2)
#
mydata$fragment.ID2 = substr(mydata$fragment.ID2,1,nchar(mydata$fragment.ID2)-2)

#step one subset data: first make subsets for all of your treatment groups 
#-----subset data--------

#by species
DLAB.df<-mydata%>%
  dplyr::filter(species=="DLAB")

MCAV.df<-mydata%>%
  filter(species=="MCAV")

PAST.df<-mydata%>%
  filter(species=="PAST")

OFRA.df<-mydata%>%
  filter(species=="OFRA")

#by Rate.Type
#-------DLAB-------
DLAB.DF.GP<-DLAB.df%>%
  filter(rate.type=="GP")

DLAB.DF.R<-DLAB.df%>%
  filter(rate.type=="R")
#-------MCAV-------
MCAV.DF.GP<-MCAV.df%>%
  filter(rate.type=="GP")

MCAV.DF.R<-MCAV.df%>%
  filter(rate.type=="R")
#-------PAST-------
PAST.DF.GP<-PAST.df%>%
  filter(rate.type=="GP")

PAST.DF.R<-PAST.df%>%
  filter(rate.type=="R")
#-------OFRA-------
OFRA.DF.GP<-OFRA.df%>%
  filter(rate.type=="GP")

OFRA.DF.R<-OFRA.df%>%
  filter(rate.type=="R")

#-finally by treatment
#-------DLAB-------
DLAB.DF.GP.deep<-DLAB.DF.GP%>%
  filter(treatment=="deep")

DLAB.DF.GP.shallow<-DLAB.DF.GP%>%
  filter(treatment=="shallow")

DLAB.DF.R.deep<-DLAB.DF.R%>%
  filter(treatment=="deep")

DLAB.DF.R.shallow<-DLAB.DF.R%>%
  filter(treatment=="shallow")
#-------MCAV-------
MCAV.DF.GP.deep<-MCAV.DF.GP%>%
  filter(treatment=="deep")

MCAV.DF.GP.shallow<-MCAV.DF.GP%>%
  filter(treatment=="shallow")

MCAV.DF.R.deep<-MCAV.DF.R%>%
  filter(treatment=="deep")

MCAV.DF.R.shallow<-MCAV.DF.R%>%
  filter(treatment=="shallow")
#-------PAST-------
PAST.DF.GP.deep<-PAST.DF.GP%>%
  filter(treatment=="deep")

PAST.DF.GP.shallow<-PAST.DF.GP%>%
  filter(treatment=="shallow")

PAST.DF.R.deep<-PAST.DF.R%>%
  filter(treatment=="deep")

PAST.DF.R.shallow<-PAST.DF.R%>%
  filter(treatment=="shallow")
#-------OFRA-------
OFRA.DF.GP.deep<-OFRA.DF.GP%>%
  filter(treatment=="deep")

OFRA.DF.GP.shallow<-OFRA.DF.GP%>%
  filter(treatment=="shallow")

OFRA.DF.R.deep<-OFRA.DF.R%>%
  filter(treatment=="deep")

OFRA.DF.R.shallow<-OFRA.DF.R%>%
  filter(treatment=="shallow")

#-------calculate fitted and CI estimates per speciers rate type and treatment-------
#Step one: make empty DF to fill with fitted values and run function
##"""it is very important to empty the dataframes before running the function multiple times"""

All.fittedT<-data.frame()#create an empty df to fill with fitted values
All.CIT<-data.frame()#create an empty df. to fill with confidence intervals

mult.fit.curves<-function(Data){
  fit2 <- nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 27),
                        data = Data,
                        iter = 500,
                        start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 285),
                        start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                        supp_errors = 'Y',
                        na.action = na.omit,
                        lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))
  #print(fit2)
  preds <- augment(fit2)
  Data%<>%
    mutate(fitted=preds$.fitted,
           residuals=preds$.resid)
  
  All.fittedT<<-rbind(All.fittedT,Data)
  
  fit_boots <- Data %>% 
    modelr::bootstrap(n = 10, id = 'boot_num') %>% #change the number for more straps
    group_by(boot_num) %>%
    mutate(fit = map(strap, ~nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 27),
                                           data = data.frame(.),
                                           iter = 100,
                                           start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 285),
                                           start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                                           lower = c(lnc=-10, E=0, Eh=0, Th=0),
                                           supp_errors = 'Y')
    ))
  fit_boots
  
  # get predictions
  preds_boot <- fit_boots %>%
    unnest(fit %>% map(augment)) %>%
    ungroup()
  
  new_preds <- Data %>%
    do(., data.frame(K = seq(min(.$K), max(.$K), length.out = 250), stringsAsFactors = FALSE))
  
  preds <- augment(fit2, newdata = new_preds)
  
  df.genotype<-(Data$genotype[1])
  df.rate.type<-(Data$rate.type[1])
  df.treatment<-(Data$treatment[1])
  df.species<-(Data$species[1])
  
  
  preds <- fit_boots %>%
    unnest(fit %>% map(augment, newdata = new_preds)) %>%
    # group by each value of K and get quantiles
    group_by(., K) %>%
    summarise(lwr_CI = quantile(.fitted, 0.025),
              upr_CI = quantile(.fitted, 0.975)) %>%
    ungroup() %>%
    merge(., preds, by = 'K')%>%
    mutate(genotype=factor(df.genotype),
           treatment=factor(df.treatment),
           rate.type=factor(df.rate.type),
           species=factor(df.species))
  All.CIT<<-rbind(All.CIT, preds)
  
}
#"""it is very important to empty the dataframes before running the function multiple times"""


#-----put all the subsetted data through the function--------
#First run all of these and get All.fitted. and All. CI
#
#DLAB
mult.fit.curves(DLAB.DF.GP.deep)
mult.fit.curves(DLAB.DF.GP.shallow)
mult.fit.curves(DLAB.DF.R.deep)
mult.fit.curves(DLAB.DF.R.shallow)
#MCAV
mult.fit.curves(MCAV.DF.GP.deep)
mult.fit.curves(MCAV.DF.GP.shallow)
mult.fit.curves(MCAV.DF.R.deep)
mult.fit.curves(MCAV.DF.R.shallow)
#PAST
mult.fit.curves(PAST.DF.GP.deep)
mult.fit.curves(PAST.DF.GP.shallow)
mult.fit.curves(PAST.DF.R.deep)
mult.fit.curves(PAST.DF.R.shallow)
#OFRA
mult.fit.curves(OFRA.DF.GP.deep)
mult.fit.curves(OFRA.DF.GP.shallow)
mult.fit.curves(OFRA.DF.R.deep)
mult.fit.curves(OFRA.DF.R.shallow)



#### Go back to Acropora data and run multcurves again to get and SAVE ALL.CIT AGAIN!! (9/13/2019)
#### 
write.csv(All.CIT,"C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/All.CI.csv")
write.csv(All.fittedT,"C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/All.fitted.csv")

#----Graph the outputs--------
#graph with your wobbly predictions

# grouped by treatment and rate.type- facetwrapped bu individual.ID
TPCall<-ggplot() +
  geom_ribbon(data=subset(All.CIT, lwr_CI>0),
              aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  theme(legend.position = "none")+
  geom_point(data=All.fittedT, aes(x=(K - 273.15), y=log.rate, shape=treatment)) +
  theme(legend.position = "none")+
  geom_line(data=All.fittedT, aes(x=(K - 273.15), y=fitted, colour=species)) +
  theme(legend.position = "none")+
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~species)
esquisse:::ggplot_to_ppt("TPCall")


#########################################
#########################################

#----Graph the outputs--------
#graph with your wobbly predictions by rate.type and genotype
All.fittedr<-All.fittedT%<>%
  mutate(group=factor(paste(rate.type,genotype)))

All.CIr<-All.CIT%<>%
  mutate(group=factor(paste(rate.type,genotype, treatment)))

#graph with your wobbly predictions by rate.type and treatment

All.fittedg<-All.fittedT%<>%
  mutate(group=factor(paste(treatment,rate.type)))

All.CIg<-All.CIT%<>%
  mutate(group=factor(paste(treatment,rate.type)))
################# bu individual. ID
All.fittedID<-All.fittedT%<>%
  mutate(group=factor(paste(treatment,rate.type, individual.ID)))
#----Graph the outputs--------
#graph with your wobbly predictions
All.fittedi<-All.fittedT%<>%
  mutate(group=factor(paste(rate.type,individual.ID)))

All.CIi<-All.CIT%<>%
  mutate(group=factor(paste(rate.type,genotype)))
####
#All.fittedT- CIT
#All.fittedr- CIr
#All.fittedg- CIg
#All.fittedID- CIT
#All.fittedi- All.CIi

#treatment
TPCallT<-ggplot() +
  geom_ribbon(data=subset(All.CIT, lwr_CI>0),
              aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedT, aes(x=(K - 273.15), y=log.rate, shape=genotype)) +
  geom_line(data=All.fittedT, aes(x=(K - 273.15), y=fitted, colour=group)) +
  theme(legend.position = "none")+
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~treatment)
esquisse:::ggplot_to_ppt("TPCallT")


#treatment
TPCallr<-ggplot() +
  geom_ribbon(data=subset(All.CIr, lwr_CI>0),
              aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedr, aes(x=(K - 273.15), y=log.rate, shape=rate.type)) +
  geom_line(data=All.fittedr, aes(x=(K - 273.15), y=fitted, colour=group)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~treatment)
esquisse:::ggplot_to_ppt("TPCallr")
###


#KEEP####
TPCalli<-ggplot() +
  geom_ribbon(data=subset(All.CIg, lwr_CI>0),
              aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedg, aes(x=(K - 273.15), y=log.rate, shape=rate.type)) +
  geom_line(data=All.fittedg, aes(x=(K - 273.15), y=fitted, colour=group)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~individual.ID)
esquisse:::ggplot_to_ppt("TPCalli")


#graph using a smoothing function
STPCall<-ggplot() +
  geom_ribbon(data=subset(All.CIT, lwr_CI>0),aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedT, aes(x=(K - 273.15), y=log.rate, shape=genotype)) +
  geom_smooth(data=All.fittedT, aes(x=(K - 273.15), y=fitted, colour=group, se=F)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~treatment)

esquisse:::ggplot_to_ppt("STPCall")

#graph using a smoothing function
STPCallg<-ggplot() +
  geom_ribbon(data=subset(All.CIg, lwr_CI>0),aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedg, aes(x=(K - 273.15), y=log.rate, shape=treatment)) +
  geom_smooth(data=All.fittedg, aes(x=(K - 273.15), y=fitted, colour=group,se=F)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~genotype)

esquisse:::ggplot_to_ppt("STPCallg")


#graph using a smoothing function and fancy up this graphical representation of biological processes

TPCTr<-ggplot() +
  geom_ribbon(data=subset(All.CIr, lwr_CI>0),aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedr, aes(x=(K - 273.15), y=log.rate, shape=treatment))+ 
  geom_smooth(data=All.fittedr, aes(x=(K - 273.15), y=fitted, colour=group,se=F)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  facet_wrap(~genotype)+
  theme_classic()+
  labs(title="Thermal Performance Curves Across treatments")
esquisse:::ggplot_to_ppt("TPCTr")

TPCr<-ggplot() +
  geom_ribbon(data=subset(All.CIr, lwr_CI>0),aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedr, aes(x=(K - 273.15), y=log.rate, shape=genotype))+ 
  geom_smooth(data=All.fittedr, aes(x=(K - 273.15), y=fitted, colour=group,se=F)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  facet_wrap(~treatment)+
  theme_classic()+
  labs(title="Thermal Performance Curves Across treatments")
esquisse:::ggplot_to_ppt("TPCr")


########################## good by genotype
STPCallg<-ggplot() +
  geom_ribbon(data=subset(All.CIg, lwr_CI>0),aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedg, aes(x=(K - 273.15), y=log.rate, shape=treatment))+ 
  geom_smooth(data=All.fittedg, aes(x=(K - 273.15), y=fitted, colour=group,se=F)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  facet_wrap(~genotype)+
  theme_classic()+
  labs(title="Thermal Performance Curves Across treatments")
esquisse:::ggplot_to_ppt("STPCallg")



###------Blank the All.parameters to get subsetted data to genotype for Topt calculation and ANOVA comparison-----

#----run this second to get parameters for each genotype
All.parameters<-data.frame()#create an empty df to fill with calculated parameters


mult.fit.curves2<-function(Data){
  fit2 <- nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 27),
                        data = Data,
                        iter = 500,
                        start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 285),
                        start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                        supp_errors = 'Y',
                        na.action = na.omit,
                        lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))
  #print(fit2)
  params <- tidy(fit2)
  
  params%<>%
    mutate(genotype=Data$genotype[1],
           treatment=Data$treatment[1],
           rate.type=Data$rate.type[1],
           temp.Cat= Data$temp.Cat[1],
           individual.ID=Data$individual.ID[1])
  print(params)
  All.parameters<<-rbind(params, All.parameters)
}

#Parameters for Topt and Pmax
#(subset
mult.fit.curves2(subset(K2.df.GP.T1, individual.ID== "K201a"))
mult.fit.curves2(subset(M5.df.GP.T1, individual.ID== "M501a"))
mult.fit.curves2(subset(M6.df.GP.T1, individual.ID== "M601a"))
mult.fit.curves2(subset(K2.df.GP.T2, individual.ID== "K202a"))
mult.fit.curves2(subset(M5.df.GP.T2, individual.ID== "M502a"))
mult.fit.curves2(subset(M6.df.GP.T2, individual.ID== "M602a"))
mult.fit.curves2(subset(K2.df.GP.C, individual.ID== "K2Ca"))
mult.fit.curves2(subset(M5.df.GP.C, individual.ID== "M5Ca"))
mult.fit.curves2(subset(M6.df.GP.C, individual.ID== "M6Ca"))
mult.fit.curves2(subset(K2.df.GP.T3, individual.ID== "K203a"))
mult.fit.curves2(subset(M5.df.GP.T3, individual.ID== "M503a"))
mult.fit.curves2(subset(M6.df.GP.T3, individual.ID== "M603a"))
mult.fit.curves2(subset(K2.df.R.T1, individual.ID==  "K201a"))
mult.fit.curves2(subset(M5.df.R.T1, individual.ID==  "M501a"))
mult.fit.curves2(subset(M6.df.R.T1, individual.ID==  "M601a"))
mult.fit.curves2(subset(K2.df.R.T2, individual.ID==  "K202a"))
mult.fit.curves2(subset(M5.df.R.T2, individual.ID==  "M502a"))
mult.fit.curves2(subset(M6.df.R.T2, individual.ID==  "M602a"))
mult.fit.curves2(subset(K2.df.R.C, individual.ID==  "K20Ca"))
mult.fit.curves2(subset(M5.df.R.C, individual.ID==  "M50Ca"))
mult.fit.curves2(subset(M6.df.R.C, individual.ID==  "M60Ca"))
mult.fit.curves2(subset(K2.df.R.T3, individual.ID==  "K203a"))
mult.fit.curves2(subset(M5.df.R.T3, individual.ID==  "M503a"))
mult.fit.curves2(subset(M6.df.R.T3, individual.ID==  "M603a"))
mult.fit.curves2(subset(K2.df.GP.T1, individual.ID== "K201b"))
mult.fit.curves2(subset(M5.df.GP.T1, individual.ID== "M501b"))
mult.fit.curves2(subset(M6.df.GP.T1, individual.ID== "M601b"))
mult.fit.curves2(subset(K2.df.GP.T2, individual.ID== "K202b"))
mult.fit.curves2(subset(M5.df.GP.T2, individual.ID== "M502b"))
mult.fit.curves2(subset(M6.df.GP.T2, individual.ID== "M602b"))
mult.fit.curves2(subset(K2.df.GP.C, individual.ID== "K20Cb"))
mult.fit.curves2(subset(M5.df.GP.C, individual.ID== "M50Cb"))
mult.fit.curves2(subset(M6.df.GP.C, individual.ID== "M60Cb"))
mult.fit.curves2(subset(K2.df.GP.T3, individual.ID== "K203b"))
mult.fit.curves2(subset(M5.df.GP.T3, individual.ID== "M503b"))
mult.fit.curves2(subset(M6.df.GP.T3, individual.ID== "M603b"))
mult.fit.curves2(subset(K2.df.R.T1, individual.ID== "K201b"))
mult.fit.curves2(subset(M5.df.R.T1, individual.ID== "M501b"))
mult.fit.curves2(subset(M6.df.R.T1, individual.ID== "M601b"))
mult.fit.curves2(subset(K2.df.R.T2, individual.ID== "K202b"))
mult.fit.curves2(subset(M5.df.R.T2, individual.ID== "M502b"))
mult.fit.curves2(subset(M6.df.R.T2, individual.ID== "M602b"))
mult.fit.curves2(subset(K2.df.R.C, individual.ID== "K20Cb"))
mult.fit.curves2(subset(M5.df.R.C, individual.ID== "M50Cb"))
mult.fit.curves2(subset(M6.df.R.C, individual.ID== "M60Cb"))
mult.fit.curves2(subset(K2.df.R.T3, individual.ID==  "K203b"))
mult.fit.curves2(subset(M5.df.R.T3, individual.ID==  "M503b"))
mult.fit.curves2(subset(M6.df.R.T3, individual.ID==  "M603b"))
mult.fit.curves2(subset(K2.df.GP.T1, individual.ID== "K201c"))
mult.fit.curves2(subset(M5.df.GP.T1, individual.ID== "M501c"))
mult.fit.curves2(subset(M6.df.GP.T1, individual.ID== "M601c"))
mult.fit.curves2(subset(K2.df.GP.T2, individual.ID== "K202c"))
mult.fit.curves2(subset(M5.df.GP.T2, individual.ID== "M502c"))
mult.fit.curves2(subset(M6.df.GP.T2, individual.ID== "M602c"))
mult.fit.curves2(subset(K2.df.GP.C, individual.ID== "K20Cc"))
mult.fit.curves2(subset(M5.df.GP.C, individual.ID== "M50Cc"))
mult.fit.curves2(subset(M6.df.GP.C, individual.ID== "M60Cc"))
mult.fit.curves2(subset(K2.df.GP.T3, individual.ID== "K203c"))
mult.fit.curves2(subset(M5.df.GP.T3, individual.ID== "M503c"))
mult.fit.curves2(subset(M6.df.GP.T3, individual.ID== "M603c"))
mult.fit.curves2(subset(K2.df.R.T1, individual.ID==  "K201c"))
mult.fit.curves2(subset(M5.df.R.T1, individual.ID==  "M501c"))
mult.fit.curves2(subset(M6.df.R.T1, individual.ID==  "M601c"))
mult.fit.curves2(subset(K2.df.R.T2, individual.ID==  "K202c"))
mult.fit.curves2(subset(M5.df.R.T2, individual.ID==  "M502c"))
mult.fit.curves2(subset(M6.df.R.T2, individual.ID==  "M602c"))
mult.fit.curves2(subset(K2.df.R.C, individual.ID== "K20Cc"))
mult.fit.curves2(subset(M5.df.R.C, individual.ID== "M50Cc"))
mult.fit.curves2(subset(M6.df.R.C, individual.ID== "M60Cc"))
mult.fit.curves2(subset(K2.df.R.T3, individual.ID==  "K203c"))
mult.fit.curves2(subset(M5.df.R.T3, individual.ID==  "M503c"))
mult.fit.curves2(subset(M6.df.R.T3, individual.ID==  "M603c"))
mult.fit.curves2(subset(K2.df.GP.T1, individual.ID== "K201d"))
mult.fit.curves2(subset(M5.df.GP.T1, individual.ID== "M501d"))
mult.fit.curves2(subset(M6.df.GP.T1, individual.ID== "M601d"))
mult.fit.curves2(subset(K2.df.GP.T2, individual.ID== "K202d"))
mult.fit.curves2(subset(M5.df.GP.T2, individual.ID== "M502d"))
mult.fit.curves2(subset(M6.df.GP.T2, individual.ID== "M602d"))
mult.fit.curves2(subset(K2.df.GP.C, individual.ID== "K20Cd"))
mult.fit.curves2(subset(M5.df.GP.C, individual.ID== "M50Cd"))
mult.fit.curves2(subset(M6.df.GP.C, individual.ID== "M60Cd"))
mult.fit.curves2(subset(K2.df.GP.T3, individual.ID== "K203d"))
mult.fit.curves2(subset(M5.df.GP.T3, individual.ID== "M503d"))
mult.fit.curves2(subset(M6.df.GP.T3, individual.ID== "M603d"))
mult.fit.curves2(subset(K2.df.R.T1, individual.ID==  "K201d"))
mult.fit.curves2(subset(M5.df.R.T1, individual.ID==  "M501d"))
mult.fit.curves2(subset(M6.df.R.T1, individual.ID==  "M601d"))
mult.fit.curves2(subset(K2.df.R.T2, individual.ID==  "K202d"))
mult.fit.curves2(subset(M5.df.R.T2, individual.ID==  "M502d"))
mult.fit.curves2(subset(M6.df.R.T2, individual.ID==  "M602d"))
mult.fit.curves2(subset(K2.df.R.C, individual.ID== "K20Cd"))
mult.fit.curves2(subset(M5.df.R.C, individual.ID== "M50Cd"))
mult.fit.curves2(subset(M6.df.R.C, individual.ID== "M60Cd"))
mult.fit.curves2(subset(K2.df.R.T3, individual.ID==  "K203d"))
mult.fit.curves2(subset(M5.df.R.T3, individual.ID==  "M503d"))
mult.fit.curves2(subset(M6.df.R.T3, individual.ID==  "M603d"))


write.csv(All.parameters, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/TPC_params/All.parameters.csv")


#mult.fit.curves2(K2.df.GP.T1)
#mult.fit.curves2(M5.df.GP.T1)
#mult.fit.curves2(M6.df.GP.T1)
#mult.fit.curves2(K2.df.GP.T2)
#mult.fit.curves2(M5.df.GP.T2)
#mult.fit.curves2(M6.df.GP.T2)
#mult.fit.curves2(K2.df.GP.C)
#mult.fit.curves2(M5.df.GP.C)
#mult.fit.curves2(M6.df.GP.C)
#mult.fit.curves2(K2.df.GP.T3)
#mult.fit.curves2(M5.df.GP.T3)
#mult.fit.curves2(M6.df.GP.T3)
#
#mult.fit.curves2(K2.df.R.T1)
#mult.fit.curves2(M5.df.R.T1)
#mult.fit.curves2(M6.df.R.T1)
#mult.fit.curves2(K2.df.R.T2)
#mult.fit.curves2(M5.df.R.T2)
#mult.fit.curves2(M6.df.R.T2)
#mult.fit.curves2(K2.df.R.C)
#mult.fit.curves2(M5.df.R.C)
#mult.fit.curves2(M6.df.R.C)
#mult.fit.curves2(K2.df.R.T3)
#mult.fit.curves2(M5.df.R.T3)
#mult.fit.curves2(M6.df.R.T3)




# plot mean of parameters  across all genotypes (not very attractive)
MeanParams<-left_join(All.CIT, All.parameters)
View(MeanParams)

write.csv(MeanParams, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/TPC_params/Mean.Params.csv")

para<-ggplot(MeanParams, aes(col = rate.type)) +
  geom_point(aes(genotype, estimate)) +
  facet_wrap(~ term, scale = 'free_x', ncol = 4) +
  geom_linerange(aes(genotype, ymin = lwr_CI, ymax = upr_CI)) +
  coord_flip() +
  scale_color_manual(values = c('green4', 'black')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  theme(legend.position = 'top') +
  xlab('curve') +
  ylab('parameter estimate')


esquisse:::ggplot_to_ppt("para")


###################### Obtaining Topt


#48 GROUPS
#192 ESTIAMTES (4 x 48)
#96 estimates/4 treatments= 24/2 (light and Dark)= 12/temperature /3 genotypes= 4 (treatments)
#------T opt from parameters-------
get_topt<-function(E, Th, Eh){
  return((Eh*Th)/(Eh+(8.62e-05*Th*(log((Eh/E))-1))))
}

topt_data<-All.parameters

topt_data<-All.parameters%>%
  select(individual.ID, genotype, rate.type, treatment, estimate, term)%>%
  group_by(individual.ID, genotype, rate.type, treatment)%>%
  spread(term,estimate)

write.csv(topt_data, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/TPC_params/Topt_data.csv")


topt_data$topt<-get_topt(E=topt_data$E, Th=topt_data$Th, Eh=topt_data$Eh)

topt_data$topt_C<-topt_data$topt-273.15

write.csv(topt_data, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/TPC_params/Top_inC.csv")

#-------linear Models does treatment or genotype affect Topt----
library(car)
library(lme4)
library(stats)
library(Rmisc)
library(Hmisc)

topt_nonegGP<-subset(topt_data, rate.type=="GP")
topt_nonegGP<-topt_nonegGP%>%
  filter(topt_C >0)
#M602a, M602c, M60Cb

topt_nonegR<-subset(topt_data, rate.type=="R")
topt_nonegR<-topt_nonegR%>%
  filter(topt_C>0)
#K201b, K20Cb




topt_noneg<-rbind(topt_nonegGP,topt_nonegR)
View(topt_noneg)

write.csv(topt_noneg, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/TPC_params/Topt_noneg.csv")


# plot distribution of Topt ##### ------ Parameter graphs--- 
Toptmeans<-ggplot(subset(topt_noneg, rate.type=="GP"), aes(topt_C)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~rate.type)+
  xlab('Optimum Temperature (ºC)') +
  ggtitle('Distribution of optimum temperatures')

esquisse:::ggplot_to_ppt("Toptmeans")

Toptmeans<-ggplot(subset(topt_noneg, rate.type=="R"), aes(topt_C)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~rate.type)+
  xlab('Optimum Temperature (ºC)') +
  ggtitle('Distribution of optimum temperatures')

esquisse:::ggplot_to_ppt("Toptmeans")


# plot distribution of estimated parameters E, Eh, lnc, Th
p1 <- ggplot(topt_noneg, aes(E)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ rate.type, scales = 'free_x')
p1
esquisse:::ggplot_to_ppt("p1")
# plot distribution of estimated parameters E, Eh, lnc, Th
p2 <- ggplot(topt_noneg, aes(Eh)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ rate.type, scales = 'free_x')
p2
esquisse:::ggplot_to_ppt("p2")

# plot distribution of estimated parameters E, Eh, lnc, Th
p3 <- ggplot(topt_noneg, aes(lnc)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ rate.type)
p3
esquisse:::ggplot_to_ppt("p3")

# plot distribution of estimated parameters E, Eh, lnc, Th
p4 <- ggplot(subset(topt_noneg, rate.type=="GP"), aes(Th)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ rate.type)
p4
esquisse:::ggplot_to_ppt("p4")

p4r <- ggplot(subset(topt_noneg, rate.type=="R"), aes(Th)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ rate.type)
esquisse:::ggplot_to_ppt("p4r")
# plot distribution of estimated parameters E, Eh, lnc, Th by genotype
# 
p6T <- ggplot(subset(topt_noneg, rate.type=="GP"), aes(topt_C)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ genotype, scales = 'free_x')
esquisse:::ggplot_to_ppt("p6T")

p6r <- ggplot(subset(topt_noneg, rate.type=="R"), aes(topt_C)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ genotype, scales = 'free_x')

esquisse:::ggplot_to_ppt("p6r")

p6Tt <- ggplot(subset(topt_noneg, rate.type=="GP"), aes(topt_C)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ treatment, scales = 'free_x')
esquisse:::ggplot_to_ppt("p6Tt")

p6rt <- ggplot(subset(topt_noneg, rate.type=="R"), aes(topt_C)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ treatment, scales = 'free_x')
esquisse:::ggplot_to_ppt("p6rt")

p6 <- ggplot(topt_noneg, aes(E)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ genotype, scales = 'free_x')
esquisse:::ggplot_to_ppt("p6")
# plot distribution of estimated parameters E, Eh, lnc, Th
p7 <- ggplot(topt_noneg, aes(Eh)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ genotype, scales = 'free_x')
esquisse:::ggplot_to_ppt("p7")

# plot distribution of estimated parameters E, Eh, lnc, Th
p8<- ggplot(topt_noneg, aes(lnc)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ genotype)
esquisse:::ggplot_to_ppt("p8")

# plot distribution of estimated parameters E, Eh, lnc, Th 
p5G <- ggplot(subset(topt_noneg, rate.type=="GP"), aes(Th)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ genotype)
esquisse:::ggplot_to_ppt("p5G")

p5r <- ggplot(subset(topt_noneg, rate.type=="R"), aes(Th)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ genotype)
esquisse:::ggplot_to_ppt("p5r")



summary(topt_noneg)
#individual.ID genotype rate.type treatment       E                   Eh               lnc        
#K201a  : 2    K2:29    GP:41     C :19     Min.   : 0.004709   Min.   :  1.588   Min.   :0.3137  
#K201b  : 2    M5:31    NP: 0     T1:24     1st Qu.: 0.127661   1st Qu.:  4.890   1st Qu.:0.3620  
#K201c  : 2    M6:27    R :46     T2:20     Median : 0.215111   Median :  8.447   Median :0.5965  
#K201d  : 2                       T3:24     Mean   : 0.747792   Mean   : 26.952   Mean   :0.6816  
#K202a  : 2                                 3rd Qu.: 0.450472   3rd Qu.: 17.128   3rd Qu.:0.6870  
#K202c  : 2                                 Max.   :15.551237   Max.   :261.107   Max.   :3.7530  
#(Other):75                                                                                       
#Th             topt           topt_C     
#Min.   :297.1   Min.   :298.6   Min.   :25.49  
#1st Qu.:306.2   1st Qu.:304.4   1st Qu.:31.21  
#Median :306.5   Median :304.8   Median :31.64  
#Mean   :306.1   Mean   :304.5   Mean   :31.33  
#3rd Qu.:307.6   3rd Qu.:305.2   3rd Qu.:32.09  
#Max.   :310.2   Max.   :307.2   Max.   :34.04  
#
### Table of Averages per genotype

topt_nonegAVG<-topt_noneg

topt_nonegAVG$avgE<-topt_nonegAVG$E

data.summarytoptE<-topt_noneg %>%
  group_by(genotype, treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(meanE=mean(E), seE=sd(E)/sqrt(n()))#%>%
data.summarytoptEh<-topt_noneg %>%
  group_by(genotype, treatment, rate.type) %>% 
  dplyr::summarise(meanEh=mean(Eh), seEh=sd(Eh)/sqrt(n()))
data.summarytoptlnc<-topt_noneg %>%
  group_by(genotype, treatment, rate.type) %>% 
  dplyr::summarise(meanlnc=mean(lnc), selnc=sd(lnc)/sqrt(n()))
data.summarytoptTh<-topt_noneg %>%
  group_by(genotype, treatment, rate.type) %>% 
  dplyr::summarise(meanTh=mean(Th), seTh=sd(Th)/sqrt(n()))
data.summarytoptTopt<-topt_noneg %>%
  group_by(genotype, treatment, rate.type) %>% 
  dplyr::summarise(meanTopt=mean(topt_C), seTopt=sd(topt_C)/sqrt(n()))

Allsum<-left_join(data.summarytoptTopt, data.summarytoptE)
Allsum2<-left_join(Allsum, data.summarytoptTh)
Allsum3<-left_join(Allsum2,data.summarytoptEh)
Allsumfinal<-left_join(Allsum3, data.summarytoptlnc)

write.csv(Allsumfinal, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/TPC_params/Allparametersums.csv")

all<-Allsumfinal


#visualize groups GP
#
GPraw<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=treatment, topt_C, colour=treatment))+
  geom_boxplot()+
  geom_point()+
  ylab("Topt (ºC)")

esquisse:::ggplot_to_ppt("GPraw")

GPoptS<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=genotype, topt_C, colour=genotype))+
  geom_boxplot()+
  geom_point()+
  ylab("Topt (ºC)")
esquisse::ggplot_to_ppt("GPoptS")


##### ugly graphs

GPtb<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=genotype, topt_C, colour=treatment, shape=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Topt ºC)")+
  labs(title="Mean Topt across Genotypes (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("GPtb")

GPtbS<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=treatment, topt_C, colour=genotype, shape=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Topt ºC)")+
  labs(title="Mean Topt across Treatments (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("GPtbS")




#linear model
gp.mod<-Anova(lm(topt_C~treatment+genotype, data=subset(topt_noneg,rate.type=="GP")))
summary(gp.mod)
gp.mod


# No significant differnece in Topt across treatment or genotype
m1N<-aov(topt_C~treatment+genotype, data=subset(topt_noneg,rate.type=="GP"))
TukeyHSD(m1N)
summary(m1N)

### No difference in Topt
### 
##visualize groups R

RTopt<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=treatment, topt_C, colour=treatment))+
  geom_boxplot()+
  geom_point()+
  ylab("Topt (ºC)")+
  title("Respiration")

esquisse:::ggplot_to_ppt("RTopt")

RtoptS<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=genotype, topt_C, colour=genotype))+
  geom_boxplot()+
  geom_point()+
  ylab("Topt (ºC)")
esquisse::ggplot_to_ppt("RtoptS")

RtbT<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=treatment, topt_C, colour=genotype, shape=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Topt (ºC)")+
  labs(title="Mean Topt across Treatment (R)")+
  coord_flip()

esquisse:::ggplot_to_ppt("RtbT")

RtbS<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=genotype, topt_C, colour=treatment, shape=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Topt (ºC)")+
  labs(title="Mean Topt across genotype (R)")+
  coord_flip()

esquisse:::ggplot_to_ppt("RtbS")




#linear model
R.mod<-Anova(lm(topt_C~treatment+genotype, data=subset(topt_noneg,rate.type=="R")))
summary(R.mod)
R.mod

m2N<-aov(topt_C~treatment+genotype, data=subset(topt_noneg,rate.type=="R")) #treatment significant (PAST) 0.0108
TukeyHSD(m2N)
#Tukey multiple comparisons of means
#95% family-wise confidence level
#
#Fit: aov(formula = topt_C ~ treatment + genotype, data = subset(topt_noneg, rate.type == "R"))
#
#$treatment
#diff       lwr       upr     p adj
#T1-C  -0.9641785 -2.884106 0.9557489 0.5397950
#T2-C  -0.4676784 -2.428899 1.4935425 0.9186728
#T3-C  -0.8143725 -2.734300 1.1055549 0.6691962
#T2-T1  0.4965002 -1.423427 2.4164276 0.8991018
#T3-T1  0.1498060 -1.727920 2.0275321 0.9964874
#T3-T2 -0.3466941 -2.266622 1.5732333 0.9622068
#
#$genotype
#diff       lwr         upr     p adj
#M5-K2 -1.5938785 -3.094899 -0.09285834 0.0352862
#M6-K2 -2.1755992 -3.700637 -0.65056120 0.0035146
#M6-M5 -0.5817207 -2.082741  0.91929941 0.6164893
######################
topt_nonegRep<-read.csv("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/Topt_estimates_noneg_repl.csv")

K2modT<-aov(topt_C~treatment, data=subset(topt_noneg,genotype=="K2")) # treatment significant across depth for PAST 0.0234
summary(K2modT)
TukeyHSD(K2modT)

K2mod<-aov(topt_C~replicas, data=subset(topt_nonegRep,genotype=="K2")) # treatment significant across depth for PAST 0.0234
summary(K2mod)
TukeyHSD(K2mod)

M5modT<-aov(topt_C~treatment, data=subset(topt_noneg,genotype=="M5")) # treatment significant across depth for PAST 0.0234
summary(M5modT)
TukeyHSD(M5modT)

M5mod<-aov(topt_C~replicas, data=subset(topt_nonegRep,genotype=="M5")) # treatment significant across depth for PAST 0.0234
summary(M5mod)
TukeyHSD(M5mod)

M6modT<-aov(topt_C~treatment, data=subset(topt_noneg,genotype=="M6")) # treatment significant across depth for PAST 0.0234
summary(M6modT)
TukeyHSD(M6modT)

M6mod<-aov(topt_C~replicas, data=subset(topt_nonegRep,genotype=="M6")) # treatment significant across depth for PAST 0.0234
summary(M6mod)
TukeyHSD(M6mod)

#-------linear models does treatment or genotype lnc----
#visualize groups GPs
lnc<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=treatment, lnc, colour=treatment))+
  geom_boxplot()+
  geom_point()+
  ylab("Topt (ºC)")

esquisse:::ggplot_to_ppt("lnc")

lnnc2<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=genotype, lnc, colour=genotype))+
  geom_boxplot()+
  geom_point()+
  ylab("Topt (ºC)")
esquisse::ggplot_to_ppt("lnnc2")

GPtLT<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=treatment, lnc, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("lnc")+
  labs(title="Mean lnc across Treatment (GP)")+
  coord_flip()

esquisse:::ggplot_to_ppt("GPtLT")

GPtbLT<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=genotype, lnc, colour=treatment, shape=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("lnc")+
  labs(title="Mean lnc across Genotype (GP)")+
  coord_flip()


esquisse:::ggplot_to_ppt("GPtbLT")


#linear model
IncGP.mod<-Anova(lm(lnc~treatment+genotype, data=subset(topt_noneg,rate.type=="GP")))
summary(IncGP.mod)
IncGP.mod

m3N<-aov(lnc~treatment+genotype, data=subset(topt_noneg,rate.type=="GP"))
TukeyHSD(m3N)
m3N

#visualize groups R
#
lncr<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=treatment, lnc, colour=treatment))+
  geom_boxplot()+
  geom_point()+
  ylab("Topt (ºC)")

esquisse:::ggplot_to_ppt("lncr")

lnnc2r<-ggplot(subset(topt_noneg,rate.type=='R'), aes(x=genotype, lnc, colour=genotype))+
  geom_boxplot()+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1)+
  geom_point()+
  ylab("lnc")+
  labs(title="Mean lnc across Genotype (R)")

esquisse::ggplot_to_ppt("lnnc2r")

RtbLT<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=treatment, lnc, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("lnc")+
  labs(title="Mean lnc across Treatment (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("RtbLT")


RtbLS<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=genotype, lnc, shape=genotype,colour=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("lnc")+
  labs(title="Mean lnc across Genotype (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("RtbLS")


#linear model
Inc.R.mod<-Anova(lm(lnc~treatment+genotype, data=subset(topt_noneg,rate.type=="R")))
Inc.R.mod
summary(Inc.R.mod)

m4N<-aov(lnc~treatment+genotype, data=subset(topt_noneg,rate.type=="R"))
TukeyHSD(m4N)
summary(m4N)




##### Statistical test for E GP

EGP<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=treatment, E, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("E")+
  labs(title="Mean E across Treatment (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("EGP")


EsGP<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=genotype, E, shape=genotype,colour=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("E")+
  labs(title="Mean E across Genotype (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("EsGP")


#linear model
E.GP.mod<-Anova(lm(E~treatment+genotype, data=subset(topt_noneg,rate.type=="GP")))
E.GP.mod
summary(E.GP.mod)

mEN<-aov(E~treatment+genotype, data=subset(topt_noneg,rate.type=="GP"))
TukeyHSD(mEN)
summary(mEN)


##### Statistical test for E R
#visualize groups R
ER<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=treatment, E, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("E")+
  labs(title="Mean E across Treatment (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("ER")


EsR<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=genotype, E, shape=genotype,colour=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("E")+
  labs(title="Mean E across Genotype (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("EsR")


#linear model
E.R.mod<-Anova(lm(E~treatment+genotype, data=subset(topt_noneg,rate.type=="R")))
E.R.mod
summary(E.R.mod)

mER<-aov(E~treatment+genotype, data=subset(topt_noneg,rate.type=="R"))
TukeyHSD(mER)
summary(mER)






##### Statistical test for Eh GP

EhGP<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=treatment, Eh, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Eh")+
  labs(title="Mean Eh across Treatment (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("EhGP")


EhsGP<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=genotype, Eh, shape=genotype,colour=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Eh")+
  labs(title="Mean Eh across Genotype (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("EhsGP")


#linear model
Eh.GP.mod<-Anova(lm(Eh~treatment+genotype, data=subset(topt_noneg,rate.type=="GP")))
Eh.GP.mod
summary(Eh.GP.mod)

mEhN<-aov(Eh~treatment+genotype, data=subset(topt_noneg,rate.type=="GP"))
TukeyHSD(mEhN)
summary(mEhN)


##### Statistical test for Eh R
#visualize groups R
EhR<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=treatment, Eh, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Eh")+
  labs(title="Mean Eh across Treatment (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("EhR")


EhsR<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=genotype, Eh, shape=genotype,colour=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Eh")+
  labs(title="Mean Eh across Genotype (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("EhsR")




#linear model
Eh.R.mod<-Anova(lm(Eh~treatment+genotype, data=subset(topt_noneg,rate.type=="R")))
Eh.R.mod
summary(Eh.R.mod)

mEhR<-aov(Eh~treatment+genotype, data=subset(topt_noneg,rate.type=="R"))
TukeyHSD(mEhR)
summary(mEhR)




##### Statistical test for Th GP

ThGP<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=treatment, Th, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Th")+
  labs(title="Mean Th across Treatment (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("ThGP")


ThsGP<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=genotype, Th, shape=genotype,colour=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Th")+
  labs(title="Mean Th across Genotype (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("ThsGP")


#linear model
Th.GP.mod<-Anova(lm(Th~treatment+genotype, data=subset(topt_noneg,rate.type=="GP")))
Th.GP.mod
summary(Th.GP.mod)

mThN<-aov(Th~treatment+genotype, data=subset(topt_noneg,rate.type=="GP"))
TukeyHSD(mThN)
summary(mThN)


##### Statistical test for Th R
#visualize groups R
ThR<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=treatment, y=Th, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=7)+
  ylab("Th")+
  labs(title="Mean Th across Treatment (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("ThR")


ThsR<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=genotype, Th, colour=treatment,shape=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Th")+
  labs(title="Mean Th across Genotype (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("ThsR")


#linear model
Th.R.mod<-Anova(lm(Th~treatment+genotype, data=subset(topt_noneg,rate.type=="R")))
Th.R.mod
summary(Th.R.mod)

mThR<-aov(Th~treatment+genotype, data=subset(topt_noneg,rate.type=="R"))
TukeyHSD(mThR)
summary(mThR)

### PMax
#calculate Pmax values between sites
## have to make dataframe long format----- FINAL FINAL
## 
## 


Pmax_datatest <- topt_noneg %>%
  #select(genotype,  term, estimate, treatment, rate.type) %>%
  #spread(term, estimate) %>%
  mutate(Pmax = schoolfield_high(lnc = lnc, E = E, Th = Th, Eh = Eh, temp = topt, Tc = 27)) %>% #add in factors that make up schoolfield function, reference topt to get pmax
  group_by(., genotype, rate.type, treatment, genotype)


write.csv(Pmax_datatest, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/TPC_params/Pmax_estimates.csv")


##### PMAX------
##### Statistical test for Pmax in GP
#visualize groups R
#
Pmaxgp<-ggplot(subset(Pmax_datatest,rate.type=='GP'), aes(x=genotype, Pmax, colour=genotype))+
  geom_boxplot()+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1)+
  geom_point()+
  ylab("Pmax")+
  xlab("Genotype")+
  labs(title="Mean Pmax across Genotype (GP)")

esquisse::ggplot_to_ppt("Pmaxgp")

Pmaxgp1<-ggplot(subset(Pmax_datatest,rate.type=='GP'), aes(x=treatment, Pmax, colour=treatment))+
  geom_boxplot()+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1)+
  geom_point()+
  ylab("Pmax")+
  xlab("Genotype")+
  labs(title="Mean Pmax across Treatments (GP)")

esquisse::ggplot_to_ppt("Pmaxgp1")

pmaxtreat<-ggplot(subset(Pmax_datatest,rate.type=="GP"), aes(x=treatment, Pmax, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Pmax")+
  labs(title="Mean Pmax across Depth (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("pmaxtreat")


pmaxsp<-ggplot(subset(Pmax_datatest,rate.type=="GP"), aes(x=genotype, Pmax, shape=treatment,colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Pmax")+
  labs(title="Mean Pmax across genotype (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("pmaxsp")


#linear model
Pmax.GP.mod<-Anova(lm(Pmax~treatment+genotype, data=subset(Pmax_datatest,rate.type=="GP")))
Pmax.GP.mod
summary(Pmax.GP.mod)

mPGP<-aov(Pmax~treatment+genotype, data=subset(Pmax_datatest,rate.type=="GP"))
TukeyHSD(mPGP)
summary(mPGP)

##### Statistical test for Pmax in R
#visualize groups R
Pmaxgprr<-ggplot(subset(Pmax_datatest,rate.type=='R'), aes(x=genotype, Pmax, colour=genotype))+
  geom_boxplot()+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1)+
  geom_point()+
  ylab("Pmax")+
  xlab("Genotype")+
  labs(title="Mean Pmax across Genotype (R)")

esquisse::ggplot_to_ppt("Pmaxgprr")

Pmaxgpr<-ggplot(subset(Pmax_datatest,rate.type=='R'), aes(x=treatment, Pmax, colour=treatment))+
  geom_boxplot()+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1)+
  geom_point()+
  ylab("Pmax")+
  xlab("Genotype")+
  labs(title="Mean Pmax across Treatments (R)")

esquisse::ggplot_to_ppt("Pmaxgpr")

Pmaxgpr2<-ggplot(subset(Pmax_datatest,rate.type=='R'), aes(x=treatment, Pmax, colour=genotype))+
  geom_boxplot()+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1)+
  geom_point()+
  ylab("Pmax")+
  xlab("Genotype")+
  labs(title="Mean Pmax across Treatments (R)")

esquisse::ggplot_to_ppt("Pmaxgpr2")


pmaxtrR<-ggplot(subset(Pmax_datatest,rate.type=="R"), aes(x=treatment, Pmax, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Pmax")+
  labs(title="Mean Pmax across Depth (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("pmaxtrR")


pmaxspR<-ggplot(subset(Pmax_datatest,rate.type=="R"), aes(x=genotype, Pmax, shape=treatment,colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Pmax")+
  labs(title="Mean Pmax across genotype (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("pmaxspR")



#linear model
Pmax.R.mod<-Anova(lm(Pmax~treatment+genotype, data=subset(Pmax_datatest,rate.type=="R")))
Pmax.R.mod
summary(Pmax.R.mod)

mPR<-aov(Pmax~treatment+genotype, data=subset(Pmax_datatest,rate.type=="R"))
TukeyHSD(mPR)
summary(mPR)



#anova function
Pmax.mod <- lm(Pmax~genotype*treatment, data=Pmax_datatest)
summary(Pmax.mod)
#Genotype M5!!! diff pmax (duh)

#check for normality, use normality plots

qqnorm(resid(Pmax.mod))
qqline(resid(Pmax.mod))

esquisse:::ggplot_to_ppt("NormPLot")

#check heteroscisity with boxplots

boxplot(resid(Pmax.mod)~Pmax_datatest$genotype*Pmax_datatest$treatment) 

#high R and low show inconsistent variances, may need to do weighted regression in the future

anova(Pmax.mod)
summary(Pmax.mod)
#intercept is significant- the intercept is different from zero
TukeyHSD(aov(Pmax.mod))
#compares genotype across 
anova(Pmax.mod)
summary(Pmax.mod)
TukeyHSD(aov(Pmax.mod))


Pmax.means<-Pmax_datatest %>%
  group_by(genotype,treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(Pmax), se=sd(Pmax)/sqrt(n())) #calculates mean and s.e.


write.csv(Pmax.means, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/TPC_params/Pmax_means.csv")



PMAX_P<-ggplot(Pmax.means, aes(x=rate.type, y=mean, col= genotype, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  labs(x="Rate Type", y="Maximum Rate of Performance (Pmax)", fill="treatment", color = "genotype") 

esquisse:::ggplot_to_ppt("PMAX_P")

PMAX_q<-ggplot(subset(Pmax.means, rate.type=="GP"), aes(x=genotype, y=mean, col= genotype, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=4) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  labs(x="genotype", y="Maximum Rate of Performance (Pmax) GP", fill="treatment", color = "genotype") 

esquisse:::ggplot_to_ppt("PMAX_q")


PMAX_R<-ggplot(subset(Pmax.means, rate.type=="R"), aes(x=genotype, y=mean, col= genotype, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=4) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  labs(x="genotype", y="Maximum Rate of Performance (Pmax) R", fill="treatment", color = "genotype") 

esquisse:::ggplot_to_ppt("PMAX_R")

PMAX_P2<-ggplot(Pmax.means, aes(x=rate.type, y=mean, col= treatment, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=3) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=1), width=0.5) +
  labs(x="Rate Type", y="Maximum Rate of Performance (Pmax)", fill="treatment", color = "treatment") 

esquisse:::ggplot_to_ppt("PMAX_P2")

ggsave(filename = "thermtol/Output/Pmax_graph.png", device = "png", width = 8, height = 5)




###### go to IndividualTPC's file to see TPC's for each genotype across each treatment and rate.type
###### 




