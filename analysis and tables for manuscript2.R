



ggline(mydataPlog.rate, "species", "log.rate",color= "treatment", add=c("mean_se", "dotplot"))
ggline(mydataRlog.rate, "species", "log.rate",color= "treatment", add=c("mean_se", "dotplot"))


boxplot(log.rate ~ species * treatment, data=mydataPlog.rate, frame = FALSE, 
        col = c("orange", "green", "purple", "blue", "orange", "green", "purple", "blue"), ylab="log.rate")

boxplot(log.rate ~ species * treatment, data=mydataRlog.rate, frame = FALSE, 
        col = c("orange", "green", "purple", "blue", "orange", "green", "purple", "blue"), ylab="log.rate")

interaction.plot(x.factor = mydataPlog.rate$species, trace.factor = mydataPlog.rate$treatment, 
                 response = mydataPlog.rate$log.rate, fun = mean, 
                 type = "b", legend = TRUE, 
                 xlab = "Species", ylab="log.rate")

interaction.plot(x.factor = mydataRlog.rate$species, trace.factor = mydataRlog.rate$treatment, 
                 response = mydataRlog.rate$log.rate, fun = mean, 
                 type = "b", legend = TRUE, 
                 xlab = "Species", ylab="log.rate")


res.aov3 <- aov(log.rate ~ species * treatment, data = mydataPlog.rate)
res.aov3r <- aov(log.rate~ species * treatment, data = mydataRlog.rate)


summary(res.aov3) 
summary(res.aov3r)

model.tables(res.aov3, type="means", se= TRUE)
model.tables(res.aov3r, type="means", se= TRUE)

#TukeyHSD(res.aov3)
#TukeyHSD(res.aov3r)

#TukeyHSD(res.aov3, which = "species") 
#TukeyHSD(res.aov3r, which = "species")
#
#summary(glht(res.aov3, linfct = mcp(species = "Tukey")))
#summary(glht(res.aov3r, linfct = mcp(species = "Tukey")))

pairwise.t.test(mydataPlog.rate$log.rate, mydataPlog.rate$species,
                p.adjust.method = "BH")
pairwise.t.test(mydataRlog.rate$log.rate, mydataRlog.rate$species,
                p.adjust.method = "BH")
##ANOVA assumes that the data are normally distributed and the variance across groups are homogeneous. We can check that with some diagnostic plots.

# 1. Homogeneity of variances
plot(res.aov3, 1) 
plot(res.aov3r, 1)
leveneTest(log.rate ~ species* treatment, data = mydataPlog.rate)
leveneTest(log.rate ~ species* treatment, data = mydataRlog.rate)
#no evidence to suggest that the variance across groups is statistically significantly different. Therefore, we can assume the homogeneity of variances in the different treatment groups.

# 2. Normality

plot(res.aov3, 2)
plot(res.aov3r, 2)

# Extract the residuals
aov_residuals3 <- residuals(object = res.aov3)
aov_residuals3r <- residuals(object = res.aov3r)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals3 )
shapiro.test(x = aov_residuals3r )

### make tables of pvalues for interaction GRoss photo

res.aov3<-aov(log.rate~species*treatment,mydataPlog.rate)
sstable<-Anova(res.aov3, type=2)
sstable
write.csv(sstable, "C:/Github/TPC_BDA_github/log.rateGPinPval.csv")

#TukeyHSD(res.aov3)
summary.aov(res.aov3, split = list(species=list(DLAB=0, MCAV=1, OFRA=2, PAST=3)))

contrasts <- summary.aov(res.aov3, split = list(species=list(DLAB=0, MCAV=1, OFRA=2, PAST=3)))[[1]][c(1,2,3, 4,5, 7), c(2, 1, 4, 5)]

# select the rows to combine
maineffects <- sstable[c(1,2), ]
me_contrasts <- contrasts[c(2,3,4,5), ]
interaction <- sstable[3, ]
resid <- sstable[4, ]

# bind the rows together in the desired order
sstableE <- rbind(maineffects, me_contrasts, interaction,  resid)
sstableE2 <- rbind(maineffects, me_contrasts, interaction,  resid)

sstableE$pes <- c(sstableE$'Sum Sq'[-nrow(sstableE)], NA)/(sstableE$'Sum Sq' + sstableE$'Sum Sq'[nrow(sstableE)]) # SS for each effect divided by the last SS (SS_residual)
View(sstableE)

options(knitr.kable.NA = '')

sstableE%>% 
  # Add a title
  kableExtra::kable( digits = 3,caption=" Anova log.rate Gross Photosynthesis ") %>% 
  # Specify how we will see the table, type ?kablestyling to see explanations and options
  kable_styling(bootstrap_options = "condensed", full_width=F) 


colnames(sstableE) <- c("SS", "df", "italic(F)", "italic(p)", "partial $\\eta^2$")

colnames(sstableE) <- c( "SS", "df", "$\\F$", "$\\alpha$p","expression(paste(*eta^2))")


rownames(sstableE) <- c("Species", "Depth", "D. labyrinthyformis", "M. cavernosa", "O. franksi", "P. astreoides", "Species:Depth", "Residuals")

kable(sstable, digits = 3, format = "pandoc", caption = "ANOVA table")

sstableE%>% 
  # Add a title
  kableExtra::kable( digits = 3,caption=" Anova log.rate Gross Photosynthesis ", col.names = c ("SS", "df", "F", "p-value","eta-squared")) %>% 
  kable_styling(bootstrap_options = "condensed", full_width=F) 

sstableE2%>% 
  # Add a title
  kableExtra::kable( digits = 3,caption=" Anova log.rate Gross Photosynthesis ", col.names = c ("SS", "df", "F", "p-value")) %>% 
  kable_styling(bootstrap_options = "condensed", full_width=F) 

### TABLES for Pvalus for respiration (E)
res.aov3<-aov(log.rate~species*treatment,mydataRlog.rate)
sstable<-Anova(res.aov3, type=2)
sstable
write.csv(sstable, "C:/Github/TPC_BDA_github/log.raterinPval.csv")

#TukeyHSD(res.aov3)
summary.aov(res.aov3, split = list(species=list(DLAB=0, MCAV=1, OFRA=2, PAST=3)))

contrasts <- summary.aov(res.aov3, split = list(species=list(DLAB=0, MCAV=1, OFRA=2, PAST=3)))[[1]][c(1,2,3, 4,5, 7), c(2, 1, 4, 5)]

# select the rows to combine
maineffects <- sstable[c(1,2), ]
me_contrasts <- contrasts[c(2,3,4,5), ]
interaction <- sstable[3, ]
resid <- sstable[4, ]

# bind the rows together in the desired order
sstableE <- rbind(maineffects, me_contrasts, interaction,  resid)
sstableE2 <- rbind(maineffects, me_contrasts, interaction,  resid)


sstableE$pes <- c(sstableE$'Sum Sq'[-nrow(sstableE)], NA)/(sstableE$'Sum Sq' + sstableE$'Sum Sq'[nrow(sstableE)]) # SS for each effect divided by the last SS (SS_residual)
sstableE

options(knitr.kable.NA = '')

sstableE%>% 
  # Add a title
  kableExtra::kable( digits = 3,caption=" Anova log.rate Gross Photosynthesis ") %>% 
  # Specify how we will see the table, type ?kablestyling to see explanations and options
  kable_styling(bootstrap_options = "condensed", full_width=F) 


colnames(sstableE) <- c("SS", "df", "F", "p.value", "(eta^2)")

rownames(sstableE) <- c("Species", "Depth", "D. labyrinthyformis", "M. cavernosa", "O. franksi", "P. astreoides", "Species:Depth", "Residuals")

kable(sstable, digits = 3, format = "pandoc", caption = "ANOVA table")

sstableE%>% 
  # Add a title
  kableExtra::kable( digits = 3,caption=" Anova log.rate Gross Photosynthesis ", col.names = c ("SS", "df", "F", "p-value","eta-squared")) %>% 
  kable_styling(bootstrap_options = "condensed", full_width=F) 

sstableE2%>% 
  # Add a title
  kableExtra::kable( digits = 3,caption=" Anova log.rate Gross Photosynthesis ", col.names = c ("SS", "df", "F", "p-value")) %>% 
  kable_styling(bootstrap_options = "condensed", full_width=F) 
```


Topt by species
```{r}
# Topt by species

# GP species


# model
mPGPsp<-aov(log.rate~species,mydataP)
summary(mPGPsp)

## check for outliers
cooksd <- cooks.distance(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")

mydataPlog.rates<-mydataP[-(15),]

#use mydataPE
mPGPsp<-aov(log.rate~species,mydataPlog.rates)
summary(mPGPsp)

# log.rate comparisons between  species for respiration

# model
mPGPsp<-aov(log.rate~species, mydataR)
summary(mPGPsp)

## check for outliers
cooksd <- cooks.distance(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")

mydataRlog.rates<-mydataR[-(15),]


mPGPsp<-aov(log.rate~species, mydataRlog.rates)
summary(mPGPsp)

#

ggline(mydataPlog.rates, "species", "log.rate",color= "treatment", add=c("mean_se", "dotplot"))
ggline(mydataRlog.rates, "species", "log.rate",color= "treatment", add=c("mean_se", "dotplot"))



boxplot(log.rate ~ species, data=mydataPlog.rates, frame = FALSE, 
        col = c("orange", "green", "purple", "blue"), ylab="log.rate")

boxplot(log.rate ~ species, data=mydataRlog.rates, frame = FALSE, 
        col = c("orange", "green", "purple", "blue"), ylab="log.rate")


res.aovagp <- aov(log.rate ~ species, data = mydataPlog.rates)
res.aovar <- aov(log.rate ~ species, data = mydataRlog.rates)

summary(res.aovagp)
summary(res.aovar)

model.tables(res.aovagp, type="means", se= TRUE)
model.tables(res.aovar, type="means", se= TRUE)



TukeyHSD(res.aovagp, which = "species")
#TukeyHSD(res.aovar, which = "species")


#summary(glht(res.aovagp, linfct = mcp(species = "Tukey")))
#summary(glht(res.aovar, linfct = mcp(species = "Tukey")))

pairwise.t.test(mydataPlog.rates$log.rate, mydataPlog.rates$species,
                p.adjust.method = "BH")
pairwise.t.test(mydataRlog.rates$log.rate, mydataRlog.rates$species,
                p.adjust.method = "BH")
##ANOVA assumes that the data are normally distributed and the variance across groups are homogeneous. We can check that with some diagnostic plots.

# 1. Homogeneity of variances
plot(res.aovagp, 1)
plot(res.aovar, 1)

leveneTest(log.rate ~ species, data = mydataPlog.rates)
leveneTest(log.rate ~ species, data = mydataRlog.rates)

#no evidence to suggest that the variance across groups is statistically significantly different. Therefore, we can assume the homogeneity of variances in the different treatment groups.

# 2. Normality
plot(res.aovagp, 2)
plot(res.aovar, 2)
# Extract the residuals
aov_residualsagp <- residuals(object = res.aovagp)
aov_residualsar <- residuals(object = res.aovar)
# Run Shapiro-Wilk test

shapiro.test(x = aov_residualsagp )
shapiro.test(x = aov_residualsar )

bartlett.test(log.rate ~ species,
              data = mydataPlog.rates) 
oneway.test(log.rate ~ species,
            data = mydataPlog.rates, var.equal=FALSE)
kruskal.test(log.rate ~ species,
             data = mydataPlog.rates)

bartlett.test(log.rate ~ species,
              data = mydataRlog.rates) 
oneway.test(log.rate ~ species,
            data = mydataRlog.rates, var.equal=FALSE)
kruskal.test(log.rate ~ species,
             data = mydataRlog.rates)
```


log.rate by treatment
```{r}
# log.rate by treatment

# GP treatment


# model
mPGPsp<-aov(log.rate~treatment,mydataP)
summary(mPGPsp)

## check for outliers
cooksd <- cooks.distance(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")
#use 
mydataPlog.ratet<-mydataP[-(15),]
mPGPsp<-aov(log.rate~treatment,mydataPlog.ratet)
summary(mPGPsp)


# log.rate comparisons between  treatment for respiration

# model
mPGPsp<-aov(log.rate~treatment, mydataR)
summary(mPGPsp)

## check for outliers
cooksd <- cooks.distance(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")

#
mydataRlog.rates2<-mydataP[-(4),]

mPGPsp<-aov(log.rate~treatment,mydataRlog.rates2)
summary(mPGPsp)

#

ggline(mydataPlog.ratet, "treatment", "E",color= "treatment", add=c("mean_se", "dotplot"))
ggline(mydataRlog.rates2, "treatment", "E",color= "treatment", add=c("mean_se", "dotplot"))



boxplot(log.rate ~ treatment, data=mydataPlog.ratet, frame = FALSE, 
        col = c("green4", "green1"), ylab="log.rate")

boxplot(log.rate ~ treatment, data=mydataRlog.rates2, frame = FALSE, 
        col = c("blue4", "blue1"), ylab="log.rate")


res.aovagp <- aov(log.rate ~ treatment, data = mydataPlog.ratet)
res.aovar <- aov(log.rate ~ treatment, data = mydataRlog.rates2)

summary(res.aovagp)
summary(res.aovar)


model.tables(res.aovagp, type="means", se= TRUE)
model.tables(res.aovar, type="means", se= TRUE)


#TukeyHSD(res.aovagp, which = "treatment")
#TukeyHSD(res.aovar, which = "treatment")


#summary(glht(res.aovagp, linfct = mcp(treatment = "Tukey")))
#summary(glht(res.aovar, linfct = mcp(treatment = "Tukey")))

pairwise.t.test(mydataPlog.ratet$log.rate, mydataPlog.ratet$treatment,
                p.adjust.method = "BH")
pairwise.t.test(mydataRlog.rates$log.rate, mydataRlog.rates2$treatment,
                p.adjust.method = "BH")

##ANOVA assumes that the data are normally distributed and the variance across groups are homogeneous. We can check that with some diagnostic plots.

# 1. Homogeneity of variances
plot(res.aovagp, 1)
plot(res.aovar, 1)

leveneTest(log.rate ~ treatment, data = mydataPlog.ratet)
leveneTest(log.rate ~ treatment, data = mydataRlog.rates2)

#no evidence to suggest that the variance across groups is statistically significantly different. Therefore, we can assume the homogeneity of variances in the different treatment groups.

# 2. Normality
plot(res.aovagp, 2)
plot(res.aovar, 2)


# Extract the residuals
aov_residualsagp <- residuals(object = res.aovagp)
aov_residualsar <- residuals(object = res.aovar)
# Run Shapiro-Wilk test

shapiro.test(x = aov_residualsagp )
shapiro.test(x = aov_residualsar )

###
kruskal.test(log.rate~treatment,mydataPlog.ratet)
kruskal.test(log.rate~treatment,mydataRlog.rates2)

#bartlett.test(log.rate ~ treatment,
#              data = mydataPlog.ratet) 
#bartlett.test(log.rate ~ treatment,
#              data = mydataRlog.rates) 
#oneway.test(log.rate ~ treatment,
#            data = mydataPlog.ratet, var.equal=FALSE)
#oneway.test(log.rate ~ treatment,
#            data = mydataRlog.rates, var.equal=FALSE)

```



old code HISTORY:
  
  
  
  
  
  
  summary(res.aovagp)
summary(res.aovar)
#TukeyHSD(res.aov2, which = "species")
#TukeyHSD(res.aov2r, which = "species")
#TukeyHSD(res.aov3, which = "species") #ShaowspR MCAV-DLAB, PAST-DLAB,OFRA-MCAV, PAST-OFRA
#TukeyHSD(res.aov3r, which = "species")
TukeyHSD(res.aovagp, which = "species")
TukeyHSD(res.aovar, which = "species")
#summary(glht(res.aov2, linfct = mcp(species = "Tukey")))
#summary(glht(res.aov2r, linfct = mcp(species = "Tukey")))
#summary(glht(res.aov3, linfct = mcp(species = "Tukey")))
#summary(glht(res.aov3r, linfct = mcp(species = "Tukey")))
summary(glht(res.aovagp, linfct = mcp(species = "Tukey")))
summary(glht(res.aovar, linfct = mcp(species = "Tukey")))
#gross photo
# model
mPGPsp<-aov(lnc~species*treatment,testGP)
summary(mPGPsp)
## check for outliers
cooksd <- cooks.distance(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")
#remove 13,15
testGPlnc<-testGP[-(15),]
mPGPsp<-aov(lnc~species*treatment, testGPlnc)
summary(mPGPsp)
# model
mPGPsp<-aov(lnc~species*treatment, testR)
summary(mPGPsp)
## check for outliers
cooksd <- cooks.distance(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")
#no outliers
testRlnc<-testR
mPGPsp<-aov(lnc~species*treatment, testRlnc)
summary(mPGPsp)
ggline(testGPlnc, "species", "lnc",color= "treatment", add=c("mean_se", "dotplot"))
ggline(testRlnc, "species", "lnc",color= "treatment", add=c("mean_se", "dotplot"))
boxplot(lnc ~ species * treatment, data=testGPlnc, frame = FALSE,
        col = c("orange", "green", "purple", "blue", "orange", "green", "purple", "blue"), ylab="lnc")
boxplot(lnc ~ species * treatment, data=testRlnc, frame = FALSE,
        col = c("orange", "green", "purple", "blue", "orange", "green", "purple", "blue"), ylab="lnc")
interaction.plot(x.factor = testGPlnc$species, trace.factor = testGPlnc$treatment,
                 response = testGPlnc$lnc, fun = mean,
                 type = "b", legend = TRUE,
                 xlab = "Species", ylab="lnc")
interaction.plot(x.factor = testRlnc$species, trace.factor = testRlnc$treatment,
                 response = testRlnc$lnc, fun = mean,
                 type = "b", legend = TRUE,
                 xlab = "Species", ylab="lnc")
res.aov3 <- aov(lnc ~ species * treatment, data = testGPlnc)
res.aov3r <- aov(lnc~ species * treatment, data = testRlnc)
summary(res.aov3)
summary(res.aov3r)
model.tables(res.aov3, type="means", se= TRUE)
model.tables(res.aov3r, type="means", se= TRUE)
pairwise.t.test(testGPlnc$lnc, testGPlnc$species,
                p.adjust.method = "BH")
pairwise.t.test(testRlnc$lnc, testRlnc$species,
                p.adjust.method = "BH")
# 1. Homogeneity of variances
plot(res.aov3, 1)
plot(res.aov3r, 1)
leveneTest(lnc ~ species* treatment, data = testGPlnc)
leveneTest(lnc ~ species* treatment, data = testRlnc)
plot(res.aov3, 2)
plot(res.aov3r, 2)
# Extract the residuals
aov_residuals3 <- residuals(object = res.aov3)
aov_residuals3r <- residuals(object = res.aov3r)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals3 )
shapiro.test(x = aov_residuals3r )
kruskal.test(lnc~species*treatment,testGPlnc)
res.aov3<-aov(lnc~species*treatment,testGPlnc)
sstable<-Anova(res.aov3, type=2)
sstable
write.csv(sstable, "C:/Github/TPC_BDA_github/lncGPinPval.csv")
#TukeyHSD(res.aov3)
summary.aov(res.aov3, split = list(species=list(DLAB=0, MCAV=1, OFRA=2, PAST=3)))
contrasts <- summary.aov(res.aov3, split = list(species=list(DLAB=0, MCAV=1, OFRA=2, PAST=3)))[[1]][c(1,2,3, 4,5, 7), c(2, 1, 4, 5)]
# select the rows to combine
maineffects <- sstable[c(1,2), ]
me_contrasts <- contrasts[c(2,3,4,5), ]
interaction <- sstable[3, ]
resid <- sstable[4, ]
# bind the rows together in the desired order
sstableE <- rbind(maineffects, me_contrasts, interaction,  resid)
sstableE2 <- rbind(maineffects, me_contrasts, interaction,  resid)
sstableE$pes <- c(sstableE$'Sum Sq'[-nrow(sstableE)], NA)/(sstableE$'Sum Sq' + sstableE$'Sum Sq'[nrow(sstableE)]) # SS for each effect divided by the last SS (SS_residual)
options(knitr.kable.NA = '')
sstableE%>%
  # Add a title
  kableExtra::kable( digits = 3,caption=" Anova lnc Gross Photosynthesis ") %>%
  # Specify how we will see the table, type ?kablestyling to see explanations and options
  kable_styling(bootstrap_options = "condensed", full_width=F)
colnames(sstableE) <- c("SS", "df", "italic(F)", "italic(p)", "partial $\\eta^2$")
colnames(sstableE) <- c( "SS", "df", "$\\F$", "$\\alpha$p","expression(paste(*eta^2))")
rownames(sstableE) <- c("Species", "Depth", "D. labyrinthyformis", "M. cavernosa", "O. franksi", "P. astreoides", "Species:Depth", "Residuals")
kable(sstable, digits = 3, format = "pandoc", caption = "ANOVA table")
sstableE%>%
  # Add a title
  kableExtra::kable( digits = 3,caption=" Anova lnc Gross Photosynthesis ", col.names = c ("SS", "df", "F", "p-value","eta-squared")) %>%
  kable_styling(bootstrap_options = "condensed", full_width=F)
sstableE2%>%
  # Add a title
  kableExtra::kable( digits = 3,caption=" Anova lnc Gross Photosynthesis ", col.names = c ("SS", "df", "F", "p-value")) %>%
  kable_styling(bootstrap_options = "condensed", full_width=F)
### TABLES for Pvalus for respiration (E)
res.aov3<-aov(lnc~species*treatment,testRlnc)
sstable<-Anova(res.aov3, type=2)
sstable
write.csv(sstable, "C:/Github/TPC_BDA_github/lncrinPval.csv")
#TukeyHSD(res.aov3)
summary.aov(res.aov3, split = list(species=list(DLAB=0, MCAV=1, OFRA=2, PAST=3)))
contrasts <- summary.aov(res.aov3, split = list(species=list(DLAB=0, MCAV=1, OFRA=2, PAST=3)))[[1]][c(1,2,3, 4,5, 7), c(2, 1, 4, 5)]
# select the rows to combine
maineffects <- sstable[c(1,2), ]
me_contrasts <- contrasts[c(2,3,4,5), ]
interaction <- sstable[3, ]
resid <- sstable[4, ]
# bind the rows together in the desired order
sstableE <- rbind(maineffects, me_contrasts, interaction,  resid)
sstableE2 <- rbind(maineffects, me_contrasts, interaction,  resid)
sstableE$pes <- c(sstableE$'Sum Sq'[-nrow(sstableE)], NA)/(sstableE$'Sum Sq' + sstableE$'Sum Sq'[nrow(sstableE)]) # SS for each effect divided by the last SS (SS_residual)
sstableE
options(knitr.kable.NA = '')
sstableE%>%
  # Add a title
  kableExtra::kable( digits = 3,caption=" Anova lnc Gross Photosynthesis ") %>%
  # Specify how we will see the table, type ?kablestyling to see explanations and options
  kable_styling(bootstrap_options = "condensed", full_width=F)
colnames(sstableE) <- c("SS", "df", "F", "p.value", "(eta^2)")
rownames(sstableE) <- c("Species", "Depth", "D. labyrinthyformis", "M. cavernosa", "O. franksi", "P. astreoides", "Species:Depth", "Residuals")
kable(sstable, digits = 3, format = "pandoc", caption = "ANOVA table")
sstableE%>%
  # Add a title
  kableExtra::kable( digits = 3,caption=" Anova lnc Gross Photosynthesis ", col.names = c ("SS", "df", "F", "p-value","eta-squared")) %>%
  kable_styling(bootstrap_options = "condensed", full_width=F)
sstableE2%>%
  # Add a title
  kableExtra::kable( digits = 3,caption=" Anova lnc Gross Photosynthesis ", col.names = c ("SS", "df", "F", "p-value")) %>%
  kable_styling(bootstrap_options = "condensed", full_width=F)
# model
mPGPsp<-aov(lnc~species,testGP)
summary(mPGPsp)
## check for outliers
cooksd <- cooks.distance(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")
testGPlncs<-testGP[-(15),]
#use testGPE
mPGPsp<-aov(lnc~species,testGPlncs)
summary(mPGPsp)
# model
mPGPsp<-aov(lnc~species, testR)
summary(mPGPsp)
# model
mPGPsp<-aov(lnc~species, testR)
summary(mPGPsp)
## check for outliers
cooksd <- cooks.distance(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")
testRlncs<-testR[-(15),]
mPGPsp<-aov(lnc~species, testRlncs)
summary(mPGPsp)
ggline(testGPlncs, "species", "lnc",color= "treatment", add=c("mean_se", "dotplot"))
ggline(testRlncs, "species", "lnc",color= "treatment", add=c("mean_se", "dotplot"))
boxplot(lnc ~ species, data=testGPlncs, frame = FALSE,
        col = c("orange", "green", "purple", "blue"), ylab="lnc")
boxplot(lnc ~ species, data=testRlncs, frame = FALSE,
        col = c("orange", "green", "purple", "blue"), ylab="lnc")
res.aovagp <- aov(lnc ~ species, data = testGPlncs)
res.aovar <- aov(lnc ~ species, data = testRlncs)
summary(res.aovagp)
summary(res.aovar)
model.tables(res.aovagp, type="means", se= TRUE)
summary(res.aovagp)
summary(res.aovar)
model.tables(res.aovagp, type="means", se= TRUE)
model.tables(res.aovar, type="means", se= TRUE)
TukeyHSD(res.aovagp, which = "species")
pairwise.t.test(testGPlncs$lnc, testGPlncs$species,
                p.adjust.method = "BH")
pairwise.t.test(testRlncs$lnc, testRlncs$species,
                p.adjust.method = "BH")
# 1. Homogeneity of variances
plot(res.aovagp, 1)
plot(res.aovar, 1)
leveneTest(lnc ~ species, data = testGPlncs)
leveneTest(lnc ~ species, data = testRlncs)
# 2. Normality
plot(res.aovagp, 2)
plot(res.aovar, 2)
# Extract the residuals
aov_residualsagp <- residuals(object = res.aovagp)
aov_residualsar <- residuals(object = res.aovar)
shapiro.test(x = aov_residualsagp )
shapiro.test(x = aov_residualsar )
bartlett.test(lnc ~ species,
              data = testGPlncs)
oneway.test(lnc ~ species,
            data = testGPlncs, var.equal=FALSE)
kruskal.test(lnc ~ species,
             data = testGPlncs)
bartlett.test(lnc ~ species,
              data = testRlncs)
oneway.test(lnc ~ species,
            data = testRlncs, var.equal=FALSE)
kruskal.test(lnc ~ species,
             data = testRlncs)





# model
mPGPsp<-aov(lnc~treatment,testGP)
summary(mPGPsp)
## check for outliers
cooksd <- cooks.distance(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")
#use
testGPlnct<-testGP[-(15),]
mPGPsp<-aov(lnc~treatment,testGPlnct)
summary(mPGPsp)
# model
mPGPsp<-aov(lnc~treatment, testR)
summary(mPGPsp)
## check for outliers
cooksd <- cooks.distance(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")
#testRlncs
testRlncs<-testGP[-(4),]
mPGPsp<-aov(lnc~treatment,testRlncs)
summary(mPGPsp)
testRlncs<-testR[-(15),]
#
testRlncs2<-testGP[-(4),]
mPGPsp<-aov(lnc~treatment,testRlncs2)
summary(mPGPsp)
ggline(testGPlnct, "treatment", "E",color= "treatment", add=c("mean_se", "dotplot"))
ggline(testRlncs2, "treatment", "E",color= "treatment", add=c("mean_se", "dotplot"))
boxplot(lnc ~ treatment, data=testGPlnct, frame = FALSE,
        col = c("green4", "green1"), ylab="lnc")
boxplot(lnc ~ treatment, data=testRlncs2, frame = FALSE,
        col = c("blue4", "blue1"), ylab="lnc")
res.aovagp <- aov(lnc ~ treatment, data = testGPlnct)
res.aovar <- aov(lnc ~ treatment, data = testRlncs2)
summary(res.aovagp)
summary(res.aovar)
model.tables(res.aovagp, type="means", se= TRUE)
model.tables(res.aovar, type="means", se= TRUE)
pairwise.t.test(testGPlnct$lnc, testGPlnct$treatment,
                p.adjust.method = "BH")
pairwise.t.test(testRlncs$lnc, testRlncs2$treatment,
                p.adjust.method = "BH")
# 1. Homogeneity of variances
plot(res.aovagp, 1)
plot(res.aovar, 1)
leveneTest(lnc ~ treatment, data = testGPlnct)
leveneTest(lnc ~ treatment, data = testRlncs2)
leveneTest(lnc ~ treatment, data = testRlncs2)
# 2. Normality
plot(res.aovagp, 2)
plot(res.aovar, 2)
# Extract the residuals
aov_residualsagp <- residuals(object = res.aovagp)
aov_residualsar <- residuals(object = res.aovar)
shapiro.test(x = aov_residualsagp )
shapiro.test(x = aov_residualsar )
###
kruskal.test(lnc~treatment,testGPlnct)
kruskal.test(lnc~treatment,testRlncs2)
###
mPGPsp<-aov(Pmax~species*treatment,testGP)
summary(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")
###
mPGPsp<-aov(Pmax~species*treatment,testGP)
summary(mPGPsp)
## check for outliers
cooksd <- cooks.distance(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")
#remove1,2
testGPP<-testGP[-c(1,2),]
mPGPsp<-aov(Pmax~species*treatment, testGPP)
summary(mPGPsp)
# model
mPGPsp<-aov(Pmax~species*treatment, testR)
summary(mPGPsp)
## check for outliers
cooksd <- cooks.distance(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")
#gross photo
# model
mPGPsp<-aov(log.rate~species*treatment,mydataP)
summary(mPGPsp)
#gross photo
# model
mPGPsp<-aov(log.rate~species*treatment,mydataP)
summary(mPGPsp)
## check for outliers
cooksd <- cooks.distance(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")
#gross photo
# model
mPGPsp<-aov(log.rate~species*treatment*temp.Cat,mydataP)
summary(mPGPsp)
## check for outliers
cooksd <- cooks.distance(mPGPsp)
plot(cooksd, pch=".", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=TRUE), col="red")  # add cutoff line--In general use, those observations that have a cook’s distance greater than 4 times the mean may be classified as influential. This is not a hard boundary.
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=TRUE),names(cooksd),""), cex=0.6, pos=2, col="red")
fit<-lm(log.rate~species*treatment*temp.Cat ,mydataR)
summary(fit)
Anova(fit)
plot_model(fit, type = "pred", terms = c("temp.Cat", "treatment", "species"))

