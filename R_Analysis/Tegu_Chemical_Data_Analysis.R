### SOIL BIOGEOCHEMICAL ANALYSIS ###

#***********************************

#Library of Packages used 
library(dplyr)
library(ggplot2)
library(Rmisc)
library(lme4)
library(emmeans)
library(ggpubr)
library(ggh4x)

#add the data file 
chemdf <- read.csv("Tegu_FinalMetadata.csv")

#time to filter so only looking at the soil data sans swabs 
Removed <- subset(chemdf, Stage != "Void ")
SoilChemData <-filter(Removed,Source == "Soil ")

#filter into either Tegu decomposition impacted or control soils 
TeguChemSoil <-filter(SoilChemData, Location == "Tegu")
ConChemSoil <-filter(SoilChemData, Location == "Control")

#change decomposition stage, location, and sample names to factors 
SoilChemData$Stage <- factor(SoilChemData$Stage, levels = c("Early", "Bloat", "Active I", "Active II", "Active III ", "Advanced I", "Advanced II", "Advanced III", "Skeletal I ", "Skeletal II"))
SoilChemData$Location <- factor(SoilChemData$Location, levels=c("Control","Tegu"))
SoilChemData$SampleName <- factor(SoilChemData$SampleName, levels=c("T0C", "T1C", "T2C", "T3C", "T4C", "T5C", "T6C", "T7C", "T8C", "T9C", "T10C", "T11C", "T13C", "T14C", "T15C", "T16C", "T17C", "T19C", "T20C", "T21C", 
                                                      "T22C", "T23C", "T24C", "T25C", "T26C", "T27C", "T28C", "T29C", "T30C", "T31C", "T0T", "T1T", "T2T", "T3T", "T4T", "T5T", "T6T", "T7T", "T8T", "T9T", "T10T", "T11T", "T13T", "T14T", "T15T", "T16T", "T17T",
                                                      "T19T", "T20T", "T21T", "T22T", "T23T", "T24T", "T25T", "T26T", "T27T", "T28T", "T29T","T30T", "T31T"))
#End of Data Staging 

##TOTAL N PER SAMPLE ANALYSIS  

tndf <- SoilChemData %>% rename_at('N', ~'Npercent') #makes a dataframe where the column name is changed 
tndf2 <- summarySE(tndf, measurevar="Npercent", groupvars=c("Location","Stage")) #this gives us the summary of the data 


svg("NpercentLine.svg", width = 9, height = 4)
ggplot(tndf2, aes(x=Stage, y=Npercent, color=Location, group=Location)) + 
  geom_errorbar(aes(ymin=Npercent-sd, ymax=Npercent+sd), width=.1) +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1)) +
  scale_color_manual(values = c("black","red")) +
  geom_line(size=1.2) +
  geom_point(aes(shape = Location), size=4) +
  scale_shape_manual(values = c(16, 15)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 20), legend.text = element_text(size=18), legend.title = element_text(size=20)) +
  ylab("%N per sample") +
  xlab("Decomposition Stage")  
dev.off()

#look at distribution of this variable with Shapiro Wilks test and visually 
shapiro.test(TeguChemSoil$N)
hist(TeguChemSoil$N)
qqnorm(TeguChemSoil$N)
qqline(TeguChemSoil$N)
shapiro.test(ConChemSoil$N)
hist(ConChemSoil$N)
qqnorm(ConChemSoil$N)
qqline(ConChemSoil$N)

#do non-parametric comparison across all time points 
kruskal.test(N ~ Location, data = SoilChemData) #significant difference (0.004421)
#also do two sample t.test comparisons at time points with ggpubr package, Holm's correction default 
compare_means(N ~ Location, data = SoilChemData, method = "t.test", paired = FALSE, group.by = "Stage")

##TOTAL C PER SAMPLE ANALYSIS 

tcdf <- SoilChemData %>% rename_at('C', ~'Cpercent') #makes a dataframe where the column name is changed 
tcdf2 <- summarySE(tcdf, measurevar="Cpercent", groupvars=c("Location","Stage")) #this gives a summary of the data 

svg("CpercentLine.svg", width = 9, height = 4)
ggplot(tcdf2, aes(x=Stage, y=Cpercent, color=Location, group=Location)) + 
  geom_errorbar(aes(ymin=Cpercent-sd, ymax=Cpercent+sd), width=.1) +
  coord_cartesian(ylim = c(0, 8)) +
  scale_y_continuous(breaks = seq(0, 8, by = 2)) +
  scale_color_manual(values = c("black","red")) +
  geom_line(size=1.2) +
  geom_point(aes(shape = Location), size=4) +
  scale_shape_manual(values = c(16, 15)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 20), legend.text = element_text(size=18), legend.title = element_text(size=20)) +
  ylab("%C per sample") +
  xlab("Decomposition Stage")
dev.off()

#look at distribution of this variable with Shapiro Wilks test and visually 
shapiro.test(TeguChemSoil$C)
hist(TeguChemSoil$C)
qqnorm(TeguChemSoil$C)
qqline(TeguChemSoil$C)
shapiro.test(ConChemSoil$C)
hist(ConChemSoil$C)
qqnorm(ConChemSoil$C)
qqline(ConChemSoil$C)

#do non-parametric comparison across all time points
kruskal.test(C ~ Location, data = SoilChemData) #no significance (0.1558) 
#also do two sample t.test comparisons at time points with ggpubr package, Holm's correction default 
compare_means(C ~ Location, data = SoilChemData, method = "t.test", paired = FALSE, group.by = "Stage")

##C:N PER SAMPLE ANALYSIS 

CNdf <- summarySE(SoilChemData, measurevar="CtoN", groupvars=c("Location","Stage")) #This gives a summery of the C:N data  

svg("CNLine.svg", width = 9, height = 4)
ggplot(CNdf, aes(x=Stage, y=CtoN, color=Location, group=Location)) + 
  geom_errorbar(aes(ymin=CtoN-sd, ymax=CtoN+sd), width=.1) +
  coord_cartesian(ylim = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 5)) +
  scale_color_manual(values = c("black","red")) +
  geom_line(size=1.2) +
  geom_point(aes(shape = Location), size=4) +
  scale_shape_manual(values = c(16, 15)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 20), legend.text = element_text(size=18), legend.title = element_text(size=20)) +
  ylab("C:N") +
  xlab("Decomposition Stage")
dev.off()

#look at distribution of this variable with Shapiro Wilks test and visually 
shapiro.test(TeguChemSoil$CtoN)
hist(TeguChemSoil$CtoN)
qqnorm(TeguChemSoil$CtoN)
qqline(TeguChemSoil$CtoN)
shapiro.test(ConChemSoil$CtoN)
hist(ConChemSoil$CtoN)
qqnorm(ConChemSoil$CtoN)
qqline(ConChemSoil$CtoN)

#do non-parametric comparison across all time points
kruskal.test(CtoN ~ Location, data = SoilChemData) #Significant difference (1.774e-9)
#also do two sample t.test comparisons at time points with ggpubr package, Holm's correction default 
compare_means(CtoN ~ Location, data = SoilChemData, method = "t.test", paired = FALSE, group.by = "Stage")

##AMMONIUM CONCENTRATION PER SAMPLE 

NH4ugdf <- summarySE(SoilChemData, measurevar="NH4ug", groupvars=c("Location","Stage")) #provides a summary of the ammonium data

svg("NH4ugLine.svg", width = 9, height = 4)
ggplot(NH4ugdf, aes(x=Stage, y=NH4ug, color=Location, group=Location)) + 
  geom_errorbar(aes(ymin=NH4ug-se, ymax=NH4ug+se), width=.1) +
  scale_y_continuous(breaks = seq(0, 800, by = 100)) + 
  scale_color_manual(values = c("black","red")) +
  geom_line(size=1.2) +
  geom_point(aes(shape = Location), size=4) +
  scale_shape_manual(values = c(16, 15)) + 
  force_panelsizes(rows = unit(2, "in"), cols = unit(6, "in")) +
  labs(color = 'Treatment', shape='Treatment') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 20), legend.text = element_text(size=18), legend.title = element_text(size=20)) +
  labs (x="Decomposition Stage", y = expression(NH[4]+ ~ (mu * g ~ gdw^-1)))
dev.off()

#look at distribution of this variable with Shapiro Wilks test and visually 
shapiro.test(TeguChemSoil$NH4ug)
hist(TeguChemSoil$NH4ug)
qqnorm(TeguChemSoil$NH4ug)
qqline(TeguChemSoil$NH4ug)
shapiro.test(ConChemSoil$NH4ug)
hist(ConChemSoil$NH4ug)
qqnorm(ConChemSoil$NH4ug)
qqline(ConChemSoil$NH4ug)

#do non-parametric comparison across all time points
kruskal.test(NH4ug ~ Location, data = SoilChemData) #Significant difference (1.369e-8)
#also do two sample t.test comparisons at time points with ggpubr package, Holm's correction default 
compare_means(NH4ug ~ Location, data = SoilChemData, method = "t.test", paired = FALSE, group.by = "Stage")

##NITRATE CONCENTRATION PER SAMPLE 

NO3ugdf <- summarySE(SoilChemData, measurevar="NO3ug", groupvars=c("Location","Stage"))

svg("NO3ugLine.svg", width = 9, height = 4)
ggplot(NO3ugdf, aes(x=Stage, y=NO3ug, color=Location, group=Location)) + 
  geom_errorbar(aes(ymin=NO3ug-se, ymax=NO3ug+se), width=.1) +
  coord_cartesian(ylim = c(0, 700)) +
  scale_y_continuous(breaks = seq(0, 700, by = 100)) +
  scale_color_manual(values = c("black","red")) +
  geom_line(size=1.2) +
  geom_point(aes(shape = Location), size=4) +
  scale_shape_manual(values = c(16, 15)) + 
  force_panelsizes(rows = unit(2, "in"), cols = unit(6, "in")) +
  labs(color = 'Treatment', shape='Treatment') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 20), legend.text = element_text(size=18), legend.title = element_text(size=20)) +
  labs (x="Decomposition Stage", y = expression(NO[3]- ~ (mu * g ~ gdw^-1)))
dev.off()

#look at distribution of this variable with Shapiro Wilks test and visually 
shapiro.test(TeguChemSoil$NO3ug)
hist(TeguChemSoil$NO3ug)
qqnorm(TeguChemSoil$NO3ug)
qqline(TeguChemSoil$NO3ug)
shapiro.test(ConChemSoil$NO3ug)
hist(ConChemSoil$NO3ug)
qqnorm(ConChemSoil$NO3ug)
qqline(ConChemSoil$NO3ug)

#do non-parametric comparison across all time points
kruskal.test(NO3ug ~ Location, data = SoilChemData) #Significant difference (1.477e-9)
#also do two sample t.test comparisons at time points with ggpubr package, Holm's correction default 
compare_means(NO3ug ~ Location, data = SoilChemData, method = "t.test", paired = FALSE, group.by = "Stage")

##pH PER SAMPLE 

pHdf <- summarySE(SoilChemData, measurevar="pH", groupvars=c("Location","Stage"))

svg("pHLine.svg", width = 9, height = 4)
ggplot(pHdf, aes(x=Stage, y=pH, color=Location, group=Location)) + 
  geom_errorbar(aes(ymin=pH-se, ymax=pH+se), width=.1) +
  coord_cartesian(ylim = c(6, 8)) +
  scale_color_manual(values = c("black","red")) +
  geom_line(size=1.2) +
  geom_point(aes(shape = Location), size=4) +
  scale_shape_manual(values = c(16, 15)) + 
  force_panelsizes(rows = unit(2, "in"), cols = unit(6, "in")) +
  labs(color = 'Treatment', shape='Treatment') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 20), legend.text = element_text(size=18), legend.title = element_text(size=20)) +
  ylab("pH") +
  xlab("Decomposition Stage")
dev.off()

#look at distribution of this variable with Shapiro Wilks test and visually 
shapiro.test(TeguChemSoil$pH)
hist(TeguChemSoil$pH)
qqnorm(TeguChemSoil$pH)
qqline(TeguChemSoil$pH)
shapiro.test(ConChemSoil$pH)
hist(ConChemSoil$pH)
qqnorm(ConChemSoil$pH)
qqline(ConChemSoil$pH)

#do non-parametric comparison across all time points
kruskal.test(pH ~ Location, data = SoilChemData) #Significant difference (0.0001725)
#also do two sample t.test comparisons at time points with ggpubr package, Holm's correction default 
compare_means(pH ~ Location, data = SoilChemData, method = "t.test", paired = FALSE, group.by = "Stage")

##EC PER SAMPLE 

ECdf <- summarySE(SoilChemData, measurevar="EC", groupvars=c("Location","Stage"))

svg("ECLine.svg", width = 9, height = 4)
ggplot(ECdf, aes(x=Stage, y=EC, color=Location, group=Location)) + 
  geom_errorbar(aes(ymin=EC-se, ymax=EC+se), width=.1) +
  scale_color_manual(values = c("black","red")) +
  geom_line(size=1.2) +
  geom_point(aes(shape = Location), size=4) +
  scale_shape_manual(values = c(16, 15)) + 
  force_panelsizes(rows = unit(2, "in"), cols = unit(6, "in")) +
  labs(color = 'Treatment', shape='Treatment') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=18), axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 20), legend.text = element_text(size=18), legend.title = element_text(size=20)) +
  labs (x="Decomposition Stage", y = expression(EC ~ (mu * S ~ cm^-1)))
dev.off()

#look at distribution of this variable with Shapiro Wilks test and visually 
shapiro.test(TeguChemSoil$EC)
hist(TeguChemSoil$EC)
qqnorm(TeguChemSoil$EC)
qqline(TeguChemSoil$EC)
shapiro.test(ConChemSoil$EC)
hist(ConChemSoil$EC)
qqnorm(ConChemSoil$EC)
qqline(ConChemSoil$EC)

#do non-parametric comparison across all time points
kruskal.test(EC ~ Location, data = SoilChemData) #Significant difference (3.879e-11)
#also do two sample t.test comparisons at time points with ggpubr package, Holm's correction default 

#END OF CHEMICAL ANALYSIS 

