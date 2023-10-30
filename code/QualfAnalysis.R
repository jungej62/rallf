dat<-read.csv("/Users/junge037/Documents/Projects/Misc/QUALF/QUALF.csv", na.strings="NA")
source.with.encoding('~/Documents/Projects/HelperFunctions.R', encoding='UTF-8') #You don't need this
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","gray30")

library(Rmisc) #need this for summarySE function.
library(ggplot2); library(lattice); library(reshape2)
library(gridExtra)

head(dat)
dat$fRep<-factor(dat$Rep)
dat$fEntry<-factor(dat$Entry)
dat$fCut<-factor(dat$Cut)

sumprotein<-summarySE(dat, "Protein", c("Loc", "fEntry", "Cut"))

ggplot(sumprotein, aes(x=Cut, y=Protein, color=fEntry))+
  facet_grid(~Loc)+
  geom_point(stat="identity") +
  geom_line()+
  ggtitle("Quality parameter changes by cut")
ggsave("ProteinTime.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")

proteinreg <- dlply(dat, .(Loc, fEntry), lm, formula = Protein ~ Cut )
proteinslope <- cbind(ldply(proteinreg, coef)[,c(1,2,4)],
                  ldply(proteinreg, function(x){coef(summary(x))[4]})[,3])
colnames(proteinslope)<-c("Loc", "fEntry", "Slope", "Slope_se")
ggplot(proteinslope, aes(x=fEntry, y=Slope))+
  facet_grid(~Loc)+
  geom_point()+
  geom_errorbar(aes(ymin=Slope-Slope_se, ymax=Slope+Slope_se))+
ggtitle("Average slopes and standard error")
ggsave("ProteinSlopes.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")

ggplot(summarySE(dat, "Protein", c("Loc","fEntry")),
       aes(x=fEntry, y=Protein))+
  facet_grid(~Loc)+
  geom_point()+
  geom_errorbar(aes(ymin=Protein-se, ymax=Protein+se))+
  ggtitle("Average over cutting dates")
ggsave("ProteinAverage.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")

### ADF
sumADF<-summarySE(dat, "ADF", c("Loc", "fEntry", "Cut"))

ggplot(sumADF, aes(x=Cut, y=ADF, color=fEntry))+
  facet_grid(~Loc)+
  geom_point(stat="identity") +
  geom_line()+
  ggtitle("Quality parameter changes by cut")
ggsave("ADFTime.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")

ADFreg <- dlply(dat, .(Loc, fEntry), lm, formula = ADF ~ Cut )
ADFslope <- cbind(ldply(ADFreg, coef)[,c(1,2,4)],
                      ldply(ADFreg, function(x){coef(summary(x))[4]})[,3])
colnames(ADFslope)<-c("Loc", "fEntry", "Slope", "Slope_se")
ggplot(ADFslope, aes(x=fEntry, y=Slope))+
  facet_grid(~Loc)+
  geom_point()+
  geom_errorbar(aes(ymin=Slope-Slope_se, ymax=Slope+Slope_se))+
  ggtitle("Average slopes and standard error")
ggsave("ADFSlopes.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")

ggplot(summarySE(dat, "ADF", c("Loc","fEntry")),
       aes(x=fEntry, y=ADF))+
  facet_grid(~Loc)+
  geom_point()+
  geom_errorbar(aes(ymin=ADF-se, ymax=ADF+se))+
  ggtitle("Average over cutting dates")
ggsave("ADFAverage.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")

### NDF
sumNDF<-summarySE(dat, "NDF", c("Loc", "fEntry", "Cut"))

ggplot(sumNDF, aes(x=Cut, y=NDF, color=fEntry))+
  facet_grid(~Loc)+
  geom_point(stat="identity") +
  geom_line()+
  ggtitle("Quality parameter changes by cut")
ggsave("NDFTime.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")

NDFreg <- dlply(dat, .(Loc, fEntry), lm, formula = NDF ~ Cut )
NDFslope <- cbind(ldply(NDFreg, coef)[,c(1,2,4)],
                  ldply(NDFreg, function(x){coef(summary(x))[4]})[,3])
colnames(NDFslope)<-c("Loc", "fEntry", "Slope", "Slope_se")
ggplot(NDFslope, aes(x=fEntry, y=Slope))+
  facet_grid(~Loc)+
  geom_point()+
  geom_errorbar(aes(ymin=Slope-Slope_se, ymax=Slope+Slope_se))+
  ggtitle("Average slopes and standard error")
ggsave("NDFSlopes.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")

ggplot(summarySE(dat, "NDF", c("Loc","fEntry")),
       aes(x=fEntry, y=NDF))+
  facet_grid(~Loc)+
  geom_point()+
  geom_errorbar(aes(ymin=NDF-se, ymax=NDF+se))+
  ggtitle("Average over cutting dates")
ggsave("NDFAverage.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")

### NDFD48
sumNDFD48<-summarySE(dat, "NDFD48", c("Loc", "fEntry", "Cut"))

ggplot(sumNDFD48, aes(x=Cut, y=NDFD48, color=fEntry))+
  facet_grid(~Loc)+
  geom_point(stat="identity") +
  geom_line()+
  ggtitle("Quality parameter changes by cut")
ggsave("NDFD48Time.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")

NDFD48reg <- dlply(dat, .(Loc, fEntry), lm, formula = NDFD48 ~ Cut )
NDFD48slope <- cbind(ldply(NDFD48reg, coef)[,c(1,2,4)],
                  ldply(NDFD48reg, function(x){coef(summary(x))[4]})[,3])
colnames(NDFD48slope)<-c("Loc", "fEntry", "Slope", "Slope_se")
ggplot(NDFD48slope, aes(x=fEntry, y=Slope))+
  facet_grid(~Loc)+
  geom_point()+
  geom_errorbar(aes(ymin=Slope-Slope_se, ymax=Slope+Slope_se))+
  ggtitle("Average slopes and standard error")
ggsave("NDFD48Slopes.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")

ggplot(summarySE(dat, "NDFD48", c("Loc","fEntry")),
       aes(x=fEntry, y=NDFD48))+
  facet_grid(~Loc)+
  geom_point()+
  geom_errorbar(aes(ymin=NDFD48-se, ymax=NDFD48+se))+
  ggtitle("Average over cutting dates")
ggsave("NDFD48Average.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")


#Lignin
sumLignin<-summarySE(dat, "Lignin", c("Loc", "fEntry", "Cut"))

ggplot(sumLignin, aes(x=Cut, y=Lignin, color=fEntry))+
  facet_grid(~Loc)+
  geom_point(stat="identity") +
  geom_line()+
  ggtitle("Quality parameter changes by cut")
ggsave("LigninTime.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")

Ligninreg <- dlply(dat, .(Loc, fEntry), lm, formula = Lignin ~ Cut )
Ligninslope <- cbind(ldply(Ligninreg, coef)[,c(1,2,4)],
                     ldply(Ligninreg, function(x){coef(summary(x))[4]})[,3])
colnames(Ligninslope)<-c("Loc", "fEntry", "Slope", "Slope_se")
ggplot(Ligninslope, aes(x=fEntry, y=Slope))+
  facet_grid(~Loc)+
  geom_point()+
  geom_errorbar(aes(ymin=Slope-Slope_se, ymax=Slope+Slope_se))+
  ggtitle("Average slopes and standard error")
ggsave("LigninSlopes.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")

ggplot(summarySE(dat, "Lignin", c("Loc","fEntry")),
       aes(x=fEntry, y=Lignin))+
  facet_grid(~Loc)+
  geom_point()+
  geom_errorbar(aes(ymin=Lignin-se, ymax=Lignin+se))+
  ggtitle("Average over cutting dates")
ggsave("LigninAverage.pdf", width=10, height=5, units="in", path="/Users/junge037/Documents/Projects/Misc/QUALF/")