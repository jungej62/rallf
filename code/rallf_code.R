#Code to analyze rallf - root and low-lignin alfalfa
#Packages ----
library(nlme);library(multcomp);library(emmeans);library(multcompView)
library(lattice); library(reshape2)
library(Rmisc); library(tidyverse); 
library(stringr); library(gtsummary); library(readxl)
#Data readin ----
zz<-getwd()
ylds<-read_excel(paste0(zz,"/data/rallf_master.xlsx"), sheet="yld_soil_2_clean", na="NA")
rics<-read_excel(paste0(zz,"/data/rallf_master.xlsx"), sheet="RIC", na="NA")
str(ylds)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Data conversions ----
ylds$location<-as.factor(ylds$location)
ylds$fyr<-as.factor(ylds$year)
ylds$plot<-as.factor(ylds$plot)
ylds$rep<-as.factor(ylds$rep)
ylds$var<-as.factor(ylds$var)
ylds$harv<-as.factor(ylds$harv)
ylds$harv_trt<-as.factor(ylds$harv_trt)
#Milk conversion ----
lab.average<-46.5
ylds$Fat<-1.96
ylds$NDFD48<-53.54+ylds$NDF*(-.178)
oldDMI<-(0.0115*1350)/0.3
ylds$EE<-ylds$Fat+1
ylds$DMI<-(120/ylds$NDF)+(ylds$NDFD48-45)*.374/1350*100
ylds$NFC<-100-((ylds$NDF*0.93)+ylds$CP+2.05+7.71)
ylds$NDFDa<-(45/lab.average)*ylds$NDFD48
ylds$NDFcp<-ylds$NDF*0.07
ylds$sTDN<-((0.93*ylds$CP)+(0.97*(ylds$Fat)*2.25)+(ylds$NFC*0.98)+((ylds$NDF-ylds$CP)*(ylds$NDFDa/100)))-7
ylds$TDN<-(ylds$NFC*.98)+(ylds$CP*.93)+(ylds$Fat*.97*2.25)+((ylds$NDFD48*.93)*(ylds$NDFD48/100))-7
ylds$RFQ<-ylds$DMI*ylds$TDN/1.23
ylds$RFQ2<-((((120/ylds$NDF)+(ylds$NDFD48-45)*0.374/1350*100)*((((100-((ylds$NDF*0.93)+ylds$CP+2.05+7.71))*0.98)+(ylds$CP*0.93)+(2.05*0.97*2.25)+((ylds$NDF*0.93)*(ylds$NDFD48/100)))-7))/1.23)
ylds$NE<-((ylds$TDN*0.0245)-0.12)/2.2
ylds$NEl<-(((((ylds$sTDN-((ylds$NDF-ylds$CP)*(ylds$NDFDa/100))+((((((lab.average+(0*12)-ylds$NDFDa)*0.374)*1.83)+ylds$NDFDa)/100)*(ylds$NDF-ylds$NDFcp)))*0.044)+0.207)*0.6741)-0.5656)/2.2
ylds$RFV<-(((88.9-(0.779*ylds$ADF))*(120/ylds$NDF))/1.29)
#base forage DMI
ylds$bDMI<-(0.0086*1350)/(ylds$NDF/100)
#adjusted forage DMI
ylds$aDMI<-((ylds$NDFD48-(45+(0*12)))*0.374)+ylds$bDMI
#adjusted total forage DMI
ylds$atDMI<-((ylds$NDFD48-(45+(0*12)))*0.374)+oldDMI
#Forage percental of total DMI
ylds$fDMI<-ylds$aDMI/ylds$atDMI
#milk yield from forage
ylds$fmilk<-((ylds$aDMI*ylds$NEl)-(0.08*(613.64^0.75)*ylds$fDMI))/0.31
#milk per ton of forage
ylds$milk_kg<-(ylds$fmilk/ylds$aDMI)*2.204#lbs milk per kg biomass
ylds$milk_kg2<-ylds$milk_kg/2.20462 #kg milk per kg biomass
ylds$myld2<-ylds$milk_kg2*ylds$dm_yld*1000 #kg milk per ha
#First reshape
wylds<-ylds %>%
  pivot_wider(names_from = harv, 
              values_from = c(dm_yld,RFV,RFQ,NDFD, NDF, myld2),
              id_cols=c(location, year, plot, rep, var, harv_trt))
wylds$fd<-as.factor(str_sub(wylds$var, 3, 4))
wylds$var2<-as.factor(str_sub(wylds$var, 1, 2))

#Summary of yields ----
wylds %>% 
  rowwise(c(location, year, plot,  rep,var, harv_trt)) %>% 
  summarise(sum_ylds = sum(c(dm_yld_1, dm_yld_2, dm_yld_3, dm_yld_4, dm_yld_5), na.rm = TRUE)) %>% 
  group_by(year, location, var, harv_trt) %>% 
  summarise(m_yld = mean(sum_ylds, na.rm=T), 
            sd_yld = sd(sum_ylds, na.rm=T), 
            n_yld = n()) %>% 
  mutate(se_yld=sd_yld/sqrt(n_yld)) %>% 
  ggplot(aes(y=m_yld, x=var, color=harv_trt))+
  geom_point(position=position_dodge(.9))+
  facet_grid(location~year)+
  geom_errorbar(aes(ymax=m_yld+se_yld, ymin=m_yld-se_yld), 
                width=0.5, position=position_dodge(.9))+
  ylab(expression("Forage yield " ~ (Mg ~ ha^{-1})))+
  coord_cartesian(ylim=c(0, 25))+
  xlab("Variety")+
  scale_color_manual(values=cbPalette)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.key.size =unit(0.75, "cm"),
        legend.text = element_text(size=12),
        axis.line = element_line(color='black'),
        legend.title=element_blank(),
        legend.position = c(.15, .9),
        axis.title.x=element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'))
ggsave("Yield.png", width=11, height=5, units="in", path="figures/")

#Summary of NDFD48 ----
ylds %>% 
  group_by(location, year, var, harv_trt) %>% 
  summarise(m_NDFD = mean(NDFD48, na.rm=T), 
            sd_NDFD = sd(NDFD48, na.rm=T), 
            n_NDFD = n()) %>% 
  mutate(se_NDFD=sd_NDFD/sqrt(n_NDFD)) %>% 
  ggplot(aes(y=m_NDFD, x=var, color=harv_trt))+
  geom_point(position=position_dodge(.9))+
  facet_grid(location~year)+
  geom_errorbar(aes(ymax=m_NDFD+se_NDFD, ymin=m_NDFD-se_NDFD), 
                width=0.5, position=position_dodge(.9))+
  ylab("NDFD (%)")+
  #coord_cartesian(ylim=c(0, 25))+
  xlab("Variety")+
  scale_color_manual(values=cbPalette)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.key.size =unit(0.75, "cm"),
        legend.text = element_text(size=12),
        axis.line = element_line(color='black'),
        legend.title=element_blank(),
        legend.position = c(.15, .9),
        axis.title.x=element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'))
ggsave("NDFD.png", width=11, height=5, units="in", path="figures/")

#Summary of NDF ----
ylds %>% 
  group_by(location, year, var, harv_trt) %>% 
  summarise(m_NDF = mean(NDF, na.rm=T), 
            sd_NDF = sd(NDF, na.rm=T), 
            n_NDF = n()) %>% 
  mutate(se_NDF=sd_NDF/sqrt(n_NDF)) %>% 
  ggplot(aes(y=m_NDF, x=var, color=harv_trt))+
  geom_point(position=position_dodge(.9))+
  facet_grid(location~year)+
  geom_errorbar(aes(ymax=m_NDF+se_NDF, ymin=m_NDF-se_NDF), 
                width=0.5, position=position_dodge(.9))+
  ylab("NDF (%)")+
  #coord_cartesian(ylim=c(0, 25))+
  xlab("Variety")+
  scale_color_manual(values=cbPalette)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.key.size =unit(0.75, "cm"),
        legend.text = element_text(size=12),
        axis.line = element_line(color='black'),
        legend.title=element_blank(),
        legend.position = c(.1, .6),
        axis.title.x=element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'))
ggsave("NDF.png", width=11, height=5, units="in", path="figures/")

#Crude protein ----
ylds %>% 
  group_by(location, year, var, harv_trt) %>% 
  summarise(m_CP = mean(CP, na.rm=T), 
            sd_CP = sd(CP, na.rm=T), 
            n_CP = n()) %>% 
  mutate(se_CP=sd_CP/sqrt(n_CP)) %>% 
  ggplot(aes(y=m_CP, x=var, color=harv_trt))+
  geom_point(position=position_dodge(.9))+
  facet_grid(location~year)+
  geom_errorbar(aes(ymax=m_CP+se_CP, ymin=m_CP-se_CP), 
                width=0.5, position=position_dodge(.9))+
  ylab("Crude protein (%)")+
  #coord_cartesian(ylim=c(0, 25))+
  xlab("Variety")+
  scale_color_manual(values=cbPalette)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.key.size =unit(0.75, "cm"),
        legend.text = element_text(size=12),
        axis.line = element_line(color='black'),
        legend.title=element_blank(),
        legend.position = c(.1, .6),
        axis.title.x=element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'))
ggsave("CrudeProtein.png", width=11, height=5, units="in", path="figures/")

#RFV ----
#There are issues with the 45 day data in 2024
ylds %>% 
  group_by(location, year, var, harv_trt) %>% 
  summarise(m_RFV = mean(RFV, na.rm=T), 
            sd_RFV = sd(RFV, na.rm=T), 
            n_RFV = n()) %>% 
  mutate(se_RFV=sd_RFV/sqrt(n_RFV)) %>% 
  ggplot(aes(y=m_RFV, x=var, color=harv_trt))+
  geom_point(position=position_dodge(.9))+
  facet_grid(location~year)+
  geom_errorbar(aes(ymax=m_RFV+se_RFV, ymin=m_RFV-se_RFV), 
                width=0.5, position=position_dodge(.9))+
  ylab("Relative feed value")+
  #coord_cartesian(ylim=c(0, 25))+
  xlab("Variety")+
  scale_color_manual(values=cbPalette)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.key.size =unit(0.75, "cm"),
        legend.text = element_text(size=12),
        axis.line = element_line(color='black'),
        legend.title=element_blank(),
        legend.position = c(.1, .6),
        axis.title.x=element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'))
ggsave("RFV.png", width=11, height=5, units="in", path="figures/")

#Milk yield HX in 2022 ----
wylds %>% 
  filter(year=="2022"&var2=="HX") %>% 
  rowwise(c(location, year, plot, rep, fd, harv_trt)) %>% 
  summarise(sum_ylds = sum(c(myld2_1, myld2_2, myld2_3, myld2_4, myld2_5), na.rm = TRUE)) %>% 
  group_by(location, fd, harv_trt) %>% 
  summarise(m_yld = mean(sum_ylds, na.rm=T), 
            sd_yld = sd(sum_ylds, na.rm=T), 
            n_yld = n()) %>% 
  mutate(se_yld=sd_yld/sqrt(n_yld)) %>% 
  filter(location=="Rosemount") %>% 
  ggplot(aes(y=m_yld, x=fd, fill=harv_trt))+
  geom_col(position="dodge")+
  geom_errorbar(aes(ymax=m_yld+se_yld, ymin=m_yld-se_yld), 
                width=0.5, position=position_dodge(.9))+
  ylab(expression("Milk yield " ~ (kg ~ ha^{-1})))+
  coord_cartesian(ylim=c(0, 28000))+
  xlab("Fall dormancy")+
  scale_fill_manual(values=cbPalette)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.key.size =unit(0.75, "cm"),
        legend.text = element_text(size=12),
        axis.line = element_line(color='black'),
        legend.title=element_blank(),
        legend.position = c(.88, .88),
        axis.title.x=element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'))
ggsave("Milk_2022_Ros.png", width=6, height=5, units="in", path="figures/")

#All Yields HX in 2022 ----
wylds %>% 
  filter(year=="2022"&var2=="HX") %>% 
  rowwise(c(location, year, plot, rep, fd, harv_trt)) %>% 
  summarise(sum_ylds = sum(c(dm_yld_1st, dm_yld_2nd, dm_yld_3rd, dm_yld_4th, dm_yld_5th), na.rm = TRUE)) %>% 
  group_by(location, fd, harv_trt) %>% 
  summarise(m_yld = mean(sum_ylds, na.rm=T), 
            sd_yld = sd(sum_ylds, na.rm=T), 
            n_yld = n()) %>% 
  mutate(se_yld=sd_yld/sqrt(n_yld)) %>% 
  filter(location=="St. Paul") %>% 
  ggplot(aes(y=m_yld, x=fd, fill=harv_trt))+
  geom_col(position="dodge")+
  geom_errorbar(aes(ymax=m_yld+se_yld, ymin=m_yld-se_yld), 
                width=0.5, position=position_dodge(.9))+
  ylab(expression("Forage yield " ~ (Mg ~ ha^{-1})))+
  coord_cartesian(ylim=c(0, 16))+
  xlab("Fall dormancy")+
  scale_fill_manual(values=cbPalette)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.key.size =unit(0.75, "cm"),
        legend.text = element_text(size=12),
        axis.line = element_line(color='black'),
        legend.title=element_blank(),
        legend.position = c(.88, .85),
        axis.title.x=element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'))
ggsave("Forage_2022_SP.png", width=6, height=5, units="in", path="figures/")

mod1<-wylds %>% 
  filter(year=="2022"&var2=="HX"&location=="St. Paul") %>% 
  rowwise(c(location, year, plot, rep, fd, harv_trt)) %>% 
  summarise(sum_ylds = sum(c(dm_yld_1st, dm_yld_2nd, dm_yld_3rd, dm_yld_4th, dm_yld_5th), na.rm = TRUE)) %>% 
  lme(sum_ylds~fd*harv_trt, random=~1|rep, data=., na.action=na.omit)
anova(mod1)
cld(emmeans(mod1, ~fd))
cld(emmeans(mod1, ~harv_trt))

#HX in 2021
wylds %>% 
  filter(year=="2022") %>% 
  rowwise(c(location, year, plot, rep, var2, fd, harv_trt)) %>% 
  summarise(sum_ylds = sum(c(dm_yld_1st, dm_yld_2nd, dm_yld_3rd, dm_yld_4th, dm_yld_5th), na.rm = TRUE)) %>% 
  group_by(location, var2, fd, harv_trt) %>% 
  summarise(m_yld = mean(sum_ylds, na.rm=T), 
            sd_yld = sd(sum_ylds, na.rm=T), 
            n_yld = n()) %>% 
  mutate(se_yld=sd_yld/sqrt(n_yld)) %>% 
  ggplot(aes(y=m_yld, x=var2, fill=harv_trt, alpha=fd))+
  facet_grid(~location)+
  geom_col(position="dodge")+
  #geom_bar(position=position_dodge(.9), stat="identity", width=.75) +
  geom_errorbar(aes(ymax=m_yld+se_yld, ymin=m_yld-se_yld), 
                width=0.5, position=position_dodge(.9))+
  ylab(expression("Forage yield " ~ (Mg ~ ha^{-1})))

#RFV plot - crazy ----
wylds %>% 
  filter(year=="2022") %>% 
  rowwise(c(location, year, plot, rep, var, harv_trt)) %>% 
  summarise(mean_rfv = mean(c(RFV_1st, RFV_2nd, RFV_3rd, RFV_4th, RFV_5th), na.rm = TRUE)) %>% 
  group_by(location, var, harv_trt) %>% 
  summarise(m_rfv= mean(mean_rfv, na.rm=T), 
            sd_rfv = sd(mean_rfv, na.rm=T), 
            n_rfv = n()) %>% 
  mutate(se_rfv=sd_rfv/sqrt(n_rfv)) %>% 
  ggplot(aes(y=m_rfv, x=var, fill=harv_trt))+
  facet_grid(~location)+
  geom_col(position="dodge")+
  #geom_bar(position=position_dodge(.9), stat="identity", width=.75) +
  geom_errorbar(aes(ymax=m_rfv+se_rfv, ymin=m_rfv-se_rfv), 
                width=0.5, position=position_dodge(.9))+
  ylab("Relative Feed Value")


#RFQ plot ----
wylds %>% 
  filter(year=="2022") %>% 
  rowwise(c(location, year, plot, rep, var, harv_trt)) %>% 
  summarise(mean_rfq = mean(c(RFQ_1st, RFQ_2nd, RFQ_3rd, RFQ_4th, RFQ_5th), na.rm = TRUE)) %>% 
  group_by(location, var, harv_trt) %>% 
  summarise(m_rfq= mean(mean_rfq, na.rm=T), 
            sd_rfq = sd(mean_rfq, na.rm=T), 
            n_rfq = n()) %>% 
  mutate(se_rfq=sd_rfq/sqrt(n_rfq)) %>% 
  ggplot(aes(y=m_rfq, x=var, fill=harv_trt))+
  facet_grid(~location)+
  geom_col(position="dodge")+
  #geom_bar(position=position_dodge(.9), stat="identity", width=.75) +
  geom_errorbar(aes(ymax=m_rfq+se_rfq, ymin=m_rfq-se_rfq), 
                width=0.5, position=position_dodge(.9))+
  ylab("Relative Feed Value")

#NDFD plot - values are super low ----
wylds %>% 
  filter(year=="2022") %>% 
  rowwise(c(location, year, plot, rep, var, harv_trt)) %>% 
  summarise(mean_ndfd = mean(c(NDFD_1st, NDFD_2nd, NDFD_3rd, NDFD_4th, NDFD_5th), na.rm = TRUE)) %>% 
  group_by(location, var, harv_trt) %>% 
  summarise(m_ndfd= mean(mean_ndfd, na.rm=T), 
            sd_ndfd = sd(mean_ndfd, na.rm=T), 
            n_ndfd = n()) %>% 
  mutate(se_ndfd=sd_ndfd/sqrt(n_ndfd)) %>% 
  ggplot(aes(y=m_ndfd, x=var, fill=harv_trt))+
  facet_grid(~location)+
  geom_col(position="dodge")+
  #geom_bar(position=position_dodge(.9), stat="identity", width=.75) +
  geom_errorbar(aes(ymax=m_ndfd+se_ndfd, ymin=m_ndfd-se_ndfd), 
                width=0.5, position=position_dodge(.9))+
  ylab("Relative Feed Value")

#NDF plot ----
wylds %>% 
  filter(year=="2022") %>% 
  rowwise(c(location, year, plot, rep, var, harv_trt)) %>% 
  summarise(mean_ndf = mean(c(NDF_1st, NDF_2nd, NDF_3rd, NDF_4th, NDF_5th), na.rm = TRUE)) %>% 
  group_by(location, var, harv_trt) %>% 
  summarise(m_ndf= mean(mean_ndf, na.rm=T), 
            sd_ndf = sd(mean_ndf, na.rm=T), 
            n_ndf = n()) %>% 
  mutate(se_ndf=sd_ndf/sqrt(n_ndf)) %>% 
  ggplot(aes(y=m_ndf, x=var, fill=harv_trt))+
  facet_grid(~location)+
  geom_col(position="dodge")+
  #geom_bar(position=position_dodge(.9), stat="identity", width=.75) +
  geom_errorbar(aes(ymax=m_ndf+se_ndf, ymin=m_ndf-se_ndf), 
                width=0.5, position=position_dodge(.9))+
  ylab("NDF (%)")

wylds %>% 
  filter(year=="2022") %>% 
  rowwise(c(location, year, plot, rep, var, harv_trt)) %>%
  summarise(mean_ndf = mean(c(NDF_1st, NDF_2nd, NDF_3rd, NDF_4th, NDF_5th), na.rm = TRUE)) %>% 
  filter(location=="Rosemount") %>% 
  lme(mean_ndf~var*harv_trt, random=~1|rep,
      na.action=na.omit, data=.) %>% 
  anova()

#Roots ----
str(rics)
rics$fyear<-as.factor(rics$year)
rics$plot<-as.factor(rics$plot)
rics$rep<-as.factor(rics$rep)
rics$var<-as.factor(rics$var)
rics$harv_trt<-as.factor(rics$harv_trt)
#GRP ----
rics %>% 
  group_by(year, location, var, harv_trt) %>% 
  summarise(m_grp= mean(GRP, na.rm=T), 
            sd_grp = sd(GRP, na.rm=T), 
            n_grp = n()) %>% 
  mutate(se_grp=sd_grp/sqrt(n_grp)) %>% 
  ggplot(aes(y=m_grp, x=var, color=harv_trt))+
  facet_grid(location~year)+
  scale_color_manual(values=cbPalette)+
  geom_point(position=position_dodge(.9))+
  geom_errorbar(aes(ymax=m_grp+se_grp, ymin=m_grp-se_grp),
                width=0.5, position=position_dodge(.9))+
  ylab("Gross Root Production")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.key.size =unit(0.75, "cm"),
        legend.text = element_text(size=12),
        axis.line = element_line(color='black'),
        legend.title=element_blank(),
        legend.position = c(.82, .88),
        axis.title.x=element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'))
ggsave("GRP.png", width=11, height=5, units="in", path="figures/")
#Root death ----
rics %>% 
  group_by(year, location, var, harv_trt) %>% 
  summarise(m_da= mean(DA, na.rm=T), 
            sd_da = sd(DA, na.rm=T), 
            n_da = n()) %>% 
  mutate(se_da=sd_da/sqrt(n_da)) %>% 
  ggplot(aes(y=m_da, x=var, color=harv_trt))+
  facet_grid(location~year)+
  geom_point(position=position_dodge(.9))+
  scale_color_manual(values=cbPalette)+
  geom_errorbar(aes(ymax=m_da+se_da, ymin=m_da-se_da),
                width=0.5, position=position_dodge(.9))+
  ylab("Root mortality (g/core/year)")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.key.size =unit(0.75, "cm"),
        legend.text = element_text(size=12),
        axis.line = element_line(color='black'),
        legend.title=element_blank(),
        legend.position = c(.82, .6),
        axis.title.x=element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'))
ggsave("Mortality.png", width=11, height=5, units="in", path="figures/")

#Root ave_rate ----
rics %>% 
  filter(fyear!="2021") %>% 
  group_by(year, location, var, harv_trt) %>% 
  summarise(m_ave_rate= mean(ave_rate, na.rm=T), 
            sd_ave_rate = sd(ave_rate, na.rm=T), 
            n_ave_rate = n()) %>% 
  mutate(se_ave_rate=sd_ave_rate/sqrt(n_ave_rate)) %>% 
  ggplot(aes(y=m_ave_rate, x=var, color=harv_trt))+
  facet_grid(location~year)+
  geom_point(position=position_dodge(.9))+
  geom_errorbar(aes(ymax=m_ave_rate+se_ave_rate, ymin=m_ave_rate-se_ave_rate),
                width=0.5, position=position_dodge(.9))+
  ylab("Daily root growth rate")

#Root ave_rate by year ----
rics %>% 
  group_by(location, fyear, var, harv_trt) %>% 
  summarise(m_ave_rate= mean(ave_rate, na.rm=T), 
            sd_ave_rate = sd(ave_rate, na.rm=T), 
            n_ave_rate = n()) %>% 
  mutate(se_ave_rate=sd_ave_rate/sqrt(n_ave_rate)) %>% 
  filter(fyear!="2021") %>% 
  ggplot(aes(y=m_ave_rate, x=var, color=harv_trt))+
  facet_grid(location~fyear)+
  geom_point(position=position_dodge(.9))+
  geom_errorbar(aes(ymax=m_ave_rate+se_ave_rate, ymin=m_ave_rate-se_ave_rate),
                width=0.5, position=position_dodge(.9))+
  scale_color_manual(values=cbPalette)+
  ylab("Daily root growth rate (g/day)")+
  xlab("Variety")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.key.size =unit(0.75, "cm"),
        legend.text = element_text(size=12),
        axis.line = element_line(color='black'),
        legend.title=element_blank(),
        legend.position = c(.82, .88),
        axis.title.x=element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'))
ggsave("RootGrowthRate.png", width=8, height=5, units="in", path="figures/")

#####Old analyses ----
#Material below this has not been edited to account for complete dataset.
#LTRI4
rics %>% 
  group_by(location, fyear, var, harv_trt) %>% 
  summarise(m_ave_rate= mean(LTRI3_mass, na.rm=T), 
            sd_ave_rate = sd(LTRI3_mass, na.rm=T), 
            n_ave_rate = n()) %>% 
  mutate(se_ave_rate=sd_ave_rate/sqrt(n_ave_rate)) %>% 
  filter(fyear=="2022"&location=="Rosemount") %>% 
  ggplot(aes(y=m_ave_rate, x=var, fill=harv_trt))+
  facet_grid(~location)+
  geom_col(position="dodge")+
  geom_errorbar(aes(ymax=m_ave_rate+se_ave_rate, ymin=m_ave_rate-se_ave_rate),
                width=0.5, position=position_dodge(.9))+
  scale_fill_manual(values=cbPalette)+
  ylab("Net root growth rate")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.key.size =unit(0.75, "cm"),
        legend.text = element_text(size=12),
        axis.line = element_line(color='black'),
        legend.title=element_blank(),
        legend.position = c(.88, .88),
        axis.title.x=element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'))
ggsave("LTRI3_2022_Ros.png", width=6, height=5, units="in", path="figures/")

#SS1 mass
rics %>% 
  group_by(location, fyear, harv_trt) %>% 
  summarise(m_SS1_mass= mean(SS1_mass, na.rm=T), 
            sd_SS1_mass = sd(SS1_mass, na.rm=T), 
            n_SS1_mass = n()) %>% 
  mutate(se_SS1_mass=sd_SS1_mass/sqrt(n_SS1_mass)) %>% 
  ggplot(aes(y=m_SS1_mass, x=fyear, fill=harv_trt))+
  facet_grid(~location)+
  geom_col(position="dodge")+
  geom_errorbar(aes(ymax=m_SS1_mass+se_SS1_mass, ymin=m_SS1_mass-se_SS1_mass),
                width=0.5, position=position_dodge(.9))+
  ylab("Standing Stock root mass")

#Standing Stock by day
ss_ric<-rics %>% 
  filter(fyear=="2022") %>% 
  select(c(location, harv_trt, var, SS1_mass, SS2_mass, SS3_mass, SS4_mass, SS5_mass)) %>% 
  pivot_longer(c(SS1_mass, SS2_mass, SS3_mass, SS4_mass, SS5_mass), names_to="SS_time", values_to="root_biomass")

ss_ric2<- rics %>% 
  filter(fyear=="2022") %>% 
  select(c(location, harv_trt, var, SS1_extract,
           SS2_extract,  SS3_extract, 
           SS4_extract,  SS5_extract)) %>% 
  pivot_longer(c(SS1_extract, SS2_extract, SS3_extract, SS4_extract, SS5_extract), names_to="SS_time", values_to="extract_date")

ss_ric$extract_date<-ss_ric2$extract_date

ss_ric %>% 
  group_by(location, harv_trt, extract_date) %>% 
  summarise(m_root_biomass= mean(root_biomass, na.rm=T), 
            sd_root_biomass = sd(root_biomass, na.rm=T), 
            n_root_biomass = n()) %>% 
  mutate(se_root_biomass=sd_root_biomass/sqrt(n_root_biomass)) %>%
ggplot(aes(x=extract_date, y=m_root_biomass, color=harv_trt))+
  geom_point()+
  geom_errorbar(aes(ymax=m_root_biomass+se_root_biomass, ymin=m_root_biomass-se_root_biomass))+
  facet_grid(~location)

#RootCN ----
#Need to finish
rics$STRI1_cn<-rics$STRI1_c/rics$STRI1_n
rics$STRI2_cn<-rics$STRI2_c/rics$STRI2_n
rics$STRI3_cn<-rics$STRI3_c/rics$STRI3_n
rics %>% 
  filter(year=="2021") %>% 
  rowwise(c(location, plot, rep, var, harv_trt)) %>% 
  summarise(mean_cn = mean(c(STRI1_cn, STRI2_cn, STRI3_cn), na.rm = TRUE)) %>% 
  group_by(location, var, harv_trt) %>% 
  summarise(m_cn = mean(mean_cn, na.rm=T), 
            sd_cn = sd(mean_cn, na.rm=T), 
            n_cn = n()) %>% 
  mutate(se_cn=sd_cn/sqrt(n_cn)) %>% 
  filter(location=="Rosemount") %>% 
  ggplot(aes(y=m_cn, x=var, fill=harv_trt))+
  facet_grid(~location)+
  geom_col(position="dodge")+
  geom_errorbar(aes(ymax=m_cn+se_cn, ymin=m_cn-se_cn),
                width=0.5, position=position_dodge(.9))+
  scale_fill_manual(values=cbPalette)+
  ylab("Mean C:N")+
  coord_cartesian(y=c(0,18))+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.key.size =unit(0.75, "cm"),
        legend.text = element_text(size=12),
        axis.line = element_line(color='black'),
        legend.title=element_blank(),
        legend.position = c(.2, .88),
        axis.title.x=element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'))
ggsave("CN_2021_Ros.png", width=6, height=5, units="in", path="figures/")
#SP
rics %>% 
  filter(year=="2021") %>% 
  rowwise(c(location, plot, rep, var, harv_trt)) %>% 
  summarise(mean_cn = mean(c(STRI1_cn, STRI2_cn, STRI3_cn), na.rm = TRUE)) %>% 
  group_by(location, var, harv_trt) %>% 
  summarise(m_cn = mean(mean_cn, na.rm=T), 
            sd_cn = sd(mean_cn, na.rm=T), 
            n_cn = n()) %>% 
  mutate(se_cn=sd_cn/sqrt(n_cn)) %>% 
  filter(location=="St. Paul") %>% 
  ggplot(aes(y=m_cn, x=var, fill=harv_trt))+
  facet_grid(~location)+
  geom_col(position="dodge")+
  geom_errorbar(aes(ymax=m_cn+se_cn, ymin=m_cn-se_cn),
                width=0.5, position=position_dodge(.9))+
  scale_fill_manual(values=cbPalette)+
  ylab("Mean C:N")+
  coord_cartesian(y=c(0,18))+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.key.size =unit(0.75, "cm"),
        legend.text = element_text(size=12),
        axis.line = element_line(color='black'),
        legend.title=element_blank(),
        legend.position = c(.88, .88),
        axis.title.x=element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.title.y = element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'))
ggsave("CN_2021_SP.png", width=6, height=5, units="in", path="figures/")
