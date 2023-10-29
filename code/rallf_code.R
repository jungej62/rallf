#Code to analyze rallf - root and low-lignin alfalfa
#Packages ----
library(nlme);library(multcomp);library(emmeans);library(multcompView)
library(lattice); library(reshape2)
library(Rmisc); library(tidyverse); 
library(stringr); library(gtsummary); library(readxl)
#Data readin ----
zz<-getwd()
ylds<-read_excel(paste0(zz,"/data/rallf_master.xlsx"), sheet="yld_soil", na="NA")
rics<-read_excel(paste0(zz,"/data/rallf_master.xlsx"), sheet="RIC", na="NA")
str(ylds)
#Data conversions ----
ylds$location<-as.factor(ylds$location)
ylds$fyr<-as.factor(ylds$year)
ylds$plot<-as.factor(ylds$plot)
ylds$rep<-as.factor(ylds$rep)
ylds$var<-as.factor(ylds$var)
ylds$harv<-as.factor(ylds$harv)
ylds$harv_trt<-as.factor(ylds$harv_trt)

wylds<-ylds %>%
  pivot_wider(names_from = harv, 
              values_from = c(dm_yld,RFV,RFQ,NDFD, NDF),
              id_cols=c(location, year, plot, rep, var, harv_trt))
wylds %>% 
  filter(year=="2022") %>% 
  rowwise(c(location, year, plot, rep, var, harv_trt)) %>% 
  summarise(sum_ylds = sum(c(dm_yld_1st, dm_yld_2nd, dm_yld_3rd, dm_yld_4th, dm_yld_5th), na.rm = TRUE)) %>% 
  group_by(location, var, harv_trt) %>% 
  summarise(m_yld = mean(sum_ylds, na.rm=T), 
            sd_yld = sd(sum_ylds, na.rm=T), 
            n_yld = n()) %>% 
  mutate(se_yld=sd_yld/sqrt(n_yld)) %>% 
  ggplot(aes(y=m_yld, x=var, fill=harv_trt))+
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
