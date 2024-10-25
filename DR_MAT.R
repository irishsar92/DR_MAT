

# Setup -------------------------------------------------------------------


suppressPackageStartupMessages({
  library(glmmTMB)
  library(car)
  library(popbio)
  library(ggplot2)
  library(ggbeeswarm)
  library(shades)
  library(RColorBrewer)
  library(DHARMa)
  library(Rmisc)
  library(dabestr)
  library(magrittr)
  library(tidyr)
  library(dplyr)
  library(ggeffects)
  library(lme4)
  library(optimx)
  library(nloptr)
  library(dfoptim)
  library(data.table)
  library(DHARMa)
  library(emmeans)
  library(ggpubr)
  library(patchwork)
  library(lubridate)
  library(performance)
  library(dplyr)#to check overdispersion in glmmTMB models
  library(survival)
  library(coxme) #for cox model
  library(purrr) #for forestplot
  library(ggpubr) #for forestplot
  library(survminer) #for forestplot
  library(AICcmodavg)
  library(forcats)
  library(fields) #for thin plate spline plots
  library(viridis)
  library(scales)
})

palette <- c('#1F77B4FF', '#FF7F0EFF', '#2CA02CFF','#D62728FF')


#edit_git_config()
#use_git()
#create_github_token()
#gitcreds_set()
#use_github()

R.version



# DRO Lifespan ------------------------------------------------------------

DRO_LS <- read.csv('DRO_LS.csv')
DRO_LS <- na.omit(DRO_LS)
DRO_LS$Treatment.ID <- as.factor(DRO_LS$Treatment.ID)
DRO_LS$Treatment.ID <- droplevels(DRO_LS$Treatment.ID)
levels(DRO_LS$Treatment.ID)

DRO_LS2 <- DRO_LS %>%
  filter(DRO_LS$Cause != 'L' & DRO_LS$Cause !='W' & DRO_LS$Cause != 'E')

surv<-survfit(Surv(Age,Event)~Treatment.ID,data=DRO_LS2)


LS_plot <-ggsurvplot(surv, ylab="Survival probability\n", data = DRO_LS2, size= 0.8, font.ylab= 18, font.xlab= 18, legend = c(0.3, 0.4), legend.title = "", title = "Matricides uncensored", censor = FALSE, xlab = "\nDay", xlim=c(0,32), break.time.by = 5, palette = palette, position= position_dodge(0.9), font.tickslab = c(14), legend.labs=c("DR" ,"DR+O", "F","F+O"), font.legend = c(16))
LS_plot

DRO_LS3 <- DRO_LS2 %>%
  filter(DRO_LS2$Cause != 'M')

surv<-survfit(Surv(Age,Event)~Treatment.ID,data=DRO_LS3)

LS_plot_matcen <-ggsurvplot(surv, ylab="Survival probability\n", data = DRO_LS3, size= 0.8, font.ylab= 18, font.xlab= 18, legend = c(0.3, 0.4), legend.title = "", palette = palette, title = "Matricides censored", censor = FALSE, xlab = "\nDay", xlim=c(0,32), break.time.by = 5, position= position_dodge(0.9), legend.labs=c("DR" ,"DR+O", "F","F+O"), font.tickslab = c(14), font.legend = c(16))
LS_plot_matcen

DRO_LS3$Treatment <- as.factor(DRO_LS3$Treatment)
DRO_LS3$Treatment <- relevel(DRO_LS3$Treatment, ref = 'F')
DRO_LS2$Treatment <- relevel(DRO_LS2$Treatment, ref = 'F')

#Function for creating forest plot
meforest <- function(cox, F){  #Eds version of forest plot
  require(AICcmodavg)
  require(ggplot2)
  require(forcats)
  store <- matrix(nrow = length(cox$coefficients) + 1, ncol = 4)
  ref <- c(paste(F), 0,NA, NA)
  store[1,] <- ref
  for (x in 1:length(cox$coefficients)){
    y = x+1
    mean<-cox$coefficients[x]
    emean <- (mean)
    CIL <- cox$coefficients[x] - (1.96*extractSE(cox)[x])
    CUL <- cox$coefficients[x] + (1.96*extractSE(cox)[x])
    eCIL <- (CIL)
    eCUL <- (CUL)
    store[y, 1:4] <- c(names(cox$coefficients[x]), emean, eCIL, eCUL)}
  colnames(store) <- c("Treatment", "mean", "CIL","CUL")
  store <- as.data.frame(store)
  store$mean <- as.numeric(as.character(store$mean))
  store$CIL <- as.numeric(as.character(store$CIL))
  store$CUL <- as.numeric(as.character(store$CUL))
  forest<-ggplot(store, aes(x = mean, xmax = CUL, xmin = CIL, y = Treatment, colour = Treatment)) +
    geom_point(size=4, shape = 19, colour = c('#1F77B4FF', '#FF7F0EFF', '#2CA02CFF','#D62728FF')) +
    geom_errorbarh(height=0.25, size=0.9,  colour = c('#1F77B4FF', '#FF7F0EFF', '#2CA02CFF','#D62728FF')) +
    theme_minimal() +
    geom_vline(xintercept=0, linetype="dotted", size=0.8) +
    xlab("Hazard Ratio")+
    ylab("")+
    theme(
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 20)
    )+
    scale_y_discrete(limits = (levels(store$Treatment)), labels = c('TreatmentDR'='DR', 'TreatmentDRO' = 'DR+O', 'TreatmentF' ='F', 'TreatmentFO' ='F+O'))+
    expand_limits(x = c(-1.5,0.5))
  return(forest)
}


cox <- coxme(Surv(Age, Event) ~ Treatment +(1|Plate.ID), data = DRO_LS2)
summary(cox)

forest <-meforest(cox, "F")
forest

cox_matcen <- coxme(Surv(Age, Event) ~ Treatment +(1|Plate.ID), data = DRO_LS3)
summary(cox_matcen)

forest_matcen <- meforest(cox_matcen, 'F')
forest_matcen

DR_all_plot <- ggarrange(LS_plot$plot, LS_plot_matcen$plot, forest, forest_matcen, nrow = 2, ncol = 2,heights = c(2, 1))
DR_all_plot
ggsave('DR_all_plot.tif', height = 8, width = 12)

#test cox models for fit
cox.zph(cox)#looks good - use cox model results
cox.zph(cox_matcen)#looks good









# DRO Reproduction --------------------------------------------------------


DR_rep <- read.csv('DR_Lab_repro_2.csv')

rep_long <- DR_rep %>% 
  pivot_longer(
    cols = `D1`:`D8`, 
    names_to = "Day",
    values_to = "value"
  )

rep_long$Treatment <- factor(rep_long$Treatment, levels = c('F','FO','DR','DRO'))

rep_long$Treatment <- revalue(rep_long$Treatment, c('F' = 'F', 'FO' = 'F+O', 'DR' = 'DR', 'DRO' = 'DR+O'))

## Age-specific reproduction plot
rep <-ggplot(data=rep_long, aes(x=factor(Day), y=value, group=Treatment, color=Treatment))+
  geom_jitter(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5))+
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", size = 1, width=0.0, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="line",  size=1, position = position_dodge(0.5)) +
  theme_classic()+
  labs(y="Offspring number", x="", size=20)+
  labs(col="")+
  theme(axis.title.y = element_text(size=12))+
  theme(axis.title.x = element_text(size=12))+
  theme(legend.key.width = unit(0.5,"cm"))+
  coord_cartesian(ylim = c(0,170))+
  theme(legend.position = c(0.8,0.9))+
  scale_colour_manual(values = palette)

rep


#Sum number of offspring over 5 days for individuals
totalrep<-na.omit(as.data.frame.table(tapply(rep_long$value,list(rep_long$Treatment, rep_long$ID),sum)))


#rename columns
names(totalrep)<-c("Treatment", "Replicate", "Totrep")

totrep_dab <-
  totalrep %>%
  load(x =Treatment, y=Totrep,
       idx = c('F','F+O','DR','DR+O'))

totrep_p <- mean_diff(totrep_dab)

tot_rep_plot <- dabest_plot(totrep_p, FALSE, raw_marker_spread = 1, swarm_label = 'LRS',custom_palette = 'd3')

tot_rep_plot

hist(totalrep$Totrep)
tot_rep_mod <- lm(Totrep ~ Treatment, data = totalrep)
summary(tot_rep_mod)

#### Lambda

#spread data back out
rep_wide <- spread(rep_long, Day, value)

rep_wide <- na.omit(rep_wide)

#calculate Lambda for each individual
L <- matrix(nrow = nrow(rep_wide), ncol = 4)

for (i in 1:nrow(rep_wide)){
  Les <- matrix(0, ncol = 10, nrow = 10) #ncol and nrow should = no. days repro +2
  diag(Les[-1,]) <- rep(1, 9) # add the 1s for survival probability diagonally
  Fert <- c(0,0, as.numeric(as.vector(rep_wide[i,][4:11]))) #the columns in data that has the reproductive counts
  Fert[is.na(Fert)] <- 0 #makes all NAs into 0s
  Les[1,] <- c(Fert)
  class(Les) <- "leslie.matrix"
  Lambda <- popbio::eigen.analysis(Les)$lambda1
  L[i, 1:4] <- c(paste0(rep_wide$Block[i]),paste0(rep_wide$ID[i]), paste0(rep_wide$Treatment[i]), Lambda)
  
}

#rename columns
colnames(L)<-c("Block", "ID", "Treatment","Lambda")

#make L a data frame
Data<-as.data.frame(L)

#make sure it's numeric
Data$Lambda<-as.numeric(as.character(Data$Lambda))

str(Data)

#create dabestr plot
Lambda_rep <-
  Data %>%
  load(Treatment, Lambda,
       idx = list(c("F", "F+O", "DR", "DR+O")))

Lambda_rep_dab <- mean_diff(Lambda_rep)

lam_plot <- dabest_plot(Lambda_rep_dab, FALSE, swarm_label = 'Lambda', raw_marker_spread = 1, custom_palette = 'd3')

lam_plot

rep/(tot_rep_plot+lam_plot)

ggsave('DRO_reproduction.tif', height = 8, width = 10)

hist(Data$Lambda)
lam_mod <- lm(Lambda ~ Treatment, data = Data)
summary(lam_mod)
sim_lam <- simulateResiduals(lam_mod, plot = T)


## Analysis of age-specific reproduction

hist(rep_long$value) #use negbin model

rep_long$Day <- revalue(rep_long$Day, c('D1' = '1', 'D2' = '2', 'D3' = '3', 'D4' = '4', 'D5' = '5', 'D6' = '6', 'D7' = '7', 'D8' = '8'))
rep_long$Day <- as.numeric(rep_long$Day)


rep_mod <- glmmTMB(value ~ Treatment * Day + (1|ID), data = rep_long, family = 'nbinom2')
summary(rep_mod)
sim1 <- simulateResiduals(rep_mod, plot = T)
testDispersion(sim1)
check_overdispersion(rep_mod)
testZeroInflation(rep_mod)

rep_mod2 <- glmmTMB(value ~ Treatment * Day + (1|ID), data = rep_long, ziformula = ~Day, family = 'nbinom2')
summary(rep_mod2)
sim2 <- simulateResiduals(rep_mod2, plot = T)
AIC(rep_mod, rep_mod2)                  

rep_mod3 <- glmmTMB(value ~ Treatment * Day + (1|ID), data = rep_long, ziformula = ~Day, disp = ~Day, family = 'nbinom2')
summary(rep_mod3)
sim3 <- simulateResiduals(rep_mod3)
AIC(rep_mod, rep_mod2, rep_mod3)




# DRO Mated Reproduction --------------------------------------------------

male_rep <- read.csv('Male_repro_DR.csv')

male_rep <- male_rep %>%
  filter(Lost != 'INF' & Lost != 'L')


male_rep$Treatment <- revalue(male_rep$Treatment, c('F' = 'F', 'FO' = 'F+O', 'DR' = 'DR', 'DRO' = 'DR+O'))

male_rep <- subset(male_rep, select = -c(D11, D12))

rep_long <- male_rep %>% 
  pivot_longer(
    cols = `D1`:`D10`, 
    names_to = "Day",
    values_to = "value"
  )

rep_long$Treatment <- factor(rep_long$Treatment, levels = c('F','F+O','DR','DR+O'))

levels(rep_long$Treatment)




rep_long$Day <- factor(rep_long$Day, levels = c('D1','D2','D3','D4','D5','D6','D7','D8','D9','D10'))

## Age-specific reproduction plot

rep <-ggplot(data=rep_long, aes(x=factor(Day), y=value, group=Treatment, color=Treatment))+
  geom_jitter(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5))+
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", size = 1, width=0.0, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="line",  size=1, position = position_dodge(0.5)) +
  theme_classic()+
  labs(y="Offspring number", x="", size=20)+
  labs(col="")+
  theme(axis.title.y = element_text(size=12))+
  theme(axis.title.x = element_text(size=12))+
  theme(legend.key.width = unit(0.5,"cm"))+
  coord_cartesian(ylim = c(0,170))+
  theme(legend.position = c(0.8,0.9))+
  scale_colour_manual(values = palette)

rep

#Sum number of offspring over 5 days for individuals
totalrep<-na.omit(as.data.frame.table(tapply(rep_long$value,list(rep_long$Treatment, rep_long$ID),sum)))


#rename columns
names(totalrep)<-c("Treatment", "Replicate", "Totrep")

totrep_dab <-
  totalrep %>%
  load(x =Treatment, y=Totrep,
       idx = c('F','F+O','DR','DR+O'))



totrep_p <- mean_diff(totrep_dab)

tot_rep_plot <- dabest_plot(totrep_p, FALSE, raw_marker_spread = 1, swarm_label = 'LRS',custom_palette = 'd3')

tot_rep_plot

Mated_tot <- lm(Totrep ~ Treatment, data = totalrep)
summary(Mated_tot)
mated_tot_sim <- simulateResiduals(Mated_tot, plot = T)


#### Lambda

#spread data back out
rep_wide <- spread(rep_long, Day, value)

rep_wide <- na.omit(rep_wide)

rep_wide <- subset(rep_wide, select = -c(Lost))

#calculate Lambda for each individual
L <- matrix(nrow = nrow(rep_wide), ncol = 4)

for (i in 1:nrow(rep_wide)){
  Les <- matrix(0, ncol = 12, nrow = 12) #ncol and nrow should = no. days repro +2
  diag(Les[-1,]) <- rep(1, 11) # add the 1s for survival probability diagonally
  Fert <- c(0,0, as.numeric(as.vector(rep_wide[i,][4:13]))) #the columns in data that has the reproductive counts
  Fert[is.na(Fert)] <- 0 #makes all NAs into 0s
  Les[1,] <- c(Fert)
  class(Les) <- "leslie.matrix"
  Lambda <- popbio::eigen.analysis(Les)$lambda1
  L[i, 1:4] <- c(paste0(rep_wide$Block[i]), paste0(rep_wide$Treatment[i]), paste0(rep_wide$ID[i]), Lambda)
  
}

#rename columns
colnames(L)<-c("Block",  "Treatment","ID","Lambda")

#make L a data frame
Data<-as.data.frame(L)

#make sure it's numeric
Data$Lambda<-as.numeric(as.character(Data$Lambda))

str(Data)


mated_lam <- lm(Lambda ~ Treatment, data = Data)
summary(mated_lam)
sim_lam_mated <- simulateResiduals(mated_lam, plot = T)

#create dabestr plot
Lambda_rep <-
  Data %>%
  load(Treatment, Lambda,
       idx = list(c("F", "F+O", "DR", "DR+O")))

Lambda_rep_dab <- mean_diff(Lambda_rep)

lam_plot <- dabest_plot(Lambda_rep_dab, FALSE, swarm_label = 'Lambda', raw_marker_spread = 1, custom_palette = 'd3')

lam_plot

rep/(tot_rep_plot + lam_plot)

ggsave('DRO_male_reproduction.tif', height = 8, width = 10)


#Analyse age-specific reproduction 

rep_long$Day <- revalue(rep_long$Day, c('D1' = '1', 'D2' = '2', 'D3' = '3', 'D4' = '4', 'D5' = '5', 'D6' = '6', 'D7' = '7', 'D8' = '8', 'D9' = '9', 'D10' = '10'))
rep_long$Day <- as.numeric(rep_long$Day)

hist(rep_long$value)

mated_mod <- glmmTMB(value ~ Treatment*Day + (1|ID), data = rep_long, family = 'nbinom2')
summary(mated_mod)
sim_mat <- simulateResiduals(mated_mod, plot = T)
testDispersion(sim_mat)
testZeroInflation(sim_mat)

mated_mod2 <- glmmTMB(value ~ Treatment*Day + (1|ID), ziformula = ~Day, data = rep_long, family = 'nbinom2')
summary(mated_mod2)
sim_mat2 <- simulateResiduals(mated_mod2, plot = T)
AIC(mated_mod2, mated_mod)

mated_mod3 <- glmmTMB(value ~ Treatment*Day + (1|ID), ziformula = ~Day, disp = ~Day, data = rep_long, family = 'nbinom2')
summary(mated_mod3)    
sim_mat3 <- simulateResiduals(mated_mod3, plot = T)
AIC(mated_mod, mated_mod2, mated_mod3)




# DRO Mated Lifespan ------------------------------------------------------

Mated_LS <- read.csv('Mated_LS.csv')
Mated_LS$Treatment.ID <- as.factor(Mated_LS$Treatment.ID)

#Combine Repro and LS plate columns to get 1 column of plate IDs
Mated_LS$Plate <- coalesce(Mated_LS$Repro.ID, Mated_LS$LS.plate)

#Remove these columns as they add unnecessary Nas to the dataset
Mated_LS = subset(Mated_LS, select = -c(Repro.ID, LS.plate))

Mated_LS <- na.omit(Mated_LS)

Mated_LS$Treatment.ID <- as.factor(Mated_LS$Treatment.ID)
Mated_LS$Treatment.ID <- droplevels(Mated_LS$Treatment.ID)
levels(Mated_LS$Treatment.ID)

Mated_LS$Treatment.ID <- revalue(Mated_LS$Treatment.ID, c('F' = 'F', 'FO' = 'F+O', 'DR' = 'DR', 'DRO' = 'DR+O'))


#Exclude lost or walled worms
Mated_LS2 <-Mated_LS %>%
  filter(Mated_LS$Cause != 'L' & Mated_LS$Cause !='W')

#Remove infected males (this infection appeared extremely pathogenic and detrimental, thus we decided to exclude these individuals)
Mated_LS2 <- Mated_LS2 %>%
  filter(Mated_LS2$Infect != 'Y')

surv<-survfit(Surv(Age,Event)~Treatment.ID,data=Mated_LS2)

Mated_LS_plot <-ggsurvplot(surv, ylab="Survival probability\n", data = Mated_LS2, size= 0.8, font.ylab= 18, font.xlab= 18, legend = c(0.2, 0.3), legend.title = "", title = "Matricides uncensored", censor = FALSE, xlab = "\nDay", xlim=c(0,25), break.time.by = 5, palette = palette, position= position_dodge(0.9), font.tickslab = c(14), legend.labs=c("DR" ,"DR+O", "F","F+O"), font.legend = c(16))
Mated_LS_plot

#Censor matricides
Mated_LS3 <- Mated_LS2 %>%
  filter(Mated_LS2$Cause != 'M')

surv<-survfit(Surv(Age,Event)~Treatment.ID,data=Mated_LS3)

Mated_LS_plot_matcen <-ggsurvplot(surv, ylab="Survival probability\n", data = Mated_LS3, size= 0.8, font.ylab= 18, font.xlab= 18, legend = c(0.2, 0.2), legend.title = "", palette = palette, title = "Matricides censored", censor = FALSE, xlab = "\nDay", xlim=c(0,25), break.time.by = 5, position= position_dodge(0.9), legend.labs=c("DR" ,"DR+O", "F","F+O"), font.tickslab = c(14), font.legend = c(16))
Mated_LS_plot_matcen

#Function for creating forest plot
meforest <- function(cox, F){  #Eds version of forest plot
  require(AICcmodavg)
  require(ggplot2)
  require(forcats)
  store <- matrix(nrow = length(cox$coefficients) + 1, ncol = 4)
  ref <- c(paste(F), 0,NA, NA)
  store[1,] <- ref
  for (x in 1:length(cox$coefficients)){
    y = x+1
    mean<-cox$coefficients[x]
    emean <- (mean)
    CIL <- cox$coefficients[x] - (1.96*extractSE(cox)[x])
    CUL <- cox$coefficients[x] + (1.96*extractSE(cox)[x])
    eCIL <- (CIL)
    eCUL <- (CUL)
    store[y, 1:4] <- c(names(cox$coefficients[x]), emean, eCIL, eCUL)}
  colnames(store) <- c("Treatment", "mean", "CIL","CUL")
  store <- as.data.frame(store)
  store$mean <- as.numeric(as.character(store$mean))
  store$CIL <- as.numeric(as.character(store$CIL))
  store$CUL <- as.numeric(as.character(store$CUL))
  forest<-ggplot(store, aes(x = mean, xmax = CUL, xmin = CIL, y = Treatment, colour = Treatment)) +
    geom_point(size=4, shape = 19, colour = c('#1F77B4FF', '#FF7F0EFF', '#2CA02CFF','#D62728FF')) +
    geom_errorbarh(height=0.25, size=0.9,  colour = c('#1F77B4FF', '#FF7F0EFF', '#2CA02CFF','#D62728FF')) +
    theme_minimal() +
    geom_vline(xintercept=0, linetype="dotted", size=0.8) +
    xlab("Hazard Ratio")+
    ylab("")+
    theme(
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 20)
    )+
    scale_y_discrete(limits = (levels(store$Treatment.ID)))+
    expand_limits(x = c(-1.5,0.5))
  return(forest)
}


cox <- coxme(Surv(Age, Event) ~ Treatment.ID + (1|Plate), data = Mated_LS2)
summary(cox)

forest <-meforest(cox, "F")
forest

cox_matcen <- coxme(Surv(Age, Event) ~ Treatment.ID +(1|Plate), data = DRO_LS3)
summary(cox_matcen)

forest_matcen <- meforest(cox_matcen, 'F')
forest_matcen

DR_all_plot <- ggarrange(LS_plot$plot, LS_plot_matcen$plot, forest, forest_matcen, nrow = 2, ncol = 2,heights = c(2, 1))
DR_all_plot
ggsave('DR_all_plot.tif', height = 8, width = 12)

#test cox models for fit
cox.zph(cox)
cox.zph(cox_matcen)







# DRO Heatshock Survival --------------------------------------------------
DR_HS <- read.csv('DR_HS.csv')
DR_HS$Age <- as.numeric(DR_HS$Age)
DR_HS <- na.omit(DR_HS)

DR_HS$Treatment <- as.factor(DR_HS$Treatment)
DR_HS$Treatment <- factor(DR_HS$Treatment, levels = c('F','FO','DR','DRO'))
#
DR_HS$Treatment <- relevel(DR_HS$Treatment, ref = 'F')



###Forest plot function
meforest <- function(cox, F){  #Eds version of forest plot
  require(AICcmodavg)
  require(ggplot2)
  require(forcats)
  store <- matrix(nrow = length(cox$coefficients) + 1, ncol = 4)
  ref <- c(paste(F), 0,NA, NA)
  store[1,] <- ref
  for (x in 1:length(cox$coefficients)){
    y = x+1
    mean<-cox$coefficients[x]
    emean <- (mean)
    CIL <- cox$coefficients[x] - (1.96*extractSE(cox)[x])
    CUL <- cox$coefficients[x] + (1.96*extractSE(cox)[x])
    eCIL <- (CIL)
    eCUL <- (CUL)
    store[y, 1:4] <- c(names(cox$coefficients[x]), emean, eCIL, eCUL)}
  colnames(store) <- c("Treatment", "mean", "CIL","CUL")
  store <- as.data.frame(store)
  store$mean <- as.numeric(as.character(store$mean))
  store$CIL <- as.numeric(as.character(store$CIL))
  store$CUL <- as.numeric(as.character(store$CUL))
  forest<-ggplot(store, aes(x = mean, xmax = CUL, xmin = CIL, y = Treatment, colour = Treatment)) +
    geom_point(size=4, shape = 19, colour = c('#1F77B4FF', '#FF7F0EFF', '#2CA02CFF','#D62728FF')) +
    geom_errorbarh(height=0.25, size=0.9,  colour = c('#1F77B4FF', '#FF7F0EFF', '#2CA02CFF','#D62728FF')) +
    theme_minimal() +
    geom_vline(xintercept=0, linetype="dotted", size=0.8) +
    xlab("Hazard Ratio")+
    ylab("")+
    theme(
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 20)
    )+
    scale_y_discrete(limits = (levels(store$Treatment)), labels = c('TreatmentDR'='DR', 'TreatmentDRO' = 'DR+O', 'TreatmentF' ='F', 'TreatmentFO' ='F+O'))+
    expand_limits(x = c(-1.5,0.5))
  return(forest)
}


#  scale_y_discrete(limits = (levels(store$Treatment)), labels = c('TreatmentF' ='F', 'TreatmentFO' ='F+O','TreatmentDR'='DR', 'TreatmentDRO' = 'DR+O'))+

cox <- coxme(Surv(Age, Event) ~ Treatment +(1|Plate.ID), data = DR_HS)
summary(cox)
cox.zph(cox)#it's OK

forest_HS <-meforest(cox, "F")
forest_HS

DR_HS$Treatment <- revalue(DR_HS$Treatment, c('F' = 'F', 'FO' = 'F+O', 'DR' = 'DR', 'DRO' = 'DR+O'))


levels(DR_HS$Treatment)


#Survival curve with line to divdie hours and days
surv<-survfit(Surv(Age,Event)~Treatment,data=DR_HS)
summary(surv)
surv <- na.omit(surv)
HS_plot <-ggsurvplot(surv, ylab="Survival probability\n", data = DR_HS, size= 0.8, font.ylab= 18, font.xlab= 18, legend = c(0.8, 0.9), legend.title = "", palette = palette, title = "", censor = FALSE, xlab = "\nHour                                                 Day", xlim=c(0,15), break.time.by = 2, position= position_dodge(0.9), font.tickslab = c(14), font.legend = c(16))

HS_plot <- HS_plot$plot + geom_vline(xintercept = 9, linetype = 'dashed', colour = 'red', size = 1)

DR_HS_plot <- ggarrange(HS_plot, forest_HS, nrow = 2, ncol = 1,heights = c(2, 1))
DR_HS_plot

ggsave('DR_HS_plot.tif', height = 12, width = 10)



# DRO Egg Size ------------------------------------------------------------

egg <- read.csv('Egg_size.csv')

egg$Treatment <- as.factor(egg$Treatment)
egg$Treatment <- factor(egg$Treatment, levels = c('F','FO','DR','DRO'))
egg$Treatment <- revalue(egg$Treatment, c('F' = 'F', 'FO' = 'F+O', 'DR' = 'DR', 'DRO' = 'DR+O'))

egg_2 <- egg %>%
  filter(Day == 2)

egg_4 <- egg %>%
  filter(Day == 4)

Egg_d4 <-
  egg_4 %>%
  load(Treatment, Egg_size,
       idx = list(c("F", "F+O", "DR", "DR+O")))

Egg_d4_dab <- mean_diff(Egg_d4)

Egg_d4_plot <- dabest_plot(Egg_d4_dab, FALSE, swarm_label = 'Area (mm^2)', raw_marker_spread = 1, custom_palette = 'd3')

Egg_d4_plot

egg_d4_mod <- lmerTest::lmer(Egg_size ~ Treatment + (1|ID), data = egg_4)
summary(egg_d4_mod)

Egg_d2 <-
  egg_2 %>%
  load(Treatment, Egg_size,
       idx = list(c("F", "F+O", "DR", "DR+O")))

Egg_d2_dab <- mean_diff(Egg_d2)

Egg_d2_plot <- dabest_plot(Egg_d2_dab, FALSE, swarm_label = 'Area (mm^2)', raw_marker_spread = 1, custom_palette = 'd3')

Egg_d2_plot

egg_plots <- ggarrange(Egg_d2_plot, Egg_d4_plot, ncol = 2, nrow = 1)

egg_plots

cowplot::plot_grid(
  plotlist = list(Egg_d2_plot, Egg_d4_plot),
  nrow = 1,
  ncol = 2,
  labels = c('.                             D2 Egg Size','.                            D4 Egg Size ')
)

ggsave('egg.tif', height = 8, width = 10)


head(egg)


egg_mod2 <- lmerTest::lmer(Egg_size ~ Treatment + (1|ID), data = egg_2)
summary(egg_mod2)

egg_mod4 <- lmerTest::lmer(Egg_size ~ Treatment + (1|ID), data = egg_4)
summary(egg_mod4)



# DRO Body Size -----------------------------------------------------------
binary <- read.csv('Binary_size.csv')

binary$Treatment <- as.factor(binary$Treatment)

binary$Treatment <- factor(binary$Treatment, levels = c('F','FO','DR','DRO'))

binary$Treatment <- revalue(binary$Treatment, c('F' = 'F', 'FO' = 'F+O', 'DR' = 'DR', 'DRO' = 'DR+O'))


Binary <-
  binary %>%
  load(Treatment, Body,
       idx = list(c("F", "F+O", "DR", "DR+O")))

Binary_dab <- mean_diff(Binary)

Binary_plot <- dabest_plot(Binary_dab, FALSE, swarm_label = 'Area (mm^2)', raw_marker_spread = 1, custom_palette = 'd3')

Binary_plot

binary_d4 <- binary %>%
  filter(Day == 4)

Binary_d4 <-
  binary_d4 %>%
  load(Treatment, Body,
       idx = list(c("F", "F+O", "DR", "DR+O")))

Binary_d4_dab <- mean_diff(Binary_d4)

Binary_d4_plot <- dabest_plot(Binary_d4_dab, FALSE, swarm_label = 'Area (mm^2)', raw_marker_spread = 1, custom_palette = 'd3')

Binary_d4_plot

binary_d2 <- binary %>%
  filter(Day == 2)

Binary_d2 <-
  binary_d2 %>%
  load(Treatment, Body,
       idx = list(c("F", "F+O", "DR", "DR+O")))

Binary_d2_dab <- mean_diff(Binary_d2)

Binary_d2_plot <- dabest_plot(Binary_d2_dab, FALSE, swarm_label = 'Area (mm^2)', raw_marker_spread = 1, custom_palette = 'd3')

Binary_d2_plot

binary = subset(binary, select = c(ID, Day, Treatment, Body))

binary_wide <- spread(binary, Day, Body)

names(binary_wide) <- c('ID', 'Treatment', 'D2','D4')
binary_wide

binary_wide$Growth <- binary_wide$D4 - binary_wide$D2

binary_growth <- 
  binary_wide %>%
  load(Treatment, Growth, 
       idx = list(c("F", "F+O", "DR", "DR+O")))

binary_growth_dab <- mean_diff(binary_growth)

binary_growth_plot <- dabest_plot(binary_growth_dab, FALSE, swarm_label = 'Area (mm^2)', raw_marker_spread = 1, custom_palette = 'd3')

binary_growth_plot

binary_means <- summarySE(measurevar = 'Body', groupvars = c('Treatment', 'Day'), data = binary)

binary$Day <- as.factor(binary$Day)
binary_means$Day <- as.factor(binary_means$Day)

growth_plot <- ggplot()+
  geom_point(data = binary_means, aes(x = Day, colour = Treatment, y = Body), position = position_dodge(width = 0.5), size = 3)+
  geom_errorbar(data = binary_means, aes(x = Day, colour = Treatment, ymin = Body - 2*se, ymax = Body + 2*se), position = position_dodge(width = 0.5), size = 1.5, width = 0.2)+
  geom_point(data = binary, aes(x = Day, colour = Treatment, y = Body), size = 2, alpha = 0.5, position = position_dodge(width = 0.5))+
  geom_line(data = binary_means, aes(x = Day, colour = Treatment, y = Body, group = Treatment), position = position_dodge(width = 0.5), size = 1.5)+
  theme_classic()+
  labs(y ='Area (mm^2)', x = '')+
  theme(axis.title.y = element_text(size=14), axis.text.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.position = 'top')+
  scale_color_manual(values = palette)+
  scale_x_discrete(labels = c('2' = 'Day 2', '4' = 'Day 4'))

growth_plot

growth_plot + binary_growth_plot

ggsave('growth.tif', height = 6, width = 10)

growth_mod <- lm(Growth ~ Treatment, data = binary_wide)
summary(growth_mod)


D2_mod <- lm(D2 ~ Treatment, data = binary_wide)
summary(D2_mod)

D4_mod <- lm(D4 ~ Treatment, data = binary_wide)
summary(D4_mod)
