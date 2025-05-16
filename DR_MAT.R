

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
  library(usethis)
})

palette <- c('#1F77B4FF' ,'#FF7F0EFF','#2CA02CFF', '#D62728FF' )


#edit_git_config()
#use_git()
#create_github_token()
#gitcreds_set()
#use_github()



#R.version


# DRO Lifespan ------------------------------------------------------------

DRO_LS <- read.csv('DRO_LS.csv')
DRO_LS <- na.omit(DRO_LS)
DRO_LS$Treatment.ID <- as.factor(DRO_LS$Treatment.ID)
DRO_LS$Treatment.ID <- droplevels(DRO_LS$Treatment.ID)
levels(DRO_LS$Treatment.ID)


DRO_LS$Treatment.ID <- relevel(DRO_LS$Treatment.ID, ref = 'F')

DRO_LS2 <- DRO_LS %>%
  filter(DRO_LS$Cause != 'L' & DRO_LS$Cause !='W' & DRO_LS$Cause != 'E')

#Survival plot
surv<-survfit(Surv(Age,Event)~Treatment.ID,data=DRO_LS2)


LS_plot <-ggsurvplot(surv, ylab="Survival probability\n", data = DRO_LS2, size= 0.8, font.ylab= 18, font.xlab= 18, legend = c(0.3, 0.4), legend.title = "", title = "Matricides uncensored", censor = FALSE, xlab = "\nDay", xlim=c(0,32), break.time.by = 5, palette = palette, position= position_dodge(0.9), legend.labs=c("F" ,"DR", "DR+O","F+O"),font.tickslab = c(14),  font.legend = c(16), font.title = c(16))
LS_plot

#Exclude matricides
DRO_LS3 <- DRO_LS2 %>%
  filter(DRO_LS2$Cause != 'M')

#Survival plot without matricides
surv<-survfit(Surv(Age,Event)~Treatment.ID,data=DRO_LS3)

LS_plot_matcen <-ggsurvplot(surv, ylab="Survival probability\n", data = DRO_LS3, size= 0.8, font.ylab= 18, font.xlab= 18, legend = c(0.3, 0.4), legend.title = "", palette = palette, title = "Matricides censored", censor = FALSE, xlab = "\nDay", xlim=c(0,32), break.time.by = 5, position= position_dodge(0.9), legend.labs=c("F" ,"DR", "DR+O","F+O"), font.tickslab = c(14), font.legend = c(16), font.title = c(16))
LS_plot_matcen

#Relevel Treatment so F control is reference level
DRO_LS3$Treatment.ID <- as.factor(DRO_LS3$Treatment.ID)
DRO_LS3$Treatment.ID <- relevel(DRO_LS3$Treatment.ID, ref = 'F')
DRO_LS2$Treatment.ID <- relevel(DRO_LS2$Treatment.ID, ref = 'F')

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
    geom_point(size=4, shape = 19, colour = c('#1F77B4FF', '#2CA02CFF', '#D62728FF','#FF7F0EFF')) +
    geom_errorbarh(height=0.25, size=0.9,  colour = c('#1F77B4FF', '#2CA02CFF', '#D62728FF','#FF7F0EFF')) +
    theme_minimal() +
    geom_vline(xintercept=0, linetype="dotted", size=0.8) +
    xlab("Hazard Ratio")+
    ylab("")+
    theme(
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 20)
    )+
    scale_y_discrete(limits = (levels(store$Treatment)), labels = c('Treatment.IDDR'='DR', 'Treatment.IDDRO' = 'DR+O', 'Treatment.IDF' ='F', 'Treatment.IDFO' ='F+O'))+
    expand_limits(x = c(-1.5,0.5))
  return(forest)
}

#Survival analysis with cox model of lifespan data with matricides
cox <- coxme(Surv(Age, Event) ~ Treatment.ID +(1|Plate.ID), data = DRO_LS2)
summary(cox)

#Pairwise comparisions
cox_emm <- emmeans(cox, 'Treatment.ID')
pairs(cox_emm)

#Create forest plot of cox model
forest <-meforest(cox, "F")
forest

#Survival analysis with cox model of lifespan data without matricides
cox_matcen <- coxme(Surv(Age, Event) ~ Treatment.ID +(1|Plate.ID), data = DRO_LS3)
summary(cox_matcen)

#Pairwise comparisons
emm_cox_matcen <- emmeans(cox_matcen, 'Treatment.ID')
pairs(emm_cox_matcen)

#create forest plot
forest_matcen <- meforest(cox_matcen, 'F')
forest_matcen

#test cox models for fit
cox.zph(cox)#looks good - use cox model results
cox.zph(cox_matcen)#looks good

#Combine plots
DR_all_plot <- ggarrange(LS_plot$plot, LS_plot_matcen$plot, forest, forest_matcen, common.legend = T, nrow = 2, ncol = 2,heights = c(2, 1), labels = c('A', 'B'))
DR_all_plot
ggsave('DR_all_plot.tif', height = 8, width = 12)



# DRO Reproduction --------------------------------------------------------


DR_rep <- read.csv('DR_Lab_repro_2.csv')

#Long form data
rep_long <- DR_rep %>% 
  pivot_longer(
    cols = `D1`:`D8`, 
    names_to = "Day",
    values_to = "value"
  )

#Change factor level order of treatment
rep_long$Treatment <- factor(rep_long$Treatment, levels = c('F','FO','DR','DRO'))

#Change factor level names for treatment
rep_long$Treatment <- revalue(rep_long$Treatment, c('F' = 'F', 'FO' = 'F+O', 'DR' = 'DR', 'DRO' = 'DR+O'))


## Age-specific reproduction plot
rep <-ggplot(data=rep_long, aes(x=factor(Day), y=value, group=Treatment, color=Treatment))+
  geom_jitter(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5))+
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", size = 1, width=0.0, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="line",  size=1, position = position_dodge(0.5)) +
  theme_classic()+
  labs(y="Offspring number", x="", size=14)+
  labs(col="")+
  theme(axis.title.y = element_text(size=14))+
  theme(axis.title.x = element_text(size=14))+
  theme(axis.text = element_text(size = 12))+
  theme(legend.text = element_text(size = 12))+
  theme(legend.key.width = unit(0.5,"cm"))+
  coord_cartesian(ylim = c(0,170))+
  theme(legend.position = c(0.8,0.9))+
  scale_colour_manual(values = palette)

rep


#Sum number of offspring over 5 days for individuals
totalrep<-na.omit(as.data.frame.table(tapply(rep_long$value,list(rep_long$Treatment, rep_long$ID),sum)))

#rename columns
names(totalrep)<-c("Treatment", "Replicate", "Totrep")

#Get means for LRS
totalrep_means <- summarySE(data = totalrep, measurevar = 'Totrep', groupvars = 'Treatment')

#Bootstrap comparison plots
totrep_dab <-
  totalrep %>%
  load(x =Treatment, y=Totrep,
       idx = c('F','F+O','DR','DR+O'))

totrep_p <- mean_diff(totrep_dab)

tot_rep_plot <- dabest_plot(totrep_p, FALSE, raw_marker_spread = 1, swarm_label = 'LRS',custom_palette = 'd3', swarm_x_text = 12, swarm_y_text = 12, contrast_y_text = 12, contrast_x_text = 12, raw_marker_alpha = 0.3, tufte_size = 1 )

tot_rep_plot

#Analyse LRS
hist(totalrep$Totrep)
tot_rep_mod <- lm(Totrep ~ Treatment, data = totalrep)
summary(tot_rep_mod)

#Pairwise comparisons
emm_LRS <- emmeans(tot_rep_mod, 'Treatment')
pairs(emm_LRS)


#Combine age-specific repro and LRS plots
DR_rep_plot <- ggarrange(rep, tot_rep_plot+theme(plot.margin = unit(c(10,30,0,10), 'pt')), ncol = 2, labels = c('A','B'), widths = c(1, 0.5))
DR_rep_plot

ggsave('DRO_rep_plot.tif', height = 5, width = 10)


#Change Day to numeric value for analysis
rep_long$Day <- revalue(rep_long$Day, c('D1' = '1', 'D2' = '2', 'D3' = '3', 'D4' = '4', 'D5' = '5', 'D6' = '6', 'D7' = '7', 'D8' = '8'))
rep_long$Day <- as.numeric(rep_long$Day)

#Add OLRE column
rep_long <- tibble::rowid_to_column(rep_long, "OLRE")


#Models
rep_mod <- glmmTMB(value ~ Treatment * poly(Day,2)+ (1|ID), data = rep_long, family = 'poisson')
summary(rep_mod)
sim1 <- simulateResiduals(rep_mod, plot = T)
testDispersion(sim1)
check_overdispersion(rep_mod)
testZeroInflation(rep_mod)

rep_mod2 <- glmmTMB(value ~ Treatment * poly(Day, 2) + (1|ID) + (1|OLRE), data = rep_long, ziformula = ~Day, family = 'poisson')
summary(rep_mod2)
sim2 <- simulateResiduals(rep_mod2, plot = T)
AIC(rep_mod, rep_mod2)                  
testDispersion(sim2)
check_overdispersion(rep_mod2)

rep_mod3 <- glmmTMB(value ~ Treatment * poly(Day, 2) + (1|ID), data = rep_long, ziformula = ~Day, dispformula = ~Day, family = 'poisson')
summary(rep_mod3)
sim3 <- simulateResiduals(rep_mod3, plot = T)
testDispersion(sim3)
testZeroInflation(sim3)

rep_mod4 <- glmmTMB(value ~ Treatment * poly(Day, 2) + (1|ID), data = rep_long, ziformula = ~Treatment*Day, family = 'poisson')
summary(rep_mod4)
sim4 <- simulateResiduals(rep_mod4, plot = T)
testDispersion(sim4)
testZeroInflation(sim4)#

rep_mod4 <- glmmTMB(value ~ Treatment * poly(Day, 2) + (1|ID) + (1|OLRE), data = rep_long, ziformula = ~Treatment+I(Day^2), family = 'poisson', control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
summary(rep_mod4)


rep_mod5 <- glmmTMB(value ~ Treatment *Day + Treatment*I(Day^2) + (1|ID), data = rep_long, family = 'nbinom2')
summary(rep_mod5)
sim5 <- simulateResiduals(rep_mod5, plot = T)
testDispersion(sim5)
check_overdispersion(sim5)
testZeroInflation(sim5)
AIC(rep_mod5, rep_mod10)

rep_mod6 <- glmmTMB(value ~ Treatment*Day + Treatment * I(Day^2) + (1|ID), data = rep_long, ziformula = ~Day, family = 'nbinom2')
summary(rep_mod6)
AIC(rep_mod5, rep_mod6)
sim6 <- simulateResiduals(rep_mod6, plot = T)
testDispersion(sim6)
check_overdispersion(rep_mod6)
testZeroInflation(rep_mod6)


rep_mod7 <- glmmTMB(value ~ Treatment* Day + Treatment*I(Day^2) + (1|ID), data = rep_long, ziformula = ~I(Day^2), family = 'nbinom2')
summary(rep_mod7)
sim7 <- simulateResiduals(rep_mod7, plot = T)
AIC(rep_mod6, rep_mod7)
testDispersion(sim7)
testZeroInflation(sim7)

rep_mod8 <- glmmTMB(value ~ Treatment*Day + Treatment*I(Day^2) + (1|ID), data = rep_long, ziformula = ~Treatment*Day, family = 'nbinom2')
summary(rep_mod8)
AIC(rep_mod7, rep_mod8)
sim8 <- simulateResiduals(rep_mod8, plot = T)
testZeroInflation(sim8)#doesn't fix zero inflation


rep_mod9 <- glmmTMB(value ~ Treatment*Day + Treatment*I(Day^2) + (1|ID), data = rep_long, ziformula = ~Treatment, family = 'nbinom2')
AIC(rep_mod7, rep_mod9)
summary(rep_mod9)
sim9 <- simulateResiduals(rep_mod9, plot = T)
testDispersion(sim9)
testZeroInflation(sim9)


Anova(rep_mod9, type = 'III')


# DRO Mated Reproduction --------------------------------------------------

male_rep <- read.csv('Male_repro_DR.csv')

#Exclude individuals lost early due to infection/human error
male_rep <- male_rep %>%
  filter(Lost != 'INF' & Lost != 'L')

#Rename treatment levels
male_rep$Treatment <- revalue(male_rep$Treatment, c('F' = 'F', 'FO' = 'F+O', 'DR' = 'DR', 'DRO' = 'DR+O'))

#Long form data
rep_long <- male_rep %>% 
  pivot_longer(
    cols = `D1`:`D10`, 
    names_to = "Day",
    values_to = "value"
  )

#Change level order of treatment
rep_long$Treatment <- factor(rep_long$Treatment, levels = c('F','F+O','DR','DR+O'))


#Change level order of Day
rep_long$Day <- factor(rep_long$Day, levels = c('D1','D2','D3','D4','D5','D6','D7','D8','D9','D10'))

## Age-specific reproduction plot
rep <-ggplot(data=rep_long, aes(x=factor(Day), y=value, group=Treatment, color=Treatment))+
  geom_jitter(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5))+
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", size = 1, width=0.0, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="line",  size=1, position = position_dodge(0.5)) +
  theme_classic()+
  labs(y="Offspring number", x="", size=16)+
  labs(col="")+
  theme(axis.title.y = element_text(size=14))+
  theme(axis.title.x = element_text(size=14))+
  theme(axis.text = element_text(size = 12))+
  theme(legend.text = element_text(size = 12))+
  theme(legend.key.width = unit(0.5,"cm"))+
  coord_cartesian(ylim = c(0,170))+
  theme(legend.position = c(0.8,0.9))+
  scale_colour_manual(values = palette)

rep

#Sum number of offspring over all days for individuals
totalrep<-na.omit(as.data.frame.table(tapply(rep_long$value,list(rep_long$Treatment, rep_long$ID),sum)))


#rename columns
names(totalrep)<-c("Treatment", "Replicate", "Totrep")


#Bootstrap estimation plots
totrep_dab <-
  totalrep %>%
  load(x =Treatment, y=Totrep,
       idx = c('F','F+O','DR','DR+O'))

totrep_p <- mean_diff(totrep_dab)

tot_rep_plot <- dabest_plot(totrep_p, FALSE, raw_marker_spread = 1, swarm_label = 'LRS',custom_palette = 'd3', swarm_x_text = 12, swarm_y_text = 12, contrast_y_text = 12, contrast_x_text = 12, raw_marker_alpha = 0.3, tufte_size = 1)

tot_rep_plot

#analyse LRS for mated nematodes
Mated_tot <- lm(Totrep ~ Treatment, data = totalrep)
summary(Mated_tot)
mated_tot_sim <- simulateResiduals(Mated_tot, plot = T)

emm_mated_LRS <- emmeans(Mated_tot, 'Treatment')
pairs(emm_mated_LRS)

#Combine plots
DR_Mrep_no_L <- ggarrange(rep, tot_rep_plot+theme(plot.margin = unit(c(10,30,0,10), 'pt')), ncol = 2, labels = c('A','B'), widths = c(1, 0.5))
DR_Mrep_no_L

ggsave('DRO_Mrep_plot_no_L.tif', height = 5, width = 10)

#Analyse age-specific reproduction 

#Make Day numeric
rep_long$Day <- revalue(rep_long$Day, c('D1' = '1', 'D2' = '2', 'D3' = '3', 'D4' = '4', 'D5' = '5', 'D6' = '6', 'D7' = '7', 'D8' = '8', 'D9' = '9', 'D10' = '10'))
rep_long$Day <- as.numeric(rep_long$Day)

hist(rep_long$value)

#Models
mated_mod <- glmmTMB(value ~ Treatment*poly(Day,2) + (1|ID), data = rep_long, family = 'poisson')
summary(mated_mod)
sim_mat <- simulateResiduals(mated_mod, plot = T)
testDispersion(sim_mat)
testZeroInflation(sim_mat)
check_overdispersion(mated_mod)

mated_mod2 <- glmmTMB(value ~ Treatment*poly(Day, 2) + (1|ID), ziformula = ~poly(Day,2), data = rep_long, family = 'poisson')
summary(mated_mod2)
sim_mat2 <- simulateResiduals(mated_mod2, plot = T)
AIC(mated_mod2, mated_mod)
testZeroInflation(sim_mat2)
testDispersion(sim_mat2)
check_overdispersion(mated_mod2)

rep_long <- tibble::rowid_to_column(rep_long, "OLRE")

mated_mod3 <- glmmTMB(value ~ Treatment*poly(Day,2) + (1|ID) + (1|OLRE), ziformula = ~poly(Day,2), data = rep_long, family = 'poisson')
summary(mated_mod3)    
sim_mat3 <- simulateResiduals(mated_mod3, plot = T)
AIC(mated_mod2, mated_mod3)
testDispersion(sim_mat3)
testZeroInflation(sim_mat3)

mated_mod4 <- glmmTMB(value ~ Treatment*Day + Treatment*I(Day^2) + (1|ID), ziformula = ~poly(Day,2), data = rep_long, family = 'nbinom2')
summary(mated_mod4)
sim4 <- simulateResiduals(mated_mod4, plot = T)
testDispersion(sim4)
testZeroInflation(sim4)


mated_mod5 <- glmmTMB(value ~ Treatment*Day + Treatment*I(Day^2) + (1|ID), ziformula = ~I(Day^2), data = rep_long, family = 'nbinom2')
summary(mated_mod5)
sim5 <- simulateResiduals(mated_mod5, plot = T)
testDispersion(sim5)
testZeroInflation(sim5)

mated_mod6 <- glmmTMB(value ~ Treatment*Day + Treatment*I(Day^2) + (1|ID), ziformula = ~Day, data = rep_long, family = 'nbinom2')
summary(mated_mod6)
sim6 <- simulateResiduals(mated_mod6, plot = T)
testDispersion(sim6)
testZeroInflation(sim6)
AIC(mated_mod4, mated_mod5, mated_mod6)

Anova(mated_mod6, type = 'III')


mated_mod7 <- glmmTMB(value ~ Treatment*Day + Treatment*I(Day^2) + (1|ID), ziformula = ~Treatment*Day, data = rep_long, family = 'nbinom2')
summary(mated_mod7)
AIC(mated_mod6, mated_mod7)
sim7 <- simulateResiduals(mated_mod7, plot = T)
testZeroInflation(sim7)
testDispersion(sim7)
# DRO Mated Lifespan
Anova(mated_mod7, type = 'III')
mated_mod8 <- glmmTMB(value ~ Treatment*Day + Treatment*I(Day^2) + (1|ID), ziformula = ~Treatment+Day, data = rep_long, family = 'nbinom2')
summary(mated_mod8)
AIC(mated_mod8, mated_mod7)
sim7 <- simulateResiduals(mated_mod7, plot = T)
testZeroInflation(sim7)
testDispersion(sim7)


mated_mod9 <- glmmTMB(value ~ Treatment*Day + Treatment*I(Day^2) + (1|ID), ziformula = ~Treatment, data = rep_long, family = 'nbinom2')
summary(mated_mod9)
AIC(mated_mod7, mated_mod9)
sim7 <- simulateResiduals(mated_mod9, plot = T)
testZeroInflation(sim9)
testDispersion(sim9)
##stick with 7


# DRO Heatshock Survival --------------------------------------------------
DR_HS <- read.csv('DR_HS.csv')

#Make age numeric
DR_HS$Age <- as.numeric(DR_HS$Age)
DR_HS <- na.omit(DR_HS)

#Treatment as factor and change level order
DR_HS$Treatment <- as.factor(DR_HS$Treatment)
DR_HS$Treatment <- factor(DR_HS$Treatment, levels = c('F','FO','DR','DRO'))
#Make control reference level
DR_HS$Treatment <- relevel(DR_HS$Treatment, ref = 'F')



###Forest plot function
meforest <- function(cox, F){  
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

#Survival model
cox <- coxme(Surv(Day, Event) ~ Treatment +(1|Plate.ID), data = DR_HS)
summary(cox)
cox.zph(cox)#it's OK

#Forest plot from results
forest_HS <-meforest(cox, "F")
forest_HS

#Change treatment factor level names
DR_HS$Treatment <- revalue(DR_HS$Treatment, c('F' = 'F', 'FO' = 'F+O', 'DR' = 'DR', 'DRO' = 'DR+O'))


#Survival curve with line to divide hours and days
surv<-survfit(Surv(Age,Event)~Treatment,data=DR_HS)
summary(surv)
surv <- na.omit(surv)
HS_plot <-ggsurvplot(surv, ylab="Survival probability\n", data = DR_HS, size= 0.8, font.ylab= 18, font.xlab= 18, legend = c(0.8, 0.9), legend.title = "", palette = palette, title = "", censor = FALSE, xlab = "\nHour                                        Day", xlim=c(0,21), break.time.by = 2, position= position_dodge(0.9), font.tickslab = c(14), font.legend = c(16), legend.labs = c('F', 'F+O', 'DR', 'DR+O'))

HS_plot <- HS_plot$plot + geom_vline(xintercept = 9, linetype = 'dashed', colour = 'red', size = 1)

HS_plot <- HS_plot + scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15,17,19,21),
                             labels = c(1,3,5,7,9,2,4,6,8,10,12))

DR_HS_plot <- ggarrange(HS_plot, forest_HS, nrow = 2, ncol = 1,heights = c(1.8, 1))
DR_HS_plot

ggsave('DR_HS_plot.tif', height = 8, width = 8)



# DRO Egg Size ------------------------------------------------------------

egg <- read.csv('Egg_size.csv')

#Change factor levels for Treatment
egg$Treatment <- as.factor(egg$Treatment)
egg$Treatment <- factor(egg$Treatment, levels = c('F','FO','DR','DRO'))
egg$Treatment <- revalue(egg$Treatment, c('F' = 'F', 'FO' = 'F+O', 'DR' = 'DR', 'DRO' = 'DR+O'))

#only day 2 egg sizes
egg_2 <- egg %>%
  filter(Day == 2)

#Only day 4 egg sizes
egg_4 <- egg %>%
  filter(Day == 4)

#Get means
egg_2_means <- summarySE(data = egg_2, measurevar = 'Egg_size', groupvars = 'Treatment')

egg_4_means <- summarySE(data = egg_4, measurevar = 'Egg_size', groupvars = 'Treatment')


#Bootstrap estimation plots
Egg_d4 <-
  egg_4 %>%
  load(Treatment, Egg_size,
       idx = list(c("F", "F+O", "DR", "DR+O")))

Egg_d4_dab <- mean_diff(Egg_d4)

Egg_d4_plot <- dabest_plot(Egg_d4_dab, FALSE, swarm_label = expression('Area mm'^2), raw_marker_spread = 1, custom_palette = 'd3', swarm_x_text = 12, swarm_y_text = 12, contrast_y_text = 12, contrast_x_text = 12, raw_marker_alpha = 0.3)

Egg_d4_plot




#Bootstrap estimation plots
Egg_d2 <-
  egg_2 %>%
  load(Treatment, Egg_size,
       idx = list(c("F", "F+O", "DR", "DR+O")))

Egg_d2_dab <- mean_diff(Egg_d2)

Egg_d2_plot <- dabest_plot(Egg_d2_dab, FALSE, swarm_label = expression('Area mm'^2), raw_marker_spread = 1.2, custom_palette = 'd3', swarm_x_text = 12, swarm_y_text = 12, contrast_y_text = 12, contrast_x_text = 12, raw_marker_alpha = 0.3)

Egg_d2_plot

#Group plots
egg_plots <- ggarrange(Egg_d2_plot+theme(plot.margin = unit(c(0,30,0,50), 'pt')), labels = c('C', 'D'), Egg_d4_plot+theme(plot.margin = unit(c(0,30,0,50), 'pt')), ncol = 2)

egg_plots

ggsave('egg.tif', height = 4, width = 8)

#Combine plots with repro plots 
all_rep_plot <- ggarrange(DR_rep_plot, egg_plots, nrow = 2)
all_rep_plot

ggsave('all_DR_rep.tif', height = 10, width = 10)


#Models 
egg_mod2 <- lmerTest::lmer(Egg_size ~ Treatment + (1|ID), data = egg_2)
summary(egg_mod2)
emm_egg2 <- emmeans(egg_mod2, 'Treatment')
pairs(emm_egg2)
Anova(egg_mod2)


egg_mod4 <- lmerTest::lmer(Egg_size ~ Treatment + (1|ID), data = egg_4)
summary(egg_mod4)
emm_egg4 <- emmeans(egg_mod4, 'Treatment')
pairs(emm_egg4)
Anova(egg_mod4)

# DRO Body Size -----------------------------------------------------------
binary <- read.csv('Binary_size.csv')

#Change factor levels for treatment
binary$Treatment <- as.factor(binary$Treatment)

binary$Treatment <- factor(binary$Treatment, levels = c('F','FO','DR','DRO'))

binary$Treatment <- revalue(binary$Treatment, c('F' = 'F', 'FO' = 'F+O', 'DR' = 'DR', 'DRO' = 'DR+O'))

#Bootstrap estimation plots
Binary <-
  binary %>%
  load(Treatment, Body,
       idx = list(c("F", "F+O", "DR", "DR+O")))

Binary_dab <- mean_diff(Binary)

Binary_plot <- dabest_plot(Binary_dab, FALSE, swarm_label = expression('Area mm'^2), raw_marker_spread = 1, custom_palette = 'd3', swarm_x_text = 12, swarm_y_text = 12, contrast_y_text = 12, contrast_x_text = 12, raw_marker_alpha = 0.3, tufte_size = 1)

Binary_plot

#Day 4 body size only
binary_d4 <- binary %>%
  filter(Day == 4)

#Day 4 size bootstrap est plots
Binary_d4 <-
  binary_d4 %>%
  load(Treatment, Body,
       idx = list(c("F", "F+O", "DR", "DR+O")))

Binary_d4_dab <- mean_diff(Binary_d4)

Binary_d4_plot <- dabest_plot(Binary_d4_dab, FALSE, swarm_label = 'Area (mm^2)', raw_marker_spread = 1, custom_palette = 'd3')

Binary_d4_plot

#Day 2 body size only
binary_d2 <- binary %>%
  filter(Day == 2)


#Day 2 bootstrap est plots
Binary_d2 <-
  binary_d2 %>%
  load(Treatment, Body,
       idx = list(c("F", "F+O", "DR", "DR+O")))

Binary_d2_dab <- mean_diff(Binary_d2)

Binary_d2_plot <- dabest_plot(Binary_d2_dab, FALSE, swarm_label = 'Area (mm^2)', raw_marker_spread = 1, custom_palette = 'd3')

Binary_d2_plot

#Keep only needed columns
binary = subset(binary, select = c(ID, Day, Treatment, Body))

#Wide form
binary_wide <- spread(binary, Day, Body)

names(binary_wide) <- c('ID', 'Treatment', 'D2','D4')
binary_wide

#Growth between D2 and D4
binary_wide$Growth <- binary_wide$D4 - binary_wide$D2

binary_wide <- binary_wide[binary_wide$Growth >= 0, ]

#Growth bootstrap est plots
binary_growth <- 
  binary_wide %>%
  load(Treatment, Growth, 
       idx = list(c("F", "F+O", "DR", "DR+O")))

binary_growth_dab <- mean_diff(binary_growth)

binary_growth_plot <- dabest_plot(binary_growth_dab, FALSE, swarm_label = expression('Area mm'^2), raw_marker_spread = 1, custom_palette = 'd3', swarm_x_text = 12, swarm_y_text = 12, contrast_y_text = 12, contrast_x_text = 12, raw_marker_alpha = 0.3, tufte_size = 1)

binary_growth_plot

#Get means
binary_means <- summarySE(measurevar = 'Body', groupvars = c('Treatment', 'Day'), data = binary)

#Day as factor for plot
binary$Day <- as.factor(binary$Day)
binary_means$Day <- as.factor(binary_means$Day)

#Growth plot
growth_plot <- ggplot()+
  geom_point(data = binary_means, aes(x = Day, colour = Treatment, y = Body), position = position_dodge(width = 0.5), size = 3)+
  geom_errorbar(data = binary_means, aes(x = Day, colour = Treatment, ymin = Body - 2*se, ymax = Body + 2*se), position = position_dodge(width = 0.5), size = 1.2, width = 0.2)+
  geom_point(data = binary, aes(x = Day, colour = Treatment, y = Body), size = 2, alpha = 0.3, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.15))+
  geom_line(data = binary_means, aes(x = Day, colour = Treatment, y = Body, group = Treatment), position = position_dodge(width = 0.5), size = 1.5)+
  theme_classic()+
  labs(y =expression('Area mm'^2), x = '')+
  theme(axis.title.y = element_text(size=14), axis.text.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.position = 'top')+
  scale_color_manual(values = palette)+
  scale_x_discrete(labels = c('2' = 'Day 2', '4' = 'Day 4'))


growth_plot

#Combine plots
growth_plots <- ggarrange(growth_plot+theme(plot.margin = unit(c(0,2,20,1), 'pt')), labels = c('A', 'B'), binary_growth_plot+theme(plot.margin = unit(c(0,2,0,1), 'pt')), ncol = 2, common.legend = TRUE, widths = c(1.2, 1))

growth_plots

ggsave('growth.tif', height = 6, width = 10)

#Models
growth_mod <- lm(Growth ~ Treatment, data = binary_wide)
summary(growth_mod)


D2_mod <- lm(D2 ~ Treatment, data = binary_wide)
summary(D2_mod)

D4_mod <- lm(D4 ~ Treatment, data = binary_wide)
summary(D4_mod)

#Post-hoc tests
emm_D2size <- emmeans(D2_mod, 'Treatment')
pairs(emm_D2size)

emm_D4size <- emmeans(D4_mod, 'Treatment')
pairs(emm_D4size)

emm_growth <- emmeans(growth_mod, 'Treatment')
pairs(emm_growth)


# Outdoor -----------------------------------------------------------------

DR <- read.csv('DR_wild2.csv')

#Change treatment factor levels
DR$Treatment <- factor(DR$Treatment, levels = c('F','FO','DR','DRO'))
DR$Treatment <- revalue(DR$Treatment, c('F' = 'F', 'FO' = 'F+O', 'DR' = 'DR', 'DRO' = 'DR+O'))




#Get means per population and per treatment
DR <- na.omit(DR)
DR_means_pop <- summarySE(DR, measurevar = 'No_worms', groupvars = c('Treatment','Population','Day'))
DR_means <- summarySE(DR, measurevar = 'No_worms', groupvars = c('Treatment', 'Day'))

DR_B1 <- DR %>%
  filter(Block == 1)
B1_means_pop <- summarySE(DR_B1, measurevar = 'No_worms', groupvars = c('Treatment','Population','Day'))
B1_means <- summarySE(DR_B1, measurevar = 'No_worms', groupvars = c('Treatment','Day'))


DR_B2 <- DR %>%
  filter(Block == 2)
B2_means_pop <- summarySE(DR_B2, measurevar = 'No_worms', groupvars = c('Treatment','Population','Day'))
B2_means <- summarySE(DR_B2, measurevar = 'No_worms', groupvars = c('Treatment','Day'))


DR_B3 <- DR %>%
  filter(Block == 3)
B3_means_pop <- summarySE(DR_B3, measurevar = 'No_worms', groupvars = c('Treatment','Population','Day'))
B3_means <- summarySE(DR_B3, measurevar = 'No_worms', groupvars = c('Treatment','Day'))

range(DR_means_pop$No_worms)
#Function for population plots including raw data
rawpop_plot_fun <- function(treat_means, pop_means){ggplot()+
    geom_point(data = treat_means, aes(x = Day, y = No_worms, colour = Treatment), size = 2.25, position = position_dodge(width = 1.2))+
    geom_line(data = treat_means, aes(x = Day, y = No_worms, colour = Treatment), size = 1.1, position = position_dodge(width = 1.2))+
    geom_point(data = pop_means, aes(x = Day, y = No_worms, colour = Treatment), size = 1.5, alpha = 0.1, position = position_jitterdodge(dodge.width = 1.2, jitter.width = 0.2))+
    geom_errorbar(data = treat_means, aes(x= Day, ymin=No_worms-se, ymax=No_worms+se, colour = Treatment), width = 0.2, size = 1.1, position = position_dodge(width = 1.2))+
    ylab('Population index (~ 95% CI)\n')+
    ylim(0, 475)+
    xlab('\nSample day')+
    theme(legend.title = element_blank(),
          legend.background = element_rect(fill = 'white', colour = 'white'),
          legend.key = element_rect(fill ='white', colour = 'white'),
          legend.text = element_text(size = (14)),
          legend.position = 'top',
          axis.title = element_text(size = (16)),
          axis.text = element_text(size = (14)),
          panel.background = element_rect(fill = 'white', colour = 'white'),   
          panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white"),
          axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"))+
    scale_color_manual( values = palette)+
    scale_x_continuous(breaks=c(7,14, 21))
  
}


#Function for population plots without raw data
pop_plot_fun <- function(means, title){ggplot()+
    geom_point(data = means, aes(x = Day, y = No_worms, colour = Treatment), size = 4, position = position_dodge(width = 0.5))+
    geom_line(data = means, aes(x = Day, y = No_worms, colour = Treatment), linewidth = 1.8, position = position_dodge(width = 0.5))+
    geom_errorbar(data = means, aes(x= Day, ymin=No_worms-(2*se), ymax=No_worms+(2*se), colour = Treatment), width = 0.2, linewidth = 1.1, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2))+
    ylab('Population index (~95% CI)')+
    xlab('Sample day')+
    labs(title = title)+
    theme(plot.title = element_text(colour = 'grey50', size = (10)),
          legend.title = element_blank(),
          legend.background = element_rect(fill = 'white', colour = 'white'),
          legend.key = element_rect(fill ='white', colour = 'white'),
          legend.text = element_text(size = (16)),
          legend.position = 'bottom',
          axis.title = element_text(size = (16)),
          axis.text = element_text(size = (14)),
          panel.background = element_rect(fill = 'white', colour = 'white'),   
          panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_line(colour = "white"),
          axis.line = element_line(linewidth =  0.5, linetype = "solid",colour = "black"))+
    scale_color_manual( values = palette)+
    scale_x_continuous(breaks=c(7,14, 21))}

#Plot including all blocks
DR_raw <- rawpop_plot_fun(DR_means, DR_means_pop)
DR_raw
ggsave('DR_pop_plot.tif', width = 7, height = 10)


#INdividual block plots
B1_raw <- rawpop_plot_fun(B1_means, B1_means_pop)
B1_raw <- B1_raw+ylim(0,450)
B1_raw
ggsave('B1_pop_plot.tif', width = 4, height = 5)


B2_raw <- rawpop_plot_fun(B2_means, B2_means_pop)
B2_raw <- B2_raw+ylim(0, 115)
B2_raw
ggsave('B2_pop_plot.tif', width = 4, height = 5)


B3_raw <- rawpop_plot_fun(B3_means, B3_means_pop)
B3_raw <- B3_raw+ylim(0,360)
B3_raw
ggsave('B3_pop_plot.tif', width = 4, height = 5)

#Combine block pkots
blocks <- ggarrange(B1_raw +rremove('axis.title')+rremove('legend')+theme(plot.margin = unit(c(0,0,0,20), 'pt')), B2_raw +rremove('legend')+rremove('axis.title')+theme(plot.margin = unit(c(0,0,0,20), 'pt')), B3_raw+rremove('legend') + rremove('axis.title')+theme(plot.margin = unit(c(0,0,0,20), 'pt')), nrow = 3, labels = c('B','C','D'), heights = c(1,1,1))

blocks

#Combine all plots
all_DR <- ggarrange(DR_raw+rremove('axis.title'), blocks, ncol = 2, labels = c('A'), common.legend = TRUE)
all_DR

all_DR <- annotate_figure(all_DR, 
                          left = text_grob('Population index [95% ci]', rot = 90, size = 16),
                          bottom = text_grob('Sample day', size = 16))

all_DR

ggsave('all_DR_outdoor.tif', height = 6, width = 8)

#Change block from numeric to factor
DR$Block <- as.factor(DR$Block)


#Models#

hist(DR$No_worms) #lots of zeroes

DR_mod <- glmmTMB(No_worms ~ Treatment * poly(Day,2) * Block + (1|Population/Rep), data = DR, family = 'poisson')
summary(DR_mod)
DR_sim1 <- simulateResiduals(DR_mod, plot = T)
testDispersion(DR_sim1)
testZeroInflation(DR_sim1)

DR_mod2 <- glmmTMB(No_worms ~ Treatment * poly(Day,2) * Block + (1|Population/Rep), data = DR, family = 'nbinom2')
summary(DR_mod2)
AIC(DR_mod, DR_mod2)

DR_sim2 <- simulateResiduals(DR_mod2, plot = T)
testDispersion(DR_sim2)
check_overdispersion(DR_sim2)
testZeroInflation(DR_sim2)

DR_mod3 <- glmmTMB(No_worms ~ Treatment * poly(Day,2) * Block + (1|Population/Rep), ziformula = ~ Day, data = DR, family = 'nbinom2')
summary(DR_mod3)
DR_sim3 <- simulateResiduals(DR_mod3, plot = T)
testZeroInflation(DR_sim3)
testDispersion(DR_sim3)
check_overdispersion(DR_mod3)

DR_mod4 <- glmmTMB(No_worms ~ Treatment * poly(Day,2) * Block + (1|Population/Rep), ziformula = ~ I(Day^2), data = DR, family = 'nbinom2', control =glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))#won't converge with optimizers or increased interations
summary(DR_mod4)
AIC(DR_mod3, DR_mod4)

DR_sim4 <- simulateResiduals(DR_mod4, plot = T)
testDispersion(DR_sim4)
check_overdispersion(DR_mod4)

AIC(DR_mod, DR_mod2, DR_mod3, DR_mod4)

DR_mod5 <- glmmTMB(No_worms ~ Treatment * poly(Day,2) * Block + (1|Population/Rep), ziformula = ~ Day*Block, data = DR, family = 'nbinom2', control =glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))#won't converge with optimizers or increased interations
summary(DR_mod5)
AIC(DR_mod4, DR_mod5)


DR_mod6 <- glmmTMB(No_worms ~ Treatment * poly(Day,2) * Block + (1|Population/Rep), ziformula = ~ Day + Block, data = DR, family = 'nbinom2')
summary(DR_mod6)
 #####Final model

DR_sim6 <- simulateResiduals(DR_mod6, plot = T)
testDispersion(DR_sim6)
testZeroInflation(DR_sim6)


DR_mod6_2 <- glmmTMB(No_worms ~ Treatment * Day * Block + Treatment *I(Day^2) * Block + (1|Population/Rep), ziformula = ~ Day + Block, data = DR, family = 'nbinom2',, control =glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
Anova(DR_mod6_2, type = 'III')


