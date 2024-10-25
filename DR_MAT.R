
#Setup

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

