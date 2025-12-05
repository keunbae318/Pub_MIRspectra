##PLS model###

##clean##
rm(list=ls()) 
graphics.off()

#library#
library(prospectr)
library(ggplot2)
library(data.table)
library(reshape2)
library(gridExtra)
library(pls)
library(ChemoSpec)
library(caret)
library(ChemoSpec)
library(hyperSpec)
library(tidyverse)
library(ggpubr)
library(spectacles)
library(MLmetrics)
library(nlme)
library(car)
library(multcomp)

#Source
source("Rcode/fun/misc.R")
source("Rcode/fun/C_Spectraaverage.R")

###################################################################################
##Average spectra based on the replicates##
###################################################################################

load("Data/Rawdata/MIRs.RData")


spectraMatrix=as.matrix(dataset$spc)
columns <- colnames(spectraMatrix) %in% c("4000":"400")

subsetMatrix <- spectraMatrix[, columns]
idMatrix=as.matrix(dataset$ID)

Averages <- SpecAver(subsetMatrix, 4, idMatrix, thr = .07)


#A list containing the average spectra's and ID's is constructed
spectraAverageslist           <- list() 
spectraAverageslist$ID        <- Averages$ID                      # The ID of the spectra
spectraAverageslist$Site      <- Averages$Site                    # The sitename
spectraAverageslist$Treatment <- Averages$Treatment               # The replica number of the soil core
spectraAverageslist$Depth     <- Averages$Depth                   # The upper depth of the sample
spectraAverages               <- data.frame(spectraAverageslist)  # The list is converted to a data frame
tmp                           <- data.frame(Averages$SpecMeans)   # The average spectra are copied
spectraAverages$spc           <- tmp                              # The average spectra are added
colnames(spectraAverages$spc) <- colnames(subsetMatrix) # The wavenumbers are used as column names for the spectra
#View(spectraAverages)

remove(list = c("spectraMatrix","tmp","idMatrix","subsetMatrix"))

spectraAverages$Depth=factor(spectraAverages$Depth,levels = c("T","M","B"),
                             labels = c("0-5 cm","15-20 cm","45-50 cm"))
spectraAverages$Treatment=factor(spectraAverages$Treatment,levels = c("U","R","D"),
                                 labels = c("Undrained","Rewetted","Drained"))

#####################################################################################
#GGPlot###
#####################################################################################

WSSL_hy <- new("hyperSpec",spc=as.matrix(spectraAverages$spc),
               wavelength = as.numeric(colnames(spectraAverages$spc)),
               data=spectraAverages[,1:4],
               label=list(.wavelength="nm",spc="Absorbance"))

WSSL_hy

###Plotting for whole data##
ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=Site,group=ID))+
  geom_line()+
  facet_grid(Depth~Treatment)+
  theme(legend.position = "top")+
  theme_bw()+
  #xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('Absorbance')+
  scale_x_reverse() ->Plot
plotly::ggplotly(Plot)

#####################################################################################
###################################################################################
# spectraAverages$spc <- savitzkyGolay(spectraAverages$spc,p=3,w=21,m=0)
# spectraAverages$spc <- Normalization100(spectraAverages$spc)
dataset <- spectraAverages
dataset$spc <- savitzkyGolay(dataset$spc,p=3,w=21,m=0)
source("Rcode/Reference/C_Areabased_Normalization.R")

#center spectra by training mean & apply to test
dataset$Depth
dataset_train <- subset(dataset, Depth %in% c("0-5 cm","45-50 cm"))
dataset_test <- subset(dataset, Depth %in% c("15-20 cm"))


dataset_train_mean <- colMeans(dataset_train$spc)
dataset_train$spc <- sweep(dataset_train$spc, 2, dataset_train_mean, "-")


#dataset4$spc2 <- binning(dataset4$spc2, bin.size = 10)

WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset_train$spc),
               wavelength = as.numeric(colnames(dataset_train$spc)),
               data=dataset_train[,1:4],
               label=list(.wavelength="nm",spc="Absorbance"))

WSSL_hy

###Plotting for whole data##
ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=Site,group=ID))+
  geom_line()+
  facet_grid(Depth~Treatment)+
  theme(legend.position = "top")+
  theme_bw()+
  xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('Absorbance')+
  scale_x_reverse()


#############################################################################################
#############################################################################################
#Labile_C_prediction##
Cummulative.dat=read.csv("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/CO2_Respiration/Rawdata/Cummulativeflux_105.csv")
EA.dat = read.csv("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Elemental_Analysis/Rdata/EA_data.csv")

Cummulative.dat$Flux <- Cummulative.dat$Flux/(EA.dat[EA.dat$Depth %in% c("Top","Bottom"), ]$C*0.01)
#save(Cummulative.dat, file =  "C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/Rcode/Github/PLSmodel_SOC.RData")
load("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/Rcode/Github/PLSmodel_SOC.RData")
########################################################################################
Cummulative.dat$Flux_centered <- Cummulative.dat$Flux - mean(Cummulative.dat$Flux)

Flux_mean <- mean(Cummulative.dat$Flux)


PLSsubset.dat <-dataset_train %>% 
  mutate(flux = Cummulative.dat$Flux_centered) %>% 
  dplyr::select(1:4,6,5) 

########################################################################################
##########################################################################################


group_cv <- function(x, k = length(unique(x))) {
  dat <- data.frame(index = seq(along = x), group = x)
  groups <- data.frame(group = unique(dat$group))
  group_folds <- createFolds(groups$group, returnTrain = TRUE, k = k)
  group_folds <- lapply(group_folds, function(x, y) y[x,,drop = FALSE], y = groups)
  dat_folds <- lapply(group_folds, function(x, y) merge(x, y), y = dat)
  lapply(dat_folds, function(x) sort(x$index))
}

# ##multilevel design#######################################################
# design <- data.frame(
#   Site=PLSsubset.dat$Site
# )
# PLSsubset.dat$spc <- mixOmics::withinVariation(PLSsubset.dat[,"spc"],design)

#PLS model index for LGOCV
plsFit <- train(
  flux ~ spc,
  data = PLSsubset.dat,
  method = "pls",
  tuneLength = 20,
  #preProcess = c("center","scale"),
  trControl = trainControl(method = "cv",
                           index = group_cv(PLSsubset.dat$Site, 13),
                           #number = 12, 
                           #repeats = 1,
                           #summaryFunction = spectroSummary,
                           verboseIter = TRUE,
                           savePredictions =TRUE,
                           #returnData =TRUE,
                           #returnResamp = "all"
  ),
  
)




###################################################################################
#-----------------------------------------------------------------------------------------------------
# One-sided SE paired t-test to select the best number of latent variables
BestTune_selection <- plsFit$pred %>%
  mutate(New_ID = paste0(ncomp, Resample)) %>%
  group_by(New_ID) %>%
  summarise(RMSE = RMSE(obs, pred), ncomp = first(ncomp), Fold =  first(Resample)) %>%
  arrange(ncomp,Fold) %>%
  group_by(ncomp) %>%
  summarise(RMSE_mean = mean(RMSE), RMSE_se = sd(RMSE)/sqrt(n()), ncomp =first(ncomp)) %>%
  mutate(RMSE_mean_plus_se = RMSE_mean+RMSE_se)

Best_nb_latents <- min(which(BestTune_selection$RMSE_mean < BestTune_selection$RMSE_mean_plus_se[plsFit$bestTune[,1]])) 

#PLOT
ggplot(BestTune_selection, aes(x=ncomp, y= RMSE_mean))+ 
  geom_point()+
  geom_line(aes(y=RMSE_mean), linetype = "dashed")+
  geom_point(aes(x=Best_nb_latents, y=RMSE_mean[Best_nb_latents]), color = "red", size = 3)+
  geom_errorbar(aes(ymin=RMSE_mean-RMSE_se , ymax=RMSE_mean+RMSE_se), width=.2,
                position=position_dodge(.9))+
  theme_bw()+
  xlab("Number of components")+
  ylab("RMSE")+
  ggtitle("One-sided paired t-test to select the best number of latent variables (Leave-OneGroup-Out)")+
  theme(legend.position = "top") -> fig1

fig1
#-----------------------------------------------------------------------------------------------------

plsFit2 <- train(
  flux ~ spc,
  data = PLSsubset.dat,
  method = "pls",
  tuneGrid = expand.grid(ncomp = Best_nb_latents),
  #preProcess = c("center","scale"),
  trControl = trainControl(method = "cv",
                           index = group_cv(PLSsubset.dat$Site, 13),
                           #number = 12, 
                           #repeats = 1,
                           #summaryFunction = spectroSummary,
                           verboseIter = TRUE,
                           savePredictions =TRUE,
                           #returnData =TRUE,
                           #returnResamp = "all"
  ),
  
)

plsFit2


###Plotting
Calibration_result  <- predict(plsFit2, newdata = PLSsubset.dat, ncomp = Best_nb_latents)


########################################################################################

stats_result <- data.frame(n=78,
                           nf=5,
                           Rsq_CV = R2_Score(plsFit$pred[plsFit$pred$ncomp == Best_nb_latents,][,1], PLSsubset.dat[,5]),#Rsqured_CV
                           RMSE_CV = RMSE(plsFit$pred[plsFit$pred$ncomp == Best_nb_latents,][,1], PLSsubset.dat[,5]),#RMSE_CV_mean
                           MAE_CV = MAE(plsFit$pred[plsFit$pred$ncomp == Best_nb_latents,][,1], PLSsubset.dat[,5]),
                           Rsq= R2_Score(Calibration_result, PLSsubset.dat[,5]),
                           RMSE= RMSE(Calibration_result, PLSsubset.dat[,5]),
                           MAE=MAE(Calibration_result, PLSsubset.dat[,5])) %>% 
  mutate_if(is.numeric,round,digits=2)


########################################################################################

plot_data <- data.frame(SampleID = row.names(PLSsubset.dat),                #Calibration results
                        LOGO_CV=plsFit$pred[plsFit$pred$ncomp == Best_nb_latents,][,1] + Flux_mean,
                        Calibration=Calibration_result+ Flux_mean, 
                        obs = PLSsubset.dat[,5]+ Flux_mean
) %>% 
  #rename(Calibration = flux.2.comps ) %>% 
  gather(key= "test", value = "pred", 2:3)  


#plot_data <- rbind(plot_data,extracted_data)
plot_data <- subset(plot_data, test == "LOGO_CV")

Pls_plot<-ggplot(aes(obs, pred, group = test, color = test, shape = test), data = plot_data) +
  stat_smooth(aes(obs, pred, group = test, color = test), 
              data = plot_data[plot_data$test %in% c("Calibration","LOGO_CV"),], 
              method='lm', 
              formula = y ~ x,
              se = F)+
  geom_point(size=2)+
  geom_abline(linetype="dashed", color="darkgray",size =1.5)+
  scale_color_manual(values=c("black","grey","red"),breaks=c('LOGO_CV', 'Calibration', 'Validation'))+
  scale_shape_discrete(breaks=c('LOGO_CV', 'Calibration', 'Validation'))+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank())+
  labs(x = expression(Measured~cumulative~flux~on~day~105*phantom(x)*(mu*g~CO[2]-C~g~SOC^{-1})),
       y = expression(Predicted~cumulative~flux~on~day~105*phantom(x)*(mu*g~CO[2]-C~g~SOC^{-1})))+
  geom_text(col = "grey",aes(x = Inf, y = -Inf, label = paste("n =", stats_result$n)), vjust = -9, hjust = 1, size = 5,show.legend = F)+
  geom_text(col = "grey",aes(x = Inf, y = -Inf, label = paste("n_latent =", stats_result$nf)), vjust = -7, hjust = 1, size = 5,show.legend = F)+
  geom_text(col = "grey",aes(x = Inf, y = -Inf, label = paste("Rsq_CV =", sprintf("%.2f", stats_result$Rsq_CV))), vjust = -5, hjust = 1, size = 5,show.legend = F)+
  geom_text(col = "grey",aes(x = Inf, y = -Inf, label = paste("RMSE_CV =", sprintf("%.2f", stats_result$RMSE_CV))), vjust = -3, hjust = 1, size = 5,show.legend = F)+
  geom_text(col = "grey",aes(x = Inf, y = -Inf, label = paste("MAE_CV =", stats_result$MAE_CV)), vjust = -1, hjust = 1, size = 5,show.legend = F)#+
geom_text(col = "grey",aes(x = Inf, y = -Inf, label = paste("Rsq_Cal =", stats_result$Rsq)), vjust = -5, hjust = 1, size = 5,show.legend = F)+
  geom_text(col = "grey",aes(x = Inf, y = -Inf, label = paste("RMSE_Cal =", stats_result$RMSE)), vjust = -3, hjust = 1, size = 5,show.legend = F)+
  geom_text(col = "grey",aes(x = Inf, y = -Inf, label = paste("MAE_Cal =", stats_result$MAE)), vjust = -1, hjust = 1, size = 5,show.legend = F)


Pls_plot

########################################################################################
# #PLS model coefficient regression
# extracted_data <- data.frame(
#   SampleID = row.names(PLSsubset.dat),
#   obs = PLSsubset.dat$flux,
#   pred = predict(plsFit2, newdata = PLSsubset.dat, comps = Best_nb_latents),
#   test = 'Validation') %>% 
#   rename(pred = "flux")

#PLS model coefficient regression
coefficent.dat<-data.frame(plsFit2$finalModel$coefficients)

coefficent2.dat <- coefficent.dat %>% 
  mutate(Wavenumer = rownames(.),
         ##########################################
         ########Change the number of component####
         ##########################################
         Coefficient = coefficent.dat$`.outcome.2.comps`) %>%  ########Change the number of component
  dplyr::select(Wavenumer,Coefficient) %>% 
  mutate(Wavenumer = as.numeric(gsub("spc","",Wavenumer)))




loading.spc <- as_tibble(coefficent2.dat[,1:2]) 

loading.spc$Wavenumer <- as.numeric(loading.spc$Wavenumer)
loading.spc$axis <- "Regression coefficient"
# loading.spc$PC_axis <- factor(loading.spc$PC_axis, levels = c("comp1","comp2"),
#                               labels = c("Latent variable 1","Latent variable 2"))


PLS_coefficients_Plot<-ggplot(loading.spc,aes(x=Wavenumer,y=Coefficient, color = axis))+ 
  geom_line(size=0.8)+
  theme_bw()+
  theme(legend.position = "top",
        legend.title =element_blank(),
        panel.grid.major.x = element_line(color = "black",
                                          size = 0.1,
                                          linetype = 2),
        panel.grid.major.y = element_line(color = "grey",
                                          size = 0.1,
                                          linetype = 2))+
  xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('PLS regression coefficients')+
  scale_x_reverse(n.breaks = 10)+
  grids(linetype = "dashed")+
  scale_color_manual(values=c("black"))


PLS_coefficients_Plot





ggarrange(Pls_plot, PLS_coefficients_Plot, ncol = 2, labels = c("A", "B"), widths = c(1,1.5))



Inspection.dat<-ggplot(as.data.frame(loading.spc),aes(x=Wavenumer,y=Coefficient, color = axis))+ 
  geom_line(size=0.8)+
  theme_bw()+
  scale_x_reverse(n.breaks = 10)+
  grids(linetype = "dashed")+
  scale_color_manual(values=c("black"))

library(plotly)
ggplotly(Inspection.dat)
Inspection.dat


############################################################################################
#Predict other layer with established model
############################################################################################

PLSsubset.dat <- dataset_train%>% 
  mutate(flux = Cummulative.dat$Flux_centered) %>% 
  dplyr::select(1:4,6,5)

# ##multilevel design#######################################################
# design <- data.frame(
#   Site=PLSsubset.dat$Site
# )
# PLSsubset.dat$spc <- mixOmics::withinVariation(PLSsubset.dat[,"spc"],design)

#PLS model index for LGOCV
group_cv(PLSsubset.dat$Site, 13)
plsFit <- train(
  flux ~ spc,
  data = PLSsubset.dat,
  method = "pls",
  tuneLength = 20,
  #preProcess = c("center","scale"),
  trControl = trainControl(method = "cv",
                           index = group_cv(PLSsubset.dat$Site, 13),
                           #number = 12, 
                           #repeats = 1,
                           #summaryFunction = spectroSummary,
                           verboseIter = TRUE,
                           savePredictions =TRUE,
                           #returnData =TRUE,
                           #returnResamp = "all"
  ),
  
)


###################################################################################
#-----------------------------------------------------------------------------------------------------
# One-sided SE paired t-test to select the best number of latent variables
BestTune_selection <- plsFit$pred %>%
  mutate(New_ID = paste0(ncomp, Resample)) %>%
  group_by(New_ID) %>%
  summarise(RMSE = RMSE(obs, pred), ncomp = first(ncomp), Fold =  first(Resample)) %>%
  arrange(ncomp,Fold) %>%
  group_by(ncomp) %>%
  summarise(RMSE_mean = mean(RMSE), RMSE_se = sd(RMSE)/sqrt(n()), ncomp =first(ncomp)) %>%
  mutate(RMSE_mean_plus_se = RMSE_mean+RMSE_se)

Best_nb_latents <- min(which(BestTune_selection$RMSE_mean < BestTune_selection$RMSE_mean_plus_se[plsFit$bestTune[,1]])) 

#PLOT
ggplot(BestTune_selection, aes(x=ncomp, y= RMSE_mean))+ 
  geom_point()+
  geom_line(aes(y=RMSE_mean), linetype = "dashed")+
  geom_point(aes(x=Best_nb_latents, y=RMSE_mean[Best_nb_latents]), color = "red", size = 3)+
  geom_errorbar(aes(ymin=RMSE_mean-RMSE_se , ymax=RMSE_mean+RMSE_se), width=.2,
                position=position_dodge(.9))+
  theme_bw()+
  xlab("Number of components")+
  ylab("RMSE")+
  ggtitle("One-sided paired t-test to select the best number of latent variables (Leave-OneGroup-Out)")+
  theme(legend.position = "top") -> fig1

fig1
#-----------------------------------------------------------------------------------------------------

plsFit2 <- train(
  flux ~ spc,
  data = PLSsubset.dat,
  method = "pls",
  tuneLength = Best_nb_latents,
  #preProcess = c("center","scale"),
  trControl = trainControl(method = "cv",
                           index = group_cv(PLSsubset.dat$Site, 13),
                           #number = 12, 
                           #repeats = 1,
                           #summaryFunction = spectroSummary,
                           verboseIter = TRUE,
                           savePredictions =TRUE,
                           #returnData =TRUE,
                           #returnResamp = "all"
  ),
  
)

plsFit2

##Predit the second layer##

dataset_train_mean <- colMeans(dataset_train$spc)
dataset_test$spc <- sweep(dataset_test$spc, 2, dataset_train_mean, "-")

Second_layer.dat <- dataset_test

Predicted_B_layer<- predict(plsFit2, newdata = Second_layer.dat, ncomp = Best_nb_latents)

Twolayer_measured.data <- data.frame(SampleID = row.names(PLSsubset.dat),
                                     Flux = PLSsubset.dat$flux + Flux_mean
)

Second_layer_predicted <- data.frame(SampleID = row.names(Second_layer.dat),
                                     Flux = Predicted_B_layer + Flux_mean
)


estimated_flux <- rbind(Twolayer_measured.data,Second_layer_predicted)
meta.dat <- dataset[,1:4] %>% 
  mutate(SampleID = row.names(spectraAverages))


plot_data <-left_join(meta.dat, estimated_flux, by = "SampleID") 
plot_data$Treatment <-factor(plot_data$Treatment, level = c("Undrained","Drained","Rewetted"))
plot_data$Depth <-factor(plot_data$Depth, labels = c("0-5 cm (measured)",
                                                     "15-20 cm (predicted)",
                                                     "45-50 cm (measured)"))

ggboxplot(plot_data, x = "Treatment", y = "Flux", width = 0.8, facet.by = "Depth",fill = "Treatment")+
  scale_fill_manual(values=c("#00AFBB","#FC4E07","#E7B800"))+
  theme_bw()+
  theme(legend.position="top",
        legend.title = element_blank())+
  xlab("")+
  ylab(expression(Cumulative~flux~on~day~105*phantom(x)*(mu*g~CO[2]-C~g~SOC^{-1})))

#saveRDS(plot_data,"C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Intergration/R_data/flux.rds")
############################################################################################
############################################################################################
######################POM C##########################                 
ANOVA.dat <- plot_data

ANOVA.dat <- subset(ANOVA.dat, Depth %in% c("0-5 cm (measured)","45-50 cm (measured)"))
ANOVA.dat <- subset(ANOVA.dat, Depth == "15-20 cm (predicted)" )

ANOVA.dat$Flux

ANOVA.dat %>%
  group_by(Depth, Treatment) %>%
  summarise(Flux = format(round(mean(Flux), 2), nsmall = 2))

max(ANOVA.dat$Flux)
min(ANOVA.dat$Flux)

ANOVA.lme=lme(Flux~Treatment*Depth,random =~1|Site,na.omit(ANOVA.dat))
#####################################################
##Assumption##
plot(ANOVA.lme) #homoscedasticity
qqnorm(ANOVA.lme,~resid(.,type="p")|Site,abline = c(0,1)) #normality of resid
qqnorm(ANOVA.lme,~ranef(.,standard=T),abline=c(0,1)) # ind. of random effects
plot(ANOVA.lme,Site~resid(.)) #normality of ranef
shapiro.test(resid(ANOVA.lme))
shapiro.test(ranef(ANOVA.lme)[,1])
#####################################################
summary(ANOVA.lme)
anova(ANOVA.lme)

#If not interaction
########################################################
ANOVA.lme=lme(Flux~Treatment+Depth,random =~1|Site,na.omit(ANOVA.dat))

summary(ANOVA.lme)
anova(ANOVA.lme)

cld(lsmeans(ANOVA.lme,~Treatment))
cld(lsmeans(ANOVA.lme,~Depth))
###If assumptions are violated, the data need to be transformd with Boxcoxfit##
ANOVA.dat$Flux <- ANOVA.dat$Flux - min(ANOVA.dat$Flux) +1

boxcoxfit(na.omit(ANOVA.dat$Flux))
ANOVA.bcx=bcPower(ANOVA.dat$Flux,lambda =  0.2672839)


ANOVA.dat$ANOVA.bcx=ANOVA.bcx
ANOVA.lme=lme(ANOVA.bcx~Treatment*Depth,random =~1|Site,na.omit(ANOVA.dat))
#Output#######################
#############################
summary(ANOVA.lme)
anova(ANOVA.lme)
#if interaction not significant###
ANOVA.lme=lme(ANOVA.bcx~Treatment+Depth,random =~1|Site,na.omit(ANOVA.dat))


summary(ANOVA.lme)
anova(ANOVA.lme)

cld(lsmeans(ANOVA.lme, ~Treatment|Depth))
cld(lsmeans(ANOVA.lme, ~Treatment|Depth))
cld(lsmeans(ANOVA.lme,~Treatment))
cld(lsmeans(ANOVA.lme,~Depth))
#########################################################

Table1=cld(lsmeans(ANOVA.lme,~Treatment))
Table2=cld(lsmeans(ANOVA.lme,~Depth))

Table_C=as.data.frame(Table1)
Table_F=as.data.frame(Table2)

Table_C <- Table_C[order(Table_C$Treatment), ]
Table_F <- Table_F[order(Table_F$Depth), ]

names(Table_C)[names(Table_C)=="Treatment"] <-"Treatment"
names(Table_F)[names(Table_F)=="Depth"] <-"Treatment"

Table=rbind(Table_C,Table_F)
Table$C_value=10^(log10(Table$lsmean* 0.2672839 +1)/ 0.2672839)
#write.csv(Table,"R:/Guillermo's Lab/POM/Research/POM/POM Data/Data statical analysis/Sequence_POM/Re_analysis_table/Table.csv")
Table


