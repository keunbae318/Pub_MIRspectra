
# The workspace is cleared
rm(list=ls()) # The workspace is cleared
graphics.off()
#initial variable
getwd()
setwd("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/")
Sys.setenv(LANGUAGE="EN")

#source("C_Areabased_Normalization.R")
load("Data/Rawdata/MIRs.RData")
#source("fun/misc.R")
source("Rcode/fun/C_Spectraaverage.R")
source("Rcode/fun/C_Spectraaverage.R")
library(tidyverse)
library(hyperSpec)
library(prospectr)
library(caret)
library(wesanderson)
library(ggplot2)
library(plotly)
library(pROC)
library(doParallel)
library(ggpubr)
#################################################################
##################################################################
# library(stats)
# library(dplyr)
# library(grDevices)
# library(ggpubr)
# library(reshape2)
# library(gridExtra)
# library(ChemoSpec)
# library(ChemoSpecUtils)
# library(ggrepel)
# library(factoextra)
# library(tidyr)
# library(mvoutlier)
# library(plotly)
# library(tidyr)

# source the file containing other functions

#------------------------------------------------------------------
#Calibratio average spectra are caluacated based on the replicates
#------------------------------------------------------------------
load("Data/Rawdata/MIRs.RData")
spectraMatrix=as.matrix(dataset$spc)
columns <- colnames(spectraMatrix) %in% c("4000":"400")
subsetMatrix <- spectraMatrix[, columns]
idMatrix=as.matrix(dataset$ID)

Averages <- SpecAver(subsetMatrix,4,idMatrix, thr = .07)


#A list containing the average spectra's and ID's is constructed
spectraAverageslist           <- list() 
spectraAverageslist$ID        <- Averages$ID                      # The ID of the spectra
spectraAverageslist$Site      <- Averages$Site                    # The sitename
spectraAverageslist$Treatment <- Averages$Treatment                  # The replica number of the soil core
spectraAverageslist$Depth     <- Averages$Depth                   # The upper depth of the sample
spectraAverages               <- data.frame(spectraAverageslist)  # The list is converted to a data frame
tmp                           <- data.frame(Averages$SpecMeans)   # The average spectra are copied
spectraAverages$spc           <- tmp                              # The average spectra are added
colnames(spectraAverages$spc) <- colnames(subsetMatrix) # The wavenumbers are used as column names for the spectra

rm(spectraMatrix,tmp)

#Data reconstruct by average

dataset2<-spectraAverages

EA.dat=read.csv("Data/EA_data.csv")
Microbial.dat=read.csv("Data/Microbial_data.csv")
Variables.dat<-left_join(EA.dat, Microbial.dat)
dataset4 <-left_join(Variables.dat,dataset2)


dataset4$Depth=factor(dataset4$Depth,levels = c("T","M","B"),
                      labels = c("0-5 cm","15-20 cm","45-50 cm"))
dataset4$Treatment=factor(dataset4$Treatment,levels = c("U","D","R"),
                          labels = c("Undrained","Drained","Rewetted"))
dataset4$TsR=factor(dataset4$TsR,levels = c("L","M","S"),
                    labels = c(">25","10-25","<10"))

WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc),
               wavelength = as.numeric(colnames(dataset4$spc)),
               data=dataset4[,1:11],
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

########################################################################
#Simple processioning - Noise removal######
#########################################################################
dataset4$spc2 <- savitzkyGolay(dataset4$spc,p=3,w=21,m=0)
source("Rcode/Reference/C_Areabased_Normalization.R")

#center spectra by training mean & apply to test
datset4_spc2_mean <- colMeans(dataset4$spc2)
dataset4$spc2 <- sweep(dataset4$spc2, 2, datset4_spc2_mean, "-")

WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc2),
               wavelength = as.numeric(colnames(dataset4$spc2)),
               data=dataset4[,1:11],
               label=list(.wavelength="nm",spc="Absorbance"))

###Plotting for whole data##
ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=Site,group=ID))+
  geom_line()+
  facet_grid(Depth~Treatment)+
  theme(legend.position = "top")+
  theme_bw()+
  xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('Absorbance')+
  scale_x_reverse()


# save(dataset4, file = "Rcode/Github/Data_PLSDAmodel.RData")
# load("Rcode/Github/Data_PLSDAmodel.RData")
# We train the PLS algorithm on the training set
source("Rcode/fun/C_Group_CV.R")

Depth_factor <- c("Full Depth","0-5 cm","15-20 cm", "45-50 cm")
Depth_ID <- c("Overall","Top","Middle","Bottom")
#### Data storage####
latent_results <- list()   # per-depth tuning results
coef_results   <- list()   # per-depth regression coefficients
class_metrics_results <- list()
cm_full_matrics <- list()

for (i in 1: length(Depth_factor)) {
  
  
  depth_name <- Depth_factor[i]
  
  if(depth_name == "Full Depth") {
    
    subset.dat <- subset(dataset4)# full depth
  } else {
    subset.dat <- subset(dataset4, Depth == depth_name)
  }
  #multilevel design
  design <- data.frame(Site=subset.dat$Site)
  subset.dat$spc2 <- mixOmics::withinVariation(subset.dat[,"spc2"],design)
  
  
  plsFit <- caret::train(Treatment~ spc2, data = subset.dat,
                         method = "pls",
                         tuneLength = 20,
                         trControl=trainControl(method = "cv",
                                                index = group_cv(subset.dat$Site,13),
                                                #classProbs = TRUE,
                                                verboseIter = TRUE,
                                                savePredictions =TRUE,
                                                summaryFunction = multiClassSummary
                                                ),
                         metric = "Accuracy"
                         )

  #Best_ncomp
  BestTune_selection <- plsFit$pred %>%
    mutate(New_ID = paste0(ncomp, Resample)) %>%
    group_by(New_ID) %>%
    summarise(Accuracy = mean(pred == obs), ncomp = first(ncomp), Fold =  first(Resample)) %>%
    arrange(ncomp,Fold) %>%
    group_by(ncomp) %>%
    summarise(Accuracy_mean = mean(Accuracy), Accuracy_se = sd(Accuracy)/sqrt(n()), ncomp =first(ncomp)) %>%
    mutate(Accuracy_mean_minus_se = Accuracy_mean-Accuracy_se)
  
    Best_nb_latents <- min(which(BestTune_selection$Accuracy_mean_minus_se[plsFit$bestTune[,1]]<= BestTune_selection$Accuracy_mean))


  plsFit_ncomp <- caret::train(Treatment~ spc2, data = subset.dat,
                         method = "pls",
                         tuneGrid = expand.grid(ncomp = as.numeric(Best_nb_latents)), 
                         trControl=trainControl(method = "cv",
                                                index = group_cv(subset.dat$Site,13),
                                                classProbs = TRUE,
                                                verboseIter = TRUE,
                                                savePredictions =TRUE,
                                                summaryFunction = multiClassSummary),
                         metric = "Accuracy"
  )

  
  
  #Metadata ------------------------------------------------------------------------------
  BestTune_selection$model <- depth_name
  BestTune_selection$Best_nb_latents <- Best_nb_latents

  #save_Result
  latent_results[[depth_name]] <- BestTune_selection
  
  #Metadata for coefficient--------------------------------------------------------------
  coefficent.dat <- data.frame(coef(plsFit_ncomp$finalModel))
  
  coefficent2.dat <- coefficent.dat %>% 
    mutate(Wavenumber = rownames(.)) %>% 
    rename(Undrained = names(.)[1],
           Drained = names(.)[2],
           Rewetted = names(.)[3],) %>%
    tidyr::gather(key = "Treatment", value = Coefficient, -Wavenumber) %>%
    mutate(Wavenumber = as.numeric(gsub("spc2", "", Wavenumber)),
           Treatment = factor(Treatment, levels = c("Undrained", "Drained","Rewetted")))
  
  loading.spc <- as_tibble(coefficent2.dat[,1:3]) 
  loading.spc$Wavenumber <- as.numeric(loading.spc$Wavenumber)
  loading.spc$axis <- "Regression coefficient"
  
  loading2.spc <- loading.spc %>%
    mutate(Depth=as.character(depth_name))
  
  #save_Result
  coef_results[[depth_name]] <- loading2.spc
  
  #Confusion matrix at best ncomp--------------------------------------------------------
  cm_full  <- confusionMatrix(plsFit_ncomp$pred[,1],plsFit_ncomp$pred[,2]) 
  cm_full_matrics[[depth_name]] <- cm_full$table
  # Sensitivity = Recall
  recall <- cm_full$byClass[,"Sensitivity"]
  precision <- cm_full$byClass[,"Pos Pred Value"]
  F1 <- 2 * (precision * recall) / (precision + recall)
  F1[is.na(F1)] <- 0
  # Combine into a data frame
  class_metrics <- data.frame(
    Class = rownames(cm_full$byClass),
    Precision = precision,
    Recall = recall,
    F1 = F1
  )
  class_metrics$macro_F1 <-mean(class_metrics$F1)
  
  class_metrics_results[[depth_name]] <- class_metrics
}


#Latent variable selection plot##################################################################
latent_variable_selection <- bind_rows(latent_results)
latent_variable_selection$model <-factor(latent_variable_selection$model, levels = c("Full Depth","0-5 cm","15-20 cm","45-50 cm"))


ggplot(latent_variable_selection,aes(x= ncomp, y= Accuracy_mean, group = model,color = model))+
  geom_point(aes(shape=model))+
  geom_point(data = latent_variable_selection[latent_variable_selection$ncomp == latent_variable_selection$Best_nb_latents,],
             aes(x=ncomp, y=Accuracy_mean), color= "red", size=3)+
  geom_errorbar(aes(ymin = Accuracy_mean - Accuracy_se, ymax = Accuracy_mean + Accuracy_se), width = 0.1)+
  geom_line()+
  geom_vline(xintercept = c(1,2,6,11), linetype = "dashed")+
  xlab("Number of Latents")+
  ylab("Accuracy")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title =element_blank())+
  scale_color_manual(values=wes_palette(n=4, name="GrandBudapest1")) -> latent.graph

Best_ncomp_result <- latent_variable_selection[latent_variable_selection$ncomp == latent_variable_selection$Best_nb_latents,]

#Coefficient plot#########################################################################################
coef_result_all <- bind_rows(coef_results)

coef_result_all$Depth
coef_result_all <- coef_result_all %>% 
  mutate(Depth = factor(Depth, levels= c("Full Depth", "0-5 cm", "15-20 cm", "45-50 cm")),
         Treatment = factor(Treatment, levels = c("Undrained", "Drained","Rewetted"))) %>% 
  subset(Depth != "Full Depth")


#saveRDS(coef_result_all, "Data/PLSDA/coef_result_all.RDS")
coef_result_all<- readRDS("Data/PLSDA/coef_result_all.RDS")

fig3<-
  ggplot(coef_result_all,aes(x=Wavenumber,y=Coefficient, color = Treatment ))+ 
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
  ylab('PLSDA coefficients')+
  scale_x_reverse(n.breaks = 10)+
  facet_grid(Depth~Treatment, scales = "free_y")+
  #grids(linetype = "dashed")+
  scale_color_manual(values=c("#00AFBB","#FC4E07","#E7B800")) 


fig3
ggplotly(fig3)

ggarrange(fig1,fig3,ncol=2,nrow=1,labels = c("A", "B"), widths = c(1,1.3))

#################################################################################
#Confusion matrix report # Table for S2
##############################
cm_full_matrics
class_metrics_results
##############################################################################33#
#################################################################################
library(pROC)


list_names <- c("Undrained","Drained","Rewetted")


for (i in seq_along(Depth_factor)) {
  depth_name <- Depth_factor[i]
  
  #subset data for depth
    if(depth_name == "Full Depth") {
    
    subset.dat <- subset(dataset4)# full depth
  } else {
    subset.dat <- subset(dataset4, Depth == depth_name)
  }
  
  #multilevel design
  design <- data.frame(Site=subset.dat$Site)
  subset.dat$spc2 <- mixOmics::withinVariation(subset.dat[,"spc2"],design)
  
  Best_nb_latents <- Best_ncomp_result$Best_nb_latents[Best_ncomp_result$model == depth_name]
  
  plsFit_pROC <- caret::train(Treatment~ spc2, data = subset.dat,
                               method = "pls",
                               tuneGrid = expand.grid(ncomp = as.numeric(Best_nb_latents)), 
                               trControl=trainControl(method = "cv",
                                                      index = group_cv(subset.dat$Site,13),
                                                      classProbs = TRUE,
                                                      verboseIter = TRUE,
                                                      savePredictions =TRUE,
                                                      summaryFunction = multiClassSummary),
                               metric = "Accuracy"
  )
  
  cl <- makeCluster(6)
  registerDoParallel(cl)
  
  ROC_HR<-foreach (i = 1:3, .combine = 'c', .packages = "pROC") %dopar% {
    
    HR <- list_names[i]
    
    #data structure#
    ROC.dat <- plsFit_pROC$pred[plsFit_pROC$pred$ncomp == Best_nb_latents,]
    ROC.dat$obs <- ifelse(ROC.dat$obs == HR, TRUE, FALSE)
    
    #ROC analysis
    formula <- as.formula(paste0("obs~", HR))
    HR.ROC <- roc(formula,data = ROC.dat,
                  #method ="bootstrap",
                  smooth=TRUE,
                  smooth.method="density",
                  # arguments for ci
                  ci=TRUE, 
                  boot.n=10000,
                  #boot.n=100,
                  boot.stratified=TRUE
    )
    
    setNames(list(HR.ROC), HR)
  }
  
  stopCluster(cl)
  registerDoSEQ()
  
  #==============================================================
  # 6) CI ribbons for ROC curves
  #==============================================================
  cl <- makeCluster(6)
  registerDoParallel(cl)
  
  ci.list<-foreach (i = 1:3, .combine = 'c', .packages = "pROC") %dopar% {
    
    HR <- names(ROC_HR)[i]
    
    
    #ROC.ci.se analysis
    ROC.ci.se <- ci.se(ROC_HR[[i]], 
                       specificities = seq(0, 1, .01),
                       conf.level=0.95, 
                       boot.n=10000,
                       #boot.n=100,
                       boot.stratified=TRUE, 
                       parallel=TRUE)
    
    setNames(list(ROC.ci.se), HR)
  }
  
  stopCluster(cl)
  registerDoSEQ()
  
  
  #==============================================================
  # 7) Save RDS output
  #==============================================================
  saveRDS(ROC_HR, file = paste0("Data/ROC/",depth_name,"_ROC_HR.Rds"))
  saveRDS(ci.list, file = paste0("Data/ROC/",depth_name,"_ci.list.Rds"))
  #==============================================================
  # 8) Save as environment variables (optional)
  #==============================================================
  
  assign(paste0("ROC_HR_", Depth_ID[i]), ROC_HR)
  assign(paste0("ci.list_", Depth_ID[i]), ci.list)
}

ROC_HR_Overall
ci.list_Overall
##################################################################################
#Graph#######
##################################################################################
ROC_HR_Overall<- readRDS("Data/ROC/Overall_ROC_HR.Rds")
ci.list_Overall <- readRDS("Data/ROC/Overall_ci.list.Rds")
ROC_HR_Top<- readRDS("Data/ROC/Top_ROC_HR.Rds")
ci.list_Top <- readRDS("Data/ROC/Top_ci.list.Rds")
ROC_HR_Middle<- readRDS("Data/ROC/Middle_ROC_HR.Rds")
ci.list_Middle <- readRDS("Data/ROC/Middle_ci.list.Rds")
ROC_HR_Bottom<- readRDS("Data/ROC/Bottom_ROC_HR.Rds")
ci.list_Bottom <- readRDS("Data/ROC/Bottom_ci.list.Rds")

ROC_HR_Overall$Undrained$auc
as.numeric(ROC_HR_Overall$Undrained$ci)

for(a in 1:4){
  
  ROC_HR<- get(paste0("ROC_HR_", Depth_ID[a]))
  ci.list <- get(paste0("ci.list_", Depth_ID[a]))
  
  p <-
    ggroc(ROC_HR, aes=c("linetype", "color"), size =1)+
    geom_abline(intercept = 1, slope = 1, color='grey',size = 0.5,linetype = "dashed") + 
    labs(x = "1-Specificity", y = "Sensitivity") + 
    scale_x_reverse(expand = c(0.01, 0.001),breaks= seq(from = 0 , to = 1 ,by = 0.25),limits = c(1,-0.01),
                    labels = c("1.00","0.75","0.5","0.25","0.0"))+
    theme_bw()+
    theme(legend.position = "top",
          legend.title =element_blank())+
    scale_color_manual(values=c("#00AFBB","#FC4E07","#E7B800"))
  p
  ############################################################################
  dat.ci.list <- lapply(ci.list, function(ciobj) 
    data.frame(x = as.numeric(rownames(ciobj)),
               lower = ciobj[, 1],
               upper = ciobj[, 3]))
  ############################################################################
  # Define a vector of colors
  fill_colors <- c("#00AFBB", "#FC4E07","#E7B800")
  for(i in 1:3) {
    p <- p + geom_ribbon(
      data = dat.ci.list[[i]],
      aes(x = x, ymin = lower, ymax = upper),
      fill = fill_colors[i],  # Set fill color from the predefined colors
      alpha = 0.2,
      inherit.aes = F
    )+
      annotate("text", 
               x = 0.3, 
               y = 0.1 - (i - 1) * 0.1,
               label = paste0(
                 sprintf("%.2f", ROC_HR[[i]]$auc),  # Always 2 digits
                 " (95 % CI: ",
                 sprintf("%.2f", as.numeric(ROC_HR[[i]]$ci)[1]), "-",
                 sprintf("%.2f", as.numeric(ROC_HR[[i]]$ci)[3]),
                 ")"
               ),
               color = fill_colors[i],
               size = 4)
    
  }
  p
  assign(paste0("Plot_", Depth_ID[a]), p)
} 

Plot_Overall
Plot_Top
Plot_Middle
Plot_Bottom

ggarrange(Plot_Overall,Plot_Top,Plot_Middle,Plot_Bottom,legend = "top",labels = c("  Full Depth","    0-5 cm", "   15-20 cm","   45-50 cm"),
          common.legend = TRUE,ncol=2,nrow=2) -> ROC.graph

ggarrange(latent.graph,ROC.graph,ncol=2,nrow=1,labels = c("A", "B"), widths = c(1,1.3))

