##------------------------------------------------
##OSC helper functions (fit on training, apply to new data)
##-----------------------------------------------
fit_osc <- function(X, Y, ncomp = 1, tol = 1e-10, max_iter = 100, center = TRUE) {
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  if(center) {
    X_center <- colMeans(X)
    Y_center <- colMeans(Y)
    X <- scale(X, center = X_center, scale = FALSE)
    Y <- scale(Y, center = Y_center, scale = FALSE)
  } else {
    X_center <- rep(0, ncol(X))
    Y_center <- rep(0, ncol(Y))
  }
  W_list <- list(); P_list <- list()
  Xcorr <- X
  for(i in 1:ncomp) {
    w <- t(Xcorr) %*% Y
    if(all(w == 0)) break
    w <- w / sqrt(sum(w^2))
    tscore <- Xcorr %*% w
    dif <- 1; iter <- 0
    while(dif > tol && iter < max_iter) {
      # orthogonalize tscore to Y (fast scalar form)
      coeff <- drop(t(Y) %*% tscore) / drop(t(Y) %*% Y)
      t_new <- tscore - Y * coeff
      # update weight
      denom <- drop(t(t_new) %*% t_new)
      if(denom == 0) break
      w <- t(Xcorr) %*% t_new / denom
      w <- w / sqrt(sum(w^2))
      t_new <- Xcorr %*% w
      dif <- sqrt(sum((t_new - tscore)^2) / sum(t_new^2))
      tscore <- t_new
      iter <- iter + 1
    }
    p <- as.vector(t(Xcorr) %*% tscore / drop(t(tscore) %*% tscore))
    Xcorr <- Xcorr - tscore %*% t(p)
    W_list[[i]] <- as.vector(w)
    P_list[[i]] <- as.vector(p)
  }
  return(list(Xcorr = Xcorr, W = W_list, P = P_list,
              X_center = X_center, Y_center = Y_center))
}

predict_osc <- function(X_new, osc_model) {
  X_new <- as.matrix(X_new)
  X_new <- sweep(X_new, 2, osc_model$X_center, FUN = "-")
  for(i in seq_along(osc_model$W)) {
    w <- osc_model$W[[i]]
    p <- osc_model$P[[i]]
    tscore <- X_new %*% w
    X_new <- X_new - tscore %*% t(p)
  }
  return(X_new)
}

##------------------------------------------------
## 2) Preprocessing helper
##    returns list(train_spc2, test_spc2, preprocName)
##------------------------------------------------
apply_preproc <- function(p_idx, Train_spc, Test_spc) {
  # Train_spc, Test_spc are numeric matrices (samples x wavelengths)
  if(p_idx == 1){
    Xtr <- movav(Train_spc, 21)
    Xte <- movav(Test_spc, 21)
    name <- "movav"
  } else if(p_idx == 2){
    Xtr <- t(diff(t(Train_spc), differences = 1))
    Xte <- t(diff(t(Test_spc), differences = 1))
    name <- "firstDer"
  } else if(p_idx == 3){
    Xtr <- t(diff(t(Train_spc), differences = 2))
    Xte <- t(diff(t(Test_spc), differences = 2))
    name <- "secondDer"
  } else if(p_idx == 4){
    Xtr <- savitzkyGolay(Train_spc, p=2, w=21, m=0)
    Xte <- savitzkyGolay(Test_spc, p=2, w=21, m=0)
    name <- "sg"
  } else if(p_idx == 5){
    Xtr <- savitzkyGolay(Train_spc, p=2, w=21, m=1)
    Xte <- savitzkyGolay(Test_spc, p=2, w=21, m=1)
    name <- "sgFirstDer"
  } else if(p_idx == 6){
    Xtr <- savitzkyGolay(Train_spc, p=2, w=21, m=2)
    Xte <- savitzkyGolay(Test_spc, p=2, w=21, m=2)
    name <- "sgSecondDer"
  } else if(p_idx == 7){
    Xtr <- msc(as.matrix(Train_spc))
    Xte <- msc(as.matrix(Test_spc))
    name <- "msc"
  } else if(p_idx == 8){
    Xtr <- standardNormalVariate(Train_spc)
    Xte <- standardNormalVariate(Test_spc)
    name <- "snv"
  } else if(p_idx == 9){
    wav_tr <- as.numeric(colnames(Train_spc))
    wav_te <- as.numeric(colnames(Test_spc))
    Xtr <- detrend(Train_spc, wav = wav_tr)
    Xte <- detrend(Test_spc, wav = wav_te)
    name <- "snv_dt"
  } else {
    Xtr <- Train_spc
    Xte <- Test_spc
    name <- "no_prepro"
  }
  return(list(Xtr = as.matrix(Xtr), Xte = as.matrix(Xte), name = name))
}
##########################################################################################################
##clean##
rm(list=ls()) 
graphics.off()

Path <- "C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/"
setwd("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/")


#PLS model##
#library
library(tools)
library(ggpubr)
library(caret)
library(tidyverse)
library(MLmetrics)
library(spectacles)
library(rstatix)
library(prospectr)
library(plsVarSel)
library(mdatools)
library(parallel)
library(doParallel)
#library(doMC)
library(signal)
library(mt)
library(prospectr)
library(ggplot2)
library(data.table)
library(reshape2)
library(gridExtra)
library(pls)
library(ChemoSpec)
library(hyperSpec)


source("Rcode/fun/C_Spectraaverage.R")
source("Rcode/fun/C_Group_CV.R")
###################################################################################
##Average spectra based on the replicates##
###################################################################################
load("Data/Rawdata/MIRs.RData")


spectraMatrix=as.matrix(dataset$spc)
#columns <- colnames(spectraMatrix) %in% c("4000":"650")
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
Spectra.dat <-spectraAverages

#####################################################################################
#GGPlot###
#####################################################################################
Spectra.dat 
#Spectra.dat$spc <- scale(Spectra.dat$spc, center = T, scale =F)
#Spectra.dat$spc <- savitzkyGolay(Spectra.dat$spc,p=3,w=21,m=0) 
#?savitzkyGolay
WSSL_hy <- new("hyperSpec",spc=as.matrix(Spectra.dat$spc),
               wavelength = as.numeric(colnames(Spectra.dat$spc)),
               data=spectraAverages[,1:4],
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

#####################################################################################
#Load the GCMS data##
GCMS.dat<-read.csv("Data/Proportion_pyrolysis_V2.csv")
GCMS.dat$Sum.aliphatic1 <- GCMS.dat$Aliphatics + GCMS.dat$Alkanes + GCMS.dat$Alkenes
GCMS.dat$Sum.aliphatic2 <- GCMS.dat$Aliphatics + GCMS.dat$Alkanes + GCMS.dat$Alkenes +GCMS.dat$Methyl.ketones
str(GCMS.dat)

Compounds_group = c("Aliphatics",
                    "Alkanes",
                    "Alkenes",
                    "Benzenes",
                    "Carbohydrates",
                    "Lignins",       
                    "Methyl.ketones",
                    "N.compounds",
                    "Phenols",
                    "Unidentified",
                    "Sum.aliphatic1",
                    "Sum.aliphatic2")

Spectra2.dat <- left_join(Spectra.dat, GCMS.dat[,c(1,7:16,22,23)], by = "ID") %>% 
  dplyr::select(c(1:4,6:17,5))

getwd()
#save(Spectra2.dat, file = "Data/PLS_model/Spectra2.RData")
#save(Spectra2.dat, file = "Rcode/Github/PLSmodel_SOM_pyGCMS.RData")

#load("Data/PLS_model/Spectra2.RData")


###################################################################################
#Partioning data from mid-infrared spectrascopy
###################################################################################
Spectra2.dat %>%
  dplyr::select(Treatment, Depth) %>%
  group_by(Treatment, Depth) %>%
  group_indices() -> indeces

set.seed(123)
P_order <- createDataPartition(as.factor(indeces), p = 0.8, list = FALSE)
Train_data <- Spectra2.dat[P_order,]; Test_data <- Spectra2.dat[-P_order,]

Org_Train_data <- Train_data; Org_Test_data <- Test_data

##------------------------------------------------
## 3)parallel setup
##------------------------------------------------
cores <- parallel::detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

cl <- makeCluster(12)
registerDoParallel(cl)

train_ctrl_template <- trainControl(method = "cv",
                                    index = group_cv(Org_Train_data$Site, 13),
                                    savePredictions = TRUE,
                                    allowParallel = FALSE,
                                    verboseIter =FALSE)


##------------------------------------------------
## 4) Prepare orginal copies of targets/spectra (do not modify orignals)
##-----------------------------------------------

Org_Train_data <- Train_data
Org_Test_data <- Test_data


##------------------------------------------------
## 5) Outer foreach (parallel over compounds)
##------------------------------------------------

preproc_names <- c("movav","firstDer","secondDer","sg","sgFirstDer","sgSecondDer","msc","snv","snv_dt","no_prepro")
osc_modes <- c("none","osc")

Compounds_group = c("Aliphatics", "Alkanes", "Alkenes", "Benzenes", "Carbohydrates", "Lignins", 
                    "Methyl.ketones", "N.compounds", "Phenols", "Unidentified", "Sum.aliphatic1", "Sum.aliphatic2")

results_list <- foreach(w = seq_along(Compounds_group),
                         .packages = c("caret","pls","prospectr","MLmetrics","tidyverse"),
                         .export = c("apply_preproc","fit_osc","predict_osc","group_cv"),
                         .combine = rbind,
                         .multicombine = TRUE) %dopar% {
                           
                          
                           index <- Compounds_group[w]
                           y_train_orig <- Org_Train_data[[index]]
                           y_test_orig <- Org_Test_data[[index]]
                           
                           train_spc_orig <-Org_Train_data$spc
                           test_spc_orig <-Org_Test_data$spc
                           
                           local_row <- list()
                           
                           for(p_idx in seq_along(preproc_names)){
                             pre <- apply_preproc(p_idx, train_spc_orig, test_spc_orig)
                             Xtr_p <- pre$Xtr
                             Xte_p <- pre$Xte
                             preprocName <- pre$name
                             
                             #center spectra by training mean & apply to test
                             Train_mean <- colMeans(Xtr_p)
                             Xtr_centered <- sweep(Xtr_p, 2, Train_mean, "-")
                             Xte_centered <- sweep(Xte_p, 2, Train_mean, "-")
                             
                             #center target using training mean 
                             ytr_centered <- as.numeric(y_train_orig - mean(y_train_orig))
                             yte_centered <- as.numeric(y_test_orig - mean(y_train_orig))
                             
                             for(osc_mode in osc_modes) {
                               if(osc_mode == "none") {
                                 Xtr_final <- Xtr_centered
                                 Xte_final <- Xte_centered
                                 osc_mode <- "none"
                               } else {
                                 osc_model <- fit_osc(Xtr_centered, matrix(ytr_centered, ncol=1), ncomp = 1, center = FALSE)
                                 Xtr_final <- osc_model$Xcorr
                                 Xte_final <- predict_osc(Xte_centered, osc_model)
                                 osc_mode <- "OSC"
                               }
                               
                               x_train <- as.matrix(Xtr_final)
                               y_train <- ytr_centered
                               
                               #Train PLS via caret
                               
                              set.seed(123)
                              plsFit <- train(x = x_train, y = y_train,
                                              method = "pls",
                                              tuneLength = 20,
                                              trControl = train_ctrl_template,
                                              preProc = NULL,
                                              metric = "RMSE")
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
                              
                              #Fit final PLS with plsr to get predict
                              df_train_pls <- data.frame(y = y_train, as.data.frame(x_train))
                              
                              modelFit <- plsr(y ~ ., data = df_train_pls, ncomp = Best_nb_latents, center = FALSE, scale = FALSE)
                              
                              df_train_new <- data.frame(as.data.frame(x_train))
                              df_test_new  <- data.frame(as.data.frame(Xte_final))
                              
                              
                              #predict on train and test
                              
                              Train_result_vec <- as.numeric(predict(modelFit, newdata = df_train_new, ncomp = Best_nb_latents))
                              Test_result_vec  <- as.numeric(predict(modelFit, newdata = df_test_new,  ncomp = Best_nb_latents))
                              
                              
                              RMSE_train <- RMSE(Train_result_vec, y_train)
                              R2_train <- R2_Score(Train_result_vec, y_train)
                              MAE_train <- MAE(Train_result_vec, y_train)
                              
                              RMSE_test <- RMSE(Test_result_vec, yte_centered)
                              R2_test <- R2_Score(Test_result_vec, yte_centered)
                              MAE_test <- MAE(Test_result_vec, yte_centered)
                              
                              # CV summary
                              
                              RMSE_CV_mean <- as.numeric(BestTune_selection[Best_nb_latents,2]) #RMSE_CV_mean
                              RMSE_CV_se <- as.numeric(BestTune_selection[Best_nb_latents,3]) #RMSE_CV_se
                              RMSE_mean_plus_se <- as.numeric(BestTune_selection[Best_nb_latents,4]) #RMSE_mean_plus_se
                              
                              R2_CV <- plsFit$results$Rsquared[Best_nb_latents]  #Rsqured_CV
                              MAE_CV <- plsFit$results$MAE[Best_nb_latents]       #MAE_CV
                              #
                              
                              local_row[[length(local_row)+1]] <- data.frame(
                                ID = paste("MIR",index, preprocName,osc_mode,Best_nb_latents, sep = "_"),
                                compound = index,
                                preprocessing = preprocName,
                                osc = osc_mode,
                                best_ncomp = Best_nb_latents,
                                RMSE_CV = RMSE_CV_mean,
                                RMSE_CV_se = RMSE_CV_se,
                                RMSE_CV_plus_se = RMSE_mean_plus_se,
                                R2_CV = R2_CV,
                                MAE_CV = MAE_CV,
                                RMSE_train = RMSE_train,
                                R2_train = R2_train,
                                MAE_train = MAE_train,
                                RMSE_test = RMSE_test,
                                R2_test = R2_test,
                                MAE_test = MAE_test,
                                stringsAsFactors = FALSE
                                )
                              }
                             }
                           do.call(rbind, local_row)
                           }
                           
          
stopCluster(cl)


#save(results_list, file = "Data/PLS_model_revision/results_list.RData")
load("Data/PLS_model_revision/results_list.RData")

results_list %>% 
  group_by(compound) %>%
  dplyr::filter(RMSE_CV < RMSE_CV_plus_se[which.min(RMSE_CV)])  %>%
  group_by(compound) %>%
  dplyr::filter(best_ncomp == min(best_ncomp)) %>% 
  group_by(compound) %>%
  dplyr::filter(RMSE_CV  == min(RMSE_CV)) %>% 
  mutate(compound = factor(compound, levels = Compounds_group)) %>%
  arrange(compound) -> best_combination
  



###########################################################################################################


cl <- makeCluster(12)
registerDoParallel(cl)


results2_list <- foreach(i = seq_len(nrow(best_combination)),
                         .combine = rbind, 
                         .packages = c("prospectr","caret","mdatools","pls","tidyverse","MLmetrics")) %dopar% {
                           
                           
                          
                           
                           # ---- read row i from best_combination ----
                           index <- as.character(best_combination$compound[i])
                           preprocName <- best_combination$preprocessing[i]
                           best_comp <- best_combination$best_ncomp[i]
                           osc_mode_in <- best_combination$osc[i]
                           
                             
                           # ---- original data (do not modify originals) ----
                           y_train_orig <- Org_Train_data[[index]]
                           y_test_orig <- Org_Test_data[[index]]
                           train_spc_orig <-Org_Train_data$spc
                           test_spc_orig <-Org_Test_data$spc
                           
                           # ---- map preprocessing name to index for apply_preproc ----
                           preproc_map <- c("movav", "firstDer", "secondDer", "sg", "sgFirstDer", 
                                            "sgSecondDer", "msc", "snv", "snv_dt", "no_prepro")
                           
                           p_idx <- match(preprocName, preproc_map)
                           
                           
                           # ---- apply preprocessing ----
                           pre <- apply_preproc(p_idx, train_spc_orig, test_spc_orig)
                           Xtr_p <- pre$Xtr
                           Xte_p <- pre$Xte
                           
            
                           #center spectra by training mean & apply to test
                           Train_mean <- colMeans(Xtr_p)
                           Xtr_centered <- sweep(Xtr_p, 2, Train_mean, "-")
                           Xte_centered <- sweep(Xte_p, 2, Train_mean, "-")
                            
                           #center target using training mean 
                           ytr_centered <- as.numeric(y_train_orig - mean(y_train_orig))
                           yte_centered <- as.numeric(y_test_orig - mean(y_train_orig))
                            
                           # ---- apply OSC if required ----
                           if(osc_mode_in == "none") {
                             Xtr_final <- Xtr_centered
                             Xte_final <- Xte_centered
                             osc_mode <- "none"
                           } else {
                             osc_model <- fit_osc(Xtr_centered, matrix(ytr_centered, ncol=1), ncomp = 1, center = FALSE)
                             Xtr_final <- osc_model$Xcorr
                             Xte_final <- predict_osc(Xte_centered, osc_model)
                             osc_mode <- "OSC"
                           }
                           
                           # ---base x_train/ y_train ----
                           x_train_full <- as.matrix(Xtr_final)
                           y_train <- ytr_centered
                          
                           
                           #Train PLS via caret
                           set.seed(123)
                           
                           ##------------------------------------------------- 
                           ## select wave numbers from variable selection and VIP values
                           ##-------------------------------------------------
                           BCM <- pls(x = x_train_full, y = y_train,
                                      cv = group_by_ipls(Org_Train_data$Site, 13),
                                      ncomp = best_comp,
                                      scale = FALSE)
                           
                           Vip_scores <- vipscores(BCM, ncomp = best_comp)
                           
                           if(is.matrix(Vip_scores)) vip_vec <- Vip_scores[,1] else vip_vec <- Vip_scores
                           vip_thresh <- median(vip_vec, na.rm = TRUE)
                           VIP_idx <- which(vip_vec > vip_thresh)
                          
                          
                           #################################################################################################
                           #extract jack-knife values
                           
                           jack_pvals <- BCM$coeffs$p.values[, best_comp, 1]
                           var_jack_idx <- which(jack_pvals < 0.05)
                         
                           # ipls# extract interval-pls values
                           ipls_fwd <- ipls(x = x_train_full, y = y_train,
                                            method = "forward",
                                            cv = group_by_ipls(Org_Train_data$Site, 13),
                                            ncomp.selcrit = "min",
                                            glob.ncomp = best_comp,
                                            int.ncomp = best_comp,
                                            int.num = 30,
                                            silent = TRUE)
                          
                           
                           values_fwd <- sort(ipls_fwd$var.selected)
                           names_vec_fwd <- as.numeric(colnames(x_train_full[,values_fwd]))
                           
                           ipls_fwd_idx <- setNames(values_fwd, names_vec_fwd)
                           
                           ipls_bwd  <- ipls(x = x_train_full, y = y_train,
                                             method = "backward",
                                             cv = group_by_ipls(Org_Train_data$Site, 13), 
                                             ncomp.selcrit = "min",
                                             glob.ncomp = best_comp,
                                             int.ncomp = best_comp,
                                             int.num = 30,
                                             silent = TRUE)
                           
                           values_bwd <- sort(ipls_bwd$var.selected)
                           names_vec_bwd <- as.numeric(colnames(x_train_full[,values_bwd]))
                           
                           ipls_bwd_idx <- setNames(values_bwd, names_vec_bwd)
                           
                           
                           raw_idx <- seq_len(ncol(x_train_full))
                          #################################################################################################
                          #Create loop
                          variable_selection_methods <- list(
                            Raw = raw_idx,
                            VIP = VIP_idx,
                            var_jack = var_jack_idx,
                            ipls_forward = ipls_fwd_idx,
                            ipls_backward = ipls_bwd_idx
                            )
                           
                          
                           
                          for(method_name in names(variable_selection_methods)) {
                            
                            variable_sel_method <- method_name
                            sel_idx <- variable_selection_methods[[method_name]]
                            if(length(sel_idx) == 0) next
                            
                            x_train_sel <- x_train_full[, sel_idx, drop = FALSE]
                            x_test_sel <- Xte_final[, sel_idx, drop = FALSE]
                            
                            
                            set.seed(123)
                            
                            ###################################################################################
                            plsFit <- train(x = x_train_sel, y = y_train,
                                            method = "pls",
                                            tuneLength = 20,
                                            trControl = train_ctrl_template,
                                            preProc = NULL,
                                            metric = "RMSE")
                            
                            #-----------------------------------------------------------------------------
                            # One-sided SE paired t-test to select the best number of latent variables
                            #-----------------------------------------------------------------------------
                            BestTune_selection <- plsFit$pred %>%
                              mutate(New_ID = paste0(ncomp, Resample)) %>%
                              group_by(New_ID) %>%
                              summarise(RMSE = RMSE(obs, pred), ncomp = first(ncomp), Fold =  first(Resample)) %>%
                              arrange(ncomp,Fold) %>%
                              group_by(ncomp) %>%
                              summarise(RMSE_mean = mean(RMSE), RMSE_se = sd(RMSE)/sqrt(n()), ncomp =first(ncomp)) %>%
                              mutate(RMSE_mean_plus_se = RMSE_mean+RMSE_se)
                            
                            Best_nb_latents <- min(which(BestTune_selection$RMSE_mean < BestTune_selection$RMSE_mean_plus_se[plsFit$bestTune[,1]])) 
                            
                            #Fit final PLS with plsr to get predict
                            df_train_pls <- data.frame(y = y_train, as.data.frame(x_train_sel))
                            
                            modelFit <- plsr(y ~ ., data = df_train_pls, ncomp = Best_nb_latents, 
                                             center = FALSE, scale = FALSE)
                            
                            df_train_new <- data.frame(as.data.frame(x_train_sel))
                            df_test_new  <- data.frame(as.data.frame(x_test_sel))
                            
                           
                            #predict on train and test
                            
                            Train_result_vec <- as.numeric(predict(modelFit, newdata = df_train_new, ncomp = Best_nb_latents))
                            Test_result_vec  <- as.numeric(predict(modelFit, newdata = df_test_new,  ncomp = Best_nb_latents))
                            #----------------------------------------------------------------------------------
                            #metrics##
              
                            RMSE_train <- RMSE(Train_result_vec, y_train)
                            R2_train <- R2_Score(Train_result_vec, y_train)
                            MAE_train <- MAE(Train_result_vec, y_train)
                            
                            RMSE_test <- RMSE(Test_result_vec, yte_centered)
                            R2_test <- R2_Score(Test_result_vec, yte_centered)
                            MAE_test <- MAE(Test_result_vec, yte_centered)
                            
                            # CV summary
                            
                            RMSE_CV_mean <- as.numeric(BestTune_selection[Best_nb_latents,2]) #RMSE_CV_mean
                            RMSE_CV_se <- as.numeric(BestTune_selection[Best_nb_latents,3]) #RMSE_CV_se
                            RMSE_mean_plus_se <- as.numeric(BestTune_selection[Best_nb_latents,4]) #RMSE_mean_plus_se
                            
                            R2_CV <- plsFit$results$Rsquared[Best_nb_latents]  #Rsqured_CV
                            MAE_CV <- plsFit$results$MAE[Best_nb_latents]       #MAE_CV
                            
                            # collect row
                            local_row[[length(local_row)+1]] <- data.frame(
                              ID = paste("MIR",index, preprocName,osc_mode,variable_sel_method, Best_nb_latents, sep = "_"),
                              compound = index,
                              preprocessing = preprocName,
                              osc = osc_mode,
                              variableSel = variable_sel_method,
                              best_ncomp = Best_nb_latents,
                              RMSE_CV = RMSE_CV_mean,
                              RMSE_CV_se = RMSE_CV_se,
                              RMSE_CV_plus_se = RMSE_mean_plus_se,
                              R2_CV = R2_CV,
                              MAE_CV = MAE_CV,
                              RMSE_train = RMSE_train,
                              R2_train = R2_train,
                              MAE_train = MAE_train,
                              RMSE_test = RMSE_test,
                              R2_test = R2_test,
                              MAE_test = MAE_test,
                              stringsAsFactors = FALSE
                              )
                          
                          
                            ID <- paste("MIR",index, preprocName,osc_mode,variable_sel_method, Best_nb_latents, sep = "_")
                              ###########################################################################
                              #Save the results for plots
                              ############################################################################
                              CV_results <- data.frame(plsFit$pred[plsFit$pred$ncomp == Best_nb_latents,]) #CV results
                              
                              Cal_results <- data.frame(SampleID = row.names((x_train)),                #Calibration results
                                                        pred=Train_result_vec + mean(y_train_orig), 
                                                        obs = y_train + mean(y_train_orig))
                              
                              Val_results <- data.frame(SampleID = row.names(Xte_final),                 #Validation results
                                                        pred= Test_result_vec + mean(y_train_orig), 
                                                        obs = yte_centered + mean(y_train_orig))
                              
                              coefficent.dat <- data.frame(modelFit$coefficients)               #Coefficient results
                              
                              Wavenumbers_results <- coefficent.dat %>% 
                                mutate(Wavenumber = rownames(.),
                                       Coefficient = coefficent.dat[,Best_nb_latents]) %>%  
                                dplyr::select(Wavenumber,Coefficient) %>%
                                mutate(Wavenumber = as.numeric(gsub("X","",Wavenumber)))

                              getwd()
                              
                              # Create a directory to store the results
                              dir.create(paste0("PLS_model_revision/loop_test/", index), showWarnings = FALSE)
                              
                              
                              # List of results and corresponding file suffixes
                              results_list <- list(CV_results = "_CV.txt", Cal_results = "_Calibration.txt", 
                                                   Val_results = "_Validation.txt", Wavenumbers_results = "_wavenumber.txt")
                              
                              
                              # Loop through the results list and write tables
                              for (result_name in names(results_list)) {
                                tableName <- paste0(Path, "PLS_model_revision/loop_test/", index, "/", ID, results_list[[result_name]])
                                write.table(get(result_name), file = tableName, row.names = FALSE)
                              }
                            }
                            do.call(rbind, local_row)
                          }


stopCluster(cl)
#save(results2_list, file = "Data/PLS_model_revision/results2_list.RData")
load("Data/PLS_model_revision/results2_list.RData")
results2_list %>% 
  group_by(compound) %>%
  dplyr::filter(RMSE_CV < RMSE_CV_plus_se[which.min(RMSE_CV)])  %>%
  group_by(compound) %>%
  dplyr::filter(best_ncomp == min(best_ncomp)) %>% 
  group_by(compound) %>%
  dplyr::filter(RMSE_CV  == min(RMSE_CV)) %>% 
  slice(1) %>% 
  mutate(compound = factor(compound, levels = Compounds_group)) %>%
  arrange(compound) -> best_combination_R2
 
View(best_combination_R2)
getwd()
#write.csv(best_combination_R2, file = "PLS_model_revision/best_combination_R22.csv", row.names = FALSE)
#####################################################################################
#VALIDATION DATA PLOT#
#####################################################################################
all_new_data <- list()
for(mirs in 1:length(Compounds_group)){
  

  #data import
  Index <- Compounds_group[mirs]
  
  #import plot
  validation.dat <-read.table(paste0("PLS_model_revision/loop_test/",Index,"/",best_combination_R2$ID[mirs],"_Validation.txt"),header = T)
  validation.dat$Index <- Compounds_group[mirs]
  #split_strings <- strsplit(best_combination$ID[7], "_")
  colnames(validation.dat)[2]       <- "pred"
  validation.dat$Best_preprocessing <- best_combination_R2$preprocessing[mirs]
  validation.dat$Best_variable      <- best_combination_R2$variableSel[mirs]
  validation.dat$Best_comp          <- best_combination_R2$best_ncomp[mirs]
  validation.dat$Rsquared           <- best_combination_R2$R2_test[mirs]
  validation.dat$RMSE               <- best_combination_R2$RMSE_test[mirs]
  validation.dat$MAE                <- best_combination_R2$MAE_test[mirs]
  
  all_new_data[[mirs]] <- validation.dat
}


final_new_data <- do.call(rbind, all_new_data) %>% 
  mutate(Rsquared = as.numeric(Rsquared),
         RMSE = as.numeric(RMSE),
         MAE = as.numeric(MAE))%>% 
  mutate_if(is.numeric,round,digits=2)


final_new_data$Index <- factor(final_new_data$Index,
                               levels = c("Alkanes","Alkenes","Benzenes","Carbohydrates",
                                          "Lignins","Methyl.ketones","N.compounds","Phenols",
                                          "Aliphatics","Unidentified","Sum.aliphatic1",
                                          "Sum.aliphatic2"),
                               labels = c("Alkanes","Alkenes","Benzenes","Carbohydrates",
                                          "Lignins","Methyl ketones","N compounds","Phenols",
                                          "Unspecified aliphatics","Unidentified","Sum.aliphatic1",
                                          "Sum.aliphatic2"))

final_new_data2 <- subset(final_new_data, !(Index %in% 
                                              c("Unidentified","Sum.aliphatic1","Sum.aliphatic2")))

ggplot(final_new_data2,aes(x = obs, y=pred))+
  geom_point()+
  geom_abline(color="grey",size = 1)+
  theme_bw()+
  facet_wrap(Index ~ ., scales = "free",  ncol= 3)+
  #geom_text(aes(x = Inf, y = -Inf, label = paste("Best_preprocess =", Best_preprocessing)), vjust = -8, hjust = 1, size = 4)+
  #geom_text(aes(x = Inf, y = -Inf, label = paste("Best_wavenumber =", Best_variable)), vjust = -6.5, hjust = 1, size = 4)+
  #geom_text(aes(x = Inf, y = -Inf, label = paste("Best_ncomp=", Best_comp)), vjust = -5, hjust = 1, size = 4)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("Rsq =", Rsquared)), vjust = -3.5, hjust = 1, size = 4)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("RMSE =", RMSE)), vjust = -2, hjust = 1, size = 4)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("MAE =", MAE)), vjust = -0.5, hjust = 1, size = 4)+
  xlab("Reference (%)")+
  ylab("Predicted (%)")



#####################################################################################
#CALIBRATION DATA PLOT#
#####################################################################################
all_new_data <- list()
for(mirs in 1:length(Compounds_group)){
  
  
  #data import
  Index <- Compounds_group[mirs]
  
  #import plot
  calibration.dat <-read.table(paste0("PLS_model_revision/loop_test/",Index,"/",best_combination_R2$ID[mirs],"_Calibration.txt"),header = T)
  calibration.dat$Index <- Compounds_group[mirs]
  #split_strings <- strsplit(best_combination_2$ID[7], "_")
  colnames(calibration.dat)[2]       <- "pred"
  calibration.dat$Best_preprocessing <- best_combination_R2$preprocessing[mirs]
  calibration.dat$Best_variable      <- best_combination_R2$variableSel[mirs]
  calibration.dat$Best_comp          <- best_combination_R2$best_ncomp[mirs]
  calibration.dat$Rsquared           <- best_combination_R2$R2_train[mirs]
  calibration.dat$RMSE               <- best_combination_R2$RMSE_train[mirs]
  calibration.dat$MAE                <- best_combination_R2$MAE_train[mirs]
  
  all_new_data[[mirs]] <- calibration.dat
}



final_new_data <- do.call(rbind, all_new_data) %>% 
  mutate(Rsquared = as.numeric(Rsquared),
         RMSE = as.numeric(RMSE),
         MAE = as.numeric(MAE))%>% 
  mutate_if(is.numeric,round,digits=2)


final_new_data$Index <- factor(final_new_data$Index,
                               levels = c("Alkanes","Alkenes","Benzenes","Carbohydrates",
                                          "Lignins","Methyl.ketones","N.compounds","Phenols",
                                          "Aliphatics","Unidentified","Sum.aliphatic1",
                                          "Sum.aliphatic2"),
                               labels = c("Alkanes","Alkenes","Benzenes","Carbohydrates",
                                          "Lignins","Methyl ketones","N compounds","Phenols",
                                          "Unspecified aliphatics","Unidentified","Sum.aliphatic1",
                                          "Sum.aliphatic2"))

final_new_data2 <- subset(final_new_data, !(Index %in% 
                                              c("Unidentified","Sum.aliphatic1","Sum.aliphatic2")))


ggplot(final_new_data2,aes(x = obs, y=pred))+
  geom_point()+
  geom_abline(color="grey",size = 1)+
  theme_bw()+
  facet_wrap(Index ~ ., scales = "free",  ncol= 3)+
  # geom_text(aes(x = Inf, y = -Inf, label = paste("Best_preprocess=", Best_preprocessing)), vjust = -8, hjust = 1, size = 4)+
  # geom_text(aes(x = Inf, y = -Inf, label = paste("Best_wavenumber=", Best_variable)), vjust = -6.5, hjust = 1, size = 4)+
  # geom_text(aes(x = Inf, y = -Inf, label = paste("Best_ncomp=", Best_comp)), vjust = -5, hjust = 1, size = 4)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("Rsq =", Rsquared)), vjust = -3.5, hjust = 1, size = 4)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("RMSE =", RMSE)), vjust = -2, hjust = 1, size = 4)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("MAE =", MAE)), vjust = -0.5, hjust = 1, size = 4)+
  xlab("Reference (%)")+
  ylab("Predicted (%)")



#####################################################################################
#CROSS-VALIDATINO DATA PLOT#
#####################################################################################
all_new_data <- list()
for(mirs in 1:length(Compounds_group)){
  
  
  #data import
  Index <- Compounds_group[mirs]
  
  #import plot
  CV.dat <-read.table(paste0("PLS_model_revision/loop_test/",Index,"/",best_combination_R2$ID[mirs],"_CV.txt"),header = T)
  CV.dat$Index <- Compounds_group[mirs]
  #split_strings <- strsplit(best_combination_R2$ID[7], "_")
  
  CV.dat$Best_preprocessing <- best_combination_R2$preprocessing[mirs]
  CV.dat$Best_variable      <- best_combination_R2$variableSel[mirs]
  CV.dat$Best_comp          <- best_combination_R2$best_ncomp[mirs]
  CV.dat$Rsquared           <- best_combination_R2$R2_CV[mirs]
  CV.dat$RMSE               <- best_combination_R2$RMSE_CV[mirs]
  CV.dat$MAE                <- best_combination_R2$MAE_CV[mirs]
  
  all_new_data[[mirs]] <- CV.dat
}
best_combination_R2

final_new_data <- do.call(rbind, all_new_data) %>% 
  mutate(Rsquared = as.numeric(Rsquared),
         RMSE = as.numeric(RMSE),
         MAE = as.numeric(MAE))%>% 
  mutate_if(is.numeric,round,digits=2)


final_new_data$Index <- factor(final_new_data$Index,
                               levels = c("Alkanes","Alkenes","Benzenes","Carbohydrates",
                                          "Lignins","Methyl.ketones","N.compounds","Phenols",
                                          "Aliphatics","Unidentified","Sum.aliphatic1",
                                          "Sum.aliphatic2"),
                               labels = c("Alkanes","Alkenes","Benzenes","Carbohydrates",
                                          "Lignins","Methyl ketones","N compounds","Phenols",
                                          "Unspecified aliphatics","Unidentified","Sum.aliphatic1",
                                          "Sum.aliphatic2"))


final_new_data2 <- subset(final_new_data, !(Index %in% 
                                              c("Unidentified","Sum.aliphatic1","Sum.aliphatic2")))

ggplot(final_new_data2,aes(x = obs, y=pred))+
  geom_point()+
  geom_abline(color="grey",size = 1)+
  theme_bw()+
  facet_wrap(Index ~ ., scales = "free",  ncol= 3)+
  #geom_text(aes(x = Inf, y = -Inf, label = paste("Best_preprocess=", Best_preprocessing)), vjust = -8, hjust = 1, size = 4)+
  #geom_text(aes(x = Inf, y = -Inf, label = paste("Best_wavenumber=", Best_variable)), vjust = -6.5, hjust = 1, size = 4)+
  #geom_text(aes(x = Inf, y = -Inf, label = paste("Best_ncomp=", Best_comp)), vjust = -5, hjust = 1, size = 4)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("Rsq =", Rsquared)), vjust = -3.5, hjust = 1, size = 4)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("RMSE =", RMSE)), vjust = -2, hjust = 1, size = 4)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("MAE =", MAE)), vjust = -0.5, hjust = 1, size = 4)+
  xlab("Reference (%)")+
  ylab("Predicted (%)")

#####################################################################################
#Coefficient data
#####################################################################################
Raw <- as.numeric(colnames(x_train_full))

all_new_data <- list()
for(mirs in 1:length(Compounds_group)){
  
  #data import
  Index <- Compounds_group[mirs]
  
  #import plot
  wavenumber.dat <-read.table(paste0("PLS_model_revision/loop_test/",Index,"/",best_combination_R2$ID[mirs],"_wavenumber.txt"),header = T)
  wavenumber.dat$Index <- Compounds_group[mirs]
  all_new_data[[mirs]] <- wavenumber.dat
}

all_new_data

final_new_data <- do.call(rbind, all_new_data)

Reference_wave_data <- data.frame(Wavenumber = rep(Raw, length(Compounds_group)),
                                  Index = rep(Compounds_group, each = length(Raw)))

variable_sel_spectra <-final_new_data
View(variable_sel_spectra)

Reference_wave_data <- left_join(Reference_wave_data, final_new_data, by = c("Wavenumber","Index")) %>% 
  mutate(Coefficient = ifelse(is.na(Coefficient),0,Coefficient))

Reference_wave_data <- subset(Reference_wave_data, Index %in%  
                                c("Alkenes","Benzenes","Carbohydrates","Lignins","N.compounds","Phenols"))

Reference_wave_data$Index <- factor(Reference_wave_data$Index,
                               levels = c("Alkenes", "Benzenes","Carbohydrates","Lignins","N.compounds","Phenols"),
                               labels = c("Alkenes", "Benzenes","Carbohydrates","Lignins","N compounds","Phenols"))

unique(Reference_wave_data$Index)

# Reference_wave_data <- Reference_wave_data %>%
#   group_by(Index) %>%
#   mutate(Scaled_Coefficient = scale(Coefficient)) %>%
#   ungroup() 

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

PLS_coefficients_Plot <- 
  ggplot(Reference_wave_data,aes(x=Wavenumber,y=Coefficient,group =Index, color = Index))+
  geom_line(size=0.8)+
  theme_bw()+
  theme(legend.position = "top",
        legend.title =element_blank(),
        panel.grid.major.x = element_line(color = "black", size = 0.1, linetype = 2),
        panel.grid.major.y = element_line(color = "grey",size = 0.1,linetype = 2))+
  xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('PLS regression coefficients')+
  scale_x_reverse(n.breaks = 10)+
  grids(linetype = "dashed")+
  scale_color_manual(values=cbPalette,guide  = guide_legend(nrow = 1, byrow = TRUE)) +
  facet_wrap(~Index, scales="free", ncol=2)

###############################################################################################################
# Permutation test##
#############################
cores <- parallel::detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

results3_list <- foreach(i = seq_len(nrow(best_combination_R2)),
                         .combine = rbind, 
                         .packages = c("prospectr","caret","mdatools","pls","tidyverse","MLmetrics")) %dopar% {
                           
                                      
                                       
                                       
                      
                        # ---- read row i from best_combination ----
                        index <- as.character(best_combination_R2$compound[i])
                        preprocName <- best_combination_R2$preprocessing[i]
                        best_comp <- best_combination_R2$best_ncomp[i]
                        osc_mode_in <- best_combination_R2$osc[i]
                        
                        
                        # ---- original data (do not modify originals) ----
                        y_train_orig <- Org_Train_data[[index]]
                        y_test_orig <- Org_Test_data[[index]]
                        train_spc_orig <-Org_Train_data$spc
                        test_spc_orig <-Org_Test_data$spc
                        
                        # ---- map preprocessing name to index for apply_preproc ----
                        preproc_map <- c("movav", "firstDer", "secondDer", "sg", "sgFirstDer", 
                                         "sgSecondDer", "msc", "snv", "snv_dt", "no_prepro")
                        
                        p_idx <- match(preprocName, preproc_map)
                        
                        # ---- apply preprocessing ----
                        pre <- apply_preproc(p_idx, train_spc_orig, test_spc_orig)
                        Xtr_p <- pre$Xtr
                        Xte_p <- pre$Xte
                        
                        #center spectra by training mean & apply to test
                        Train_mean <- colMeans(Xtr_p)
                        Xtr_centered <- sweep(Xtr_p, 2, Train_mean, "-")
                        Xte_centered <- sweep(Xte_p, 2, Train_mean, "-")
                        
                        #center target using training mean 
                        ytr_centered <- as.numeric(y_train_orig - mean(y_train_orig))
                        yte_centered <- as.numeric(y_test_orig - mean(y_train_orig))
                        
                        # ---- apply OSC if required ----
                        if(osc_mode_in == "none") {
                          Xtr_final <- Xtr_centered
                          Xte_final <- Xte_centered
                          osc_mode <- "none"
                        } else {
                          osc_model <- fit_osc(Xtr_centered, matrix(ytr_centered, ncol=1), ncomp = 1, center = FALSE)
                          Xtr_final <- osc_model$Xcorr
                          Xte_final <- predict_osc(Xte_centered, osc_model)
                          osc_mode <- "OSC"
                        }
                        
                        # ---base x_train/ y_train ----
                        x_train_full <- as.matrix(Xtr_final)
                        y_train <- ytr_centered
                        
                        #-------------------------------------------------
                        sel_wavenumber <- variable_sel_spectra$Wavenumber[variable_sel_spectra$Index == index]
                        sel_wavenumber <- as.numeric(sel_wavenumber)
                        
                        x_train_sel <- x_train_full[, as.character(sel_wavenumber), drop = FALSE]
                        x_test_sel <- Xte_final[, as.character(sel_wavenumber), drop = FALSE]
                        
                        get_cv_rmse <- function(x,y) {
                          
                          plsFit <- train(x = x, y = y,
                                        method = "pls",
                                        tuneGrid = expand.grid(ncomp = best_comp) ,
                                        trControl = train_ctrl_template,
                                        preProc = NULL,
                                        metric = "RMSE")
                          plsFit$results[["RMSE"]]
                          
                        }
                        
                        observed <- get_cv_rmse(x_train_sel,y_train)
                        
                        
                        ###permutation scores ##########
                        n_perm <- 999
                        
                        perm_scores <- numeric(n_perm) 
                        
                        for(p in seq_len(n_perm)) {
                          perm_scores[p] <- get_cv_rmse(x_train_sel, sample(y_train))
                        }
                        
                        pvalue <- (1 + sum(perm_scores <= observed)) / (1 + n_perm)
                        
                        data.frame(
                          compound = index,
                          p_value = pvalue,
                          stringsAsFactors = FALSE)
                        }



stopCluster(cl)

#save(results3_list, file = "Data/PLS_model_revision/results3_list.RData")
load("Data/PLS_model_revision/results3_list.RData")
results3_list
  







