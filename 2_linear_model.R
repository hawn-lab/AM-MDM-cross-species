library(tidyverse)
library(limma)
library(kimma)

#### Human Data ####
load("data_clean/human_dat_clean.RData")

#### Human linear model ####
model_hum <- kmFit(dat = voom_hum, 
      model="~cell*mtb + (1|ptID)",
      run_lme = TRUE, 
      run_contrast = TRUE, contrast_var = "cell:mtb",
      metrics = TRUE, use_weights = TRUE)

save(model_hum, file = "results/lme_human.RData")

#### Murine Data ####
load("data_clean/mouse_dat_clean.RData")

#### Murine linear model ####
model_mur <- kmFit(dat = voom_mur, 
                   model="~treatment*mtb",
                   run_lm = TRUE, 
                   run_contrast = TRUE, contrast_var = "treatment:mtb",
                   metrics = TRUE, use_weights = TRUE)

save(model_mur, file = "results/lm_mouse.RData")
