# Load packages -----------------------------------------------------------

(start.time <- Sys.time())

package.list <- c("jagsUI", "coda",
                  'dplyr', 'stringr',
                  'magrittr', 'tidyr',
                  'mcmcplots','ggplot2') 

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load model --------------------------------------------------------------

mod <- readRDS("/scratch/userID/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model2_checkN.RDS")



# Check Rhat --------------------------------------------------------------

Rhat <- mod$Rhat

saveRDS(Rhat, file ="/scratch/userID/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model2_Rhat.RDS")


# Get summary -------------------------------------------------------------

sum <- summary(mod$samples)

saveRDS(sum, "/scratch/userID/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model2_summary.RDS")
