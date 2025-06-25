#################################################################################
## Compute effective (R) and angular reproduction numbers (Omega) on simulalted
# data with no variations in generation time
#################################################################################

# Set working directory to source
setwd("C:/Users/esthe/OneDrive - Nanyang Technological University/Rt/")
setwd("/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/Rt/")
# Folder path for results
folres = paste0("./plots/")

# Main functions to run EpiFilter
files.sources = list.files(path = "./omegaFunctions")
for (i in 1:length(files.sources)) {
  source(paste0(c("./omegaFunctions/", files.sources[i]), collapse = ''))
}

# Load workspace
load("./out/rtdata.Rdata")
load("./data/temp.Rdata")

#########################Not Temp-Dependent########################################
################## Input daily dengue data and length of time series ##############

nday <- length(cases)
tday <- 1:nday
Lday <- totalinf
Iday <- cases

# For window delta compute rmse infections
Irms = rep(0, nday); delta = 28 #set delta to 2x expected mean generation time
for(i in 2:nday){
  # Edge effect cases
  if(i-1 < delta){
    Irms[i] = sqrt(sum(Iday[seq(i-1, 1, -1)]^2))/sqrt(i-1)
  }else{
    # Omega only truly valid from this point
    Irms[i] = sqrt(sum(Iday[seq(i-1, i-delta, -1)]^2))/sqrt(delta)
  }
}

#################################################################################
################ Estimate R and Omega using EpiFilter ###########################

# Setup grid [Rmin Rmax] noise eta and CIs conf
Rmin = 0.01; Rmax = 10; eta = 0.1
# Uniform prior over grid of size m
m = 1000; pR0 = (1/m)*rep(1, m)

# Delimited grid defining space of R
Rgrid = seq(Rmin, Rmax, length.out = m)
# Maximum of prediction grid and confidence (1-conf)
Imax = 2*max(Iday); conf = 0.025

# Filtered (causal) estimates as list [Rmed, Rci, Rmean, pR, pRup, pstate]
Rfilt = epiFilter(Rgrid, m, eta, pR0, nday, Lday[tday], Iday[tday], conf)
# Causal predictions from filtered estimates [pred predci]
#Ifilt = recursPredict(Rgrid, Rfilt$pR, Lday[tday], Rfilt$Rmean, conf, Iday, Imax) 
#// just to compare filter vs smooth predictions, Ifilt doesnt go to smoothing step

# Smoothed estimates as list of [Rmed, Rhatci, Rmean, qR]
Rsmooth = epiSmoother(Rgrid, m, Rfilt$pR, Rfilt$pRup, nday, Rfilt$pstate, conf)
# Smoothed predictions from filtered estimates [pred predci]
Ismooth = recursPredict(Rgrid, Rsmooth$qR, Lday[tday], Rsmooth$Rmean, conf, Iday, Imax)

# Prob of R > 1 (resurgence)
pR1 = rep(0, nday); id1 = which(Rgrid >= 1); id1 = id1[1]
for (i in 1:nday){
  pR1[i] = sum(Rsmooth$qR[i, id1:m]) #sum of probs when R>1 for each timepoint i
}

# Filtered (causal) estimates as list [Ommed, Omci, Ommean, pOm, pOmup, pOmstate]
Omfilt = epiFilter(Rgrid, m, eta, pR0, nday, Irms[tday], Iday[tday], conf)
# Causal predictions from filtered estimates [pred predci]
IfiltOm = recursPredict(Rgrid, Omfilt$pR, Irms[tday], Omfilt$Rmean, conf, Iday, Imax)

# Smoothed estimates as list of [Ommed, Omhatci, Ommean, qOm]
Omsmooth = epiSmoother(Rgrid, m, Omfilt$pR, Omfilt$pRup, nday, Omfilt$pstate, conf)
# Smoothed predictions from filtered estimates [pred predci]
IsmoothOm = recursPredict(Rgrid, Omsmooth$qR, Irms[tday], Omsmooth$Rmean, conf, Iday, Imax)
# Prob of R > 1 (resurgence)
pOm1 = rep(0, nday); id1 = which(Rgrid >= 1); id1 = id1[1]
for (i in 1:nday){
  pOm1[i] = sum(Omsmooth$qR[i, id1:m])
}

###########################Temp-Dependent######################################
############## Input daily dengue data and Length of time series ##############

nday <- length(cases)
tday <- 1:nday
Lday_td <- totalinf_td
Iday <- cases

# For window delta compute rmse infections
Irms = rep(0, nday); delta = 28 #set delta to 2x expected mean generation time
for(i in 2:nday){
  # Edge effect cases
  if(i-1 < delta){
    Irms[i] = sqrt(sum(Iday[seq(i-1, 1, -1)]^2))/sqrt(i-1)
  }else{
    # Omega only truly valid from this point
    Irms[i] = sqrt(sum(Iday[seq(i-1, i-delta, -1)]^2))/sqrt(delta)
  }
}

#################################################################################
################ Estimate R and Omega using EpiFilter ###########################

# Setup grid [Rmin Rmax] noise eta and CIs conf
Rmin = 0.01; Rmax = 10; eta = 0.1
# Uniform prior over grid of size m
m = 1000; pR0 = (1/m)*rep(1, m)

# Delimited grid defining space of R
Rgrid = seq(Rmin, Rmax, length.out = m)
# Maximum of prediction grid and confidence (1-conf)
Imax = 2*max(Iday); conf = 0.025

# Filtered (causal) estimates as list [Rmed, Rci, Rmean, pR, pRup, pstate]
Rfilt_td = epiFilter(Rgrid, m, eta, pR0, nday, Lday_td[tday], Iday[tday], conf)
# Causal predictions from filtered estimates [pred predci]
Ifilt_td = recursPredict(Rgrid, Rfilt_td$pR, Lday_td[tday], Rfilt_td$Rmean, conf, Iday, Imax) 
#// just to compare filter vs smooth predictions, Ifilt doesnt go to smoothing step

# Smoothed estimates as list of [Rmed, Rhatci, Rmean, qR]
Rsmooth_td = epiSmoother(Rgrid, m, Rfilt_td$pR, Rfilt_td$pRup, nday, Rfilt_td$pstate, conf)
# Smoothed predictions from filtered estimates [pred predci]
Ismooth_td = recursPredict(Rgrid, Rsmooth_td$qR, Lday_td[tday], Rsmooth_td$Rmean, conf, Iday, Imax)

# Prob of R > 1 (resurgence)
pR1_td = rep(0, nday); id1_td = which(Rgrid >= 1); id1_td = id1_td[1]
for (i in 1:nday){
  pR1_td[i] = sum(Rsmooth_td$qR[i, id1_td:m]) #sum of probs when R>1 for each timepoint i
}

# # Filtered (causal) estimates as list [Ommed, Omci, Ommean, pOm, pOmup, pOmstate]
# Omfilt_td = epiFilter(Rgrid, m, eta, pR0, nday, Irms[tday], Iday[tday], conf)
# # Causal predictions from filtered estimates [pred predci]
# IfiltOm_td = recursPredict(Rgrid, Omfilt_td$pR, Irms[tday], Omfilt_td$Rmean, conf, Iday, Imax)
# 
# # Smoothed estimates as list of [Ommed, Omhatci, Ommean, qOm]
# Omsmooth_td = epiSmoother(Rgrid, m, Omfilt_td$pR, Omfilt_td$pRup, nday, Omfilt_td$pstate, conf)
# # Smoothed predictions from filtered estimates [pred predci]
# IsmoothOm_td = recursPredict(Rgrid, Omsmooth_td$qR, Irms[tday], Omsmooth_td$Rmean, conf, Iday, Imax)
# # Prob of R > 1 (resurgence)
# pOm1_td = rep(0, nday); id1_td = which(Rgrid >= 1); id1_td = id1_td[1]
# for (i in 1:nday){
#   pOm1_td[i] = sum(Omsmooth_td$qR[i, id1_td:m])
# }

save.image("./real_data_results.Rdata")


