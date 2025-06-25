library("EpiEstim")
library("caTools")
library(xts)
library(coga)

#start with real temp data
#generate temp-dep Rt
#generate simulated cases based on known Rt
#estimate omega using simulated cases

##set up
setwd("C:/Users/esthe/OneDrive - Nanyang Technological University/Rt/")
setwd("/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/Rt/")
files.sources = list.files(path = "./omegaFunctions")
for (i in 1:length(files.sources)) {
  source(paste0(c("./omegaFunctions/", files.sources[i]), collapse = ''))
}

##load in real temperature data
load(file="./data/temp.Rdata")
n <- 200
temp1 <- temp[1:n]
temp2 <- rnorm(n, 28, 1.2)
temp3 <- rnorm(n, 28, 2)

period <- c(0.05, 0.1, 0.15)
ampli <- c(1, 2, 3)
smooth <- c("T", "F")
temp <- c("temp1", "temp2", "temp3") #increasing in variation

vars <- expand.grid(period=period, ampli=ampli, smooth=smooth, temp=temp)

simGenerator <- function(ampli, period, smooth, temp0, newlam){ #generate different temperature and Rt data set to apply to epifilter
  ##generate Rt using sine curve and random noise
  
  temp <- get(as.character(temp0))
  
  if(smooth == "F"){
    Rt <- ampli*sin(period*(1:n)) + rnorm(n, 0, 0.2) + 1.5
    Rt <- (Rt + abs(min(Rt))) / 2 + 0.2 #floor 0 and make median around 1
    #plot(Rt, type="l")
  } else {
    Rt <- ampli*sin(period*(1:n)) + 1.5
    Rt <- (Rt + abs(min(Rt))) / 2 + 0.2
    #plot(Rt, type="l")
  }
  
  ##generate cases using Rt (poisson dist)
  
  #calc current total infectiousness
  aH <- 16 ; sH <- 2.7 #shape and rate param of intrinsic incub (gamma dist)
  aM <- 4.3 ; b0 <- 7.9 ; b1 <- 0.21; #shape and rate params of extrinsic incub (gamma dist)
  lamTM <- newlam; lamTH <- newlam #say mosquito population reduces by 75% or triples
  
  real_inc <- round(runif(36, 0, 20)) #initiate case values so that we can initiate total infectiousness
  for(i in 1:(n-1)){
    s <- length(real_inc)
    t <- temp[i] #contemperaneous daily temp at t=36
    sM <- 1/(b0-b1*t) 
    w <- dcoga(1:35, shape=c(aM, aH, 1, 1), rate=c(sM, sH, lamTM, lamTH))
    # w_real <- dcoga(1:35, shape=c(aM, aH, 1, 1), rate=c(sM, sH, 1, 1))
    # gi <- c(gi, which.max(w))
    # gi_real <- c(gi_real, which.max(w_real))
    totalinf <- w %*% rev(real_inc[(s-34):(s)])
    real_inc <- c(real_inc, rpois(1, ((totalinf) * Rt[i])))
    #print(paste0("Iteration ", i, " of ", n-1, ":", totalinf))
    #print(paste0("Iteration ", i, " of ", n-1, ":", real_inc[i+36]))
    #plot(real_inc, type="l", ylab="Cases")
    
  }
  
return(list("Amp"=ampli, "Period"=period, "Smooth"=as.character(smooth), "Rt"=Rt, "Temp"=temp, "True Inc"=real_inc))  
}

params <- list()
lamset <- c(4, 3, 0.33, 0.25)
#for each set of params, output true incidence and true Rt
for(newlam in lamset){
  for(row in 1:nrow(vars)){
    params[[row]] <- simGenerator(vars$ampli[row], vars$period[row], vars$smooth[row], vars$temp[row], newlam)
  }
  
  newlam <- gsub("\\.", "", as.character(newlam))
  
  save(params, file=paste0("./out/VCpar", newlam, ".Rdata"))
}


#####train everything on each set of temp, Rt and true inc#####
#return list() of each metric vs real Rt agreement areas - plot this as bar plot
#return list() of estimated metric appended to params

for(newlam in lamset){
  lamname <- gsub("\\.", "", as.character(newlam))
  load(paste0("./out/VCpar", lamname, ".Rdata"))
  
  result <- list()
  for(index in 1:length(params)){
    print(index)
    cases <- params[[index]][["True Inc"]]
    temp <- params[[index]][["Temp"]]
    Rt <- params[[index]][["Rt"]]
    
    #total infectiousness
    aH <- 16 ; sH <- 2.7 #shape and rate param of intrinsic incub (gamma dist)
    aM <- 4.3 ; b0 <- 7.9 ; b1 <- 0.21; #shape and rate params of extrinsic incub (gamma dist)
    lamTM <- 1 ; lamTH <- 1
    
    #calc total infectiousness based on simulated cases
    totalinf_td <- c()
    for(s in 36:length(cases)){
      t <- temp[(s-35)]
      sM <- 1/(b0-b1*t)
      w <- dcoga(1:35, shape=c(aM, aH, 1, 1), rate=c(sM, sH, lamTM, lamTH))
      totalinf_td <- c(totalinf_td, w %*% rev(cases[(s-35):(s-1)]))
    }
    
    cases <- cases[36:length(cases)]
    if(index==20){cases[1] <- 5} #cannot start with 0 else omega estimation breaks
    
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
    
    # Smoothed estimates as list of [Rmed, Rhatci, Rmean, qR]
    Rsmooth_td = epiSmoother(Rgrid, m, Rfilt_td$pR, Rfilt_td$pRup, nday, Rfilt_td$pstate, conf)
    #Ismooth_td = recursPredict(Rgrid, Rsmooth_td$qR, Lday_td[tday], Rsmooth_td$Rmean, conf, Iday, Imax)
    print("td estimated")
    
    # Filtered (causal) estimates as list [Ommed, Omci, Ommean, pOm, pOmup, pOmstate]
    Omfilt = epiFilter(Rgrid, m, eta, pR0, nday, Irms[tday], Iday[tday], conf)
    
    # Smoothed estimates as list of [Ommed, Omhatci, Ommean, qOm]
    Omsmooth = epiSmoother(Rgrid, m, Omfilt$pR, Omfilt$pRup, nday, Omfilt$pstate, conf)
    #IsmoothOm = recursPredict(Rgrid, Omsmooth$qR, Irms[tday], Omsmooth$Rmean, conf, Iday, Imax)
    print("om estimated")
    
    # Prob of R > 1 (resurgence)
    # pOm1 = rep(0, nday); id1 = which(Rgrid >= 1); id1 = id1[1]
    # for (i in 1:nday){
    #   pOm1[i] = sum(Omsmooth$qR[i, id1:m])
    # }
    
    # actual Rt > 1
    # pRt_real <- ifelse(Rt>1, 1, 0)
    
    #########################Not Temp-Dependent########################################
    ################## Input daily dengue data and length of time series ##############
    
    #calc total infectiousness based on simulated cases
    
    w <- dcoga(1:35, shape=c(1, aH, 1, 1), rate=c(0.067, sH, lamTM, lamTH))
    
    cases <- params[[index]][["True Inc"]]
    totalinf <- c()
    for(s in 36:(length(cases))){
      totalinf <- c(totalinf, w %*% rev(cases[(s-35):(s-1)]))
    }
    cases <- cases[36:length(cases)]
    
    Lday <- totalinf
    
    #################################################################################
    ################ Estimate R using EpiFilter ###########################
    
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
    #Ismooth = recursPredict(Rgrid, Rsmooth$qR, Lday[tday], Rsmooth$Rmean, conf, Iday, Imax)
    print("ti estimated")
    
    # # Prob of R > 1 (resurgence)
    # pR1 = rep(0, nday); id1 = which(Rgrid >= 1); id1 = id1[1]
    # for (i in 1:nday){
    #   pR1[i] = sum(Rsmooth$qR[i, id1:m]) #sum of probs when R>1 for each timepoint i
    # }
    
    result[[index]] <- list("Omsmooth"=Omsmooth, "Rsmooth"=Rsmooth, "Rsmooth_td"=Rsmooth_td)
  }
  
  om_score <- c()
  td_score <- c()
  ti_score <- c()
  #calc agreement with true rt
  for(i in 1:length(result)){
    Rt <- params[[i]]$Rt
    Omsmooth <- result[[i]]$Omsmooth
    Rsmooth <- result[[i]]$Rsmooth
    Rsmooth_td <- result[[i]]$Rsmooth_td
    
    ##omega
    #R>1 sig
    OmOver1 <- Omsmooth$Rci[1,31:200] > 1 & Omsmooth$Rci[2,31:200] > 1 & Rt[31:200]>1 #17
    #R<=1 sig
    OmUnder1 <- Omsmooth$Rci[1,31:200] <= 1 & Omsmooth$Rci[2,31:200] <= 1 & Rt[31:200] <=1 #64
    om_score <- c(om_score, (sum(OmOver1) + sum(OmUnder1))/170)
    
    ##td-Rt
    #R>1 sig
    tdOver1 <- Rsmooth_td$Rci[1,31:200] > 1 & Rsmooth_td$Rci[2,31:200] > 1 & Rt[31:200]>1 #66
    #R<=1 sig
    tdUnder1 <- Rsmooth_td$Rci[1,31:200] <= 1 & Rsmooth_td$Rci[2,31:200] <= 1 & Rt[31:200] <=1 #63
    td_score <- c(td_score, (sum(tdOver1) + sum(tdUnder1))/170)
    
    ##ti-Rt
    #R>1 sig
    tiOver1 <- Rsmooth$Rci[1,31:200] > 1 & Rsmooth$Rci[2,31:200] > 1 & Rt[31:200]>1 #78
    #R<=1 sig
    tiUnder1 <- Rsmooth$Rci[1,31:200] <= 1 & Rsmooth$Rci[2,31:200] <= 1 & Rt[31:200] <=1 #38
    ti_score <- c(ti_score, (sum(tiOver1) + sum(tiUnder1))/170)
    
  }
  
  vars$om_score <- om_score
  vars$td_score <- td_score
  vars$ti_score <- ti_score
  
  om_score <- c()
  td_score <- c()
  ti_score <- c()
  #calc agreement with true rt
  for(i in 1:length(result)){
    Rt <- params[[i]]$Rt
    Omsmooth <- result[[i]]$Omsmooth
    Rsmooth <- result[[i]]$Rsmooth
    Rsmooth_td <- result[[i]]$Rsmooth_td
    
    ##omega
    #R>1 sig
    OmOver1 <- Omsmooth$Rmean[31:200] > 1 & Rt[31:200]>1
    #R<=1 sig
    OmUnder1 <- Omsmooth$Rmean[31:200] <= 1 & Rt[31:200] <=1
    om_score <- c(om_score, (sum(OmOver1) + sum(OmUnder1))/170)
    
    ##td-Rt
    #R>1 sig
    tdOver1 <- Rsmooth_td$Rmean[31:200] > 1 & Rt[31:200]>1 #66
    #R<=1 sig
    tdUnder1 <- Rsmooth_td$Rmean[31:200] <= 1 & Rt[31:200] <=1
    td_score <- c(td_score, (sum(tdOver1) + sum(tdUnder1))/170)
    
    ##ti-Rt
    #R>1 sig
    tiOver1 <- Rsmooth$Rmean[31:200] > 1 & Rt[31:200]>1 
    #R<=1 sig
    tiUnder1 <- Rsmooth$Rci[1,31:200] <= 1 & Rt[31:200] <=1
    ti_score <- c(ti_score, (sum(tiOver1) + sum(tiUnder1))/170)
    
  }
  
  vars$om_score_ns <- om_score
  vars$td_score_ns <- td_score
  vars$ti_score_ns <- ti_score
  
  varsVC <- vars
  paramsVC <- params
  resultVC <- result
  
  save(list=c("varsVC", "paramsVC", "resultVC"), file=paste0("./out/VC", lamname, "_sim_output.Rdata"))
  
  write.csv(file=paste0("./out/tables/vars", lamname, ".csv"), varsVC, row.names=F)  
}

lamnames <- c("4", "3", "033", "025")

for(lamname in lamnames){
  files <- paste0("./out/VC", lamname, "_sim_output.Rdata")
  df <- as.data.frame(matrix(nrow=54, ncol=3))
  load(files)
  
  # Setup grid [Rmin Rmax] noise eta and CIs conf
  Rmin = 0.01; Rmax = 10; eta = 0.1
  # Uniform prior over grid of size m
  m = 1000; pR0 = (1/m)*rep(1, m)
  
  # Delimited grid defining space of R
  Rgrid = seq(Rmin, Rmax, length.out = m)
  nday = 200
  m = 1000
  
  for(s in 1:54){
    pOm1 = rep(0, nday); id1 = which(Rgrid >= 1); id1 = id1[1]
    pR1 = pOm1; pR1_td = pOm1
    for (i in 1:nday){
      pOm1[i] = sum(result[[s]]$Omsmooth$qR[i, id1:m])
      pR1[i] = sum(result[[s]]$Rsmooth$qR[i, id1:m])
      pR1_td[i] = sum(result[[s]]$Rsmooth_td$qR[i, id1:m])
    }
    
    
    roc_td <- roc(as.numeric(params[[s]]$Rt > 1)[36:nday], pR1_td[36:nday])
    if(is.na(pOm1[1])){print(s)} else {roc_om <- roc(as.numeric(params[[s]]$Rt > 1)[36:nday], pOm1[36:nday])}
    roc_ti <- roc(as.numeric(params[[s]]$Rt > 1)[36:nday], pR1[36:nday])
    
    df[s, ] <- c(roc_td$auc, roc_om$auc, roc_ti$auc)
  }
  
  colnames(df) <- c("td-Rt", "Omega", "ti-Rt")
  df$OmDiff <- df$`td-Rt` - df$Omega
  df$tiDiff <- df$`td-Rt` - df$`ti-Rt`
  
  df <- round(df, 3)
  write.csv(df, paste0("./out/tables/VC", lamname, "aucroc.csv"), row.names=F)
}



