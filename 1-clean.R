library(EpiEstim)
library(caTools)
library(coga)
library(tidyverse)

setwd("/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/Rt/")
folres = paste0("./out") # Folder path for results

##aggregate national daily dengue cases
cases <- read.csv("./data/sgdenguedaily.csv")
natlcases <- rowSums(cases[2:ncol(cases)])
tscases <- ts(data=natlcases, start=c(2012, 2), frequency=365)
cases <- tscases

save(cases, file="./data/cases.Rdata")
load("./data/cases.Rdata")

####clean daily sg temperature
temp <- read.csv("./data/all_districts_weather_data.csv")
temp$Mean.Temperature...C. <- as.numeric(temp$Mean.Temperature...C.)
temp <- temp %>% 
  group_by(Year, Month, Day) %>%
  summarise(MeanTemp=mean(`Mean.Temperature...C.`, na.rm=T))
temp <- temp[-1,] #dengue data starts from 2 jan 2012
temp <- temp[1:length(cases),]
temp <- temp$MeanTemp
save(temp, file="./data/temp.Rdata")
load("./data/temp.Rdata")

##parameters
aH <- 16 ; sH <- 2.7 #shape and rate param of intrinsic incub (gamma dist)
aM <- 4.3 ; b0 <- 7.9 ; b1 <- 0.21; #shape and rate params of extrinsic incub (gamma dist)
lamTM <- 1 ; lamTH <- 1
mgt_td <- c() #modal gen time

totalinf_td <- c()
for(s in 36:length(cases)){
  print(s)
  t <- temp[s]
  sM <- 1/(b0-b1*t)
  w_td <- dcoga(1:35, shape=c(aM, aH, 1, 1), rate=c(sM, sH, lamTM, lamTH))
  mgt_td <- c(mgt_td, which.max(w_td))
  totalinf_td <- c(totalinf_td, w_td %*% rev(cases[(s-35):(s-1)]))
}

totalinf <- c()
w <- dcoga(1:35, shape=c(1, aH, 1, 1), rate=c(0.23, sH, lamTM, lamTH))
mgt <- which.max(w)
for(s in 36:length(cases)){
  print(s)
  totalinf <- c(totalinf, w %*% rev(cases[(s-35):(s-1)]))
}

cases <- cases[36:length(cases)]
save(list=c("totalinf", "totalinf_td", "cases"), file="./rtdata.Rdata")
save(list=c("mgt", "mgt_td"), file="./out/gentime.Rdata")

