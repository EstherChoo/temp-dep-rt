library(ggplot2)
library(tidyverse)
library(cowplot)
library(coga)

setwd("C:/Users/esthe/OneDrive - Nanyang Technological University/Rt/")
setwd("/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/Rt/")

######################plot 1: gen int dist ######################
w_u <- function(){ ###serial interval dist for temp-indep
  aH <- 16 ; sH <- 2.7 #shape and rate param of intrinsic incub (gamma dist)
  lamTM <- 1 ; lamTH <- 1
  
  w <- dcoga(seq(1,35,0.5), shape=c(1, aH, 1, 1), rate=c(0.23, sH, lamTM, lamTH))
  return(w)
}
seq(1,35,0.5)[which.max(w)]

w_u_td <- function(tp){
  aH <- 16 ; sH <- 2.7 #shape and rate param of intrinsic incub (gamma dist)
  aM <- 4.3 ; b0 <- 7.9 ; b1 <- 0.21; #shape and rate params of extrinsic incub (gamma dist)
  lamTM <- 1 ; lamTH <- 1
  
  sM <- 1/(b0-b1*tp)
  w <- dcoga(seq(1,35,0.5), shape=c(aM, aH, 1, 1), rate=c(sM, sH, lamTM, lamTH))
  
  return(w)
}
seq(1,35,0.5)[which.max(w)]

serint <- data.frame(matrix(ncol=4, nrow=69))
colnames(serint) <- c("Temp-Independent", "Temp-Dep (33\u00b0C)", "Temp-Dep (35\u00b0C)", "Temp-Dep (37\u00b0C)")
serint$`Temp-Dep (33°C)` <- w_u_td(33)
serint$`Temp-Dep (35°C)` <- w_u_td(35)
serint$`Temp-Dep (37°C)` <- w_u_td(37)
serint$`Temp-Independent` <- w_u()
serint$n <- seq(1,35,0.5)
serint <- pivot_longer(serint, cols=`Temp-Independent`:`Temp-Dep (37°C)`, names_to="Type", values_to="Temperature")

ggplot(serint) +
  geom_line(aes(x=n, y=Temperature, color=Type)) +
  labs(x="Days", y="Serial Interval Distribution") +
  scale_colour_manual(values=c("Temp-Dep (33°C)" = "cadetblue3", "Temp-Dep (35°C)" = "dodgerblue", 
                               "Temp-Dep (37°C)" = "navy", "Temp-Independent" = "magenta4")) +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside= c(0.8,0.8))
  

##################### plot 2: agreement over time ######################
load("./out/real_data_results.Rdata")
load("./data/cases.Rdata")
cases <- cases[36:4564]
dates <- seq(from=as.Date("2012-01-03"), by="day", length.out=length(cases))
year <- data.frame(dates) %>%
  mutate(year = year(dates)) %>%
  group_by(year) %>%
  summarise(first_date = min(dates)) %>%
  arrange(year)

p0 <- ggplot() +
  geom_line(aes(x=dates, y=cases)) +
  scale_x_continuous(labels=year$year, breaks=c(year$first_date)) +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(hjust = -1)) +
  labs(title="A. Dengue Cases in Singapore", y="Cases", x="")

##Omega vs unspecified Rt##
Rt_sig <- (Rsmooth$Rci[1,2:nday] < 1 & Rsmooth$Rci[2,2:nday] < 1) | (Rsmooth$Rci[1,2:nday] > 1 & Rsmooth$Rci[2,2:nday] > 1)
Om_sig <- (Omsmooth$Rci[1,2:nday] < 1 & Omsmooth$Rci[2,2:nday] < 1) | (Omsmooth$Rci[1,2:nday] > 1 & Omsmooth$Rci[2,2:nday] > 1)
td_Rt_sig <- (Rsmooth_td$Rci[1,2:nday] < 1 & Rsmooth_td$Rci[2,2:nday] < 1) | (Rsmooth_td$Rci[1,2:nday] > 1 & Rsmooth_td$Rci[2,2:nday] > 1)

Om_Rt <- (Rsmooth$Rmean < 1 & Omsmooth$Rmean < 1) | (Rsmooth$Rmean > 1 & Omsmooth$Rmean > 1) 
Om_Rt <- Rt_sig & Om_sig & Om_Rt[2:nday] 
Om_Rt_ns <- !Rt_sig & !Om_sig
Om_Rt <- Om_Rt | Om_Rt_ns
mean(Om_Rt) #48.0%
rle1 <- rle(Om_Rt)
from1 <- as.Date(c())
to1 <- from1
for(i in 1:length(rle1$lengths)){
  if(rle1$values[i]==F) next
  start <- sum(rle1$lengths[1:(i-1)]) + 1
  end <- sum(rle1$lengths[1:i])
  from1 <- c(from1, dates[start])
  to1 <- c(to1, dates[end])
}

p1 <- ggplot() +
  geom_ribbon(aes(ymin=Rsmooth$Rci[1,], ymax=Rsmooth$Rci[2,], x=dates, fill="ti-Rt"), alpha=0.5) +
  geom_ribbon(aes(ymin=Omsmooth$Rci[1,], ymax=Omsmooth$Rci[2,], x=dates, fill="Omega"), alpha=0.5) +
  geom_rect(aes(xmin=from1, xmax=to1, ymin=0, ymax=3), fill="slategray", alpha=0.3) +
  geom_line(aes(x=dates, y=Rsmooth$Rmean, color="ti-Rt")) +
  geom_line(aes(x=dates, y=Omsmooth$Rmean, color="Omega")) +
  geom_line(aes(x=dates, y=1), linetype="dashed", colour="black") +
  theme_bw() +
  scale_x_continuous(labels=year$year, breaks=c(year$first_date)) +
  scale_color_manual(values=c("ti-Rt" = "navyblue", "Omega"="chartreuse4")) +
  scale_fill_manual(values=c("ti-Rt" = "navyblue", "Omega"="chartreuse4")) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(hjust = -1)) +
  labs(title="B. Angular Reproduction Number (Omega) against ti-Rt", y="Rt", x="")

##Omega vs temp-dependent Rt##
Om_td_Rt <- (Rsmooth_td$Rmean < 1 & Omsmooth$Rmean < 1) | (Rsmooth_td$Rmean > 1 & Omsmooth$Rmean > 1) #92.7% ns
Om_td_Rt <- Om_td_Rt[2:nday] & Om_sig & td_Rt_sig #22.3% s
Om_td_Rt_ns <- !td_Rt_sig & !Om_sig
Om_td_Rt <- Om_td_Rt | Om_td_Rt_ns
rle2 <- rle(Om_td_Rt)
from2 <- as.Date(c())
to2 <- from2
for(i in 1:length(rle2$lengths)){
  if(rle2$value[i]==F) next
  start <- sum(rle2$lengths[1:(i-1)]) + 1
  end <- sum(rle2$lengths[1:i])
  from2 <- c(from2, dates[start])
  to2 <- c(to2, dates[end])
}

p2 <- ggplot() +
  geom_ribbon(aes(ymin=Rsmooth_td$Rci[1,], ymax=Rsmooth_td$Rci[2,], x=dates, fill="td-Rt"), alpha=0.5) +
  geom_ribbon(aes(ymin=Omsmooth$Rci[1,], ymax=Omsmooth$Rci[2,], x=dates, fill="Omega"), alpha=0.5) + 
  geom_rect(aes(xmin=from2, xmax=to2, ymin=0, ymax=3), fill="slategray", alpha=0.3) +
  scale_x_continuous(labels=year$year, breaks=c(year$first_date)) +
  geom_line(aes(x=dates, y=Rsmooth_td$Rmean, color="td-Rt")) +
  geom_line(aes(x=dates, y=Omsmooth$Rmean, color="Omega")) +
  geom_line(aes(x=dates, y=1), linetype="dashed", colour="black") +
  theme_bw() +
  scale_color_manual(values=c("td-Rt" = "darkorchid", "Omega"="chartreuse4")) +
  scale_fill_manual(values=c("td-Rt" = "darkorchid", "Omega"="chartreuse4")) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(hjust = -1)) +
  labs(title="C. Angular Reproduction Number (Omega) against td-Rt", y="Rt", x="")



png("./plots/2-realdataest.png", width=1200, height=830)
print(plot_grid(p0, p1, p2, nrow=3))
dev.off()

#######AUC-ROC plot (real data)#############
# library(pROC)
# td_R1 <- as.numeric([2:nday] >= 0.5)
# R1 <- as.numeric([2:nday] >= 0.5)
# Om1 <- as.numeric(pOm1[2:nday] >= 0.5)
# 
# load("./out/real_data_results.Rdata")
# roc_rt <- roc(pR1, pOm1)
# auc(roc_rt) #0.77
# roc_om <- roc(pR1_td, pOm1)
# auc(roc_om) #0.88
# 
# par(mfrow=c(2,1))
# plot(roc_rt, main="ti-Rt vs Omega")
# plot(roc_om, main="td-Rt vs Omega")

##################### plot 3: array of plots from simulation results ######################
#temp / true rt / td-rt / ti-rt / om
#smooth case+amp1, period 123, temp123
load("./out/par.Rdata")
load("./out/sim_output.Rdata")

letter <-c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
ind <- c(1,2,3,19,20,21,37,38,39)

###temperature###
templist <- list()
for(i in 1:length(ind)){
  if(i==1){
    templist[[i]] <- local({
      i <- i
      ggplot() +
        geom_line(aes(x=1:200, y=params[[ind[i]]][["Temp"]])) +
        theme_classic() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        labs(y="", title="Temperature (\u00B0C)") +
        scale_y_continuous(limits=c(20,33))
    })
  } else {
    templist[[i]] <- local({
      i <- i
      ggplot() +
        geom_line(aes(x=1:200, y=params[[ind[i]]][["Temp"]])) +
        theme_classic() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        labs(y="") +
        scale_y_continuous(limits=c(20,33))
    })
  }
}

tempplot <- plot_grid(plotlist=templist, ncol=1, labels=LETTERS[1:9])

####true Rt#####
rtlist <- list()
for(i in 1:length(ind)){
  if(i==1){
    rtlist[[i]] <- local({
      i <- i
      ggplot() +
        geom_line(aes(x=1:200, y=params[[ind[i]]][["Rt"]])) +
        theme_classic() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        labs(y="", title="True Rt") +
        scale_y_continuous(limits=c(0,2))
    })
  } else {
    rtlist[[i]] <- local({
      i <- i
      ggplot() +
        geom_line(aes(x=1:200, y=params[[ind[i]]][["Rt"]])) +
        theme_classic() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        labs(y="") +
        scale_y_continuous(limits=c(0,2))
    })
  }
}

rtplot <- plot_grid(plotlist=rtlist, ncol=1)

####cases#####
caseslist <- list()
for(i in 1:length(ind)){
  if(i==1){
    caseslist[[i]] <- local({
      i <- i
      ggplot() +
        geom_line(aes(x=1:200, y=params[[ind[i]]][["True Inc"]][36:235])) +
        theme_classic() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        labs(y="", title="Cases") +
        scale_y_continuous(limits=c(0,280))
    })
  } else {
    caseslist[[i]] <- local({
      i <- i
      ggplot() +
        geom_line(aes(x=1:200, y=params[[ind[i]]][["True Inc"]][36:235])) +
        theme_classic() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        labs(y="") +
        scale_y_continuous(limits=c(0,280))
    })
  }
}

casesplot <- plot_grid(plotlist=caseslist, ncol=1)

#####temp-dependent#####
tdlist <- list()
for(i in 1:length(ind)){
  if(i==1){
    tdlist[[i]] <- local({
      i <- i
      ggplot() +
        geom_line(aes(x=1:200, y=params[[ind[i]]][["Rt"]], color="red")) +
        geom_line(aes(x=1:200, y=result[[ind[i]]][["Rsmooth_td"]][["Rmean"]])) +
        geom_ribbon(aes(x=1:200, ymin=result[[ind[i]]][["Rsmooth_td"]][["Rci"]][1,], ymax=result[[ind[i]]][["Rsmooth_td"]][["Rci"]][2,]), alpha=0.3) +
        geom_label(aes(x=150, y=2.5, label=paste0(round(vars[ind[i], "td_score"]*100, 1), "%")), size=4) +
        theme_classic() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              legend.position="none") +
        labs(y="", title="td-Rt") +
        scale_y_continuous(limits=c(0,3))
    })
  } else {
    tdlist[[i]] <- local({
      i <- i
      ggplot() +
        geom_line(aes(x=1:200, y=params[[ind[i]]][["Rt"]], color="red")) +
        geom_line(aes(x=1:200, y=result[[ind[i]]][["Rsmooth_td"]][["Rmean"]])) +
        geom_ribbon(aes(x=1:200, ymin=result[[ind[i]]][["Rsmooth_td"]][["Rci"]][1,], ymax=result[[ind[i]]][["Rsmooth_td"]][["Rci"]][2,]), alpha=0.3) +
        geom_label(aes(x=150, y=2.5, label=paste0(round(vars[ind[i], "td_score"]*100, 1), "%")), size=4) +
        theme_classic() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              legend.position="none") +
        labs(y="") +
        scale_y_continuous(limits=c(0,3))
    })
  }
}

tdplot <- plot_grid(plotlist=tdlist, ncol=1, hjust=0.3)

#####temp-indep#####
tilist <- list()
for(i in 1:length(ind)){
  if(i==1){
    tilist[[i]] <- local({
      i <- i
      ggplot() +
        geom_line(aes(x=1:200, y=params[[ind[i]]][["Rt"]], color="red")) +
        geom_line(aes(x=1:200, y=result[[ind[i]]][["Rsmooth"]][["Rmean"]])) +
        geom_ribbon(aes(x=1:200, ymin=result[[ind[i]]][["Rsmooth"]][["Rci"]][1,], ymax=result[[ind[i]]][["Rsmooth"]][["Rci"]][2,]), alpha=0.3) +
        geom_label(aes(x=150, y=2.5, label=paste0(round(vars[ind[i], "ti_score"]*100, 1), "%")), size=4) +
        theme_classic() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              legend.position="none") +
        labs(y="", title="ti-Rt") +
        scale_y_continuous(limits=c(0,3))
    })
  } else {
    tilist[[i]] <- local({
      i <- i
      ggplot() +
        geom_line(aes(x=1:200, y=params[[ind[i]]][["Rt"]], color="red")) +
        geom_line(aes(x=1:200, y=result[[ind[i]]][["Rsmooth"]][["Rmean"]])) +
        geom_ribbon(aes(x=1:200, ymin=result[[ind[i]]][["Rsmooth"]][["Rci"]][1,], ymax=result[[ind[i]]][["Rsmooth"]][["Rci"]][2,]), alpha=0.3) +
        geom_label(aes(x=150, y=2.5, label=paste0(round(vars[ind[i], "ti_score"]*100, 1), "%")), size=4) +
        theme_classic() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              legend.position="none") +
        labs(y="") +
        scale_y_continuous(limits=c(0,3))
    })
  }
}

tiplot <- plot_grid(plotlist=tilist, ncol=1, hjust=0.3)

#####omega#####
omlist <- list()
for(i in 1:length(ind)){
  if(i==1){
    omlist[[i]] <- local({
      i <- i
      ggplot() +
        geom_line(aes(x=1:200, y=params[[ind[i]]][["Rt"]], color="red")) +
        geom_line(aes(x=1:200, y=result[[ind[i]]][["Omsmooth"]][["Rmean"]])) +
        geom_ribbon(aes(x=1:200, ymin=result[[ind[i]]][["Omsmooth"]][["Rci"]][1,], ymax=result[[ind[i]]][["Omsmooth"]][["Rci"]][2,]), alpha=0.3) +
        geom_label(aes(x=150, y=2.5, label=paste0(round(vars[ind[i], "om_score"]*100, 1), "%")), size=4) +
        theme_classic() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              legend.position="none") +
        labs(y="", title="Omega") +
        scale_y_continuous(limits=c(0,3))
    })
  } else {
    omlist[[i]] <- local({
      i <- i
      ggplot() +
        geom_line(aes(x=1:200, y=params[[ind[i]]][["Rt"]], color="red")) +
        geom_line(aes(x=1:200, y=result[[ind[i]]][["Omsmooth"]][["Rmean"]])) +
        geom_ribbon(aes(x=1:200, ymin=result[[ind[i]]][["Omsmooth"]][["Rci"]][1,], ymax=result[[ind[i]]][["Omsmooth"]][["Rci"]][2,]), alpha=0.3) +
        geom_label(aes(x=150, y=2.5, label=paste0(round(vars[ind[i], "om_score"]*100, 1), "%")), size=4) +
        
        theme_classic() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              legend.position="none") +
        labs(y="") +
        scale_y_continuous(limits=c(0,3))
    })
  }
}

omplot <- plot_grid(plotlist=omlist, ncol=1, hjust=0.3)

finalplot <- plot_grid(tempplot, rtplot, casesplot, tdplot, tiplot, omplot, nrow=1)
finalplot1 <- plot_grid(finalplot, ggplot()+theme_void(), nrow=2, rel_heights = c(1,0.03))


png("./plots/3-simarray.png", res=300, units="in", width=12, height=9)
ggdraw(finalplot1) +
  draw_label("Time (Days)", x = 0.1, y = 0.025, size = 10) +
  draw_label("Time (Days)", x = 0.27, y = 0.025, size = 10) +
  draw_label("Time (Days)", x = 0.44, y = 0.025, size = 10) +
  draw_label("Time (Days)", x = 0.6, y = 0.025, size = 10) +
  draw_label("Time (Days)", x = 0.76, y = 0.025, size = 10) +
  draw_label("Time (Days)", x = 0.93, y = 0.025, size = 10)
dev.off()

 


