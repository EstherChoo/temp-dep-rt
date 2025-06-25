library(ggplot2)
library(cowplot)

##################### prediction plot ######################
#for each rt, real vs pred + ci + mape % label at side
load("./out/real_data_results.Rdata")
load("./data/cases.Rdata")
cases <- cases[37:4564]
dates <- seq(from=as.Date("2012-01-04"), by="day", length.out=length(cases))
year <- data.frame(dates) %>%
  mutate(year = year(dates)) %>%
  group_by(year) %>%
  summarise(first_date = min(dates)) %>%
  arrange(year)

mape <- function(real, pred){
  round(mean(abs(pred-real)/pmax(real,0.1)) * 100,1)
}


mape(cases, Ismooth$pred)
mape(cases, IsmoothOm$pred)

td <- ggplot() +
  geom_point(aes(y=cases, x=dates, color="Real")) +
  geom_line(aes(y=Ismooth_td$pred, x=dates, color="Predictions")) +
  geom_ribbon(aes(ymin=Ismooth_td$predci[1,], ymax=Ismooth_td$predci[2,], x=dates),alpha=0.7, fill="tomato") +
  geom_text(aes(x=dates[600],y=200,label=paste0("MAPE: ", mape(cases, Ismooth_td$pred), "%"))) +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(hjust = -1)) +
  scale_color_manual(values=c("tomato", "black"))+
  scale_x_continuous(labels=year$year, breaks=c(year$first_date)) +
  labs(x="",y="Cases",title="Temperature-dependent Reproduction Number")

ti <- ggplot() +
  geom_point(aes(y=cases, x=dates, color="Real")) +
  geom_line(aes(y=Ismooth$pred, x=dates, color="Predictions")) +
  geom_ribbon(aes(ymin=Ismooth$predci[1,], ymax=Ismooth$predci[2,], x=dates),alpha=0.7, fill="pink") +
  geom_text(aes(x=dates[600],y=200,label=paste0("MAPE: ", mape(cases, Ismooth$pred), "%"))) +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(hjust = -1)) +
  scale_color_manual(values=c("pink", "black"))+
  scale_x_continuous(labels=year$year, breaks=c(year$first_date)) +
  labs(x="",y="Cases",title="Temperature-independent Reproduction Number")

om <-ggplot() +
  geom_point(aes(y=cases, x=dates, color="Real")) +
  geom_line(aes(y=IsmoothOm$pred, x=dates, color="Predictions")) +
  geom_ribbon(aes(ymin=IsmoothOm$predci[1,], ymax=IsmoothOm$predci[2,], x=dates),alpha=0.7, fill="goldenrod") +
  geom_text(aes(x=dates[600],y=200,label=paste0("MAPE: ", mape(cases, IsmoothOm$pred), "%"))) +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(hjust = -1)) +
  scale_color_manual(values=c("goldenrod", "black"))+
  scale_x_continuous(labels=year$year, breaks=c(year$first_date)) +
  labs(x="",y="Cases",title="Angular Reproduction Number")

png("./plots/s1_pred.png", width=1200,height=800)
plot_grid(td, ti, om, ncol=1)
dev.off()

##################### temp plot ######################
load("./out/par.Rdata")

ggplot() +
  geom_line(aes(x=1:200, y=params[[1]][["Temp"]], colour="Singapore Data"),linewidth=1) +
  geom_line(aes(x=1:200, y=params[[19]][["Temp"]], colour="More Variation")) +
  geom_line(aes(x=1:200, y=params[[37]][["Temp"]], colour="Most Variation")) +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank()) +
  labs(x="Days", y="Temperature (\u00B0C)") +
  scale_color_manual(values=c("purple3", "forestgreen", "tomato"))

##################### array of plots from simulation results ######################
#temp / true rt / td-rt / ti-rt / om
#smooth case+amp1, period 123, temp123
load("./out/par.Rdata")
load("./out/sim_output.Rdata")

ind <- 1:54

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

tempplot <- plot_grid(plotlist=templist, ncol=1, labels=1:54, label_size=11)

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
        scale_y_continuous(limits=c(0,4))
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
        scale_y_continuous(limits=c(0,4))
    })
  }
}

rtplot <- plot_grid(plotlist=rtlist, ncol=1, labels=1:54, label_size=11)

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
        labs(y="", title="Cases")
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
        labs(y="")
    })
  }
}

casesplot <- plot_grid(plotlist=caseslist, ncol=1, labels=1:54, label_size=11)

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
        scale_y_continuous(limits=c(0,4))
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
        scale_y_continuous(limits=c(0,4))
    })
  }
}

tdplot <- plot_grid(plotlist=tdlist, ncol=1, labels=1:54, label_size=11)

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
        scale_y_continuous(limits=c(0,4))
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
        scale_y_continuous(limits=c(0,4))
    })
  }
}

tiplot <- plot_grid(plotlist=tilist, ncol=1, labels=1:54, label_size=11)

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
        scale_y_continuous(limits=c(0,4))
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
        scale_y_continuous(limits=c(0,4))
    })
  }
}

omplot <- plot_grid(plotlist=omlist, ncol=1, labels=1:54, label_size=11)

finalplot <- plot_grid(tempplot, rtplot, casesplot, tdplot, tiplot, omplot, nrow=1)

pdf("./plots/simarrayfull.pdf", width=11.5, height=54)
print(finalplot)
dev.off()

#######AUC-ROC plot (sim data)#############
# 3 x 54?
library(pROC)
load("./out/sim_output.Rdata")

#load in auc value

# Setup grid [Rmin Rmax] noise eta and CIs conf
Rmin = 0.01; Rmax = 10; eta = 0.1
# Uniform prior over grid of size m
m = 1000; pR0 = (1/m)*rep(1, m)

# Delimited grid defining space of R
Rgrid = seq(Rmin, Rmax, length.out = m)
nday = 200
m = 1000
aucroc_list <-  list()
for(s in 1:54){
  aucroc_list[[s]] <- list()
  pOm1 = rep(0, nday); id1 = which(Rgrid >= 1); id1 = id1[1]
  pR1 = pOm1; pR1_td = pOm1
  for (i in 1:nday){
    pOm1[i] = sum(result[[s]]$Omsmooth$qR[i, id1:m])
    pR1[i] = sum(result[[s]]$Rsmooth$qR[i, id1:m])
    pR1_td[i] = sum(result[[s]]$Rsmooth_td$qR[i, id1:m])
  }
  
  aucroc_list[[s]]$td_R1 <- pR1_td[36:nday] 
  aucroc_list[[s]]$R1 <- pR1[36:nday]
  aucroc_list[[s]]$Om1 <- pOm1[36:nday] 
  aucroc_list[[s]]$real_R1 <- as.numeric(params[[s]]$Rt > 1)[36:nday] 
}

letter <-c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))

pdf("./plots/test-sim_auc.pdf", width=11, height=7, paper="a4r")
par(mfrow=c(4,3))


for(s in 1:54){
  roc_td <- roc(aucroc_list[[s]]$real_R1, aucroc_list[[s]]$td_R1)
  roc_om <- roc(aucroc_list[[s]]$real_R1, aucroc_list[[s]]$Om1)
  roc_ti <- roc(aucroc_list[[s]]$real_R1, aucroc_list[[s]]$R1)
  plot(roc_ti, main=paste(s), col="navyblue")
  text(x=0, y=0.6,labels=round(roc_ti$auc, 3), cex=1.2)
  plot(roc_om, col="chartreuse4")
  text(x=0, y=0.6,labels=round(roc_om$auc, 3), cex=1.2)
  plot(roc_td, col="magenta4")
  text(x=0, y=0.6,labels=round(roc_td$auc, 3), cex=1.2)
}

dev.off()

####AUC-ROC table###
library(pROC)

files <- c("./out/sim_output.Rdata")
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
  roc_om <- roc(as.numeric(params[[s]]$Rt > 1)[36:nday], pOm1[36:nday])
  roc_ti <- roc(as.numeric(params[[s]]$Rt > 1)[36:nday], pR1[36:nday])

  df[s, ] <- c(roc_td$auc, roc_om$auc, roc_ti$auc)
}

colnames(df) <- c("td-Rt", "Omega", "ti-Rt")
df$OmDiff <- df$`td-Rt` - df$Omega
df$tiDiff <- df$`td-Rt` - df$`ti-Rt`

df <- round(df, 3)
write.csv(df, "./out/tables/aucroc.csv", row.names=F)


###edit table
csv <- read.csv("./out/tables/vars4.csv")
csv <- csv[1:7]
csv <- arrange(csv, smooth)
colnames(csv) <- c("Period", "Amplitude", "Smooth", "Temp", "Omega", "td-Rt", "ti-Rt")
csv <- csv[-3]
csv[4:6] <- round(csv[4:6], 3)
csv$Temp <- gsub("temp1", "SG Data", csv$Temp)
csv$Temp <- gsub("temp2", "More Variance", csv$Temp)
csv$Temp <- gsub("temp3", "Most Variance", csv$Temp)
csv2 <- csv[28:54,]
csv <- csv[1:27,] #smooth on left, noisy on right
csv <- cbind(csv, csv2)
write.csv(csv, "./out/tables/vars4.csv", row.names=F)



###gen time with real data
load("./out/gentime.Rdata")
ggplot() +
  geom_line(aes(x=1:4529, y=mgt, color="Temperature-Independent"), linewidth=1) +
  geom_line(aes(x=1:4529, y=mgt_td, color="Temperature-Dependent")) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  labs(y="Mean Generation Time", x="Time (Days)") +
  scale_color_manual(values=c("navy", "black"))

auc <- read.csv("./out/tables/aucroc.csv")[1:3]
auc <- auc - 0.9
auc$ind <- 1:54
auc <- pivot_longer(auc, cols=1:3, names_to="Rt")
ggplot(auc) +
  facet_wrap(~ind) +
  geom_col(aes(x=Rt, y=value))
