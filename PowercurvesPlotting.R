library(ggplot2)

load("powercurves1.RData")
dfplot <- data.frame(t = (0:20)/20,
                     means = c(powercurves[,2,]),
                     lowers = c(powercurves[,1,]),
                     uppers = c(powercurves[,3,]),
                     type = as.factor(rep(c("T1","T2","T3","T4","T5"),each=21)))
ggplot(dfplot,aes(x=t,y=means,color=type,fill=type,group=type))+
  geom_ribbon(aes(x=t,y=means,ymin=lowers,ymax=uppers),lwd=0,alpha=0.2)+
  geom_line()+
  labs(title="Diffuse noise, dimension 3",fill="Statistic",color="Statistic")+xlab("t")+ylab("Power")
rm(list=ls())

load("powercurves2.RData")
dfplot <- data.frame(t = (0:20)/20,
                     means = c(powercurves[,2,]),
                     lowers = c(powercurves[,1,]),
                     uppers = c(powercurves[,3,]),
                     type = as.factor(rep(c("T1","T2","T3","T4","T5"),each=21)))
ggplot(dfplot,aes(x=t,y=means,color=type,fill=type,group=type))+
  geom_ribbon(aes(x=t,y=means,ymin=lowers,ymax=uppers),lwd=0,alpha=0.2)+
  geom_line()+
  labs(title="Rank-1 with theta=0, dimension 3",fill="Statistic",color="Statistic")+xlab("t")+ylab("Power")
rm(list=ls())

load("powercurves3.RData")
dfplot <- data.frame(t = (0:20)/20,
                     means = c(powercurves[,2,]),
                     lowers = c(powercurves[,1,]),
                     uppers = c(powercurves[,3,]),
                     type = as.factor(rep(c("T1","T2","T3","T4","T5"),each=21)))
ggplot(dfplot,aes(x=t,y=means,color=type,fill=type,group=type))+
  geom_ribbon(aes(x=t,y=means,ymin=lowers,ymax=uppers),lwd=0,alpha=0.2)+
  geom_line()+
  labs(title="Rank-1 with theta=pi/4, dimension 3",fill="Statistic",color="Statistic")+xlab("t")+ylab("Power")
rm(list=ls())

load("powercurves4.RData")
dfplot <- data.frame(t = (0:20)/20,
                     means = c(powercurves[,2,]),
                     lowers = c(powercurves[,1,]),
                     uppers = c(powercurves[,3,]),
                     type = as.factor(rep(c("T1","T2","T3","T4","T5"),each=21)))
ggplot(dfplot,aes(x=t,y=means,color=type,fill=type,group=type))+
  geom_ribbon(aes(x=t,y=means,ymin=lowers,ymax=uppers),lwd=0,alpha=0.2)+
  geom_line()+
  labs(title="Rank-1 with theta=pi/2, dimension 3",fill="Statistic",color="Statistic")+xlab("t")+ylab("Power")
rm(list=ls())

load("powercurves5.RData")
dfplot <- data.frame(t = (0:20)/20,
                     means = c(powercurves[,2,]),
                     lowers = c(powercurves[,1,]),
                     uppers = c(powercurves[,3,]),
                     type = as.factor(rep(c("T1","T2","T3","T4","T5"),each=21)))
ggplot(dfplot,aes(x=t,y=means,color=type,fill=type,group=type))+
  geom_ribbon(aes(x=t,y=means,ymin=lowers,ymax=uppers),lwd=0,alpha=0.2)+
  geom_line()+
  labs(title="Salt-and-pepper, dimension 3",fill="Statistic",color="Statistic")+xlab("t")+ylab("Power")
rm(list=ls())



