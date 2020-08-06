

#CLEAR WORKSPACE
cat("\014")
rm(list=ls())

# LOAD LIBRARIES OF ALL REQUIRED PACKAGES
library(deSolve)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(RColorBrewer)
library(rstudioapi)

#SET WORKING DIRECTORY 
setwd("/Users/sandraadele/Desktop/PLACEMENT")
getwd()


#SET THE DURATION OF THE MODEL  
#MODEL DURATION (in years) 
#MODEL DURATION (in years)
duration <- 1 # years
time.seq <- seq(0, duration*365, by = 1/30) # in months, by 1 day intervals


#MODEL PARAMETERS FOR HUMAN AND MOSQUITO 
#MODEL PARAMETERS HUMAN 

mu_s= 0.08 #case fatality rate for severe malaria 
mu_h=0.012 #human death rate 
mu_bh= 0.01225 #human birth rate
alpha=0.13 #rate of progression from liver to blood stage 
beta=0.09 #rate of progression from blood stage to mild malaria 
delta_one=0.2 #rate of progression from mild malaria to severe malaria 
delta_two=0.38 #proportion of mild malaria that receive treatment 
delta_3=0.06 #rate of recovery from mild disease without treatment 
sigma=0.7 #rate of hospitalization for severe malaria 
gamma=0.8 #rate of recovery after treatment 
tau=0.14 #rate of recovery after hospitalization
epsilon=0.0019 #rate of recovery from infectiousness in clinically immune
omega=0.00055 #rate of loss of immunity 
chi=0.05 #proportion of mild malaria that become severe 
theta=0.11 #proportion of mild cases that receive treatment 
kappa=0.9 #proportion of clinically immune that develop latent infection 

#MODEL PARAMETERS MOSQUITO
mu_bm=0.13 #mosquito birth rate 
mu_m=0.13 #mosquito death rate 
gamma_m=0.091 #rate of progression from exposed to infectious 

#FORCE OF INFECTION FOR HUMAN AND MOSQUITO COMPARTMENTS 
bee=6 #mosquito biting rate 
amp=1.01 #relative amplitude of seasonal forcing
phi=0.98 #seasonal offset parameter
ai=0.05 #per bite probability of transmission from infected mosquito to human
ai_m=0.47 #per bite probability of transmission from infectious human to mosquito
t=1 #placeholder


#MODEL INITIAL CONDITIONS FOR HUMAN AND MOSQUITO POPULATION 
#HUMAN POPULATION
initPh<-201000000 #population size
initLs<-0 #liver stage non-immune
initBs<-0 #blood stage non-infectious non-immune
initIc<- 0 #gametocyte stage mild malaria
initIs<-0 #gametocyte stage severe malaria 
initTr<-0 #treatment
initTh<-0 #hospitalization
initR<-0 #recovered
initLr<-0 #liver stage clinical immunity
initBr<-0 #blood stage clinical immunity
initIr<-0 #gametocyte stage clinical immunity 

#INITIAL SUSCEPTIBLE COMPARTMENT-NO CLINICAL IMMMUNITY 
initSh<-initPh-initLs-initBs-initIc-initIs-initTr-initTh-initR-initLr-initBr-initIr

#MOSQUITO POPULATION
initPm<-201000000 #mosquito population size
initEm<-1 #initial exposed mosquito population
initIm<-4562700 #initial infected mosquitoes

#INITIAL MOSQUITO SUSCEPTIBLE COMPARTMENT 
initSm<-initPm-initEm-initIm 
state<-c(Sh=initSh,Ls=initLs,Bs=initBs, Ic=initIc, Is=initIs, Tr=initTr, Th=initTh, R=initR, Lr=initLr, Br=initBr, Ir=initIr, Sm = initSm, Em=initEm, Im = initIm)


#MALARIA MODEL FUNCTION
malaria<-function(t, state, parameters) 
{
  parameters <- c(
    #MODEL PARAMETERS FOR HUMAN AND MOSQUITO 
    #MODEL PARAMETERS HUMAN 
    
    mu_bh=mu_bh, #human birth rate
    mu_h=mu_h, #human death rate 
    mu_s=mu_s, #case fatality rate for severe malaria 
    alpha=alpha, #rate of progression from liver to blood stage 
    beta=beta, #rate of progression from blood stage to mild malaria 
    delta_one=delta_one, #rate of progression from mild malaria to severe malaria 
    delta_two=delta_two, #proportion of mild malaria that receive treatment 
    delta_3=delta_3, #rate of recovery from mild disease without treatment 
    sigma=sigma, #rate of hospitalization for severe malaria 
    gamma=gamma, #rate of recovery after treatment 
    tau=tau, #rate of recovery after hospitalization
    epsilon=epsilon, #rate of recovery from infectiousness in clinically immune
    omega=omega, #rate of loss of immunity 
    chi=chi, #proportion of mild malaria that become severe 
    theta=theta, #proportion of mild cases that receive treatment 
    kappa=kappa,#proportion of clinically immune that develop latent infection 
    
    #MODEL PARAMETERS MOSQUITO
    mu_bm=mu_bm, #mosquito birth rate 
    mu_m=mu_m, #mosquito death rate 
    gamma_m=gamma_m, #rate of progression from exposed to infectious 
    
    #FORCE OF INFECTION FOR HUMAN AND MOSQUITO COMPARTMENTS 
    bee=bee, #mosquito biting rate 
    amp=amp, #relative amplitude of seasonal forcing
    phi=phi, #seasonal offset parameter
    ai=ai, #per bite probability of transmission from infected mosquito to human
    ai_m=ai_m, #per bite probability of transmission from infectious human to mosquito
    t=t #placeholder
    
  )
  with(as.list(c(state, parameters)), 
       {
         #baseline parameters 
         seas<- 1+amp*cos(2*pi*(t-phi)/300)
         t=1
         
         
         
         # define variables for mosquito and human population 
         Ph <- (Sh+Ls+Bs+Ic+Is+Tr+Th+R+Lr+Br+Ir)
         Pm <- (Sm+Em+Im)
         
         # force of infection
         lambda <- seas*bee*ai*(Im/Pm)
         lambda_m <- seas*bee*ai_m*(Ic/Ph)
         
         #rate of change human and mosquito
         dSh<- mu_bh *Ph + omega*(1-kappa)*R-mu_h*Sh-lambda*Sh
         dLs<-lambda*Sh-alpha*Ls-mu_h*Ls
         dBs<-alpha*Ls-beta*Bs-mu_h*Bs
         dIc<-beta*Bs-chi*(delta_one)*Ic-theta*(delta_two)*Ic-(1-chi)*(1-theta)*delta_3*Ic-mu_h*Ic
         dIs<-chi*(delta_one)*Ic-sigma*Is-(mu_s+mu_h)*Is
         dTr<-theta*(delta_two)*Ic-gamma*Tr-mu_h*Tr
         dTh<-sigma*Is-tau*Th-(mu_s+mu_h)*Th
         dR<-gamma*Tr+tau*Th+(1-chi)*(1-theta)*delta_3*Ic*(t)+epsilon*Ir-omega*(1-kappa)*R-kappa*(lambda)*R-mu_h*R
         dLr<-kappa*(lambda)*R-alpha*Lr-mu_h*Lr
         dBr<-alpha*Lr-beta*Br-mu_h*Br
         dIr<-beta*Br-epsilon*Ir-mu_h*Ir
         dMt<-Is*mu_s + Th*mu_s 
         
         
         dSm<-mu_bm*Pm-mu_m*Sm-lambda_m*Sm
         dEm<-lambda_m*Sm-mu_m*Em-gamma_m*Em
         dIm<-gamma_m*Em-mu_m*Im
         
         list(c(dSh, dLs, dBs, dIc, dIs, dTr, dTh, dR, dLr, dBr, dIr,dSm,dEm,dIm))
       })}

out  <- ode(y = state, times = time.seq, func = malaria, parms = parameters)
pop  <- out[,"Sh"]+out[,"Ls"]+out[,"Bs"]+out[,"Ic"]+out[,"Is"]+out[,"Tr"]+out[,"Th"]+out[,"R"]+out[,"Lr"]+out[,"Br"]+out[,"Ir"]+out[,"Sm"]+out[,"Em"]+out[,"Im"]
time <- out[,"time"]
OUT  <- data.frame(out)

# total human population 
OUT$Ph_tot <- OUT$Sh + OUT$Ls + OUT$Bs + OUT$Ic + OUT$Is + OUT$Tr + OUT$Th + OUT$R +OUT$Lr + OUT$Br + OUT$Ir 

# total mosquito population
OUT$Pm_tot <- OUT$Sm + OUT$Em + OUT$Im
OUT$year <- time/365
# Percentage population prevalence (asymptomatic and symptomatic people in blood stage)
OUT$Prevalence.perc <- ((OUT$Ic + OUT$Is + OUT$Tr + OUT$Th + OUT$Ir)/OUT$Ph_tot)*100 
#  Population prevalence (asymptomatic and symptomatic people in blood stage)
OUT$Prevalence.pop <- (OUT$Ic + OUT$Is + OUT$Tr + OUT$Th + OUT$Ir) 
# yearly reported incidence = everyone in blood stages who is symptomatic and received treatment)
OUT$Incidence <- beta*OUT$Ic + (chi*delta_one)*OUT$Is 
#Proportion of cases with no immunity / partial immunity
OUT$Clinical <- (OUT$Ic + OUT$Is)/(OUT$Ir) 


# function to round up to neaerest hundred
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

#MODEL PLOTS

#Change format of dataframe so that all compartments/states are in one columnn, 
#count in 2nd column, year in 3rd column
cbPalette <- c("light blue", "dark blue", "red", "maroon","green","dark green",
               "brown","yellow","pink","orange")
Shtab<-cbind(cbind(rep('Sh',length(time.seq)),as.numeric(OUT$Sh)),OUT$year)
Lstab<-cbind(cbind(rep('Ls',length(time.seq)),as.numeric(OUT$Ls)),OUT$year)
Bstab<-cbind(cbind(rep('Bs',length(time.seq)),as.numeric(OUT$Bs)),OUT$year)
Ictab<-cbind(cbind(rep('Ic',length(time.seq)),as.numeric(OUT$Ic)),OUT$year)
Istab<-cbind(cbind(rep('Is',length(time.seq)),as.numeric(OUT$Is)),OUT$year)
Trtab<-cbind(cbind(rep('Tr',length(time.seq)),as.numeric(OUT$Tr)),OUT$year)
Thtab<-cbind(cbind(rep('Th',length(time.seq)),as.numeric(OUT$Th)),OUT$year)
Rtab<-cbind(cbind(rep('R',length(time.seq)),as.numeric(OUT$R)),OUT$year)
Lrtab<-cbind(cbind(rep('Lr',length(time.seq)),as.numeric(OUT$Lr)),OUT$year)
Brtab<-cbind(cbind(rep('Br',length(time.seq)),as.numeric(OUT$Br)),OUT$year)
Irtab<-cbind(cbind(rep('Ir',length(time.seq)),as.numeric(OUT$Ir)),OUT$year)

Prev.perc.tab<-cbind(cbind(rep('Prevalence',length(time.seq)),
                           as.numeric(OUT$Prevalence.perc)),OUT$year)
Prev.pop.tab<-cbind(cbind(rep('Prevalence',length(time.seq)),
                          as.numeric(OUT$Prevalence.pop)),OUT$year)
Inctab<-cbind(cbind(rep('Incidence',length(time.seq)),
                    as.numeric(OUT$Incidence)),OUT$year)
MMWInctab<-cbind(cbind(rep('Incidence',length(time.seq)),
                       as.numeric(OUT$Incidence)),OUT$year)
MMWInctab[,2] <- as.integer(MMWInctab[,2]) # proportion of people treated by the MMW
MMWInctab[,2] <- (as.integer(MMWInctab[,2])/initPh) # incidence per 1000 population

# Plot 1: All compartments
dd<-as.data.frame(rbind(Shtab,Lstab,Bstab,Ictab,Istab,Trtab,Trtab,Rtab,Lrtab,Brtab,Irtab))
colnames(dd)<-c('State','Count','year')
dd$Count<-as.numeric(as.character(dd$Count))
dd$year<-as.numeric(as.character(dd$year))
p1 <- ggplot(data=dd,aes(x=year,y=Count,group=State,color=State)) +
  geom_line(size=1) + 
  theme_bw()+
  scale_colour_manual(values=cbPalette,breaks=c("Sh","Ls","Bs","Ic","Is","Tr","Th","R","Lr","Br","Ir"))+
  labs(title = "1Compartmental model of malaria transmission", x="year", 
       y="Population") +
  theme(legend.title=element_blank())  +
  theme(plot.title = element_text(lineheight=.10, face="bold")) +
  theme(legend.key = element_blank())
p1

# Plot 2: All compartments, except susceptible population
dd2<-as.data.frame(rbind(Lstab,Bstab,Ictab,Istab,Trtab,Trtab,Rtab,Lrtab,Brtab,Irtab))
colnames(dd2)<-c('State','Count','year')
dd2$Count<-as.numeric(as.character(dd2$Count))
dd2$year<-as.numeric(as.character(dd2$year))  

p2 <- ggplot(data=dd2,aes(x=year,y=Count,group=State,color=State)) +
  geom_line(size=1) +
  theme_bw()+
  
  scale_colour_manual(values=cbPalette,breaks=c("Ls","Bs","Ic","Is","Tr","Th","R","Lr","Br","Ir")) +
  labs(title = "2 Compartmental model of malaria transmission", 
       x="year", y="population") +
  theme(legend.title=element_blank())  +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  theme(legend.key = element_blank())

p2

# Plot 3: Infected compartments
dd3<-as.data.frame(rbind(Lstab,Bstab,Ictab,Istab,Trtab,Trtab,Rtab,Lrtab,Brtab,Irtab))
colnames(dd3)<-c('State','Count','year')
dd3$Count<-as.numeric(as.character(dd3$Count))
dd3$year<-as.numeric(as.character(dd3$year))  

p3 <- ggplot(data=dd3,aes(x=year,y=Count,group=State,color=State)) +
  geom_line(size=1) +
  theme_bw()+
  
  scale_colour_manual(values=cbPalette,breaks=c("Ls","Bs","Ic","Is","Tr","Th","R","Lr","Br","Ir")) +
  labs(title = "3 Compartmental model of malaria transmission", 
       x="year", y="Population") +
  theme(legend.title=element_blank())  +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  theme(legend.key = element_blank())
p3

# Plot 4: Prevalence
dd4<-as.data.frame(rbind(Prev.perc.tab))
colnames(dd4)<-c('State','Count','year')
dd4$Count<-as.numeric(as.character(dd4$Count))
dd4$year<-as.numeric(as.character(dd4$year))  

p4 <- ggplot(data=dd4,aes(x=year,y=Count, color="orange")) +
  geom_line(size=1) +
  theme_bw()+
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  labs(title = "4 Malaria prevalence", x="year", y="Prevalence (%)") +
  theme(legend.position="none")
p4



# Plot 5: Incidence
MMWInctab <- as.data.frame(rbind(MMWInctab))
colnames(MMWInctab)<-c('State','Count','year')
MMWInctab$Count<-as.numeric(as.character(MMWInctab$Count))
MMWInctab$year<-as.numeric(as.character(MMWInctab$year)) 

p5 <- ggplot(data=MMWInctab,aes(x=year,y=Count,color="orange")) +
  geom_line(size=1) +
  theme_bw()+
  labs(title = "5 Malaria incidence in Nigeria", x="year", 
       y="Population") +
  theme(legend.title=element_blank())  +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  theme(legend.position="none")
p5

MMWInctab$Date <- seq.Date(as.Date("01/01/20","%d/%m/%y"),
                           by = "day", length.out = 1826)
MMWInctab$CountX1000 <- MMWInctab$Count*1000
plot(MMWInctab[1:1862,4],MMWInctab[1:1862,5],type = 'l', 
     col="purple",lwd=2, ylab = "Incidence per 1000 population", xlab="Year",
     main="Main title")

# Plot 6: Phase plot
plot(MMWInctab[,2],Prev.perc.tab[,2], type = 'l', 
     xlab = "Incidence of symptomatic cases per 1000 population", 
     ylab = "Prevalence (%)", main = "6 Phase plot of incidence vs prevalance", 
     col="light blue", lwd =2)



#######Model fit###########
##line below exports data frame to csv file; after I loaded it into excel I modified it there to add columns for month 
write.csv(dd4,"/Users/sandraadele/Desktop/PLACEMENT/dd4.csv", row.names = FALSE)
##now load in new modified dd4 file
dd4_mod <- read.csv("dd4_mod.csv")
##below I summed total incident cases for each month
month_2020 <- aggregate(dd4_mod$incidence, by=list(Category=dd4_mod$month), FUN=sum)
month_2020
##export summed values by month to csv
write.csv(month_2020,"/Users/sandraadele/Desktop/PLACEMENT/month2020.csv", row.names = FALSE)
cori <- read.csv("cori.csv")
cori$incidence <- (cori$denom/201000000)*1000
##now plot 2017 datapoints
plot(cori[1:12,5],cori[1:12,18], type="l",
     col="purple",lwd=2, ylab = "Incidence per 1000 population", xlab="Month",
     main="Yearly malaria trends (baseline)",xaxt='n')
xtick<-seq(1, 12, by=1)
axis(side=1, at=xtick, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                  "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))
par(new=T)
##2018 datapoints
plot(cori[13:24,5],cori[13:24,18],type="l", ylab = "", xlab="",
     col="green",lwd=2, 
     main="", xaxt='n', yaxt='n')
par(new=T)
##2019 datapoints
plot(cori[25:36,5],cori[25:36,18],type="l",ylab = "", xlab="",
     col="red",lwd=2, 
     main="", xaxt='n', yaxt='n')
par(new=T)
##2020 data
plot(cori[37:48,5],cori[37:48,18],,type="l",ylab = "", xlab="",
     col="black",lwd=2, 
     main="", xaxt='n', yaxt='n')
##make legend 
legend("bottomright", inset=.05,
       c("2017","2018", "2019", "2020"), fill=c("purple","green", "red", "black"), horiz=TRUE)

########### NOT TO BE RUN####################################################







##model fit  
dd4$Date <- seq.Date(as.Date("01/01/15","%d/%m/%y"),by = "day", length.out = 3651)
plot(dd4[1463:1827,4],dd4[1463:1827,2],type = 'l',
     main="Modelfit" ,
     col="black",lwd=2, ylim=c(0,5), xlab="Month", ylab ="Prevalence %")
points(dd4[2709,4],0,col= "red", pch = 19) #MK 0-2.84, ST 0.53-1.68
arrows(dd4[2709,4],0,dd4[2709,4],0.67,length=0.05, angle=90, code=3,col="red")
points(dd4[2831,4],2.84,col= "red", pch = 19)
arrows(dd4[2831,4],1.59,dd4[2831,4],5,length=0.05, angle=90, code=3,col="red")
legend("bottomleft", inset=.05,
       c("Model   ","Data   "), fill=c("purple","green"), horiz=TRUE)

dd4$Date <- seq.Date(as.Date("01/01/15","%d/%m/%y"),by = "day", length.out = 3651)
plot(dd4[1463:1827,4],dd4[1463:1827,2],type = 'l',
     main="Modelfit",
     col="black",lwd=3, ylim=c(0,3.7), xlab="Month", ylab ="Prevalence %")
points(dd4[1500,4],0.53,col= "red", pch = 19) #MK 0-2.84, ST 0.53-1.68
arrows(dd4[1500,4],0.18,dd4[1500,4],1.56,length=0.05, angle=90, code=3,col="red")
points(dd4[1731,4],1.68,col= "red", pch = 19)
arrows(dd4[1731,4],0.77,dd4[1731,4],3.61,length=0.05, angle=90, code=3,col="red")
legend("bottomleft", inset=.05,
       c("Model   ","Data   "), fill=c("blue","orange"), horiz=TRUE)

### Mortality 
initMt<-0
Mttab<-cbind(cbind(rep('Mt',length(time.seq)),as.numeric(OUT$Mt)),OUT$year)
plot(initMt)
month_2020
dd4
