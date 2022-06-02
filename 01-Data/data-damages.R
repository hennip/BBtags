library(runjags)
require(rjags)
require(sp)
require(openxlsx)
library(readxl)
library(tidyverse)
library(lubridate)

# gear n_tagged n_return      p
# <dbl>    <int>    <dbl>  <dbl>
#  1     1      201       24 0.119 #Kaukalo
#  2     2      211       35 0.166 #Sukka

dat2=read_xlsx("C:/Users/03080932/OneDrive - Valtion/Projects/BBTags/dat/orig/Vauriot_Data_2021_ver2.xlsx")%>%
  filter(pituus>60) #Jätetään pois 1 mahd. 1SW kala kuten JAGS-mallissa

t1<-dat2%>%group_by(rysätyyppi,suomuvaurio)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=ifelse(suomuvaurio==1,1,0), type="suomuvaurio1")

t2<-dat2%>%group_by(rysätyyppi,suomuvaurio)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=ifelse(suomuvaurio==2,1,0), type="suomuvaurio2")

t3<-dat2%>%group_by(rysätyyppi,silmävaurio)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=silmävaurio, type="silmävaurio")

t4<-dat2%>%group_by(rysätyyppi,`haavoja/viiltoja/naarmuja`)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=`haavoja/viiltoja/naarmuja`, type="HaaViiNaarm")

t5<-dat2%>%group_by(rysätyyppi,havaspyydys)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=havaspyydys, type="havaspyydys")

t6<-dat2%>%group_by(rysätyyppi,suuvaurio)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=suuvaurio, type="suuvaurio")

t7<-dat2%>%group_by(rysätyyppi,Evävaurio)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=Evävaurio, type="Evavaurio")


# Suomuvauriot
# t30<-dat2%>%
#   mutate(suomuvaurio12=ifelse(suomuvaurio==2, 1, suomuvaurio))%>%
#   group_by(rysätyyppi,suomuvaurio12)%>%
#   summarise(n=n())%>%
#   mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
#   mutate(p=n/Ntot, x=suomuvaurio12, type="suomuvaurio_1tai2")



dfd<-full_join(t1,t2)%>%
  full_join(t3)%>%
  full_join(t4)%>%
  full_join(t5)%>%
  full_join(t6)%>%
   full_join(t7)%>%
  filter(x==1)%>%
  select(rysätyyppi, n, Ntot, p,x, type)%>%
  arrange(x,rysätyyppi)
#View(dfd)


dfd%>%filter(rysätyyppi=="Kaukalo")%>%select(n, type)
dfd%>%filter(rysätyyppi=="Sukka")%>%select(n, type)


Mdamage<-"
model{  
for(i in 1:7){
  # Kaukalo
  
  xK[i]~dbin(pK[i],NK)
  pK[i]~dbeta(1,1)
  
  # Sukka
  xS[i]~dbin(pS[i],NS)
  pS[i]~dbeta(1,1)

  prob[i]<-step(pK[i]-pS[i])
}



}"

cat(Mdamage,file="Mdamage.txt")

data=list(xK=c(134,35,4,15,3,11,141),
           xS=c(108,87,1,10,0,3,186),
           NK=201,NS=211)

 run0 <- run.jags(Mdamage,
                  monitor= c("prob", "pK", "pS"),data=data, #inits = inits,
                  n.chains = 2, method = 'parallel', thin=10, burnin = 1000,
                  modules = "mix",keep.jags.files=F,sample =10000, adapt = 1000,
                  progress.bar=TRUE)
 
 run<-run0
 
 round(summary(run),3)

chains<-as.mcmc(run) 
round(t(hdi(chains, credMass=0.9)),4)

par(mfrow=c(3,3))
plot(density(chains[,"pK[1]"]), main="A", lwd=2, col="red", xlim=c(0,1), ylim=c(0,15), xlab="Probability")
lines(density(chains[,"pS[1]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[2]"]), main="B", ylim=c(0,15), lwd=2, col="red", xlim=c(0,1), xlab="Probability")
lines(density(chains[,"pS[2]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[3]"]), main="C", ylim=c(0,90), lwd=2, col="red", xlim=c(0,0.2), xlab="Probability")
lines(density(chains[,"pS[3]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[4]"]), main="D", ylim=c(0,30), lwd=2, col="red", xlim=c(0,0.2), xlab="Probability")
lines(density(chains[,"pS[4]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[5]"]), main="E", ylim=c(0,200), lwd=2, col="red", xlim=c(0,0.2), xlab="Probability")
lines(density(chains[,"pS[5]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[6]"]), main="F", ylim=c(0,50), lwd=2, col="red", xlim=c(0,0.2), xlab="Probability")
lines(density(chains[,"pS[6]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[7]"]), main="G", ylim=c(0,20), lwd=2, col="red", xlim=c(0.5,1), xlab="Probability")
lines(density(chains[,"pS[7]"]), lwd=2, col="blue")  
