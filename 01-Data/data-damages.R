library(runjags)
require(rjags)
require(sp)
require(openxlsx)
library(readxl)
library(tidyverse)
library(lubridate)

# gear n_tagged n_return      p
# <dbl>    <int>    <dbl>  <dbl>
# 1     1      201       20 0.0995 #Kaukalo
# 2     2      211       31 0.147 #Sukka

dat2=read_xlsx("C:/Users/03080932/OneDrive - Valtion/Projects/BBTags/dat/orig/Vauriot_Data_2021_ver2.xlsx")%>%
  filter(pituus>60) #Jätetään pois 1 mahd. 1SW kala kuten JAGS-mallissa

t1<-dat2%>%group_by(rysätyyppi,vaurio)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=vaurio, type="vaurio")

t2<-dat2%>%group_by(rysätyyppi,vaurio_ei_suomu)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=vaurio_ei_suomu, type="vaurio_ei_suomu")

t4<-dat2%>%group_by(rysätyyppi,Silmävaurio)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=Silmävaurio, type="silmävaurio")

t5<-dat2%>%group_by(rysätyyppi,`Tuoreita purema/viiltohaavoja`)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=`Tuoreita purema/viiltohaavoja`, type="Tuor_pur_vi")

t6<-dat2%>%group_by(rysätyyppi,`Tuoreita pantamaisia haavoja (havaspyydys)`)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=`Tuoreita pantamaisia haavoja (havaspyydys)`, type="Tuor_panta")

t7<-dat2%>%group_by(rysätyyppi,`Kroonisia purema/viiltohaavoja`)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=`Kroonisia purema/viiltohaavoja`, type="Kroo_pur_vi")

t8<-dat2%>%group_by(rysätyyppi,`Kroonisia pantamaisia haavoja (havaspyydys)`)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=`Kroonisia pantamaisia haavoja (havaspyydys)`, type="Kroo_panta")

t9<-dat2%>%group_by(rysätyyppi,`Muita tuoreita haavoja ihon läpi`)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=`Muita tuoreita haavoja ihon läpi`, type="Muu_tuore")

t10<-dat2%>%group_by(rysätyyppi,suuvaurio)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=suuvaurio, type="suuvaurio")

t11<-dat2%>%group_by(rysätyyppi,Evävaurio)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=Evävaurio, type="evavaurio")


# Suomuvauriot
t30<-dat2%>%
  mutate(suomuvaurio12=ifelse(suomuvaurio==2, 1, suomuvaurio))%>%
  group_by(rysätyyppi,suomuvaurio12)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=suomuvaurio12, type="suomuvaurio_1tai2")

t31<-dat2%>%group_by(rysätyyppi,suomuvaurio)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=ifelse(suomuvaurio==1,1,0), type="suomuvaurio1")

t32<-dat2%>%group_by(rysätyyppi,suomuvaurio)%>%
  summarise(n=n())%>%
  mutate(Ntot=if_else(rysätyyppi=="Kaukalo", 201, 211))%>%
  mutate(p=n/Ntot, x=ifelse(suomuvaurio==2,1,0), type="suomuvaurio2")


dfd<-full_join(t1,t2)%>%
  full_join(t4)%>%
  full_join(t5)%>%
  full_join(t6)%>%
  full_join(t7)%>%
  full_join(t8)%>%
  full_join(t9)%>%
  full_join(t10)%>%
full_join(t11)%>%
  full_join(t30)%>%
  full_join(t31)%>%
  full_join(t32)%>%
  filter(x==1)%>%
  select(rysätyyppi, n, Ntot, p,x, type)%>%
  arrange(x,rysätyyppi)
View(dfd)


dfd%>%filter(rysätyyppi=="Kaukalo")%>%select(n, type)
dfd%>%filter(rysätyyppi=="Sukka")%>%select(n, type)


Mdamage<-"
model{  
for(i in 1:13){
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

data=list(xK=c(184,151,4,6,1,7,2,2,11,141,169,134,35),
          xS=c(207,188,1,5,0,4,0,1,3,186,195,108,87),
          NK=201,NS=211)

 run0 <- run.jags(Mdamage,
                  monitor= c("prob", "pK", "pS"),data=data, #inits = inits,
                  n.chains = 2, method = 'parallel', thin=10, burnin = 1000,
                  modules = "mix",keep.jags.files=F,sample =10000, adapt = 1000,
                  progress.bar=TRUE)
 
 run<-run0
 
 round(summary(run),3)

chains<-as.mcmc(run) 

par(mfrow=c(3,3))
plot(density(chains[,"pK[1]"]), main="vaurio", lwd=2, col="red", xlim=c(0,1), ylim=c(0,40))
lines(density(chains[,"pS[1]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[2]"]), main="vaurio ei suomu", ylim=c(0,20), lwd=2, col="red", xlim=c(0,1))
lines(density(chains[,"pS[2]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[3]"]), main="silmävaurio", ylim=c(0,100), lwd=2, col="red", xlim=c(0,1))
lines(density(chains[,"pS[3]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[4]"]), main="Tuoreita purema/viiltohaavoja", ylim=c(0,50), lwd=2, col="red", xlim=c(0,1))
lines(density(chains[,"pS[4]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[5]"]), main="Tuoreita pantamaisia haavoja", ylim=c(0,200), lwd=2, col="red", xlim=c(0,1))
lines(density(chains[,"pS[5]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[6]"]), main="Kroonisia purema/viiltohaavoja", ylim=c(0,50), lwd=2, col="red", xlim=c(0,1))
lines(density(chains[,"pS[6]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[7]"]), main="Kroonisia pantamaisia haavoja", ylim=c(0,200), lwd=2, col="red", xlim=c(0,1))
lines(density(chains[,"pS[7]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[8]"]), main="Muita tuoreita haavoja ihon läpi", ylim=c(0,80), lwd=2, col="red", xlim=c(0,1))
lines(density(chains[,"pS[8]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[9]"]), main="suuvaurio", ylim=c(0,60), lwd=2, col="red", xlim=c(0,1))
lines(density(chains[,"pS[9]"]), lwd=2, col="blue")  

par(mfrow=c(2,2))
plot(density(chains[,"pK[10]"]), main="evävaurio", ylim=c(0,20), lwd=2, col="red", xlim=c(0,1))
lines(density(chains[,"pS[10]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[11]"]), main="suomuvaurio vähäinen tai huomattava", ylim=c(0,30), lwd=2, col="red", xlim=c(0,1))
lines(density(chains[,"pS[11]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[12]"]), main="suomuvaurio vähäinen", ylim=c(0,20), lwd=2, col="red", xlim=c(0,1))
lines(density(chains[,"pS[12]"]), lwd=2, col="blue")  
plot(density(chains[,"pK[13]"]), main="suomuvaurio huomattava", ylim=c(0,20), lwd=2, col="red", xlim=c(0,1))
lines(density(chains[,"pS[13]"]), lwd=2, col="blue")  

      suomuvaurio (vähäinen tai huomattava
                   suomuvaurio vähäinen
                   suomuvaurio huomattava
                   
