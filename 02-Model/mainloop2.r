#setwd("c:/models/bbtags")
require(rjags)
require(sp)
require(openxlsx)
library(readxl)
library(tidyverse)
library(lubridate)

#### Data preparation

# Pooling observations to weekly counts of releases and recaptures
# Stratfying by defined areas: Selk?meri, Merenkurkku, Pohjanlahti, Joki

dat=read_xlsx("H:/Projects/BBTags/mark_recapture.xlsx", na="")
#View(dat[1:3,])
coord<-dat%>%mutate(x=`E ETRS-TM35FIN`, y=`N ETRS-TM35FIN`,
                    mark2=`paikka...2`,
                    recap=palautuspaikka,markday2=merkintapvm,
                    recapday2=saantipvm)%>%
  mutate(mark=ifelse(mark2=="Merikarvia" | mark2=="Pori",1,
                      ifelse(mark2=="Luoto",2,
                             ifelse(mark2=="Haukipudas" | mark2=="Ajos",3,4))))%>%
  mutate(recap=ifelse(y<7000000,1,
                            ifelse(y>700000 & y<7200000,2,
                                   ifelse(y>7200000 & recap=="meri",3,4)))
  )%>%
  select(x,y,mark,markday2, recapday2)%>%
#mutate(markday=markday2-min(markday2)+1)%>%
  #mutate(recapday=recapday2-min(recapday2,na.rm=T)+1)%>%
  #select(x,y,mark,markday, markday2, recapday, recapday2)
  #coord=data.frame(x=as.numeric(dat$`E.ETRS-TM35FIN`),y=as.numeric(dat$`N.ETRS-TM35FIN`),
  #                mark=dat$paikka,
  #                recap=dat$palautuspaikka,markday=dat$merkint?pvm,recapday=dat$saantipvm)
  # coord$mark=ifelse(coord$mark=="Merikarvia" | coord$mark=="Pori",1,
  #                   ifelse(coord$mark=="Luoto",2,
  #                          ifelse(coord$mark=="Haukipudas" | coord$mark=="Ajos",3,4)))
  # 
  # coord$recap=ifelse(coord$y<7000000,1,
#                    ifelse(coord$y>700000 & coord$y<7200000,2,
#                           ifelse(coord$y>7200000 & coord$recap=="meri",3,4)))
mutate(mark_yday=yday(markday2))%>%
  mutate(recap_yday=yday(recapday2))
  
min(coord$mark_yday)
min(coord$recap_yday, na.rm=T)
152-145


# recapday:n määrittely vähän jännä, miksi tämäkin alkaa 1:stä?
# Toimisiko ydayn avulla määrittely paremmin?
tmp<-coord%>%mutate(markday=mark_yday-min(mark_yday)+1)%>%
  mutate(recapday=recap_yday-min(recap_yday, na.rm = T)+1)
View(tmp)                    
View(coord)
plot(coord$x,coord$y,col=coord$recap)

#coord$markday=coord$markday-min(coord$markday)+1
#coord$recapday=coord$recapday-min(coord$recapday,na.rm=TRUE)+1

markweek<-c()
recapweek<-c()
for(i in 1:20){
  beg=(i-1)*7+1
  end=(i-1)*7+7
  for(k in 1:length(coord$x)){
  if(coord$markday[k]<=end && coord$markday[k] >=beg){ markweek[k]=i}
    if(!is.na(coord$recapday[k]) && coord$recapday[k]<=end && coord$recapday[k] >=beg){ recapweek[k]=i}
  if(is.na(coord$recapday[k])) recapweek[k]=NA
    }
}
#coord<-
  coord%>%mutate(recapweek=recapweek)
View(coord)

Tagged=array(0,dim=c(20,4,3))  #nrow=20,ncol=4)
Recap=array(0,dim=c(20,4,3)) #nrow=20,ncol=4)
r=coord[!is.na(coord$x),]

for(t in 1:3){
for(w in 1:20){
  for(a in 1:4){
    Tagged[w,a,t]=length(coord$x[coord$markweek==w & coord$mark==a & coord$mark==t])
    Recap[w,a,t]=length(r$x[r$recapweek==w & r$recap==a & r$mark==t])
  }
}
}

data=list(T=19,A=3,K=3,Tagged_o=Tagged,R=Recap,
          M_sea=0.1,           # Assuming known natural mortalities for now, should use a prior later
          M_river=0.2,
         
          reporting_river=0.8, # Assuming known reporting rates for now, should use a prior later
          reporting_sea=0.5)   

jm=jags.model("jags_model.r",data,n.chains=4)  # compiling the jags model

# Specifying variables to be monitored

monitor=c("h","F_sea","F_river","move","harvest_river","harvest_sea","keep_tag","survive_sea","survive_river")
chains=coda.samples(jm,monitor,n.iter=1000) # Running the MCMC simulation

# Dumping all MCMC plots into a PDF file

pdf(file="results.pdf")
plot(chains)
dev.off()

# Comparing the trap-related mortality components

windows(14,7)

d=as.matrix(chains)

par(mfrow=c(1,2))

plot(density(d[,"h[1]"]),lwd=4, type="l",xlab="Inst. mortality rate/week",ylab="Probability density",main="",xlim=c(0,1.3))
points(density(d[,"h[2]"]),lwd=4,type="l",col="red")
points(density(d[,"h[3]"]),lwd=4,type="l",col="blue")

p1=sum(d[,"h[1]"]>d[,"h[2]"])/4000
p2=sum(d[,"h[1]"]>d[,"h[3]"])/4000
p3=sum(d[,"h[2]"]>d[,"h[3]"])/4000

text(1,6,paste("P(P-U > Paunet)=",p1))
text(1,5.5,paste("P(P-U > P-U sock)=",p2))
text(1,5,paste("P(Paunet > P-U sock)=",p3))

legend("topright",legend=c("Selk?meri : P-U ","Merenkurkku : Paunet","Per?meri : P-U sock"),col=c("black","red","blue"),lwd=c(4,4,4))

plot(density(exp(-d[,"h[1]"])),lwd=4, type="l",xlab="Survival rate/week",ylab="Probability density",main="",xlim=c(0,1))
points(density(exp(-d[,"h[2]"])),lwd=4,type="l",col="red")
points(density(exp(-d[,"h[3]"])),lwd=4,type="l",col="blue")

legend("topleft",legend=c("Selk?meri : P-U ","Merenkurkku : Paunet","Per?meri : P-U sock"),col=c("black","red","blue"),lwd=c(4,4,4))

text(0.3,6,paste("P(P-U > Paunet)=",1-p1))
text(0.3,5.5,paste("P(P-U > P-U sock)=",1-p2))
text(0.3,5,paste("P(Paunet > P-U sock)=",1-p3))
