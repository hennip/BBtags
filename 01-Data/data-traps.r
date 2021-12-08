require(rjags)
require(sp)
require(openxlsx)
library(readxl)
library(tidyverse)
library(lubridate)

#### Data preparation


dat=read_xlsx("H:/Projects/BBTags/dat/der/Lohirysä_data_2021_der.xlsx")
#              options = "ENCODING=WINDOWS-1252")

dat%>%filter(pituus<70)%>% select(pituus, everything())
dim(dat)

df<-dat%>%filter(pituus>60)%>% # Take only MSW salmon
  select(pvm, fisher, rysätyyppi, palautettu, saantipvm, joki, meri)%>%
  mutate(gear=ifelse(rysätyyppi=="Kaukalo", 1, ifelse(rysätyyppi=="Sukka", 2, NA)))%>%
  #mutate(recap=ifelse(palautettu==1,1,ifelse(is.na(palautettu)==T, 0, NA)))
mutate(recap=ifelse(is.na(palautettu)==T, 0, 1))%>%
  mutate(lag=saantipvm-pvm)%>%
  mutate(week1=week(pvm)-20)%>%
  mutate(week2=week(saantipvm)-20)%>%
  mutate(recap_area= case_when(meri==1 & is.na(joki)==T ~1,
                         joki==1 & is.na(meri)==T~2))%>%
  select(week1, week2, pvm, recap_area, saantipvm, everything())

dim(df)
head(df)
#View(df)


df%>%group_by(gear)%>%
  summarise(n_tagged=n(),
            n_return=sum(palautettu, na.rm = T),
            p=n_return/n_tagged)

# Simple input
(input_dat<-df%>%group_by(gear, fisher)%>%
    summarise(n_tagged=n(),
              n_return=sum(palautettu, na.rm = T),
              p=n_return/n_tagged))


#View(filter(df, palautettu==1))

# 3 palautettua kalaa joille tieto saantipäivästä puuttuu:
# Viikolla 4 merkitty sukka-kala
# Viikolla 5 merkitty kaukalo- ja sukkakala

filter(df, week1==4, rysätyyppi=="Sukka", palautettu==1, is.na(week2)==F)%>%
  summarise(w=round(mean(week2))) #8
filter(df, week1==5, rysätyyppi=="Sukka", palautettu==1, is.na(week2)==F)%>%
  summarise(w=round(mean(week2))) #8
filter(df, week1==5, rysätyyppi=="Kaukalo", palautettu==1, is.na(week2)==F)%>%
  summarise(w=round(mean(week2))) #9

# imputoidaan keskimääräiset palautusviikot puuttuviin

df<-df%>%
  mutate(week2=ifelse((week1==4 & rysätyyppi=="Sukka" & palautettu==1 & is.na(week2)==T), 8, week2 ))%>%
  mutate(week2=ifelse((week1==5 & rysätyyppi=="Sukka" & palautettu==1 & is.na(week2)==T), 8, week2 ))%>%
  mutate(week2=ifelse((week1==5 & rysätyyppi=="Kaukalo" & palautettu==1 & is.na(week2)==T), 9, week2 ))

#View(filter(df, palautettu==1))

# Tsekataan ne jotka palautettu, mutta palautuspaikka tuntematon
filter(df, palautettu==1, is.na(joki)==T & is.na(meri)==T) # 1kpl
filter(df, week1==4, rysätyyppi=="Sukka", palautettu==1, is.na(week2)==F)%>%
  select(joki, everything()) # joki 7 kpl, meri 8 kpl

# imputoidaan puuttuva saantitieto mereltä tulleeksi
df<-df%>%
  mutate(recap_area=ifelse((rysätyyppi=="Sukka" & palautettu==1 & is.na(meri)==T &is.na(joki)==T), 
                     1, recap_area ))
  

#View(filter(df, palautettu==1, rysätyyppi=="Sukka"))
#View(filter(df, palautettu==1, rysätyyppi=="Kaukalo"))
max(df$week1) #6
max(df$week2, na.rm=T) #15


dfT_S<-filter(df, rysätyyppi=="Sukka")
dfT_K<-filter(df, rysätyyppi=="Kaukalo")

dfRsea_S<-filter(df, palautettu==1, recap_area==1, rysätyyppi=="Sukka")
dfRsea_K<-filter(df, palautettu==1, recap_area==1, rysätyyppi=="Kaukalo")
dfRriv_S<-filter(df, palautettu==1, recap_area==2, rysätyyppi=="Sukka")
dfRriv_K<-filter(df, palautettu==1, recap_area==2, rysätyyppi=="Kaukalo")



# Tagged and recaptured per week 1:15
tagged<-array(NA, dim=c(15,2)) # week,  gear 
for(j in 1:15){
  tmpS<-tmpK<-0
  for(i in 1:dim(dfT_S)[1]){
    if(dfT_S$week1[i]==j){tmpS<-tmpS+1}
  }
  tagged[j,2]<-tmpS

  for(i in 1:dim(dfT_K)[1]){
    if(dfT_K$week1[i]==j){tmpK<-tmpK+1}
  }
  tagged[j,1]<-tmpK
}
tagged

recap<-array(NA, dim=c(15,2,2)) # week, recap_area, gear
for(j in 1:15){
  tmp1S<-tmp1K<-0 # recap sea
  tmp2S<-tmp2K<-0 # recap river
  
  # sea
  for(i in 1:dim(dfRsea_S)[1]){
    if(dfRsea_S$week2[i]==j){tmp1S<-tmp1S+1}
  }
  recap[j,1,2]<-tmp1S
  
  for(i in 1:dim(dfRsea_K)[1]){
    if(dfRsea_K$week2[i]==j){tmp1K<-tmp1K+1}
  }
  recap[j,1,1]<-tmp1K
  
  # river
  for(i in 1:dim(dfRriv_S)[1]){
    if(dfRriv_S$week2[i]==j){tmp2S<-tmp2S+1}
  }
  recap[j,2,2]<-tmp2S
  
  for(i in 1:dim(dfRriv_K)[1]){
    if(dfRriv_K$week2[i]==j){tmp2K<-tmp2K+1}
  }
  recap[j,2,1]<-tmp2K
  
}
recap

sum(tagged)
sum(recap)







# Huomioita: 
# Kalastaja 1:llä vähiten merkittyjä kaloja, merkinnät keskittyvät kauden alkuun
# ja merkkipalautuksia on vain 1/(21+20) (sukka)
# kalastajilla 1 ja 2 hyvin vähän eroa takaisinsaannissa 
# pyydystyyppien välillä. 
#Kalastajilla 3 ja 4 erot selkeitä, tosin
# kalastaja 4:llä merkinnät keskittyneet voimakkaasti tietylle päivälle


df%>%group_by(fisher)%>%
  summarise(n_tagged=n(),
            n_return=sum(palautettu, na.rm = T),
            p=n_return/n_tagged)


(sum_tagged<-df%>%group_by(gear)%>%
  summarise(n_tagged=n()))

print("Kaukalo")
df%>%filter(gear==1)%>%group_by(pvm)%>%
            summarise(n_tagged=n(),
                      osuus_pvm=round(n_tagged/sum_tagged[1,2],2),
                      n_return=sum(palautettu, na.rm = T),
                      p=round(n_return/n_tagged,2))
print("Sukka")
df%>%filter(gear==2)%>%group_by(pvm)%>%
  summarise(n_tagged=n(),
            osuus_pvm=round(n_tagged/sum_tagged[2,2],2),
            n_return=sum(palautettu, na.rm = T),
            p=round(n_return/n_tagged,2))


for(k in 1:4){#Kalastaja
  #k<-4
  df2<-filter(df, fisher==k); print (paste0("Fisher ",k))
  sum_tagged<-df2%>%group_by(gear)%>%
    summarise(n_tagged=n())
  
    print(
    df2%>%filter(gear==1)%>%group_by(pvm)%>%
    summarise(n_tagged=n(),
              osuus_pvm=round(n_tagged/sum_tagged[1,2],2),
              n_return=sum(palautettu, na.rm = T),
              p=round(n_return/n_tagged,2)))
    
  print(
    df2%>%filter(gear==2)%>%group_by(pvm)%>%
    summarise(n_tagged=n(),
              osuus_pvm=round(n_tagged/sum_tagged[2,2],2),
              n_return=sum(palautettu, na.rm = T),
              p=round(n_return/n_tagged,2))
  )
  }
# Huomioita: Kalastaja 4 on merkinnyt 22.6. 65% kaukalomerkeistään ja 16.6. 53% sukkamerkeistään
# N'iden isojen merkintäerien takaisinsaanti näyttäisi olevan ehkä hieman alakanttiin 
# verrattuna omiin ja kollegoiden merkintöihin. Katsotaan miten tilanne muuttuu kun saadaan lisää dataa.




df%>%
  summarise(min_lag=min(lag, na.rm=T), 
            max_lag=max(lag, na.rm=T), 
            mean_lag=mean(lag, na.rm=T),
            sd_lag=sd(lag, na.rm=T))


df%>%group_by(fisher,pvm)%>%summarise(n=n())


par(mfrow=c(1,2))
tmp<-filter(df, gear==1)
plot(jitter(tmp$fisher,2)~tmp$pvm, xlab="date", ylab="fisher", main="Kaukalo", col=tmp$fisher)
tmp<-filter(df, gear==2)
plot(jitter(tmp$fisher,2)~tmp$pvm, xlab="date", ylab="fisher", main="Sukka", col=tmp$fisher)


# 
# 
# 
# 
# # Pooling observations to weekly counts of releases and recaptures
# # Stratfying by defined areas: Selk?meri, Merenkurkku, Pohjanlahti, Joki
# 
# dat=read.xlsx("H:/Projects/BBTags/mark_recapture.xlsx")
# coord=data.frame(x=as.numeric(dat$`E.ETRS-TM35FIN`),y=as.numeric(dat$`N.ETRS-TM35FIN`),mark=dat$paikka
#                  ,recap=dat$palautuspaikka,markday=dat$merkintapvm,recapday=dat$saantipvm)
# coord$mark=ifelse(coord$mark=="Merikarvia" | coord$mark=="Pori",1,
#                   ifelse(coord$mark=="Luoto",2,
#                          ifelse(coord$mark=="Haukipudas" | coord$mark=="Ajos",3,4)))
# 
# coord$recap=ifelse(coord$y<7000000,1,
#                    ifelse(coord$y>700000 & coord$y<7200000,2,
#                           ifelse(coord$y>7200000 & coord$recap=="meri",3,4)))
# plot(coord$x,coord$y,col=coord$recap)
# 
# coord$markday=coord$markday-min(coord$markday)+1
# # Miksi recapday alkaa 1:stä eikä mene samassa suhteessa markdayn kanssa?
# coord$recapday=coord$recapday-min(coord$recapday,na.rm=TRUE)+1
# # > min(coord$markday)
# #[1] 43975
# #> min(coord$recapday,na.rm=TRUE)
# #[1] 43982
# 43982-43975
# 
# 
# # Mitä tämä tismalleen ottaen tekee?
# # Määrittää viikon kullekin merkinnälle ja palautukselle?
# # Perustuu markday:hin ja recapday:hin, siis markweek ja recapweek eivät ole yhteydessä toisiinsa?
# for(i in 1:20){
#   beg=(i-1)*7+1
#   end=(i-1)*7+7
#   for(k in 1:length(coord$x)){
#   if(coord$markday[k]<=end && coord$markday[k] >=beg){ coord$markweek[k]=i}
#     if(!is.na(coord$recapday[k]) && coord$recapday[k]<=end && coord$recapday[k] >=beg){ coord$recapweek[k]=i}
#   if(is.na(coord$recapday[k])) coord$recapweek[k]=NA
#     }
# }
# 
# 
# 
# Tagged=array(0,dim=c(20,4,3))  #nrow=20,ncol=4)
# Recap=array(0,dim=c(20,4,3)) #nrow=20,ncol=4)
# r=coord[!is.na(coord$x),]
# 
# # Mikä on a:n ja t:n ero? Molemmat perustuvat merkintäpaikkaan (mark)
# for(t in 1:3){
# for(w in 1:20){ # viikko
#   for(a in 1:4){
#     Tagged[w,a,t]=length(coord$x[coord$markweek==w & coord$mark==a & coord$mark==t])
#     Recap[w,a,t]=length(r$x[r$recapweek==w & r$recap==a & r$mark==t])
#   }
# }
# }
# 
# data=list(T=19,A=3,K=3,Tagged_o=Tagged,R=Recap,
#           M_sea=0.1,           # Assuming known natural mortalities for now, should use a prior later
#           M_river=0.2,
#          
#           reporting_river=0.8, # Assuming known reporting rates for now, should use a prior later
#           reporting_sea=0.5)   
# 
# jm=jags.model("jags_model.r",data,n.chains=4)  # compiling the jags model
# 
# # Specifying variables to be monitored
# 
# monitor=c("h","F_sea","F_river","move","harvest_river","harvest_sea","keep_tag","survive_sea","survive_river")
# chains=coda.samples(jm,monitor,n.iter=1000) # Running the MCMC simulation
# 
# # Dumping all MCMC plots into a PDF file
# 
# pdf(file="results.pdf")
# plot(chains)
# dev.off()
# 
# # Comparing the trap-related mortality components
# 
# windows(14,7)
# 
# d=as.matrix(chains)
# 
# par(mfrow=c(1,2))
# 
# plot(density(d[,"h[1]"]),lwd=4, type="l",xlab="Inst. mortality rate/week",ylab="Probability density",main="",xlim=c(0,1.3))
# points(density(d[,"h[2]"]),lwd=4,type="l",col="red")
# points(density(d[,"h[3]"]),lwd=4,type="l",col="blue")
# 
# p1=sum(d[,"h[1]"]>d[,"h[2]"])/4000
# p2=sum(d[,"h[1]"]>d[,"h[3]"])/4000
# p3=sum(d[,"h[2]"]>d[,"h[3]"])/4000
# 
# text(1,6,paste("P(P-U > Paunet)=",p1))
# text(1,5.5,paste("P(P-U > P-U sock)=",p2))
# text(1,5,paste("P(Paunet > P-U sock)=",p3))
# 
# legend("topright",legend=c("Selk?meri : P-U ","Merenkurkku : Paunet","Per?meri : P-U sock"),col=c("black","red","blue"),lwd=c(4,4,4))
# 
# plot(density(exp(-d[,"h[1]"])),lwd=4, type="l",xlab="Survival rate/week",ylab="Probability density",main="",xlim=c(0,1))
# points(density(exp(-d[,"h[2]"])),lwd=4,type="l",col="red")
# points(density(exp(-d[,"h[3]"])),lwd=4,type="l",col="blue")
# 
# legend("topleft",legend=c("Selk?meri : P-U ","Merenkurkku : Paunet","Per?meri : P-U sock"),col=c("black","red","blue"),lwd=c(4,4,4))
# 
# text(0.3,6,paste("P(P-U > Paunet)=",1-p1))
# text(0.3,5.5,paste("P(P-U > P-U sock)=",1-p2))
# text(0.3,5,paste("P(Paunet > P-U sock)=",1-p3))
