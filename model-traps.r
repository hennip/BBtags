
library(runjags)



  # yhdistetty rannikko + jokikalastuskuolevuus, tälle priori elämänkiertomallista?
  # Otetaanko esim. viimeiseltä kolmelta vuodelta, ts 2018-2020?
  
  # Kalastuskuolevuus hetkelliseksi, otetaan kuitenkin huomioon koko kauden ajalta
  # Suurin osa kaloista merkitty kesäkuun lopussa, kalastuskausi kuitenkin käynnistyy 
  # kunnolla vasta samoihin aikoihin

#  s[j]=survive_sea[j]*keep_tag
#  survive_sea[j]=exp(-(M_sea)+h[j]+F_sea)) 
  ## survival based on instantenous mortality rates per 4 weeks
  # M: natural mortality, h: trap induced mortality
  # F: fishing mortality
  

  # pitäisikö tässä kuitenkin ottaa huomioon aika? Eli kuinka
  # monta viikkoa kala selviää hengissä merkinnän/käsitttelyn jälkeen?
  

M2<-"model{
  for(j in 1:2){ # number of gear types 1:Kaukalo; 2:Sukka
  for(t in 1:T){
    #Tagged[t,j]~dbin((1-handling_mort),Tagged_o[t,j])
    Tagged[t,j]<-Tagged_o[t,j]

    C[t,1,j]=harv_sea[j]*N[t,1,j]              # number of tagged fish caught at sea
    R[t,1,j]~dpois(C[t,1,j]*reporting_sea+0.001)    # number of tagged fish reported = data

    C[t,2,j]=harv_river[j]*N[t,2,j] # number of tagged fish caught at river
    R[t,2,j]~dpois(C[t,2,j]*reporting_river+0.001)
  }
  }
  
    
  # First time step: the initial release, all areas except river
  for(j in 1:2){
    N[1,1,j]=Tagged[1,j]    
    N[1,2,j]=0                # no releases in the river (A+1)
  
  # Number of fish at  sea through time 
  # s:surviving that stays at the area
  for(t in 1:(T-1)){
    N[t+1,1,j]=N[t,1,j]*s[t,1,j]+Tagged[t+1,j]
  }
 
  # Number of fish at river through time 
  # m:surviving that moves from sea to river
  # s:surviving that stays at the area
  for(t in 1:T){
    N[t+1,2,j]=N[t,1,j]*m[t,j]+N[t,2,j]*s[t,2,j]
  }
 
  
  for(t in 1:T){
    # surviving fractions, these should be made stochastic in future, using dbeta()
    s[t,1,j]=surv_sea[j]*keep_tag*(1-move) # surviving fraction that stays at sea
    s[t,2,j]=surv_river[j]*keep_tag        # surviving fraction that stays at river
    m[t,j]=surv_sea[j]*keep_tag*move     # surviving fraction that moves from sea to river
  }
  
  Pdie_sea[j]=(h[j]+hand_inst*haav_prop)/(M_sea/12+h[j]+F_sea/11+hand_inst*haav_prop)*
  (1-exp(-(M_sea/12+h[j]+F_sea/11+hand_inst*haav_prop)))

  Pdie_river[j]=(h[j]+hand_inst*haav_prop)/(M_river/52+h[j]+F_river/12+hand_inst*haav_prop)*
  (1-exp(-(M_river/52+h[j]+F_river/12+hand_inst*haav_prop)))
  
    surv_sea[j]=exp(-(M_sea/12+h[j]+F_sea/11+hand_inst))  # survival based on instantenous mortality rates per week
                                                # M: natural mortality, h: trap induced mortality
                                                # F: fishing mortality
                                                
    harv_sea[j]=keep_tag*(1-surv_sea[j])*(F_sea/11)/(M_sea/12+F_sea/11+h[j]+hand_inst) # harvest based on inst. mortalities
    
    surv_river[j]=exp(-(M_river/52+F_river/12+h[j]+hand_inst))
    harv_river[j]=keep_tag*(1-surv_river[j])*(F_river/12)/(M_river/52+F_river/12+h[j]+hand_inst)
  
    h[j]~dunif(0,12)#trap induced mortality
  }

  F_sea~dlnorm(log(0.1566)-0.5/T_sea, T_sea) # sea F based on WGBAST
  T_sea=1/log(pow(0.02062/0.1566,2)+1)

  F_river~dlnorm(log(0.2119)-0.5/T_river, T_river) # river F based on WGBAST
  T_river=1/log(pow(0.07043/0.2119,2)+1)

  M_river~dlnorm(log(0.0975)-0.5/T_Mr, T_Mr) # normal inst M based on WGBAST, 52 weeks
  T_Mr=1/log(pow(0.0312/0.0975,2)+1)
  
  M_sea~dlnorm(log(0.1622)-0.5/T_Mc, T_Mc) #coastal inst M during 2 months in WGBAST model
  T_Mc<-1/log(pow(0.052/0.1622,2)+1)


  move~dunif(0,1)#dbeta(2,2)             # weakly informative prior for movement probability 
  moveP~dunif(0,1)#dbeta(2,2)             # weakly informative prior for movement probability 

  #handling_mort~dbeta(27,132) # from Siira et al 
  handling_mort~dbeta(1.8,7.2) # Ruokonen et al 2021, 12 weeks
  hand_inst=-log(1-handling_mort/12)
  haav_prop~dbeta(4.53, 5.47) # haavinnan osuus käsittelykuolevuudesta Siira et al allaskokeen perusteella
haav_propP~dbeta(4.53, 5.47) # haavinnan osuus käsittelykuolevuudesta Siira et al allaskokeen perusteella

  loose_tag~dbeta(12,68)  # from Siira et al, over three months
  keep_tag_inst=-log(1-loose_tag/12) # inst for one week
  keep_tag=exp(-keep_tag_inst)  # probability for keepin the tag for 1 week

  loose_tagP~dbeta(12,68) 
  handling_mortP~dbeta(1.8,7.2) 
  
  F_seaP~dlnorm(log(0.1566)-0.5/T_seaP, T_seaP) 
  T_seaP=1/log(pow(0.02062/0.1566,2)+1)

  F_riverP~dlnorm(log(0.2119)-0.5/T_riverP, T_riverP) 
  T_riverP=1/log(pow(0.07043/0.2119,2)+1)

  M_riverP~dlnorm(log(0.0975)-0.5/T_MrP, T_MrP) # normal inst M based on WGBAST, 52 weeks
  T_MrP=1/log(pow(0.0312/0.0975,2)+1)
  
  M_seaP~dlnorm(log(0.1622)-0.5/T_McP, T_McP) #coastal inst M during 2 months in WGBAST model
  T_McP<-1/log(pow(0.052/0.1622,2)+1)

reporting_sea=C_X[Y]
reporting_river=R_X[Y]

Y~dcat(p[1:3])

for(i in 1:3){
  p[i]=1/3
}

# Timo
C_X[1]~dbeta(0.69*25,(1-0.69)*25)
R_X[1]~dbeta(0.77*30,(1-0.77)*30)

# Petri
C_X[2]~dbeta(0.71*60, (1-0.71)*60)
R_X[2]~dbeta(0.755*50, (1-0.755)*50)

#Tapsa
C_X[3]~dbeta(0.78*13, (1-0.78)*13)
R_X[3]~dbeta(0.56*13, (1-0.56)*13)


}"


data=list(
  Tagged_o=tagged,
  R=recap,  
  T=15
)  


var_names<- c(
  "haav_prop",  "haav_propP",
  "Pdie_sea", "Pdie_river",
  #"Y", "p",
  "reporting_sea", "reporting_river",
  "reporting_seaP", "reporting_riverP",
  "h", "F_sea", "F_river", 
              "M_sea",             "M_river",
  "M_seaP",             "M_riverP",
  "harv_sea", "harv_river", "move","moveP",
              "handling_mort", "handling_mortP","loose_tag","loose_tagP",
              "F_seaP", "F_riverP")

run0 <- run.jags(M2,
                 monitor= var_names,data=data, #inits = inits,
                 n.chains = 2, method = 'parallel', thin=10, burnin = 1000,
                 modules = "mix",keep.jags.files=F,sample =10000, adapt = 1000,
                 progress.bar=TRUE)


run<-run0

round(summary(run),3)
plot(run)

chains<-as.mcmc(run)

# Rysäkohtainen kuolevuus
#===========================
summary(chains[,"h[1]"], quantiles=c(0.05,0.5,0.95))
summary(chains[,"h[2]"], quantiles=c(0.05,0.5,0.95))

# Absoluuttiseksi kääntäminen mielekästä vain jos ajatellaan että muita kuolevuuksia ei olisi
summary(1-exp(-chains[,"h[1]"]), quantiles=c(0.05,0.5,0.95))
summary(1-exp(-chains[,"h[2]"]), quantiles=c(0.05,0.5,0.95))

summary((1-exp(-chains[,"h[1]"]))/(1-exp(-chains[,"h[2]"])), quantiles=c(0.05,0.5,0.95))





# Osuus kokonaiskuolevuudesta
##################################

summary(chains[,"Pdie_sea[1]"], quantiles=c(0.05,0.5,0.95))
summary(chains[,"Pdie_sea[2]"], quantiles=c(0.05,0.5,0.95))

summary(chains[,"Pdie_sea[1]"]/chains[,"Pdie_sea[2]"], quantiles=c(0.12,0.25,0.5,0.75,0.95))
#summary(chains[,"Pdie_river[1]"]/chains[,"Pdie_river[2]"], quantiles=c(0.11,0.25,0.5,0.75,0.95))


par(mfrow=c(2,2))
plot(density(chains[,"reporting_seaP"]), main="reporting at coastal fishery")
lines(density(chains[,"reporting_sea"]), lwd=2)

plot(density(chains[,"reporting_riverP"]), main="reporting at river fishery")
lines(density(chains[,"reporting_river"]), lwd=2)

plot(density(chains[,"F_seaP"]), main="F at coastal fishery")
lines(density(chains[,"F_sea"]), lty=2)

plot(density(chains[,"F_riverP"]), main="F at river fishery")
lines(density(chains[,"F_river"]), lty=2)


# Instantaneous mortalities
par(mfrow=c(3,3))
plot(density(chains[,"F_seaP"]), main="Kalastuskuol., rannikko", xlab="Hetkellinen kuolleisuus", ylab="Tiheys", xlim=c(0,0.5))
lines(density(chains[,"F_sea"]),lwd=2)
plot(density(chains[,"F_riverP"]), main="Kalastuskuol., jokialue", xlab="Hetkellinen kuolleisuus", ylab="Tiheys", xlim=c(0,0.6))
lines(density(chains[,"F_river"]),lwd=2)
plot(density(chains[,"M_seaP"]), main="Luonnollinen kuol., rannikko", xlab="Hetkellinen kuolleisuus", ylab="Tiheys", xlim=c(0,0.5), ylim=c(0,10))
lines(density(chains[,"M_sea"]),lwd=2)
plot(density(chains[,"M_riverP"]), main="Luonnollinen kuol., jokialue", xlab="Hetkellinen kuolleisuus", ylab="Tiheys", xlim=c(0,0.5))
lines(density(chains[,"M_river"]),lwd=2)

#par(mfrow=c(2,2))
plot(density(chains[,"handling_mortP"]), main="Käsittelykuolleisuus", xlab="Hetkellinen kuolleisuus", ylab="Tiheys", xlim=c(0,1), ylim=c(0,4))
lines(density(chains[,"handling_mort"]),lwd=2)
plot(density(chains[,"haav_propP"]), main="Haavinnan osuus käsittelykuolleisuudesta", xlab="Osuus", ylab="Tiheys", xlim=c(0,1))
lines(density(chains[,"haav_prop"]),lwd=2)
plot(density(chains[,"loose_tagP"]), main="Merkin irtoamisen tn", xlab="Todennäköisyys", ylab="Tiheys", xlim=c(0,0.4), ylim=c(0,12))
lines(density(chains[,"loose_tag"]),lwd=2)

plot(density(chainsP[,"C_A"][[1]]), main="Raportointiaktiivisuus, rannikko", xlab="Todennäköisyys", ylim=c(0,6))
lines(density(chains[,"reporting_sea"]), lwd=2)

plot(density(chainsP[,"R_A"][[1]]), main="Raportointiaktiivisuus, jokialueet", xlab="Todennäköisyys", ylim=c(0,6))
lines(density(chains[,"reporting_river"]), lwd=2)



# Raportoinnin priorit
prior<-"
model{

C_A<-C_X[Y]
R_A<-R_X[Y]

Y~dcat(p[1:3])

for(i in 1:3){
  p[i]<-1/3
}

# Timo
C_X[1]~dbeta(0.69*25,(1-0.69)*25)
R_X[1]~dbeta(0.77*30,(1-0.77)*30)

# Petri
C_X[2]~dbeta(0.71*60, (1-0.71)*60)
R_X[2]~dbeta(0.755*50, (1-0.755)*50)

#Tapsa
C_X[3]~dbeta(0.78*13, (1-0.78)*13)
R_X[3]~dbeta(0.56*13, (1-0.56)*13)

}"


cat(prior,file="priorTot.txt")
system.time(jm<-jags.model('priorTot.txt',n.adapt=100,n.chains=1))

system.time(chainsP<-coda.samples(jm,variable.names=c(
  "C_A","R_A","C_X","R_X", "p", "Y"
),n.iter=1000000,thin=1))

#summary(chainsP, quantiles=c(0.025,0.5,0.975))

par(mfrow=c(1,2))
plot(density(chainsP[,"C_A"][[1]]), ylim=c(0,7), main="Rannikkokalastus")
lines(density(chainsP[,"C_X[1]"][[1]]), lty=2, col="red")
lines(density(chainsP[,"C_X[2]"][[1]]), lty=3, col="green")
lines(density(chainsP[,"C_X[3]"][[1]]), lty=4, col=4)
lines(density(chains[,"reporting_sea"]), lty=1, lwd=2)


plot(density(chainsP[,"R_A"][[1]]), ylim=c(0,8), main="Jokikalastus")
lines(density(chainsP[,"R_X[1]"][[1]]), lty=2, col="red")
lines(density(chainsP[,"R_X[2]"][[1]]), lty=3, col="green")
lines(density(chainsP[,"R_X[3]"][[1]]), lty=4, col=4)
lines(density(chains[,"reporting_river"]), lty=1, lwd=2)

