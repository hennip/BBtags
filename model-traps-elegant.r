
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
  for(t in 1:T){
    for(i in 1:2){ # Area (1: sea, 2: river)  
      for(j in 1:2){ # Gear (1: pontoon trap with an emptying chute, 2: PU trap with a lifting bag)

        C[t,i,j]=HR[i,j]*N[t,i,j]              
        R[t,i,j]~dpois(C[t,i,j]*p_rep[i]+0.001)

      }
    }
  }  
    
  # First time step: the initial release
  for(j in 1:2){
    N[1,1,j]=Tagged[1,j]    
    N[1,2,j]=0        # no releases in the river
  
    # Number of fish at sea through time 
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
      s[t,1,j]=surv[1,j]*k*(1-move) # surviving fraction that stays at sea
      s[t,2,j]=surv[2,j]*k        # surviving fraction that stays at river
      m[t,j]=surv[1,j]*k*move     # surviving fraction that moves from sea to river
    }
  
    for(i in 1:2){
      Pdie[i,j]=(h[j]+inst_hand*q_net)/(M[i]+h[j]+F[i]+inst_hand*q_net)*
                (1-exp(-(M[i]+h[j]+F[i]+inst_hand*q_net)))
      surv[i,j]=exp(-(M[i]+h[j]+F[i]+inst_hand))  # survival based on instantenous mortality rates per week
      HR[i,j]=k*(1-surv[i,j])*F[i]/(M[i]+F[i]+h[j]+inst_hand) # harvest based on inst. mortalities
    }
  }
  
  # Uninformative priors
  move~dunif(0,1)
  for(j in 1:2){
    h[j]~dunif(0,12)#trap induced mortality
  }

  # Informative priors based on Baltic salmon stock assessment model
  F_sea~dlnorm(log(0.1566)-0.5/T_sea, T_sea) # sea F, 11 weeks 
  T_sea=1/log(pow(0.02062/0.1566,2)+1)

  F_river~dlnorm(log(0.2119)-0.5/T_river, T_river) # river F, 12 weeks 
  T_river=1/log(pow(0.07043/0.2119,2)+1)

  M_river~dlnorm(log(0.0975)-0.5/T_Mr, T_Mr) # river M, 52 weeks
  T_Mr=1/log(pow(0.0312/0.0975,2)+1)
  
  M_sea~dlnorm(log(0.1622)-0.5/T_Mc, T_Mc) # sea  M, 8 weeks
  T_Mc<-1/log(pow(0.052/0.1622,2)+1)

  # Instantaneous mortalities per week
  M[1]=M_sea/8
  M[2]=M_river/52
  F[1]=F_sea/11
  F[2]=F_river/12

  # Other informative priors

  handling_mort~dbeta(1.8,7.2) # Ruokonen et al 2021, 12 weeks
  inst_hand=-log(1-handling_mort/12)

  q_net~dlnorm(-0.723,7.2485)T(,1)# haavinnan osuus käsittelykuolevuudesta Siira et al allaskokeen perusteella
  q_netP~dlnorm(-0.723,7.2485)T(,1)

  loose_tag~dbeta(12,68)  # from Siira et al, over three months
  keep_tag_inst=-log(1-loose_tag/12) # inst for one week
  k=exp(-keep_tag_inst)  # probability for keepin the tag for 1 week
  loose_tagP~dbeta(12,68) 

  # Priors for reporting rates of tags based on expert elicitation
  p_rep[1]=S_X[Y]
  p_rep[2]=R_X[Y]
  
  Y~dcat(w[1:3])
  
  for(i in 1:3){
    w[i]=1/3 # Equal weights given to each expert
  }
  
  # Expert 1
  S_X[1]~dbeta(0.69*25,(1-0.69)*25)
  R_X[1]~dbeta(0.77*30,(1-0.77)*30)
  
  # Expert 2
  S_X[2]~dbeta(0.71*60, (1-0.71)*60)
  R_X[2]~dbeta(0.755*50, (1-0.755)*50)
  
  # Expert 2
  S_X[3]~dbeta(0.78*13, (1-0.78)*13)
  R_X[3]~dbeta(0.56*13, (1-0.56)*13)


}"


data=list(
  Tagged=tagged,
  R=recap,  
  T=15
)  


var_names<- c(
  #"haav_prop",  "haav_propP",
  "q_net",
  "Pdie", 
  "p_rep",
  "h", "F", "M",
  "F_sea", "F_river", 
  "M_sea",  "M_river",
  "HR", "move",
  "handling_mort", "loose_tag","loose_tagP")

# run0 <- run.jags(M2,
#                  monitor= var_names,data=data, #inits = inits,
#                  n.chains = 2, method = 'parallel', thin=10, burnin = 1000,
#                  modules = "mix",keep.jags.files=F,sample =10000, adapt = 1000,
#                  progress.bar=TRUE)
# 

run1 <- autorun.jags(M2,
                 monitor= var_names,data=data, #inits = inits,
                 n.chains = 2, method = 'parallel', thin=10, 
                 modules = "mix",progress.bar=TRUE)


run<-run1

sum_run<-round(summary(run),3);sum_run
filter(as.data.frame(sum_run), `MC%ofSD`>5)
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





# Tn kuolla tietystä rysätyypistä vapauttamisen vuoksi
##################################

summary(chains[,"Pdie_sea[1]"], quantiles=c(0.05,0.5,0.95)) #Kaukalo
summary(chains[,"Pdie_sea[2]"], quantiles=c(0.05,0.5,0.95)) #Sukka

# Tn kaukalo > sukka
summary(chains[,"Pdie_sea[1]"]/chains[,"Pdie_sea[2]"], quantiles=c(0.05,0.25,0.5,0.75,0.95))
#summary(chains[,"Pdie_river[1]"]/chains[,"Pdie_river[2]"], quantiles=c(0.11,0.25,0.5,0.75,0.95))

tmp<-c()
for(i in 1:length(chains[,"Pdie_sea[1]"])){
  tmp[i]<-ifelse(chains[,"Pdie_sea[1]"][i]>chains[,"Pdie_sea[2]"][i], 1,0)
}
mean(tmp)

par(mfrow=c(1,2))
plot(density(chains[,"Pdie_sea[1]"]),lwd=2, main="Vapautuskuolleisuus", ylim=c(0,17), xlab="Todennäköisyys kuolla")
lines(density(chains[,"Pdie_sea[2]"]), lty=2, lwd=2)
#plot(density(chains[,"Pdie_seaP"]), lwd=1) #Priori piikkaa ykköseen, mutta tasainen hännässä
legend("topright", lty=c(1,2), lwd=c(2,2), legend=c("Kaukalo", "Sukka"))

plot(density(chains[,"Pdie_sea[1]"]/chains[,"Pdie_sea[2]"]), main="Kuoll. kaukalo / kuoll. sukka", xlim=c(0,20), ylim=c(0,0.35), xlab="Osamäärä")

# eng ###########################
par(mfrow=c(1,2))
plot(density(chains[,"Pdie_sea[1]"]),lwd=2, main="Release mortality", ylim=c(0,17), xlab="Probability to die")
lines(density(chains[,"Pdie_sea[2]"]), lty=2, lwd=2)
legend("topright", lty=c(1,2), lwd=c(2,2), legend=c("Emptying chute", "Lifting bag"))
plot(density(chains[,"Pdie_sea[1]"]/chains[,"Pdie_sea[2]"]), main="Mortality by emptying chute / Mortality by lifting bag", xlim=c(0,20), ylim=c(0,0.35), xlab="Quotient")
########################################



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

# eng ###########################
par(mfrow=c(3,3))
plot(density(chains[,"F_seaP"]), main="F at coastal fisheries", xlab="Instantaneous mortality", ylab="Density", xlim=c(0,0.5))
lines(density(chains[,"F_sea"]),lwd=2)
plot(density(chains[,"F_riverP"]), main="F at river fisheries", xlab="Instantaneous mortality", ylab="Density", xlim=c(0,0.6))
lines(density(chains[,"F_river"]),lwd=2)
plot(density(chains[,"M_seaP"]), main="M at coastal fisheries", xlab="Instantaneous mortality", ylab="Density", xlim=c(0,0.5), ylim=c(0,10))
lines(density(chains[,"M_sea"]),lwd=2)
plot(density(chains[,"M_riverP"]), main="M at river fisheries", xlab="Instantaneous mortality", ylab="Density", xlim=c(0,0.5))
lines(density(chains[,"M_river"]),lwd=2)
plot(density(chains[,"handling_mortP"]), main="Handling mortality", xlab="Instantaneous mortality", ylab="Density", xlim=c(0,1), ylim=c(0,4))
lines(density(chains[,"handling_mort"]),lwd=2)
plot(density(chains[,"haav_propP"]), main="(Handling mortality) / (tagging + handling mortality)", xlab="Proportion", ylab="Density", xlim=c(0,1))
lines(density(chains[,"haav_prop"]),lwd=2)
plot(density(chains[,"loose_tagP"]), main="Prob. to loose a tag", xlab="Probability", ylab="Density", xlim=c(0,0.4), ylim=c(0,12))
lines(density(chains[,"loose_tag"]),lwd=2)
plot(density(chainsP[,"C_A"][[1]]), main="Prob. to return a tag at coastal fisheries", xlab="Probability", ylim=c(0,6))
lines(density(chains[,"reporting_sea"]), lwd=2)
plot(density(chainsP[,"R_A"][[1]]), main="Prob. to return a tag at river fisheries", xlab="Probability", ylim=c(0,6))
lines(density(chains[,"reporting_river"]), lwd=2)
########################################


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

