
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
  
# Priors:

pow<-function(x,z){
  y=x^z
  return(y)
}


T_Fs=1/log(pow(0.02062/0.1566,2)+1);1/sqrt(T_Fs)
log(0.1566)-0.5/T_Fs
#F_sea~dlnorm(log(0.1566)-0.5/T_Fs, T_Fs) # sea F, 11 weeks 
F_sea~logN(-1.863,0.1311^2)

T_Fr=1/log(pow(0.07043/0.2119,2)+1);1/sqrt(T_Fr)
log(0.2119)-0.5/T_Fr
#F_river~dlnorm(log(0.2119)-0.5/T_Fr, T_Fr) # river F, 12 weeks 
F_river~logN(-1.604,0.3237^2)
  
  
T_Mr=1/log(pow(0.0312/0.0975,2)+1);1/sqrt(T_Mr)
log(0.0975)-0.5/T_Mr
#M_river~dlnorm(log(0.0975)-0.5/T_Mr, T_Mr) # river M, 52 weeks
M_river~logN(-2.377,0.3122^2)
  
  
T_Ms<-1/log(pow(0.052/0.1622,2)+1);1/sqrt(T_Ms)
log(0.1622)-0.5/T_Ms
#M_sea~dlnorm(log(0.1622)-0.5/T_Ms, T_Ms) # sea  M, 8 weeks
M_sea~logN(-1.868,0.3128^2)

M3<-"model{
  for(t in 1:T){
      for(j in 1:2){ # Gear (1: pontoon trap with an emptying chute, 2: PU trap with a lifting bag)
        N_alive[t,j]~dbin(p_alive[j],N_tagged[t,j])

      for(i in 1:2){ # Area (1: sea, 2: river)  
        N_recap[t,i,j]~dpois(HR[i]*n[t,i,j]*r_rep[i]+0.001)

      }
    }
  }  
    
  
  for(j in 1:2){

    p_alive[j]=exp(-(h[j]+ip_mark))
    p_die_tmp[j]=1-exp(-(h[j]+ip_mark))
    p_die[j]=1-exp(-(h[j]+ip_mark*q_hand))

    n[1,1,j]=N_alive[1,j]    
    n[1,2,j]=0        # no releases in the river
  
    # Number of fish at sea through time 
    # s:surviving that stays at the area
    for(t in 1:(T-1)){
      n[t+1,1,j]=n[t,1,j]*q_SS[1]+N_alive[t+1,j]
    }
 
    # Number of fish at river through time 
    # m:surviving that moves from sea to river
    # s:surviving that stays at the area
    for(t in 1:T){
      n[t+1,2,j]=n[t,1,j]*q_SM+n[t,2,j]*q_SS[2]
    }
  } 

  # Surviving fractions  
  q_SS[1]=p_surv[1]*p_keep*(1-p_move) # survives and stays at sea
  q_SS[2]=p_surv[2]*p_keep        # survives and that stays at river
  q_SM=p_surv[1]*p_keep*p_move   # survives and moves from sea to river
  
  for(i in 1:2){
    p_surv[i]=exp(-(M[i]+F[i]))  # survival based on instantenous mortality rates per week
    HR[i]=p_keep*(1-p_surv[i])*F[i]/(M[i]+F[i]) # harvest based on inst. mortalities
  }
  
  # Uninformative priors
  p_move~dunif(0,1) # probability to move from sea to river, 1 week

  for(j in 1:2){
    h[j]~dunif(0,12) # trap induced mortality
  }

  # Informative priors based on Baltic salmon stock assessment model
  F_sea~dlnorm(log(0.1566)-0.5/T_Fs, T_Fs) # sea F, 11 weeks 
  T_Fs=1/log(pow(0.02062/0.1566,2)+1)

  F_river~dlnorm(log(0.2119)-0.5/T_Fr, T_Fr) # river F, 12 weeks 
  T_Fr=1/log(pow(0.07043/0.2119,2)+1)

  M_river~dlnorm(log(0.0975)-0.5/T_Mr, T_Mr) # river M, 52 weeks
  T_Mr=1/log(pow(0.0312/0.0975,2)+1)
  
  M_sea~dlnorm(log(0.1622)-0.5/T_Ms, T_Ms) # sea  M, 8 weeks
  T_Ms<-1/log(pow(0.052/0.1622,2)+1)

  # Instantaneous mortalities per week
  M[1]=M_sea/8
  M[2]=M_river/52
  F[1]=F_sea/11
  F[2]=F_river/12

  # Other informative priors

  p_mark~dbeta(1.8,7.2) # mortality due to tagging + handling, 12 weeks
  ip_mark=-log(1-p_mark) # instantaneous tagging + handling mortality, 12 weeks

# Share of handling mortality in total tagging + handling
  q_hand~dlnorm(-0.723,7.2485)T(,1)
  
  p_lose~dbeta(12,68)  # probability to loose a tag, 12 weeks
  ip_keep=-log(1-p_lose/12) # prob to keep tag for one week, instantaneous scale
  p_keep=exp(-ip_keep)  # probability to keep a tag for 1 week
  
  # Priors for reporting rates of tags based on expert elicitation
  r_rep[1]=E[Y,1] # sea fisheries
  r_rep[2]=E[Y,2] # river fisheries
  
  Y~dcat(w[1:3])
  
  for(i in 1:3){
    w[i]=1/3 # Equal weights given to each expert
  }
  
  # Expert 1
  E[1,1]~dbeta(17.25,7.75)
  E[1,2]~dbeta(23.1,6.9)

  # Expert 2
  E[2,1]~dbeta(42.6,17.4)
  E[2,2]~dbeta(37.75,12.25)

  # Expert 3
  E[3,1]~dbeta(10.14,2.86)
  E[3,2]~dbeta(7.28,5.72)



}"


data=list(
  N_tagged=tagged,
  N_recap=recap,  
  T=15
)  


var_names<- c(
  "q_hand",
  "p_alive", 
  "p_die", 
  #"p_die_tmp", 
  "r_rep",
  "h", "F", "M",
  "F_sea", "F_river", 
  "M_sea",  "M_river",
  "HR", "p_move",
  "p_mark", "p_lose", "p_keep")

run1 <- run.jags(M3,
                 monitor= var_names,data=data, #inits = inits,
                 n.chains = 2, method = 'parallel', thin=10, sample=30000,
                 modules = "mix",progress.bar=TRUE)


run<-run1
save(run, file=str_c("chains-traps-final.RData"))

sum_run<-round(summary(run, confidence=c(0.9)),3);sum_run
filter(as.data.frame(sum_run), `MC%ofSD`>5)
plot(run)
chains<-as.mcmc(run)
summary(chains, quantiles=c(0.05,0.5,0.95) )





# Priors
prior<-"
model{

 p_move~dunif(0,1) # probability to move from sea to river, 1 week

  for(j in 1:2){
    h[j]~dunif(0,12) # trap induced mortality
  }

  # Informative priors based on Baltic salmon stock assessment model
  F_sea~dlnorm(log(0.1566)-0.5/T_Fs, T_Fs) # sea F, 11 weeks 
  T_Fs=1/log(pow(0.02062/0.1566,2)+1)

  F_river~dlnorm(log(0.2119)-0.5/T_Fr, T_Fr) # river F, 12 weeks 
  T_Fr=1/log(pow(0.07043/0.2119,2)+1)

  M_river~dlnorm(log(0.0975)-0.5/T_Mr, T_Mr) # river M, 52 weeks
  T_Mr=1/log(pow(0.0312/0.0975,2)+1)
  
  M_sea~dlnorm(log(0.1622)-0.5/T_Ms, T_Ms) # sea  M, 8 weeks
  T_Ms<-1/log(pow(0.052/0.1622,2)+1)

  p_mark~dbeta(1.8,7.2) # mortality due to tagging + handling, 12 weeks
  q_hand~dlnorm(-0.723,7.2485)T(,1)
  p_lose~dbeta(12,68)  # probability to loose a tag, 12 weeks

  # Priors for reporting rates of tags based on expert elicitation
  r_rep[1]=E[Y,1] # sea fisheries
  r_rep[2]=E[Y,2] # river fisheries
  
  Y~dcat(w[1:3])
  
  for(i in 1:3){
    w[i]=1/3 # Equal weights given to each expert
  }
  

 # Expert 1
  E[1,1]~dbeta(17.25,7.75)
  E[1,2]~dbeta(23.1,6.9)

  # Expert 2
  E[2,1]~dbeta(42.6,17.4)
  E[2,2]~dbeta(37.75,12.25)

  # Expert 3
  E[3,1]~dbeta(10.14,2.86)
  E[3,2]~dbeta(7.28,5.72)


}"



var_namesP<-c(
  "q_hand",
  "r_rep",
  "h", 
  "F_sea", "F_river", 
  "M_sea",  "M_river",
  "p_move",
  "p_mark", "p_lose", "p_keep",
  "r_rep", "E")


runP <- autorun.jags(prior,
                     monitor= var_namesP,#data=data, #inits = inits,
                     n.chains = 2, method = 'parallel', modules = "mix",
                     progress.bar=TRUE)


chainsP<-as.mcmc(runP)


#summary(chainsP, quantiles=c(0.025,0.5,0.975))

par(mfrow=c(1,2))
plot(density(chainsP[,"r_rep[1]"]), ylim=c(0,7), main="Sea fisheries", xlab="Reporting rate")
lines(density(chainsP[,"E[1,1]"]), lty=2, col="red")
lines(density(chainsP[,"E[2,1]"]), lty=3, col="green")
lines(density(chainsP[,"E[3,1]"]), lty=4, col=4)
lines(density(chains[,"r_rep[1]"]), lty=1, lwd=2)
legend("topleft", legend=c("Expert 1","Expert 2","Expert 3", "Average"), lty=c(2,3,4, 1), col=c("red", "green", 4, 1))


plot(density(chainsP[,"r_rep[2]"]), ylim=c(0,8), main="River fisheries", xlab="Reporting rate")
lines(density(chainsP[,"E[1,2]"]), lty=2, col="red")
lines(density(chainsP[,"E[2,2]"]), lty=3, col="green")
lines(density(chainsP[,"E[3,2]"]), lty=4, col=4)
lines(density(chains[,"r_rep[2]"]), lty=1, lwd=2)
legend("topleft", legend=c("Expert 1","Expert 2","Expert 3", "Average"), lty=c(2,3,4, 1), col=c("red", "green", 4, 1))

# eng ###########################
par(mfrow=c(4,3))
plot(density(chainsP[,"F_sea"]), main=expression(F[sea]), xlab="Instantaneous mortality", ylab="Density", xlim=c(0,0.5))
lines(density(chains[,"F_sea"]),lwd=2)
plot(density(chainsP[,"F_river"]), main=expression(F[river]), xlab="Instantaneous mortality", ylab="Density", xlim=c(0,0.6), ylim=c(0,7))
lines(density(chains[,"F_river"]),lwd=2)
plot(density(chainsP[,"M_sea"]), main=expression(M[sea]), xlab="Instantaneous mortality", ylab="Density", xlim=c(0,0.5), ylim=c(0,10))
lines(density(chains[,"M_sea"]),lwd=2)
plot(density(chainsP[,"M_river"]), main=expression(M[river]), xlab="Instantaneous mortality", ylab="Density", xlim=c(0,0.5))
lines(density(chains[,"M_river"]),lwd=2)
plot(density(chainsP[,"h[1]"]), main=expression(h[1]), xlab="Instantaneous mortality", ylim=c(0,3), xlim=c(0,1))
lines(density(chains[,"h[1]"]), lwd=2)
plot(density(chainsP[,"h[2]"]), main=expression(h[2]), xlab="Instantaneous mortality", ylim=c(0,8), xlim=c(0,1))
lines(density(chains[,"h[2]"]), lwd=2)
plot(density(chainsP[,"p_mark"]), main=expression(p^mark), xlab="Instantaneous mortality", ylab="Density", xlim=c(0,0.6), ylim=c(0,10))
lines(density(chains[,"p_mark"]),lwd=2)
plot(density(chainsP[,"p_move"]), main=expression(p^move), xlab="Probability", ylim=c(0,8), xlim=c(0,0.6))
lines(density(chains[,"p_move"]), lwd=2)
plot(density(chainsP[,"p_lose"]), main=expression(p^lose), xlab="Probability", ylab="Density", xlim=c(0,0.4), ylim=c(0,12))
lines(density(chains[,"p_lose"]),lwd=2)
plot(density(chainsP[,"q_hand"]), main=expression(q^hand), xlab="Proportion", ylab="Density", xlim=c(0,1))
lines(density(chains[,"q_hand"]),lwd=2)
plot(density(chainsP[,"r_rep[1]"]), main=expression(r^rep), xlab="Probability", ylim=c(0,5), xlim=c(0.2,1))
lines(density(chains[,"r_rep[1]"]), lwd=2)
plot(density(chainsP[,"r_rep[2]"]), main=expression(r^rep), xlab="Probability", ylim=c(0,5), xlim=c(0.2,1))
lines(density(chains[,"r_rep[2]"]), lwd=2)




########################################


# Rysäkohtainen kuolevuus
#===========================
summary(chains[,"p_die[1]"]/chains[,"p_die[2]"], quantiles=c(0.05, 0.22,0.5,0.95))
# 78% tn:llä sukkaeloonjäänti yhtä hyvä tai parempi kuin kaukaloeloonjäänti
# 50% tn:llä vähintään 1.9 kertainen

par(mfrow=c(1,1))
plot(density(chains[,"p_die[1]"]),lwd=2, main="Vapautuskuolleisuus", ylim=c(0,8), xlab="Todennäköisyys")
lines(density(chains[,"p_die[2]"]), lty=2, lwd=2)
legend("topright", lty=c(1,2), lwd=c(2,2), legend=c("Kaukalo", "Sukka"))

plot(density(chains[,"p_die[1]"]),lwd=2, main="Release mortality", ylim=c(0,8), xlab="Probability")
lines(density(chains[,"p_die[2]"]), lty=2, lwd=2)
legend("topright", lty=c(1,2), lwd=c(2,2), legend=c("Emptying chute", "Lifting bag"))


abline(v=mean(chains[,"p_die[1]"]))
abline(v=mean(chains[,"p_die[2]"]))


summary((1-chains[,"p_die[2]"])/(1-chains[,"p_die[1]"]), quantiles=c(0.05, 0.20,0.5,0.95))


summary(chains[,"h[1]"], quantiles=c(0.05,0.5,0.95))
summary(chains[,"h[2]"], quantiles=c(0.05,0.5,0.95))

# Absoluuttiseksi kääntäminen mielekästä vain jos ajatellaan että muita kuolevuuksia ei olisi
summary(1-exp(-chains[,"h[1]"]), quantiles=c(0.05,0.5,0.95))
summary(1-exp(-chains[,"h[2]"]), quantiles=c(0.05,0.5,0.95))

summary((1-exp(-chains[,"h[1]"]))/(1-exp(-chains[,"h[2]"])), quantiles=c(0.05,0.5,0.95))





par(mfrow=c(2,2))
plot(density(chains[,"reporting_seaP"]), main="reporting at coastal fishery")
lines(density(chains[,"reporting_sea"]), lwd=2)

plot(density(chains[,"reporting_riverP"]), main="reporting at river fishery")
lines(density(chains[,"reporting_river"]), lwd=2)

plot(density(chains[,"F_seaP"]), main="F at coastal fishery")
lines(density(chains[,"F_sea"]), lty=2)

plot(density(chains[,"F_riverP"]), main="F at river fishery")
lines(density(chains[,"F_river"]), lty=2)







# Priors vs posteriors
par(mfrow=c(3,3))
plot(density(chainsP[,"F_seaP"]), main="Kalastuskuol., rannikko", xlab="Hetkellinen kuolleisuus", ylab="Tiheys", xlim=c(0,0.5))
lines(density(chains[,"F_sea"]),lwd=2)
plot(density(chainsP[,"F_riverP"]), main="Kalastuskuol., jokialue", xlab="Hetkellinen kuolleisuus", ylab="Tiheys", xlim=c(0,0.6))
lines(density(chains[,"F_river"]),lwd=2)
plot(density(chainsP[,"M_seaP"]), main="Luonnollinen kuol., rannikko", xlab="Hetkellinen kuolleisuus", ylab="Tiheys", xlim=c(0,0.5), ylim=c(0,10))
lines(density(chains[,"M_sea"]),lwd=2)
plot(density(chainsP[,"M_riverP"]), main="Luonnollinen kuol., jokialue", xlab="Hetkellinen kuolleisuus", ylab="Tiheys", xlim=c(0,0.5))
lines(density(chains[,"M_river"]),lwd=2)

#par(mfrow=c(2,2))
plot(density(chainsP[,"handling_mortP"]), main="Käsittelykuolleisuus", xlab="Hetkellinen kuolleisuus", ylab="Tiheys", xlim=c(0,1), ylim=c(0,4))
lines(density(chains[,"handling_mort"]),lwd=2)
plot(density(chainsP[,"haav_propP"]), main="Haavinnan osuus käsittelykuolleisuudesta", xlab="Osuus", ylab="Tiheys", xlim=c(0,1))
lines(density(chains[,"haav_prop"]),lwd=2)
plot(density(chainsP[,"loose_tagP"]), main="Merkin irtoamisen tn", xlab="Todennäköisyys", ylab="Tiheys", xlim=c(0,0.4), ylim=c(0,12))
lines(density(chains[,"loose_tag"]),lwd=2)

plot(density(chainsP[,"C_A"]), main="Raportointiaktiivisuus, rannikko", xlab="Todennäköisyys", ylim=c(0,6))
lines(density(chains[,"reporting_sea"]), lwd=2)

plot(density(chainsP[,"R_A"]), main="Raportointiaktiivisuus, jokialueet", xlab="Todennäköisyys", ylim=c(0,6))
lines(density(chains[,"reporting_river"]), lwd=2)

