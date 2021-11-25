
library(runjags)

# Lisää aika, onko 20 viikkoa riittävä? Rannikkokalastus ja jokikalastus erikseen


M1<-"model{
  

  for(j in 1:2){ # number of gear types 1:Kaukalo; 2:Sukka
  for(i in 1:4){ #number of fishermen
    Tagged[i,j]~dbin((1-handling_mort),Tagged_o[i,j])
  
    C[i,j]=harv_tot[j]*Tagged[i,j]              # number of tagged fish caught
    R[i,j]~dpois(C[i,j]*reporting+0.001)    # number of tagged fish reported = data
  }

 harv_tot[j]=keep_tag*(1-surv[j])*F/(M/4+F+h[j]) # harvest based on inst. mortalities
 #harvest_sea[a,k]=keep_tag*(1-survive_sea[a,k])*F_sea[a]/(M_sea/52+F_sea[a]+h[k]) # harvest based on inst. mortalities
 
    surv[j]=exp(-(M+F+h[j])) 
    #HR[j]=1-surv[j]
    # surv_sea[j]=exp(-(M_sea+F_sea+h[j])) 
    # surv_river[j]=exp(-(M_river+F_river+h[j]))

    #harv_sea[j]=keep_tag*(1-surv_sea[j]) * F_sea/(M_sea+F_sea+h[j]) # harvest based on inst. mortalities
    #harv_river[j]=keep_tag*(1-surv_river[j])*F_river/(M_river/52+F_river+h[k])
    
    #h[j]~dunif(0,12)#trap induced mortality
    h[j]~dunif(0,20)#trap induced mortality
  }

  #F~dunif(0,100)
  
  F~dlnorm(log(0.319)-0.5/T_F, T_F) 
  T_F<-1/log(pow(0.0716/0.319,2)+1)
  #F_sea~dunif(0,100) # gulf of bothnia    
  #F_river~dunif(0,100) # river fishery
  M~dlnorm(log(0.0975)-0.5/T_M, T_M) #inst M based on WGBAST
  T_M<-1/log(pow(0.0312/0.0975,2)+1)

  handling_mort~dbeta(27,132) # from Siira et al 
  loose_tag~dbeta(12,68)  # from Siira et al, over three months
  keep_tag_inst=-log(1-loose_tag) # inst for three months
  keep_tag=exp(-keep_tag_inst)  # probability for keepin the tag for 3 months

  loose_tagP~dbeta(12,68) 
  handling_mortP~dbeta(27,132)


}"

data=list(#Tagged_o=c(202, 211),R=c(18,27),
Tagged_o=cbind(
  filter(input_dat, gear==1)$n_tagged,
  filter(input_dat, gear==2)$n_tagged),
R=cbind(
  filter(input_dat, gear==1)$n_return,
  filter(input_dat, gear==2)$n_return),

#M=0.1,
          #M_sea=0.1,           # Assuming known natural mortalities for now, should use a prior later
          #M_river=0.2,
          
reporting=0.8          
#reporting_river=0.8, # Assuming known reporting rates for now, should use a prior later
          #reporting_sea=0.5
)  


var_names<- c("h", #"F_sea", "F_river", 
              #"HR",
              "F","M", "harv_tot",
              "handling_mort", "handling_mortP","loose_tag","loose_tagP")

run0 <- run.jags(M1,
                 monitor= var_names,data=data, #inits = inits,
                 n.chains = 2, method = 'parallel', thin=10, burnin = 1000,
                 modules = "mix",keep.jags.files=F,sample =1000, adapt = 1000,
                 progress.bar=TRUE)


run<-run0

summary(run)
plot(run)

chains<-as.mcmc(run)

summary(chains[,"h[1]"])
summary(chains[,"h[1]"]/chains[,"h[2]"], quantiles=c(0.1,0.5,0.95))



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
    Tagged[t,j]~dbin((1-handling_mort),Tagged_o[t,j])

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
  
    surv_sea[j]=exp(-(M/52+h[j]+F_sea))  # survival based on instantenous mortality rates per week
                                                # M: natural mortality, h: trap induced mortality
                                                # F: fishing mortality
                                                
    harv_sea[j]=keep_tag*(1-surv_sea[j])*F_sea/(M/52+F_sea+h[j]) # harvest based on inst. mortalities
    
    surv_river[j]=exp(-(M/52+F_river+h[j]))
    harv_river[j]=keep_tag*(1-surv_river[j])*F_river/(M/52+F_river+h[j])
  
    h[j]~dunif(0,12)#trap induced mortality
  }

  F_sea~dlnorm(log(0.1566)-0.5/T_sea, T_sea) # sea F based on WGBAST
  T_sea<-1/log(pow(0.02062/0.1566,2)+1)

  F_river~dlnorm(log(0.2119)-0.5/T_river, T_river) # river F based on WGBAST
  T_river<-1/log(pow(0.07043/0.2119,2)+1)

  M~dlnorm(log(0.0975)-0.5/T_M, T_M) #inst M based on WGBAST
  T_M<-1/log(pow(0.0312/0.0975,2)+1)

  move~dbeta(4,1)             # weakly informative prior for movement probability 

  handling_mort~dbeta(27,132) # from Siira et al 
  loose_tag~dbeta(12,68)  # from Siira et al, over three months
  keep_tag_inst=-log(1-loose_tag) # inst for three months
  keep_tag=exp(-keep_tag_inst)  # probability for keepin the tag for 3 months

  loose_tagP~dbeta(12,68) 
  handling_mortP~dbeta(27,132)
  
  F_seaP~dlnorm(log(0.1566)-0.5/T_seaP, T_seaP) 
  T_seaP<-1/log(pow(0.02062/0.1566,2)+1)

  F_riverP~dlnorm(log(0.2119)-0.5/T_riverP, T_riverP) 
  T_riverP<-1/log(pow(0.07043/0.2119,2)+1)

  MP~dlnorm(log(0.0975)-0.5/T_MP, T_MP) 
  T_MP<-1/log(pow(0.0312/0.0975,2)+1)


}"


data=list(
  Tagged_o=tagged,
  R=recap,  
  reporting_river=0.8, # Assuming known reporting rates for now, should use a prior later
  reporting_sea=0.5,
  T=15
)  


var_names<- c("h", "F_sea", "F_river", 
              "M", "harv_sea", "harv_river", "move",
              "handling_mort", "handling_mortP","loose_tag","loose_tagP",
              "F_seaP", "F_riverP")

run0 <- run.jags(M2,
                 monitor= var_names,data=data, #inits = inits,
                 n.chains = 2, method = 'parallel', thin=10, burnin = 1000,
                 modules = "mix",keep.jags.files=F,sample =1000, adapt = 1000,
                 progress.bar=TRUE)


run<-run0

round(summary(run),3)
plot(run)

chains<-as.mcmc(run)

summary(chains[,"h[1]"])
summary(chains[,"h[1]"]/chains[,"h[2]"], quantiles=c(0.05,0.5,0.75,0.80))







  
# Original model, somwhat
model{   
  ## Movement of fish through space & time
  #==========================================
  
  # Account for handling mortality
  for(t in 1:(T+1)){
    for(a in 1:(A+1)){
      
      # Number of tagged fish released: Tagged_o
      # Number of tagged fish alive after handling mortality
      
      Tagged[t,a,j]~dbin((1-handling_mort),Tagged_o[t,a,k])
    }
  }
  
  
  # First time step: the initial release, all areas except river
  for(a in 1:A){
    N[1,a,k]=Tagged[1,a,k]    
  }
  N[1,A+1,k]=0                # no releases in the river (A+1)
  
  # Number of fish at area 1 through time 
  for(t in 1:T){
    N[t+1,1,k]=N[t,1,k]*s[t,1,k]+Tagged[t+1,1,k]
  }
  
  # Number of fish at areas 2:A+1 through time
  for(k in 1:K){       # K different tagging places
  for(t in 1:T){       # T+1 time steps = 20 weeks
    for(a in 1:A){     # A+1 different recapture locations
      
      # N: number of fish alive in the beginning of time step t, at recap area a, 
      # originated from tagging place k
      
      N[t+1,a+1,k]=N[t,a,k]*m[t,a,k]+N[t,a+1,k]*s[t,a+1,k]+Tagged[t+1,a+1,k]
      }
  }

  
  for(t in 1:T){
    
    # surviving fractions, these should be made stochastic in future, using dbeta()
    
    for(a in 1:A){
      s[t,a,k]=survive_sea[a,k]*keep_tag*(1-move) # surviving fraction that stays in area a
      m[t,a,k]=survive_sea[a,k]*keep_tag*move     # surviing fraction that moves to next area
    }
    s[t,(A+1),k]=survive_river[k]*keep_tag
    m[t,(A+1),k]=0
  }
  }
  for(k in 1:K){
  for(a in 1:A){
  survive_sea[a,k]=exp(-(M_sea/52+h[k]+F_sea[a]))   # survival based on instantenous mortality rates per week
                                                    # M: natural mortality, h: trap induced mortality
                                                    # F: fishing mortality
  
  harvest_sea[a,k]=keep_tag*(1-survive_sea[a,k])*F_sea[a]/(M_sea/52+F_sea[a]+h[k]) # harvest based on inst. mortalities
  }
  survive_river[k]=exp(-(M_river/52+F_river+h[k]))
  

  for(t in 1:(T+1)){
    for(a in 1:A){
      C[t,a,k]=harvest_sea[a,k]*N[t,a,k]              # number of tagged fish caught
      R[t,a,k]~dpois(C[t,a,k]*reporting_sea+0.001)    # number of tagged fish reported = data
    }
    C[t,(A+1),k]=harvest_river[k]*N[t,(A+1),k]
    R[t,(A+1),k]~dpois(C[t,(A+1),k]*reporting_river+0.001)
  }
  
 
  harvest_river[k]=keep_tag*(1-survive_river[k])*F_river/(M_river/52+F_river+h[k])
  
    h[k]~dunif(0,52)
  }
  
  F_sea[1]~dunif(0,100)      # vague priors for fishing mortality at different stages
  F_sea[2]=F_sea[1]          # assuming equal fishing mortality across the coast, except
  F_sea[3]~dunif(0,100)      # in the Bothnian Bay
  F_river~dunif(0,100)
  move~dbeta(4,1)             # weakly informative prior for movement probability 
  handling_mort~dbeta(27,132) # from Siira et al 
  loose_tag~dbeta(12,68)  # from Siira et al, over three months = 12 weeks
  keep_tag_inst=-log(1-loose_tag)/12 # inst for one week
  keep_tag=exp(-keep_tag_inst)  # weekly pronbability for keepin the tag
}

