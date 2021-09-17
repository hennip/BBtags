model{
  
  for(i in 1:4){ # number of fishermen
    for(j in 1:2){ # number of gear types
      Tagged[i,j]~dbin((1-handling_mort),Tagged_o[i,j])
    }
  }

  s[i,j]=survive_sea[j]*keep_tag
  
  ## Movement of fish through space & time
  #==========================================
  
  # Account for handling mortality
  for(t in 1:(T+1)){
    for(a in 1:(A+1)){
      
      # Number of tagged fish released: Tagged_o
      # Number of tagged fish alive after handling mortality
      
      Tagged[t,a,k]~dbin((1-handling_mort),Tagged_o[t,a,k])
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

