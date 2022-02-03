
model{
  
  # x_1sw[i]: Number of dead 1SW salmon in group i 
  # N_1sw[i]: Total number of 1SW salmon in group i
  # x_Msw[i]: Number of dead MSW salmon in group i 
  # N_Msw[i]: Total number of MSW salmon in group i
  for(i in 1:2){ # 1: Tagged and handled; 2: handled (not tagged)
    x_1sw[i]~dbin(p_1sw[i],N_1sw[i])
    p_1sw[i]<-1-exp(-inst_1sw[i]) # finite mortality of tagged 1SW salmon
    x_msw[i]~dbin(p_msw[i],N_msw[i])
    p_msw[i]<-1-exp(-inst_msw[i])	# finite mortality of tagged MSW salmon
    
  }
  # Instantaneous mortalities:
  inst_1sw[1]=instTag_1sw+instHan_1sw
  inst_1sw[2]=instHan_1sw
  inst_msw[1]=instTag_msw+instHan_msw
  inst_msw[2]=instHan_msw
  
  instTag_1sw<- -log(survTag_1sw)
  instHan_1sw<- -log(survHan_1sw) 
  instTag_msw<- -log(survTag_msw) 
  instHan_msw<- -log(survHan_msw)
  
  # Prior distributions
  survTag_1sw~dbeta(1,1) 
  survHan_1sw~dbeta(1,1)  
  survTag_msw~dbeta(1,1) 
  survHan_msw~dbeta(1,1) 
  
  # The estimated proportion of instantaneous handling out 
  # of combined tagging and handling 
  q_1sw<-inst_1sw[2]/inst_1sw[1] 
  q_msw<-inst_msw[2]/inst_msw[1] 
  
  
}

# Data
data<-list(x_1sw=c(14,20), N_1sw=c(40,221),
           x_msw=c(8,24), N_msw=c(40,221))
