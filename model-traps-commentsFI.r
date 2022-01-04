model{
  for(j in 1:2){ # Rysätyyppi 1:Kaukalo; 2:Sukka
    for(t in 1:15){ # Aika viikkoina
      C[t,1,j]=harv_sea[j]*N[t,1,j]                  # C[,1,]: Rannikkopyynnissä saaliiksi saatu määrä merkittyjä kaloja
      R[t,1,j]~dpois(C[t,1,j]*reporting_sea+0.001)   # R[,1,]: Merkkipalautusten lukumäärä rannikkopyynnistä 
      
      C[t,2,j]=harv_river[j]*N[t,2,j]                # C[,2,]: Jokipyynnissä saaliiksi saatu määrä merkittyjä kaloja
      R[t,2,j]~dpois(C[t,2,j]*reporting_river+0.001) # R[,2,]: Merkkipalautusten lukumäärä jokipyynnistä 
      
    }
  }
  
  # Ensimmäinen viikko
  for(j in 1:2){
    N[1,1,j]=Tagged[1,j]  # Ensimmäinen merkintäerä (rannikko)
    N[1,2,j]=0            # Ei merkintää jokialueilla
    
    # N[,1,]: Kalojen lukumäärä rannikkoalueella viikoittain
    for(t in 1:(T-1)){
      N[t+1,1,j]=N[t,1,j]*s[t,1,j]+Tagged[t+1,j]
    }
    
    # N[,2,]: Kalojen lukumäärä jokialueilla viikoittain
    for(t in 1:T){
      N[t+1,2,j]=N[t,1,j]*m[t,j]+N[t,2,j]*s[t,2,j]
    }
    
    
    # s ja m: osuudet jotka selviytyvät hengissä seuraavaan viikkoon
    for(t in 1:T){
      s[t,1,j]=surv_sea[j]*keep_tag*(1-move) # osuus joka selviää hengissä, pitää merkkinsä ja pysyy rannikolla
      s[t,2,j]=surv_river[j]*keep_tag        # osuus joka selviää hengissä, pitää merkkinsä ja pysyy jokialueella
      m[t,j]=surv_sea[j]*keep_tag*move       # osuus joka selviää hengissä, pitää merkkinsä ja siirtyy rannikolta jokialueelle
    }
    
    # Pdie_sea: Vapautuskuolevuus rysätyypistä j
    # Huom! Kalastus- ja luonnollisten kuolevuuksien nimittäjissä olevat viikkomäärät (11, 12 tai 52) ovat seurausta
    # kunkin kalastuksen oletetusta kestosta viikkoina WGBAST-mallissa (ks. ICES 2021) 
    Pdie_sea[j]=(h[j]+hand_inst*haav_prop)/(M_sea/12+h[j]+F_sea/11+hand_inst*haav_prop)*
      (1-exp(-(M_sea/12+h[j]+F_sea/11+hand_inst*haav_prop)))
    
    Pdie_river[j]=(h[j]+hand_inst*haav_prop)/(M_river/52+h[j]+F_river/12+hand_inst*haav_prop)*
      (1-exp(-(M_river/52+h[j]+F_river/12+hand_inst*haav_prop)))
    
    # surv_sea: Eloonjäänti rannikoilla, kaikki kuolinsyyt hetkellisinä kuolevuuksina
    # M_sea/12: Luonnollinen kuolevuus rannikolla viikon aikana
    # F_sea/11: Kalastuskuolevuus rannikolla viikon aikana
    # hand_inst: Käsittelykuolevuus (ml. merkintä) viikon aikana
    # h: Kustakin rysätyypistä aiheutuva kuolevuus viikon aikana
    surv_sea[j]=exp(-(M_sea/12+h[j]+F_sea/11+hand_inst))
    
    # harv_sea: Rannikolla saaliiksi jäävä osuus merkittyjä kaloja viikon aikana
    # keep_tag: Merkin pysymisen todennäköisyys viikon aikana
    harv_sea[j]=keep_tag*(1-surv_sea[j])*(F_sea/11)/(M_sea/12+F_sea/11+h[j]+hand_inst) # harvest based on inst. mortalities
    
    # surv_river: Eloonjäänti jokialueilla, kaikki kuolinsyyt hetkellisinä kuolevuuksina
    # M_river/52: Luonnollinen kuolevuus jokialueilla viikon aikana
    # F_river/12: Kalastuskuolevuus jokialueilla viikon aikana
    surv_river[j]=exp(-(M_river/52+F_river/12+h[j]+hand_inst))
    
    # harv_river: Jokialueilla saaliiksi jäävä osuus merkittyjä kaloja viikon aikana
    harv_river[j]=keep_tag*(1-surv_river[j])*(F_river/12)/(M_river/52+F_river/12+h[j]+hand_inst)
    
    h[j]~dunif(0,12) # Rysäkohtaisen kuolevuuden priorijakauma
  }
  
  # WGBAST-malliin perustuvat priorijakaumat:
  F_sea~dlnorm(log(0.1566)-0.5/T_sea, T_sea) # Hetkellinen kalastuskuolevuus rannikkokalastuksessa 11 viikon aikana
  T_sea=1/log(pow(0.02062/0.1566,2)+1) # Hajontamuuttuja
  F_river~dlnorm(log(0.2119)-0.5/T_river, T_river) # Hetkellinen kalastuskuolevuus jokialueilla 12 viikon aikana
  T_river=1/log(pow(0.07043/0.2119,2)+1) # Hajontamuuttuja
  M_river~dlnorm(log(0.0975)-0.5/T_Mr, T_Mr) # (Yleinen) hetkellinen luonnollinen kuolevuus vuoden aikana 
  T_Mr=1/log(pow(0.0312/0.0975,2)+1) # Hajontamuuttuja
  M_sea~dlnorm(log(0.1622)-0.5/T_Mc, T_Mc) # Hetkellinen luonnollinen kuolevuus rannikkokalastuksen aikana, coastal inst M during 2 months in WGBAST model
  T_Mc<-1/log(pow(0.052/0.1622,2)+1) # Hajontamuuttuja
  
  
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
  
  
}