require(rjags)
require(coda)


# instantaneous M should have a prior roughly
#Mean             SD     
#0.0975      0.0312       
#5%     25%     50%     75%     90%     95% 
#0.049 0.073 0.0973    0.121   0.141   0.153 


M1<-"
model{  
M~dlnorm(log(0.0975)-0.5/T, T) #inst M
T<-1/log(pow(0.0312/0.0975,2)+1)


}"

cat(M1,file="prior.txt")


system.time(jm<-jags.model('prior.txt',n.adapt=100,n.chains=1))


system.time(chains1<-coda.samples(jm,
                                  variable.names=c(
                                    "M"
                                  ),
                                  n.iter=10000,
                                  thin=1))

chainsM<-chains1
summary(chainsM, quantiles=c(0.05,0.5,0.95))


# F (=-log(1-HR)) in coastal + river fishery together
#Mean             SD        
#0.318964       0.071594       

#  5%    25%    50%    75%    90%    95% 
#  0.2182 0.2708 0.3103 0.3579 0.4157 0.4481 

M1<-"
model{  
F~dlnorm(log(0.319)-0.5/T, T) 
T<-1/log(pow(0.0716/0.319,2)+1)


}"

cat(M1,file="prior.txt")


system.time(jm<-jags.model('prior.txt',n.adapt=100,n.chains=1))


system.time(chains1<-coda.samples(jm,
                                  variable.names=c(
                                    "F"
                                  ),
                                  n.iter=10000,
                                  thin=1))

chainsM<-chains1
summary(chainsM, quantiles=c(0.05,0.5,0.95))


# Joki- ja rannikko HR erikseen
# Aja WGBAST-projektin F42310_HR-skripti riville 71 saakka
# Vuodet 1:32 == 1989-2020, jolloin 
#hrW[2,30:32,]; hr_coastW[2,30:32,]
# on MSW HR vuosina 2018-2020 joelle ja rannikolle 


# Aja seuraavat summary-komennot em. projektissa
#summary(as.mcmc(apply(hrW[2,30:32,],2, mean)), quantiles=c(0.05,0.25,0.5,0.75,0.95))
# Mean             SD      
# 0.188997       0.056039     
#     5%    25%    50%    75%    95% 
# 0.1053 0.1500 0.1846 0.2240 0.2863 

#summary(as.mcmc(apply(hr_coastW[2,30:32,],2, mean)), quantiles=c(0.05,0.25,0.5,0.75,0.95))
# Mean           SD        
# 0.1448100      0.0175908       
#5%    25%    50%    75%    95% 
#0.1164 0.1326 0.1451 0.1554 0.1762

# Hetkellisenä  F (=-log(1-HR)):

# Jokikalastus
#F_R<--log(1-as.mcmc(apply(hrW[2,30:32,],2, mean)))
#summary(F_R, quantiles=c(0.05,0.25,0.5,0.75,0.95))
#Mean             SD       
#0.211929       0.070434       0.001819       0.002590 
#  5%    25%    50%    75%    95% 
#0.1112 0.1626 0.2041 0.2536 0.3373 

# Rannikkokalastus
#F_C<--log(1-as.mcmc(apply(hr_coastW[2,30:32,],2, mean)))
#summary(F_C, quantiles=c(0.05,0.25,0.5,0.75,0.95))
#Mean             SD      
#0.1566437      0.0206177      
#  5%    25%    50%    75%    95% 
#0.1238 0.1423 0.1568 0.1688 0.1938 



# Jokikalastus
M3<-"
model{  
F~dlnorm(log(0.2119)-0.5/T, T) 
T<-1/log(pow(0.07043/0.2119,2)+1)


}"

cat(M3,file="prior.txt")


system.time(jm<-jags.model('prior.txt',n.adapt=100,n.chains=1))


system.time(chains1<-coda.samples(jm,
                                  variable.names=c(
                                    "F"
                                  ),
                                  n.iter=10000,
                                  thin=1))

chainsM<-chains1
summary(chainsM, quantiles=c(0.05,0.5,0.95))

# Rannikkokalastus
M3<-"
model{  
F~dlnorm(log(0.1566)-0.5/T, T) 
T<-1/log(pow(0.02062/0.1566,2)+1)


}"

cat(M3,file="prior.txt")


system.time(jm<-jags.model('prior.txt',n.adapt=100,n.chains=1))


system.time(chains1<-coda.samples(jm,
                                  variable.names=c(
                                    "F"
                                  ),
                                  n.iter=10000,
                                  thin=1))

chainsM<-chains1
summary(chainsM, quantiles=c(0.05,0.5,0.95))




