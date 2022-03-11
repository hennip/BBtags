require(rjags)
require(coda)


# instantaneous M at the coastal fishery (2 months in WGBAST-model)
#should have a prior roughly
# Mean             SD        
# 0.162237       0.051923 
#     5%     25%     50%     75%     90%     95% 
#  0.08216 0.12066 0.16195 0.20050 0.23469 0.25385 
0.051923/0.162237


M1<-"
model{  
#M~dlnorm(log(0.0975)-0.5/T, T) #inst M during river fishery
#T<-1/log(pow(0.0312/0.0975,2)+1)


M~dlnorm(log(0.1622)-0.5/T, T) #inst M during coastal fishery
T<-1/log(pow(0.052/0.1622,2)+1)


}"

cat(M1,file="prior.txt")


system.time(jm<-jags.model('prior.txt',n.adapt=100,n.chains=1))


system.time(chains1<-coda.samples(jm,
                                  variable.names=c(
                                    "M"
                                  ),
                                  n.iter=10000,
                                  thin=1))

summary(chains1, quantiles=c(0.05,0.5,0.95))

plot(density(as.mcmc(chains1[,"M"])))


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


# Käsittelykuolevuus, ml. merkintä 
# Aikaisempi analyysi (Ruokonen et al)

#Mean             SD       Naive SE Time-series SE 
#0.198994       0.180809       0.004668       0.009740 
#  5%      25%      50%      75%      90%      95% 
#  -0.11920  0.09339  0.21339  0.32926  0.41562  0.46088 

#handling_mort~dbeta(27,132) # from Siira et al 

# Haavinnan osuus kokonaiskäsittelykuolevuudesta (haavinta+merkintä =1)
# Ks. snapshot_BB_allaskoe
#> summary(as.mcmc(q_markMSW), quantiles=c(0.05,0.25,0.5,0.75,0.95))
#Mean             SD       Naive SE Time-series SE 
#0.509953       0.173782       0.004487       0.004099 

#  5%    25%    50%    75%    95% 
#0.2642 0.3767 0.4862 0.6188 0.8347 
# Obs! Beta(1,1) prior used in the one above!




M4<-"
model{  

hand_siira~dbeta(27,132) # from Siira et al 

handlingM~dbeta(a, b)
a<-mu*eta
b<-(1-mu)*eta

mu<-0.2
eta<-9

handlingM2~dbeta(1.8,7.2)

inst_hand<--log(1-handlingM2)

# #haav_prop[1]~dbeta(4.53, 5.47)
# haav_prop[1]~dbeta(a2, b2)
# a2<-mu2*eta2
# b2<-(1-mu2)*eta2
# 
# mu2<-0.453
# eta2<-10

# haav_prop~dlnorm(M,T)I(,1)
# M<-log(mu3)-0.5/T
# T<-1/log(cv*cv+1)
# cv<-0.2/mu3
# mu3<-0.52
haav_prop~dlnorm(-0.723,7.2485)I(,1)


}"

cat(M4,file="priorM4.txt")

# 
# mu3<-0.52
# cv<-0.2/mu3
#  Tau<-1/log(cv*cv+1)
#  M<-log(mu3)-0.5/Tau
#  M;Tau

system.time(jm<-jags.model('priorM4.txt',n.adapt=100,n.chains=1))


system.time(chainsM<-coda.samples(jm,
                                  variable.names=c(
                                    "hand_siira",
                                    "handlingM", "handlingM2",
                                    "haav_prop"
                                  ),
                                  n.iter=100000,
                                  thin=1))
summary(chainsM, quantiles=c(0.05,0.25,0.5,0.75,0.95))


# Raportointiaktiivisuudet rannikko- ja jokikalastuksessa:
# Elisitointi kolmelta ekspertiltä (Petri, Timo, Tapani)
# min max mode/median, min-max vastaa 95% tn-väliä 

#Rannikkokalastus:

# Timo: 0.7 (0.5,0.85)
# Petri: 0.7 (0.55,0.80)
# Tapsa: moodi=0.85 95%PI=(0.5,0.95) 



M3<-"
model{  
x[1]~dbeta(mu*eta, (1-mu)*eta) 
mu<-0.69
eta<-25

x[2]~dbeta(0.71*60, (1-0.71)*60)


x[3]~dbeta(0.78*13, (1-0.78)*13)


}"

cat(M3,file="prior.txt")


system.time(jm<-jags.model('prior.txt',n.adapt=100,n.chains=1))

system.time(chains1<-coda.samples(jm,variable.names=c("x"),n.iter=1000000,thin=1))

summary(chains1, quantiles=c(0.025,0.05,0.5,0.95,0.975))
par(mfrow=c(2,2))
plot(density(chains1[,"x[1]"][[1]]), main="Rannikkokalastus, Timo")
abline(v=c(0.5,0.7,0.85))
plot(density(chains1[,"x[2]"][[1]]), main="Rannikkokalastus, Petri")
abline(v=c(0.55,0.72,0.8))
plot(density(chains1[,"x[3]"][[1]]), main="Rannikkokalastus, Tapsa")
abline(v=c(0.5,0.85,0.95))

summary(chains1[,"x[2]"], quantiles=c(0.025,0.5,0.975))


#Jokikalastus:

# Timo: 0.8 (0.6,0.9)
# Petri: 0.75 (0.6,0.85)
# Tapsa: moodi=0.6 95% PI= (0.3,0.8)



M3<-"
model{  
x[1]~dbeta(mu*eta, (1-mu)*eta) 
mu<-0.77
eta<-30

x[2]~dbeta(0.755*50, (1-0.755)*50)

x[3]~dbeta(0.56*13, (1-0.56)*13)


}"

cat(M3,file="prior.txt")

system.time(jm<-jags.model('prior.txt',n.adapt=100,n.chains=1))

system.time(chains1<-coda.samples(jm,variable.names=c("x"),n.iter=1000000,thin=1))

summary(chains1, quantiles=c(0.025,0.05,0.5,0.95,0.975))
par(mfrow=c(2,2))
plot(density(chains1[,"x[1]"][[1]]), main="Jokikalastus, Timo")
abline(v=c(0.6,0.8,0.9))
plot(density(chains1[,"x[2]"][[1]]), main="Jokikalastus, Petri")
abline(v=c(0.6,0.75,0.85))
plot(density(chains1[,"x[3]"][[1]]), main="Jokikalastus, Tapsa")
abline(v=c(0.3,0.6,0.8))

summary(chains1[,"x[2]"], quantiles=c(0.025,0.5,0.975))

library(runjags)

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

# var_names<-c("C_X")
# 
# run0 <- run.jags(prior,
#                  monitor= var_names,#data=data, #inits = inits,
#                  n.chains = 2, #method = 'parallel', thin=10, burnin = 1000,
#                  #modules = "mix",
#                  keep.jags.files=F,sample =1000, adapt = 1000,
#                  progress.bar=TRUE)


cat(prior,file="priorTot.txt")
system.time(jm<-jags.model('priorTot.txt',n.adapt=100,n.chains=1))

system.time(chains<-coda.samples(jm,variable.names=c(
  "C_A","R_A","C_X","R_X", "p", "Y"
  ),n.iter=1000000,thin=1))

summary(chains, quantiles=c(0.025,0.5,0.975))

par(mfrow=c(1,2))
plot(density(chains[,"C_A"][[1]]), ylim=c(0,7), main="Rannikkokalastus")
lines(density(chains[,"C_X[1]"][[1]]), lty=2, col="red")
lines(density(chains[,"C_X[2]"][[1]]), lty=3, col="green")
lines(density(chains[,"C_X[3]"][[1]]), lty=4, col=4)


plot(density(chains[,"R_A"][[1]]), ylim=c(0,8), main="Jokikalastus")
lines(density(chains[,"R_X[1]"][[1]]), lty=2, col="red")
lines(density(chains[,"R_X[2]"][[1]]), lty=3, col="green")
lines(density(chains[,"R_X[3]"][[1]]), lty=4, col=4)








