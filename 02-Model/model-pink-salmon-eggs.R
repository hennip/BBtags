library(runjags)
library(readxl)
require(rjags)
require(sp)
require(openxlsx)
library(readxl)
library(tidyverse)
library(lubridate)



dat=read_xlsx("C:/Users/03080932/OneDrive - Valtion/Projects/BBTags/dat/Kyttyrämunadata.xlsx",
              range = "A2:E34")

df1<-dat%>%filter(Week==36, Site==1)%>%select(-Week, -Site)
df2<-dat%>%filter(Week==36, Site==2)%>%select(-Week, -Site)
df3<-dat%>%filter(Week==36, Site==3)%>%select(-Week, -Site)

x36<-array(NA, dim=c(3,3,3))
x36[,,1]<-as.matrix(df1)
x36[,,2]<-as.matrix(df2)
x36[,,3]<-as.matrix(df3)

N36<-array(NA, dim=c(3,3))
N36[,1]=apply(x36[,,1], 1, sum)
N36[,2]=apply(x36[,,2], 1, sum)
N36[,3]=apply(x36[,,3], 1, sum)
# ==============================

df2<-dat%>%filter(Week==37, Site==2)%>%select(-Week, -Site)
df3<-dat%>%filter(Week==37, Site==3)%>%select(-Week, -Site)

x37<-array(NA, dim=c(3,3,2))
#x37[,,1]<-as.matrix(df1)
x37[,,1]<-as.matrix(df2)
x37[,,2]<-as.matrix(df3)

N37<-array(NA, dim=c(3,2))
N37[,1]=apply(x37[,,1], 1, sum)
N37[,2]=apply(x37[,,2], 1, sum)

# ==============================
df1<-dat%>%filter(Week==39, Site==1)%>%select(-Week, -Site)
x39<-array(NA, dim=c(3,3))
x39[,]<-as.matrix(df1)

N39<-array(NA, dim=c(3,1))
N39=apply(x39[,], 1, sum)

# ==============================
df1<-dat%>%filter(Week==40, Site==1)%>%select(-Week, -Site)
df2<-dat%>%filter(Week==40, Site==2)%>%select(-Week, -Site)
df3<-dat%>%filter(Week==40, Site==3)%>%select(-Week, -Site)

x40<-array(NA, dim=c(3,3,3))
x40[,,1]<-as.matrix(df1)
x40[,,2]<-as.matrix(df2)
x40[,,3]<-as.matrix(df3)

N40<-array(NA, dim=c(3,3))
N40[,1]=apply(x40[,,1], 1, sum)
N40[,2]=apply(x40[,,2], 1, sum)
N40[,3]=apply(x40[,,3], 1, sum)


# ==============================
df1<-dat%>%filter(Week==43, Site==1)%>%select(-Week, -Site)
df2<-dat%>%filter(Week==43, Site==2)%>%select(-Week, -Site)
df3<-dat%>%filter(Week==43, Site==3)%>%select(-Week, -Site)

x43<-array(NA, dim=c(3,3))
x43[,1]<-as.matrix(df1)
x43[,2]<-as.matrix(df2)
x43[,3]<-as.matrix(df3)

N43<-c()
N43[1]=sum(x43[,1])
N43[2]=sum(x43[,2])
N43[3]=sum(x43[,3])

# ==============================
df2<-dat%>%filter(Week==45, Site==2)%>%select(-Week, -Site)
x45<-array(NA, dim=c(2,3))
x45[,]<-as.matrix(df2)

N45<-c()
N45[1]=sum(x45[1,])
N45[2]=sum(x45[2,])



M1<-"model{

# Week 36
for(s in 1:3){ #site
  for(i in 1:3){ #toisto
    x36[i,1:3,s] ~ dmulti(p36[1:3,s],N36[i,s])
  }
  p36[1:3,s]~ddirich(alpha36[1:3,s])
}

# Week 37
for(s in 1:2){ #site
  for(i in 1:3){ #toisto
    x37[i,1:3,s] ~ dmulti(p37[1:3,s],N37[i,s])
  }
  p37[1:3,s]~ddirich(alpha37[1:3,s])
}

# Week 39
for(i in 1:3){ #toisto
  x39[i,1:3] ~ dmulti(p39[1:3],N39[i])
}
p39[1:3]~ddirich(alpha39[1:3])

# Week 40
for(s in 1:3){ #site
  for(i in 1:3){ #toisto
    x40[i,1:3,s] ~ dmulti(p40[1:3,s],N40[i,s])
  }
  p40[1:3,s]~ddirich(alpha40[1:3,s])
}

# Week 43
for(s in 1:3){ #site
 # for(i in 1:3){ #toisto
    x43[1:3,s] ~ dmulti(p43[1:3,s],N43[s])
 # }
  p43[1:3,s]~ddirich(alpha43[1:3,s])
}

# Week 45
for(i in 1:2){ #toisto
  x45[i,1:3] ~ dmulti(p45[1:3],N45[i])
}
p45[1:3]~ddirich(alpha45[1:3])



}"


data=list(
  alpha36=array(1, dim=c(3,3)), x36=x36, N36=N36,
  alpha37=array(1, dim=c(3,2)), x37=x37, N37=N37,
  alpha39=c(1,1,1), x39=array(c(0,0,0,20,20,20,1,1,0), dim=c(3,3)), N39=c(21,21,20),
  alpha40=array(1, dim=c(3,3)), x40=x40, N40=N40,
  alpha43=array(1, dim=c(3,3)), x43=x43, N43=N43,
  alpha45=c(1,1,1), x45=array(c(0,0,2,0,30,20), dim=c(2,3)), N45=c(32,20)
)  
data


var_names<- c(
  "p36", "p37", "p39", "p40", "p43", "p45")

run1 <- autorun.jags(M1,
                     monitor= var_names,data=data, #inits = inits,
                     n.chains = 2, method = 'parallel', thin=10, 
                     modules = "mix",progress.bar=TRUE)


run<-run1

sum_run<-round(summary(run, quantiles=c(0.05,0.5,0.95)),3);sum_run

chains<-as.mcmc(run)
summary(chains, quantiles=c(0.05,0.5,0.95))

# x[i,1:3,s]~dmulti (p[1:3,s], N[i,s])
# p[1:3,s]~ddirich(alpha[1:3,s])
# p[1:3,s]~ddirich(1,1,1)

# alpha[1:3,s]=c(1,1,1)