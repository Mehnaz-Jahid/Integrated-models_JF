# Blanc et al. 

# Supplementary Information for Blanc et al.
# Improving abundance estimation by combining capture-recapture and occupancy data: 
# example with a large carnivore
# This is the model Mh using Schofield and Barker 2014 method see Appendix 1 of their paper for details

# Load packages
library(coda)
library(rjags)
#library(mcmcplots)

# Setup directory

#--------------------- capture-recapture data -----------------------------#
# Capture-recapture data structure : rows = individuals; columns = capture occasions
mydata <- read.delim('BlancEtAl/HTcapthist.txt', sep = "\t")

extra = 250 # define large number of extra individual capture histories
n = nrow(mydata) # number of observed individuals
M = extra + n
xn = rowSums(mydata[,2:9])
x = c(xn,rep(0,extra))
k = ncol(mydata)
zerouse <- 0
#-------------------- specify model ------------------------------#
sink("CRandPO_HET.txt")
cat("
model{
    
    
    # Capture-recapture (CR) likelihood
 for(i in 1:M){
    x[i] ~ dbin(pee[i],k)
    pee[i] <- pcr[i]*step(N-i)  # pcr = individual detection probability pi from CR model
    logit(pcr[i]) <- mumup + xi*etap[i] # integration of an individual random effect ??i
    etap[i] ~ dnorm(0,sigmeta)
 }
    
    N ~ dpois(lamn)T(0,M)
  #  zerouse ~ dpois(comb)           						#P(zerouse=0)=exp(-comb)
  #  comb <- logfact(M)-logfact(N)+logfact((N-n)*step(N-n)). #ln(M!)-ln(N!)+ln[(N-n)! I(N>n)]        
    
    
   
     # priors for detection probability in CR likelihood
    mumup ~ dlogis(0,1)
    sigmeta ~ dgamma(1.5,37.5)
    xi ~ dnorm(0,1) 
    lamn ~ dgamma(10,1/3.7)
} 
    ",fill=TRUE)
sink()

#-------------------- Jags stuff -----------------------#
# List of data
mydatax <- list (x=x,M=M,n=n,k=k,zerouse=zerouse,y=y,nsites=nsites,nsurvs=nsurvs)

# Parameters monitored
parameters = c ('N', 'ppo', 'mumup0', 'psi0', 'mupsi', 'mean.p0', 'sdeps', 'taups', 'betapsi', 'etap')


# Initial values
init1 <- list(mupsi=runif(1),truocc=as.numeric(apply(y,1,sum)>0),xi = rnorm(1,sd=2), N = sample(n:M,1),mumup = rnorm(1), taup = rlnorm(1))
init2 <- list(mupsi=runif(1),truocc=as.numeric(apply(y,1,sum)>0),xi = rnorm(1,sd=2), N = sample(n:M,1),mumup = rnorm(1), taup = rlnorm(1))
init3 <- list(mupsi=runif(1),truocc=as.numeric(apply(y,1,sum)>0),xi = rnorm(1,sd=2), N = sample(n:M,1),mumup = rnorm(1), taup = rlnorm(1))
inits <- list(init1,init2,init3)

#-------------------- Call jags from R -----------------------#
jmodel <- jags.model("CRandPO_HET.txt", mydatax, inits, n.chains = 3,n.adapt = 2500)
jsample <- coda.samples(jmodel, parameters, n.iter=15000, thin = 1)
mcmcplot(jsample, dir = "C:/Users/mehnazjahid/Desktop/UVic/CRT/CRTmodelPaper/trash", filename = "jsample_mcmcplot")
summary(jsample)

# capture recapture
install.packages("Rcapture")
library(Rcapture)
cr<- mydata[,2:9]
crmod<- closedp(cr)  
