rm(list = ls())

library(nimble)
library(sf)

#dir <- "C:/Users/mehnazjahid/Desktop/UVic/CRT/CRTmodelPaper"
dir <- "./"
setwd(dir)

M<- 200
J_1<- 50 # total number of  dna traps
J_2<- 50 # total number of camera traps

mask <- read_sf(dsn = file.path(dir,"2012_mask_utm11n.shp"))
traps_hair<- read_sf(dsn = file.path(dir,"dna_traps_utm11n.shp"))
traps_ct<- read_sf(dsn = file.path(dir,"ct_traps_utm11n.shp"))

y_1<- read.csv("capthist_tourani.csv")
ydot_2<- read.csv("CTcapthist_tourani.csv", header = FALSE)

Model_3 <- nimbleCode({
  ##-----------------------------------------------------------------
  ## INDIVIDUAL INCLUSION
  ## AC LOCATIONS
  psi ~ dunif(0, 1)
  
  for (i in 1:M) {
    z[i] ~ dbern(psi) ## equation (2)
  }
  
  N <- sum(z[1:M]) ## equation (3)
  sigma ~ dunif(0, 100)
  p0_1 ~ dunif(0, 1) # for survey type 1 (HT in our case)
  p0_2 ~ dunif(0, 1) # for survey type 2 (CT in our case)
  ## IDENTIFIED DETECTIONS, SURVEY TYPE 1
  for (i in 1:M) {
    d_squared_1[i, 1:J_1] <- (mask[i, 2] - traps_hair[1:J_1,1])^2 +
      (sxy[i, 3] - traps_hair[1:J_1,2])^2
    p_1[i, 1:J_1] <- p0_1 * exp(-d_squared_1[i,1:J_1]/(2*sigma*sigma))
## equation (4)
    ##p_1 in model 1 and 3 is equivalent to p_1*alpha in model 2 and 4
    for(j in 1:J_1){
      y_1[i, j] ~ dbern(p_1[i, j] * z[i]) ## equation (6)
    }
    ## y_1[i, 1:J_1] ~ dbern_vector(p_1[i, 1:J_1], z[i]) ## equation (6)
  }
  ## UNIDENTIFIED DETECTIONS, SURVEY TYPE 2
  for (i in 1:M) {
    d_squared_2[i, 1:J_2] <- (mask[i, 2] - traps_ct[1:J_2,1])^2 +
      (mask[i, 3] - traps_ct[1:J_2,2])^2
    p_2[i, 1:J_2] <- p0_2 * exp(-d_squared_2[i,1:J_2]/(2*sigma*sigma)) *z[i] ## equation (4)
  }
  for (j in 1:J_2) {
    pdot_2[j] <- 1 - prod(1-p_2[1:M, j]) ## equation (5)
    ydot_2[j] ~ dbern(pdot_2[j]) ## equation (7)
  }
  ## ydot_2[1:J_2] ~ dbern_vector(pdot_2[1:J_2], 1) ## equation (7)
})

data <- list(mask = cbind(mask$id,
                          mask$x,
                          mask$y),
             traps_hair = cbind(traps_hair$x,
                                traps_hair$y),
             traps_ct =  cbind(traps_ct$x,
                               traps_ct$y),
             y_1 = y_1,
             ydot_2 = as.numeric(ydot_2$V1))

constants <- list(M = nrow(mask),
                  J_1 = J_1,
                  J_2 = J_2)

inits <- list(z = rep(z,constants$M),
              sigma = 1,
              p0_1 = .5,
              p0_2 = .5,
              psi = .5)

model <- nimbleModel(Model_3, 
                     data = data,
                     constants = constants,
                     inits = inits)

sample_model3<-nimbleMCMC(model,
                          data=data,
                          niter=100, #10000,
                          nburnin=10, #1000,
                          thin= 1, #20,
                          nchains=2,
                          samplesAsCodaMCMC = T)                                  
summary(sample_model3)

mcmcplot(sample_model3, dir = "C:/Users/mehnazjahid/Desktop/trash", filename = "tourani_mcmcplot")
