#### Monte Carlo Sensitivity Test ####

## Marie Zahn
## 06-18-2021

## Code for Monte Carlo sensitivity analysis to examine banter model 
## classification rate variation using data from independent acoustic encounters

## load required packages
library(PAMpal)
library(banter)
library(rfPermute)

## function to run banter model -----------------------------------------------------------
# x is Acoustic Study object

runBanter <- function(x) {
  ## export banter model
  banterAll <- export_banter(x)
  ## initialize banter model
  bant.mdl <- initBanterModel(banterAll$events)
  ## run RF models for each Detector added
  bant.mdl <- addBanterDetector(
    bant.mdl, 
    data = banterAll$detectors[c(2,3,4,5)], # does not include Det 1 and 6
    ntree = 10000, 
    sampsize = 50,
    importance = FALSE # FALSE because we don't need imp scores for Monte Carlo test
  )
  
  ## run Event model
  bant.mdl <- runBanterModel(bant.mdl, ntree = 10000, sampsize = 9)
  
}

## use banter function to run Monte Carlo sensitivity test --------------------------------

## load data - narwhal and beluga echolocation parameters obtained from PAMGuard and PAMpal
load(file='Data/Narwhal-Beluga-data_043021.rdata')

## specify number of monte carlo simulations
nsims <- 1000

## set up matrix to store model output
output <- matrix(data = NA, nrow = nsims, ncol = 1)
colnames(output) <- "correct.class"

## loop that draws random subsamples from beluga and narwhal encounters and runs banter function

for (i in 1:nsims) {
  ## draw random 2 events from beluga and narwhal encounters
  
  ## runBanter()
  
  ## extract classification rate
  ## get RF data from banter model
  event.rf <- getBanterModel(bant.mdl, "event")
  confusionMat <- confusionMatrix(event.rf)
  class.rate <- confusionMat[3,3] # extract overall pct.correct
  
  ## assign result to output matrix
  output[i,] <- class.rate
}