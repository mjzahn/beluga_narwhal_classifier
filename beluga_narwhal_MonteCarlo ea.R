rm(list = ls())
library(PAMpal)
library(banter)
library(rfPermute)
library(tidyverse)

# initialize, run banter model, and return full model object
runMyBanter <- function(events, detectors, ntree) {
  initBanterModel(events) %>% 
    addBanterDetector(detectors, ntree = ntree) %>% 
    runBanterModel(ntree = ntree)
}

# return a random n.event event.ids sampled from a random n.encounters for each species
sampleEvents <- function(n.encounters, n.events, df) {
  unname(unlist(lapply(split(df, df$species), function(spp.df) {
    enc.ids <- unique(spp.df$encounter.id)
    ran.encounters <- sample(enc.ids, min(length(enc.ids), n.encounters))
    unlist(lapply(ran.encounters, function(en) {
      en.rows <- which(spp.df$encounter.id == en)
      spp.df$event.id[sample(en.rows, min(length(en.rows), n.events))]
    }))
  })))
}

# load data
load("Data/sensitivity_test_data.rdata")
ev.det.data$detectors[c(1, 6)] <- NULL

# run sensitivity replicates
rep.result <- replicate(
  1000, 
  {
    # randomly sample encounters and events for this replicate
    rep.events <- sampleEvents(3, 5, encounters)
    # extract event data
    events <- filter(ev.det.data$events, event.id %in% rep.events)
    # extract detector data
    detectors <- lapply(ev.det.data$detectors, function(x) {
      x.df <- filter(x, event.id %in% rep.events)
      if(nrow(x.df) == 0) return(NULL) 
      x.ev <- filter(events, event.id %in% x.df$event.id)
      sp.freq <- table(x.ev$species)
      if(length(sp.freq) < 2 | any(sp.freq == 1)) return(NULL) else x.df
    })
    detectors <- detectors[!sapply(detectors, is.null)]
    # run banter model for this replicate
    runMyBanter(events, detectors, ntree = 10000)
  },
  simplify = FALSE
)

save(ev.det.data, rep.result, runMyBanter, file = "sensitivity_test_ws.rdata")
