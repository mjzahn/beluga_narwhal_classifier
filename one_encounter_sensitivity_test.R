rm(list = ls())
library(banter)
library(rfPermute)
library(tidyverse)

# initialize, run banter model, and return full model object
runMyBanter <- function(event.data, ntree) {
  initBanterModel(event.data$events) %>% 
    addBanterDetector(
      event.data$detectors, ntree = ntree, sampsize = 0.5, num.cores = 7
    ) %>% 
    runBanterModel(ntree = ntree, sampsize = 0.5)
}

extractEvents <- function(event.data, event.ids) {
  events = filter(event.data$events, event.id %in% event.ids)
  list(
    events = events,
    detectors = {
      det <- lapply(event.data$detectors, function(x) {
        x.df <- filter(x, event.id %in% event.ids)
        if(nrow(x.df) == 0) return(NULL) 
        x.ev <- filter(events, event.id %in% x.df$event.id)
        sp.freq <- table(x.ev$species)
        if(length(sp.freq) < 2 | any(sp.freq == 1)) return(NULL) else x.df
      })
      det[!sapply(det, is.null)]
    }
  )
}

# load data
load("Data/sensitivity_test_data.rdata")
ev.det.data$detectors[c(1, 6)] <- NULL

enc.pairs <- do.call(
  expand.grid,
  lapply(split(encounters$encounter.id, encounters$species), unique)
)
colnames(enc.pairs) <- paste0("enc.id.", colnames(enc.pairs))

rep.result <- lapply(
  1:nrow(enc.pairs),
  function(i) {
    enc.ids <- as.numeric(enc.pairs[i, ])
    event.ids <- encounters$event.id[encounters$encounter.id %in% enc.ids]
    rep.data <- extractEvents(ev.det.data, event.ids)
    rep.banter <- runMyBanter(rep.data, ntree = 10000)
    pred.events <- setdiff(encounters$event.id, event.ids)
    pred.data <- extractEvents(ev.det.data, pred.events)
    list(model = rep.banter, event.pred = predict(rep.banter, pred.data))
  }
)

obs.mdl <- runMyBanter(ev.det.data, 10000)
obs.cm <- confusionMatrix(getBanterModel(obs.mdl))

smry <- sapply(rep.result, function(x) {
  rf <- getBanterModel(x$model)
  cm <- confusionMatrix(rf)
  pred <- x$event.pred$validation.matrix
  c(
    n.events.45 = sum(cm[1, 1:2]),
    n.events.85 = sum(cm[2, 1:2]),
    pred.cor = cor(as.vector(cm[1:2, 1:2]), as.vector(pred)),
    pct.correct = cm[3, 3]
  )
}) 

pair.smry <- arrange(cbind(enc.pairs, t(smry)), pred.cor)
pair.smry

mean(pair.smry$pred.cor)
mean(pair.smry$pred.cor > 0)
mean(pair.smry$pred.cor > 0.5)


save(
  enc.pairs, rep.result, runMyBanter, extractEvents,
  file = "one_encounter_test_ws.rdata"
)