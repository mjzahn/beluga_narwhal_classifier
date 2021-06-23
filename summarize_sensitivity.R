rm(list = ls())
library(rfPermute)
library(tidyverse)
library(banter)

load("sensitivity_test_ws.rdata")

obs.mdl <- runMyBanter(ev.det.data$events, ev.det.data$detectors, 10000)
obs.cm <- confusionMatrix(getBanterModel(obs.mdl))

cm <- sapply(
  rep.result, 
  function(x) confusionMatrix(getBanterModel(x)), 
  simplify = "array"
)

plot(density(cm["Overall", "pct.correct", ]))
abline(v = obs.cm["Overall", "pct.correct"])

obs.cm
median(cm["Overall", "pct.correct", ])
quantile(cm["Overall", "pct.correct", ], c(0.025, 0.975))
mean(cm["Overall", "pct.correct", ] > obs.cm["Overall", "pct.correct"]) # average of how many monte carlo sims are greater than observed
