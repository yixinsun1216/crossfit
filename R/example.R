library(tidyverse)
library(Formula)
library(grf)
library(modelr)
library(furrr)
library(splines)

library("MASS")
data(Boston)

rfseed <- 19310821
set.seed(rfseed)

crime_dml <-
  "crim ~ zn + indus + nox" %>%
  as.formula() %>%
  dml(Boston, psi_plr, psi_plr_grad, psi_plr_op, n = 101)

