## I believe this script will be used to train some markov models.. looking into how individuals transition through the various HH versus LL on the
## eeg ID task (Amplitude & latency will be used; median splits will determine the group assignment for these two)


## Load libraries
source("~/adroseHelperScripts/R/afgrHelpFunc.R")
library(markovchain)

## Create functions
## This function will claculate the transition probabilites given a matrix
trans.matrix <- function(X, prob=T)
{
  tt <- table( c(X[,-ncol(X)]), c(X[,-1]) )
  if(prob) tt <- tt / rowSums(tt)
  tt
}

## load data
load("./fname.gz")
ds_eegPred <- out.data[[9]]

## Identify the states
mod.plot <- ds_eegPred[complete.cases(ds_eegPred$P2crymaxIDLatpost),]
mod.plot <- mod.plot[!duplicated(mod.plot$record_id),]
mod.plot$latTertilesCry <- cut(mod.plot$P2crymaxIDLatpost,     breaks = quantile(mod.plot$P2crymaxIDLatpost, probs = c(0,.5, 1)), labels = c("Q","S"), include.lowest = T)
mod.plot$ampTertilesCry <- cut(mod.plot$P2crymaxIDpost,        breaks = quantile(mod.plot$P2crymaxIDpost   , probs = c(0, .5,1)), labels = c("L","H"), include.lowest = T)
mod.plot$latTertilesNeu <- cut(mod.plot$P2neutralmaxIDLatpost, breaks = quantile(mod.plot$P2neutralmaxIDLatpost, probs = c(0, .5,1)), labels = c("Q","S"), include.lowest = T)
mod.plot$ampTertilesNeu <- cut(mod.plot$P2neutralmaxIDpost,    breaks = quantile(mod.plot$P2neutralmaxIDpost   , probs = c(0, .5,1)), labels = c("L","H"), include.lowest = T)
mod.plot$latTertilesHap <- cut(mod.plot$P2happymaxIDLatpost,   breaks = quantile(mod.plot$P2happymaxIDLatpost, probs = c(0, .5,1)), labels = c("Q","S"), include.lowest = T)
mod.plot$ampTertilesHap <- cut(mod.plot$P2happymaxIDpost,      breaks = quantile(mod.plot$P2happymaxIDpost   , probs = c(0, .5,1)), labels = c("L","H"), include.lowest = T)
mod.plot$latTertilesUnh <- cut(mod.plot$P2unhappymaxIDLatpost, breaks = quantile(mod.plot$P2unhappymaxIDLatpost, probs = c(0, .5,1)), labels = c("Q","S"), include.lowest = T)
mod.plot$ampTertilesUnh <- cut(mod.plot$P2unhappymaxIDpost,    breaks = quantile(mod.plot$P2unhappymaxIDpost, probs = c(0, .5,1)), labels = c("L","H"), include.lowest = T)
## Now make the combination of these groups
mod.plot$plot.vals <- paste(mod.plot$latTertilesCry, mod.plot$ampTertilesCry)
mod.plot$plot.valsCry <- paste(mod.plot$latTertilesCry, mod.plot$ampTertilesCry)
mod.plot$plot.valsHap <- paste(mod.plot$latTertilesHap, mod.plot$ampTertilesHap)
mod.plot$plot.valsNeu <- paste(mod.plot$latTertilesNeu, mod.plot$ampTertilesNeu)
mod.plot$plot.valsUnh <- paste(mod.plot$latTertilesUnh, mod.plot$ampTertilesUnh)

## Now create the first transition matrix going from Crying --> Unhappy
transition.probs <- trans.matrix(as.matrix(cbind(mod.plot$plot.valsCry, mod.plot$plot.valsUnh, mod.plot$plot.valsNeu, mod.plot$plot.valsHap)))

## Now create the MC object
brainStates <- c("Q H", "Q L", "S H", "S L")
bsMat <- new("markovchain", states=brainStates, transitionMatrix=matrix(transition.probs, byrow = F, nrow=4), name="cryToUnh")
## Now visualize this model
plot(bsMat)
