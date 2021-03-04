## Thislibrary will be used to run a simulation study fir the impact that DIF has on factor estimates
## In order to do this I am going to have to simulate two different populations that have differences in:
## specific items difficulty parameters
## discrimination params
## the difficult thing will to be simulate 1 covariate that continulsy relates with either difficulty & discrimination
## and pull from a single population 

# ---- load-sources -----------------------------------------------------------------
library(psych)
library(parallel)
library(doParallel)
library(foreach)
library(ggplot2)
library(reshape2)
library(caret)
library(mirt)
library(utils)
source("~/adroseHelperScripts/R/afgrHelpFunc.R")

# ---- set-parallel-env -----------------------------------------------------------------
cl <- makeCluster(3)
registerDoParallel(cl)


# ---- declare-sim-uniform-dif -----------------------------------------------------------------
## First identify the number of items with DIF
dif.items   <- c(.1, .3, .5)
total.items <- c(20)
dif.in.dif <- c(.3,.6,1)
floor.dif <- c(-1, 1)
floor.dis <- c(.3, 1.5)
sample.size <- c(1000)
sample.dist <- c(.25,.5)
sim.total <- 100
## Now create this combination of permutations
all.folds <- expand.grid(dif.items, total.items, dif.in.dif, floor.dif, floor.dis, sample.size, sample.dist)

# ---- run-sim-uniform-dif -----------------------------------------------------------------
# Create the progress bar.
pb <- txtProgressBar(min = 1, max = iterations, style=3)

## Now loop through each of these and return a vector of difficulty parameters
## The output of this process will be a list with the referent group and the focal group difficulty params
## Each total permutation (row in all.folds) will have 500 samples drawn; so the output list will have dim(all.folds)[1] * 500 observations
all.dif.values <- list()
vals <- foreach(i=1:dim(all.folds)[1], .export = "results") %do%{
  ## Create some variables which will be static across all samples of this permutation
  total.item.val <- all.folds[i,2]
  total.dif.items <- floor(total.item.val * all.folds[i,1])
  floor.difficulty.ref <- all.folds[i,4]
  floor.discrim.ref <- all.folds[i,5]
  ceiling.difficulty.ref <- floor.difficulty.ref + 2
  ceiling.discrim.ref <- floor.discrim.ref + 1.5
  floor.difficulty.foc <- floor.difficulty.ref - all.folds[i,3]
  ceiling.difficulty.foc <- ceiling.difficulty.ref - all.folds[i,3]
  sample.size.ref <- floor(all.folds[i,6] * (1-all.folds[i,7]))
  sample.size.foc <- floor(all.folds[i,6] * (all.folds[i,7]))
  results <- foreach(W = 1:sim.total, .errorhandling = "pass", .packages = c("psych")) %dopar% {
      # sample the difficulty estimates
      difficulty.no.dif <- runif(total.item.val, floor.difficulty.ref, ceiling.difficulty.ref)
      difficulty.with.dif.ref <- runif(total.dif.items, floor.difficulty.foc, ceiling.difficulty.ref)
      difficulty.with.dif <- difficulty.no.dif
      difficulty.with.dif[1:total.dif.items] <- difficulty.with.dif.ref
      # Now sample the discrimination values
      discrim.no.dif <- runif(total.item.val, floor.discrim.ref, ceiling.discrim.ref)
      # Now sample the true theta
      true.theta <- rnorm(all.folds[i,6])
      theta.ref <- true.theta[1:sample.size.ref]
      theta.foc <- true.theta[-c(1:sample.size.ref)]
      ## Now simulate these data
      sim.ref <- sim.irt(nvar=total.item.val, n=sample.size.ref, a = discrim.no.dif, d = difficulty.no.dif, theta = theta.ref)
      sim.foc <- sim.irt(nvar=total.item.val, n=sample.size.foc, a = discrim.no.dif, d = difficulty.with.dif, theta = theta.foc)
      ## Now combine these data and estimate the item params
      all.vals <- rbind(sim.ref$items, sim.foc$items)
      mod.one <- irt.fa(all.vals, plot=FALSE)
      ## Now score the estimates
      IRT_scores <- scoreIrt(mod.one, all.vals)$theta1
      ## Now output all of these values
      output.list <- list(ref.dif = difficulty.no.dif, foc.dif = difficulty.with.dif.ref, IRTmod = mod.one, true.scores = true.theta, est.vals = IRT_scores)
  }
  setTxtProgressBar(pb, i)
  results
}
stopCluster(cl)


# ---- create-an-anova-modeling-differences-incor -----------------------------------------------------------------
all.anova.vals <- NULL
for(i in 1:dim(all.folds)[1]){
  dif.in.dif <-  all.folds[i,3]
  floor.dif  <-  all.folds[i,4]
  floor.dis  <- all.folds[i,5]
  num.items <- all.folds[i,2]
  num.dif <- all.folds[i,1]
  samp.size <- all.folds[i,6]
  ## Now grab the cor values
  for(q in 1:sim.total){
    ## Get the estimated values
    est.vals <- vals[[i]][[q]]$est.vals
    ## Get the true values
    tru.vals <- vals[[i]][[q]]$true.scores
    ## Cor these values
    cor.val <- cor(est.vals, tru.vals)
    ## NOw write this row
    out.row <- c(dif.in.dif, floor.dif, floor.dis, num.items, num.dif, samp.size, cor.val)
    all.anova.vals <- rbind(all.anova.vals, out.row)
  }
}
## Now train the model
all.anova.vals <- as.data.frame(all.anova.vals)
all.anova.vals[,1:6] <- apply(all.anova.vals[,1:6], 2, as.factor)
mod <- lm(V7 ~ (V1 + V2 + V3 + V5 )^4, data=all.anova.vals)

# ---- plot-uni-dif-results -----------------------------------------------------------------
library(visreg)
## Most interesting interactions will be plotted here
p1 <- visreg(mod, "V2", by="V1", cond=list(V5="0.1"), overlay=T, gg=TRUE)
p2 <- visreg(mod, "V2", by="V1", cond=list(V5="0.3"), overlay=T, gg=TRUE)
p3 <- visreg(mod, "V2", by="V1", cond=list(V5="0.5"), overlay=T, gg=TRUE)
multiplot(p1, p2, p3, cols = 3)

p1 <- visreg(mod, "V3", by="V2", cond=list(V5="0.1"), overlay=T, gg=TRUE)
p2 <- visreg(mod, "V3", by="V2", cond=list(V5="0.3"), overlay=T, gg=TRUE)
p3 <- visreg(mod, "V3", by="V2", cond=list(V5="0.5"), overlay=T, gg=TRUE)
multiplot(p1, p2, p3, cols = 3)

# ---- declare-sim-nonuniform-dif -----------------------------------------------------------------
dif.items   <- c(.1, .3, .5)
total.items <- c(20)
dif.in.dif <- c(.3,.6,1)
dif.in.dis <- c(.1, .3, .5)
floor.dif <- c(-1, 1)
floor.dis <- c(.3, 1.5)
sample.size <- c(1000)
sample.dist <- c(.25,.5)
sim.total <- 100