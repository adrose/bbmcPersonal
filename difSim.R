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
total.items <- c(20) ## Keep this static - it is not built into the directory creation structure ATM!!!
dif.in.dif <- c(.3,.6,1)
floor.dif <- c(-3, -1, 1)
floor.dis <- c(.3, 1.5)
sample.size <- c(200, 1000)
sample.dist <- c(.5)
sim.total <- 100
## Now create this combination of permutations
all.folds <- expand.grid(dif.items, total.items, dif.in.dif, floor.dif, floor.dis, sample.size, sample.dist)

## Now create an out directory to store all of the sim data sets
base.dir <- "./data/difSimData/uniData/"

# ---- run-sim-uniform-dif -----------------------------------------------------------------
# Create the progress bar.
pb <- txtProgressBar(min = 1, max = dim(all.folds)[1], style=3)

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
  ceiling.discrim.ref <- floor.discrim.ref + 1.25
  floor.difficulty.foc <- floor.difficulty.ref - all.folds[i,3]
  ceiling.difficulty.foc <- ceiling.difficulty.ref - all.folds[i,3]
  sample.size.ref <- floor(all.folds[i,6] * (1-all.folds[i,7]))
  sample.size.foc <- floor(all.folds[i,6] * (all.folds[i,7]))
  # Now create the output directory to store all of the values
  out.dir <- paste(base.dir, "/sampleSize_", all.folds[i,6], "sampleProp", all.folds[i,7],"/totalItems_", 
                   total.item.val,"/propDif_", total.dif.items,"/floorDifRef_", floor.difficulty.ref, "/floorDifFoc_", 
                   floor.difficulty.foc, "/floorDisRef_", floor.discrim.ref ,sep='')
  # Create directory if it does not exist
  if(!dir.exists(out.dir)){
    dir.create(out.dir, recursive = T)
  }
  results <- foreach(W = 1:sim.total, .errorhandling = "pass", .packages = c("psych")) %dopar% {
      # Create an output file name
      out.file <- paste(out.dir,"/",W, ".RData", sep='')
      ## Now check if the file exists or not
      if(!file.exists(out.file)){
        # sample the difficulty estimates
        difficulty.no.dif <- runif(total.item.val, floor.difficulty.ref, ceiling.difficulty.ref)
        difficulty.with.dif.ref <- runif(total.dif.items, floor.difficulty.foc, max(difficulty.no.dif[1:total.dif.items]))
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
        # Now do a model without any of the DIF items
        mod.two <- irt.fa(all.vals[,-c(1:total.dif.items)])
        ## Now score the estimates
        IRT_scores <- scoreIrt(mod.one, all.vals)$theta1
        IRT_scores_nd <- scoreIrt(mod.two, all.vals[,-c(1:total.dif.items)])$theta1
        ## Now output all of these values
        output.list <- list(ref.dif = difficulty.no.dif, foc.dif = difficulty.with.dif.ref, IRTmod = mod.one, IRTmodND = mod.two, true.scores = true.theta, est.vals = IRT_scores, est.valsND = IRT_scores_nd)
        ## Now write these to disk
        save(output.list, file=out.file)
      }
      ## Now read the output
      load(out.file, verbose = T)
      print(output.list)
      
  }
  setTxtProgressBar(pb, i)
  results
}
vals.uni <- vals

# ---- create-an-anova-modeling-differences-incor -----------------------------------------------------------------
all.anova.vals <- NULL
for(i in 1:dim(all.folds)[1]){
  dif.in.dif <-  all.folds[i,3]
  floor.dif  <-  all.folds[i,4]
  floor.dis  <- all.folds[i,5]
  num.items <- all.folds[i,2]
  num.dif <- all.folds[i,1]
  samp.size <- all.folds[i,6]
  prop.dif <- all.folds[i,7]
  ## Now grab the cor values
  for(q in 1:sim.total){
    ## Get the estimated values
    if(!is.na(class(vals[[i]][[q]])[2])){
      print(paste("error", "I==", i, "q==", q))
      next
    }
    # Grab the estimated values
    est.vals <- vals[[i]][[q]]$est.vals
    est.vals.ND <- vals[[i]][[q]]$est.valsND
    ## Get the true values
    tru.vals <- vals[[i]][[q]]$true.scores
    ## Cor these values
    cor.val <- cor(est.vals, tru.vals)
    cor.val2 <- cor(est.vals.ND, tru.vals)
    cor.val3 <- cor(est.vals, est.vals.ND)
    ## NOw write this row
    out.row <- c(dif.in.dif, floor.dif, floor.dis, num.items, num.dif, samp.size, prop.dif, cor.val, cor.val2, cor.val3)
    all.anova.vals <- rbind(all.anova.vals, out.row)
  }
}
## Now train the model
all.anova.vals <- as.data.frame(all.anova.vals)
colnames(all.anova.vals) <- c("dif.in.dif",
                              "floor.dif",
                              "floor.dis",
                              "num.items",
                              "num.dif",
                              "samp.size",
                              "prop.dif",
                              "cor.val",
                              "cor.val2",
                              "cor.val3")
all.anova.vals[,1:7] <- apply(all.anova.vals[,1:7], 2, as.factor)
# Identify any variable with multiple values
length.vals <- apply(all.anova.vals, 2, function(x) length(table(x)))
## Now create a formula based on any varaible that has multiple values
outcome.var <- "cor.val"

mod <- lm(cor.val ~ (dif.in.dif+floor.dif+floor.dis+num.dif+samp.size)^3, data=all.anova.vals)
mod2 <- lm(cor.val2 ~ (dif.in.dif+floor.dif+floor.dis+num.dif+samp.size)^3, data=all.anova.vals)
mod3 <- lm(cor.val3 ~ (dif.in.dif+floor.dif+floor.dis+num.dif+samp.size)^3, data=all.anova.vals)

# ---- plot-uni-dif-results -----------------------------------------------------------------
library(visreg)
library(tidyverse)
## Most interesting interactions will be plotted here
p1 <- visreg(mod, "dif.in.dif", by="floor.dif", cond=list(num.dif="0.1"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1")
p1$layers[[4]] <- NULL
p2 <- visreg(mod, "dif.in.dif", by="floor.dif", cond=list(num.dif="0.3"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1")
p2$layers[[4]] <- NULL
p3 <- visreg(mod, "dif.in.dif", by="floor.dif", cond=list(num.dif="0.5"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1")
p3$layers[[4]] <- NULL
multiplot(p1, p2, p3, cols = 3)

p1 <- visreg(mod, "floor.dis", by="floor.dif", cond=list(samp.size="200"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1")
p1$layers[[3]] <- NULL
p2 <- visreg(mod, "floor.dis", by="floor.dif", cond=list(samp.size="1000"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1")
p2$layers[[3]] <- NULL
multiplot(p1, p2, cols = 2)


p1 <- visreg(mod, "dif.in.dif", by="samp.size", cond=list(num.dif="0.1"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1")
p1$layers[[4]] <- NULL
p2 <- visreg(mod, "dif.in.dif", by="samp.size", cond=list(num.dif="0.3"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1")
p2$layers[[4]] <- NULL
p3 <- visreg(mod, "dif.in.dif", by="samp.size", cond=list(num.dif="0.5"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1")
p3$layers[[4]] <- NULL
multiplot(p1, p2, p3, cols = 3)


p1 <- visreg(mod, "dif.in.dif", by="floor.dif", cond=list(floor.dis="1.5"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1") + ggtitle("Floor Dis = 1.5")
p1$layers[[4]] <- NULL
p2 <- visreg(mod, "dif.in.dif", by="floor.dif", cond=list(floor.dis="0.3"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1") + ggtitle("Floor Dis = 0.3")
p2$layers[[4]] <- NULL
multiplot(p1, p2, cols = 2)

  

# ---- declare-sim-nonuniform-dif -----------------------------------------------------------------
dif.items   <- c(.1, .3, .5)
total.items <- c(20)
dif.in.dif <- c(0, .3,.6, 1)
dif.in.dis <- c(0, .3, .5, .7)
floor.dif <- c(-1, 1)
floor.dis <- c(.3, 1.5)
sample.size <- c(200, 1000)
sample.dist <- c(.5)
sim.total <- 100

## Now create this combination of permutations
all.folds <- expand.grid(dif.items, total.items, dif.in.dif, floor.dif, floor.dis, sample.size, sample.dist, dif.in.dis)
## Now add a null row



## Now create an out directory to store all of the sim data sets
base.dir <- "./data/difSimData/nuniData/"


# ---- run-sim-nonuniform-dif -----------------------------------------------------------------
# Create the progress bar.
pb <- txtProgressBar(min = 1, max = dim(all.folds)[1], style=3)

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
  floor.discrim.foc <- floor.discrim.ref - all.folds[i,8]
  ceiling.discrim.foc <- ceiling.discrim.ref - all.folds[i,8]
  sample.size.ref <- floor(all.folds[i,6] * (1-all.folds[i,7]))
  sample.size.foc <- floor(all.folds[i,6] * (all.folds[i,7]))
  # Now create the output directory to store all of the values
  out.dir <- paste(base.dir, "/sampleSize_", all.folds[i,6], "sampleProp", all.folds[i,7],"/totalItems_", total.item.val,"/propDif_", 
                  total.dif.items,"/floorDifRef_", floor.difficulty.ref, "/floorDifFoc_", floor.difficulty.foc, "/floorDisRef_",floor.discrim.ref,
                  "/floorDisFoc_", floor.discrim.foc, sep='')
  # Create directory if it does not exist
  if(!dir.exists(out.dir)){
    dir.create(out.dir, recursive = T)
  }
  
  results <- foreach(W = 1:sim.total, .errorhandling = "pass", .packages = c("psych")) %dopar% {
    # Create an output file name
    out.file <- paste(out.dir,"/",W, ".RData", sep='')
    ## Now check if the file exists or not
    if(!file.exists(out.file)){
      # sample the difficulty estimates
      difficulty.no.dif <- runif(total.item.val, floor.difficulty.ref, ceiling.difficulty.ref)
      difficulty.with.dif.ref <- runif(total.dif.items, floor.difficulty.foc, max(difficulty.no.dif[1:total.dif.items]))
      difficulty.with.dif <- difficulty.no.dif
      difficulty.with.dif[1:total.dif.items] <- difficulty.with.dif.ref
      # Now sample the discrimination values
      discrim.no.dif <- runif(total.item.val, floor.discrim.ref, ceiling.discrim.ref)
      discrim.with.dif.foc <- runif(total.dif.items, floor.discrim.foc, max(discrim.no.dif[1:total.dif.items]))
      discrim.with.dif <- discrim.no.dif
      discrim.with.dif[1:total.dif.items] <- discrim.with.dif.foc
      # Now sample the true theta
      true.theta <- rnorm(all.folds[i,6])
      theta.ref <- true.theta[1:sample.size.ref]
      theta.foc <- true.theta[-c(1:sample.size.ref)]
      ## Now simulate these data
      sim.ref <- sim.irt(nvar=total.item.val, n=sample.size.ref, a = discrim.no.dif, d = difficulty.no.dif, theta = theta.ref)
      sim.foc <- sim.irt(nvar=total.item.val, n=sample.size.foc, a = discrim.with.dif, d = difficulty.with.dif, theta = theta.foc)
      ## Now combine these data and estimate the item params
      all.vals <- rbind(sim.ref$items, sim.foc$items)
      mod.one <- irt.fa(all.vals, plot=FALSE)
      # Now do a model without any of the DIF items
      mod.two <- irt.fa(all.vals[,-c(1:total.dif.items)])
      ## Now score the estimates
      IRT_scores <- scoreIrt(mod.one, all.vals)$theta1
      IRT_scores_nd <- scoreIrt(mod.two, all.vals[,-c(1:total.dif.items)])$theta1
      ## Now output all of these values
      output.list <- list(ref.dif = difficulty.no.dif, foc.dif = difficulty.with.dif.ref, IRTmod = mod.one, IRTmodND = mod.two, true.scores = true.theta, est.vals = IRT_scores, est.valsND = IRT_scores_nd)
      ## Now write these to disk
      save(output.list, file=out.file)
    }
    ## Now read the output
    load(out.file, verbose = T)
    print(output.list)
   }
  setTxtProgressBar(pb, i)
  results
}

stopCluster(cl)

vals.nuni <- vals
# ---- create-an-anova-modeling-differences-incor-nonuni -----------------------------------------------------------------
all.anova.vals <- NULL
for(i in 1:dim(all.folds)[1]){
  dif.in.dif <-  all.folds[i,3]
  floor.dif  <-  all.folds[i,4]
  floor.dis  <- all.folds[i,5]
  num.items <- all.folds[i,2]
  num.dif <- all.folds[i,1]
  samp.size <- all.folds[i,6]
  dif.in.dis <- all.folds[i,8]
  ## Now grab the cor values
  for(q in 1:sim.total){
    ## Get the estimated values
    if(!is.na(class(vals[[i]][[q]])[2])){
      print(paste("error", "I==", i, "q==", q))
      next
    }
    ## Get the estimated values
    est.vals <- vals[[i]][[q]]$est.vals
    # Grab the estimated values
    est.vals <- vals[[i]][[q]]$est.vals
    est.vals.ND <- vals[[i]][[q]]$est.valsND
    ## Get the true values
    tru.vals <- vals[[i]][[q]]$true.scores
    ## Cor these values
    cor.val <- cor(est.vals, tru.vals)
    cor.val2 <- cor(est.vals.ND, tru.vals)
    cor.val3 <- cor(est.vals, est.vals.ND)
    ## NOw write this row
    out.row <- c(dif.in.dif, floor.dif, floor.dis, num.items, num.dif, samp.size, cor.val, dif.in.dis, cor.val2, cor.val3)
    all.anova.vals <- rbind(all.anova.vals, out.row)
  }
}
## Now train the model
all.anova.vals <- as.data.frame(all.anova.vals)
all.anova.vals[,c(1:6, 8)] <- apply(all.anova.vals[,c(1:6, 8)], 2, as.factor)
colnames(all.anova.vals) <- c("dif.in.dif", "floor.dif", "floor.dis", "num.items", "prop.dif", "samp.size", "cor.val", "dif.in.dis", "cor.val2", "cor.val3")
mod1 <- lm(cor.val ~ (dif.in.dif + floor.dif + floor.dis + prop.dif + samp.size + dif.in.dis)^4, data=all.anova.vals)
mod2 <- lm(cor.val2 ~ (dif.in.dif + floor.dif + floor.dis + prop.dif + samp.size + dif.in.dis)^3, data=all.anova.vals)
mod3 <- lm(cor.val3 ~ (dif.in.dif + floor.dif + floor.dis + prop.dif + samp.size + dif.in.dis)^3, data=all.anova.vals)


# ---- plot-nonuni-dif-results -----------------------------------------------------------------
library(visreg)
## floor.dif:floor.dis:prop.dif:dif.in.dis 
p1 <- visreg(mod1, "prop.dif", by="dif.in.dis", cond=list(floor.dis=.3, floor.dif="-1"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.6, 1)) + scale_color_brewer(palette="Set1") + ggtitle("Floor Dis = .3 ; Floor Dif = -1")
p2 <- visreg(mod1, "prop.dif", by="dif.in.dis", cond=list(floor.dis=.3, floor.dif="1"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.6, 1)) + scale_color_brewer(palette="Set1") + ggtitle("Floor Dis = .3 ; Floor Dif = 1")
p3 <- visreg(mod1, "prop.dif", by="dif.in.dis", cond=list(floor.dis=1.5, floor.dif="-1"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.6, 1)) + scale_color_brewer(palette="Set1") + ggtitle("Floor Dis = 1.5 ; Floor Dif = -1")
p4 <- visreg(mod1, "prop.dif", by="dif.in.dis", cond=list(floor.dis=1.5, floor.dif="1"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.6, 1)) + scale_color_brewer(palette="Set1") + ggtitle("Floor Dis = 1.5 ; Floor Dif = 1")
p1$layers[[4]] <- NULL
p2$layers[[4]] <- NULL
p3$layers[[4]] <- NULL
p4$layers[[4]] <- NULL
multiplot(p1, p2 , p3, p4, cols=2)

# dif.items   <- c(.1, .3, .5)
# total.items <- c(20)
# dif.in.dif <- c(0, .3,.6, 1)
# dif.in.dis <- c(0, .3, .5, .7)
# floor.dif <- c(-1, 1)
# floor.dis <- c(.3, 1.5)
# sample.size <- c(200, 1000)
# sample.dist <- c(.5)
# sim.total <- 100


## floor.dis:prop.dif:dif.in.dis
dif.in.dis <- c(0, .3, .5, .7)
p1 <- visreg(mod1, "floor.dis", by="prop.dif", cond=list(dif.in.dis="0"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set1") + ggtitle("Dif in Dis = 0")
p1$layers[[3]] <- NULL
p2 <- visreg(mod1, "floor.dis", by="prop.dif", cond=list(dif.in.dis="0.3"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set1") + ggtitle("Dif in Dis = 0.3")
p2$layers[[3]] <- NULL
p3 <- visreg(mod1, "floor.dis", by="prop.dif", cond=list(dif.in.dis="0.5"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set1") + ggtitle("Dif in Dis = 0.5")
p3$layers[[3]] <- NULL
p4 <- visreg(mod1, "floor.dis", by="prop.dif", cond=list(dif.in.dis="0.7"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set1") + ggtitle("Dif in Dis = 0.7")
p4$layers[[3]] <- NULL
multiplot(p1, p2, p3, p4, cols=4)

## floor.dis:prop.dif:samp.size
p1 <- visreg(mod1, "floor.dis", by="prop.dif", cond=list(samp.size="200"), overlay=T, gg=T) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set1") + ggtitle("Samp size = 200")
p1$layers[[3]] <- NULL
p2 <- visreg(mod1, "floor.dis", by="prop.dif", cond=list(samp.size="1000"), overlay=T, gg=T) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set1") + ggtitle("Samp size = 1000")
p2$layers[[3]] <- NULL
multiplot(p1, p2, cols=2)

## floor.dif:prop.dif:dif.in.dis 
p1 <- visreg(mod1, "floor.dif", by="prop.dif", cond=list(dif.in.dis="0"), overlay=T, gg=T) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set2") + ggtitle("Dif in Dis = 0")
p1$layers[[3]] <- NULL
p2 <- visreg(mod1, "floor.dif", by="prop.dif", cond=list(dif.in.dis="0.3"), overlay=T, gg=T) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set2") + ggtitle("Dif in Dis = 0.3")
p2$layers[[3]] <- NULL
p3 <- visreg(mod1, "floor.dif", by="prop.dif", cond=list(dif.in.dis="0.5"), overlay=T, gg=T) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set2") + ggtitle("Dif in Dis = 0.5")
p3$layers[[3]] <- NULL
p4 <- visreg(mod1, "floor.dif", by="prop.dif", cond=list(dif.in.dis="0.7"), overlay=T, gg=T) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set2") + ggtitle("Dif in Dis = 0.7")
p4$layers[[3]] <- NULL
multiplot(p1, p2 ,p3, p4, cols=4)

#floor.dif:floor.dis:samp.size
p1 <- visreg(mod1, "floor.dif", by="floor.dis", cond=list(samp.size="200"), overlay=T, gg=T) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set1")
p1$layers[[3]] <- NULL
p2 <- visreg(mod1, "floor.dif", by="floor.dis", cond=list(samp.size="1000"), overlay=T, gg=T) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set1")
p2$layers[[3]] <- NULL
multiplot(p1, p2, cols=2)

## dif.in.dif:floor.dif:floor.dis
p1 <- visreg(mod1, "dif.in.dif", by="floor.dis", cond=list(floor.dif="-1"), overlay=T, gg=T) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set2") + ggtitle("Floor Dif = -1")
p1$layers[[5]] <- NULL
p2 <- visreg(mod1, "dif.in.dif", by="floor.dis", cond=list(floor.dif="1"), overlay=T, gg=T) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set2") + ggtitle("Floor Dif = 1")
p2$layers[[5]] <- NULL
multiplot(p1, p2, cols=2)


# ---- plot-nonuni-dif-results-wd-vs-nd -----------------------------------------------------------------
p1 <- visreg(mod3, "floor.dif", by="prop.dif", cond=list(samp.size="1000", floor.dis=.3), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set1") + ggtitle("n=1000 - floor.dis=.3") + coord_cartesian(ylim=c(.85, 1)) 
p1$layers[[4]] <- NULL
p2 <- visreg(mod3, "floor.dif", by="prop.dif", cond=list(samp.size="1000", floor.dis=1.5), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set1") + ggtitle("n=1000 - floor.dis=1.3") + coord_cartesian(ylim=c(.85, 1))
p2$layers[[4]] <- NULL
p4 <- visreg(mod3, "floor.dif", by="prop.dif", cond=list(samp.size="200", floor.dis=.3), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set1") + ggtitle("n=200 - floor.dis=.3") + coord_cartesian(ylim=c(.85, 1))
p4$layers[[4]] <- NULL
p5 <- visreg(mod3, "floor.dif", by="prop.dif", cond=list(samp.size="200", floor.dis=1.5), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.7, 1)) + scale_color_brewer(palette="Set1") + ggtitle("n=200 - floor.dis=1.3") + coord_cartesian(ylim=c(.85, 1))
p5$layers[[4]] <- NULL

multiplot(p1, p2, p4, p5, cols=2)
