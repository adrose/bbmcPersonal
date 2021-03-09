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
cl <- makeCluster(4)
registerDoParallel(cl)


# ---- declare-sim-uniform-dif -----------------------------------------------------------------
## First identify the number of items with DIF
dif.items   <- c(.1, .3, .5)
total.items <- c(20) ## Keep this static - it is not built into the directory creation structure ATM!!!
dif.in.dif <- c(.3,.6,1)
floor.dif <- c(-1, 1)
floor.dis <- c(.3, 1.5)
sample.size <- c(200, 1000)
sample.dist <- c(.25,.5)
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
  out.dir <- paste(base.dir, "/sampleSize_", all.folds[i,6], "sampleProp", all.folds[i,7],"/totalItems_", total.item.val,"/propDif_", total.dif.items,"/floorDifRef_", floor.difficulty.ref, "/floorDifFoc_", floor.difficulty.foc, sep='')
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
    est.vals <- vals[[i]][[q]]$est.vals
    ## Get the true values
    tru.vals <- vals[[i]][[q]]$true.scores
    ## Cor these values
    cor.val <- cor(est.vals, tru.vals)
    ## NOw write this row
    out.row <- c(dif.in.dif, floor.dif, floor.dis, num.items, num.dif, samp.size, prop.dif, cor.val)
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
                              "cor.val")
all.anova.vals[,1:7] <- apply(all.anova.vals[,1:7], 2, as.factor)
# Identify any variable with multiple values
length.vals <- apply(all.anova.vals, 2, function(x) length(table(x)))
## Now create a formula based on any varaible that has multiple values
outcome.var <- "cor.val"

mod <- lm(cor.val ~ (dif.in.dif+floor.dif+floor.dis+num.dif+samp.size+prop.dif)^3, data=all.anova.vals)

# ---- plot-uni-dif-results -----------------------------------------------------------------
library(visreg)
## Most interesting interactions will be plotted here
p1 <- visreg(mod, "dif.in.dif", by="floor.dif", cond=list(num.dif="0.1"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1")
p2 <- visreg(mod, "dif.in.dif", by="floor.dif", cond=list(num.dif="0.3"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1")
p3 <- visreg(mod, "dif.in.dif", by="floor.dif", cond=list(num.dif="0.5"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1")
multiplot(p1, p2, p3, cols = 3)

p1 <- visreg(mod, "dif.in.dif", by="samp.size", cond=list(num.dif="0.1"), overlay=F, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1")
p2 <- visreg(mod, "dif.in.dif", by="samp.size", cond=list(num.dif="0.3"), overlay=F, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1")
p3 <- visreg(mod, "dif.in.dif", by="samp.size", cond=list(num.dif="0.5"), overlay=F, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, .95)) + scale_color_brewer(palette="Set1")
multiplot(p1, p2, p3, cols = 3)

## Now visualize all contrasts
all.contrasts <- summarySE(all.anova.vals, measurevar = "cor.val", groupvars = c("dif.in.dif", "floor.dif", "floor.dis", "num.dif", "samp.size", "prop.dif"), na.rm=T)

all.contrasts[which(all.contrasts$floor.dis==.3),] %>% ggplot(., aes(x=dif.in.dif, y=cor.val, fill=samp.size)) +
  geom_bar(stat='identity', position = 'dodge') +
  geom_errorbar(aes(ymin=cor.val-se, ymax=cor.val+se)) +
  facet_grid(floor.dif~num.dif)
  

# ---- declare-sim-nonuniform-dif -----------------------------------------------------------------
dif.items   <- c(.1, .3, .5)
total.items <- c(20)
dif.in.dif <- c(.3,.6, 1)
dif.in.dis <- c(.1, .3, .5)
floor.dif <- c(-1, 1)
floor.dis <- c(.3, 1.5)
sample.size <- c(200, 1000)
sample.dist <- c(.25,.5)
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
      difficulty.with.dif.ref <- runif(total.dif.items, floor.difficulty.foc, ceiling.difficulty.ref)
      difficulty.with.dif <- difficulty.no.dif
      difficulty.with.dif[1:total.dif.items] <- difficulty.with.dif.ref
      # Now sample the discrimination values
      discrim.no.dif <- runif(total.item.val, floor.discrim.ref, ceiling.discrim.ref)
      discrim.with.dif.foc <- runif(total.dif.items, floor.discrim.foc, ceiling.discrim.foc)
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
      ## Now score the estimates
      IRT_scores <- scoreIrt(mod.one, all.vals)$theta1
      ## Now output all of these values
      output.list <- list(ref.dif = difficulty.no.dif, foc.dif = difficulty.with.dif.ref, IRTmod = mod.one, true.scores = true.theta, est.vals = IRT_scores)
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
    ## Get the true values
    tru.vals <- vals[[i]][[q]]$true.scores
    ## Cor these values
    cor.val <- cor(est.vals, tru.vals)
    ## NOw write this row
    out.row <- c(dif.in.dif, floor.dif, floor.dis, num.items, num.dif, samp.size, cor.val, dif.in.dis)
    all.anova.vals <- rbind(all.anova.vals, out.row)
  }
}
## Now train the model
all.anova.vals <- as.data.frame(all.anova.vals)
all.anova.vals[,c(1:6, 8)] <- apply(all.anova.vals[,c(1:6, 8)], 2, as.factor)
mod <- lm(V7 ~ (V1 + V2 + V3 + V5 + V8)^4, data=all.anova.vals)


# ---- plot-nonuni-dif-results -----------------------------------------------------------------
library(visreg)
## Most interesting interactions will be plotted here
p1 <- visreg(mod, "V3", by="V8", cond=list(V5="0.1"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, 1)) + scale_color_brewer(palette="Set1")
p2 <- visreg(mod, "V3", by="V8", cond=list(V5="0.3"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, 1)) + scale_color_brewer(palette="Set1")
p3 <- visreg(mod, "V3", by="V8", cond=list(V5="0.5"), overlay=T, gg=TRUE) + theme_minimal() + coord_cartesian(ylim=c(.75, 1)) + scale_color_brewer(palette="Set1")
multiplot(p1, p2, p3, cols = 3)

p1 <- visreg(mod, "V3", by="V2", cond=list(V5="0.1"), overlay=T, gg=TRUE) + ggtitle("Prop DIF = .1") + theme(legend.position = "bottom") + scale_color_brewer(palette="PRGn")
p2 <- visreg(mod, "V3", by="V2", cond=list(V5="0.3"), overlay=T, gg=TRUE) + ggtitle("Prop DIF = .3") + theme(legend.position = "none") + scale_color_brewer(palette="PRGn")
p3 <- visreg(mod, "V3", by="V2", cond=list(V5="0.5"), overlay=T, gg=TRUE) + ggtitle("Prop DIF = .5") + theme(legend.position = "none") + scale_color_brewer(palette="PRGn")
multiplot(p1, p2, p3, cols = 3)

## DO all of the main effects here
visreg(mod, "V1", gg=TRUE) + ggtitle("ME: Magnitude of Dif in Dif")
visreg(mod, "V2", gg=TRUE) + ggtitle("ME: Floor Difficulty")
visreg(mod, "V3", gg=TRUE) + ggtitle("ME: Floor Discrimination")
visreg(mod, "V5", gg=TRUE) + ggtitle("ME: Prop of items with DIF")
visreg(mod, "V6", gg=TRUE) + ggtitle("ME: Sample Size")
visreg(mod, "V8", gg=TRUE) + ggtitle("ME: Magnitude of Dif in Dis")

## Now do the bar plots
plot.vals <- summarySE(data=all.anova.vals, measurevar = "V7", groupvars = c("V2", "V3", "V5", "V8"))
tmp.plot <- ggplot(plot.vals, aes(x=V8, y=V7, fill=V3, color=V3)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=V7-se, ymax=V7+se),position="dodge", stat="identity") +
  facet_grid(V2~V5) +
  coord_cartesian(ylim=c(.7, 1))


## FInd the range of values
tmp.vals <- summarySE(data = all.anova.vals, measurevar = "V7", groupvars = c("V1","V2","V3","V4","V5","V6","V8"))
