## This file is going to be used to explore  patterns with the items that show either uDIF or nuDIF
## this is toatlly exploratory... and I need to find a good way to analyze it, but I think it might be interesting?

## Declare/load libraries
source("~/adroseHelperScripts/R/afgrHelpFunc.R")
install_load("reshape2", "progress", "mirt", "psych", "ggplot2", "polycor", "corrplot")

## Declare any functions
binary.flip <- function (x) {
  x * -1 + 1
}

## load data
dif.vals <- read.csv("./eegBehavioralData/eegDIFItems.csv", na.strings = ".")
load("./fname.gz")
ds_eegPred <- out.data[[9]]

## RT data will take a little management to obtain.. this will all be done below
load("./eegBehavioralData/IRTCareData.RData")
real.all.wide <- as.data.frame(real.all.wide.rt)
real.all.wide.1 <- real.all.wide[,c(1, grep("1_", names(real.all.wide)))]
real.all.wide.2 <- real.all.wide[,c(1, grep("2_", names(real.all.wide)))]
real.all.wide.3 <- real.all.wide[,c(1, grep("3_", names(real.all.wide)))]
real.all.wide.4 <- real.all.wide[,c(1, grep("4_", names(real.all.wide)))]

names(real.all.wide.1) <- gsub(names(real.all.wide.1), pattern = '1_', replacement = '')
names(real.all.wide.2) <- gsub(names(real.all.wide.2), pattern = '2_', replacement = '')
names(real.all.wide.3) <- gsub(names(real.all.wide.3), pattern = '3_', replacement = '')
names(real.all.wide.4) <- gsub(names(real.all.wide.4), pattern = '4_', replacement = '')

## Now combine em
for.irt.rt <- rbind(real.all.wide.1, real.all.wide.2, real.all.wide.3, real.all.wide.4)
for.irt.rt[,2:97] <- apply(for.irt.rt[,2:97], c(1,2), function(x) as.numeric(as.character(x)))
frt.long <- melt(for.irt.rt, id.vars = "V1")

## Clean up the item names
frt.long$variable <- strSplitMatrixReturn(frt.long$variable, "_")[,2]
frt.long$variable <- gsub(frt.long$variable, pattern = ',', replacement = ' ')
## Now combine all of the data
all.data <- merge(frt.long, dif.vals, by.x="variable", by.y="itemFull", all=T)
## Now clean up the NA for the DIF and nuDIF
all.data$uDIF[is.na(all.data$uDIF)] <- 0
all.data$nuDIF[is.na(all.data$nuDIF)] <- 0

## Now add the eeg data
all.data.eeg <- merge(all.data, ds_eegPred, by.x="V1", by.y="record_id", all=T)

## DO any modeling down here
mod.rt.diff <- lmerTest::lmer(value ~ uDIF + (1|variable / V1), data=all.data)
mod.rt.diff <- lmerTest::lmer(value ~ nuDIF + (1|variable / V1), data=all.data)
mod.rt.diff.eeg <- lmerTest::lmer(value ~ (uDIF + P2crymaxCarepost + P2crymaxCareLatpost)^3 + (1|variable/V1), data=all.data.eeg)
