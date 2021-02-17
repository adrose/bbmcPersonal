## This will be a small script which will be used to look at 
## interrealtions amongst the P200 varaibility across the emotions
## It is essentially just going to make a correlation matrix
## I will probablly throw this script into the organizeIDItemData.R code

# Clear ram
rm(list=ls())

## Load library(s)
library(ggplot2)
library(GGally)
library(reshape2)
library(ggalluvial)
source("~/adroseHelperScripts/R/afgrHelpFunc.R")

## Load data
load("./fname.gz")
ds_eegPred <- out.data[[9]]



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

## Now plot the amount of movement across these groups; this will be done with a river plot
## This is going to be complicated... I need to create all combinations of data?
alluvial.data <- melt(table(mod.plot$plot.valsCry, mod.plot$plot.valsHap, mod.plot$plot.valsNeu, mod.plot$plot.valsUnh))
## Now prepare a cry specific river data
alluvial.data.cry1 <- melt(table(mod.plot$plot.valsCry, mod.plot$plot.valsHap))
alluvial.data.cry1$target <- "Hap"
alluvial.data.cry2 <- melt(table(mod.plot$plot.valsCry, mod.plot$plot.valsNeu))
alluvial.data.cry2$target <- "Neu"
alluvial.data.cry3 <- melt(table(mod.plot$plot.valsCry, mod.plot$plot.valsUnh))
alluvial.data.cry3$target <- "Unh"
alluvial.data.cry <- rbind(alluvial.data.cry1, alluvial.data.cry2, alluvial.data.cry3)
## Now clean all empty cells
alluvial.data <- alluvial.data[which(alluvial.data$value!=0),]
alluvial.data.cry <- alluvial.data.cry[which(alluvial.data.cry$value!=0),]
colnames(alluvial.data)[1:4] <- c("Cry", "Happy", "Neutral", "Unhappy")
ggplot(alluvial.data, aes(y = value, axis1=Cry, axis2=Unhappy, axis3=Neutral, axis4=Happy)) +
  geom_alluvium(width=1/12) +
  geom_stratum(width=1/12, fill='black', color='gray') +
  geom_label(stat="stratum", aes(label=after_stat(stratum))) +
  scale_x_discrete(limits=c("Cry", "Happy", "Neutral", "Unhappy"))+
  theme_bw() +
  ggtitle("Group EEG Pattern Progression")
ggplot(alluvial.data.cry, aes(y=value, axis1=Var1, axis2=Var2)) +
  geom_alluvium(width=1/12) +
  geom_stratum(width=1/12, fill='black', color='gray') +
  geom_label(stat="stratum", aes(label=after_stat(stratum))) +
  facet_grid(. ~ target) +
  scale_x_discrete(limits=c("1", "2"))


## Now prepare a correlation matrix down here
## First isolate the data of interest
## I also need to create the interaction varaibles and add these to the matrix
cor.mat.dat <- grep("^P2", names(ds_eegPred), value = T, perl = T)
cor.mat.dat <- grep("post$", cor.mat.dat, value = T, perl = T)
## Now go through and create the interaction terms for all of these
emo.vals <- c("all", "happy", "neutral", "unhappy", "cry")
task.vals <- c("ID", "Care")
for(t in task.vals){
  for(e in emo.vals){
    str.one <- paste("P2", e, "max", t, "post", sep='')
    #P2allmaxIDLatfront
    str.two <- paste("P2", e,"max", t, "Lat","post", sep='')
    str.out <- paste("P2", e, "max", t, "Int","post", sep='')
    val.out <- scale(ds_eegPred[,str.one]) * scale(ds_eegPred[,str.two])
    ds_eegPred[,str.out] <- val.out
    cor.mat.dat <- c(cor.mat.dat, str.out)
  }
}

to.plot <- scale(ds_eegPred[,cor.mat.dat])[,]
# Now grab the cor vals
cor.vals <- melt(cor(to.plot, use='complete'))
# Now clean them
# Start with the first col
cor.vals$Var1.emo <- strSplitMatrixReturn(strSplitMatrixReturn(charactersToSplit = cor.vals$Var1, splitCharacter = "P2")[,2], "max")[,1]
cor.vals$Var1.val <- "Amp"
cor.vals$Var1.val[grep("Lat", cor.vals$Var1)] <- "Lat"
cor.vals$Var1.val[grep("Int", cor.vals$Var1)] <- "Int"
cor.vals$Val1.task <- "I"
cor.vals$Val1.task[grep("Care", cor.vals$Var1)] <- "C"
# Now do the second col
cor.vals$Var2.emo <- strSplitMatrixReturn(strSplitMatrixReturn(charactersToSplit = cor.vals$Var2, splitCharacter = "P2")[,2], "max")[,1]
cor.vals$Var2.val <- "Amp"
cor.vals$Var2.val[grep("Lat", cor.vals$Var2)] <- "Lat"
cor.vals$Var2.val[grep("Int", cor.vals$Var2)] <- "Int"
cor.vals$Val2.task <- "I"
cor.vals$Val2.task[grep("Care", cor.vals$Var2)] <- "C"

# Now prepare the column names
cor.vals$Var1CN <- paste(cor.vals$Val1.task, cor.vals$Var1.emo, cor.vals$Var1.val, sep="_")
cor.vals$Var2CN <- paste(cor.vals$Val2.task, cor.vals$Var2.emo, cor.vals$Var2.val, sep="_")
cor.vals$value2 <- round(cor.vals$value, 2)
cor.vals$value3 <- round(cor.vals$value, 1)

## Now add a facet
cor.vals$facet <- "Pure ID"
cor.vals$facet[which(cor.vals$Val1.task=="C" & cor.vals$Val2.task=="C")] <- "Pure Care"
cor.vals$facet[which(cor.vals$Val1.task=="I" & cor.vals$Val2.task=="C")] <- "Mixed"
cor.vals$facet[which(cor.vals$Val1.task=="C" & cor.vals$Val2.task=="I")] <- "Mixed"


## Now plot
pal <- wesanderson::wes_palette("Zissou1", type = "continuous")
tmp <- ggplot(data=cor.vals) +
  theme_minimal() +
  geom_tile(aes(x=Var2CN, y=Var1CN, fill=value)) +
  geom_text(data=cor.vals,aes(y=Var1CN, x=Var2CN, label=value3), check_overlap = F) +
  coord_fixed() +
  scale_fill_gradientn(limits=c(-1,1), colours = pal) + 
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  xlab("") +
  ylab("")


tmp2 <- ggplot(data=cor.vals[-which(cor.vals$facet=="Mixed"),]) +
  theme_minimal() +
  geom_tile(aes(x=Var2CN, y=Var1CN, fill=value)) +
  geom_text(aes(y=Var1CN, x=Var2CN, label=value2), check_overlap = T) +
  scale_fill_gradientn(limits=c(-1,1), colors = pal) + 
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  facet_wrap(facet~., scales="free") +
  xlab("") +
  ylab("")
tmp
tmp2


