## This script will be used to calculate a variable efficiency value taking into account which of the questions have acceptable IRT parameters
## It will go through calculation of the IRT params within each emotion for the ID task - and then throw out the faces with unaccpetable discrimination

## Clean everything
rm(list=ls())

## Load library(s) Declare function(s)
source("~/adroseHelperScripts/R/afgrHelpFunc.R")
install_load("reshape2", "progress", "mirt", "psych", "ggplot2", "polycor", "corrplot")

binary.flip <- function (x) {
  x * -1 + 1
}

confusionMatrixLeft <- function(df, Right=FALSE, timeMax=FALSE){
  ## Convert the values to the theoretical responses
  ## Here are the response keys
  # A|2 <- Happy
  # S|3 <- Neutral
  # K|4 <- Unhappy
  # L|5 <- Crying
  ## Here are the stim keys
  # 1 <- Crying
  # 6 <- Unhappy
  # 7 <- Neutral
  # 8 <- Happy
  
  ## Define some vars
  in.response <- as.character(df$TriNoAns)
  in.stim <- as.character(df$TriNoStim)
  resp.time <- df$timeToAnswer
  
  ## Modify the null responses
  resp.time[resp.time==0] <- NA
  resp.time <- resp.time * 10e-7
  
  ## Now check for time max flag
  if(timeMax==TRUE){
    # Now change all of the na values to the maximum allowed time
    index <- which(is.na(resp.time))
    # Now modify these values
    resp.time[index] <- ((df$TimeStim[index+1] - df$TimeStim[index]) - 1) * 10e-7
  }
  
  ## First convert the stimuli
  in.key <- c(1,6,7,8)
  out.key <- c("Crying", "Unhappy", "Neutral", "Happy")
  out.stim <- in.stim
  index <- 1
  for(i in in.key){
    out.stim[which(out.stim==i)] <- out.key[index]
    index <- index + 1
  }
  ## Now do the response
  # Now check to see if it is a right or left happy entry
  in.key <- c(2,3,4,5)
  out.key <- c("Happy", "Neutral", "Unhappy", "Crying")
  if(Right){
    out.key <- rev(out.key)
  }
  out.resp <- in.response
  index <- 1
  for(i in in.key){
    out.resp[which(out.resp==i)] <- out.key[index]
    index <- index + 1
  }
  ## Now make sure the out response have the right number of levels
  flag.check <- setdiff(out.resp, out.key)
  if(length(flag.check)>1){
    ## Set unwanted values to ""
    # Find which flag.check is equal to ""
    index <- which(flag.check=="")
    ## Now change all of those that are flagged to ""
    flag.check <- flag.check[-index]
    out.resp[out.resp %in% flag.check] <- ""
  }
  
  ## Now create the table
  out.table <- table(out.stim, out.resp, useNA = 'ifany')
  
  ## Now I need to calculate the efficiency
  # I am going to do this overall and also within each emotion
  total.correct <- sum(out.stim == out.resp, na.rm=T)
  total.percent <- total.correct / length(out.stim)
  total.time <- mean(resp.time, na.rm = T)
  out.vals <- data.frame("contrast" = c("All"), "correct" = c(total.correct), "time" = c(total.time), "percent" = c(total.percent))
  out.vals$contrast <- as.character(out.vals$contrast)
  # Now do the emotion specific sums
  for(i in out.key){
    tmp.index <- which(out.stim == i)
    tmp.correct <- sum(out.stim[tmp.index] == out.resp[tmp.index], na.rm=T)
    tmp.time <- mean(resp.time[tmp.index], na.rm=T)
    tmp.percent <- tmp.correct / length(tmp.index)
    out.vals <- rbind(out.vals, list(i, tmp.correct, tmp.time, tmp.percent))
  }
  ## Now calculate the efficiency
  out.vals$efficiency <- out.vals$correct / out.vals$time
  
  ## Now return a vector with the correct incorrect indicies
  out.vector <- which(out.stim == out.resp)
  
  ## Now return the vector of number of responses
  
  ## Now perpare a list with all of the output
  return(out.list <- list(out.table=out.table, out.vals=out.vals, out.vec=out.vector, out.resp=df$noOfResposesStim, out.stim=out.stim))
}

timeForError <- function(df){
  ## First grab the difference between the stim and the ans
  indexStim <- 2:dim(df)[1]
  ## Now calc the difference
  downTime <- NULL
  downTime <- df$TimeStim[indexStim] - df$TimeAns[indexStim-1]
  downTime <- c(NA, downTime)
  ## Now flag for NA's
  # If a participant has a no response then it is impossible to have an error judgment
  # So these error judgments need to be flagged as NA's
  # When an individual has a null response @ lower index (answer) than the upper index (stim) is invalid
  null.answer <- which(is.na(df$TriNoAns))
  downTime[null.answer] <- NA
  
  ## Now attach to the df
  df$downTime <- downTime
  
  ## Now return the new df
  return(df)
  
}



## read data
load("./eegBehavioralData/IRTIDData.RData")
load("./fname.gz")
ds_eegPred <- out.data[[9]]

## Check for the 

# first orgainze the data
real.all.wide <- as.data.frame(real.all.wide2)
real.all.wide.1 <- real.all.wide[,c(1, grep("1_", names(real.all.wide)))]
real.all.wide.2 <- real.all.wide[,c(1, grep("2_", names(real.all.wide)))]
real.all.wide.3 <- real.all.wide[,c(1, grep("3_", names(real.all.wide)))]
real.all.wide.4 <- real.all.wide[,c(1, grep("4_", names(real.all.wide)))]

names(real.all.wide.1) <- gsub(names(real.all.wide.1), pattern = '1_', replacement = '')
names(real.all.wide.2) <- gsub(names(real.all.wide.2), pattern = '2_', replacement = '')
names(real.all.wide.3) <- gsub(names(real.all.wide.3), pattern = '3_', replacement = '')
names(real.all.wide.4) <- gsub(names(real.all.wide.4), pattern = '4_', replacement = '')

## Now combine these and run IRT
for.irt <- rbind(real.all.wide.1, real.all.wide.2, real.all.wide.3, real.all.wide.4)
# fix the factor levels
for.irt[,2:97] <- apply(for.irt[,2:97], 2, function(x) factor(x, levels=c("happy", "neutral", "unhappy", "crying")))
for.irt2 <- for.irt
## Now go through and change all of these characters into numerics
# Cry == 1; Unhappy == 2; Neutral ==3; Happy == 4
new.val <- 1
for(i in c("crying", "unhappy", "neutral", "happy")){
  ## Grab the indicies
  index.vals <- which(for.irt2==i)
  index.vals <- matrix(rc.ind(for.irt2, index.vals), byrow=F, ncol=2)
  for(L in 1:dim(index.vals)[1]){
    for.irt2[index.vals[L,1], index.vals[L,2]] <- new.val
  }
  new.val <- new.val + 1
}
for.irt2[,2:97] <- apply(for.irt2[,2:97], 2, as.numeric)
## Now train the model
mod.1 <- mirt(data=for.irt2[,2:97], itemtype = 'nominal', SE=T, model=1)

## Now do a binary () correct vs incorrect) IRT model
for.irt3 <- matrix(NA, ncol = 96, nrow = dim(for.irt)[1])
# I need to go through and find the correct string for each column and then identify the correct responses within each cell
for(i in 2:97){
  string.index <- for.irt[,i]
  irt_three_col <- i -1
  ## First grab the correct index
  # Grab the character from the colname
  char.val <- tolower(substr(strSplitMatrixReturn(colnames(for.irt[2]), "_")[,2][1], 1, 1))
  if(char.val == 'c') cor.index <- which(string.index=="crying")
  if(char.val == 'n') cor.index <- which(string.index=="neutral")
  if(char.val == 'u') cor.index <- which(string.index=="unhappy")
  if(char.val == 'h') cor.index <- which(string.index=="happy")
  ## Now change the values
  for.irt3[which(!is.na(for.irt[,i])),irt_three_col] <- 0
  for.irt3[cor.index,irt_three_col] <- 1
}
for.irt3 <- as.data.frame(for.irt3)
colnames(for.irt3) <- colnames(for.irt2)[2:97]
## Now run the IRT with the binary data
# I will run this within the intensity bins as this is where we se the best fit
mod.11 <- mirt(for.irt3[,1:48], 1, SE=T, itemtype='2PL')
plot(mod.11, type='trace')
mod.13 <- mirt(for.irt3[,49:96], 1, SE=T, itemtype='2PL')
plot(mod.13, type='trace')


## Now isolate the discrimination for these items
high.Inten.Values <- coef(mod.11, IRTpars=T, simplify=T, as.data.frame=F)$items
low.Inten.Values <- coef(mod.13, IRTpars=T, simplify=T, as.data.frame=F)$items


## Now find the items to remove
# I am going to start with discrimination values less than 1
# The majority of these come from the low I items
bad.items <- names(c(which(abs(low.Inten.Values[,"a"]) < 1), which(abs(high.Inten.Values[,"a"]) < 1)))

## Now go back and recalc the total number of correct responses while removing these items
# First find all of the columns to remove
bad.cols <- NULL
for(i in 1:length(bad.items)){
  bad.cols <- c(bad.cols, grep(bad.items[i], colnames(real.all.wide2)))
}
recalc.vals <- real.all.wide2[,-bad.cols]
## Now find the correct vs incorrect responses
calc.vals.cor <- recalc.vals
for(i in 2:dim(recalc.vals)[2]){
  ## First find the correct string
  char.val <- tolower(substr(strSplitMatrixReturn(colnames(recalc.vals[i]), "_")[,3][1], 1, 1))
  if(char.val == 'c') cor.index <- "crying"
  if(char.val == 'n') cor.index <- "neutral"
  if(char.val == 'u') cor.index <- "unhappy"
  if(char.val == 'h') cor.index <- "happy"
  ## Now change all values == the string index to 1; and not equal -> 0; while preserving the NA values
  calc.vals.cor[,i] <- as.character(calc.vals.cor[,i])
  calc.vals.cor[which(recalc.vals[,i]==cor.index),i] <- 1
  calc.vals.cor[which(recalc.vals[,i]!=cor.index),i] <- 0
  calc.vals.cor[which(is.na(recalc.vals[,i])),i] <- NA
  calc.vals.cor[,i] <- as.numeric(calc.vals.cor[,i])
}

## Now calculate the sums for each emotion and overall
overallSum <- rowSums(calc.vals.cor[,2:313], na.rm = T)
## Now get the sums within emotion
crying.sum <- rowSums(calc.vals.cor[,grep("_C", colnames(calc.vals.cor))], na.rm = T)
happy.sum <- rowSums(calc.vals.cor[,grep("_H", colnames(calc.vals.cor))], na.rm = T)
neutral.sum <- rowSums(calc.vals.cor[,grep("_N", colnames(calc.vals.cor))], na.rm = T)
unhappy.sum <- rowSums(calc.vals.cor[,grep("_U", colnames(calc.vals.cor))], na.rm = T)
# now do the percentages
percent.overall <- overallSum /  apply(!is.na(calc.vals.cor[,2:313]), 1, sum)
percent.crying <- crying.sum /   apply(!is.na(calc.vals.cor[,grep("_C", colnames(calc.vals.cor))]), 1, sum)
percent.happy <- happy.sum /     apply(!is.na(calc.vals.cor[,grep("_H", colnames(calc.vals.cor))]), 1, sum)
percent.neutral <- neutral.sum / apply(!is.na(calc.vals.cor[,grep("_N", colnames(calc.vals.cor))]), 1, sum)
percent.unhappy <- unhappy.sum / apply(!is.na(calc.vals.cor[,grep("_U", colnames(calc.vals.cor))]), 1, sum)
percent.data <- cbind(as.character(calc.vals.cor[,1]), percent.overall, percent.crying, percent.happy, percent.neutral, percent.unhappy)
colnames(percent.data)[1] <- "record_id"
percent.data <- as.data.frame(percent.data)

## Now run the lmer with the updated values
lm.data <- merge(ds_eegPred, percent.data, by="record_id", all=T)
id.vars <- c("record_id", "ace_score", "dose_hv_visit_count", "income_h4")
vars.of.int <- colnames(percent.data)[-c(1,2)]
data.iso <- lm.data[,c(id.vars, vars.of.int)]
data.iso[,5:8] <- scale(data.iso[,5:8])[,]
x <- reshape2::melt(data.iso, id.vars = id.vars)
x$value <- as.numeric(x$value)
## Add valence and intensity values
x$valence <- "+"
x$valence[grep("unhappy", x$variable)] <- "-"
x$valence[grep("cry", x$variable)] <- "-"

## Now do intensity
x$intensity <- "strong"
x$intensity[grep("unhappy", x$variable)] <- "weak"
x$intensity[grep("neutral", x$variable)] <- "weak"

# Train the model
mod.continous <- lmerTest::lmer(value ~ (ace_score+dose_hv_visit_count+intensity+valence)^4 + (1|record_id), data=x)
mod.continous.2 <- lmerTest::lmer(value ~ (ace_score+dose_hv_visit_count+variable)^3 + (1|record_id), data=x)

mod.one <- lm(as.numeric(as.character(percent.crying)) ~ ace_score * dose_hv_visit_count,  data=data.iso)
mod.two <- lm(as.numeric(as.character(percent.happy)) ~ ace_score * dose_hv_visit_count,   data=data.iso)
mod.thr <- lm(as.numeric(as.character(percent.neutral)) ~ ace_score * dose_hv_visit_count, data=data.iso)
mod.fou <- lm(as.numeric(as.character(percent.unhappy)) ~ ace_score * dose_hv_visit_count, data=data.iso)
  
## Now get the reaction time values
recalc.rt.time <- as.data.frame(real.all.wide.rt)



