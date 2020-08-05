## This script will be used to produce the response probabilities for the care task
## This is what nick sent me though
#Stim Trigs (Care)
#3 (Crying)
#4 (Unhappy)
#5 (Neutral)
#6 (Happy)

#Response Trigs (Care)
#1 ("L Shift") --> "Yes"
#2 ("R Shift") --> "No" ** I believe, I need to double check to see if this is correct 20200723

## Turns out the yes and right flip....
## So I need to go find which ever one they endorse crying the most with?


## Clean everything
rm(list=ls())

## Load library(s)
source("~/adroseHelperScripts/R/afgrHelpFunc.R")
install_load("tidyverse", "reshape2", "gganimate", 'plotly', 'sjPlot')

## Create function(s)
convertToResponse <- function(df, rightYes=TRUE){
  ## First convert the stimuli
  #Stim Trigs (Care)
  #3 (Crying)
  #4 (Unhappy)
  #5 (Neutral)
  #6 (Happy)
  
  ## Define some vars
  in.response <- as.character(df$TriNoAns)
  in.stim <- as.character(df$TriNoStim)
  resp.time <- df$timeToAnswer
  
  ## Modify the null responses
  resp.time[resp.time==0] <- NA
  resp.time <- resp.time * 10e-7
  
  ## First convert the stimuli
  in.key <- c(3,4,5,6)
  out.key <- c("Crying", "Unhappy", "Neutral", "Happy")
  out.stim <- in.stim
  index <- 1
  for(i in in.key){
    out.stim[which(out.stim==i)] <- out.key[index]
    index <- index + 1
  }
 
  # Now check to see if it is a right or left happy entry
  in.key <- c(1,2)
  out.key <- c("Yes", "No")
  if(rightYes){
    out.key <- rev(out.key)
  }
  out.resp <- in.response
  index <- 1
  for(i in in.key){
    out.resp[which(out.resp==i)] <- out.key[index]
    index <- index + 1
  }
 
  ## Now prepare an out table of the responses
  out.table <- table(out.resp, out.stim, useNA = "ifany")
  #out.table <- melt(out.table)
  
  ## Now preapre mean resp times within each yes/no permutation
  out.time <- summarySE(data.frame(out.stim, out.resp, resp.time), measurevar = "resp.time", groupvars = c("out.stim", "out.resp"), na.rm = TRUE)
  # Now quickly remove any NA rows if there are any
  if(sum(is.na(out.time$out.resp))>0){
    # Now remove the rows w/ NA
    out.time <- out.time[complete.cases(out.time),]
  }
  
  ## Now return all of these values
  return(out.list <- list(out.table=out.table, out.time=out.time))
}


## Now create a function which will grab the time between stimuli and previous answer
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

## Load data
in.path <- "/home/arosen/Documents/bbmcPersonal/eegBehavioralData/CARE_mod/"
in.example <- "/home/arosen/Documents/bbmcPersonal/eegBehavioralData/CARE_mod/117-114_OutResponse.csv"
in.data <- read.csv(in.example)

## Now in a loop go through and find the response pattern for each emotion for each person
all.id <- system(paste("ls", in.path), intern = TRUE)
out.check <- NULL
for(q in all.id){
  in.file <- paste(in.path, q, sep='')
  in.data <- read.csv(in.file)
  left.vals <- table(in.data$TriNoStim, in.data$TriNoAns)[1,]
  out.row <- c(q, left.vals)
  out.check <- rbind(out.check, out.row)
}
out.check <- data.frame(out.check)
out.check$X1 <-as.numeric(as.character(out.check$X1))
out.check$X2 <-as.numeric(as.character(out.check$X2))
out.check$Left <- "No"
out.check$Left[which(out.check$X1 - out.check$X2 > 0)] <- "Yes"

## Now load all of the data
all.data <- list()
all.response <- list()
all.time <- NULL
for(i in 1:dim(out.check)[1]){
  q <- all.id[i]
  in.file <- paste(in.path, q, sep='')
  in.data <- read.csv(in.file)
  right.flag <- FALSE
  if(out.check$Left[i]=="No"){
    right.flag <- TRUE
  }
  resp.table <- convertToResponse(in.data, rightYes = right.flag)
  check <- melt(resp.table$out.table)
  all.response[[i]] <- resp.table$out.table
  tmp.time <- resp.table$out.time
  tmp.time$ID <- q
  all.time <- rbind(all.time, tmp.time)
}
## Now organize the data
all.resp <- dcast(melt(lapply(all.response, as.list)), L1 ~L2)
# Now fix the column names
colnames(all.resp)[1] <- "ID"
colnames(all.resp)[2:4] <- paste(colnames(resp.table$out.table)[1], rownames(resp.table$out.table), sep = '_')
colnames(all.resp)[5:7] <- paste(colnames(resp.table$out.table)[2], rownames(resp.table$out.table), sep = '_')
colnames(all.resp)[8:10] <- paste(colnames(resp.table$out.table)[3], rownames(resp.table$out.table), sep = '_')
colnames(all.resp)[11:13] <- paste(colnames(resp.table$out.table)[4], rownames(resp.table$out.table), sep = '_')
all.resp$ID <- all.id

## Now organize the time vars
all.time$perm <- paste(all.time$out.stim, all.time$out.resp, sep='_')
all.time.mean <- spread(all.time[,c("ID", "perm", "resp.time")], perm, resp.time)
all.time.n <- spread(all.time[,c("ID", "perm", "N")], perm, N)
# Convert all NAs in the n df to 0 values
all.time.n[is.na(all.time.n)] <- 0
all.time.var <- spread(all.time[,c("ID", "perm", "sd")], perm, sd)
colnames(all.time.var)[2:9] <- paste(colnames(all.time.var)[2:9], "_var", sep='')

## Now merge these all together
all.time.out <- merge(all.time.mean, all.time.n, by="ID", suffixes = c("_mean", "_n"), all=T)
all.time.out <- merge(all.time.out, all.time.var, by="ID", all=T)

all.out <- merge(all.resp, all.time.out, by="ID")
## Now fix the id var
all.out$ID <- strSplitMatrixReturn(all.out$ID, "_")[,1]

## Now write the csv
write.csv(all.out, "careData-RAW-FirstResponse.csv", quote=F, row.names=F)

#### Of note - ID 117-114 appears to be wildly broken!!!

## Now plot and explore these values
# First plot the happy response counts
pdf("careResponsePatterns.pdf")
## Crying Response
melt(all.out[,c("ID", "Crying_Yes", "Crying_No", "Crying_NA")], id.vars = "ID") %>% ggplot(., aes(x=value)) +
  geom_histogram() +
  facet_grid(.~variable, scales = "free") +
  theme_bw()
## Crying speed
melt(all.out[,c("ID", "Crying_Yes_mean", "Crying_No_mean")], id.vars = "ID") %>% ggplot(., aes(x=value)) +
  geom_histogram() +
  facet_grid(.~variable, scales = "free") +
  theme_bw()

## Happy Response
melt(all.out[,c("ID", "Happy_Yes", "Happy_No", "Happy_NA")], id.vars = "ID") %>% ggplot(., aes(x=value)) +
  geom_histogram() +
  facet_grid(.~variable, scales = "free") +
  theme_bw()
## Happy speed
melt(all.out[,c("ID", "Happy_Yes_mean", "Happy_No_mean")], id.vars = "ID") %>% ggplot(., aes(x=value)) +
  geom_histogram() +
  facet_grid(.~variable, scales = "free") +
  theme_bw()

## Neutral Response
melt(all.out[,c("ID", "Neutral_Yes", "Neutral_No", "Neutral_NA")], id.vars = "ID") %>% ggplot(., aes(x=value)) +
  geom_histogram() +
  facet_grid(.~variable, scales = "free") +
  theme_bw()
## Neutral speed
melt(all.out[,c("ID", "Neutral_Yes_mean", "Neutral_No_mean")], id.vars = "ID") %>% ggplot(., aes(x=value)) +
  geom_histogram() +
  facet_grid(.~variable, scales = "free") +
  theme_bw()

## Unhappy Response
melt(all.out[,c("ID", "Unhappy_Yes", "Unhappy_No", "Unhappy_NA")], id.vars = "ID") %>% ggplot(., aes(x=value)) +
  geom_histogram() +
  facet_grid(.~variable, scales = "free") +
  theme_bw()
## Unhappy speed
melt(all.out[,c("ID", "Unhappy_Yes_mean", "Unhappy_No_mean")], id.vars = "ID") %>% ggplot(., aes(x=value)) +
  geom_histogram() +
  facet_grid(.~variable, scales = "free") +
  theme_bw()

dev.off()
