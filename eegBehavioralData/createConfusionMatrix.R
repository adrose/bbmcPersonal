## This script will be used to produce the confusion matricies for the ID task
## It will also give efficiency estimates for individuals
## Important to note that the order of answers varies between individuals
## There is a left presenation where happy is on the right and one where it is on the left
## I am still unsure what this means though
## This is what nick sent me though
## Stim Trigs (ID)
## 1 (Crying)
## 6 (Unhappy)
## 7 (Neutral)
## 8 (Happy)
## 
## Response Trigs (ID)
## 2 ("A")
## 3 ("S")
## 4 ("K")
## 5 ("L")

## Load library(s)
source("~/adroseHelperScripts/R/afgrHelpFunc.R")
install_load("tidyverse", "reshape2", "gganimate", 'plotly', 'sjPlot')

## Create functions
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
in.path <- "/home/arosen/Documents/bbmcPersonal/eegBehavioralData/ID_mod/"
in.example <- "/home/arosen/Documents/bbmcPersonal/eegBehavioralData/ID_mod/118-101_OutResponse.csv"
in.data <- read.csv(in.example)

## Now run an example
foo <- confusionMatrixLeft(in.data)

## Now in a loop go through and try to identify if happy right or happy left based on the overall correct responses
all.id <- system(paste("ls", in.path), intern = TRUE)
out.check <- NULL
for(q in all.id){
  in.file <- paste(in.path, q, sep='')
  in.data <- read.csv(in.file)
  left.vals <- confusionMatrixLeft(in.data, Right = FALSE)$out.vals[1,'correct']
  right.vals <- confusionMatrixLeft(in.data, Right = TRUE)$out.vals[1,'correct']
  out.row <- c(q, left.vals, right.vals)
  out.check <- rbind(out.check, out.row)
}
out.check <- as.data.frame(out.check)
out.check$V1 <- as.character(out.check$V1)
out.check$V2 <- as.numeric(as.character(out.check$V2))
out.check$V3 <- as.numeric(as.character(out.check$V3))
out.check$Guess <- "Right"
out.check$Guess[which(out.check$V2 - out.check$V3 > 0)] <- "Left"

## Now re run this with the guess as the input direction
all.data <- list()
all.response <- NULL
all.vals <- NULL
downTime <- NULL
for(i in 1:dim(out.check)[1]){
  q <- all.id[i]
  in.file <- paste(in.path, q, sep='')
  in.data <- read.csv(in.file)
  right.flag <- FALSE
  if(out.check$Guess[i]=="Right"){
    right.flag <- TRUE
  }
  in.data.new <- timeForError(in.data)
  if(identical(in.data$downTime, NULL)){
    write.csv(in.data.new, in.file, quote=F, row.names=F)
    print("Writing")
  }
  cor.table <- confusionMatrixLeft(in.data, Right = right.flag, timeMax = TRUE)
  all.data[[i]] <- cor.table
  all.response.tmp <- melt(cor.table$out.table)
  all.response.tmp$id <- q
  all.response <- rbind(all.response, all.response.tmp)
  all.vals.tmp <- melt(cor.table$out.vals)
  all.vals.tmp$id <- q
  all.vals <- rbind(all.vals, all.vals.tmp)
  ## Now create a bool array with the correct incorrect questions
  out.bool <- rep(0, length(in.data.new$downTime))
  out.bool[cor.table$out.vec] <- 1
  downTime.tmp <- cbind(q, in.data.new$downTime, out.bool, cor.table$out.resp, cor.table$out.stim)
  downTime <- rbind(downTime, downTime.tmp)
}
## Clean up the ID for all.vals
all.vals$id <- strSplitMatrixReturn(all.vals$id, "_")[,1]

## Now create an error time value explore
downTime <- data.frame(downTime)
downTime$q <- strSplitMatrixReturn(downTime$q, "_")[,1]
downTime$V2 <- as.numeric(as.character(downTime$V2))

## Now create a plot of all of these values
downTime$V2 <- downTime$V2 * 10e-7
downTime$q <- as.factor(downTime$q)
## Quick clean
downTime$V2[downTime$V2>2.5] <- NA
q <- downTime[which(downTime$out.bool==0),] %>% ggplot(., aes(x = V2, frame=q)) +
  geom_histogram(stat='bin') +
  theme_bw() +
  ggtitle("Reponse and Presentation Interval") +
  facet_wrap(vars(q), scales = "free")

to.write <- summarySE(data = downTime[which(downTime$out.bool==0),], measurevar = "V2", groupvars = c("q", "V4"), na.rm=T)
to.write <- to.write[-(which(to.write$V4==0 & to.write$N==0)),]
write.csv(to.write, "downTimeIncorrectAnswersWithMultiple.csv", quote=F, row.names=F)

## Now without # of responses
to.write <- summarySE(data = downTime[which(downTime$out.bool==0),], measurevar = "V2", groupvars = c("q"), na.rm=T)
write.csv(to.write, "downTimeIncorrectAnswers.csv", quote=F, row.names=F)

## Now check for differences between correct and incorrect response times within subjects
mem.mod <- lmerTest::lmer(V2 ~ out.bool + (1|q), data=downTime)
## Plot our random effect

pdf("exploreHisto.pdf")
## Now z score some of these things
## First isolate the efficiency variables
data.to.score <- all.vals[which(all.vals$variable=="efficiency"),]
data.to.score %>% ggplot(., aes(x=value, color=contrast, fill=contrast)) +
  geom_histogram() +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(. ~contrast, scales = 'free_x') +
  ggtitle("Efficiency histogram (raw)")
data.to.score <- spread(data = data.to.score, value = value, key = contrast)[,-1]
data.to.score[,2:6] <- scale(data.to.score[,2:6])[,]
data.to.score <- melt(data.to.score)
## Now plot it
data.to.score %>% ggplot(., aes(x=value, color=variable, fill=variable)) +
  geom_histogram() +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(. ~variable) +
  ggtitle("Efficiency histogram (z-score)")

## Now do the same for the response time
data.to.score <- all.vals[which(all.vals$variable=="time"),]
data.to.score %>% ggplot(., aes(x=value, color=contrast, fill=contrast)) +
  geom_histogram() +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(. ~contrast) +
  ggtitle("RT histogram (raw)")

data.to.score <- spread(data = data.to.score, value = value, key = contrast)[,-1]
data.to.score[,2:6] <- scale(data.to.score[,2:6])[,]
data.to.score <- melt(data.to.score)
## Now plot it
data.to.score %>% ggplot(., aes(x=value, color=variable, fill=variable)) +
  geom_histogram() +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(. ~variable) +
  ggtitle("RT histogram (z-score)")

## Now do correct responses
data.to.score <- all.vals[which(all.vals$variable=="correct"),]
data.to.score %>% ggplot(., aes(x=value, color=contrast, fill=contrast)) +
  geom_histogram() +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(. ~contrast) +
  ggtitle("Correct responses histogram (raw)")

data.to.score <- spread(data = data.to.score, value = value, key = contrast)[,-1]
data.to.score[,2:6] <- scale(data.to.score[,2:6])[,]
data.to.score <- melt(data.to.score)
## Now plot it
data.to.score %>% ggplot(., aes(x=value, color=variable, fill=variable)) +
  geom_histogram() +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(. ~variable) +
  ggtitle("Correct responses histogram (z-score)")

## Now do percent correct
data.to.score <- all.vals[which(all.vals$variable=="percent"),]
data.to.score %>% ggplot(., aes(x=value, color=contrast, fill=contrast)) +
  geom_histogram() +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(. ~contrast) +
  ggtitle("% correct responses histogram (raw)")

data.to.score <- spread(data = data.to.score, value = value, key = contrast)[,-1]
data.to.score[,2:6] <- scale(data.to.score[,2:6])[,]
data.to.score <- melt(data.to.score)
## Now plot it
data.to.score %>% ggplot(., aes(x=value, color=variable, fill=variable)) +
  geom_histogram() +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(. ~variable) +
  ggtitle("% correct responses histogram (z-score)")
dev.off()
## Now figure some method to plot the confusion matrix
all.response$id <- strSplitMatrixReturn(all.response$id, "_")[,1]
table(all.response$out.resp)
#all.response <- all.response[-which(all.response$out.resp==0),]
all.response$out.resp <- as.character(all.response$out.resp)
all.response$out.resp[is.na(all.response$out.resp)] <- "No Answer"
all.response$out.resp <- factor(all.response$out.resp, levels = c("Crying", "Happy", "Neutral", "Unhappy", "No Answer"))
p <- ggplot(all.response, aes(y=out.stim, x=out.resp, frame=id, fill=value)) +
  geom_tile() +
  coord_equal() +
  theme_bw() +
  ggtitle("ID confusion matrix")
ggplotly(p)

## Now write these data out
all.response.tw <- recast(all.response, id ~ out.stim + out.resp, measure.var = 'value')
all.vals.tw <- recast(all.vals, id ~ variable + contrast, measure.var = "value")
# Now merge and write em
all.out <- merge(all.response.tw, all.vals.tw, by='id')
## Now remove all na columns
#all.out <- all.out[,-grep("NA", names(all.out))]
write.csv(all.out, "behaviorData-RAW-FirstResponse.csv", quote=F, row.names=F)