## Clean everything
rm(list=ls())

## Test script to organize data for IRT analysis using care task faces
source("~/adroseHelperScripts/R/afgrHelpFunc.R")
install_load("reshape2", "progress", "mirt", "psych", "ggplot2")

## Declare any functions
binary.flip <- function (x) {
  x * -1 + 1
}

## Run the bash script to organize the files as desired
system("/home/arosen/Documents/bbmcPersonal/eegBehavioralData/organizeIRTData.sh")

## This is going to have to be done in a loop
all.files.pic <- system("ls ~/Documents/bbmcPersonal/eegBehavioralData/careFace/careMod/pic*", intern = T)

## Now go through this in a loop and merge em
all.out <- NULL
for(i in all.files.pic){
  ## First grab the id
  tmp <- strSplitMatrixReturn(charactersToSplit = i, splitCharacter = '/')[,9]
  ## Now isolate the id
  p1 <- strSplitMatrixReturn(charactersToSplit = tmp, splitCharacter = 'c')[,2]
  p2 <- strSplitMatrixReturn(charactersToSplit = p1, splitCharacter = '\\.')[,1]
  print(p2)
  ## Now decalre the files
  file.1 <- paste("~/Documents/bbmcPersonal/eegBehavioralData/careFace/careMod/pic", p2,".csv", sep='')
  file.2 <- paste("~/Documents/bbmcPersonal/eegBehavioralData/careFace/careMod/resp", p2,".csv", sep='')
  in.file.1 <- read.csv(file.1, header=F, sep='\t')
  in.file.2 <- read.csv(file.2, header=F, sep='\t')
  all.data <- merge(in.file.1, in.file.2, by=c("V1", "V2"))
  #print(dim(all.data))
  ## Now I need to go through and find all of the instances with 2 responses and take the later responses
  ## First find if they have any two response counts
  if(dim(table(table(all.data$V2)))>1){
    print("Duplicate")
    ## Now find which rows are duplicates
    row.vals <- names(which(table(all.data$V2)>1))
    ## Now go through each of the row vals and take the latest response
    # First create an index for the rows to keep
    rows.to.lose <- NULL
    for(w in row.vals){
      # Find max response time
      rt.max <- max(as.numeric(as.character(all.data[which(all.data$V2==w),"V6.y"])))
      row.index <- which(as.numeric(as.character(all.data[which(all.data$V2==w),"V6.y"]))==rt.max)
      row.to.vom <- which(all.data$V2==w)[-row.index]
      rows.to.lose <- c(rows.to.lose, row.to.vom) 
    }
    all.data <- all.data[-rows.to.lose,]
  }
  ## Now fix column names
  colnames(all.data) <- c("fileID", "itemVal", "Pre", "PictureDisplayed", "timeShown", "Null", "Null2", "respVal", "ResponseGiven", "timeGiven", "responseTime", "Null3")
  all.data$participantID <- p2
  all.out <- rbind(all.out, all.data)
}

## Now remove all NULL columns
null.col <- grep("Null", names(all.out))
all.out <- all.out[,-null.col]


## Now go through each of the ids and turn each item from long to wide data
id.vals <- names(table(all.out$participantID))
picture.vals <- names(table(all.out$PictureDisplayed))
all.out.wide <- NULL
real.all.wide <- NULL
question.vals <- names(table(all.out.wide$PictureDisplayed))
pb <- progress_bar$new(total = length(id.vals))
for(i in id.vals){
  # Isolate the data of interest
  data.to.use <- all.out[which(all.out$participantID==i),]
  ## Now go through an isolate the questions of interest
  # First create the participant specific output
  participant.wide <- NULL
  index <- 1
  irt.tmp <- i
  for(p in picture.vals){
    pic.data.to.use <- data.to.use[which(data.to.use$PictureDisplayed==p),]
    ## Now check to see if our dim is greater than 0
    if(dim(pic.data.to.use)[1]>0){
      ## Add a order column
      pic.data.to.use$order <- rank(as.numeric(as.character(pic.data.to.use$timeGiven)))
      ## Now reshape using the order as the time val
      resp.given <- dcast(pic.data.to.use, formula = PictureDisplayed + participantID ~ order, value.var = c('ResponseGiven'))
      resp.time <- dcast(pic.data.to.use, formula = PictureDisplayed + participantID ~ order, value.var = 'responseTime')
      # Now merge these
      row.out <- merge(resp.time, resp.given, by=c("PictureDisplayed", "participantID"), suffixes = c("_responseTime", "_responseGiven"))
      # Now make sure the out row is 1x10
      if(!identical(dim(row.out), c(1, 10))){
        #print("short")
        ## Find which column names are not in here
        all.names <- c("PictureDisplayed", "participantID", "1_responseGiven", "2_responseGiven", "3_responseGiven", "4_responseGiven", "1_responseTime", "2_responseTime", "3_responseTime", "4_responseTime")
        names.noin <- which(all.names %in% names(row.out)=="FALSE")
        ## Now create the noin names columns
        for(z in names.noin){
          row.out[all.names[z]] <- NA
        }
      }
      ## Now merge
      if(index ==1){
        participant.wide <- rbind(participant.wide, row.out)
      }else{
        participant.wide <- merge(row.out, participant.wide, by=c("PictureDisplayed", "participantID"), all=T) 
      }
      ## Now prepare the response given IRT data
      row.out.irt <- row.out[, c("PictureDisplayed", "participantID", "1_responseGiven", "2_responseGiven", "3_responseGiven", "4_responseGiven")]
      colnames(row.out.irt)[3:6] <- paste(colnames(row.out.irt)[3:6], row.out.irt$PictureDisplayed, sep='_')
      row.out.irt <- row.out.irt[,3:6]
      irt.tmp <- unlist(as.array(c(irt.tmp, row.out.irt)))
    }
  }
  all.out.wide <- rbind(all.out.wide, participant.wide)
  real.all.wide <- rbind(real.all.wide, irt.tmp)
  pb$tick()
}

## Now I want to go through and calculate internal consitency per person
output.person.con <- NULL
for(i in id.vals){
  ## Isolate the individual
  data.to.use <- all.out.wide[which(all.out.wide$participantID==i),]
  ## Now convert the data to binary responses
  columns.to.loop <- grep("Given", names(data.to.use))
  data.to.calc <- matrix(NA, nrow = dim(data.to.use)[1], ncol = 4)
  ## Now loop through
  index <- 1
  for(q in columns.to.loop){
   c.index <- which(data.to.use[,q]==10)
   data.to.calc[c.index,index] <- 0
   c.index <- which(data.to.use[,q]==18)
   data.to.calc[c.index,index] <- 0
   index <- index + 1
  }
  index <- 1
  for(q in columns.to.loop){
    c.index <- which(data.to.use[,q]==12)
    data.to.calc[c.index,index] <- 1
    c.index <- which(data.to.use[,q]==24)
    data.to.calc[c.index,index] <- 1
    index <- index + 1
  }
  ## Now estimate the cronbach's alpha for each emotion
  alpha.value           <- try(psych::alpha(data.to.calc, check.keys = T)$total[1])
  unh.alpha             <- try(psych::alpha(data.to.calc[grep("^U", data.to.use$PictureDisplayed),], check.keys = T, use='complete.obs', impute = 'mean')$total[1])
  cry.alpha             <- try(psych::alpha(data.to.calc[grep("^C", data.to.use$PictureDisplayed),], check.keys = T, use='complete.obs', impute = 'mean')$total[1])
  hap.alpha             <- try(psych::alpha(data.to.calc[grep("^H", data.to.use$PictureDisplayed),], check.keys = T, use='complete.obs', impute = 'mean')$total[1])
  neu.alpha             <- try(psych::alpha(data.to.calc[grep("^N", data.to.use$PictureDisplayed),], check.keys = T, use='complete.obs', impute = 'mean')$total[1])

  ## Now check to see if it is a numeric, and if not add a flag
  alpha.value <- ifelse(is.na(as.numeric(alpha.value)), NA, as.numeric(alpha.value))
  unh.alpha <- ifelse(is.na(as.numeric(unh.alpha)), NA, as.numeric(unh.alpha))
  cry.alpha <- ifelse(is.na(as.numeric(cry.alpha)), NA, as.numeric(cry.alpha))
  hap.alpha <- ifelse(is.na(as.numeric(hap.alpha)), NA, as.numeric(hap.alpha))
  neu.alpha <- ifelse(is.na(as.numeric(neu.alpha)), NA, as.numeric(neu.alpha))
  
  # Now prepare the output row
  out.row <- cbind(i, alpha.value[1], unh.alpha[1], cry.alpha[1], hap.alpha[1], neu.alpha[1])
  colnames(out.row) <- c("record_id", "overall", "unhappy", "crying", "happy", "neutral")
  output.person.con <- rbind(output.person.con, out.row)
  colnames(output.person.con) <- c("record_id", "overall", "unhappy", "crying", "happy", "neutral")
}
## It is kind of interesting -- the emotion specific values are very different than the overall???
## I wonder why this is??
## Also looks like they don't average out - the emtoin specific are consistently lower than the overall
output.person.con[,2:6] <- apply(output.person.con[,2:6], c(1,2), as.numeric)

## Now go through and clauclate the item veraiability
question.vals <- names(table(all.out.wide$PictureDisplayed))
## Change the items named: "UX,64" --> "UX,32"
all.out.wide$PictureDisplayed[which(all.out.wide$PictureDisplayed=="UX,64")] <- "UX,32"
question.vals <- names(table(all.out.wide$PictureDisplayed))[1:96]
output.con <- NULL
#for(i in question.vals){
#  data.to.use <- all.out.wide[which(all.out.wide$PictureDisplayed==i),]
#  #print(dim(data.to.use))
#  alpha.value <- psych::alpha(apply(data.to.use[,c("1_responseGiven","2_responseGiven","3_responseGiven","4_responseGiven")], c(1,2), function(x) as.numeric(as.character(x))))
#  alpha.value <- alpha.value$total[1]
#  out.row <- cbind(i, alpha.value)
#  output.con <- rbind(output.con, out.row)
#}


### Now run the IRT? and see if it'll work
#real.all.wide[real.all.wide=="10"] <- 0
#real.all.wide[real.all.wide=="18"] <- 0
#real.all.wide[real.all.wide=="12"] <- 1
#real.all.wide[real.all.wide=="24"] <- 1
#real.all.wide[,2:385] <- apply(real.all.wide[,2:385], c(1,2), function(x) as.numeric(as.character(x)))
#for.irt <- apply(real.all.wide[,2:385], c(1,2), function(x) as.numeric(as.character(x)))
#mod.1 <- mirt(for.irt, 1, IRTpars=T)
## Now print the coef
#out.mod.1 <- coef(mod.1)
#
## Now we are going to go rhtough and grab the difficulty paramter for each item over time


## Now go through and use the long version to run the IRT to get more stable difficulty estimates
# I am going to collapse the real.all.wide and stack them via the order of presentation
real.all.wide <- as.data.frame(real.all.wide)
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
for.irt[,2:97] <- apply(for.irt[,2:97], c(1,2), function(x) as.numeric(as.character(x)))

## Now remove any fully empty rows
index <- which(apply(for.irt[,2:97], 1, function(x) sum(is.na(x)))==96)
for.irt <- for.irt[-index,]

## I am going to go through and see which response is provided more frequently for the crying faces - this will be the 1 variable or the yes care flag
# This array gives us the values of yes response codes
yes.vals <- apply(for.irt[,c(2:25)], 1, function(x) names(which(table(x)==max(table(x)))))
no.vals <- apply(for.irt[,c(2:25)], 1, function(x) names(which(table(x)==min(table(x)))))
for(i in 1:length(yes.vals)){
  ## Go thorugh each row and change the yes to the yes val and the no to the no vals
  yes.val <- yes.vals[i]
  if(length(yes.val[[1]])>1){yes.val <- yes.val[[1]][1]}
  no.val <- no.vals[i]
  if(length(no.val[[1]])>1){no.val <- no.val[[1]][2]}
  yes.index <- which(for.irt[i, 2:97]==yes.val) +1
  no.index <-  which(for.irt[i, 2:97]==no.val)+1
  for.irt[i,yes.index] <- 0
  for.irt[i,no.index] <- 1
  ## Now print the dim of the responses
  print(dim(table(t(for.irt[i,2:97]))))
}

## Now convert the stragglers.. not sure why I have them
for.irt[,2:97] <- apply(for.irt[,2:97], 2, function(x) ifelse(x > 1, 0, x))
for.irt[,2:97] <- binary.flip(for.irt[,2:97])
## Calc ICC


## Now run IRT
mod.2 <- mirt(for.irt[,c(2:97)], 1, IRTpars=T)
mod.2.unhappy <- mirt(for.irt[,c(74:97)], 1, IRTpars=T)
mod.2.crying <-  mirt(for.irt[,c(2:25)], 1, IRTpars=T)
mod.2.neutral <-  mirt(for.irt[,c(50:73)], 1, IRTpars=T)
mod.2.happy <-  mirt(for.irt[,c(26:49)], 1, IRTpars=T)

## Now run a 2 param model
mod.3 <- mirt(for.irt[,c(2:97)], 2, IRTpars=T)

## Now calc some ICC for the factor scores
icc.data <- cbind(as.character(for.irt$V1), fscores(mod.3))
# Now break it down into 1 - 4 episodes
icc.data.1 <- as.data.frame(icc.data[1:61,])
icc.data.2 <- as.data.frame(icc.data[62:122,])
icc.data.3 <- as.data.frame(icc.data[123:184,])
icc.data.4 <- as.data.frame(icc.data[185:244,])
## Now merge em
icc.data <- merge(icc.data.1, icc.data.2, by='V1')
icc.data <- merge(icc.data, icc.data.3, by="V1")
icc.data <- merge(icc.data, icc.data.4, by="V1")
colnames(icc.data) <- c("record_id", 'fOnePredOne', "fTwoPredOne", "fOnePredTwo", "fTwoPredTwo", "fOnePredThree", "fTwoPredThree", "fOnePredFour", "fTwoPredFour")
icc.data[,2:9] <- apply(icc.data[,2:9], 2, function(x) as.numeric(as.character(x)))

## Now do the ICC for factor one
psych::ICC(icc.data[,c(2,4,6,8)]) # pretty decent!
psych::ICC(icc.data[,c(3,5,7,9)]) # less decent but acceptable



vals <- coef(mod.2)
output=NULL
for(i in 1:96) {
  output=rbind(output,
               matrix(vals[[i]],ncol=4,byrow=TRUE))
  rownames(output)[i] <- names(vals[i])
}
## Now get the difficulty values
output <- cbind(output, rownames(output))
output <- data.frame(output)
colnames(output)[2] <- "Difficulty"
output$Emotion <- "Happy"
output$Emotion[grep("_U", output$X5)] <- "Unhappy"
output$Emotion[grep("_C", output$X5)] <- "Cry"
output$Emotion[grep("_N", output$X5)] <- "Neutral"

## Now plot the IRT characteristics
plot(mod.2, type='trace', which.items=c(1:96))
plot(mod.2, type='infotrace', which.items=c(1:96))
plot(mod.2, type='info')
plot(mod.2, type='SE')
plot(mod.2.neutral, type='trace', which.items=c(1:24))
plot(mod.2.happy, type='trace', which.items=c(1:24))
plot(mod.2.crying, type='trace', which.items=c(1:24))
plot(mod.2.unhappy, type='trace', which.items=c(1:24))

## Now explore the dimenionslaity of these models following methods of:
# https://doi.org/10.3389/feduc.2019.00045


## Write the data
write.csv(output, "itemDifficultyCare.csv", quote=F, row.names=F)
## Now plot these
output$Difficulty <- as.numeric(as.character(output$Difficulty))
tmp.plot <- ggplot(output, aes(x=Difficulty)) +
  geom_histogram() +
  facet_grid(Emotion ~.) +
  theme_light()

colnames(output)[1] <- "Discrimination"
output$Discrimination <- as.numeric(as.character(output$Discrimination))
tmp.plot2 <- ggplot(output, aes(x=Discrimination)) +
  geom_histogram() +
  facet_grid(Emotion ~.) +
  theme_light()

## Now plot the consistency
output.person.con <- as.data.frame(output.person.con)
output.person.con[,2:6] <- apply(output.person.con[,2:6], c(1,2), function(x) as.numeric(as.character(x)))
tmp.plot3 <- ggplot(output.person.con, aes(x=overall)) +
  geom_histogram() +
  theme_light()

## Now plot emotion specific ish
output.person.con.2 <- melt(output.person.con, id.vars = 'record_id')
tmp.plot4 <- ggplot(output.person.con.2, aes(x=abs(value))) +
  geom_histogram() +
  facet_grid(variable ~.) +
  theme_light()

## Now write the data
write.csv(output.person.con, "./consistencyDataCARE.csv", quote=F, row.names=F)
write.csv(icc.data, "./careFactorScores.csv", quote=F, row.names=F)
