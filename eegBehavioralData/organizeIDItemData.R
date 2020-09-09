## Test script to organize data for IRT analysis using care task faces
source("~/adroseHelperScripts/R/afgrHelpFunc.R")
install_load("reshape2", "progress", "mirt", "psych", "ggplot2")

## Run the bash script to organize the files as desired
system("/home/arosen/Documents/bbmcPersonal/eegBehavioralData/organizeIDItemData.sh")

## This is going to have to be done in a loop
all.files.pic <- system("ls ~/Documents/bbmcPersonal/eegBehavioralData/careFace/idMod/pic*", intern = T)
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
  file.1 <- paste("~/Documents/bbmcPersonal/eegBehavioralData/careFace/idMod/pic", p2,".csv", sep='')
  file.2 <- paste("~/Documents/bbmcPersonal/eegBehavioralData/careFace/idMod/resp", p2,".csv", sep='')
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
out.best.guesses <- NULL
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
  ## Now find the best estimate for the emotion values
  # First collapse the emotions into the four possible emotions
  data.to.use$emotion <- "unhappy"
  data.to.use$emotion[grep("^N", data.to.use$PictureDisplayed)] <- "neutral"
  data.to.use$emotion[grep("^H", data.to.use$PictureDisplayed)] <- "happy"
  data.to.use$emotion[grep("^C", data.to.use$PictureDisplayed)] <- "crying"
  ## now create an output with the best guess
  out.best.guess <- matrix(NA, nrow=1, ncol = 5)
  colnames(out.best.guess) <- c("record_id", rownames(table(data.to.use$emotion, as.character(data.to.use$ResponseGiven))))
  response.patterns <- table(data.to.use$emotion, data.to.use$ResponseGiven)
  ## Now go through each row and find the maximum value -- and return a flag if double!
  for(W in row.names(response.patterns)){
    row.to.check <- response.patterns[W,]
    val.max <- names(which(row.to.check == max(row.to.check)))
    if(length(val.max)==1){
      out.best.guess[,W] <- val.max
    }
    if(length(val.max)>1){
      out.best.guess[,W] <- "Tie"
    }
  }
  out.best.guess[,"record_id"] <- i
  out.best.guesses <- rbind(out.best.guesses, out.best.guess)
  pb$tick()
}

## It looks like this method failed and there is far too much confusion -- I am going to grab the left vs right estimates from the
## createConfusionMatrix.R script
## Here if the value is right the emotions will be:
# happy == 18 ; neutral == 20 ; unhappy == 22 ; crying == 24
## if emtoins are left then:
# happy == 24 ; unhappy == 20 ; neutral == 22 ; crying == 18 
direction.vals <- read.csv('./leftRightValues.csv')
direction.vals$record_id <- strSplitMatrixReturn(direction.vals$V1, "_")[,1]
direction.vals <- merge(direction.vals, out.best.guesses, all=T)
### Looks like I am missing data for 118-81 and 118-103 from the original 
### Looks like I am issing data for 118-62 from the new

## Now create new datasets using the best guesses as the item answers
## Now go through and change to the correct emotion triggers
direction.vals$crying <- 24
direction.vals$crying[direction.vals$Guess=="Left"] <- 18
direction.vals$unhappy <- 22
direction.vals$unhappy[direction.vals$Guess=="Left"] <- 20
direction.vals$neutral <- 20
direction.vals$neutral[direction.vals$Guess=="Left"] <- 22
direction.vals$happy <- 18
direction.vals$happy[direction.vals$Guess=="Left"] <- 24

## Now change the 


## Go through and make the propoer changes
# Adding a tmp id vals as we are missing some data somewhere along this line
id.vals <- direction.vals$record_id[complete.cases(direction.vals)]
all.out.wide2 <- NULL
for(i in id.vals){
  ## Isolate the individual
  data.to.use <- all.out.wide[which(all.out.wide$participantID==i),]
  if(dim(data.to.use)[1]==0){next}
  ## answer values
  mod.vals <- direction.vals[which(direction.vals$record_id==i),]
  ## Now go through each of the emotions and change each of the response values appropriatly
  orig.vals <- c(18, 20, 22, 24)
  for(W in orig.vals){
    # Find the emotion to use
    new.val <- colnames(mod.vals)[which(mod.vals==W)]
    # Now find all of the indices
    orig.indices <- which(data.to.use==W)
    if(length(orig.indices)==0){
      print(paste("Participtant:", i, "Missing:", W,";", new.val, sep = " "))
      next
    }
    orig.indices <- matrix(rc.ind(data.to.use, orig.indices), byrow=F, ncol=2)
    ## Now loop through the indices
    for(L in 1:dim(orig.indices)[1]){
      data.to.use[orig.indices[L,1], orig.indices[L,2]] <- new.val
    }
  }
  ## Now attach this to the new output
  all.out.wide2 <- rbind(all.out.wide2, data.to.use)
}
## Now I can check to see if both the reaction time & response are NA --
## If they are not then I can go back and add the response depending on missing response values
## Thats going to be for a later time though


## Make the variables factors so we can compute the alpha for them
all.out.wide2[,7] <- factor(all.out.wide2[,7], c("crying", "unhappy", "neutral", "happy"))
all.out.wide2[,8] <- factor(all.out.wide2[,8], c("crying", "unhappy", "neutral", "happy"))
all.out.wide2[,9] <- factor(all.out.wide2[,9], c("crying", "unhappy", "neutral", "happy"))
all.out.wide2[,10] <- factor(all.out.wide2[,10], c("crying", "unhappy", "neutral", "happy"))

## Now I want to go through and calculate internal consistency per person
output.person.con <- NULL
for(i in id.vals){
  ## Isolate the individual
  data.to.use <- all.out.wide2[which(all.out.wide2$participantID==i),]
  if(dim(data.to.use)[1]==0){next}
  data.to.calc <- matrix(NA, nrow = dim(data.to.use)[1], ncol = 4)
  data.to.calc[,1] <- as.numeric(data.to.use[,7])
  data.to.calc[,2] <- as.numeric(data.to.use[,8])
  data.to.calc[,3] <- as.numeric(data.to.use[,9])
  data.to.calc[,4] <- as.numeric(data.to.use[,10])
  
  ## Now estimate the cronbach's alpha
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

## Now go through and calculate the item variability
question.vals <- names(table(all.out.wide$PictureDisplayed))
## Change the items named: "UX,64" --> "UX,32"
all.out.wide$PictureDisplayed[which(all.out.wide$PictureDisplayed=="UX,64")] <- "UX,32"
question.vals <- names(table(all.out.wide$PictureDisplayed))[1:96]
output.con <- NULL
for(i in question.vals){
  data.to.use <- all.out.wide[which(all.out.wide$PictureDisplayed==i),]
  #print(dim(data.to.use))
  alpha.value <- psych::alpha(apply(data.to.use[,c("1_responseGiven","2_responseGiven","3_responseGiven","4_responseGiven")], c(1,2), function(x) as.numeric(as.character(x))))
  alpha.value <- alpha.value$total[1]
  out.row <- cbind(i, alpha.value)
  output.con <- rbind(output.con, out.row)
}


## Now try to train a LMER mod for response time ~ picture displayed
# First load the IRT data
item.diff <- read.csv("itemDifficultyCare.csv", row.names = NULL)
item.diff$itemName <- paste(strSplitMatrixReturn(item.diff$X4, "_")[,2], item.diff$X5, sep=',')
# Merge em
for.lmer <- merge(all.out, item.diff, by.x = "PictureDisplayed", by.y="itemName", all=T)
for.lmer$row.names <- as.numeric(for.lmer$row.names)
for.lmer$responseTime <- as.numeric(as.character(for.lmer$responseTime))
tmp <- lmerTest::lmer(responseTime ~   Emotion+ X1+(1|PictureDisplayed /participantID), data=for.lmer)
write.csv(for.lmer, "responseTimeQuestion.csv", quote=F, row.names=F)


## Now check for endurance effects
tmp.2 <- lmerTest::lmer(responseTime~ timeShown + (1|participantID), data=for.lmer)


## Now plot some of these values
output.person.con <- as.data.frame(output.person.con)
output.person.con.2 <- reshape2::melt(output.person.con, id.vars = 'record_id')
output.person.con.2$value <- as.numeric(as.character(output.person.con.2$value))
tmp.plot4 <- ggplot(output.person.con.2, aes(x=abs(value))) +
  geom_histogram() +
  facet_grid(variable ~.) +
  theme_light()
