# ## Clean everything
rm(list=ls())
 
# ## Test script to organize data for IRT analysis using ID task faces
source("~/adroseHelperScripts/R/afgrHelpFunc.R")
install_load("reshape2", "progress", "mirt", "psych", "ggplot2", "polycor", "corrplot", "ggalluvial")

## Run the bash script to organize the files as desired
system("/home/arosen/Documents/bbmcPersonal/eegBehavioralData/organizeIDItemData.sh")

## Declare any functions
binary.flip <- function (x) {
  x * -1 + 1
}

# ## This is going to have to be done in a loop
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
real.all.wide.rt <- NULL
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
  irt.tmp2 <- i
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
      ## Now do the same on RT
      row.out.irt2 <- row.out[,c("PictureDisplayed", "participantID", "1_responseTime", "2_responseTime", "3_responseTime", "4_responseTime")]
      colnames(row.out.irt2)[3:6] <- paste(colnames(row.out.irt2)[3:6], row.out.irt2$PictureDisplayed, sep='_')
      row.out.irt2 <- row.out.irt2[,3:6]
      irt.tmp2 <- unlist(as.array(c(irt.tmp2, row.out.irt2)))
    }
  }
  all.out.wide <- rbind(all.out.wide, participant.wide)
  real.all.wide <- rbind(real.all.wide, irt.tmp)
  real.all.wide.rt <- rbind(real.all.wide.rt, irt.tmp2)

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
### Looks like I am missing data for 118-62 from the new

## Now create new datasets using the best guesses as the item answers
## Now go through and change to the correct emotion triggers
## Manually flip the 117-114 patterns!
direction.vals$Guess[which(direction.vals$record_id=="117-114")] <- "Left"
## Also need to manually flip: 118-286; I am not sure why this individual failed?
direction.vals$Guess[which(direction.vals$record_id=="118-286")] <- "Right"
direction.vals$Guess[which(direction.vals$record_id=="118-300")] <- "Right"
direction.vals$crying <- 24
direction.vals$crying[direction.vals$Guess=="Right"] <- 18
direction.vals$unhappy <- 22
direction.vals$unhappy[direction.vals$Guess=="Right"] <- 20
direction.vals$neutral <- 20
direction.vals$neutral[direction.vals$Guess=="Right"] <- 22
direction.vals$happy <- 18
direction.vals$happy[direction.vals$Guess=="Right"] <- 24

## Go through and make the proper changes
# Adding a tmp id vals as we are missing some data somewhere along this line
id.vals <- direction.vals$record_id[complete.cases(direction.vals)]
all.out.wide2 <- NULL
for(i in id.vals){
  ## Isolate the individual
  data.to.use <- all.out.wide[which(all.out.wide$participantID==i),]
  if(dim(data.to.use)[1]==0){next}
  ## answer values
  mod.vals <- direction.vals[which(direction.vals$record_id==i),]
  ## Now go through each of the emotions and change each of the response values appropriately
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

## This code was added 2020-10-06 because ADROSE found some issues with the correct vs incorrect coding!!!
## Now go through and double check the direction vals
all.out.wide2$emotion <- "crying"
all.out.wide2$emotion[grep("^H", all.out.wide2$PictureDisplayed)] <- "happy"
all.out.wide2$emotion[grep("^U", all.out.wide2$PictureDisplayed)] <- "unhappy"
all.out.wide2$emotion[grep("^N", all.out.wide2$PictureDisplayed)] <- "neutral"

## Now double check the scores by participant
crying.response <- table(all.out.wide2$emotion, all.out.wide2$participantID, all.out.wide2$"1_responseGiven")[,,1]
happy.response <- table(all.out.wide2$emotion, all.out.wide2$participantID, all.out.wide2$"1_responseGiven")[,,2]
neutral.response <- table(all.out.wide2$emotion, all.out.wide2$participantID, all.out.wide2$"1_responseGiven")[,,3]
unhappy.response <- table(all.out.wide2$emotion, all.out.wide2$participantID, all.out.wide2$"1_responseGiven")[,,4]

## Now do the real.all.wide
real.all.wide2 <- NULL
real.all.wide <- as.data.frame(real.all.wide)
for(i in id.vals){
  ## Isolate the individual
  data.to.use <- real.all.wide[which(real.all.wide[,1]==i),]
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
    ## Now change the values to the new value
    data.to.use[,orig.indices] <- new.val
  }
  ## Now attach this to the new output
  real.all.wide2 <- rbind(real.all.wide2, data.to.use)
}

## Now I can check to see if both the reaction time & response are NA --
## If they are not then I can go back and add the response depending on missing response values
## Thats going to be for a later time though

## Make the variables factors so we can compute the alpha for them
all.out.wide2[,7] <- factor(all.out.wide2[,7], c("crying", "unhappy", "neutral", "happy"))
all.out.wide2[,8] <- factor(all.out.wide2[,8], c("crying", "unhappy", "neutral", "happy"))
## Add in a checkpoint here!!!
all.out.wide2[,9] <- factor(all.out.wide2[,9], c("crying", "unhappy", "neutral", "happy"))
all.out.wide2[,10] <- factor(all.out.wide2[,10], c("crying", "unhappy", "neutral", "happy"))

save.image(file = "./eegBehavioralData/IRTIDData.RData")
load("./eegBehavioralData/IRTIDData.RData")


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

# ## Now try to train a LMER mod for response time ~ picture displayed
# # First load the IRT data
# item.diff <- read.csv("./eegBehavioralData/itemDifficultyID.csv", row.names = NULL)
# item.diff$itemName <- paste(strSplitMatrixReturn(item.diff$X4, "_")[,2], item.diff$X5, sep=',')
# # Merge em
# for.lmer <- merge(all.out, item.diff, by.x = "PictureDisplayed", by.y="itemName", all=T)
# for.lmer$row.names <- as.numeric(for.lmer$row.names)
# for.lmer$responseTime <- as.numeric(as.character(for.lmer$responseTime))
# tmp <- lmerTest::lmer(responseTime ~   Emotion+ X1+(1|PictureDisplayed /participantID), data=for.lmer)
# write.csv(for.lmer, "responseTimeQuestion.csv", quote=F, row.names=F)


# ## Now check for endurance effects
# tmp.2 <- lmerTest::lmer(responseTime~ timeShown + (1|participantID), data=for.lmer)

## Now play with the IRT models here
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
## Now go through and change all of these characters into numeric
# Cry == 1; Unhappy == 2; Neutral ==3; Happy == 4
new.val <- 1
for.irt2 <- for.irt
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
plot(mod.1, type='trace', which.items=c(1:96))
plot(mod.1, type='infotrace', which.items=c(1:96))

## Now do a binary () correct vs incorrect) IRT model
for.irt3 <- matrix(NA, ncol = 96, nrow = dim(for.irt)[1])
# I need to go through and find the correct string for each column and then identify the correct responses within each cell
for(i in 2:97){
  string.index <- for.irt[,i]
  irt_three_col <- i -1
  ## First grab the correct index
  # Grab the character from the colname
  char.val <- tolower(substr(strSplitMatrixReturn(colnames(for.irt)[i], "_")[,2][1], 1, 1))
  print(char.val)
  if(char.val == 'c') cor.index <- which(string.index=="crying")
  if(char.val == 'n') cor.index <- which(string.index=="neutral")
  if(char.val == 'u') cor.index <- which(string.index=="unhappy")
  if(char.val == 'h') cor.index <- which(string.index=="happy")
  ## Now change the values
  for.irt3[which(!is.na(for.irt[,i])),irt_three_col] <- 0
  for.irt3[cor.index,irt_three_col] <- 1
}
for.irt3 <- as.data.frame(for.irt3)
## Now run the IRT with the binary data
#mod.2 <- mirt(for.irt3, 1, SE=T)
#extract.mirt(mod.2, "BIC") ## 12803.46
mod.12 <- mirt(for.irt3, 1, SE=T, itemtype = "2PL")
plot(mod.12, type='trace', which.items=c(49:96))
mod.22 <- mirt(for.irt3, 2, SE=T, itemtype = "2PL")

## Now prepare for the item level covariates analysis
for.irt3.cv <- for.irt3
colnames(for.irt3.cv) <- strSplitMatrixReturn(colnames(for.irt)[2:97], "_")[,2]
for.irt3.cv$record_id <- for.irt$V1
for.irt3.cv <- melt(for.irt3.cv, id.vars = "record_id")
## Now attach item level covariates
item.covs <- read.csv("./eegBehavioralData/itemDemographics.csv", na.strings = ".")
## Now prep the item names
item.covs$variable <- paste(item.covs$itemCode, item.covs$Val2, sep=',')

## Now load the participants demographics
load("./fname.gz")
ds_eegPred <- out.data[[9]]
for.irt3.cv <- merge(for.irt3.cv, ds_eegPred, by="record_id")
## Now merge the item covariates
for.irt3.cv <- merge(for.irt3.cv, item.covs, by="variable")
## Now create a race agreement column
for.irt3.cv$raceAgree <- 0
for.irt3.cv$raceAgree[which(for.irt3.cv$race_ethnicity_h2 == "White" & for.irt3.cv$Race == "CA")] <- 1
for.irt3.cv$raceAgree[which(for.irt3.cv$race_ethnicity_h2 == "Black" & for.irt3.cv$Race == "AA")] <- 1
for.irt3.cv$raceAgree[which(for.irt3.cv$race_ethnicity_h2 == "Latino" & for.irt3.cv$Race == "LA")] <- 1
for.irt3.cv$raceAgree[which(for.irt3.cv$race_ethnicity_h2 == "American Indian" & for.irt3.cv$Race == "NA")] <- 1
## Now train a model
item.cov.explore <- lme4::glmer(value ~ (Gender + Race + Emotion)^2 + (1|record_id), data=for.irt3.cv, family = binomial)
visreg::visreg(item.cov.explore, "Gender", by="Emotion", overlay=F)
visreg::visreg(item.cov.explore, "Race", by="Emotion", overlay=F)
for.irt3.cv$raceAgree <- factor(for.irt3.cv$raceAgree)
item.cov.explore2 <- lme4::glmer(value ~ (Gender + raceAgree + Emotion)^2 + (1|record_id), data=for.irt3.cv, family = binomial)
visreg::visreg(item.cov.explore2, "Gender", by="raceAgree")
#item.cov.explore <- glm(value ~ Gender + Race * raceAgree + variable, data=for.irt3.cv, family = binomial)

## Now explore variance across trials within a participant
# First grab the individual ID values
for.irt3.cv <- for.irt3
colnames(for.irt3.cv) <- strSplitMatrixReturn(colnames(for.irt)[2:97], "_")[,2]
for.irt3.cv$record_id <- for.irt$V1
# Now grab the correct values within an individual across time
# Start with overall
for.irt3.cv$allCor <- rowSums(for.irt3.cv[,1:96], na.rm = T)
for.irt3.cv$allPer <- for.irt3.cv$allCor / rowSums(!is.na(for.irt3.cv[,1:96]))
# Now do the emotion specific values
for.irt3.cv$cryCor <- rowSums(for.irt3.cv[,1:24], na.rm = T)
for.irt3.cv$cryPer <- for.irt3.cv$cryCor / rowSums(!is.na(for.irt3.cv[,1:24]))
for.irt3.cv$hapCor <- rowSums(for.irt3.cv[,25:48], na.rm = T)
for.irt3.cv$hapPer <- for.irt3.cv$hapCor / rowSums(!is.na(for.irt3.cv[,25:48]))
for.irt3.cv$neuCor <- rowSums(for.irt3.cv[,49:72], na.rm = T)
for.irt3.cv$neuPer <- for.irt3.cv$neuCor / rowSums(!is.na(for.irt3.cv[,49:72]))
for.irt3.cv$unhCor <- rowSums(for.irt3.cv[,73:96], na.rm = T)
for.irt3.cv$unhPer <- for.irt3.cv$unhCor / rowSums(!is.na(for.irt3.cv[,73:96]))

# Now add a marker for the round each response pattern is associated with rounds go from 1--> 4
for.irt3.cv$round <- 1
for.irt3.cv$round[62:122] <- 2
for.irt3.cv$round[123:183] <- 3
for.irt3.cv$round[184:244] <- 4

## Now plot these performances across time within each individual
library(tidyverse)
tmp.plot <- for.irt3.cv %>% ggplot(., aes(x=round, y=allPer, group=record_id)) +
  geom_point() +
  geom_line()
## Now do the same for each emotion
for.irt3.cv.tmp <- for.irt3.cv[,c("record_id", "round", "cryPer", "hapPer", "neuPer", "unhPer")]
for.irt3.cv.tmp <- melt(for.irt3.cv.tmp, id.vars = c("record_id", "round"))
tmp.plot2 <- for.irt3.cv.tmp %>% ggplot(., aes(y=record_id, x=value, group=record_id, color=variable)) +
  geom_point() +
  #geom_line() +
  theme_bw()
  #facet_wrap( ~ variable)

## Now try to model practice effects?
mod.practice <- lmerTest::lmer(value ~ variable * round + (1|record_id), data=for.irt3.cv.tmp)

## Now fit 2 1-factor models using the high intensity and low intensity data seperate
mod.11 <- mirt(for.irt3[,1:48], 1, SE=T, itemtype='2PL')
mod.13 <- mirt(for.irt3[,49:96], 1, SE=T, itemtype='2PL')
mod.14 <- mirt(for.irt3[,49:96], 2, SE=T, itemtype='2PL')

plot(mod.11, type='trace', which.items=c(1:48))
plot(mod.13, type='trace', which.items=c(1:48))

mod.2.unhappy <- mirt(for.irt3[,c(73:96)], 1, IRTpars=T, itemtype = "2PL")
mod.2.crying <-  mirt(for.irt3[,c(1:24)], 1, IRTpars=T, itemtype = "2PL")
mod.2.neutral <-  mirt(for.irt3[,c(49:72)], 1, IRTpars=T, itemtype = "2PL")
mod.2.happy <-  mirt(for.irt3[,c(25:48)], 1, IRTpars=T, itemtype = "2PL")


extract.mirt(mod.2, "BIC") ## 12608.97
#mod.2 <- mirt(for.irt3, 3, SE=T)
#extract.mirt(mod.2, "BIC") ## 12955.42

## Now prep the RT data
for.irt <- rbind(real.all.wide.1, real.all.wide.2, real.all.wide.3, real.all.wide.4)
for.irt[,2:97] <- apply(for.irt[,2:97], c(1,2), function(x) as.numeric(as.character(x)))

## Now do the RT
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



## Now try a FA approach
for.irt.mixed <- for.irt3
for.irt3 <- as.data.frame(apply(for.irt3, 2, factor))
colnames(for.irt3) <- colnames(for.irt)[2:97]
cor.mat <- polycor::hetcor(for.irt3, use=c("pairwise.complete.obs"))
tmp.dat <- cor.mat$cor
colnames(tmp.dat) <- gsub(colnames(tmp.dat), pattern = "responseGiven_", replacement = "")
rownames(tmp.dat) <- gsub(rownames(tmp.dat), pattern = "responseGiven_", replacement = "")
## Now correct the na values
na.index <- matrix(rc.ind(tmp.dat, which(is.na(tmp.dat))), byrow=F, ncol=2)
for(i in 1:dim(na.index)[1]){
  x.val <- na.index[i,1]
  y.val <- na.index[i,2]
  ## Now get the variable names
  x.val.char <- rownames(tmp.dat)[x.val]
  y.val.char <- colnames(tmp.dat)[y.val]
  ## Now compute the correlation
  tmp.dat[x.val,y.val] <- polycor::polychor(for.irt3[,x.val], for.irt3[,y.val])
}
corrplot(tmp.dat)
fa.2 <- fa(r = tmp.dat, n.obs = nrow(for.irt3), rotate = "varimax", nfactors = 2) ## FA suggests a 2 factor solution is the most optimal

## Now run a fa on the rt data
fa.3 <- fa(r = for.irt.rt[,2:97], rotate='varimax', nfactors=6)

## Now run a FA with the rt data too
for.irt.rt <- for.irt.rt[-which(for.irt.rt[,1] %in% for.irt$V1 == FALSE),]
for.irt4 <- as.data.frame(cbind(for.irt3, for.irt.rt[,2:97]))
cor.mat <- polycor::hetcor(for.irt4)
fa.3 <- fa(r = cor.mat$cor, n.obs = nrow(for.irt4), rotate = "varimax", nfactors = 2) ## FA suggests a 2 factor solution is the most optimal

## Try a bifactor IRT solution here
mod.7 <- bfactor(data.matrix(for.irt3), c(rep(1, 48), rep(2, 48)))


#plot(mod.2, type='trace', which.items=c(1:96))
#plot(mod.2, type='infotrace', which.items=c(1:96))

## Now calc some ICC for the factor scores
icc.data <- cbind(as.character(for.irt$V1), fscores(mod.1))
# Now break it down into 1 - 4 episodes
icc.data.1 <- as.data.frame(icc.data[1:61,])
icc.data.2 <- as.data.frame(icc.data[62:122,])
icc.data.3 <- as.data.frame(icc.data[123:184,])
icc.data.4 <- as.data.frame(icc.data[185:244,])
## Now merge em
icc.data <- merge(icc.data.1, icc.data.2, by='V1')
icc.data <- merge(icc.data, icc.data.3, by="V1")
icc.data <- merge(icc.data, icc.data.4, by="V1")
#colnames(icc.data) <- c("record_id", 'fOnePredOne', "fTwoPredOne", "fOnePredTwo", "fTwoPredTwo", "fOnePredThree", "fTwoPredThree", "fOnePredFour", "fTwoPredFour")
icc.data[,2:5] <- apply(icc.data[,2:5], 2, function(x) as.numeric(as.character(x)))

## Now do the ICC for factor one
psych::ICC(icc.data[,c(2,3,4,5)]) # pretty decent!
#psych::ICC(icc.data[,c(3,5,7,9)]) # less decent but acceptable

## Now plot some of these values
output.person.con <- as.data.frame(output.person.con)
output.person.con.2 <- reshape2::melt(output.person.con, id.vars = 'record_id')
output.person.con.2$value <- as.numeric(as.character(output.person.con.2$value))
tmp.plot4 <- ggplot(output.person.con.2, aes(x=abs(value))) +
  geom_histogram() +
  facet_grid(variable ~.) +
  theme_light()

## Now write the data
write.csv(output.person.con, "./consistencyDataID.csv", quote=F, row.names=F)
write.csv(icc.data, "./idFactorScores.csv", quote=F, row.names=F)


## Now explore some question specific effects
#load covars
load("./fname.gz")
ds_eegPred <- out.data[[9]]

## First try to run a IRT w/ covariates from the MIRT package
## see help here: https://rdrr.io/cran/mirt/man/mixedmirt.html

# Prepare the data
library(tidyverse)
source('~/adroseHelperScripts/R/afgrHelpFunc.R')
for.irt3 <- for.irt.mixed
for.irt3$record_id <- for.irt2$V1
# Fix broken ID values
for.irt3 <- for.irt3 %>% mutate(record_id = as.character(record_id)) %>% mutate(record_id = if_else(record_id %in% "117-115", "118-64", record_id))
## Add factor scores to the for.irt3 data
for.irt3$mod2Crying <- fscores(mod.2.crying)
for.irt.mixed <- merge(for.irt3, ds_eegPred, all=F)
# iter.vals <- c("ace_score", "dose_hv_visit_count", names(for.irt.mixed)[grep("^P2",names(for.irt.mixed), perl=T)])
# for(i in iter.vals){
#   for.irt.mixed[!is.na(for.irt.mixed[,i]),i] <- range01(for.irt.mixed[!is.na(for.irt.mixed[,i]),i])+1
# }
## Now do some models with the factor scores
lmer.mod.crying <- lmerTest::lmer(mod2Crying ~ P2crymaxCarepost * P2crymaxCareLatpost + (1|record_id), data=for.irt.mixed)
lmer.mod.crying.ace <- lmerTest::lmer(mod2Crying ~ ace_score + (1|record_id), data=for.irt.mixed)
lmer.mod.crying.hv <- lmerTest::lmer(mod2Crying ~ dose_hv_visit_count + (1|record_id), data=for.irt.mixed)
lmer.mod.crying.hvBace <- lmerTest::lmer(mod2Crying ~ dose_hv_visit_count*ace_score + (1|record_id), data=for.irt.mixed)
## Now for s & g try to predict the product of amp & lat w/ ace * hv
for.irt.mixed.tmp <- for.irt.mixed
for.irt.mixed.tmp$productOutcome <- for.irt.mixed.tmp$P2crymaxCarepost * for.irt.mixed.tmp$P2crymaxCareLatpost
## Now limit it to one per individual because there will be no variability in record_id outcomes
for.irt.mixed.tmp[!duplicated(for.irt.mixed.tmp$record_id),]
for.irt.mixed.tmp <- for.irt.mixed.tmp[,c("record_id", "ace_score", "")]

lmer.mod.crying.prod <- lm(productOutcome ~ dose_hv_visit_count * ace_score, data=for.irt.mixed.tmp)


## Now look at these effects
visreg::visreg(lmer.mod.crying, "P2crymaxCarepost", by="P2crymaxCareLatpost", overlay=T, gg=T) + theme_bw()
visreg::visreg(lmer.mod.crying.hvBace, "dose_hv_visit_count", by="ace_score", overlay=T, gg=T) + theme_bw()
visreg::visreg(lmer.mod.crying.prod, "dose_hv_visit_count", by="ace_score", overlay=T, gg=T) + theme_bw()

## Prepare multiple cores
mirtCluster(spec = 4)

## Declare the model
# This model takes far too long to run; I will try at a later date; but I am goin got try an emotion specific model
mod.mixed.cry1 <- mixedmirt(for.irt.mixed[,2:25], covdata = for.irt.mixed[,c("record_id", "ace_score")], model=1, fixed=~ace_score,  random=~1|record_id,itemtype = "2PL", technical = list(removeEmptyRows=TRUE))
mod.mixed.cry2 <- mixedmirt(for.irt.mixed[,2:25], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count")], model=1, fixed=~dose_hv_visit_count,  random=~1|record_id,itemtype = "2PL")
# Dose is a highly sig effect!
mod.mixed.cry3 <- mixedmirt(for.irt.mixed[,2:25], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "ace_score")], model=1, fixed=~dose_hv_visit_count*ace_score, random=~1|record_id, itemtype = "2PL")
# Interaction is crazy sig!? --> not with the random effect for individual =/

## Now try some eeg outcome
mod.mixed.cry4 <- mixedmirt(for.irt.mixed[,2:25], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2crymaxCarefront")], model=1, fixed=~P2crymaxCarefront, random=~1|record_id,itemtype = "2PL")
mod.mixed.cry5 <- mixedmirt(for.irt.mixed[,2:25], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2crymaxCarepost")], model=1, fixed=~P2crymaxCarepost, random=~1|record_id,itemtype = "2PL")
mod.mixed.cry6 <- mixedmirt(for.irt.mixed[,2:25], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2crymaxCareLatfront")], model=1, fixed=~P2crymaxCareLatfront, random=~1|record_id,itemtype = "2PL")
mod.mixed.cry7 <- mixedmirt(for.irt.mixed[,2:25], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2crymaxCareLatpost")], model=1, fixed=~P2crymaxCareLatpost, random=~1|record_id,itemtype = "2PL")
mod.mixed.cry8 <- mixedmirt(for.irt.mixed[,2:25], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2crymaxCarefront", "P2crymaxCareLatfront")], model=1, fixed=~P2crymaxCarefront*P2crymaxCareLatfront, random=~1|record_id,itemtype = "2PL")
mod.mixed.cry8d <- mixedmirt(for.irt.mixed[,2:25], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2crymaxCarefront", "P2crymaxCareLatfront")], model=1, fixed=~P2crymaxCarefront*P2crymaxCareLatfront, random=~1|items,itemtype = "2PL")
# Interactions are very sig
mod.mixed.cry9 <- mixedmirt(for.irt.mixed[,2:25], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2crymaxCarepost", "P2crymaxCareLatpost")], model=1, fixed=~P2crymaxCarepost*P2crymaxCareLatpost, random=list(~1|record_id, ~1|items),itemtype = "2PL", return.design = FALSE)
mod.mixed.cry9D <- mixedmirt(for.irt.mixed[,2:25], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2crymaxCarepost", "P2crymaxCareLatpost")], model=1, fixed=~P2crymaxCarepost*P2crymaxCareLatpost, random=list(~1|record_id, ~1|items),itemtype = "2PL", return.design = TRUE)

## Now try to plot these difficulty values for the rear crying questions
## First thing I will have to identify the tertiles for the latency variable
mod.plot <- for.irt.mixed[complete.cases(for.irt.mixed$P2crymaxIDLatpost),]
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



## Now estimate the difficult params
mod.plot$step1 <- mod.plot$P2crymaxIDpost * -1.719
mod.plot$step2 <- mod.plot$P2crymaxIDLatpost * -0.016
mod.plot$step3 <- mod.plot$P2crymaxIDpost * mod.plot$P2crymaxIDLatpost
mod.plot$step4 <- mod.plot$step3 * 0.006
mod.plot$step5 <- mod.plot$step1 + mod.plot$step2 + mod.plot$step4 + 2.936
## Now plot a histogram of these vlaues
mod.plot <- mod.plot[complete.cases(mod.plot$P2happymaxIDfront),]
p1 <- ggplot(mod.plot, aes(x=step5, group=plot.vals, fill=plot.vals)) +
  geom_density(position = "dodge") +
  theme_bw()
mod.plot$Theta <- randef(mod.mixed.cry9)$Theta
p2 <- ggplot(mod.plot, aes(x=Theta, group=plot.vals, fill=plot.vals)) +
  geom_density(position = "dodge") +
  theme_bw()
multiplot(p1, p2, cols = 2)

## Try to plot some ICC values here
# First obtain all of the fixed effect values
vals.one <- mod2values(mod.mixed.cry9)
denomFixed <- function(intercept=2.935601570, P2crymaxIDpost=-1.719401835, P2crymaxIDLatpost=-0.015931629, int=0.005759492, subjLat=NULL, subjAmp=NULL, record_id=0, items=2.818588, theta=0, itemDiscrim=0){
  ## Return a vector of values
  val <-theta * itemDiscrim + intercept + P2crymaxIDpost * subjAmp + P2crymaxIDLatpost * subjLat + int*(subjLat*subjAmp) + record_id + items
  val <- val * -1
  val <- exp(val)
  val <- 1 + val
  return(val)
}

## Now in order to create the ICC for each item I need to loop through:
## Each ITEM random effect; each subjects AMP & LAT
## First prepare the output data
out.vals <- NULL
item.re <- randef(mod.mixed.cry9)$items
item.disc <- unlist(lapply(coef(mod.mixed.cry9, IRTpars=T, as.data.frame=TRUE), function(x)x[13]))[1:24]
subject.amplitude <- quantile(for.irt.mixed$P2crymaxCarepost, c(0, .75), na.rm = T)
subject.lat <- quantile(for.irt.mixed$P2crymaxCareLatpost, c(0, .75), na.rm = T)
all.perm <- expand.grid(subject.amplitude, subject.lat)
for(i in 1:length(item.re)){
  for(l in 1:dim(all.perm)[1]){
    ## Calc the prob of correct response across the theta board
    vals.tmp <- 1/denomFixed(theta = seq(-4,4,.01), record_id = 0, itemDiscrim = item.disc[i], items = item.re[i], subjLat = all.perm$Var2[l], subjAmp = all.perm$Var1[l])
    ## Now prep the output data
    out.mat <- matrix(NA, nrow=length(vals.tmp), ncol=5)
    out.mat[,1] <- seq(-4,4,.01)
    out.mat[,2] <- all.perm$Var2[l]
    out.mat[,3] <- all.perm$Var1[l]
    out.mat[,4] <- i
    out.mat[,5] <- vals.tmp
    ## Now combine it
    out.vals <- rbind(out.vals, out.mat)
  }
}
## Now plot these
out.vals <- as.data.frame(out.vals)
out.vals$V4 <- factor(out.vals$V4)
out.vals$brainPerm <- paste(out.vals$V2, out.vals$V3, sep='_')
## Now plot it
ggplot(out.vals, aes(x=V1, y=V5, group=brainPerm, color=brainPerm, fill=brainPerm)) +
  geom_line() +
  facet_wrap(.~V4) +
  theme_bw()

## Do the happy questions here
mod.mixed.hap1 <- mixedmirt(for.irt.mixed[,26:49], covdata = for.irt.mixed[,c("record_id", "ace_score")], model=1, fixed=~ace_score,  random=~1|record_id,itemtype = "2PL")
mod.mixed.hap2 <- mixedmirt(for.irt.mixed[,26:49], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count")], model=1, fixed=~dose_hv_visit_count,  random=~1|record_id,itemtype = "2PL")
mod.mixed.hap3 <- mixedmirt(for.irt.mixed[,26:49], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "ace_score")], model=1, fixed=~dose_hv_visit_count*ace_score, random=~1|record_id, itemtype = "2PL")

## Now try some eeg outcome
mod.mixed.hap4 <- mixedmirt(for.irt.mixed[,26:49], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2happymaxIDfront")], model=1, fixed=~P2happymaxIDfront, random=~1|record_id,itemtype = "2PL")
mod.mixed.hap5 <- mixedmirt(for.irt.mixed[,26:49], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2happymaxIDpost")], model=1, fixed=~P2happymaxIDpost, random=~1|record_id,itemtype = "2PL")
mod.mixed.hap6 <- mixedmirt(for.irt.mixed[,26:49], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2happymaxIDLatfront")], model=1, fixed=~P2happymaxIDLatfront, random=~1|record_id,itemtype = "2PL")
mod.mixed.hap7 <- mixedmirt(for.irt.mixed[,26:49], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2happymaxIDLatpost")], model=1, fixed=~P2happymaxIDLatpost, random=~1|record_id,itemtype = "2PL")
mod.mixed.hap8 <- mixedmirt(for.irt.mixed[,26:49], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2happymaxIDfront", "P2happymaxIDLatfront")], model=1, fixed=~P2happymaxIDfront*P2happymaxIDLatfront, random=list(~1|record_id, ~1|items),itemtype = "2PL")
mod.mixed.hap9 <- mixedmirt(for.irt.mixed[,26:49], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2happymaxIDpost", "P2happymaxIDLatpost")], model=1, fixed=~P2happymaxIDpost*P2happymaxIDLatpost, random=list(~1|record_id, ~1|items),itemtype = "2PL")

## Now onto neutral faces
mod.mixed.neu1 <- mixedmirt(for.irt.mixed[-c(38),50:73], covdata = for.irt.mixed[-c(38),c("record_id", "ace_score")], model=1, fixed=~ace_score,  random=~1|record_id,itemtype = "2PL")
mod.mixed.neu2 <- mixedmirt(for.irt.mixed[-c(38),50:73], covdata = for.irt.mixed[-c(38),c("record_id", "dose_hv_visit_count")], model=1, fixed=~dose_hv_visit_count,  random=~1|record_id,itemtype = "2PL")
mod.mixed.neu3 <- mixedmirt(for.irt.mixed[-c(38),50:73], covdata = for.irt.mixed[-c(38),c("record_id", "dose_hv_visit_count", "ace_score")], model=1, fixed=~dose_hv_visit_count*ace_score, random=~1|record_id, itemtype = "2PL")

## Now try some eeg outcome
mod.mixed.neu4 <- mixedmirt(for.irt.mixed[,50:73], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2neutralmaxIDfront")], model=1, fixed=~P2neutralmaxIDfront, random=~1|record_id,itemtype = "2PL")
mod.mixed.neu5 <- mixedmirt(for.irt.mixed[,50:73], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2neutralmaxIDpost")], model=1, fixed=~P2neutralmaxIDpost, random=~1|record_id,itemtype = "2PL")
mod.mixed.neu6 <- mixedmirt(for.irt.mixed[,50:73], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2neutralmaxIDLatfront")], model=1, fixed=~P2neutralmaxIDLatfront, random=~1|record_id,itemtype = "2PL")
mod.mixed.neu7 <- mixedmirt(for.irt.mixed[,50:73], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2neutralmaxIDLatpost")], model=1, fixed=~P2neutralmaxIDLatpost, random=~1|record_id,itemtype = "2PL")
mod.mixed.neu8 <- mixedmirt(for.irt.mixed[,50:73], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2neutralmaxIDfront", "P2neutralmaxIDLatfront")], model=1, fixed=~P2neutralmaxIDfront*P2neutralmaxIDLatfront, random=~1|record_id,itemtype = "2PL")
mod.mixed.neu9 <- mixedmirt(for.irt.mixed[,50:73], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2neutralmaxIDpost", "P2neutralmaxIDLatpost")], model=1, fixed=~P2neutralmaxIDpost*P2neutralmaxIDLatpost, random=~1|record_id,itemtype = "2PL")

## Now onto unhappy
mod.mixed.unh1 <- mixedmirt(for.irt.mixed[,74:97], covdata = for.irt.mixed[,c("record_id", "ace_score")], model=1, fixed=~ace_score,  random=~1|record_id,itemtype = "2PL")
mod.mixed.unh2 <- mixedmirt(for.irt.mixed[,74:97], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count")], model=1, fixed=~dose_hv_visit_count,  random=~1|record_id,itemtype = "2PL")
mod.mixed.unh3 <- mixedmirt(for.irt.mixed[,74:97], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "ace_score")], model=1, fixed=~dose_hv_visit_count*ace_score, random=~1|record_id, itemtype = "2PL")

## Now try some eeg outcome
mod.mixed.unh4 <- mixedmirt(for.irt.mixed[,74:97], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2unhappymaxIDfront")], model=1, fixed=~P2unhappymaxIDfront, random=~1|record_id,itemtype = "2PL")
mod.mixed.unh5 <- mixedmirt(for.irt.mixed[,74:97], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2unhappymaxIDpost")], model=1, fixed=~P2unhappymaxIDpost, random=~1|record_id,itemtype = "2PL")
mod.mixed.unh6 <- mixedmirt(for.irt.mixed[,74:97], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2unhappymaxIDLatfront")], model=1, fixed=~P2unhappymaxIDLatfront, random=~1|record_id,itemtype = "2PL")
mod.mixed.unh7 <- mixedmirt(for.irt.mixed[,74:97], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2unhappymaxIDLatpost")], model=1, fixed=~P2unhappymaxIDLatpost, random=~1|record_id,itemtype = "2PL")
mod.mixed.unh8 <- mixedmirt(for.irt.mixed[,74:97], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2unhappymaxIDfront", "P2unhappymaxIDLatfront")], model=1, fixed=~P2unhappymaxIDfront*P2unhappymaxIDLatfront, random=~1|record_id,itemtype = "2PL")
mod.mixed.unh9 <- mixedmirt(for.irt.mixed[,74:97], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2unhappymaxIDpost", "P2unhappymaxIDLatpost")], model=1, fixed=~P2unhappymaxIDpost*P2unhappymaxIDLatpost, random=~1|record_id,itemtype = "2PL")

## Now try to collapse high intensity emotions
mod.high.int1 <- mixedmirt(for.irt.mixed[,2:49], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "ace_score")],  model=1, fixed=~ace_score * dose_hv_visit_count,  random=~1|record_id,itemtype = "2PL")

mod.all.em1 <- mixedmirt(for.irt.mixed[,2:97], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "ace_score")],  model=1, fixed=~ace_score * dose_hv_visit_count,  random=~1|record_id,itemtype = "2PL")
mod.all.em2 <- mixedmirt(for.irt.mixed[,2:97], covdata = for.irt.mixed[,c("record_id", "dose_hv_visit_count", "ace_score")],  model=2, fixed=~ace_score * dose_hv_visit_count,  random=list(~1|record_id,~1|items), itemtype = "2PL")

save.image(file="allMixedMods.RData")

## Write the files for the MIMIC models here
## These will be exported to a machine that can run MPlus -- my worktop
cry.dat <- for.irt.mixed[,c(2:25)]
cry.dat <- cbind(cry.dat, for.irt.mixed[,c("record_id", "dose_hv_visit_count", "","P2crymaxCarepost","P2crymaxCareLatpost")])
write.csv(cry.dat, "./forMIMICCry.csv", quote=F, row.names=F)
hap.dat <- for.irt.mixed[,c(26:49)]
hap.dat <- cbind(hap.dat, for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2happymaxCarepost","P2happymaxCareLatpost")])
write.csv(hap.dat, "./forMIMICHap.csv", quote=F, row.names=F)
neu.dat <- for.irt.mixed[,c(50:73)]
neu.dat <- cbind(neu.dat, for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2neutralmaxCarepost","P2neutralmaxCareLatpost")])
write.csv(neu.dat, "./forMIMICNeu.csv", quote=F, row.names=F)
unh.dat <- for.irt.mixed[,c(74:97)]
unh.dat <- cbind(unh.dat, for.irt.mixed[,c("record_id", "dose_hv_visit_count", "P2unhappymaxCarepost","P2unhappymaxCareLatpost")])
write.csv(unh.dat, "./forMIMICUnh.csv", quote=F, row.names=F)

## To do plot out the random effect as a function of brain 
## Two directions: ICC curves;  bar graphs for the 4 groups and plot the difference in item excluding random effect for the item (potentially include facet)
## Amplitude as x axis (high versus low) plot difficulties for 25 percentile, 50th percentile and 75th percentile