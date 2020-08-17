## Test script to organize data for IRT analysis using care task faces
source("~/adroseHelperScripts/R/afgrHelpFunc.R")
install_load("reshape2", "progress", "mirt", "psych")
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
  ## Now estimate the cronbach's alpha
  alpha.value <- psych::alpha(data.to.calc, check.keys = T)
  alpha.value <- alpha.value$total[1]
  # Now prepare the output row
  out.row <- cbind(i, alpha.value)
  output.person.con <- rbind(output.person.con, out.row)
}

## Now go through and clauclate the item veraiability
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


## Now run the IRT? and see if it'll work
real.all.wide[real.all.wide=="10"] <- 0
real.all.wide[real.all.wide=="18"] <- 0
real.all.wide[real.all.wide=="12"] <- 1
real.all.wide[real.all.wide=="24"] <- 1
real.all.wide[,2:385] <- apply(real.all.wide[,2:385], c(1,2), function(x) as.numeric(as.character(x)))
for.irt <- apply(real.all.wide[,2:385], c(1,2), function(x) as.numeric(as.character(x)))
mod.1 <- mirt(for.irt, 1, IRTpars=T)
# Now print the coef
out.mod.1 <- coef(mod.1)

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


## Now run IRT
mod.2 <- mirt(for.irt[,-c(1)], 1, IRTpars=T)
