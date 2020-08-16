## Test script to organize data for IRT analysis using care task faces
source("~/adroseHelperScripts/R/afgrHelpFunc.R")
install_load("reshape2", "progress")
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
pb <- progress_bar$new(total = length(id.vals))
for(i in id.vals){
  # Isolate the data of interest
  data.to.use <- all.out[which(all.out$participantID==i),]
  ## Now go through an isolate the questions of interest
  # First create the participant specific output
  participant.wide <- NULL
  index <- 1
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
    }
  }
  all.out.wide <- rbind(all.out.wide, participant.wide)
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
  alpha.value <- psych::alpha(data.to.calc)
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
for(i in question.vals){
  data.to.use <- all.out.wide[which(all.out.wide$PictureDisplayed==i),]
  print(dim(data.to.use))
}
