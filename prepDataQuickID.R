# --- prep-environ ---------------------
## Clean everything
rm(list=ls())

## Test script to organize data for IRT analysis using ID task faces
source("~/adroseHelperScripts/R/afgrHelpFunc.R")
install_load("reshape2", "progress", "mirt", "psych", "ggplot2", "polycor", "corrplot", "ggalluvial", "GGally", "tidyverse")

## Declare any functions
binary.flip <- function (x) {
  x * -1 + 1
}

# ---- read-data -------------------------
# Load behavioral data
load("./eegBehavioralData/IRTIDData.RData")

## Now load the participants demographics
load("./fname.gz")
ds_eegPred <- out.data[[9]]

## Now load the latest two stage data from Emma
amp.vals <- read.csv("./data/eegData/ampVals.csv")
lat.vals <- read.csv("./data/eegData/latVals.csv")

## Now load the old two stage data from Emma
amp.vals.old <- read.csv("./data/eegData/febDump/ampVals.csv")
lat.vals.old <- read.csv("./data/eegData/febDump/latVals.csv")

## Now load the latest 4:1 data from emma
four.to.one <- read.csv("./data/eegData/fourToOne.csv")

# ---- prep-IRT-data --------------------------
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
for.irt3$record_id <- for.irt2$V1
for.irt3 <- for.irt3 %>% mutate(., record_id = as.character(record_id)) %>% 
  mutate(., record_id = if_else(record_id %in% "117-115", "118-64", record_id))

# 118-103 118-62 118-81

# ---- prep-EEG-data ---------------------------
amp.vals$ERPset <- gsub(strSplitMatrixReturn(as.character(amp.vals$ERPset), "_")[,1], pattern = "#", replacement = '')
lat.vals$ERPset <- gsub(strSplitMatrixReturn(as.character(lat.vals$ERPset), "_")[,1], pattern = "#", replacement = '')


amp.vals <- tidyr::pivot_wider(amp.vals, id_cols="ERPset", names_from="binlabel", values_from="value") %>% 
  mutate(record_id = as.character(ERPset)) %>% 
  mutate(record_id = if_else(record_id %in% "117-115", "118-64", record_id)) %>%
  melt(.) %>% 
  mutate(cycle = substring(sapply(lapply(strsplit(as.character(variable), NULL), rev), paste, collapse=""), 1, 1)) %>%
  mutate(variable = gsub('.{1}$', '', variable)) %>% 
  select(-ERPset) %>% 
  tidyr::pivot_wider(., id_cols=c(record_id, cycle), names_from=c("variable"), values_from=c(value))

amp.vals.old <- tidyr::pivot_wider(amp.vals.old, id_cols="ERPset", names_from="binlabel", values_from="value") %>% 
  mutate(record_id = as.character(ERPset)) %>% 
  mutate(record_id = if_else(record_id %in% "117-115", "118-64", record_id)) %>%
  melt(.) %>% 
  mutate(cycle = substring(sapply(lapply(strsplit(as.character(variable), NULL), rev), paste, collapse=""), 1, 1)) %>%
  mutate(variable = gsub('.{1}$', '', variable)) %>% 
  select(-ERPset) %>% 
  tidyr::pivot_wider(., id_cols=c(record_id, cycle), names_from=c("variable"), values_from=c(value))
colnames(amp.vals.old)[3:6] <- paste(gsub("\\s","",names(amp.vals.old)[3:6]), "_Amp_Old", sep='')

lat.vals <- tidyr::pivot_wider(lat.vals, id_cols="ERPset", names_from="binlabel", values_from="value") %>% 
  mutate(record_id = as.character(ERPset)) %>% 
  mutate(record_id = if_else(record_id %in% "117-115", "118-64", record_id)) %>%
  melt(.) %>% 
  mutate(cycle = substring(sapply(lapply(strsplit(as.character(variable), NULL), rev), paste, collapse=""), 1, 1)) %>%
  mutate(variable = gsub('.{1}$', '', variable)) %>% 
  select(-ERPset) %>% 
  tidyr::pivot_wider(., id_cols=c(record_id, cycle), names_from=c("variable"), values_from=c(value))

lat.vals.old <- tidyr::pivot_wider(lat.vals.old, id_cols="ERPset", names_from="binlabel", values_from="value") %>% 
  mutate(record_id = as.character(ERPset)) %>% 
  mutate(record_id = if_else(record_id %in% "117-115", "118-64", record_id)) %>%
  melt(.) %>% 
  mutate(cycle = substring(sapply(lapply(strsplit(as.character(variable), NULL), rev), paste, collapse=""), 1, 1)) %>%
  mutate(variable = gsub('.{1}$', '', variable)) %>% 
  select(-ERPset) %>% 
  tidyr::pivot_wider(., id_cols=c(record_id, cycle), names_from=c("variable"), values_from=c(value))
colnames(lat.vals.old)[3:6] <- paste(gsub("\\s","",names(lat.vals.old)[3:6]), "_Lat_Old", sep='')

multi.vals <- full_join(amp.vals, lat.vals, by=c("record_id", "cycle"), suffix=c("_Amp", "_Lat")) %>% 
  full_join(., amp.vals.old, by=c("record_id", "cycle")) %>% 
  full_join(., lat.vals.old, by=c("record_id", "cycle"))

names(multi.vals) <- make.names(names(multi.vals))
colnames(multi.vals) <- gsub(names(multi.vals), pattern = "X", replacement = "")
colnames(multi.vals) <- gsub(colnames(multi.vals), pattern = "\\.", replacement = "", perl=F)

## Now widen the multi.vals
wide.multi <- pivot_wider(multi.vals, id_cols = c("record_id"), values_from = c("Crying_Amp", "Unhappy_Amp","Neutral_Amp","Happy_Amp","Crying_Lat","Unhappy_Lat","Neutral_Lat","Happy_Lat","Crying_Amp_Old","Unhappy_Amp_Old","Neutral_Amp_Old","Happy_Amp_Old","Crying_Lat_Old","Unhappy_Lat_Old","Neutral_Lat_Old","Happy_Lat_Old"), names_from = c("cycle"))

## Now add the new four to one ratio data
# First fix the record  id column
four.to.one$record_id <- gsub(strSplitMatrixReturn(as.character(four.to.one$record_id), "_")[,1], pattern = "#", replacement = '')
## Now spread the data
four.to.one <- four.to.one %>% pivot_wider(., id_cols=c("record_id"), values_from=c("amplitude", "latency"), names_from=c("emotion")) %>% 
  mutate(record_id = if_else(record_id %in% "117-115", "118-64", record_id))

## Now combine all of these values
all.eeg.wide <- four.to.one %>% full_join(., wide.multi, by="record_id") 
all.eeg.long <- four.to.one %>% full_join(., multi.vals, by="record_id") %>% 
  mutate(cycle=as.numeric(cycle))


# --- combine-values ----------------------------------
# The following section will combine the following data formats:
# for.irt3 ==== the binary correct versus incorrect response data from the ID task
# all.eeg.wide === eeg values with cycle data as the repeated value
# ds_eegPred === original data w/ the ace & dose values
original.target <- for.irt3

## Add a cycle to the behavioral data 
original.target <- original.target %>% group_by(record_id) %>% 
  mutate(round=1:n()) %>% 
  mutate(cycle=1) %>% 
  mutate(cycle=(if_else(round>2, 2, 1)))

## Now add the eeg data w/ cycle
original.target <- original.target %>% full_join(., all.eeg.long, by=c("record_id", "cycle"))

## Now add the additional demographic.vars 
final.out <- original.target %>% full_join(., ds_eegPred, by=c("record_id"))

# --- write-data ----------------------------------
write.csv(final.out, "./data/forMPlus.csv", quote=F, row.names=F)
