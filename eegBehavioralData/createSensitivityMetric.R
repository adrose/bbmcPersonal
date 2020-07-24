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


## Load library(s)
source("~/adroseHelperScripts/R/afgrHelpFunc.R")
install_load("tidyverse", "reshape2", "gganimate", 'plotly', 'sjPlot')

## Create function(s)

## Load data
in.path <- "/home/arosen/Documents/bbmcPersonal/eegBehavioralData/CARE_mod/"
in.example <- "/home/arosen/Documents/bbmcPersonal/eegBehavioralData/CARE_mod/118-101_OutResponse.csv"
in.data <- read.csv(in.example)

## Convert emotions 


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
out.check$Left <- "Yes"
