### This script will be used to plot the DIF effects when the interaction of AMP & LAT are significant
### This script will not run the IRT models; but will take the output from MPlus and plot the ICC in ggplot2
### The models will be writeen on the PC and transfered over to the linux workstation



# load-library--------------------------------------------------------
library("tidyverse")
library("GGally")
library("reshape2")

# ----- load-data--------------------------------------------------------
## Now load the latest two stage data from Emma
amp.vals <- read.csv("./data/eegData/ampVals.csv")
lat.vals <- read.csv("./data/eegData/latVals.csv")
## Now read in the coefficients
load("./data/outParamsNU.RData")
out.list.nu <- out.list 
load("./data/outParams.RData")

# ---- mod-data--------------------------------------------------------
## NOw add these back into the eeg Pred values
amp.vals <- tidyr::pivot_wider(amp.vals, id_cols="ERPset", names_from="binlabel", values_from="value") %>% 
  mutate(record_id = as.character(ERPset)) %>% 
  mutate(record_id = if_else(record_id %in% "117-115", "118-64", record_id)) %>%
  melt(.) %>% 
  mutate(cycle = substring(sapply(lapply(strsplit(as.character(variable), NULL), rev), paste, collapse=""), 1, 1)) %>%
  mutate(variable = gsub('.{1}$', '', variable)) %>% 
  select(-ERPset) %>% 
  tidyr::pivot_wider(., id_cols=c(record_id, cycle), names_from=c("variable"), values_from=c(value))
colnames(amp.vals) <- gsub(colnames(amp.vals), pattern = ' ', replacement = "")

lat.vals <- tidyr::pivot_wider(lat.vals, id_cols="ERPset", names_from="binlabel", values_from="value") %>% 
  mutate(record_id = as.character(ERPset)) %>% 
  mutate(record_id = if_else(record_id %in% "117-115", "118-64", record_id)) %>%
  melt(.) %>% 
  mutate(cycle = substring(sapply(lapply(strsplit(as.character(variable), NULL), rev), paste, collapse=""), 1, 1)) %>%
  mutate(variable = gsub('.{1}$', '', variable)) %>% 
  select(-ERPset) %>% 
  tidyr::pivot_wider(., id_cols=c(record_id, cycle), names_from=c("variable"), values_from=c(value))
colnames(lat.vals) <- gsub(colnames(amp.vals), pattern = ' ', replacement = "")

## Create the z score values for all of the imaging data
all.eeg <- merge(amp.vals, lat.vals, by=c("record_id", "cycle"), suffixes = c("_Amp", "_Lat"))
all.eeg[,3:10] <- scale(all.eeg[,3:10])[,]

## Now prepare the interaction terms
emo.vals <- c("Crying", "Unhappy", "Neutral", "Happy")
for(e in emo.vals){
  ## First declare the column names
  amp.col <- paste(e, "_Amp", sep='')
  lat.col <- paste(e, "_Lat", sep='')
  int.col <- paste(e, "_Int", sep='')
  # Now estimate the cross prod
  all.eeg[,int.col] <- all.eeg[,amp.col] * all.eeg[,lat.col]
}
# prep-functions--------------------------------------------------------
## Prepare functions
## Create a function which will operate on unifrom dif items
returnUnifromICC <- function(baseDIF=NULL, baseDIS=NULL, p2AmpMod=NULL, p2LatMod=NULL, p2IntMod=NULL, p2AmpSubj=NULL, p2LatSubj=NULL, p2IntSubj=NULL){
  ## First create a range of values for which the ICC will be calulcated over
  theta.values <- seq(-4, 4, by=.11)
  ### Now create the formula which will be used to return proabability of endorsment across the entire range of theta values
  all.vals <- -(baseDIF + p2LatMod * p2LatSubj + p2AmpMod * p2AmpMod + p2IntMod * p2IntSubj) + baseDIS * theta.values / 1 + -(baseDIF + p2LatMod * p2LatSubj + p2AmpMod * p2AmpMod + p2IntMod * p2IntSubj) + baseDIS * theta.values
  all.vals <- exp(all.vals) / (1 + exp(all.vals))
  all.vals <-  cbind(theta.values, all.vals)
  all.vals <- data.frame(all.vals)
  colnames(all.vals) <- c("theta", "probEndorse")
  ### Now return these values
  return(all.vals)
}

## Now create f unction which will perform the same steps but for nonunifrom DIF
returnNunifromICC <- function(baseDIF=NULL, baseDIS=NULL, p2AmpMod=NULL, p2LatMod=NULL, p2IntMod=NULL, p2AmpSubj=NULL, p2LatSubj=NULL, p2IntSubj=NULL, p2Ampn=NULL, p2Latn=NULL, p2Intn=NULL){
  ## First create a range of values for which the ICC will be calulcated over
  theta.values <- seq(-6, 6, by=.11)
  ### Now create the formula which will be used to return proabability of endorsment across the entire range of theta values
  all.vals <- -(baseDIF + p2LatMod * p2LatSubj + p2AmpMod * p2AmpMod + p2IntMod * p2IntSubj) + (baseDIS + p2Ampn*p2AmpSubj + p2Intn*p2Intn + p2Intn*p2IntSubj ) * theta.values / 1 + -(baseDIF + p2LatMod * p2LatSubj + p2AmpMod * p2AmpMod + p2IntMod * p2IntSubj) + (baseDIS + p2Ampn*p2AmpSubj + p2Intn*p2Intn + p2Intn*p2IntSubj ) * theta.values
  all.vals <- exp(all.vals) / (1 + exp(all.vals))
  all.vals <-  cbind(theta.values, all.vals)
  all.vals <- data.frame(all.vals)
  colnames(all.vals) <- c("theta", "probEndorse")
  ### Now return these values
  return(all.vals)
}

## Now in a loop prepare all of the coefficient values
emotion.value <- NULL
base.dif.vals <- NULL
base.dis.vals <- NULL
p2AmpModvals <- NULL
p2LatModvals <- NULL
p2IntModvals <- NULL
### Now work with the mplus output
for(i in 1:4){
  tmp.data <- out.list[[i]]
  ## Now identify significant interaction DIF
  uni.dif.t <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(paramHeader %in% paste0("V",1:24,".ON")) %>%
    dplyr::filter(param=="P2INT") %>% 
    dplyr::filter(pval<0.05) %>% 
    dplyr::select(item)
  ## Now isolate the params for this item
  base.dif.vals.t <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>%
    dplyr::filter(paramHeader=="Thresholds") %>%
    dplyr::filter(param %in% paste("V", uni.dif.t$item, "$1", sep='')) %>%
    ## Now isolate the item the model was trained on
    dplyr::mutate(itemTrain=paste("V", item, "$1", sep='')) %>% 
    dplyr::filter(param==itemTrain) %>% 
    dplyr::select(est)
  base.dis.vals.t <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>%
    dplyr::filter(paramHeader=="F.BY") %>%
    dplyr::filter(param %in% paste("V", uni.dif.t$item, sep='')) %>%
    ## Now isolate the item the model was trained on
    dplyr::mutate(itemTrain=paste("V", item, sep='')) %>% 
    dplyr::filter(param==itemTrain) %>% 
    dplyr::select(est)
  p2AmpModvals.t <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>%
    dplyr::filter(paramHeader %in% paste("V", uni.dif.t$item,".ON", sep='')) %>%
    dplyr::filter(param=="P2AMP") %>% 
    dplyr::select(est)  
  p2LatModvals.t <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>%
    dplyr::filter(paramHeader %in% paste("V", uni.dif.t$item,".ON", sep='')) %>%
    dplyr::filter(param=="P2LAT") %>% 
    dplyr::select(est)
  p2IntModvals.t <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>%
    dplyr::filter(paramHeader %in% paste("V", uni.dif.t$item,".ON", sep='')) %>%
    dplyr::filter(param=="P2INT") %>% 
    dplyr::select(est)
  emotion.sub <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>% 
    dplyr::select(emotion)
  base.dif.vals <- append(base.dif.vals, base.dif.vals.t$est)
  base.dis.vals <- append(base.dis.vals, base.dis.vals.t$est)
  p2AmpModvals <- append(p2AmpModvals, p2AmpModvals.t$est)
  p2LatModvals <- append(p2LatModvals, p2LatModvals.t$est)
  p2IntModvals <- append(p2IntModvals, p2IntModvals.t$est)
  emotion.value <- append(emotion.value, rep(unique(emotion.sub$emotion), length(uni.dif.t$item)))  
}
all.plot.vals <- NULL
colnames(all.eeg) <- tolower(colnames(all.eeg))
## Now loop through all of the base values
for(i in 1:length(base.dif.vals)){
  base.dif <- base.dif.vals[i]
  base.dis <- base.dis.vals[i]
  p2.amp <- p2AmpModvals[i]
  p2.lat <- p2LatModvals[i]
  p2.int <- p2IntModvals[i]
  e <- emotion.value[i]
  ## Now declare the column names so we can isolate the subject values
  amp.col <- paste(e, "_amp", sep='')
  lat.col <- paste(e, "_lat", sep='')
  int.col <- paste(e, "_int", sep='')
  all.subj.vals <- unique(all.eeg[,c(amp.col, lat.col, int.col)])
  ## Now loop through all of the subject values and estimate the ICC for each eeg pattern
  for(s in 1:dim(all.subj.vals)[1]){
    tmp.vals <- returnUnifromICC(baseDIF = base.dif, baseDIS = base.dis, p2AmpMod = p2.amp, p2LatMod = p2.lat, p2IntMod = p2.int, p2AmpSubj = all.subj.vals[s,1], p2LatSubj = all.subj.vals[s,2], p2IntSubj = all.subj.vals[s,3])
    # Now prepare these for the output dataframe
    subj.amp <- all.subj.vals[s,1]
    subj.lat <- all.subj.vals[s,2]
    subj.int <- all.subj.vals[s,3]
    out.vals <- cbind(i, base.dif, base.dis, p2.amp, p2.lat, p2.int, e, s, subj.amp, subj.lat, subj.int, tmp.vals)
    all.plot.vals <- rbind(all.plot.vals, out.vals)
  }
}

## Now clean these data
all.plot.vals$subjMagnitude <- abs(all.plot.vals$subj.int)

## Now plot these
all.plot.vals %>% ggplot(., aes(x=theta, y=probEndorse, color=subj.int, group=s)) +
  geom_line() +
  theme_bw() +
  facet_wrap(i~ e ) +
  ggtitle("uDIF ICCs") +
  ylab("Probability of Correct Response") +
  xlab("Theta")
## Now loop through each plot and return an individual png for each file so I can prep it for some vector art
for(q in unique(all.plot.vals$i)){
  tmp.plot <- all.plot.vals[which(all.plot.vals$i==q),] %>% ggplot(., aes(x=theta, y=probEndorse, color=subj.int, group=s)) +
    geom_line(size=2) +
    theme_bw() +
    ylab("") +
    xlab("") +
    theme(legend.position = "none", text = element_text(size=32))
  file.name <- paste("./paperImages/", unique(all.plot.vals[which(all.plot.vals$i==q),"e"]), "_", q, ".png", sep='')
  png(filename = file.name)
  print(tmp.plot)
  dev.off()
}

## Now make a density plot of the interaction term
x <- all.plot.vals$subj.int[which(all.plot.vals$i==q)]
y <- density(x, n = 2^12, kernel = "triangular", adjust = 8)
png("./paperImages/interactionHistogram.png")
ggplot(data.frame(x = y$x, y = y$y), aes(x, y)) + geom_line() + 
  geom_segment(aes(xend = x, yend = 0, colour = x)) +
  xlab("Interaction Magnitude") +
  ylab("Density") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size=32))
dev.off()

##############################################################
#### Work with nonunifrom DIF here
##############################################################
## Now in a loop prepare all of the coefficient values
emotion.value <- c("Crying", "Crying", "Crying", "Neutral", "Neutral", "Neutral", "Happy")
base.dif.vals <- c(-2.025, -1.529, -2.258, -0.905, -1.009, -0.911, -2.329)
base.dis.vals <- c(1.318, 1.476, 1.107, 1.569, 1.534, 0.768, 1.393)
p2AmpModvals <-  c(0.225, 0.424, 0.060, 0.193, -0.372, -0.136, 0.014)
p2LatModvals <-  c(0.613, -0.540, 0.355, -0.138, 0.386, 0.008, -0.064)
p2IntModvals <-  c(-0.526, -0.841, 0.087, -0.132, -0.758, -0.043, -0.263)
p2AmpModDvals <- c(-0.059, 0.508, -0.338, -0.043, -0.741, -0.580, -0.480)
p2LatModDvals <- c(0.034, -0.693, 0.233, -0.129, -0.153, 0.007, 0.389)
p2IntModDvals <- c(-1.239, -1.127, -0.873, -0.801, -0.809, -0.603, 0.622)
all.plot.vals <- NULL

## Now in a loop prepare all of the coefficient values
emotion.value <- NULL
base.dif.vals <- NULL
base.dis.vals <- NULL
p2AmpModvals <- NULL
p2LatModvals <- NULL
p2IntModvals <- NULL
p2AmpModDvals <- NULL
p2LatModDvals <- NULL
p2IntModDvals <- NULL
### Now work with the mplus output
for(i in 1:4){
  tmp.data <- out.list.nu[[i]]
  ## Now identify significant interaction DIF
  uni.dif.t <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(paramHeader %in% paste0("V",1:24,".ON")) %>%
    dplyr::filter(param=="INTERACT4") %>% 
    dplyr::filter(pval<0.05) %>% 
    dplyr::select(item)
  ## Now isolate the params for this item
  base.dif.vals.t <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>%
    dplyr::filter(paramHeader=="Thresholds") %>%
    dplyr::filter(param %in% paste("V", uni.dif.t$item, "$1", sep='')) %>%
    ## Now isolate the item the model was trained on
    dplyr::mutate(itemTrain=paste("V", item, "$1", sep='')) %>% 
    dplyr::filter(param==itemTrain) %>% 
    dplyr::select(est)
  base.dis.vals.t <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>%
    dplyr::filter(paramHeader=="F.BY") %>%
    dplyr::filter(param %in% paste("V", uni.dif.t$item, sep='')) %>%
    ## Now isolate the item the model was trained on
    dplyr::mutate(itemTrain=paste("V", item, sep='')) %>% 
    dplyr::filter(param==itemTrain) %>% 
    dplyr::select(est)
  p2AmpModvals.t <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>%
    dplyr::filter(paramHeader %in% paste("V", uni.dif.t$item,".ON", sep='')) %>%
    dplyr::filter(param=="P2AMP") %>% 
    dplyr::select(est)  
  p2LatModvals.t <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>%
    dplyr::filter(paramHeader %in% paste("V", uni.dif.t$item,".ON", sep='')) %>%
    dplyr::filter(param=="P2LAT") %>% 
    dplyr::select(est)
  p2IntModvals.t <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>%
    dplyr::filter(paramHeader %in% paste("V", uni.dif.t$item,".ON", sep='')) %>%
    dplyr::filter(param=="P2INT") %>% 
    dplyr::select(est)
  p2AmpModDvals.t <-tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>%
    dplyr::filter(paramHeader %in% paste("V", uni.dif.t$item,".ON", sep='')) %>%
    dplyr::filter(param=="INTERACT2") %>% 
    dplyr::select(est)  
  p2LatModDvals.t <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>%
    dplyr::filter(paramHeader %in% paste("V", uni.dif.t$item,".ON", sep='')) %>%
    dplyr::filter(param=="INTERACT3") %>% 
    dplyr::select(est)
  p2IntModDvals.t <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>%
    dplyr::filter(paramHeader %in% paste("V", uni.dif.t$item,".ON", sep='')) %>%
    dplyr::filter(param=="INTERACT4") %>% 
    dplyr::select(est)
  
  emotion.sub <- tmp.data %>% bind_rows() %>% 
    dplyr::filter(item %in% c(uni.dif.t$item)) %>% 
    dplyr::select(emotion)
  base.dif.vals <- append(base.dif.vals, base.dif.vals.t$est)
  base.dis.vals <- append(base.dis.vals, base.dis.vals.t$est)
  p2AmpModvals <- append(p2AmpModvals, p2AmpModvals.t$est)
  p2LatModvals <- append(p2LatModvals, p2LatModvals.t$est)
  p2IntModvals <- append(p2IntModvals, p2IntModvals.t$est)
  p2AmpModDvals <- append(p2AmpModDvals, p2AmpModDvals.t$est)
  p2LatModDvals <- append(p2LatModDvals, p2LatModDvals.t$est)
  p2IntModDvals <- append(p2IntModDvals, p2IntModDvals.t$est)
  emotion.value <- append(emotion.value, rep(unique(emotion.sub$emotion), length(uni.dif.t$item)))  
}

## Now loop through all of the base values
for(i in 1:length(base.dif.vals)){
  base.dif <- base.dif.vals[i]
  base.dis <- base.dis.vals[i]
  p2.amp <- p2AmpModvals[i]
  p2.lat <- p2LatModvals[i]
  p2.int <- p2IntModvals[i]
  p2.ampn <- p2AmpModDvals[i]
  p2.latn <- p2LatModDvals[i]
  p2.intn <- p2IntModDvals[i]
  e <- emotion.value[i]
  ## Now declare the column names so we can isolate the subject values
  amp.col <- paste(e, "_amp", sep='')
  lat.col <- paste(e, "_lat", sep='')
  int.col <- paste(e, "_int", sep='')
  all.subj.vals <- unique(all.eeg[,c(amp.col, lat.col, int.col)])
  ## Now loop through all of the subject values and estimate the ICC for each eeg pattern
  for(s in 1:dim(all.subj.vals)[1]){
    tmp.vals <- returnNunifromICC(baseDIF = base.dif, baseDIS = base.dis, p2AmpMod = p2.amp, p2LatMod = p2.lat, 
                                  p2IntMod = p2.int, p2AmpSubj = all.subj.vals[s,1], p2LatSubj = all.subj.vals[s,2], p2IntSubj = all.subj.vals[s,3],
                                  p2Ampn = p2.ampn, p2Latn = p2.latn, p2Intn = p2.intn)
    # Now prepare these for the output dataframe
    subj.amp <- all.subj.vals[s,1]
    subj.lat <- all.subj.vals[s,2]
    subj.int <- all.subj.vals[s,3]
    out.vals <- cbind(i, base.dif, base.dis, p2.amp, p2.lat, p2.int, e, s, subj.amp, subj.lat, subj.int, tmp.vals)
    all.plot.vals <- rbind(all.plot.vals, out.vals)
  }
}
## Now clean these data
all.plot.vals$subjMagnitude <- abs(all.plot.vals$subj.int)
## Remove subj 50 from q 8
all.plot.vals <- all.plot.vals[-which(all.plot.vals$s==50 & all.plot.vals$i==8),]
## Now plot these
tmp <- all.plot.vals %>% ggplot(., aes(x=theta, y=probEndorse, color=subj.int, group=s)) +
  geom_line() +
  theme_bw() +
  facet_wrap(i ~ e ) +
  ggtitle("nuDIF ICCs") +
  ylab("Probability of Correct Response") +
  xlab("Theta")
ggplotly(tmp)
## Now loop through each plot and return an individual png for each file so I can prep it for some vector art
for(q in unique(all.plot.vals$i)){
  tmp.plot <- all.plot.vals[which(all.plot.vals$i==q),] %>% ggplot(., aes(x=theta, y=probEndorse, color=subj.int, group=s)) +
    geom_line(size=2) +
    theme_bw() +
    ylab("") +
    xlab("") +
    theme(legend.position = "none", text = element_text(size=32))
  file.name <- paste("./paperImages/nonUni/", unique(all.plot.vals[which(all.plot.vals$i==q),"e"]), "_", q, ".png", sep='')
  png(filename = file.name)
  print(tmp.plot)
  dev.off()
}
