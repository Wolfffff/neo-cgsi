library(purrr)
library(tidyverse)
library(edgeR)
library(tools)

curr_time <- function(){
  return(format(Sys.time(), "%d%m%Y_%H%M"))
}

plots_dir <- "results/plots/"
narrowPeak_cols <- c("chr", "startCoord", "endCoord", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")

files <- list.files(pattern = "nar", full.names = T,recursive=T) #List of files to concatenate

samples = list()

for (file in files) {
  temp_read <- read_tsv(file,col_names = F)
  colnames(temp_read) <- narrowPeak_cols
  temp_read$nearTSS <- nearestTSS(temp_read$chr, 0.5*(temp_read$startCoord + temp_read$endCoord), species="Hs")
  samples[[file_path_sans_ext(basename(file))]]  <- temp_read
}

for (name in names(samples)){
  sample = samples[[name]]
  print(ggplot(sample,aes(x=nearTSS$distance)) +
    geom_histogram(aes(y=..density..),binwidth=1000,col="darkblue") +
    xlim(-1e5,1e5) +
    xlab("Distance") +
    ylab("Count") + 
    ggtitle(paste0("Nearest TSS Distribution - ",name," - ",curr_time())) +
    theme_minimal())
  ggsave(file=paste0(plots_dir,"NearestTSS_",name,"_",curr_time(),".png"))
}
ggplot(samples[[1]],aes(x=nearTSS$distance)) +
          geom_histogram(aes(y=..density..,fill="red"),binwidth=2000,alpha=.4) +
          geom_histogram(data=samples[[2]],aes(nearTSS$distance,y=..density..,,fill="green"),binwidth=2000,alpha=.4) +
          geom_histogram(data=samples[[3]],aes(x=nearTSS$distance,y=..density..,fill="blue"),binwidth=2000,alpha=.4) +
          xlim(-1e5,1e5) +
          xlab("Distance") +
          ylab("Density") +
          ggtitle("Nearest TSS Distribution") +
          scale_fill_discrete(name = "Sample",labels=c(names(samples))) +
  theme_minimal()
ggsave(file=paste0(plots_dir,"NearestTSS_GroupedPlot_",curr_time(),".png"))

ggplot(samples[[1]],aes(x=nearTSS$distance)) +
  geom_density(aes(y=..density..,col="red")) +
  geom_density(data=samples[[2]],aes(nearTSS$distance,y=..density..,,col="green")) +
  geom_density(data=samples[[3]],aes(x=nearTSS$distance,y=..density..,col="blue")) +
  xlim(-1e5,1e5) +
  xlab("Distance") +
  ylab("Density") +
  ggtitle("Nearest TSS Distribution") +
  scale_color_discrete(name = "Sample",labels=c(names(samples))) +
  scale_fill_manual(values = c("red","green","blue"), name = "text",labels=c(names(samples))) +
  theme_minimal() 
ggsave(file=paste0(plots_dir,"NearestTSS_GroupedPlot_KDE_",curr_time(),".png"))
