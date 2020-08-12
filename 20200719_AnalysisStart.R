library(purrr)
library(tidyverse)
library(edgeR)
library(tools)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

library(scales)

mylog10_trans <- function (base = 10) 
{
  trans <- function(x) log(x + 1, base)
  inv <- function(x) base^x
  trans_new(paste0("log-", format(base)), trans, inv, log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

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
  temp_TSS <- nearestTSS(temp_read$chr, 0.5*(temp_read$startCoord + temp_read$endCoord), species="Hs")
  colnames(temp_TSS) = paste0("TSS_Nearest.",colnames(temp_TSS))
  temp_read <- cbind(temp_read,temp_TSS)
  samples[[file_path_sans_ext(basename(file))]]  <- temp_read
}

for (name in names(samples)){
  sample = samples[[name]]
  print(ggplot(sample,aes(x=TSS_Nearest.distance)) +
    geom_histogram(,binwidth=1000,col="darkblue") +
    xlim(-1e5,1e5) +
    xlab("Distance") +
    ylab("Count") + 
    ggtitle(paste0("Nearest TSS Distribution - ",name," - ",curr_time())) +
    theme_minimal()) +
    scale_y_log10()
   ggsave(file=paste0(plots_dir,"NearestTSS_LogCount_",name,"_",curr_time(),".png"))
}
ggplot(samples[[1]],aes(x=TSS_Nearest.distance)) +
          geom_histogram(aes(y=..density..,fill="red"),binwidth=2000,alpha=.4) +
          geom_histogram(data=samples[[2]],aes(TSS_Nearest.distance,y=..density..,,fill="green"),binwidth=2000,alpha=.4) +
          geom_histogram(data=samples[[3]],aes(TSS_Nearest.distance,y=..density..,fill="blue"),binwidth=2000,alpha=.4) +
          xlim(-1e5,1e5) +
          xlab("Distance") +
          ylab("Density") +
          ggtitle("Nearest TSS Distribution") +
          scale_fill_discrete(name = "Sample",labels=c(names(samples))) +
          theme_minimal()
# ggsave(file=paste0(plots_dir,"NearestTSS_GroupedPlot_",curr_time(),".png"))

ggplot(samples[[1]],aes(x=TSS_Nearest.distance)) +
  geom_density(aes(y=..density..,col="red")) +
  geom_density(data=samples[[2]],aes(TSS_Nearest.distance,y=..density..,,col="green")) +
  geom_density(data=samples[[3]],aes(x=TSS_Nearest.distance,y=..density..,col="blue")) +
  xlim(-1e5,1e5) +
  xlab("Distance") +
  ylab("Density") +
  ggtitle("Nearest TSS Distribution") +
  scale_color_discrete(name = "Sample",labels=c(names(samples))) +
  scale_fill_manual(values = c("red","green","blue"), name = "text",labels=c(names(samples))) +
  theme_minimal() 
# ggsave(file=paste0(plots_dir,"NearestTSS_GroupedPlot_KDE_",curr_time(),".png"))





# covplot(peaks, weightCol="V5")

peaks = readPeakFile(files[[1]])
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peaks, windows=promoter)
sample_title = "S1"
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red",title=sample_title,xlab="Distance from TSS",ylab="Peak #")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency",title=sample_title)
peakAnno <- annotatePeak(files[[1]], tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)

peaks = readPeakFile(files[[2]])
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peaks, windows=promoter)
sample_title = "S4"
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red",title=sample_title,xlab="Distance from TSS",ylab="Peak #")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency",title=sample_title)
peakAnno <- annotatePeak(files[[2]], tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)


peaks = readPeakFile(files[[3]])
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peaks, windows=promoter)
sample_title = "S7"
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red",title=sample_title,xlab="Distance from TSS",ylab="Peak #")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency",title=sample_title)
peakAnno <- annotatePeak(files[[3]], tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)


et <- exactTest(y)
topTags(peaks)

