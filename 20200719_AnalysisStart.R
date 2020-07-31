library(purrr)
library(tidyverse)
library(edgeR)
f_files<- list.files("neo-cgsi/macs2/S4", pattern = "nar", full.names = T) #List of files to concatenate
read_in_feature_counts<- function(file){
  cnt<- read_tsv(file, col_names = F, comment = "#")
  colnames(cnt) <- c("chr", "startCoord", "endCoord", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
  cnt<- cnt %>% dplyr::select(1, 2, 3)
  return(cnt) #Function to concat the files, but only adds the Chr, Start, End, Str, Length from the first
}
TM3Seq_counts<- map(f_files, read_in_feature_counts) #Does the function above for the list given
TM3Seq_counts_df<- purrr::reduce(TM3Seq_counts, inner_join) #Makes into a data frame 
BiocManager::install("org.Hs.eg.db")
nearestTSS("chr5", 385724, species="Hs")
TM3Seq_counts_df$nearTSS <- nearestTSS(TM3Seq_counts_df$chr, 0.5*(TM3Seq_counts_df$startCoord + TM3Seq_counts_df$endCoord), species="Hs")
TM3Seq_counts_df <- TM3Seq_counts_df[TM3Seq_counts_df$nearTSS$distance > -10000 & TM3Seq_counts_df$nearTSS$distance < 10000,]
p <- ggplot(TM3Seq_counts_df, aes(x=nearTSS$distance)) + geom_histogram(bins=1000)
