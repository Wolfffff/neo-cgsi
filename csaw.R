# Ref:
# https://bioconductor.org/packages/3.11/workflows/vignettes/chipseqDB/inst/doc/h3k9ac.html

library(csaw)
library(Rsamtools)
library(edgeR)
library(ChIPpeakAnno)
library(statmod)
library(tidyverse)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(hash)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


files <- list.files(pattern = ".bam$", full.names = T,recursive=T) #List of files to concatenate
acdata <- data.frame(Path=files,Name=basename(files),Description = c("S1","S4","S7"),Lipid=c("No","No","Yes"))

diagnostics <- list()
for (bam in files) {
  total <- countBam(bam)$records
  mapped <- countBam(bam, param=ScanBamParam(
    flag=scanBamFlag(isUnmapped=FALSE)))$records
  marked <- countBam(bam, param=ScanBamParam(
    flag=scanBamFlag(isUnmapped=FALSE, isDuplicate=TRUE)))$records
  diagnostics[[basename(bam)]] <- c(Total=total, Mapped=mapped, Marked=marked)
}
diag.stats <- data.frame(do.call(rbind, diagnostics))
diag.stats$Prop.mapped <- diag.stats$Mapped/diag.stats$Total*100
diag.stats$Prop.marked <- diag.stats$Marked/diag.stats$Mapped*100
diag.stats

standard.chr <- paste0("chr", c(1:23, "X", "Y"))
param <- readParam(minq=20, restrict=standard.chr)

x <- correlateReads(files, param=reform(param, dedup=TRUE))
frag.len <- maximizeCcf(x)
frag.len

plot(1:length(x)-1, x, xlab="Delay (bp)", ylab="CCF", type="l")
abline(v=frag.len, col="red")
text(x=frag.len, y=min(x), paste(frag.len, "bp"), pos=4, col="red")


win.data <- windowCounts(files, param=param, width=105, ext=frag.len)
win.data

bins <- windowCounts(files, bin=TRUE, width=2000, param=param)
filter.stat <- filterWindowsGlobal(win.data, bins)
min.fc <- 3
keep <- filter.stat$filter > log2(min.fc)
summary(keep)

hist(filter.stat$back.abundances, main="", breaks=50,
     xlab="Background abundance (log2-CPM)")
threshold <- filter.stat$abundances[1] - filter.stat$filter[1] + log2(min.fc)
abline(v=threshold, col="red")

filtered.data <- win.data[keep,]


win.ab <- scaledAverage(filtered.data)
adjc <- calculateCPM(filtered.data, use.offsets=FALSE)
logfc <- adjc[,3] - adjc[,2]
smoothScatter(win.ab, logfc, ylim=c(-6, 6), xlim=c(0, 5),
              xlab="Average abundance", ylab="Log-fold change")

lfit <- smooth.spline(logfc~win.ab, df=5)
o <- order(win.ab)
lines(win.ab[o], fitted(lfit)[o], col="red", lty=2)

filtered.data <- normOffsets(filtered.data)
offsets <- assay(filtered.data, "offset")
head(offsets)

norm.adjc <- calculateCPM(filtered.data, use.offsets=TRUE)
norm.fc <- norm.adjc[,3]-norm.adjc[,2]
smoothScatter(win.ab, norm.fc, ylim=c(-6, 6), xlim=c(0, 5),
              xlab="Average abundance", ylab="Log-fold change")

lfit <- smooth.spline(norm.fc~win.ab, df=5)
lines(win.ab[o], fitted(lfit)[o], col="red", lty=2)

celltype <- acdata$Lipid
celltype <- factor(celltype)
design <- model.matrix(~0+celltype)
colnames(design) <- levels(celltype)
design


y <- asDGEList(filtered.data)
str(y)

y <- estimateDisp(y, design)
summary(y$trended.dispersion)

plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

summary(fit$df.prior)

plotMDS(norm.adjc, labels=celltype,
        col=c("red", "blue")[as.integer(celltype)])

contrast <- makeContrasts(No-Yes, levels=design)
res <- glmQLFTest(fit, contrast=contrast)
head(res$table)

merged <- mergeResults(filtered.data, res$table, tol=100, 
                       merge.args=list(max.width=5000))
merged$regions

tabcom <- merged$combined
tabcom

is.sig <- tabcom$FDR <= 0.05
summary(is.sig)

table(tabcom$direction[is.sig])

tabbest <- merged$best
tabbest

is.sig.pos <- (tabbest$rep.logFC > 0)[is.sig]
summary(is.sig.pos)

out.ranges <- merged$regions
mcols(out.ranges) <- DataFrame(tabcom,
                               best.pos=mid(ranges(rowRanges(filtered.data[tabbest$rep.test]))),
                               best.logFC=tabbest$rep.logFC)
saveRDS(file="brant_results.rds", out.ranges)

simplified <- out.ranges[is.sig]
simplified$score <- -10*log10(simplified$FDR)
export(con="brant_results.bed", object=simplified)



out.ranges = FDR_filtered

readRDS(file="brant_results.rds", out.ranges)



# Interpretation

anno <- detailRanges(out.ranges, orgdb=org.Hs.eg.db, txdb=txdb)

head(anno$overlap)
head(anno$left)
head(anno$right)

meta <- mcols(out.ranges)
mcols(out.ranges) <- data.frame(meta, anno)


data(TSS.human.GRCh38)
minimal <- out.ranges
elementMetadata(minimal) <- NULL
anno.regions <- annotatePeakInBatch(minimal, AnnotationData=TSS.human.GRCh38)
colnames(elementMetadata(anno.regions))



prom <- suppressWarnings(promoters(txdb,upstream=3000, downstream=1000, columns=c("tx_name", "gene_id")))
entrez.ids <- sapply(prom$gene_id, FUN=function(x) x[1]) # Using the first Entrez ID.
gene.name <- select(org.Hs.eg.db, keys=entrez.ids, keytype="ENTREZID", column="SYMBOL")
prom$gene_name <- gene.name$SYMBOL[match(entrez.ids, gene.name$ENTREZID)]


olap.out <- overlapResults(filtered.data, regions=prom, res$table)
olap.out

simple <- DataFrame(ID=prom$tx_name, Gene=prom$gene_name, olap.out$combined)
simple[!is.na(simple$PValue),]

PRGs <- read_csv("PRGs-2.csv")

PRG_mapping = hash()
for (gene in PRGs$Gene) {
  PRG_mapping[[gene]] =simple[which(simple$Gene == toupper(gene)),]
}

PRG_mapping = simple[0,]
for (gene in PRGs$Gene) {
  PRG_mapping = rbind(PRG_mapping,simple[which(simple$Gene == toupper(gene)),])
}
PRG_mapping_FDR10 <- PRG_mapping[which(PRG_mapping$FDR < 0.10),]
PRG_mapping_FDR10 <- PRG_mapping_FDR10[order(PRG_mapping_FDR10$FDR),]

PRG_mapping_FDR01 <- PRG_mapping[which(PRG_mapping$FDR < 0.01),]
PRG_mapping_FDR01 <- PRG_mapping_FDR01[order(PRG_mapping_FDR01$FDR),]




sorted = simple[order(simple$FDR),]
FDR_filtered_10 = sorted[which(sorted$FDR < 0.10),]
FDR_filtered_01 = sorted[which(sorted$FDR < 0.01),]


nrow(PRG_mapping_FDR10)/length(unique(PRG_mapping$Gene))
nrow(FDR_filtered_10)/length(unique(sorted$Gene))



gax <- GenomeAxisTrack(col="black", fontsize=15, size=2)
greg <- GeneRegionTrack(txdb, showId=TRUE,
                        geneSymbol=TRUE, name="", background.title="transparent")
symbols <- unlist(mapIds(org.Hs.eg.db, gene(greg), "SYMBOL",
                         "ENTREZID", multiVals = "first"))
symbol(greg) <- symbols[gene(greg)]

o <- order(out.ranges$PValue)
sorted.ranges <- out.ranges[o]
sorted.ranges

cur.region <- sorted.ranges[1]
cur.region


collected <- list()
lib.sizes <- filtered.data$totals/1e6
for (i in seq_along(acdata$Path)) {
  reads <- extractReads(bam.file=acdata$Path[i], cur.region, param=param)
  cov <- as(coverage(reads)/lib.sizes[i], "GRanges")
  collected[[i]] <- DataTrack(cov, type="histogram", lwd=0, ylim=c(0,10),
                              name=acdata$Description[i], col.axis="black", col.title="black",
                              fill="darkgray", col.histogram=NA)
}
plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
           from=start(cur.region), to=end(cur.region))



# No complex differential binding

complex <- sorted.ranges$num.up.logFC > 0 & sorted.ranges$num.down.logFC > 0
cur.region <- sorted.ranges[complex][2]
cur.region

collected <- list()
for (i in seq_along(acdata$Path)) {
  reads <- extractReads(bam.file=acdata$Path[i], cur.region, param=param)
  cov <- as(coverage(reads)/lib.sizes[i], "GRanges")
  collected[[i]] <- DataTrack(cov, type="histogram", lwd=0, ylim=c(0,3),
                              name=acdata$Description[i], col.axis="black", col.title="black",
                              fill="darkgray", col.histogram=NA)
}
plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
           from=start(cur.region), to=end(cur.region))


# Simple DB across small region

sharp <- sorted.ranges$num.tests < 20
cur.region <- sorted.ranges[sharp][1]
cur.region

collected <- list()
for (i in seq_along(acdata$Path)) {
  reads <- extractReads(bam.file=acdata$Path[i], cur.region, param=param)
  cov <- as(coverage(reads)/lib.sizes[i], "GRanges")
  collected[[i]] <- DataTrack(cov, type="histogram", lwd=0, ylim=c(0,3),
                              name=acdata$Description[i], col.axis="black", col.title="black",
                              fill="darkgray", col.histogram=NA)
}
plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
           from=start(cur.region), to=end(cur.region))


sessionInfo()
