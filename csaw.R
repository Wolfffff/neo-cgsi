# Ref:
# https://bioconductor.org/packages/3.11/workflows/vignettes/chipseqDB/inst/doc/h3k9ac.html


BiocManager::install('org.Hs.eg.db')
library(csaw)
library(Rsamtools)
library(edgeR)
library(ChIPpeakAnno)
library(statmod)
library(tidyverse)
library(Gviz)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
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
export(con="brant_results.bed", object=simplified, format="csv")




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



prom <- suppressWarnings(promoters(txdb,
                                   upstream=3000, downstream=1000, columns=c("tx_name", "gene_id")))
entrez.ids <- sapply(prom$gene_id, FUN=function(x) x[1]) # Using the first Entrez ID.
gene.name <- select(org.Hs.eg.db, keys=entrez.ids, keytype="ENTREZID", column="SYMBOL")
prom$gene_name <- gene.name$SYMBOL[match(entrez.ids, gene.name$ENTREZID)]
head(prom)

olap.out <- overlapResults(filtered.data, regions=prom, res$table)
olap.out

simple <- DataFrame(ID=prom$tx_name, Gene=prom$gene_name, olap.out$combined)
simple[!is.na(simple$PValue),]




gax <- GenomeAxisTrack(col="black", fontsize=15, size=2)
greg <- GeneRegionTrack(txdb, showId=TRUE,
                        geneSymbol=TRUE, name="", background.title="transparent")
symbols <- unlist(mapIds(org.Hs.eg.db, gene(greg), "SYMBOL",
                         "ENTREZID", multiVals = "first"))
symbol(greg) <- symbols[gene(greg)]