########################################
## CAGE data analysis script ###########
#
# Script to generate Transcript and Genewise count tables given CTSS files from CAGE-seq data
# The script is mostly following the same steps as the CAGE pipeline (Rscripts 4 and 5)

#
# Required input:
# - the bam files of interest
# - a GTF file for Gene annotation
# - directory where to perform the analysis
# The genome is hardcoded to HG39 (BSgenome.Hsapiens.UCSC.hg38)
#
# Author: Kevin Menden

#======================#


#==== Hard coded section =======================================#
script_name = "CAGE_gene_counts_from_ctss_v1.0.R"
date = Sys.Date()
current_time = gsub(":", ".", gsub(" ", "_", Sys.time()))
# The genome
bsgenome = "BSgenome.Hsapiens.UCSC.hg38"
#================================================================#

# Load libraries
library(CAGEr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(ChIPseeker)
library(org.Hs.eg.db)
options(stringsAsFactors = F)

setwd("/data/kevin/rimod_cageseq/")
# Parse arguments
ctss_directory = "All_ctss_files"
gtf_file = "gencode.v34.annotation.gtf"
analysis_dir = "/data/kevin/rimod_cageseq/"
tss_range = 3000
annot_level = "gene"
keep_singletons = "TRUE"


# Save parameters in config file
params <- c(current_time, script_name, ctss_directory, gtf_file, analysis_dir, as.character(tss_range), annot_level, keep_singletons)
param.names <- c("Time", "Script_name", "CTSS_directory", "GTF_fiile", "Analysis_dir", "TSS_range", "Annotation_level", "Keep_singletons")
params <- data.frame(param.name = param.names, param = params)
write.table(params, paste("config_file", current_time, sep="_"), quote = F, row.names = F, sep="\t")


### Data preparation
### List all bam files and extract sample labels
paths_to_ctss<- list.files(ctss_directory, full.names=TRUE, pattern = "*.ctss$")
sample_labels <- list.files(ctss_directory, pattern = "*ctss$")
sample_labels <- gsub(".ctss", "", sample_labels)

#==== Data loading and analysis ==================#

# Creat TxDB genome
gencode_TxDb <- makeTxDbFromGFF(gtf_file, format = "gtf", dataSource = "GENCODE",  organism = "Homo sapiens")

annoPeakList <- list()

### Aggregrate CTSS files
# parse it
# generate peak bed file
# call the peaks
# annotate the peaks

# make cage-seq object
myCAGEset <- new("CAGEset",
                 genomeName = bsgenome,
                 inputFiles = paths_to_ctss,
                 inputFilesType = "ctss",
                 sampleLabels = sample_labels)
# parsing
getCTSS(myCAGEset, removeFirstG = TRUE, correctSystematicG = TRUE)

# setting ctss object
ctss <- CTSStagCount(myCAGEset)

# normalization
normalizeTagCount(myCAGEset, method="none")

# clustering
print("Clustering...")
clusterCTSS(object = myCAGEset, threshold = 10, thresholdIsTpm = FALSE,
            nrPassThreshold = 2, method = "distclu", maxDist =20,
            removeSingletons = TRUE)


print("aggreagate cluters")
aggregateTagClusters(myCAGEset, tpmThreshold = 0, qLow = NULL, qUp = NULL, maxDist = 100)

countTable <- myCAGEset@consensusClustersTpmMatrix

for_deseq_consensus_cluster <- consensusClusters(myCAGEset)

for_deseq_consensus_cluster$id <- paste(for_deseq_consensus_cluster$chr,
                                        for_deseq_consensus_cluster$start, for_deseq_consensus_cluster$end,
                                        for_deseq_consensus_cluster$strand, sep="_")

rownames(countTable) <- for_deseq_consensus_cluster$id
write.table(countTable, "RiMod_CAGEseq_cluster_count_table_all.txt", quote=F, sep="\t")

#=============================================#


##############
# Annotation as genes
##############

# Bring in bed-file format
countBed <- myCAGEset@consensusClustersTpmMatrix
consClus <- consensusClusters(myCAGEset)
chr <- consClus$chr
start <- consClus$start
end <- consClus$end
df <- data.frame(chr, start, end)
countBed <- cbind(df, countBed)

peak_file <- "cage_consClusters_peak_file.bed"
write.table(countBed, peak_file, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
# Generate peak object
print("Peak calling ...")
peaks <- readPeakFile(peak_file)
# Transform to GenomicRange object
peaks_gr <- GRanges(peaks, strand = Rle(strand(consClus$strand)))

# define promoters around 3K;
promoter_3k <- getPromoters(TxDb=gencode_TxDb, upstream=tss_range, downstream=tss_range, by="gene")
length(promoter_3k)

# Get a tag-matrix around that promoter region;
peaks_gr_tagmat_3k <- getTagMatrix(peaks_gr, windows=promoter_3k)

# Plots Tags with heatmap
print("Plotting tag matrix ...")
pdf("tag_heatmap_peaks_3k.pdf")
tagHeatmap(peaks_gr_tagmat_3k,  xlim=c(-tss_range, tss_range), xlab="Genomic range", ylab="Transcripts", color="red")
dev.off()
# Plot average profile heatmaps
pdf("tag_avgprf_heatmap_peaks_3k.pdf")
plotAvgProf(peaks_gr_tagmat_3k,  xlim=c(-tss_range, tss_range), xlab="Genomic range", ylab="Read Count Frequency")
dev.off()


## Peak annotation

#=============== Genewise =======================#
# Gene Annotation within the range of +/- 3kbases;
print("Peak annotation ...")
peaks_gr_3k <- annotatePeak(peaks_gr, tssRegion=c(-tss_range, tss_range),
                            TxDb=gencode_TxDb,
                            level = annot_level,
                            assignGenomicAnnotation = TRUE,
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                            annoDb="org.Hs.eg.db",
                            addFlankGeneInfo = FALSE, flankDistance = 5000, verbose = TRUE)


### Plot the results of the annotation
pdf("annotation_gene_piechart.pdf")
plotAnnoPie(peaks_gr_3k)
dev.off()
pdf("annotation_gene_barplot.pdf")
plotAnnoBar(peaks_gr_3k)
dev.off()
pdf("annotation_gene_upsetplot_.pdf")
upsetplot(peaks_gr_3k)
dev.off()
pdf("annotation_gene_vennpie.pdf")
vennpie(peaks_gr_3k)
dev.off()
pdf("annotation_gene_disttoTSS.pdf")
plotDistToTSS(peaks_gr_3k, title="Distribution of CAGE clusters relative to TSS")
dev.off()


# Get annotates dataframe
annoDF <- as.data.frame(peaks_gr_3k)
# get correct colnames back
colnames(annoDF)[6:253] <- colnames(countTable)


# test if order is correct
count_test <- as.numeric(countTable[1,])
anno_test <- as.numeric(annoDF[1,])[6:253]
all(count_test == anno_test)
# [1] TRUE
# --> order is still correct and we can happily save and have our data

write.table(annoDF, "RiMod_genewise_annotated_CAGE_clusters.txt", sep="\t", quote=F)








