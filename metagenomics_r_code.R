setwd("/cloud/project/fastq_files")

#### Importing files from SRA Explorer ####
urls <- c('ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/042/SRR22993542/SRR22993542_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/042/SRR22993542/SRR22993542_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/032/SRR22993532/SRR22993532_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/032/SRR22993532/SRR22993532_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/043/SRR22993543/SRR22993543_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/043/SRR22993543/SRR22993543_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/041/SRR22993541/SRR22993541_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/041/SRR22993541/SRR22993541_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/044/SRR22993544/SRR22993544_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/044/SRR22993544/SRR22993544_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/033/SRR22993533/SRR22993533_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/033/SRR22993533/SRR22993533_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/047/SRR22993547/SRR22993547_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/047/SRR22993547/SRR22993547_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/046/SRR22993546/SRR22993546_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/046/SRR22993546/SRR22993546_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/031/SRR22993531/SRR22993531_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/031/SRR22993531/SRR22993531_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/045/SRR22993545/SRR22993545_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/045/SRR22993545/SRR22993545_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/049/SRR22993549/SRR22993549_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/049/SRR22993549/SRR22993549_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/030/SRR22993530/SRR22993530_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/030/SRR22993530/SRR22993530_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/051/SRR22993551/SRR22993551_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/051/SRR22993551/SRR22993551_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/050/SRR22993550/SRR22993550_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/050/SRR22993550/SRR22993550_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/029/SRR22993529/SRR22993529_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/029/SRR22993529/SRR22993529_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/028/SRR22993528/SRR22993528_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/028/SRR22993528/SRR22993528_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/040/SRR22993540/SRR22993540_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/040/SRR22993540/SRR22993540_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/048/SRR22993548/SRR22993548_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/048/SRR22993548/SRR22993548_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/039/SRR22993539/SRR22993539_1.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/039/SRR22993539/SRR22993539_2.fastq.gz',
          'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/038/SRR22993538/SRR22993538_1.fastq.gz',
          "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/038/SRR22993538/SRR22993538_2.fastq.gz",
          "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/037/SRR22993537/SRR22993537_1.fastq.gz",
          "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/037/SRR22993537/SRR22993537_2.fastq.gz",
          "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/036/SRR22993536/SRR22993536_1.fastq.gz",
          "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/036/SRR22993536/SRR22993536_2.fastq.gz",
          "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/035/SRR22993535/SRR22993535_1.fastq.gz",
          "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/035/SRR22993535/SRR22993535_2.fastq.gz",
          "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/034/SRR22993534/SRR22993534_1.fastq.gz",
          "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/034/SRR22993534/SRR22993534_2.fastq.gz")
for (url in urls) {
  command <- paste("wget", url, "&")
  system(command)
}

#### Unzipping Files ####
file_list <- list.files(pattern = "\\.gz$")
for (file in file_list) {
  command <- paste("gzip -d", file)
  system(command)
}

#### importing trimmed seq files ####

setwd("/cloud/project/trimmed_seqs")
#uploaded zip file of all fastq files

#### Getting Ready ####

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
library(dada2); packageVersion("dada2")

path <- setwd("/cloud/project/trimmed_seqs") 
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)


#### Inspect read quality profiles ####

plotQualityProfile(fnRs[1:2])
plotQualityProfile(fnFs[1:2])

#### Filter and Trim ####

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                           maxN=0, maxEE=c(3,5), truncQ=2,rm.phix=TRUE,
                           compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

#getting an error that it can't find 34_1.fastq file so trying to figure out why
readLines("trimmed_SRR22993534_1.fastq")
readLines("trimmed_SRR22993534_2.fastq")
#FIXED for some reason fastq file for 34_1 did not make it to my folder so when i trimmed everything in cutadapt, it produced no output, thus empty files

#error bc vectors are not same length in out <- filterAndTrim command
length(fnFs)
length(filtFs)
length(fnRs)
length(filtRs)
length(unique(sample.names))
#FIXED had to change multithread to false

#### Learn Error Rates ####

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)

#### Sample Inference ####

dadaFs <- dada(filtFs, err=errF, multithread = TRUE)
dadaRs <- dada(filtRs, err=errR, multithread = TRUE)
dadaFs[[1]]

#### Merge Paired Reads ####

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, maxMismatch = 1, verbose=TRUE) #allowing for 1 mismatch to ensure adequate merging
head(mergers[[1]])

#### Construct sequence table ####

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#### Remove Chimeras ####

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) 
#chimeras account for ~10% of merged seq reads when considering abundance of variants

#### Track reads through the pipeline ####

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")
rownames(track) <- sample.names
head(track)
track

#### Assign Taxonomy ####

#restart R first
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("DECIPHER", force = TRUE)
library(DECIPHER)

dna <- DNAStringSet(getSequences(seqtab.nochim))

system("wget http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData")

load("SILVA_SSU_r138_2019.RData") 
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

#### Handoff to phyloseq ####

#installing necessary packages
BiocManager::install("phyloseq")
BiocManager::install("Biostrings")
install.packages("ggplot2")
install.packages("readxl")

#loading libraries
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(readxl)
library(viridisLite)
library(RColorBrewer)
library(vegan)
theme_set(theme_bw())

meta_data_xl <- read_excel("TakinMetadata.xlsx")
head(meta_data_xl)

#converting to data frame
meta_data <- as.data.frame(meta_data_xl)
rownames(meta_data) <- sample_names(otu_table(seqtab.nochim, taxa_are_rows = FALSE))

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               sample_data(meta_data),
               tax_table(taxid))
ps <- prune_samples(sample_names(ps) != "Mock", ps)

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#### plots ####

#abundance
plot_richness(ps, x="Season", measures = c("Shannon", "Simpson"), color = "TakinID")
plot_richness(ps, x="Season", measures = c("Observed", "Chao1"), color = "TakinID")

#alpha diversity
alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Observed", "Chao1"))
alpha_div

#ordinate

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Season", title="Bray NMDS")

# sorting taxa into top 20
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
tax_table(ps.top20)
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

#this is to sort it by Takin_# rather than sample
get_numeric_part <- function(x) as.numeric(gsub("Takin_", "", x))
sorted_levels <- unique(sample_data(ps.top20)$TakinID[order(get_numeric_part(sample_data(ps.top20)$TakinID))])
sample_data(ps.top20)$TakinID <- factor(sample_data(ps.top20)$TakinID, levels = sorted_levels)

plot_bar(ps.top20, x="TakinID", fill="phylum") + facet_wrap(~Season, scales="free_x")

# Create distance matrix
distance_matrix <- dist(otu_table(ps))
# Perform NMDS ordination
nmds <- ordinate(ps, method = "NMDS", distance = distance_matrix)
plot_ordination(ps, nmds, color = "Season", type = "samples", title = "NMDS Ordination of Microbial Community Samples by Season")

#### Get the data out of phyloseq! ####
  
install.packages("remotes")
remotes::install_github("vmikk/metagMisc")

# Export the taxonomy and abundance data
library(metagMisc); packageVersion("metagMisc")
export_ps <- phyloseq_to_df(ps)
setwd("/cloud/project/ASV_CVS")
write.table(export_ps, file = "16SASVData.csv", col.names=TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

#export metadata
export_MD <- as.matrix(meta_data)
write.table(export_MD, file = "16S_Metadata.csv", col.names=TRUE, row.names = TRUE, sep = "\t", quote = FALSE)


