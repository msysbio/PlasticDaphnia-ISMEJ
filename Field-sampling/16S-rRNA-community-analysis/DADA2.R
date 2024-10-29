"
DADA2 piepline for daphnia field sequences: bacteriolplankton & Daphnia
"

library(dada2)
path <- '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Miseq/Raw_reads/labelled'

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- list.files(path, pattern="_R1.fastq.gz", full.names = TRUE)
fnRs <- list.files(path, pattern="_R2.fastq.gz", full.names = TRUE)

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #extract sample names = forward reads
sample.names2 <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)#extract sample names = reverse reads

#Checking if all names match:
all(sample.names == sample.names2)

#plot quality
plotQualityProfile(fnFs[1:2]) #Forward reads quality
plotQualityProfile(fnRs[1:2]) #Reverse reads quality

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Check if forward and reverse contain the same number of sequences: TAKES LONG TIME
library(ShortRead)
count.reads <- function(fn) {
  length(readFastq(fn))
}
nreadF <- sapply(fnFs, count.reads) 
nreadR <- sapply(fnRs, count.reads) 
table(nreadF == nreadR) #true: 157

#Reads filtering

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200), trimLeft=c(25,25),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)


#Quality after trimming:
plotQualityProfile(filtFs[1:10]) #Forward reads quality
plotQualityProfile(filtRs[1:10]) #Reverse reads quality

#learning the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plot forward
plotErrors(errF, nominalQ=TRUE)
#plot reverse
plotErrors(errR, nominalQ=TRUE)

#Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# MERGING READS
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Removing chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#0.90 remains

#Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


for(i in 1:length(track[,1])){
  for(j in 2:length(track[1,])){
    track[i,j] <- round((track[i,j]/track[i,1])*100, 2)
  }
}

head(track)

#Taxonomy assigment
##Try both databases and see where the annotations are better == fewer unknowns 

write.csv(seqtab.nochim,"/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/Miscellanious/Natasha/RenamedData/tables/seqtab.nochim")

#SILVA
taxa <- assignTaxonomy(seqtab.nochim, "/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Miseq/Databases/Database-updated/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Miseq/Databases/Database-updated/silva_species_assignment_v138.1.fa")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

write.csv(taxa,"/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Miseq/output/dada2/new/taxa")

#Remove mitochondria/chloroplasts and unidentified sequences
#Chloroplast
is.chloroplast <- taxa[,"Order"] %in% "Chloroplast"
seqtable.nochloro <- seqtab.nochim[,!is.chloroplast]
dim(seqtable.nochloro) 
taxa_table.nochloro <- taxa[!is.chloroplast,]
dim(taxa_table.nochloro)

#Mitchondria
is.mitochondria <- taxa_table.nochloro[,"Family"] %in% "Mitochondria"
seqtable.nomito <- seqtable.nochloro[,!is.mitochondria]
taxonomy.nomito <- taxa_table.nochloro[!is.mitochondria,]
dim(seqtable.nomito)

#Unidentified bacteria: Anything with unknown Phylum will be chucked out
is.unknown = taxonomy.nomito[,"Phylum"] %in% NA
seqtable.clean = seqtable.nomito[,!is.unknown]
taxonomy.clean = taxonomy.nomito[!is.unknown,]

write.csv(seqtable.clean, file="/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Miseq/output/dada2/new/otus_clean.csv")
write.csv(taxonomy.clean, file="/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Miseq/output/dada2/new/tax_clean.csv")

##end of pipeline