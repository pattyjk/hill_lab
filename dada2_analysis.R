#from: https://benjjneb.github.io/dada2/tutorial.html

#install if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.22")

#load DADA2
library(dada2)

#set the path where your raw fastq is. This should be changed specifically for where your data is in your computer
path<-"~/Documents/GitHub/hill_lab/demux_fastq/"

#make sure the files are there, would be embarssing if they weren't and you ran all this code...
list.files(path)

#get a list of forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

#plot quality scores to figure out where to chop
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#pull out sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#make a lsit of file names for writing the files to
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#quality filter and trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 13,
                     truncLen=c(100,100),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

#look at the output and see what we're losing with QC
head(out)

#learn error rate for forward and then rverse reads
errF <- learnErrors(filtFs, multithread=TRUE)
#100424448 total bases in 1154304 reads from 39 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
#100424448 total bases in 1154304 reads from 39 samples will be used for learning the error rates.

#plot error rates
plotErrors(errF, nominalQ=TRUE)

#run dada
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#merge paired ends reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap = 5)

#make ASV table
seqtab_merge <- makeSequenceTable(mergers)
seqtab <- makeSequenceTable(dadaFs)
seqtab_rev <- makeSequenceTable(dadaRs)

#look at depth
rowSums(seqtab_merge)
rowSums(seqtab)
rowSums(seqtab_rev)

#the reverse reads didn't do a great job at being merged, so let's jsut focus on the formward reads goign forward, no big deal

#remove chimeras by denovo clustering
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#get dimensions of ASV table
dim(seqtab.nochim)
#[1]   48 1114

#write table
write.table(seqtab.nochim, '~/Documents/GitHub/hill_lab/asv_table.txt', row.names=T, quote=F, sep='\t')

#if you want to read in the ASV table again
library(dada2)
seqtab.nochim<-read.delim('~/Documents/GitHub/hill_lab/pat_analysis/asv_table.txt', header=T, row.names=1)

#if want to use it in DADA2 need as a matrix, so
seqtab.nochim<-as.matrix(seqtab.nochim)

#assign taxonomy with SILVA 138 (newest as of 12/29/25) and make into a data frame
taxa<-assignTaxonomy(seqtab.nochim, "~/Downloads/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)
taxa2<-as.data.frame(taxa)
taxa2$asv<-row.names(taxa2)

#save taxonomy so don't have to rerun
write.table(taxa2, '~/Documents/GitHub/hill_lab/pat_analysis/taxonomy.txt', sep='\t', quote=F, row.names=F)

#read in taxonomy if necessary
taxa2<-read.delim('~/Documents/GitHub/hill_lab/pat_analysis/taxonomy.txt', header=T)

#prune eukaryotes from the dataset
which(taxa2$Kingdom =='Eukaryota')
#[1]  124  191  315  331  340  360  550  577  689  732  756  760  762  798  873  916  934  985 1025 1068 1076 1102
#22 ASVs are eukaryotes

euk_asvs <- rownames(taxa)[taxa[, "Kingdom"] == "Eukaryota"]
seqtab.noeuk <- seqtab.nochim[, !colnames(seqtab.nochim) %in% euk_asvs]

#double check make sure nothing weird happened
dim(seqtab.nochim)
#[1]   48 1114

dim(seqtab.noeuk)
#[1]   48 1092

#double check math
1114-1092
#[1] 22
#math checks out!

#write asv table to file for safe keeping
write.table(seqtab.noeuk, '~/Documents/GitHub/hill_lab/pat_analysis/asv_table_no_euks.txt', quote=F, sep='\t', row.names = T)

##Use the negative controls to 'decontaminate' the real samples
#install new packages if needed, https://joey711.github.io/phyloseq/ & https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
library(phyloseq) 
library(decontam)

#read in metadata
meta<-read.delim('~/Documents/GitHub/hill_lab/hill_lab_map', header=T)

#if you need to re-read in ASV table
seqtab.noeuk<-read.delim('~/Documents/GitHub/hill_lab/pat_analysis/asv_table_no_euks.txt', header=T, row.names=1)

#transpose the dataframe to make things easier to work with
seqtab.noeuk<-as.data.frame(t(seqtab.noeuk))

# Keep only samples that exist in both files
##not necessary here but in case you get a dataset that needs that (e.g. a sample didn't sequence)
shared.samples <- intersect(colnames(seqtab.noeuk), meta$SampleID)
seqtab.noeuk <- seqtab.noeuk[, shared.samples]
meta <- meta[match(shared.samples, meta$SampleID), ]

# Convert to phyloseq objects, makes things easier for decontam package
OTU <- otu_table(as.matrix(seqtab.noeuk), taxa_are_rows = TRUE)
SAM <- sample_data(meta)
sample_names(SAM) <- meta$SampleID  # ensure they match
sample_names(OTU) <- colnames(seqtab.noeuk)
ps <- phyloseq(OTU, SAM)

# Identify negative controls
# Adjust this line if your controls are named differently, you can change 'NegativeControl' or add more thigns with a comma
sample_data(ps)$is.neg <- sample_data(ps)$Type %in% c("NegativeControl")

# Run decontam based on prevalance method
#check the other methods in the help file, but the prevalance method is most commonly used
contamdf.prev <- isContaminant(ps, method = "prevalence", neg = "is.neg")

# Output contaminants
contaminants <- as.data.frame(rownames(contamdf.prev[contamdf.prev$contaminant, ]))
names(contaminants)<-'OTUS'
write.table(contaminants, "~/Documents/GitHub/hill_lab/pat_analysis/contaminants.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#how many? (n=14)
dim(contaminants)
#[1] 14  1

#remove contaminants
dim(seqtab.noeuk)
#[1] 1094   48

otu_no_contam<-seqtab.noeuk[-which(row.names(seqtab.noeuk) %in% contaminants$OTUS), ]
dim(otu_no_contam)
#[1] 1080   48

#the dim isn't necessary, I'm jsut paranoid about not removing what I say I want removed

#write to file
write.table(otu_no_contam, '~/Documents/GitHub/hill_lab/pat_analysis/no_contam_asv_table.txt', row.names = T, sep='\t', quote=F)

##Data analysis time, after all that processing stuffs
#load the needed libraries
library(vegan)
library(ggplot2)

#re-read in ASV table if needed (aka you start here instead of rerunning everything)
otu_no_contam<-read.delim('~/Documents/GitHub/hill_lab/pat_analysis/no_contam_asv_table.txt', header=T, row.names = 1)

#look at sequencing depth
colSums(otu_no_contam)
min(colSums(otu_no_contam))
#[1] 6018

#plot rarefaction curve
rarecurve(seqtab.nochim, step=75)

#rarefy table
asv_rare<-rrarefy(seqtab.nochim, sample = 6010)

#calculate alpha diversity
div_asv<-as.data.frame(specnumber(asv_rare))
asv_shan<-as.data.frame(diversity(asv_rare, index = 'shannon'))

#calculate beta diversity with bray-curtis similarity
pcoa1<-capscale(asv_rare~1, distance = 'bray')

#extract coordinates for plotting
ko.scores<-scores(pcoa1)
str(ko.scores)
ko.scores2<-as.data.frame(ko.scores$sites)
ko.scores2$SampleID<-row.names(ko.scores2)

#plot it
ggplot(ko.scores2, aes(MDS1, MDS2, label=SampleID))+
  geom_text()+
  theme_bw()
