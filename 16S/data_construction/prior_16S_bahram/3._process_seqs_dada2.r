#Processing Bahram sequence data.
#dada2 tutorial here: https://benjjneb.github.io/dada2/tutorial.html
#NOTE: this is definitely a job for the cluster. Get a node with as much RAM and as many cores as you can.
#pair ends, quality filter, generate ESV tables using dada2. You will need the cluster for (1) memory and (2) parallel computing.
#remove junk seqs using the positive and negative controls.
#assign taxonomy to ESVs using RDP via dada2.
#save ESV (otu) and taxonomy tables. Make sure you can link these to mapping file of environmental data.
#you may need to trim primers off fastq files first. This seems useful: http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/Biostrings/html/trimLRPatterns.html

#load dada2 (installed via the bioconductor package. Instruction in link above.)
rm(list=ls())
library(dada2)
source('paths.r')
source('NEFI_functions/tic_toc.r')


#start with a test directory that includes only two samples, forward and reverse reads.
#specify forward and reverse read motifs. 
#in an ideal world, this is all you would need to input.
#Should have this output all figures and summary of which reads passed as separate files to check.
path <- "/projectnb/talbot-lab-data/NEFI_data/big_data/bahram_test"
#path <- bahram.seq.dir
forward.read.motif <- '_1.fastq'
reverse.read.motif <- '_2.fastq'

#set output for ESV table and tracking path.
output.dir <- paste0(path,'/dada2_output/')
cmd <- paste0('mkdir -p ',output.dir)
system(cmd)
esv.table.path <- paste0(output.dir,'esv_table.rds')
    track.path <- paste0(output.dir,    'track.rds')

#### read in the names of the fastq files. ####
fnFs <- sort(list.files(path, pattern = forward.read.motif, full.names = T))
fnRs <- sort(list.files(path, pattern = reverse.read.motif, full.names = T))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#### quality filtering and trucantion. ####
#We're going to perform some quality filtering and truncation to clip off the parts of the reads where the uality scores get gnarly.
#Colin wishes there was an algorithm to decide this, rather than just eyeballing it.
#If Zoey is reading this in the future, she should create this algorithm and call it as a function!
#Just find where median line crosses some threshold on average or minimum for a bunch of samples, then use that.

#setup filtered files in a filtered sub directory.
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#Inspect quality score patterns.
#You can get buy on doing this with like two samples. Its just a visual check.
#quality scores will drop off at the end of the read. This will happen sooner for reverse reads.
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])

#choose truncation length.
#Colin is calling this whenever the green line drops below 30, according to his eyeball.
truncation.length.forward <- 240
truncation.length.reverse <- 150
#filter
tic()
cat('Begin quality filtering...\n')
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncation.length.forward,truncation.length.reverse),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
cat('quality filtering complete.')
toc()

#check whether you lost all your reads or not.
#head(out)


#### "Learn" the error rates. ####
#this takes a while, even with not a lot of reads.
#this works with a subset of the data. So it takes about the same time with big or small data sets.
#default number of reads is 1,000,000. You could speed up by dropping the nreads parameter.
#This is a machine learning algorithm to dial in the error rate model of your reads.
#Colin isn't really sure what an error rate model is. ¯\_(ツ)_/¯
tic() #start timer loop.
cat('Learning error rates...\n')
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
cat('Error models fitted!')
toc() #end timer loop.

#sanity check plot. Errors should drop as quality scores increase.
#plotErrors(errF, nominalQ = T)


####Dereplicate the sequences. ####
#This combines identical reads in unique sequences.
#dada2 does a little more than this: it keeps quality information for downstream inference.
#there is a mod to do this in parallel described here: https://benjjneb.github.io/dada2/bigdata.html
tic()
cat('Dereplicating sequences...\n')
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
cat('Sequences dereplicated.')
toc()
# Name the derep objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


#### Sample Inference ####
#RAM requirements go way down if you mod this as described in the bigdata tutorial for dada2.
tic()
cat('Performing sample inference with the dada2 algorithm...\n')
dadaFs <- dada(derepFs, err=errF, multithread = T)
dadaRs <- dada(derepRs, err=errR, multithread = T)
cat('Sample inference complete.')
toc()


#### Merge paired reads. ####
tic()
cat('Merging paired ends...\n')
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
cat('Paired ends merged.')
toc()

#### Construct sequence table. ####
#rownames are the samples. columnnames are sequences of the ESVs.
seqtab <- makeSequenceTable(mergers)

#### Remove chimeras. ####
tic()
cat('Removing chimeric sequences...\n')
seqtab.nochim <- removeBimeraDenovo(seqtab, method = 'consensus', multithread = T)
#check proportion of reads that pass chimeras filter.
#lots of ESVs can be chimeras, but they usually represent less than 5-10% of total reads.
chim.retain <- sum(seqtab.nochim / sum(seqtab))
cat('Chimeric sequences removed.',chim.retain*100,'% of sequences retained.')
toc()


#### Track reads through pipeline. ####
#This is a useful check. If you are losing all your reads at some point this will give you an idea of where.
#this dataframe saves to dada2_output within the sequence folder.
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

#### save output path, tracking file and some plots. ####
saveRDS        (track,     track.path)
saveRDS(seqtab.nochim, esv.table.path)

#quality score plot.
#plotQualityProfile(fnFs[1:2])

#error rate plot.
#plotErrors(errF, nominalQ = T)

cat('Analysis complete.')
#end script.
