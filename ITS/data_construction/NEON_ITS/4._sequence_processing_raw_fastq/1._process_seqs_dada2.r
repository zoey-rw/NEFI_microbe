#processing new raw fastq files passed along by Lee Stanish October 2018.
#Raw sequences were joined using fastq-join.
#Seqs demultiplexed using split_libraries_fastq.py, minimal filtering as suggested by dada2 tutorial.
#Details in evernote notebook. Unfortunately, processing these required a lot of "by-hand" steps.
#From here we have per-sample fastq files to be processed in dada2.
#clear environment, source paths.
rm(list=ls())
source('paths.r')
library(dada2)
library(Biostrings)
library(ShortRead)
source('NEFI_functions/tic_toc.r')

#Set path to sequences.----
path <- NEON_ITS_fastq.dir
path <- substr(path, 1, nchar(path)-1) #for compatibility with dada2

#set output paths for ESV table and tracking path.----
output.dir <- paste0(path,'/dada2_output/')
cmd <- paste0('mkdir -p ',output.dir)
system(cmd)
#output file path.
output_filepath1 <- paste0(path,'SV_table.rds')
      track.path <- paste0(path,   'track.rds')
output_filepath2 <- NEON_ITS_fastq_SV.table.path


#Grab all file paths.----
seqs <- sort(list.files(path, pattern = '.fastq', full.names = T))
#subset to test
testing = F
if(testing == T){
  seqs <- seqs[1:3]
}

#sample names.
sample.names <- sapply(strsplit(basename(seqs), ".fastq"), `[`, 1)

#Identify ITS primers.----
#We don't need to do this. Primers already removed when sequences delivered.
#Correct primers were ITS1f/ITS2
FWD <- 'CTTGGTCATTTAGAGGAAGTAA' #ITS1F
#FWD <- 'TCCGTAGGTGAACCTGCGG'    #ITS1
#FWD <- 'GCATCGATGAAGAACGCAGC'   #ITS3
REV <- 'GCTGCGTTCTTCATCGATGC'   #ITS2
#REV <- 'TCCTCCGCTTATTGATATGC'   #ITS4
#FWD <- 'CGGCTGCGTTCTTCATCGATGC' #linker-primer sequence

#We are going to check for all possible orientations of the primers.
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

#Remove primers.----
#Not necessary today since they wre already removed.

#Finish filter and trim.----
tic()
cat('Begin quality filtering...\n')
filtered_path <- file.path(path, 'filtered/')
filtered <- file.path(filtered_path, basename(seqs))
out <- filterAndTrim(seqs, filtered, maxN = 0, maxEE = 2,  #make maxEE = c(2,2) if processing multiple files at once (like fwd/rev reads separately).
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = T)
cat('quality filtering complete.')
toc()

#Some reads don't make it through. Update paths and sample.names.
filtered <- list.files(filtered_path)
filtered <- file.path(filtered_path,filtered)
sample.names <- sapply(strsplit(basename(filtered), ".fastq"), `[`, 1)

#Learn the error rates.----
#This is a machine learning algorithm to dial in the error rate model of your reads.
#Colin isn't really sure what an error rate model is. ¯\_(ツ)_/¯
tic() #start timer loop.
cat('Learning error rates...\n')
err <- learnErrors(filtered, multithread = TRUE)
cat('Error models fitted!')
toc() #end timer loop.

#Dereplicate the reads.----
tic()
cat('Dereplicating sequences...\n')
derep <- derepFastq(filtered, verbose = T)
cat('Sequences dereplicated.')
toc()
#add sample names.
names(derep) <- sample.names
      
#Sample Inference. Send it through dada2!-----
tic()
cat('Performing sample inference with the dada2 algorithm...\n')
dada_out <- dada(derep, err = err, multithread = TRUE)
cat('Sample inference complete.')
toc()

#Construct sequence table.----
seqtab <- makeSequenceTable(dada_out)

#Remove Chimeras.----
tic()
cat('Removing chimeric sequences...\n')
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#check proportion of reads that pass chimeras filter.
#lots of ESVs can be chimeras, but they usually represent less than 5-10% of total reads.
chim.retain <- sum(seqtab.nochim / sum(seqtab))
cat('Chimeric sequences removed.',chim.retain*100,'% of sequences retained.')
toc()

#Track reads through pipeline.----
#This is a useful check. If you are losing all your reads at some point this will give you an idea of where.
#this dataframe saves to dada2_output within the sequence folder.
getN <- function(x) sum(getUniques(x))
out_sub <- out[rownames(out) %in% paste0(sample.names,'.fastq'),]
track <- cbind(out_sub, sapply(dada_out, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised","nonchim")
rownames(track) <- sample.names

#Save output.----
saveRDS(seqtab.nochim, output_filepath1)
saveRDS(seqtab.nochim, output_filepath2) #also save to scc_gen
saveRDS(track, track.path)
cat('Analysis complete.')
#end script.
