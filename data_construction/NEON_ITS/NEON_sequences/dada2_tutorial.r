#testing dada2
#based on tutorial found here: https://benjjneb.github.io/dada2/tutorial.html
cur.dir <- '/fs/data3/caverill/NEON_ASVs/ITS/test/MiSeq_SOP'

#get file paths and name
F.seqs <- sort(list.files(cur.dir, pattern = '_R1_001.fastq', full.names = T))
R.seqs <- sort(list.files(cur.dir, pattern = '_R2_001.fastq', full.names = T))
sample.names <- sapply(strsplit(basename(F.seqs), "_"), `[`, 1)

filt.F.seqs <- file.path(cur.dir, "filtered", paste0(sample.names,"_F_filt.fastq.gz"))
filt.R.seqs <- file.path(cur.dir, "filtered", paste0(sample.names,"_R_filt.fastq.gz"))
out <- filterAndTrim(F.seqs, filt.F.seqs, R.seqs, filt.R.seqs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

#learn the error rates.
err.F <- learnErrors(filt.F.seqs, multithread = T)
err.R <- learnErrors(filt.R.seqs, multithread = T)

#dereplicate identical sequences to unique sequences with corresponding abundances.
derep.F <- derepFastq(filt.F.seqs, verbose = T)
derep.R <- derepFastq(filt.R.seqs, verbose = T)
names(derep.F) <- sample.names
names(derep.R) <- sample.names

#apply dada2 inference algorithm. whatever.
dada.F <- dada(derep.F, err=err.F, multithread = T)
dada.R <- dada(derep.R, err=err.R, multithread = T)

#merge paired reads. This generates a list of unique sequences and counts for each sample.
mergers <- mergePairs(dada.F, derep.F, dada.R, derep.R, verbose = T)

#make a ASV table.
#columns are unique sequences. Column names are the actual unique sequences.
#rows are unique samples. rownames are the acutal sample names.
seqtab <- makeSequenceTable(mergers)
output.path <- paste0(cur.dir,'/seqtab.rds')
saveRDS(seqtab,output.path)

#remove chimeras
seqtab <- readRDS(output.path)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
