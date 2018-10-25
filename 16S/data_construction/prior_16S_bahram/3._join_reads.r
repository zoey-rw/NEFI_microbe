# Joins forward and reverse reads using fastq-join.
# input is a folder of forward and reverse reads, output is a new folder with joined reads
rm(list=ls())
source('paths.r')

# Input path to sequences:
#seq.dir <- bahram.seq.dir
seq.dir <-  "/projectnb/talbot-lab-data/NEFI_data/big_data/bahram_test/"

# Path to trimmed sequences
trimmed_seqs <- paste0(seq.dir, "q.trim/")

# Join forward and reverse reads, output into new joined_reads directory
cmd <- paste0("multiple_join_paired_ends.py -i ", trimmed_seqs, " -o ", seq.dir, "joined_seqs --read1_indicator '_1' --read2_indicator '_2'")
system(cmd)

# Import script to clean up the output of fastq-join; goes into each folder, extracts the joined file, renames as the sample name, and puts it in the enclosing directory ("joined_seqs")
cmd <- paste0("cp /usr3/graduate/zrwerbin/NEFI_microbe/16S/data_construction/prior_16S_bahram/condenseJoinedFiles.sh ", seq.dir, "joined_seqs/")
system(cmd)

# execute that script
cmd <- paste0("cd ", seq.dir, "joined_seqs/", ";", seq.dir, "joined_seqs/", "condenseJoinedFiles.sh")
system(cmd)


# removes everything other than the .fastq files

cmd <- paste0("cd ", seq.dir, "joined_seqs/", ";", 'find . ! -name "*.fastq" -exec rm -r {} \\;')
system(cmd)
