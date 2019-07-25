#OTU table can be loaded into R now, its sufficiently small.
#This assigns taxonomy in parallel using dada2.

#clear environment, source paths.
rm(list=ls())
source('paths.r')
library(doParallel)

#how many cores to run on and therefore how many subsets to break taxonomy string into.
n <- 16 #important to request a node with 16 (or fewer) cores.
registerDoParallel(cores=n)

#load otu table where rownames are unique sequences.
otu <- readRDS(emp_esv_clean.path)
to_assign <- rownames(otu)

#set breakpoints for subsetting taxonomy list.
brk <- round(length(to_assign) / n)

#foreach loops are 
output.list <-
foreach(i = 1:n) %dopar% {
  #tell loop where i of n taxonomy subset starts and ends.
  start <- (brk*i - brk) + 1
    end <- brk*i
    #if you on the last subset go to end.
    if(i == n){end = length(to_assign)}
    
    #assign taxa
    tax.out <- dada2::assignTaxonomy(to_assign[start:end],greengenes.path)
    #return output to list
    return(tax.out)
}

#merge together output of parallel assignment.
out <- data.frame(do.call('rbind',output.list))

#save output as your taxonomy file.
saveRDS(out, emp_tax.path)
