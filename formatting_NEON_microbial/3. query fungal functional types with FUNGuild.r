#Code to assign fungal guilds using the FUNGuild database
#clear environment, load packages.
rm(list=ls())
library(rvest)
library(jsonlite)

#set taxonomy path and output path
#Colin is setting it to be the same, just appending a column to my taxonomy file.
tax.path <- '/fs/data3/caverill/NEFI_microbial/map_otu/ITS_tax_clean.rds'
out.path <- '/fs/data3/caverill/NEFI_microbial/map_otu/ITS_tax_clean.rds'

#path to FUNGuild database
url <- "http://www.stbates.org/funguild_db.php"

#download database, convert it to something R interpretable.
fg <- url %>% 
  read_html() %>%
  html_text() 
fg <- jsonlite::fromJSON(gsub("funguild_db", "", fg))

#load your fungal taxonomy data.
tax <- readRDS(tax.path)

#first, I need to remove leading `n__` from each of my taxonomic vectors.
for(i in 1:7){
tax[,i] <- substring(tax[,i],4)
}

#There are 9 unique levels of taxonomic resolution actually in FUNGuild (though 24 potential levels)
#0-keyword, 3-Phylum, 7-Order, 9-Family, 13-genus, 20species, 21-subspecies, 24-Form
#I only have data on k/c/p/o/f/g/s, so only going to deal with levels 3,7,9,13,20
#basically write a series of if statements to assign function
#start with highest level of taxonomy and go down.

#series of ifelse statements to match each level.
#there is almost certainly a better and faster way to do this. But this works and is fast enough.
tax$guild <- NA
#phylum level match.
tax$guild <- ifelse(tax$taxonomy2 %in% fg$taxon,
                    fg[match(tax$taxonomy2, fg$taxon),5],
                    tax$guild)
#class level match.
tax$guild <- ifelse(tax$taxonomy3 %in% fg$taxon,
                    fg[match(tax$taxonomy3, fg$taxon),5],
                    tax$guild)
#order level match.
tax$guild <- ifelse(tax$taxonomy4 %in% fg$taxon,
                    fg[match(tax$taxonomy4, fg$taxon),5],
                    tax$guild)
#family level match.
tax$guild <- ifelse(tax$taxonomy5 %in% fg$taxon,
                    fg[match(tax$taxonomy5, fg$taxon),5],
                    tax$guild)
#genus level match.
tax$guild <- ifelse(tax$taxonomy6 %in% fg$taxon,
                    fg[match(tax$taxonomy6, fg$taxon),5],
                    tax$guild)

#species level match.
tax$guild <- ifelse(tax$taxonomy7 %in% fg$taxon,
                    fg[match(tax$taxonomy7, fg$taxon),5],
                    tax$guild)

#save output.
saveRDS(tax, out.path)

#print fraction of fungal taxa assigned to a functional guild.
cat(sum(!is.na(tax$guild))/(nrow(tax))*100,'% of fungal taxa assigned a functional guild.', sep = '')