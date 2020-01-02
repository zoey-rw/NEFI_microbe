#' fg_assign_parallel.r
#' assigns functional groups to a taxonomy table based on the FUNGuild database.
#' This function was built by Colin Averill, however the FUNGuild database was developed by Nguyen et al.
#' PLEASE CITE: Nguyen NH, Song Z, Bates ST, Branco S, Tedersoo L, Menke J, Schilling JS, Kennedy PG. 2016. FUNGuild: An open annotation tool for parsing fungal community datasets by ecological guild. Fungal Ecology, 20: 241-248. doi:10.1016/j.funeco.2015.06.006
#' depends on R packages rvest, jsonlite.
#' I load this function at the top of a script using: source("/path/to/fg_assign.r")
#'
#' @param tax_table  #dataframe of taxonomy to assign. kingdom/phylum/class/order/family/genus/species as column names.
#' @param url        #URL where FUNGuild database lives. currently default to: http://www.stbates.org/funguild_db.php.
#' @param n.cores    #number of processing cores. If not specified script will automatically detect the number of cores avaialble and use all of them.
#'
#' @return           #returns the taxonomy table with FUNGuild assignments appended.
#' @export
#'
#' @examples
#' Testing the function - generate some artificial data to assign guild to.
#' tax_table <- structure(list(kingdom = "Fungi", phylum = "Ascomycota", class = "Pezizomycetes", order = "Pezizales", family = "Tuberaceae", genus = "Tuber", species = "Tuber melosporum"), 
#'                             .Names = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), 
#'                             row.names = 4L, class = "data.frame")
#' k <- tax_table
#' tax_table <- rbind(k,k,k,k,k,k,k,c('Fungi',NA,NA,NA,NA,NA,NA))
#' test <- fg_assign(tax_table)
#' 
fg_assign_parallel <- function(tax_table, url = "http://www.stbates.org/funguild_db.php", n.cores = NA){
  #check if dependencies are installed. If not, stop.----
  if (!require('rvest'   ,character.only = TRUE)){
    stop("please install the rvest package.")
  }
  if (!require('jsonlite',character.only = TRUE)){
    stop("please install the jsonlite package.")
  }
  if(!require('doParallel', character.only = TRUE)){
    stop('please install the doParallel package.')
  }
  
  #check that the input is formatted right. If not, stop, throw an error.
  if (!is.data.frame(tax_table)){
    stop('Your taxonomy table needs to be a data.frame. Try again.')
  }
  
  #make sure column names are lower case.
  colnames(tax_table) <- tolower(colnames(tax_table))
  
  #setup parallel.----
    library(doParallel)
    if(is.na(n.cores)){
      n.cores <- detectCores()
    }
    #if number of cores is still NA for some reason set to 1.
    if(is.na(n.cores)){
      n.cores <- 1
      cat('detectCores() returned NA. Setting n.cores to 1.\n')
    }
    registerDoParallel(n.cores)

  #download FUNGuild database, convert it to something R interpretable.----
  fg <- url %>% 
    xml2::read_html() %>%
    rvest::html_text() 
  fg <- jsonlite::fromJSON(gsub("funguild_db", "", fg))
  
  #There are 9 unique levels of taxonomic resolution actually in FUNGuild (though 24 potential levels)
  #0-keyword, 3-Phylum, 7-Order, 9-Family, 13-genus, 20-Species, 21-Subspecies, 24-Form
  #This function requires data on k/c/p/o/f/g/s, so only deals with levels 3,7,9,13,20
  #What follows is a series of if statements to assign function within a parallel loop.
  #start with highest level of taxonomy and go down.

  fg$taxon <- tolower(fg$taxon)
  #setup output list.
  out <- list()
  #run parallel fg assign loop.
  out <-
  foreach(i = 1:nrow(tax_table)) %dopar% {
    to_return <- NA
    #phylum level match.
    if(tax_table$phylum[i] %in% fg$taxon){
      to_return <- fg[match(tax_table$phylum[i], fg$taxon),4:10]
    }
    #class level match.
    if(tax_table$class[i] %in% fg$taxon){
      to_return <- fg[match(tax_table$class[i], fg$taxon),4:10]
    }
    #order level match.
    if(tax_table$order[i] %in% fg$taxon){
      to_return <- fg[match(tax_table$order[i], fg$taxon),4:10]
    }
    #family level match.
    if(tax_table$family[i] %in% fg$taxon){
      to_return <- fg[match(tax_table$family[i], fg$taxon),4:10]
    }
    #genus level match.
    if(tax_table$genus[i] %in% fg$taxon){
      to_return <- fg[match(tax_table$genus[i], fg$taxon),4:10]
    }
    #species level match.
    if(tax_table$species[i] %in% fg$taxon){
      to_return <- fg[match(tax_table$species[i], fg$taxon),4:10]
    }
    #return output.
    return(to_return)
  } #end parallel loop.
  
  #bind up output, append to taxonomy table.
  out <- do.call(rbind,out)
  tax_table <- cbind(tax_table, out)

  #report and return output.
  cat(sum(!is.na(tax_table$guild))/(nrow(tax_table))*100,'% of taxa assigned a functional guild.\n', sep = '')
  return(tax_table)
} #end function.
