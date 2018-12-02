#' common_group_quantification.r
#' Supply a sv/otu table, a taxonomy file, and a phylogenetic level.
#' Return sample abundances, relative abundances and sequence depth for all taxa at that level 
#' found in greater than X% of samples (default 50%).
#'
#' @param sv           #sv or otu table where columns are unique taxonomic assignments, rows are samples.
#' @param tax          #taxonomy for sv/otu table. 7 columns for k/p/c/o/f/g/s.
#' @param groups       #unique groups you are looking for. ex: c('Ascomycota','Basidiomycota')
#' @param phyla_level  #level of phylogeny in quotes. i.e. 'phylum'.
#'
#' @return
#' @export
#'
#' @examples
common_group_quantification <- function(sv, tax, groups, tax_level, samp_freq = 0.5){
  #some tests.
  if(ncol(sv) != nrow(tax)){
    stop('Number of columns of sv table does not match number of rows in taxonomy table.')
  }
  tax_level <- tolower(tax_level)
  #make sure taxonomy column names are lower case.
  tax_names <- c('kingdom','phylum','class','order','family','genus','species')
  colnames(tax) <- tax_names
  
  #merge taxonomy and sv file.
  k <- cbind(tax, t(sv))
  k <- data.table(k)
  to_change <- colnames(k)[colnames(k) %in% tax_level]
  setnames(k,to_change,'tax_level')
  #count sequence abundance in each group by sample, and number of samples in which it occurs.
  abundance <- list()
  frequency <- list()
  for(i in 1:length(groups)){
    z <- k[tax_level == groups[i],]
    start <- ncol(tax) + 1
    seq.count <- colSums(z[,start:ncol(z)])
    seq.freq  <- length(seq.count[seq.count > 0]) / nrow(sv)
    abundance[[i]] <- seq.count
    frequency[[i]] <- seq.freq
  }
  abundance <- do.call(rbind,abundance)
  frequency <- do.call(rbind,frequency)
  
  #name some stuff.
  rownames(abundance) <- groups
  colnames(abundance) <- rownames(sv)
  
  #put together frequency table and order.
  frequency <- data.frame(cbind(groups,frequency))
  colnames(frequency)[2] <- 'sample_frequency'
  frequency$sample_frequency <- as.character(frequency$sample_frequency)
  frequency$sample_frequency <- as.numeric(frequency$sample_frequency)
  frequency <- frequency[order(-frequency$sample_frequency),]
  #subset to those that are found in > sam_freq of samples (default 50%).
  ref.frequency <- frequency[frequency$sample_frequency > samp_freq,]
  #kill unknown and anything unknown
  ref.frequency <- ref.frequency[!(ref.frequency$groups %in% c('unknown','Unknown')),]
  
  #get sequence depth per sample.
  seq_total <- rowSums(sv)
  
  #only keep groups that are found in >samp_freq% of samples.
  abundance <- abundance[rownames(abundance) %in% ref.frequency$groups,]
  abundance <- t(abundance)
  
  #get abundances and relative abundances.
          other <- seq_total - rowSums(abundance)
      abundance <- cbind(other, abundance)
  rel.abundance <- abundance / colSums(abundance)
  
  #return output: abundances, relative abundances, sequence depth, group sample frequencies.
         to_return  <- list(abundance,rel.abundance,seq_total,frequency)
   names(to_return) <- c('abundances','rel.abundances','seq_total','group_frequencies')
  return(to_return)
}