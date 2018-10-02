#' Set truncation length
#' 
#' Decides on truncation length for trimmed reads based on quality score means. Default cutoff is a score of 30. 
#' Warns you if the beginning (first 10) bases are low-quality, but returns the first low-quality base after the 10th base.
#' If no bases are below score, returns the last base.
#' 
#' @param fl input files (prints suggested length for each; if used in script, could take the first input, or the lowest)
#' @param qscore default = 30.
#' @param n default = 5e+05 don't know exactly why this matters, but it's part of how 'qa' is called within the dada2 scripts...
#' 
#' @export
#' 
#' @examples

get_truncation_length <-function (fl, qscore = 30, n = 5e+05, quiet = TRUE){
  trunc_lengths <- data.frame(file = character(0), early_lowqual = numeric(0),
                              trunc_length = numeric(0))
  for (f in fl) {
    srqa <- ShortRead::qa(f, n = n)
    df <- srqa[["perCycle"]]$quality
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count, 
                                                          df$Cycle)
    lowqual <- which(means < qscore)
    lowqual_in_first_10 <- length(which(lowqual <= 10))
    lowqual_after_10 <- lowqual[lowqual > 10][1]
    
    trunc_lengths <- rbind(trunc_lengths, data.frame(file = f, early_lowqual = lowqual_in_first_10,
                                                    trunc_length = ifelse(is.na(lowqual_after_10), length(means), lowqual_after_10))) 

if (quiet == FALSE) {
  if(lowqual_in_first_10 > 0){
    cat(paste0(basename(f),': ', lowqual_in_first_10, ' bases within first 10 bases are below your quality score. Consider trimming left.\n'))
    } 
  if (is.na(lowqual_after_10)){
    cat(paste0(basename(f), ': After first 10 bases, no bases with a mean under your quality score found. Truncate at end of read, base: ', length(means),'\n'))
  
    } else if (!is.na(lowqual_after_10)){
    cat(paste0(basename(f),': After first 10 bases, first mean below your quality score is base: ',lowqual_after_10,'\n'))  
  
    } else "Something's wrong. Inspect function/reads."
} # end printing to console  
} # end loop
return(trunc_lengths$trunc_length)  
} # end function
