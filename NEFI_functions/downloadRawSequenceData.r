# Function downloads raw sequence data and metadata for NEON ITS and 16S sequencing
# Written by Lee Stanish, NEON, 2019
# Doesn't download provisional raw sequence data (2013-2014)
# Provisional data was obtained directly from Lee Stanish

downloadRawSequenceData <- function(site="all", startdate="YYYY-MM", enddate="YYYY-MM", 
                                    outdir="", check.size=TRUE, taxa=NULL) {
  # outdir - path to directory to ouput the data. Defaults to the R default directory if none provided
  # change checkFileSize to FALSE to override file size checks
  
  library(neonUtilities)
  library(utils)
  
  dat <- loadByProduct(dpID="DP1.10108.001", site=site, startdate=startdate, enddate=enddate,
                       package="expanded", check.size=check.size)
  
  u.urls <- unique(dat$mmg_soilRawDataFiles$rawDataFilePath)
  
  if (!is.null(taxa)) {
    if (taxa == "16S") {
      u.urls <- u.urls[grep("16S", u.urls)]
    } else if (taxa == "ITS") {
      u.urls <- u.urls[grep("ITS", u.urls)]
    }
  }
  
  fileNms <- gsub('^.*\\/', "", u.urls)
  print(paste("There are", length(u.urls), " unique files to download.") )
  
  for(i in 1:length(u.urls)) {
    download.file(url=as.character(u.urls[i]), destfile = ifelse(dir.exists(outdir), 
                                                                 paste(outdir, fileNms[i], sep="/"), 
                                                                 paste(getwd(), fileNms[i], sep="/" )) ) 
  }
  
  return(dat)
  ## END FUNCTION ##
}