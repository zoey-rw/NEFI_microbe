#downloading a NEON site-month set of sequences, process into ASV tables.
#This merely gets unique sequences from the already paired ends and quality filtered files available.
#The format is compatible with downstream chimera checking and RDP assignment with dada2.
#This downloads and builds unqiue ASV tables for each site-date combination.
#There are functions in dada2 to merge ASV tables.
source('paths.r')

#ITS ASV master path. Where all site by date ASV tables are stored.
main.path <- ASV_ITS.dir

#Load NEON data that we know has sequence data
dna_ids <- readRDS(dp1.10108.00_output.path)
site_dates <- readRDS(ITS_site_dates.path)
#this mg-rast file was passed to CA by Lee Stanish so we could rip from MG-RAST.
mg_rast_link <- read.csv(mg_rast.key)
#merge these together.
dna_ids <- merge(dna_ids, mg_rast_link[,c('dnaSampleID','metagenome_id')], all.x = T)
dna_ids$metagenome_id <- as.character(dna_ids$metagenome_id)


#subset to a single site for testing.
site_dates <- site_dates[1]

#test function after loading the function below.
test <- rip_NEON_seqs(site_dates, dna_ids)


rip_NEON_seqs <- function(site_dates, dna_ids){
  for(i in 1:length(site_dates)){
    site <- names(site_dates)[i]
    for(k in 1:length(site_dates[[site]])){
      ##### Specify date, check if you already have this data, if not. go for it. ####
      date <- site_dates[[site]][k]
      
      #make site-date directory if it doesnt already exist.
      #if upstream site directory doesnt exist, it will also be created.
      cur.dir <- paste0(main.path,site,'/',date)
      pre <- paste0('mkdir -p ', cur.dir)
      system(pre)
      
      #grab mgm_ids for a particular site_date
      mgm.id <- dna_ids[dna_ids$siteID == site & dna_ids$dateID == date,]$metagenome_id
      mgm.id <- mgm.id[!is.na(mgm.id)]
      
      #output path for final asv table.
      asv.output.path <- paste0(cur.dir,'/ASV_',site,'.',date,'.rds')
      if(file.exists(asv.output.path) == T){
        cat(paste0(site,' ',date,' ASV table already present. you good! moving on.\n'))
        next
      }
      #### Download the sequences. ####
      cat(paste0('Downloading sequences for ',site,' ',date,'...\n'))
      for(i in 1:length(mgm.id)){
        mgrast.link <- paste0('http://api.metagenomics.anl.gov/1/download/',mgm.id[i],'?file=050.1')
        file <- paste0(cur.dir,'/',mgm.id[i],'.fasta')
        pre <- paste0('curl ',mgrast.link,' > ',file)
        system(pre)
      }
      cat('Download complete.\n')
      
      #get names of all sequence files in cur.dir.
      files <- sort(list.files(cur.dir, pattern = '.fasta', full.names = T))
      sample.names <- sapply(strsplit(basename(files), ".fasta"), `[`, 1)
      
      #### Build an ASV table from all sequence files. ####
      cat(paste0('Building ASV table for ',site,' ',date,'...\n'))
      for(j in 1:length(sample.names)){
        #get just the sequences into a separate file.
        file <- paste0(cur.dir,'/'    ,sample.names[j],'.fasta')
        s.file <- paste0(cur.dir,'/seq.',sample.names[j],'.fasta')
        pre <- paste0("sed -n '0~2p' ",file,' > ',s.file)
        system(pre)
        
        #Convert sequences to a sequence table. Write table to a file. No need to load all sequences into R memory.
        c.file <- paste0(cur.dir,'/counts.',sample.names[j],'.txt')   #sample-specific ASV table.
        pre <- paste0('cat ',s.file,' | sort | uniq -c > ', c.file)
        system(pre)
        asv <- data.table::data.table(read.table(c.file))
        colnames(asv) <- c(sample.names[j],'seq')
        data.table::setkey(asv,seq)
        
        #merge multiple sample-specific ASV tables.
        if(j == 1){ out = asv}
        if(j > 1){ out <- merge(out, asv, all = T)} #this merge really requires data.table be loaded into the environment.
        
        #remove duplicated seq and counts files.
        system(paste0('rm ',c.file))
        system(paste0('rm ',s.file))
      }
      #convert back to dataframe, replace NA values with zeros.
      out <- as.data.frame(out)
      out[is.na(out)] <- 0
      
      #transpose to be consistent with dada2
      #this is being weird.
      t.out <- t(out[,-1])
      colnames(t.out) <- out$seq
      
      #convert from numeric dataframe to integer matrix. this is important for dada2 commands downstream.
      t.out <- as.matrix(t.out)
      t.out <- apply (t.out, c (1, 2), function (x) {(as.integer(x))})
      
      #save ASV table, remove downloaded sequence data.
      saveRDS(t.out,asv.output.path)
      system(paste0('rm ',cur.dir,'/*.fasta'))
      
      cat(paste0(site,' ',date,' ASV table constructed.\n'))
      
    } #end date loop.
  } #end site loop.
} #end function loop.

