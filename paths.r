#all paths for NEFI fungal amplicon analysis and forecasts.
#The path to the data directory will depend if you are on your local computer, the scc or somewhere else.
#To deal with this Colin has condiitonally setup the path to data based on the hostname of the computer.
#Forst instance, if the hostname is 'pecan2' colin tells the omputer that all the data are in /fs/data3/caverill/NEFI_16S_data/
#You can your local or remote machine to this list, following the instructions below. You just need your hostname.
#To get your hostname open up terminal and type 'hostname'.
#I often save data objects as ".rds" files, which is just an R data storage type.
#.rds files can be loaded with "readRDS()", and saved with "saveRDS(R_object, path/to/file.rds)".
host <- system('hostname', intern=T)
#data directory conditional to which computer you are working on.
#defaults to scc directory.
data.dir <- '/projectnb/talbot-lab-data/caverill/NEFI_microbial_data/'
#conditional data directory assignment.
if(host == 'pecan2'){data.dir <- '/fs/data3/caverill/NEFI_microbial_data/'}
#if(host == 'scc1'  ){data.dir <- '/projectnb/talbot-lab-data/caverill/NEFI_16S_data/'} #hashed out as I haven't set this up on scc yet.
#if(host == 'scc2'  ){data.dir <- '/projectnb/talbot-lab-data/caverill/NEFI_16S_data/'}
#if(host == 'Colins-MacBook-Pro.local'){data.dir <- '/Users/colin/Documents/rstudio_data/NEFI_16S_data'}
#make directory if it doesn't exist
cmd <- paste0('mkdir -p ',data.dir)
system(cmd)

#JAGS output.
dir <- paste0(data.dir,'JAGS_output/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
ted_ITS.prior_fg_JAGSfit    <- paste0(dir,'ted_ITS.prior_fg_JAGSfit.rds')
ted_ITS.prior_20gen_JAGSfit <- paste0(dir,'ted_ITS.prior_20gen_JAGSfit.rds')
ITS.prior_linear_fg_cov.selection_JAGS <- paste0(dir,'ITS.prior_linear_fg_cov.selection_JAGS.rds')

#Bahram sequence data.
bahram_dir <- paste0(data.dir,'bahram_2018_seqs')
cmd <- paste0('mkdir -p ',bahram_dir)
system(cmd)

#Raw NEON sequence data from MG-rast
mg_rast.key <- paste0(data.dir,'reference_data/MG-RAST_mapped_identifiers.csv')
#setup place to save ITS data.
ASV_ITS.dir <- paste0(data.dir,'NEON_ASVs/ITS/')
ASV_16S.dir <- paste0(data.dir,'NEON_ASVs/16S/')
cmd <- paste0('mkdir -p ',ASV_ITS.dir)
cmd <- paste0('mkdir -p ',ASV_16S.dir)


#Tedersoo ITS prior paths
#raw sequence directory.
ted.seq.dir <- paste0(data.dir,'tedersoo_SRA_seqs/')

dir <- paste0(data.dir,'prior_data/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
ted.ITSprior_data <- paste0(dir,'ted_all_prior_data.rds')
      ted_map_raw <- "/fs/data3/caverill/Microbial_Space_Time_data/tedersoo_2014.data/merging seq and site data.csv"
      ted_otu_raw <- "/fs/data3/caverill/Microbial_Space_Time_data/tedersoo_2014.data/alldata.biom_rdp_tedersoo_otu.txt"
      ted_srr_codes <- paste0(dir,'tedersoo_srr_codes.csv')
      
#downloading Tedersoo sequence data from SRA.
ted2014_seq_dir <- paste0(data.dir,'tedersoo_SRA_seqs')
cmd <- paste0('mkdir -p ',ted2014_seq_dir)
system(cmd)


#product ouput paths.
dir <- paste0(data.dir,'NEON_data_aggregation/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
         ITS_site_dates.path <- paste0(dir,'ITS_site_dates.rds')
        site_level_data.path <- paste0(dir,'site_level_data.rds')
    dp1.10108.00_output.path <- paste0(dir,'dp1.10108.00_output.rds')
    dp1.10109.00_output.path <- paste0(dir,'dp1.10109.00_output.rds')
    dp1.10086.00_output.path <- paste0(dir,'dp1.10086.00_output.rds')
    dp1.10078.00_output.path <- paste0(dir,'dp1.10078.00_output.rds')
    dp1.10098.00_output.path <- paste0(dir,'dp1.10098.00_output.rds')
dp1.10098.00_plot.level.path <- paste0(dir,'dp1.10098.00_plot.level.rds')

#Paths to spatial x and y variables for modeling.
 obs.table.path <- paste0(dir, 'obs.table.rds') 
core.table.path <- paste0(dir,'core.table.rds') 
plot.table.path <- paste0(dir,'plot.table.rds') 
site.table.path <- paste0(dir,'site.table.rds') 

#Get NEON map, tax and OTU tables.
#This wille ventually be repalced by ASV and tax tables genereated by CA pipeline, which will pull and process raw data NEON sequence data from MG-RAST.
#mapping file will be constructed from DP1.10801.001 product, linking it to core and site data.
dir <- paste0(data.dir,'map_otu/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
its_otu.path <- paste0(dir,'ITS_otu_clean.rds')
its_map.path <- paste0(dir,'ITS_map_clean.rds')
its_tax.path <- paste0(dir,'ITS_tax_clean.rds')
its_fun.path <- paste0(dir,'ITS_fun_clean.rds')

#reference data product paths
dir <- paste0(data.dir,'reference_data/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
              em_traits.path <- paste0(dir,'ecto_genus_traits_hobbie_Jan2018.csv')
             em_species.path <- paste0(dir,'myc_assignments.rds')
             poa_genera.path <- paste0(dir,'poaceae_genera_wikipedia.rds')
            #em_species.path <- paste0(dir,'mycorrhizal_SPCD_data.csv')
              em_genera.path <- paste0(dir,'tedersoo_2017_genera.csv')
 NEON_plantStatus_codes.path <- paste0(dir,'NEON_DP1.10098.plantStatus_decode.csv')
 

#Forecast covariates
#Aggregated NEON site level covariates and global level uncertainty for predictors
dir <- paste0(data.dir,'NEON_covariates/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)

NEON_site_covs.path <- paste0(dir,'NEON_site_covs.rds')
NEON_glob_covs.path <- paste0(dir,'NEON_glob_covs.rds')