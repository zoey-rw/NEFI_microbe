#all paths for NEFI fungal amplicon analysis and forecasts.
#The path to the data directory will depend if you are on pecan2 or the scc.
#To deal with this Colin has condiitonally setup the path to data based on the hostname of the computer.
#Forst instance, if the hostname is 'pecan2' colin tells the computer that all the data are in /fs/data3/caverill/NEFI_data
#I often save data objects as ".rds" files, which is just an R data storage type.
#.rds files can be loaded with "readRDS()", and saved with "saveRDS(R_object, path/to/file.rds)".

#### High level directory structure. ####
#NEFI_data - master data directory.
host <- system('hostname', intern=T)
#data directory conditional to which computer you are working on.
data.dir <- '/projectnb/talbot-lab-data/NEFI_data/'
#conditional data directory assignment.
if(host == 'pecan2'){data.dir <- '/fs/data3/caverill/NEFI_data/'}
#make directory if it doesn't exist
cmd <- paste0('mkdir -p ',data.dir)
system(cmd)

#NEFI_data has 3 sub-directories. big_data, ITS and 16S.
#Each ITS and 16S directory has sub directories for data generated in scc or pecan (scc_gen or pecan_gen)
#big_data only lives in the scc.
#scc_gen and pecan_gen live on both, but updating the directory from one computer to the other only goes one way.
#scc_gen only goes scc->pecan. pecan_gen only goes pecan->scc. Setup is in an rsync.sh script.
#This keeps data synced and nothing goes missing. There are transfer scripts that make sure this is done correctly.
big_data_dir <- paste0(data.dir,'big_data/')
data_ITS_dir <- paste0(data.dir,'ITS/')
data_16S_dir <- paste0(data.dir,'16S/')
  ITS_scc_gen_dir <- paste0(data_ITS_dir,  'scc_gen/')
ITS_pecan_gen_dir <- paste0(data_ITS_dir,'pecan_gen/')
  scc_gen_16S_dir <- paste0(data_16S_dir,  'scc_gen/')
pecan_gen_16S_dir <- paste0(data_16S_dir,'pecan_gen/')

#make the directories if they don't exist.
system(paste0('mkdir -p ',big_data_dir))
system(paste0('mkdir -p ',data_ITS_dir))
system(paste0('mkdir -p ',  ITS_scc_gen_dir))
system(paste0('mkdir -p ',ITS_pecan_gen_dir))
system(paste0('mkdir -p ', scc_gen_16S_dir))
system(paste0('mkdir -p ', pecan_gen_16S_dir))


#### big_data file paths. ####
#Bahram 2018 sequences.
bahram.seq.dir <- paste0(big_data_dir,'bahram_2018_seqs/')
cmd <- paste0('mkdir -p ',bahram.seq.dir)
system(cmd)
# NEON 16S sequences.
NEON.seq.dir <- paste0(big_data_dir,'NEON_raw_16S_fastq/')
cmd <- paste0('mkdir -p ',NEON.seq.dir)
system(cmd)

#Tedersoo ITS prior paths
#raw sequence directory.
ted.seq.dir <- paste0(big_data_dir,'tedersoo_SRA_seqs/')
cmd <- paste0('mkdir -p ',ted.seq.dir)
system(cmd)

#Raw NEON ITS sequence data from custom links provided by L.Stanish.
      NEON_ITS.dir <- paste0(big_data_dir,'NEON_raw_ITS_seqs/')
NEON_ITS_link_file <- paste0(NEON_ITS.dir,'NEON_rawFilesList.csv')
  NEON_ITS_SV.path <- paste0(NEON_ITS.dir,'NEON_ITS_sv.rds')
 NEON_ITS_tax.path <- paste0(NEON_ITS.dir,'NEON_ITS_tax.rds')

#Raw NEON sequence data from MG-rast. only works for 16S currently.
mg_rast.key <- paste0(data.dir,'reference_data/MG-RAST_mapped_identifiers.csv')
#setup place to save MG-RAST sequence data.
ASV_ITS.dir <- paste0(big_data_dir,'NEON_mgRAST_ASVs/ITS/')
ASV_16S.dir <- paste0(big_data_dir,'NEON_mgRAST_ASVs/16S/')
cmd <- paste0('mkdir -p ',ASV_ITS.dir)
system(cmd)
cmd <- paste0('mkdir -p ',ASV_16S.dir)
system(cmd)


#### ITS/scc_gen: JAGS output paths. ####
dir <- paste0(ITS_scc_gen_dir,'JAGS_output/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
ted_ITS.prior_fg_JAGSfit    <- paste0(dir,'ted_ITS.prior_fg_JAGSfit.rds')
ted_ITS.prior_20gen_JAGSfit <- paste0(dir,'ted_ITS.prior_20gen_JAGSfit.rds')
ITS.prior_linear_fg_cov.selection_JAGS <- paste0(dir,'ITS.prior_linear_fg_cov.selection_JAGS.rds')

#### ITS/scc_gen: tedersoo 2014 SV and taxonomy paths. ####
ted_2014_SV.table.path <- paste0(ITS_scc_gen_dir,'ted_2014_SV.table.rds')
     ted_2014_tax.path <- paste0(ITS_scc_gen_dir,'ted_2014_tax.rds') 
     
#### ITS/scc_gen: NEON to forecast SV and taxonomy paths. ####
NEON_SV.table.path <- paste0(ITS_scc_gen_dir,'NEON_SV.table.rds')
     NEON_tax.path <- paste0(ITS_scc_gen_dir,'NEON_tax.rds') 
     
#### ITS/pecan_gen: Uncertainty product paths. ####
dir <- paste0(ITS_pecan_gen_dir,'uncertainty_products/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
#worldclim2
wc2.dir <- paste0(dir,'worldclim2_uncertainty/')
cmd <- paste0('mkdir -p ',wc2.dir)
system(cmd)
wc_prec_raw_data.path <- paste0(wc2.dir,'worldclim2_prec_raw_data.rds')
wc_temp_raw_data.path <- paste0(wc2.dir,'worldclim2_temp_raw_data.rds')
wc_prec_JAGS.path     <- paste0(wc2.dir,'precipitation_JAGS_model.rds')
wc_temp_JAGS.path     <- paste0(wc2.dir,'temperature_JAGS_model.rds')

#NEON soil uncertainty products.
neon_uncertainty_dir <- paste0(dir,'NEON_soil_products/')
cmd <- paste0('mkdir -p ',neon_uncertainty_dir)
system(cmd)
pC_NEON_uncertainty_model.path <- paste0(neon_uncertainty_dir,'pC_NEON_uncertainty_model.rds')
cn_NEON_uncertainty_model.path <- paste0(neon_uncertainty_dir,'cn_NEON_uncertainty_model.rds')

#### ITS/pecan_gen: Tedersoo 2014 prior data. ####
dir <- paste0(ITS_pecan_gen_dir,'tedersoo_2014_data/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
tedersoo_ITS.prior_for_analysis.path <- paste0(dir,'tedersoo_ITS.prior_for_analysis.rds')
tedersoo_ITS.prior_fromSV_analysis.path <- paste0(dir,'tedersoo_ITS.prior_fromSV_analysis.rds')

#### ITS/pecan_gen: NEON data aggregation paths. ####
dir <- paste0(ITS_pecan_gen_dir,'NEON_data_aggregation/')
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

#paths for hierarchical co-variate aggregation.
core_obs.path  <- paste0(dir,'core_obs.rds')
core_core.path <- paste0(dir,'core_core.rds')
core_plot.path <- paste0(dir,'core_plot.rds')
core_site.path <- paste0(dir,'core_site.rds')
core_glob.path <- paste0(dir,'core_glob.rds')
plot_plot.path <- paste0(dir,'plot_plot.rds')
plot_site.path <- paste0(dir,'plot_site.rds')
plot_glob.path <- paste0(dir,'plot_glob.rds')
site_site.path <- paste0(dir,'site_site.rds')
site_glob.path <- paste0(dir,'site_glob.rds')
hierarch_filled.path <- paste0(dir,'hierarch_filled.path')

#### ITS/pecan_gen: Forecast covariate paths.----
#Aggregated NEON site level covariates and global level uncertainty for predictors
dir <- paste0(ITS_pecan_gen_dir,'NEON_covariates/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
NEON_site_covs.path <- paste0(dir,'NEON_site_covs.rds')
NEON_glob_covs.path <- paste0(dir,'NEON_glob_covs.rds')

#### ITS/pecan_gen: Forecast output paths.----
dir <- paste0(ITS_pecan_gen_dir,'NEON_forecasts/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
NEON_site_fcast_fg.path <- paste0(dir,'NEON_fcast_site.level_fg.rds')
NEON_site_fcast_20gen.path <- paste0(dir,'NEON_fcast_site.level_20gen.rds')

#### ITS/pecan_gen: NEON map, tax and OTU table paths. ####
#This will eventually be repalced by ASV and tax tables genereated by CA pipeline, which will pull and process raw data NEON sequence data from MG-RAST.
#mapping file will be constructed from DP1.10801.001 product, linking it to core and site data.
dir <- paste0(ITS_pecan_gen_dir,'NEON_map_otu/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
its_otu.path <- paste0(dir,'NEON_ITS_otu_clean.rds')
its_map.path <- paste0(dir,'NEON_ITS_map_clean.rds')
its_tax.path <- paste0(dir,'NEON_ITS_tax_clean.rds')
its_fun.path <- paste0(dir,'NEON_ITS_fun_clean.rds')

#### ITS/pecan_gen: Figure paths.
#all figures generated from project.
dir <- paste0(ITS_pecan_gen_dir,'figures/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
         ecto_prior_prediction.fig.path <- paste0(dir,'ecto_prior_prediction.png')
           fg_prior_prediction.fig.path <- paste0(dir,'fg_prior_prediction.png')
        top20_prior_prediction.fig.path <- paste0(dir,'top20_prior_prediction.png')
NEON_site.level_top20_forecast.fig.path <- paste0(dir,'NEON_site.level_20gen_forecast.png')
   NEON_site.level_fg_forecast.fig.path <- paste0(dir,'NEON_site.level_fg_forecast.png')

#### ITS/pecan_gen: reference data product paths. ####
#these are not actually generated by scripts. Maybe this should move to a "small_data" higher level directory.
#All stuff that was downloaded by "hand".
dir <- paste0(ITS_pecan_gen_dir,'reference_data/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
              em_traits.path <- paste0(dir,'ecto_genus_traits_hobbie_Jan2018.csv')
             em_species.path <- paste0(dir,'myc_assignments.rds')
             poa_genera.path <- paste0(dir,'poaceae_genera_wikipedia.rds')
            #em_species.path <- paste0(dir,'mycorrhizal_SPCD_data.csv')
              em_genera.path <- paste0(dir,'tedersoo_2017_genera.csv')
 NEON_plantStatus_codes.path <- paste0(dir,'NEON_DP1.10098.plantStatus_decode.csv')
          ted_srr_codes.path <- paste0(dir,'tedersoo_srr_codes.csv')
            ted_map_raw.path <- paste0(dir,'tedersoo_2014_mapping_file.csv')
            ted_otu_raw.path <- paste0(dir,'tedersoo_2014_otu_file.txt')
     ted_sampling_dates.path <- paste0(dir,'tedersoo2014_dates.csv')
 
#### NEON pre-release: otu tables from L. Stanish. ####
#This is a one off, will be replaced by our own pipeline eventually.
#These data live outside the NEFI_data directory.
NEON_pre_release.dir <- '/fs/data3/caverill/NEON_pre-release_ITS-16S_otus/'
#outputs from pre-release processing.
neon_pre_release_otu.out_16S <- paste0(NEON_pre_release.dir,'16S_otu_clean.rds')
neon_pre_release_map.out_16S <- paste0(NEON_pre_release.dir,'16S_map_clean.rds')
neon_pre_release_tax.out_16S <- paste0(NEON_pre_release.dir,'16S_tax_clean.rds')
neon_pre_release_otu.out_ITS <- paste0(NEON_pre_release.dir,'ITS_otu_clean.rds')
neon_pre_release_map.out_ITS <- paste0(NEON_pre_release.dir,'ITS_map_clean.rds')
neon_pre_release_tax.out_ITS <- paste0(NEON_pre_release.dir,'ITS_tax_clean.rds')



########## ---------- 16S file paths ---------- ###########

#### --------- Bahram prior paths --------- ####

#### 16S/scc_gen: bahram processed sequences ####
dir <- paste0(scc_gen_16S_dir,'bahram_processed/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
bahram_dada2_SV_table.path <- paste0(dir,'bahram_dada2_SV_table.rds')
bahram_dada2_tax_table.path <- paste0(dir,'bahram_dada2_tax_table.rds')
bahram_dada2_track_table.path <- paste0(dir,'bahram_dada2_track_table.rds')

#### 16S/scc_gen: JAGS output ####
dir <- paste0(scc_gen_16S_dir, 'JAGS_output/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
bahram_16S.prior_12gen_JAGSfit <- paste0(dir,'bahram_16S.prior_12gen_JAGSfit.rds')
bahram_16S.prior_12gen_JAGSfit_moisture <- paste0(dir,'bahram_16S.prior_12gen_JAGSfit_moisture.rds')
bahram_16S.prior_12gen_JAGSfit_south_lat <- paste0(dir,'bahram_16S.prior_12gen_JAGSfit_south_lat.rds')
bahram_16S.prior_15phyla_JAGSfit <- paste0(dir,'bahram_16S.prior_15phyla_JAGSfit.rds')
  
#### 16S/scc_gen: raw data for prior metadata ####
dir <-  scc_gen_16S_dir
SRA.path <- paste0(dir, "SraRunTable.csv")
metadata_bahram_raw.path <- paste0(dir, "metadata_bahram_raw.csv")

#### 16S/pecan_gen: Bahram 2018 prior data. ####
dir <- pecan_gen_16S_dir
cmd <- paste0('mkdir -p ',dir)
system(cmd)
bahram_prior_gen.path <- paste0(dir,'bahram_prior_gen.rds') # 12 cosmopolitan genera
cosmo_output_16S.path <- paste0(dir, "cosmo_output_16S.rds") # 20 cosmopolitan genera - use this one, it's the same approach as for ITS
cosmo_output_16S_south_lat.path <- paste0(dir, "cosmo_output_16S_south_lat.rds") # 20 cosmopolitan genera incl southern lats
phyla_output_16S.path <- paste0(dir,"phyla_output_16S.rds") # 15 cosmopolitan phyla
phyla_output_16S_south_lat.path <- paste0(dir,"phyla_output_16S_south_lat.rds") # 15 cosmopolitan phyla incl southern lats
bahram_metadata.path <-  paste0(dir, "bahram_metadata.rds") # prior metadata - northern temperate latitudes
bahram_metadata_south_lat.path <-  paste0(dir, "bahram_metadata_south_lat.rds") # prior metadata - northern and souther temperate latitudes
bahram_prior_metadata.path <-  paste0(dir, "bahram_prior_metadata.rds")


#### ------------ NEON paths ------------- ####

### 16S/scc_gen:  NEON processed data ###
dir <- paste0(scc_gen_16S_dir,'NEON_processed/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
NEON_dada2_SV_table.path <- paste0(dir,'NEON_dada2_SV_table.rds')
NEON_dada2_tax_table.path <- paste0(dir,'NEON_dada2_tax_table.rds')
NEON_dada2_track_table.path <- paste0(dir,'NEON_dada2_track_table.rds')
NEON_plot.level_genera_obs_16S.path <- paste0(dir, 'NEON_plot.level_genera_obs_16S.rds')
NEON_site.level_genera_obs_16S.path <- paste0(dir, 'NEON_site.level_genera_obs_16S.rds')

#### 16S/pecan_gen ####

# paths to NEON cosmo forecasts and variable importance data
dir <- pecan_gen_16S_dir
NEON_cosmo_abundances_16S.path <- paste0(dir, "NEON_cosmo_abundances_16S.rds")
NEON_site_fcast_genera_16S.path <- paste0(dir, "NEON_site_fcast_cosmo_16S.rds")
hierarch_filled_16S.path <- paste0(dir, "hierarch_filled_16S.rds")
NEON_genera_variable_importance_data_16S.path <- paste0(dir, "NEON_genera_variable_importance_data_16S.rds")

#### 16S/pecan_gen: NEON data aggregation paths. ####
dir <- paste0(pecan_gen_16S_dir,'NEON_data_aggregation/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
    site_dates_16S.path <- paste0(dir,'site_dates_16S.rds')
    site_level_data_16S.path <- paste0(dir,'site_level_data_16S.rds')
    dp1.10108.00_output_16S.path <- paste0(dir,'dp1.10108.00_output_16S.rds')
    dp1.10109.00_output_16S.path <- paste0(dir,'dp1.10109.00_output_16S.rds')
    dp1.10086.00_output_16S.path <- paste0(dir,'dp1.10086.00_output_16S.rds')
    dp1.10078.00_output_16S.path <- paste0(dir,'dp1.10078.00_output_16S.rds')
    dp1.10098.00_output_16S.path <- paste0(dir,'dp1.10098.00_output_16S.rds')
    dp1.10098.00_plot.level_16S.path <- paste0(dir,'dp1.10098.00_plot.level_16S.rds')
    core_obs_16S.path <- paste0(dir, "core_obs_16S.rds")
    core_core_16S.path <- paste0(dir, "core_core_16S.rds")
    core_plot_16S.path <- paste0(dir, "core_plot_16S.rds")
    core_site_16S.path <- paste0(dir, "core_site_16S.rds")
    core_glob_16S.path <- paste0(dir, "core_glob_16S.rds")
    plot_plot_16S.path <- paste0(dir, "plot_plot_16S.rds")
    plot_site_16S.path <- paste0(dir, "plot_site_16S.rds")
    plot_glob_16S.path <- paste0(dir, "plot_glob_16S.rds")
    site_site_16S.path <- paste0(dir, "site_site_16S.rds")
    site_glob_16S.path <- paste0(dir, "site_glob_16S.rds")

#Paths to spatial x and y variables for modeling.
obs.table_16S.path <- paste0(dir, 'obs.table_16S.rds') 
core.table_16S.path <- paste0(dir,'core.table_16S.rds') 
plot.table_16S.path <- paste0(dir,'plot.table_16S.rds') 
site.table_16S.path <- paste0(dir,'site.table_16S.rds') 

# path to NEON 16S top 20 genera
NEON_gen_abundances.path <- paste0(dir, "NEON_gen_abundances.rds")
# path to NEON 16S all genera
NEON_all_gen_abundances.path <- paste0(dir, "NEON_all_gen_abundances.rds")
 
#### 16S/pecan_gen: paths to NEON covariates ####
dir <- paste0(pecan_gen_16S_dir,'NEON_covariates/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
NEON_site_covs_16S.path <- paste0(dir,'NEON_site_covs_16S.rds')
NEON_glob_covs_16S.path <- paste0(dir,'NEON_glob_covs_16S.rds')

#### 16S/pecan_gen: NEON forecast figures #### 
dir <- paste0(pecan_gen_16S_dir,'figures/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
NEON_cps.fcast_cosmo_16S.path <- paste0(dir, "NEON_cps.fcast_cosmo_16S.pdf")
NEON_core.fcast_cosmo_fig_16S.path <- paste0(dir, "NEON_core.fcast_cosmo_16S.pdf")
NEON_plot.fcast_cosmo_fig_16S.path <- paste0(dir, "NEON_plot.fcast_cosmo_fig_16S.pdf")
NEON_site.fcast_cosmo_fig_16S.path <- paste0(dir, "NEON_site.fcast_cosmo_fig_16S.pdf")

