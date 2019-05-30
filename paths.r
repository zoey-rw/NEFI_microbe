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

#Raw NEON ITS provided in fastq format, both reads!
NEON_ITS_fastq.dir <- paste0(big_data_dir,'NEON_raw_ITS_fastq/ITS_run150922/per_sample_demux/')
#broken out per sequencing run, fwd/rev reads.
NEON_ITS_run150225_r1_fastq.dir <- paste0(big_data_dir,'NEON_raw_ITS_fastq/ITS_run150225/r1_per_sample_demux/')
NEON_ITS_run150225_r2_fastq.dir <- paste0(big_data_dir,'NEON_raw_ITS_fastq/ITS_run150225/r2_per_sample_demux/')
NEON_ITS_run150922_r1_fastq.dir <- paste0(big_data_dir,'NEON_raw_ITS_fastq/ITS_run150922/r1_per_sample_demux/')
NEON_ITS_run150922_r2_fastq.dir <- paste0(big_data_dir,'NEON_raw_ITS_fastq/ITS_run150922/r2_per_sample_demux/')
#dada2 output- make these.
NEON_ITS_run150225_dada2_out.dir <- paste0(big_data_dir,'NEON_raw_ITS_fastq/ITS_run150225/dada2_output/')
NEON_ITS_run150922_dada2_out.dir <- paste0(big_data_dir,'NEON_raw_ITS_fastq/ITS_run150922/dada2_output/')
system(paste0('mkdir -p ',NEON_ITS_run150225_dada2_out.dir))
system(paste0('mkdir -p ',NEON_ITS_run150922_dada2_out.dir))
     
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
ted_ITS.prior_phyla_JAGSfit <- paste0(dir,'ted_ITS.prior_phyla_JAGSfit.rds')
ITS.prior_linear_fg_cov.selection_JAGS <- paste0(dir,'ITS.prior_linear_fg_cov.selection_JAGS.rds')
ted_ITS_prior_phylo.group_JAGSfits <- paste0(dir,'ted_ITS.pror_phylo.groups_JAGSfits.rds')

#dmulti-ddirch output, fit with fastq-derived data.
ted_ITS.prior_dmulti.ddirch_fg_JAGSfit <- paste0(dir,'ted_ITS.prior_dmulti-ddirch_fg_JAGSfit.rds')
ted_ITS.prior_dmulti.ddirch_cosmo_JAGSfit <- paste0(dir,'ted_ITS.prior_dmulti-ddirch_cosmo_JAGSfit.rds')
ted_ITS.prior_dmulti.ddirch_yeast_JAGSfit <- paste0(dir,'ted_ITS.prior_dmulti-ddirch_yeast_JAGSfit.rds')

#### ITS/scc_gen: tedersoo 2014 SV and taxonomy paths. ####
ted_2014_SV.table.path <- paste0(ITS_scc_gen_dir,'ted_2014_SV.table.rds')
     ted_2014_tax.path <- paste0(ITS_scc_gen_dir,'ted_2014_tax.rds') 
     
#### ITS/scc_gen: NEON to forecast SV and taxonomy paths. ####
NEON_SV.table.path <- paste0(ITS_scc_gen_dir,'NEON_SV.table.rds')
     NEON_tax.path <- paste0(ITS_scc_gen_dir,'NEON_tax.rds') 
     NEON_fun.path <- paste0(ITS_scc_gen_dir,'NEON_fun.rds')
#table derived from .fastq raw files.
     NEON_ITS_fastq_SV.table.path <- paste0(ITS_scc_gen_dir,'NEON_fastq_SV.table.rds')
     NEON_ITS_fastq_tax.path <- paste0(ITS_scc_gen_dir,'NEON_fastq_tax.rds')
     NEON_ITS_fastq_fun.path <- paste0(ITS_scc_gen_dir,'NEON_fastq_fun.rds')
          
#### ITS/pecan_gen: Figures.----
dir <- paste0(ITS_pecan_gen_dir,'figures/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)

#ddirch forecast-validation.
NEON_cps_fg_forecast_figure.path <- paste0(dir,'NEON_CPS_fg_forecast_figure.png')
NEON_cps_cosmo_forecast_figure.path <- paste0(dir,'NEON_CPS_cosmo_forecast_figure.png')
NEON_cps_phyla_forecast_figure.path <- paste0(dir,'NEON_CPS_phyla_forecast_figure.png')
NEON_ddirch_var.decomp_fg_figure.path <- paste0(dir,'NEON_ddirch_var.decomp_fg.png')
NEON_ddirch_var.decomp_cosmo_figure.path <- paste0(dir,'NEON_ddirch_var.decomp_cosmo.png')
NEON_ddirch_var.importance_fg_figure.path <- paste0(dir,'NEON_ddirch_var.importance_fg.png')
NEON_ddirch_var.importance_cosmo_figure.path <- paste0(dir,'NEON_ddirch_var.importance_fg.png')

#dmulti-ddirch forecast-validation.
NEON_dmulti.ddirch_fg_forecast_validation.path <- paste0(dir,'NEON_dmulti.ddirch_fg_forecast_validation.png')
NEON_dmulti.ddirch_cosmo_forecast_validation.path <- paste0(dir,'NEON_dmulti.ddirch_cosmo_forecast_validation.png')
NEON_dmulti.ddirch_yeast_forecast_validation.path <- paste0(dir,'NEON_dmulti.ddirch_yeast_forecast_validation.png')

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
            tedersoo_ITS_clean_map.path <- paste0(dir,'tedersoo_ITS_clean_map.rds')
              tedersoo_ITS_fg_list.path <- paste0(dir,'tedersoo_ITS_fg_list.rds')
           tedersoo_ITS_hydro_list.path <- paste0(dir,'tedersoo_ITS_hydro_list.rds')
           tedersoo_ITS_yeast_list.path <- paste0(dir,'tedersoo_ITS_yeast_list.rds')
    tedersoo_ITS_cosmo_genera_list.path <- paste0(dir,'tedersoo_ITS_cosmo_genera_list.rds')
           tedersoo_ITS_phyla_list.path <- paste0(dir,'tedersoo_ITS_phyla_list.rds')
tedersoo_ITS_common_phylo_groups_list.path <- paste0(dir,'tedersoo_ITS_common_phylo_groups_list.rds')

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

#cleanup of NEON fun (which also has taxonomy) and SV tables.
NEON_SV.table_clean.path <- paste0(dir,'NEON_SV.table_clean.rds')
     NEON_fun_clean.path <- paste0(dir,'NEON_fun_clean.rds')
     
#cleanup of NEON fun (which also has taxonomy) and SV tables - DERIVED FROM RAW FASTQ.
NEON_ITS_fastq_SV.table_clean.path <- paste0(dir,'NEON_ITS_fastq_SV.table_clean.rds')
     NEON_ITS_fastq_fun_clean.path <- paste0(dir,'NEON_ITS_fastq_fun_clean.rds')

#### ITS/pecan_gen: Forecast covariate paths.----
#Aggregated NEON site level covariates and global level uncertainty for predictors
dir <- paste0(ITS_pecan_gen_dir,'NEON_covariates/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
        NEON_site_covs.path <- paste0(dir,'NEON_site_covs.rds')
        NEON_glob_covs.path <- paste0(dir,'NEON_glob_covs.rds')
          NEON_taxa_fg.path <- paste0(dir,'NEON_taxa_fg.rds')
     NEON_cosmo_genera.path <- paste0(dir,'NEON_cosmo_genera.rds')
NEON_site.level_fg_obs.path <- paste0(dir,'NEON_site.level_fg_obs.rds')
NEON_plot.level_fg_obs.path <- paste0(dir,'NEON_plot.level_fg_obs.rds')
NEON_site.level_genera_obs.path <- paste0(dir,'NEON_site.level_genera_obs.rds')
NEON_plot.level_genera_obs.path <- paste0(dir,'NEON_plot.level_genera_obs.rds')

#Aggregated fg and cosmo_genera for fastq-dervied data.
     NEON_ITS_fastq_taxa_fg.path <- paste0(dir,'NEON_ITS_fastq_taxa_fg.rds')
NEON_ITS_fastq_yeast_taxa.path <- paste0(dir,'NEON_ITS_fastq_yeast_taxa.rds')
NEON_ITS_fastq_cosmo_genera.path <- paste0(dir,'NEON_ITS_fastq_cosmo_genera.rds')
NEON_ITS_fastq_all_cosmo_phylo_groups.path <- paste0(dir,'NEON_ITS_fastq_all_cosmo_phylo_groups.rds')
NEON_ITS_fastq_phyla.path <- paste0(dir,'NEON_ITS_fastq_phyla.rds')
NEON_plot.level_fg_obs_fastq.path <- paste0(dir,'NEON_plot.level_fg_obs_fastq.rds')
NEON_site.level_fg_obs_fastq.path <- paste0(dir,'NEON_site.level_fg_obs_fastq.rds')
NEON_plot.level_genera_obs_fastq.path <- paste0(dir,'NEON_plot.level_genera_obs_fastq.rds')
NEON_site.level_genera_obs_fastq.path <- paste0(dir,'NEON_site.level_genera_obs_fastq.rds')
NEON_plot.level_yeast_obs_fastq.path <- paste0(dir,'NEON_plot.level_yeast_obs_fastq.rds')
NEON_site.level_yeast_obs_fastq.path <- paste0(dir,'NEON_site.level_yeast_obs_fastq.rds')
NEON_plot.level_phyla_obs_fastq.path <- paste0(dir,'NEON_plot.level_phyla_obs_fastq.rds')
NEON_site.level_phyla_obs_fastq.path <- paste0(dir,'NEON_site.level_phyla_obs_fastq.rds')


#### ITS/pecan_gen: Forecasts, variance decomposition and variable importance output paths.----
dir <- paste0(ITS_pecan_gen_dir,'NEON_forecasts/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)

#ddirch forecasts and variable importance output.
NEON_site_fcast_fg.path <- paste0(dir,'NEON_fcast_site.level_fg.rds')
NEON_site_fcast_genera.path <- paste0(dir,'NEON_fcast_site.level_genera.rds')
NEON_site_fcast_phyla.path <- paste0(dir,'NEON_fcast_phyla.rds')
NEON_site_fcast_all_phylo_levels.path <- paste0(dir,'NEON_fcast_all_phylo_levels.rds')
NEON_fg_variable_importance_data.path <- paste0(dir,'NEON_fg_variable_importance_data.rds')
NEON_genera_variable_importance_data.path <- paste0(dir,'NEON_genera_variable_importance_data.rds')
NEON_phyla_variable_importance_data.path <- paste0(dir,'NEON_phyla_variable_importance_data.rds')
NEON_ddirch_var.decomp_fg.path <- paste0(dir,'NEON_ddirch_var.decomp_fg.rds')
NEON_ddirch_var.decomp_cosmo.path <- paste0(dir,'NEON_ddirch_var.decomp_cosmo.rds')
NEON_ddirch_var.decomp_phyla.path <- paste0(dir,'NEON_ddirch_var.decomp_phyla.rds')

#dmulti-ddirch forecasts, variance decomposition and variable importance output.
NEON_dmulti.ddirch_fcast_fg.path <- paste0(dir,'NEON_dmulti.ddirch_fcast_fg.rds')
NEON_dmulti.ddirch_fcast_cosmo.path <- paste0(dir,'NEON_dmulti.ddirch_fcast_cosmo.rds')
NEON_dmulti.ddirch_fcast_yeast.path <- paste0(dir,'NEON_dmulti.ddirch_fcast_yeast.rds')
NEON_dmulti.ddirch_var.importance_fg.path <- paste0(dir,'NEON_dmulti.ddirch_var.importance_fg.rds')
NEON_dmulti.ddirch_var.importance_cosmo.path <- paste0(dir,'NEON_dmulti.ddirch_var.importance_cosmo.rds')
NEON_dmulti.ddirch_var.importance_yeast.path <- paste0(dir,'NEON_dmulti.ddirch_var.importance_yeast.rds')
NEON_dmulti.ddirch_var.decomposition_fg.path <- paste0(dir,'NEON_dmulti.ddirch_var.decomposition_fg.rds')
NEON_dmulti.ddirch_var.decomposition_cosmo.path <- paste0(dir,'NEON_dmulti.ddirch_var.decomposition_cosmo.rds')
NEON_dmulti.ddirch_var.decomposition_yeast.path <- paste0(dir,'NEON_dmulti.ddirch_var.decomposition_yeast.rds')

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
          NEON_ITS_link_file <- paste0(dir,'NEON_rawFilesList.csv')
            field_sites.path <- paste0(dir,'field-sites.csv')
     
 
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

#####ITS/scc_gen: seq depth analysis paths.----
dir <- paste0(ITS_scc_gen_dir,'seq_depth_analysis/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
ted_fg_ddirch_fit_seqdepth.path <- paste0(dir,'ted_fg_ddirch_fit.rds')
ted_cosmo_ddirch_fit_seqdepth.path <- paste0(dir,'ted_fg_ddirch_fit.rds')

####ITS/pecan_gen: seq depth analysis paths.----
dir <- paste0(ITS_pecan_gen_dir,'seq_depth_analysis/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
ted_fg_seq.depth_ddirch_foreacsts.path <- paste0(dir,'ted_fg_seq.depth_ddirch_foreacsts.rds')
ted_cosmo_seq.depth_ddirch_foreacsts.path <- paste0(dir,'ted_cosmo_seq.depth_ddirch_foreacsts.rds')

########## ---------- 16S file paths ---------- ###########

#### --------- Bahram prior paths --------- ####

#### 16S/scc_gen: bahram processed sequences ####
dir <- paste0(scc_gen_16S_dir,'bahram_processed/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
bahram_dada2_SV_table.path <- paste0(dir,'bahram_dada2_SV_table.rds')
bahram_dada2_SV_table_rare_not_subset.path <- paste0(dir,'bahram_dada2_SV_table_rare_not_subset.rds')
bahram_dada2_SV_table_rare_all.samples.path <- paste0(dir, 'bahram_dada2_SV_table_rare_all.samples.rds')
bahram_dada2_SV_table_rare.path <- paste0(dir,'bahram_dada2_SV_table_rare.rds')
bahram_dada2_tax_table.path <- paste0(dir,'bahram_dada2_tax_table.rds')
bahram_dada2_track_table.path <- paste0(dir,'bahram_dada2_track_table.rds')

#### 16S/scc_gen: JAGS output ####
dir <- paste0(scc_gen_16S_dir, 'JAGS_output/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
bahram_16S.prior_15phyla_JAGSfit <- paste0(dir,'bahram_16S.prior_15phyla_JAGSfit.rds')
bahram_16S.prior_cop_olig_JAGSfit <-   paste0(dir,'bahram_16S.prior_cop_olig_JAGSfit.rds')
bahram_16S.prior_cop_olig_all_nutr_JAGSfit <- paste0(dir,'bahram_16S.prior_cop_olig_all_nutr_JAGSfit.rds')
bahram_16S.prior_cop_olig_all_nutr_no_moist_JAGSfit <- paste0(dir,'bahram_16S.prior_cop_olig_all_nutr_no_moist_JAGSfit.rds')
bahram_16S.prior_cop_olig_all_var_JAGSfit <- paste0(dir,'bahram_16S.prior_cop_olig_all_JAGSfit.rds')
bahram_16S_prior_phylo.group_JAGSfits <- paste0(dir,'bahram_16S.prior_phylo_JAGSfit.rds')
bahram_16S_prior_N_cycle_JAGSfits <- paste0(dir,'bahram_16S.prior_N_cycle_JAGSfit.rds')
bahram_16S_prior_C_cycle_JAGSfits <- paste0(dir,'bahram_16S.prior_C_cycle_JAGSfit.rds')
bahram_16S_prior_fg_cov.selection_JAGS <- paste0(dir,'bahram_16S_prior_fg_cov.selection_JAGS')
prior_16S_all.fg.groups_JAGSfits.path <- paste0(dir,'prior_16S_all.fg.groups_JAGSfits')
bahram_16S_prior_dmulti.ddirch_all.group_JAGSfits <- paste0(dir,'bahram_16S_prior_dmulti.ddirch_all.group_JAGSfits')

#### 16S/scc_gen: raw sequence and metadata for prior metadata ####
dir <-  scc_gen_16S_dir
SRA.path <- paste0(dir, "SraRunTable.csv")
metadata_bahram_raw.path <- paste0(dir, "metadata_bahram_raw.csv")

#### 16S/pecan_gen: Bahram 2018 prior abundance data. ####
dir <- paste0(pecan_gen_16S_dir, 'prior_data/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
# group abundances
bahram_16S_common_phylo_fg_abun.path <- paste0(dir,"bahram_16S_common_phylo_fg_abun.rds") # all cosmopolitan taxa and functional groups
phyla_output_16S.path <- paste0(dir,"phyla_output_16S.rds") # 15 cosmopolitan phyla
prior_cop_olig_abundances.path <- paste0(dir,"prior_cop_olig_abundances.rds") # copiotrophic/oligotrophic relative abundances - Bahram
prior_N_cyclers_abundances.path <- paste0(dir, "prior_N_cyclers_abundances.rds")
prior_C_cyclers_abundances.path <- paste0(dir, "prior_C_cyclers_abundances.rds")
prior_fg_abundances_16S.path <- paste0(dir, "prior_fg_abundances_16S.rds")
bahram_tax_fg_16S.path <- paste0(dir, "bahram_tax_fg_16S.rds")
bahram_16S_phylo_lineage_abun_16S.path <- paste0(dir, "bahram_16S_phylo_lineage_abun_16S.rds")

# metadata
bahram_metadata.path <-  paste0(dir, "bahram_metadata.rds") # processed prior metadata


#### ------------ NEON paths ------------- ####

### 16S/scc_gen:  NEON processed abundance data ###
dir <- paste0(scc_gen_16S_dir,'NEON_processed/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
NEON_dada2_SV_table.path <- paste0(dir,'NEON_dada2_SV_table.rds')
NEON_dada2_SV_table_rare.path <- paste0(dir,'NEON_dada2_SV_table_rare.rds')
NEON_dada2_tax_table.path <- paste0(dir,'NEON_dada2_tax_table.rds')
NEON_dada2_track_table.path <- paste0(dir,'NEON_dada2_track_table.rds')
NEON_plot.level_genera_obs_16S.path <- paste0(dir, 'NEON_plot.level_genera_obs_16S.rds')
NEON_site.level_genera_obs_16S.path <- paste0(dir, 'NEON_site.level_genera_obs_16S.rds')
NEON_plot.level_phyla_obs_16S.path <- paste0(dir,'NEON_plot.level_phyla_obs_16S.rds')
NEON_site.level_phyla_obs_16S.path <- paste0(dir,'NEON_site.level_phyla_obs_16S.rds')
NEON_all.fg_plot.site_obs_16S.path <- paste0(dir,'NEON_all.fg_plot.site_obs_16S.rds')
NEON_all.phylo.levels_plot.site_obs_16S.path <- paste0(dir, 'NEON_all.phylo.levels_plot.site_obs_16S.rds')

#### 16S/pecan_gen ####

# paths to NEON cosmo forecasts and variable importance data
dir <- paste0(pecan_gen_16S_dir, 'NEON_forecast_data/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
NEON_cosmo_abundances_16S.path <- paste0(dir, "NEON_cosmo_abundances_16S.rds")
NEON_site_fcast_genera_16S.path <- paste0(dir, "NEON_site_fcast_cosmo_16S.rds")
NEON_genera_variable_importance_data_16S.path <- paste0(dir, "NEON_genera_variable_importance_data_16S.rds")
NEON_cps_fcast_all_phylo_16S.path <- paste0(dir, "NEON_cps_fcast_all_phylo_16S.rds")
NEON_cps_fcast_phyla_16S.path <- paste0(dir, "NEON_cps_fcast_phyla_16S.rds")
NEON_cps_fcast_N.cycler_16S.path  <- paste0(dir, "NEON_cps_fcast_N.cycler_16S.rds")
NEON_cps_fcast_fg_16S.path <- paste0(dir, "NEON_cps_fcast_fg_16S.rds")
NEON_cps_fcast_fg_16S_no.nutr.path <- paste0(dir, "NEON_cps_fcast_fg_16S_no.nutr.rds")
NEON_ddirch_var.decomp_fg_16S.path <- paste0(dir, "NEON_ddirch_var.decomp_fg_16S.rds")
NEON_fg_variable_importance_data_16S.path <- paste0(dir, "NEON_fg_variable_importance_data_16S.rds")
NEON_ddirch_var.decomp_fg_16S.path <- paste0(dir, "NEON_ddirch_var.decomp_fg_16S.rds")
NEON_ddirch_var.decomp_phylo.group_16S.path <- paste0(dir, "NEON_ddirch_var.decomp_phylo.group_16S.rds")

#### 16S/pecan_gen: NEON covariate data aggregation paths. ####
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
    dp1.10008.00_output_16S.path <- paste0(dir,'dp1.10008.00_output_16S.rds')
    dp1.10047.00_output_16S.path <- paste0(dir,'dp1.10047.00_output_16S.rds')
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
    hierarch_filled_16S.path <- paste0(dir, "hierarch_filled_16S.rds")
    missing_data_removed_16S.path <- paste0(dir, "missing_data_removed_16S.rds")
    
#Paths to spatial x and y variables for modeling.
obs.table_16S.path <- paste0(dir, 'obs.table_16S.rds') 
core.table_16S.path <- paste0(dir,'core.table_16S.rds') 
plot.table_16S.path <- paste0(dir,'plot.table_16S.rds') 
site.table_16S.path <- paste0(dir,'site.table_16S.rds') 

# path to NEON group abundances
NEON_16S_phylo_groups_abundances.path <- paste0(dir, "NEON_16S_phylo_groups_abundances.rds")
NEON_fg_abundances_16S.path <- paste0(dir, "NEON_fg_abundances_16S.rds")


# NEON covariate data (16S/ITS)
dir <- paste0(pecan_gen_16S_dir,'NEON_covariates_16S_ITS/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)

# raw data from NEON, WorldClim, etc
dir <- paste0(pecan_gen_16S_dir,'NEON_covariates_16S_ITS/NEON_covariates_raw/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
site_dates.path <- paste0(dir, "sites_dates.rds")
dp1.10108.001_output.path <- paste0(dir, "dp1.10108.001_output.rds")
dp1.10086.001_output.path <- paste0(dir, "dp1.10086.001_output.rds")
dp1.10109.001_output.path <- paste0(dir, "dp1.10109.001_output.rds")
dp1.10078.001_output.path <- paste0(dir, "dp1.10078.001_output.rds")
dp1.10098.001_output.path <- paste0(dir, "dp1.10098.001_output.rds")
dp1.10098.001_plot.level_output.path <- paste0(dir, "dp1.10098.001_plot.level_output.rds")
site_level_covariates.path <- paste0(dir, "site_level_covariates.rds")

# aggregated data for forecasts
dir <- paste0(pecan_gen_16S_dir,'NEON_covariates_16S_ITS/NEON_covariates_aggregated/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
core_obs_data.path <- paste0(dir, "core_obs_data.rds")
core_core_data.path <- paste0(dir, "core_core_data.rds")
core_plot_data.path <- paste0(dir, "core_plot_data.rds")
core_site_data.path <- paste0(dir, "core_site_data.rds")
core_glob_data.path <- paste0(dir, "core_glob_data.rds")
plot_plot_data.path <- paste0(dir, "plot_plot_data.rds")
plot_site_data.path <- paste0(dir, "plot_site_data.rds")
plot_glob_data.path <- paste0(dir, "plot_glob_data.rds")
site_site_data.path <- paste0(dir, "site_site_data.rds")
site_glob_data.path <- paste0(dir, "site_glob_data.rds")
hierarch_filled_data.path <- paste0(dir, "hierarch_filled_data.rds")













#to delete:
NEON_gen_abundances.path <- paste0(dir, "NEON_gen_abundances.rds")
NEON_all_gen_abundances.path <- paste0(dir, "NEON_all_gen_abundances.rds")
NEON_phyla_abundances_16S.path  <- paste0(dir, "NEON_phyla_abundances.rds")
NEON_N_cyclers_abundances.path <- paste0(dir, "NEON_N_cyclers_abundances.rds")
NEON_C_cyclers_abundances.path <- paste0(dir, "NEON_C_cyclers_abundances.rds")
NEON_cop_olig_abundances.path <- paste0(dir,"NEON_cop_olig_abundances.rds") # copiotrophic/oligotrophic relative abundances - NEON


#### 16S/pecan_gen: paths to NEON covariates ####
dir <- paste0(pecan_gen_16S_dir,'NEON_covariates/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
NEON_site_covs_16S.path <- paste0(dir,'NEON_site_covs_16S.rds')
NEON_glob_covs_16S.path <- paste0(dir,'NEON_glob_covs_16S.rds')

#### 16S/pecan_gen: figures
dir <- paste0(pecan_gen_16S_dir,'figures/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
# NEON forecast figures
NEON_cps.fcast_cosmo_16S.path <- paste0(dir, "NEON_cps.fcast_cosmo_16S.pdf")
NEON_core.fcast_cosmo_fig_16S.path <- paste0(dir, "NEON_core.fcast_cosmo_16S.pdf")
NEON_plot.fcast_cosmo_fig_16S.path <- paste0(dir, "NEON_plot.fcast_cosmo_fig_16S.pdf")
NEON_site.fcast_cosmo_fig_16S.path <- paste0(dir, "NEON_site.fcast_cosmo_fig_16S.pdf")
NEON_cps.fcast_phyla_16S.path <- paste0(dir, "NEON_cps.fcast_phyla_16S.pdf")
NEON_cps.fcast_fg_fig_16S.path <- paste0(dir, "NEON_cps.fcast_fg_fig_16S.pdf")
# prior calibration figures 
prior_16S_r2_distribution.density_figure.path <- paste0(dir, "prior_16S_r2_distribution.density_fig.pdf")

#### 16S/pecan_gen: reference_data
dir <- paste0(pecan_gen_16S_dir,'reference_data/')
cmd <- paste0('mkdir -p ',dir)
system(cmd)
pH_conversion.path <- paste0(dir, "pH_conversion_data.rds")
  