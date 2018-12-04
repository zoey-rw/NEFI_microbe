#plotting phyla taxonomic group forecasts vs. validation data.
#clear environment, source paths, functions and packages.
rm(list=ls())
source('paths.r')

#set output path.----
output.path <- NEON_cps_phyla_forecast_figure.path
output.path <- 'test.png'

#load data.----
d <- readRDS(NEON_site_fcast_phyla.path)
core.truth <- readRDS(NEON_ITS_fastq_phyla.path)
core.truth <- core.truth$rel.abundances
rownames(core.truth) <- gsub('-GEN','',core.truth$geneticSampleID)
plot.truth <- readRDS(NEON_plot.level_phyla_obs_fastq.path)
site.truth <- readRDS(NEON_site.level_phyla_obs_fastq.path)

#DEFINE OUTLIER SITES- Russula follow the ECM fg, so sample outliers.----
out_sites <- c()
out_spp   <- c()

#setup output spec.----
png(filename=output.path,width=12,height=12,units='in',res=300)

#end plot.----
dev.off()
