#assigning mycorrhizal status and plot level relative EM abundance, forest cover, and conifer presence.
#clear environment, source paths.
rm(list=ls())
source('paths.r')

#load tree data
dat <- readRDS(dp1.10098.00_output.path)

#load lookup table for legacy plantStatus codes provided by Katie Jones at NEON.
p.codes <- read.csv(NEON_plantStatus_codes.path)

#load mycorrhizal data
myc.spp <- readRDS(em_species.path)
myc.gen <- read.csv(em_genera.path)
myc.spp$genus_spp <- paste(myc.spp$GENUS,myc.spp$SPECIES,sep = ' ')
poa.gen <- readRDS(poa_genera.path)

#known US AM genera.
#AM genera based on searching Colin's super mycorrhizal database. 
#If there are at least 5 recrods and all are AM, then it gets assigned AM at the genus level.
#if a ton of AM and one AM_ECM, still count as AM.
am_genera <- c('Thuja','Fraxinus','Nyssa','Celtis','Cornus','Diospyros','Ilex','Lonicera','Magnolia','Viburnum',poa.gen)
erm_genera <- c('Rhododendron','Vaccinium')

#Everything is ID'd at least to genus. Species has lots of qualifiers. Lets clear these up.
dat$species <- sub("^(\\S*\\s+\\S+).*", "\\1", dat$scientificName)
dat$genus   <- sub(" .*", ""                 , dat$scientificName)

#assign mycorrhizal status
dat <- merge(dat,myc.spp[,c('Species','MYCO_ASSO')], by.x = 'species', by.y = 'Species', all.x=T)
dat[dat$genus %in% myc.gen$genus,]$MYCO_ASSO <- 'ECM'
dat[dat$genus %in% am_genera,    ]$MYCO_ASSO <- 'AM'
dat[dat$genus %in% erm_genera,   ]$MYCO_ASSO <- 'ERM'

#subset to trees that have a stem diameter measurement.
dat <- dat[grep('tree',dat$growthForm),]
dat <- dat[!(is.na(dat$stemDiameter)),]

#how much of the basal area is assigned a mycorrhizal association? 96%. We good.
dat$basal <- pi * (dat$stemDiameter/2)^2
(1 - (sum(dat[is.na(dat$MYCO_ASSO),]$basal) / sum(dat$basal)))*100

#We now need to account for dead trees, insect damaged trees.
#deal with legacy codes in data using key.
for(i in 1:nrow(dat)){
  if(dat$plantStatus[i] %in% p.codes$lovElementName){
     dat$plantStatus[i] <- as.character(p.codes[p.codes$lovElementName == dat$plantStatus[i],]$lovElementCode)
  }
}

#Assign live, live_ECM and dead basal area.
dat$basal_live <- ifelse(grepl('Live', dat$plantStatus) == T, dat$basal, 0)
dat$basal_dead <- ifelse(grepl('Dead', dat$plantStatus) == T | grepl('dead', dat$plantStatus) == T, dat$basal, 0)
dat$basal_ECM  <- ifelse(dat$MYCO_ASSO == 'ECM', dat$basal_live, 0)

#aggregate to plot scale
plot.level            <- aggregate(basal_live ~ plotID, data = dat, FUN = sum, na.rm=T, na.action = na.pass)
plot.level$basal_ECM  <- aggregate(basal_ECM  ~ plotID, data = dat, FUN = sum, na.rm=T, na.action = na.pass)[,2]
plot.level$basal_dead <- aggregate(basal_dead ~ plotID, data = dat, FUN = sum, na.rm=T, na.action = na.pass)[,2]
#three disney plots have zero basal area. drop these.
plot.level <- plot.level[!(plot.level$basal_live == 0),]
plot.level$relEM <- plot.level$basal_ECM / plot.level$basal_live
plot.level$live_fraction <- plot.level$basal_live / (plot.level$basal_live + plot.level$basal_dead)
plot.level$siteID <- substring(plot.level$plotID,1,4)

#save plot-level EM relative abundance.
saveRDS(plot.level, dp1.10098.00_plot.level.path)
