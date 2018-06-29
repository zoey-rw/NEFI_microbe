#assigning mycorrhizal status and plot level relative EM abundance, forest cover, and conifer presence.
#clear environment, source paths.
rm(list=ls())
source('paths.r')

#load tree data
dat <- readRDS(dp1.10098.00_output.path)
#load mycorrhizal data
myc.spp <- readRDS(em_species.path)
myc.gen <- read.csv(em_genera.path)
myc.spp$genus_spp <- paste(myc.spp$GENUS,myc.spp$SPECIES,sep = ' ')
poa.gen <- readRDS(poa_genera.path)

#known US AM genera.
#AM genera based on searching super database, if there are at least 5 recrods and all are AM, then it gets assigned AM at the genus level.
#if a ton of AM and one AM_ECM, still count as AM.
am_genera <- c('Thuja','Fraxinus','Nyssa','Celtis','Cornus','Diospyros','Ilex','Lonicera','Magnolia','Viburnum',poa.gen)
erm_genera <- c('Rhododendron','Vaccinium')

#Everything is ID'd at least to genus. Species has lots of qualifiers. Lets clear these up.
dat$species <- sub("^(\\S*\\s+\\S+).*", "\\1", dat$scientificName)
dat$genus   <- sub(" .*", ""                 , dat$scientificName)

#assign mycorrhizal status
dat <- merge(dat,myc.spp[,c('Species','MYCO_ASSO')], by.x = 'species', by.y = 'Species', all.x=T)
dat[dat$genus %in% myc.gen$genus,]$MYCO_ASSO <- 'ECM'
dat[dat$genus %in% am_genera,]$MYCO_ASSO <- 'AM'
dat[dat$genus %in% erm_genera,]$MYCO_ASSO <- 'ERM'

#subset to trees that have a stem diameter measurement.
dat <- dat[grep('tree',dat$growthForm),]
dat <- dat[!(is.na(dat$stemDiameter)),]

#how much of the basal area is assigned? 96%. We good.
dat$basal <- pi * (dat$stemDiameter/2)^2
unassigned <- dat[is.na(dat$MYCO_ASSO),]
(1 - (sum(unassigned$basal) / sum(dat$basal)))*100

#We now need to account for dead trees, insect damaged trees. FIgure out plantStatus codes.
#Need to pick a time point (if there are even multiple time points) for each site.
unique(dat$plantStatus)
dat[dat$plantStatus == '4',]

#Aggregate by plot.
#exclude 50% damaged/dead?