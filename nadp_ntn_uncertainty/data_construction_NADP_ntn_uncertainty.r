#quantifying NAPP_NTN wet N deposition uncertainity.
#Need to get annual NH4, NO3 data extracted.
#sum together both, do a regression. See if variance scales with elevation.
#Getting NADP_NTN wet deposition uncertainty
ntn_data <- read.csv('/fs/data3/caverill/NEFI_microbial/NADP_ntn_uncertainty/NTN-All-cy.csv')
ntn_site <- read.csv('/fs/data3/caverill/NEFI_microbial/NADP_ntn_uncertainty/NTN_site_data.csv')

#Get sd of precited vs observed N deposition.
##NH4 and NO3 are in mg / L (ppm)
##ppt is the total centimeters of rain
##svol is the total mL of precipitation collected.
##svol (cm3) / area (cm2) = centimeters (cm1)

d <- merge(ntn_data, ntn_site[,c('Site.ID','Latitude','Longitude')], by.x = 'siteID', by.y = 'Site.ID')