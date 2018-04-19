#this is the neon product lookup table. 
#Each neon product query returns a column of file paths.
#only one of those filepaths is actually the data you want for a given site by year_month combination.
#The position of that path within the column varies even within product.
#the unique portion of that name is unique to each product.
#this table will get you where you need to go. Would love to not have to use it.

#specify product by grep_string combons
lookup <- data.frame(rbind(
      c('DP1.10086.001','soilCore','soil physical properties'),
      c('DP1.10078.001','soilChemistry','soil chemical properties'),
      c('DP1.00097.001','mgc_perbiogeosample', 'soil chemical properties (Megapit)'),
      c('DP1.10033.001','ltr_massdata','litterfall and fine woody debris sampling'),
      c('DP1.10066.001','','root sampling (megapit)'),
      c('DP1.10080.001','','soil inorganic nitrogen pools and transformations'),
      c('DP1.10104.001','','soil microbe biomass'),
      c('DP1.10026.001','','plant foliar physical and chemical properties'),
      c('DP1.00098.001','','coarse down wood log survey'),
      c('DP1.10014.001','','Coarse downed wood bulk density sampling')
      ))

#name the columns
colnames(lookup) <- c('neon_product','grep_string','name')

#save the file
write.csv(lookup,'/home/caverill/NEFI_microbe/neon_product_grep_lookup.csv')

#check out a neon product to figure out grep string.
nneo::nneo_product('DP1.10080.001')
test <- nneo::nneo_data(product_code = 'DP1.10080.001',site = 'SCBI', year_month = '2017-03')
test$data$files$name
prod <- read.delim(test$data$files$url
                   [grep('NEON.D02.SCBI.DP1.10080.001.variables.20171221T212823Z.csv',test$data$files$name)], sep = ',')
prod <- read.delim(test$data$files$url[16], sep = ',')
