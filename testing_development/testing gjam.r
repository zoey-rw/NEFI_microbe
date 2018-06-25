#Fitting a model with gjam.
library(gjam)
library(data.table)

#load data
d <- data.table(readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/ted_all_prior_data.rds'))

#fit a gjam, intercept only, to ecto abundance.

#problem- sum of ecto + spro is grater than 1, because some are listed as both. Gotta deal with this in data cleaning ahead of this.
#specify your y dependent variable matrix of species abundances
spp.y <- data.frame(d$Ectomycorrhizal,d$Saprotroph)
#get an other column, cause things need to sum to 1.
spp.y$S3other <- 1 - (rowSums(spp.y))

#make x data.
spp.x <- data.frame(1, d$pH)
colnames(spp.x) <- c('intercept','x1')

#get a vector declaring the species abundance data types. Here we want fractional composition, FC.
#for some reason this work swith the sim genenrated stuff, but not my list...
types <- c('FC','FC','FC')
sim <- gjamSimData(S = 3, typeNames = 'FC')

#get a formula.
g.formula <- as.formula(~ x1)

#make model
m1 <- list(ng = 50, burnin = 5, typeNames = sim$typeNames)

#works. Now how do we extract species interactions from this?
test <- gjam(g.formula, spp.x, spp.y, modelList = m1)




