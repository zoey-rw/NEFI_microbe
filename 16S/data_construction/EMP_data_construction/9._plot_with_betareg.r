#looping models of bacterial relative abundance.
#clear environment, load packages.
rm(list=ls())
source('paths.r')
library(betareg)

bahram_tax <- c("Proteobacteria", "Verrucomicrobia", "Actinobacteria", "Planctomycetes",
                "Chloroflexi", "Acidobacteria", "Firmicutes", "Bacteroidetes",
                "Gemmatimonadetes", "Armatimonadetes", "Alphaproteobacteria",
                "[Spartobacteria]", "Actinobacteria", "Planctomycetia", "Ellin6529",
                "Thermoleophilia", "Acidobacteria-6", "Bacilli", "Gammaproteobacteria",
                "Betaproteobacteria", "[Chloracidobacteria]", "Deltaproteobacteria",
                "Solibacteres", "[Saprospirae]", "[Pedosphaerae]", "Acidimicrobiia",
                "Opitutae", "Verrucomicrobiae", "Sva0725", "TK10", "Sphingobacteriia",
                "Gemmatimonadetes", "Phycisphaerae", "Rhizobiales", "[Chthoniobacterales]",
                "Actinomycetales", "Gemmatales", "Solirubrobacterales", "iii1-15",
                "Sphingomonadales", "Bacillales", "Xanthomonadales", "Burkholderiales",
                "Rhodospirillales", "Myxococcales", "Caulobacterales", "Planctomycetales",
                "Syntrophobacterales", "Gaiellales", "Pirellulales", "Solibacterales",
                "[Saprospirales]", "[Pedosphaerales]", "Acidimicrobiales", "Ellin329",
                "Opitutales", "Verrucomicrobiales", "Sva0725", "SC-I-84", "B07_WMSP1",
                "Sphingobacteriales", "Ellin5290", "WD2101", "Bradyrhizobiaceae",
                "[Chthoniobacteraceae]", "Hyphomicrobiaceae", "Mycobacteriaceae",
                "Isosphaeraceae", "Conexibacteraceae", "Sphingomonadaceae", "Sinobacteraceae",
                "Comamonadaceae", "Rhodospirillaceae", "Haliangiaceae", "Caulobacteraceae",
                "Phyllobacteriaceae", "Streptomycetaceae", "Planctomycetaceae",
                "Syntrophobacteraceae", "Gaiellaceae", "Acetobacteraceae", "Pirellulaceae",
                "Pseudonocardiaceae", "Chitinophagaceae", "Solibacteraceae",
                "auto67_4W", "Gemmataceae", "EB1017", "Microbacteriaceae", "Xanthomonadaceae",
                "Opitutaceae", "Verrucomicrobiaceae", "Paenibacillaceae", "Bradyrhizobium",
                "DA101", "Rhodoplanes", "Mycobacterium", "Methylibium", "Hyphomicrobium",
                "Phenylobacterium", "Devosia", "Planctomyces", "Candidatus Xiphinematobacter",
                "Sphingomonas", "Pseudonocardia", "Chthoniobacter", "Candidatus Solibacter",
                "A17", "Opitutus", "Gemmata", "Paenibacillus")

#load some bacterial data
map <- readRDS(emp_map_clean.path) #environmental metadata
rownames(map) <- map$`#SampleID`
data.list <- readRDS(emp_phylo.level.list_esv.path) #relative abundance 

r2_lev <- list()
for(k in 1:length(data.list)){
  taxa <- t(data.list[[k]])
  n <- ncol(taxa) - 1
  #match the sample names and merge data into one data frame
  d <- data.frame()
  d <- merge(map, taxa, by="row.names")

  #isolate y data (relative abundance)
  y <- d[,(ncol(d)-n):ncol(d)]
  
  y <- y[,colnames(y) %in% bahram_tax]  
  #we need to transform y data to not have any hard zeros or 1s. This is a limitation of beta regression.
  #function to transform [0,1] to (0,1) a la Cribari-Neto & Zeileis 2010, statisticians who have written the most on this topic.
  crib.fun <- function(x){(x * (length(x) - 1) + 0.5) / length(x)}
  y <- as.data.frame(lapply(y,crib.fun))

  #lets fit bacterial models as a function of ph, map, mat, map/mat cv, mdr, and aridity
  preds <- c('map','mat','map_CV','mat_CV','mdr','aridity')
  preds <- c('map','mat')
  x <- d[,preds] #just the metadata we want to look at

  #setup some empty lists.
  model.out <- list()
  predicted.out <- list()
  observed.out <- list()
  #loop across the columns of y.
  for(i in 1:ncol(y)){ 
    #setup data object
    print(i)
    data <- cbind(y[,i],x) 
    data <- data[complete.cases(data),]
    group <- data[,1]
    if (colnames(y)[i]=="NA.") next()
    #setup regression formula.
    form <- as.formula(paste0('group ~ ',paste(colnames(x), collapse="+")))
    #fit model.
    mod <- betareg(form, data = data)
    #save fitted and predicted values.
    model.out[[i]] <- mod
    predicted.out[[i]] <- fitted(mod)
    observed.out[[i]] <- data[,1]
  }
  names(model.out) <- colnames(y)
  names(predicted.out) <- colnames(y)
  names(observed.out) <- colnames(y)

  #collapse predicted and observed data to a data.frame 
  predicted.out <- data.frame(do.call('cbind',predicted.out))
  observed.out <- data.frame(do.call('cbind', observed.out))

  #plot them all!
  #fname <- paste0('/Users/student/Desktop/NEFI_16S_results/',names(data.list)[k],'_',paste0(preds, collapse="_"),'_betareg.png')
  #png(filename = fname, width = 16, height = 7, units = 'in',res=300)
  par(mfrow = c(3,3)) # 5,6 for phyla; 5,10 for all other levels of taxonomy

  par(mar=c(1,1,1,1))
  
  r2_out <- list()
  for(i in 1:ncol(observed.out)){
    mod <- lm(observed.out[,i] ~ predicted.out[,i])
    rsq <- round(summary(mod)$r.squared, 2)
    r2_out[[i]] <- rsq
    plot(observed.out[,i] ~ predicted.out[,i], cex = 0.3, col=d$study_id)
    mtext(colnames(observed.out)[i], side = 3, line = -1, adj = 0.05, cex = 0.8, col = 'purple')
    mtext(paste0('R2 = ',rsq), side = 3, line = -2.4, adj = 0.05, cex = 0.8, col = 'purple')
    #drop 1:1 line
    abline(0,1, lwd = 2)
  }
  title(paste(names(data.list)[k], 'predicted with variables', paste(preds, collapse=', ')), outer=TRUE, line = -0.75)
  r2_lev[[k]] <- mean(unlist(r2_out))
  #dev.off()
}
