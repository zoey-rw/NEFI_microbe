#splits apart a semi-comma separated taxonomy list into a table of k/p/o/c/f/g/s
linaeus<-function(tax.vector){
  list<- strsplit(tax.vector,split=";")
  matrix<-stringi::stri_list2matrix(list,byrow=T)
  colnames(matrix) <- c("kingdom","phylum","class","order","family","genus","species")
  table<-as.data.frame(matrix)
  return(table)
}