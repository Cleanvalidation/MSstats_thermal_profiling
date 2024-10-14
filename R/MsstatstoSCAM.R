#convert MSstats output to fit scam models
MSstatstoSCAM<-function(Result){
  Result<-Result %>% separate(Condition,into=c("temp_ref","treatment"),sep="_") %>% dplyr::mutate(I=Abundance,
                                                                                                  replicate=TechRepMixture,
                                                                                                  Accession=Protein,
                                                                                                  Spectrum.File=Run,
                                                                                                  sample_id=Subject)
  stopifnot(names(Result) %in% c("I","replicate","Accession","Spectrum.File","sample_id","treatment"))
  return(Result)
}
