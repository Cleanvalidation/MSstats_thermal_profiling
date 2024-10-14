Annotation2MSstatsTMT<-function(input,solvent="DMSO",temps,reference){
  #pivot_longer
  input_long<-input %>%
    tidyr::pivot_longer(.,cols=names(input)[stringr::str_detect(names(input),"126|[[:digit:]]+N|[[:digit:]]+C")],
                        names_to="id",
                        values_to="value")

  input_long<-input_long %>%
    dplyr::mutate(temp_ref=stringr::str_extract(id,"126|[[:digit:]]+N|[[:digit:]]+C|[[:digit:]]+|131"))
  names(input_long)<-stringr::str_replace_all(names(input_long)," ","_")

  input_long$replicate<-stringr::str_extract(input_long$Spectrum_File,"[:upper:]+[[:digit:]]+_|[:lower:][[:digit:]]+_")
  input_long$replicate<-stringr::str_extract(input_long$replicate,"[[:digit:]]+")
  input_long$TechRepMixture<-input_long$replicate
  input_long$Mixture<-input_long$File_ID
  if(any(stringr::str_detect(input_long$File_ID,"."))){
    input_long$Fraction<-stringr::str_remove(input_long$File_ID,"F[[:digit:]]+.")
  }



  if(any(names(input_long)=="Intensity")){
    input_long<-input_long %>% dplyr::select(-Intensity)
  }
  if(!is.na(reference)){
  input_long<-input_long %>%
    dplyr::rename("Protein"="Master_Protein_Accessions",
                  "Run"="Spectrum_File",
                  "Abundance"="value",
                  "Channel"="temp_ref") %>%
    dplyr::mutate(PrecursorCharge=Charge,
                  Mixture=Run,
                  BioReplicate=paste0(Run,"_",Channel),
                  Condition=ifelse(Channel==reference,"Norm",paste0(Channel,"_",ifelse(stringr::str_detect(Run,solvent),"vehicle","treated"))),
                  Group=Condition,
                  Channel=as.factor(Channel))
  }else{
    input_long<-input_long %>%
      dplyr::rename("Protein"="Master_Protein_Accessions",
                    "Run"="Spectrum_File",
                    "Abundance"="value",
                    "Channel"="temp_ref") %>%
      dplyr::mutate(PrecursorCharge=Charge,
                    Mixture=Run,
                    BioReplicate=paste0(Run,"_",Channel),
                    Condition=paste0(Channel,"_",ifelse(stringr::str_detect(Run,solvent),"vehicle","treated")),
                    Group=Condition,
                    Channel=as.factor(Channel))
}
  output<-input_long %>%
    dplyr::inner_join(temps) %>%
    dplyr::select(Run,Fraction,TechRepMixture,Channel,Condition,Mixture,BioReplicate) %>%
    dplyr::distinct()
  return(output)
}
