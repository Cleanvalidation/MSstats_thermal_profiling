
Sim2MSstatsTMT<-function(Sim){
  if(any(names(Sim)=="ProteinLevelData")){
    Sim<-Sim$ProteinLevelData |> as.data.frame()
  }
  x<-Sim |> dplyr::mutate(
    BioReplicate=paste0(Protein,"_",Channel),
                           Mixture=Protein,
                           Run="Simulation",
                           Experiment=Subject,
                           Fraction=1,
                           Condition=paste0(Channel,"_",Condition))

  TechReps<-data.frame(Experiment=as.factor(unique(x$Subject)),TechRepMixture=rep(c(1,2),2))

  if(any(names(x)=="TechRepMixture")){
    x <-x |> dplyr::select(-TechRepMixture)
  }
  x<-x |> dplyr::inner_join(TechReps,by="Experiment")|> as.data.frame()
  return(x)
}
