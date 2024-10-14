
#prepares MSstats groupComparisons output to TPP
MSstats_format_TO_TPP<-function(Sim,df.temps,CARRIER=TRUE){
  if(class(Sim)=="list"&any(stringr::str_detect(names(Sim),"TPPdata"))){
  x<-Sim$TPPdata
  }else if(any(names(Sim)=="ProteinLevelData")){
    x<-Sim$ProteinLevelData
  }else{
    x<-Sim
  }
  if(any(names(df.temps)=="Channel")){
    df.temps<-df.temps |> dplyr::rename("temp_ref"="Channel")
  }
  if(isTRUE(CARRIER)){
    if(any(stringr::str_detect(df.temps$temp_ref,"131C"))){
      df.temps<-df.temps[!stringr::str_detect(df.temps$temp_ref,"131C"),]
    }

  }

  df.temps<-df.temps |> dplyr::filter(df.temps$temp_ref %in% unique(x$Channel))
  x<-x |> dplyr::filter(Channel %in% unique(df.temps$temp_ref))
  #rename_TPP
  if(any(names(x)=="Protein")&!any(names(x)=="uniqueID")){
    x$uniqueID<-x$Protein
  }
  if(any(names(x)=="uniqueID")&!any(names(x)=="gene_name")){
    x$gene_name<-x$uniqueID
  }
  if(any(names(x)=="Protein")&!any(names(x)=="gene_name")){
    x$gene_name<-x$Protein
  }
  if(any(names(x)=="treatment")&!any(names(x)=="Condition")){
    x$Condition<-x$treatment
  }
  if(any(names(x)=="dataset")&!any(names(x)=="Condition")){
    x$Condition<-x$dataset
  }
  if(!any(names(x)=="uniqueID")&any(names(x)=="Protein")){
    x$uniqueID<-x$Protein
  }
  if(any(names(x)=="Channel")&!any(names(x)=="temp_ref")){
    x$temp_ref<-x$Channel
  }
  if(any(names(x)=="value")&!any(names(x)=="I")){
    x$I<-x$value
  }
  if(any(names(x)=="Abundance")&!any(names(x)=="I")){
    x$I<-x$Abundance
  }
  if(any(stringr::str_detect(x$Condition,"_"))){
    x$Condition<-stringr::str_extract(stringr::str_to_lower(x$Condition),"[[:lower:]]+")
  }
  if(any(names(x)=="Subject")&!any(names(x)=="sample_id")){
    x$sample_id<-x$Subject
  }
  if(any(names(x)=="Mixture")&!any(names(x)=="Experiment")){
    x$Experiment<-x$Mixture
  }
  if(any(names(x)=="Mixture")&!any(names(x)=="sample_id")){
    x$sample_id<-x$Mixture
  }
  if(any(names(x)=="Experiment")){
    TechReps<-data.frame(Experiment=as.factor(unique(x$Experiment)),TechRepMixture=rep(c(1,2),2))
    x$Experiment<-as.factor(x$Experiment)
    if(any(names(x)=="TechRepMixture")){
      x <-x |> dplyr::select(-TechRepMixture)
    }
    x<-x |> dplyr::inner_join(TechReps,by="Experiment")
  }else if(any(names(x)=="Subject")){
    TechReps<-data.frame(Experiment=as.factor(unique(x$Subject)),TechRepMixture=rep(c(1,2),2))
    x<-x |> dplyr::rename(Experiment=Subject)
    if(any(names(x)=="TechRepMixture")){
      x <-x |> dplyr::select(-TechRepMixture)
    }
    x<-x |> dplyr::inner_join(TechReps,by="Experiment")
  }else if(any(names(x)=="Condition")){
    if(any(stringr::str_detect(x$Condition,"_"))){
      x$Experiment<-stringr::str_extract(stringr::str_to_lower(x$Condition),"[[:lower:]]+")
    }
  }


  #stopifnot(all(stringr::str_detect(names(x),c("uniqueID","gene_name","sample_id","TechRepMixture","I","temp_ref"))))
  TPP_Cliff<-dplyr::bind_rows(x) |> dplyr::select(uniqueID,gene_name,sample_id,TechRepMixture,I,temp_ref,Experiment) |> dplyr::filter(!is.na(I))
  TPP_Cliff_pivot<-tidyr::pivot_wider(
    TPP_Cliff,
    id_cols = NULL,
    names_from = temp_ref,
    names_prefix = "rel_fc_",
    names_sep = "_",
    names_repair = "minimal",
    values_from = I,
    values_fn=unique

  )

  TPP_Cliff<-TPP_Cliff_pivot |> dplyr::distinct()
  check<-names(TPP_Cliff)
  check1<-check[stringr::str_detect(check,"[:digit:][:upper:]")]
  #column numbers that have reporter ion data
  data2<-which(check %in% check1)
  #replace C or N with L and H
  check1<-stringr::str_replace(check1,"C","H")
  check1<-stringr::str_replace(check1,"N","L")
  #replace names
  check[data2]<-check1
  names(TPP_Cliff)<-check
  x1<-x |> dplyr::select(uniqueID,gene_name,Condition,sample_id,TechRepMixture)
  TPP_Cliff<-TPP_Cliff |> dplyr::inner_join(x1)
  TPP_Cliff$treatment<-stringr::str_extract(stringr::str_to_lower(TPP_Cliff$Condition),"[[:lower:]]+")
  TPP_Cliff$Condition<-ifelse(TPP_Cliff$treatment=="vehicle","Vehicle","Treatment")
  TPP_Cliff$Experiment<-ifelse(TPP_Cliff$Condition=="Treatment",
                               paste0("Treatment_",TPP_Cliff$TechRepMixture),paste0(TPP_Cliff$Condition,"_",TPP_Cliff$TechRepMixture))

  TPP_Cliff$ComparisonVT1<-NA
  TPP_Cliff$ComparisonVT2<-NA

  TPP_Cliff$ComparisonVT1<-ifelse(TPP_Cliff$TechRepMixture==1,"x","")
  TPP_Cliff$ComparisonVT2<-ifelse(TPP_Cliff$TechRepMixture==2,"x","")

  check1<-check[stringr::str_detect(check,"rel_fc_[[:digit:]]+|rel_fc_[[:digit:]]+[:upper:]")]
  #column numbers that have reporter ion data
  data2<-which(check %in% check1)

  config<-dplyr::bind_rows(TPP_Cliff)
  names(config)<-stringr::str_replace(names(config),"rel_fc_","")
  check<-c(config |> dplyr::select(Experiment,Condition,ComparisonVT1,ComparisonVT2),config[data2])
  temp_ref<-stringr::str_replace(df.temps$temp_ref,"C","H")
  temp_ref<-stringr::str_replace(temp_ref,"N","L")

  temps<-df.temps |> dplyr::mutate(temp_ref=temp_ref) |> dplyr::distinct()
  temps<-tidyr::pivot_wider(temps,names_from=temp_ref,values_from=temperature)
  temps<-purrr::map_dfr(seq_len(nrow(TPP_Cliff)), ~temps)
  temp_ref<-stringr::str_replace(names(temps),"C","H")
  temp_ref<-stringr::str_replace(temp_ref,"N","L")

  names(temps)<-temp_ref
  TPP_Cliff<-cbind(TPP_Cliff,temps)
  #keep two replicates
  TPP_Cliff<-TPP_Cliff |> dplyr::filter(ComparisonVT1=="x" | ComparisonVT2=="x")
  TPP_Cliff$qssm<-as.integer(5)
  TPP_Cliff$qupm<-as.integer(10)

  temp_ref<-stringr::str_replace(df.temps$temp_ref,"C","H")
  temp_ref<-stringr::str_replace(temp_ref,"N","L")

  TPPconfig<-TPP_Cliff |>
    dplyr::select(Experiment,Condition,ComparisonVT1,ComparisonVT2,tidyr::starts_with("1")) |>
    dplyr::distinct()|>
    dplyr::mutate(Experiment=as.character(Experiment)) |>
    dplyr::arrange(Experiment)

  TPPdata<-TPP_Cliff |>
    dplyr::select(uniqueID,gene_name,Experiment,qssm,qupm,tidyr::starts_with("rel_fc")) |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(Experiment)) |>
    dplyr::arrange(Experiment) |>
    dplyr::mutate(Experiment=as.character(Experiment)) |>
    dplyr::group_split(gene_name)
  TPPdata<-TPPdata |> purrr::keep(function(x) length(unique(x$Experiment))==4)

  TPPdata<-dplyr::bind_rows(TPPdata) |> dplyr::group_split(Experiment)
  names(TPPdata)<-unique(TPPconfig$Experiment)
  resultPath<-file.path(getwd())

  TPPdata<-purrr::map(TPPdata,function(x) x |> dplyr::filter(!is.na(rel_fc_126)))
  TPPdata1<-purrr::map(TPPdata,function(x) x |> dplyr::filter(!is.na(gene_name)) |> reshape2::melt())
  check<-purrr::map2(TPPdata,TPPdata1,function(x,y)y |> dplyr::right_join(x))
  check<-lapply(check,function(x) x |> dplyr::mutate(uniqueID=as.factor(uniqueID),
                                                      gene_name=uniqueID))
  out<-list(TPPconfig,check)
  names(out)<-c("TPPconfig","TPPdata")
  return(out)

}
