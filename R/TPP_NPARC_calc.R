#Run TPP or NPARC calculation pipeline, renames data from MSstats to TPP format and then uses method = TPP or NPARC to apply the pipeline
TPP_NPARC_calc<-function(Sim,method="NPARC",DF=5,df.temps,CARRIER=TRUE){
  if(any(names(df.temps)=="Channel")){
    df.temps<-df.temps |> dplyr::rename("temp_ref"="Channel")
  }
  #rename data to MSStatsTMT format
  start=proc.time()
  Sim_All_Data<-list(ProteinLevelData=suppressWarnings(MSstats_format_TO_TPP(Sim,df.temps,CARRIER=CARRIER)))
  end=proc.time()
  print(paste0("Renamed data to match",as.character(method)," format in  ",as.numeric(signif((end-start)[1],2))," seconds"))

  if(method=="NPARC"){

    stopifnot(any(DF %in% c(3,4,5,6,7)))
    #simulate data for spline DF=5
    start=proc.time()
    Sim_NPARC<-runTPP_splineFtest(Sim_All_Data$ProteinLevelData,DF=DF)
    end=proc.time()
    print(paste0("NPARC results calculated in ",as.numeric(signif((end-start)[1],2)),"seconds"))
    Sim_NPARC$unmoderatedFp_val<-1-pf(Sim_NPARC$F_statistic, Sim_NPARC$df1,Sim_NPARC$df2)

    return(Sim_NPARC)
  }else if(method=="TPP"){

    stopifnot(any(DF %in% c(3,4,5,6,7)))
    #simulate data for sigmoid DF=5
    start=proc.time()
    Sim_TPP<-runTPP_sigmoid(Sim_All_Data$ProteinLevelData,DF=DF)
    end=proc.time()
    print(paste0("TPP results calculated in ",as.numeric(signif((end-start)[1],2))," seconds"))

    return(Sim_TPP)
  }
}
