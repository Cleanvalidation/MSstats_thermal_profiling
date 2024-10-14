#Run TPP sigmoids on normalized data
runTPP_sigmoid<-function(Data,DF=future::availableCores()*0.5,NORM=FALSE){
  if(any(names(Data$TPPdata[[1]])=="uniqueID")){
  Data$TPPdata <- furrr::future_map(Data$TPPdata,function(x) x %>% dplyr::select(-dplyr::starts_with("uniqueID"),-dplyr::starts_with("variable"),-dplyr::starts_with("value")))
  }
  tpptrData<-TPP::tpptrImport(configTable = Data$TPPconfig,data=Data$TPPdata)

  if(isTRUE(NORM)){
    TRreqs<-TPP::tpptrDefaultNormReqs()
    TRreqs$fcRequirements<-TRreqs$fcRequirements[1,]
    TRreqs$fcRequirements$thresholdLower<-c(-1)
    TRreqs$fcRequirements$thresholdUpper<- Inf
    normData<-TPP::tpptrNormalize(data=tpptrData,normReqs = TRreqs)
    Fit_sig <- TPP::tpptrCurveFit(data=normData$normData,
                                  resultPath=getwd(),
                                  nCores = DF,
                                  verbose=FALSE,
                                  maxAttempts = 20,
                                  doPlot=FALSE)
  }else{
    Fit_sig <- TPP::tpptrCurveFit(data=tpptrData,
                                  resultPath=getwd(),
                                  nCores = DF,
                                  verbose=FALSE,
                                  maxAttempts = 20,
                                  doPlot=FALSE)
  }
  pval_results<-TPP::tpptrAnalyzeMeltingCurves(Fit_sig)

  return(pval_results)
}
