#Run NPARC implementation F test from TPP package
runTPP_splineFtest<-function(normData,DF=5,returnModels=FALSE){

  tpptrData<-TPP::tpptrImport(configTable = normData$TPPconfig,data=normData$TPPdata)


  check_longTable <- TPP::tpptrTidyUpESets(tpptrData)
  #TPP fit splines

  SimSplineFits<- TPP::tpptrFitSplines(data = check_longTable,
                                        factorsH1 = c("condition"),
                                        splineDF = DF,
                                        nCores = 6,
                                        returnModels = returnModels)

  SimFtest<-TPP::tpptrFTest(SimSplineFits)

  return(SimFtest)
}
