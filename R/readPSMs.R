
readPSMs<-function(wd){
  h<-list.files(wd,pattern="PSMs.xls") |> as.list()

  PSMs<-purrr::map(h,function(x) readxl::read_xlsx(path=paste0(wd,"/",x))) |> dplyr::bind_rows()

  return(PSMs)
}
