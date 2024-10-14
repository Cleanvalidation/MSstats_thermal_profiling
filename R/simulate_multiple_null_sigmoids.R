
simulate_multiple_null_sigmoids = function(runs=1000,n=4,error_sd=0.05,rho=0){
  simdata_list = vector(mode="list",length=runs)
  simdata_params_list = vector(mode="list",length=runs)

  for(i in 1:runs){
    Protein = paste("Null_Corr",rho,"Sim",i,sep="_")

    sim = simulate_null_sigmoid(n=n,error_sd=error_sd,rho=rho)
    simdata_list[[i]] = sim$simdata
    simdata_list[[i]]$Protein = Protein

    simdata_params_list[[i]] = sim$sample_params
    simdata_params_list[[i]]$Protein = Protein
  }

  simdata = dplyr::bind_rows(simdata_list)
  simdata_params = dplyr::bind_rows(simdata_params_list)

  return(list(simdata=simdata,simdata_params=simdata_params))
}
