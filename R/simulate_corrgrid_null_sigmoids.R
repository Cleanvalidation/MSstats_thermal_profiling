simulate_corrgrid_null_sigmoids = function(runs=1000,n=4,error_sd=0.05,rho=seq(0,0.9,0.1)){
  simdata_list = vector(mode="list",length=length(rho))
  simdata_params_list = vector(mode="list",length=length(rho))

  for(i in 1:length(rho)){
    r = rho[i]
    sim = simulate_multiple_null_sigmoids(runs=runs,n=n,error_sd=error_sd,rho=r)
    simdata_list[[i]] = sim$simdata
    simdata_params_list[[i]] = sim$simdata_params
    print(paste0("Generated Null Data for Corr of ",r))
  }

  simdata = dplyr::bind_rows(simdata_list)
  simdata_params = dplyr::bind_rows(simdata_params_list)

  return(list(ProteinLevelData=simdata,params=simdata_params))
}
