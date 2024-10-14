
simulate_null_sigmoid = function(n=4,error_sd=0.05,rho=0){
  stopifnot(n%%2==0)
  temperature = c(37,41,44,47,50,53,56,59,63,67)
  Ntemps = length(temperature)


  error_var = error_sd^2
  error_var_vec = rep(error_var,Ntemps)

  R = diag(Ntemps)
  R[R==0] = rho
  S = cor2cov(R,error_var_vec)

  sample_idx = sample(1:nrow(null_sigmoid_posterior_df),1)
  sample_params = null_sigmoid_posterior_df[sample_idx,] #get posterior params

  abundance = Cetsa(sample_params$p,sample_params$k,sample_params$m,t=temperature)

  simdata_wide = MASS::mvrnorm(n=n,mu=abundance,Sigma=S)
  colnames(simdata_wide) = temperature
  Subject = paste0("F",seq(1,n,1))
  simdata_wide = dplyr::as_tibble(simdata_wide) |> dplyr::mutate(Subject = Subject,Condition=c(rep("treated",n/2),rep("vehicle",n/2)))

  simdata = simdata_wide |> tidyr::pivot_longer(cols=matches("[[:digit:]]+"),names_to="temperature",values_to="Abundance") |>
    dplyr::mutate(temperature=as.numeric(temperature)) |>
    dplyr::inner_join(Channel2Temps_HumanData,by="temperature") |>
    dplyr::mutate(Group = paste(Channel,Condition,sep="_"),
                  Subject = forcats::as_factor(Subject),
                  Condition = forcats::as_factor(Condition),
                  Condition = forcats::fct_relevel(Condition,"vehicle","treated")) |>
    dplyr::arrange(temperature)

  sample_params = sample_params |> dplyr::select(p,k,m) |> dplyr::mutate(sigma=error_sd,rho=rho)

  return(list(simdata=simdata,sample_params=sample_params))


  return(simdata)
}
