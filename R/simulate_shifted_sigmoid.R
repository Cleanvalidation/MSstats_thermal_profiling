simulate_shifted_sigmoid = function(n=4,error_sd=0.08,rho=0){
  stopifnot(n%%2==0) #we restrict to balanced scenarios only
  temperature = c(37,41,44,47,50,53,56,59,63,67)
  Ntemps = length(temperature)

  error_var = error_sd^2
  error_var_vec = rep(error_var,Ntemps)

  R = diag(Ntemps)
  R[R==0] = rho
  S = cor2cov(R,error_var_vec) #assume same cov matrix for treated/vehicle

  sample_idx = sample(1:nrow(truepositive_sample_sigmoid_posterior_df),1)
  sample_params = truepositive_sample_sigmoid_posterior_df[sample_idx,] #get posterior params

  abundance_treated = Cetsa(p=sample_params$p_Tmt,k=sample_params$k_Tmt,m=sample_params$m_Tmt,t=temperature)
  abundance_vehicle = Cetsa(p=sample_params$p_Ctrl,k=sample_params$k_Ctrl,m=sample_params$m_Ctrl,t=temperature)

  simdata_wide_treated = MASS::mvrnorm(n=n/2,mu=abundance_treated,Sigma=S)
  colnames(simdata_wide_treated) = temperature


  simdata_wide_vehicle = MASS::mvrnorm(n=n/2,mu=abundance_vehicle,Sigma=S)
  colnames(simdata_wide_vehicle) = temperature

  Subject = paste0("F",seq(1,n,1))

  simdata_wide = rbind(simdata_wide_treated,simdata_wide_vehicle) |> dplyr::as_tibble() |>
    dplyr::mutate(Subject = Subject,Condition=c(rep("treated",n/2),rep("vehicle",n/2)))

  simdata = simdata_wide |> tidyr::pivot_longer(cols=matches("[[:digit:]]+"),names_to="temperature",values_to="Abundance") |>
    dplyr::mutate(temperature=as.numeric(temperature)) |>
    dplyr::inner_join(Channel2Temps_HumanData,by="temperature") |>
    dplyr::mutate(Group = paste(Channel,Condition,sep="_"),
                  Subject = forcats::as_factor(Subject),
                  Condition = forcats::as_factor(Condition),
                  Condition = forcats::fct_relevel(Condition,"vehicle","treated")) |>
    dplyr::arrange(temperature)

  sample_params = sample_params |> dplyr::select(p_Tmt,k_Tmt,m_Tmt,p_Ctrl,k_Ctrl,m_Ctrl,ATE) |> dplyr::mutate(sigma=error_sd,rho=rho)

  return(list(simdata=simdata,sample_params=sample_params))
}
