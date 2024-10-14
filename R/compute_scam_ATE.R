compute_scam_ATE = function(data){
  model = scam::scam(Abundance ~ s(temperature,by=Condition,bs="mpd") + Condition + s(Subject,bs="re"),data = data)
  ATE = marginaleffects::tidy(marginaleffects::comparisons(model,var="Condition"))
  return(ATE)
}
