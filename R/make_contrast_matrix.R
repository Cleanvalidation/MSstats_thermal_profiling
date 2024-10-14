
make_contrast_matrix = function(data){


  condition_levels = levels(data$ProteinLevelData$Condition)
  if(any(stringr::str_detect(condition_levels,"Norm"))){
    stop("Remove Norm condition")
  }
  if("131C" %in% unique(data$ProteinLevelData$Channel)){
    stop("Remove Reference Channel 131C")
  }

  #Channel 127C => Temp=44, Channel 130C =>Temp = 63

  #Treated - Vehicle
  levels_Condition = levels(data$ProteinLevelData$Condition)
  ATE = 2/length(levels_Condition)*ifelse(stringr::str_detect(levels_Condition,"treated"),1,-1)

  contrast_matrix = t(as.matrix(ATE))

  colnames(contrast_matrix) = condition_levels
  rownames(contrast_matrix) = c("ATE")
  return(contrast_matrix)
}
