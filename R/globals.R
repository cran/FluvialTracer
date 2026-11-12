# R/globals.R
utils::globalVariables(c(
  # ggplot2 aesthetics and plot helpers
  "LD1", "LD2", "LD3", "Sources", "variable", "value", "n", "Group",
  "id", "prop", "x", "y", "z", "groups", "xvar", "yvar", "varname",
  "angle", "hjust", "..scaled..", "red", "ID", "GOF", "tracer",
  "Dw1", "Dw2", "Dw3", "Dw4", "Dw5", "mean_val", "sd_val", "s_max",
  "t_max", "s_min", "t_min", "row_num", "after_stat",
  
  # ggbiplot internals (third-party embedded code - scales/plyr)
  "muted", "ddply",
  
  # Package-internal functions called from other files
  # (R CMD check sees these as globals in some call contexts)
  "is_averaged", "unmix",
  # cluster package
  "cluster", "dist", "label",
  #Ternary
  "complete.cases",
  #pull analysis
  "id", "Tracer", "Source", "Weight", "Median_Weight", "SD_Weight",
  
  # Package dataset used unqualified inside plotting functions
  "fluvial_palette",
  
  # Diff mix variables (added to clear check notes)
  "Proportion", "Difference", "Difference_Plot", "Label_HTML", 
  "Sort_Val", "Direction", "Label_Orig", ":="
))