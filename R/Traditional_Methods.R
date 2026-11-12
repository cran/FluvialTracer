#' Range test for tracer pre-selection
#'
#' Removes tracer properties whose mixture value falls outside the envelope
#' defined by the minimum and maximum values observed across all source groups.
#' Tracers failing this test cannot be explained by any linear combination of
#' the measured sources and should be excluded before unmixing.
#'
#' @param data Data frame in the standard FluvialTracer format (sources + mixture).
#'
#' @return A data frame retaining only tracer columns that pass the range test.
#'
#' @export
#'
rangeTest <- function(data) {
  # reorder factor levels in order of appearance
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))

  # read groups (second column)
  groups <- data[, 2]
  
  # replace all groups with 'source' except the last one
  levels(groups)[1:nlevels(groups) - 1] <- "source"
  
  # compute min values for mixture and source groups
  min <- aggregate(data[, 3:ncol(data)], list(groups), min)
  n_min <- c("s_min", "t_min")
  mint <- as.data.frame(t(min[, -1]))
  colnames(mint) <- n_min
  
  # compute max values for mixture and source groups
  max <- aggregate(data[, 3:ncol(data)], list(groups), max)
  n_max <- c("s_max", "t_max")
  maxt <- as.data.frame(t(max[, -1]))
  colnames(maxt) <- n_max
  
  # merge min and max values
  ranges1 <- merge(maxt, mint, by.x = "row.names", by.y = "row.names")
  
  # filter properties out of range
  ranges <- subset(ranges1, s_max < t_max | s_min > t_min)
  rows <- as.vector(ranges$Row.names)
  data <- data[, !(names(data) %in% rows)]
  
  message("Attention-> ", nrow(ranges), " variables from a total of ", nrow(ranges1), 
          " were removed: ", crayon::red(crayon::bold(paste(ranges[, 1], collapse = " "))), " .")
  message(" The variable/variables that remains in your dataset is/are: ", 
          crayon::green(crayon::bold(paste(names(data[, 3:ncol(data)]), collapse = " "))), " .")
  
  return(data)
}



#' Kruskal-Wallis test for tracer pre-selection
#'
#' Retains only tracer properties that show statistically significant
#' differences across source groups according to the non-parametric
#' Kruskal-Wallis rank-sum test. Tracers with no discriminant power add noise
#' to the unmixing model and should be removed before unmixing.
#'
#' @param data Data frame in the standard FluvialTracer format (sources + mixture).
#' @param pvalue Numeric. Significance threshold. Default 0.05.
#'
#' @return A data frame retaining only tracer columns with p-value below pvalue.
#'
#' @export
#'
#'
KWTest <- function(data, pvalue = 0.05) {
  # reorder factor levels in order of appearance
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  
  # read groups (second column)
  groups <- data[, 2]
  
  # asume last group is mixtures
  mixture <- levels(groups)[nlevels(groups)]
  
  # read sources
  sources <- data[!groups == mixture, ]
  
  # remove mixture level
  s_groups <- droplevels(sources[, 2])
  
  # extract properties
  var_s <- sources[, 3:ncol(sources)]
  var <- colnames(var_s)
  
  # perform KW test
  KW_p <- c()
  for (i in 1:ncol(var_s)) {
    KW_p <- append(KW_p, kruskal.test(var_s[, i], s_groups)$p.value)
  }
  KW <- data.frame(var, KW_p)
  
  # remove properties under pvalue
  KWOFF <- subset(KW, KW_p > pvalue)
  varOFFKW <- as.vector(KWOFF$var)
  total_KW <- data[, !(names(data) %in% varOFFKW)]
  
  message("Attention-> ", nrow(KWOFF), " variables from a total of ", ncol(data[3:ncol(data)]), 
          " were removed: ", crayon::red(crayon::bold(paste(varOFFKW, collapse = " "))), " .")
  message(" The variable/variables that remains in your dataset is/are: ", 
          crayon::green(crayon::bold(paste(names(total_KW[3:ncol(total_KW)]), collapse = " "))), " .")
  
  return(total_KW)
}


#' Discriminant Function Analysis for tracer pre-selection
#'
#' Performs stepwise forward variable selection using the Wilks Lambda
#' criterion. Tracers are added one at a time when they improve inter-source
#' discrimination beyond the specified significance threshold.
#'
#' @param data Data frame in the standard FluvialTracer format (sources + mixture).
#' @param niveau Numeric. Significance level for the F-test at each step. Default 0.1..
#'
#' @return A list with two elements: (1) a data frame retaining only the tracer
#'   columns selected by the stepwise procedure, and (2) the object returned by
#'   \code{klaR::greedy.wilks} with the details of the selection.
#'
#' @export
#'
DFATest <- function(data, niveau = 0.1) {
  # reorder factor levels in order of appearance
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  
  # read groups (second column)
  groups <- data[, 2]
  
  # asume last group is mixtures
  mixture <- levels(groups)[nlevels(groups)]
  
  # read sources
  sources <- data[!groups == mixture, ]
  
  # remove mixture level
  s_groups <- droplevels(sources[, 2])
  
  # extract properties
  var_A <- sources[3:ncol(sources)]
  
  # assign groups
  var_A$groups <- as.factor(s_groups)
  
  # perform DFA
  DFA <- greedy.wilks(groups ~ ., data = var_A, niveau = niveau)
  
  # filter selected properties
  var <- as.vector(DFA[[1]][, 1])
  DFA_OFF1 <- data[, 3:ncol(data)]
  DFA_OFF1 <- DFA_OFF1[, !(names(DFA_OFF1) %in% var)]
  var_OFF.DFA <- names(DFA_OFF1)
  total_DFA <- data[, !(names(data) %in% var_OFF.DFA)]
  
  
  message("Attention-> ", ncol(DFA_OFF1), " variables were removed from a total of ", 
          ncol(data[3:ncol(data)]), ": ", crayon::red(crayon::bold(paste(var_OFF.DFA, collapse = " "))), " .")
  message(" The variable/variables that remains in your dataset is/are: ", 
          crayon::green(crayon::bold(paste(names(total_DFA[3:ncol(total_DFA)]), collapse = " "))), " .")
  
  return(list(total_DFA, DFA))
  
}

