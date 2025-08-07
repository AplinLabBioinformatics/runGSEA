
#' Check if Java is running
#'
#' Check if Java is running to determine if GSEA should be performed.
#'
#' @return An R logical value (`TRUE` or `FALSE`)
#' @export
#'
#' @examples
#' checkjava()

checkjava <- function(){
  temp_sysCheck <- system("wmic path Win32_PerfFormattedData_PerfProc_Process get Name,PercentProcessorTime",
                          intern = TRUE)
  temp_df_check <- do.call(rbind, lapply(strsplit(temp_sysCheck, " "), function(x) {
    x <- x[x != ""]
    data.frame(process = x[1], cpu = x[2])
  }))
  temp_javaVal <- temp_df_check[grepl("java", temp_df_check$process), "cpu"]
  return(temp_javaVal)
}

#' Convert Filename from R Format to Command Prompt Format
#'
#' Converts an R string representing a file path to a format compatible with Command Prompt handling.
#'
#' @param strIn A character string representing a file path.
#'
#' @return A Command Prompt path string.
#' @export
#'
#' @examples
#' rtocppath("C:/Users/Example/Documents/file.txt")

rtocppath <- function(strIn) {
  if ((grepl("[:][/]", strIn) == TRUE) & (grepl("\\\\", strIn) == FALSE)) {
    strIn <- gsub("[/]", "\\\\\\", strIn)
  }
  paste0("\"", strIn, "\"")
  return(strIn)
}

#' Convert R Logical to Command Prompt Logical
#'
#' Converts an R logical value (`TRUE` or `FALSE`) to a format compatible with command prompt.
#'
#' @param boolIn A logical value in R (`TRUE` or `FALSE`).
#'
#' @return A logical value in Command Prompt (`true` or `false`).
#' @export
#'
#' @examples
#' rtocpbool(TRUE)
#' rtocpbool(FALSE)

rtocpbool <- function(boolIn) {
  if (is.logical(boolIn)) {
    if (boolIn == T) {
      boolIn <- "true"
    } else {
      boolIn <- "false"
    }
  }
  return(boolIn)
}

#' Run Gene Set Enrichment Analysis (GSEA) via Command Prompt
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) using expression data and class labels.
#' It supports various customization options including scoring schemes, permutation settings, and 
#' output formats. Defaults are set to current GSEA defaults.
#'
#' @param gctExprFile Path to the GCT file containing expression data.
#' @param clsLabelsFile Path to the CLS file containing class labels.
#' @param clsComparison Comparison string specifying the phenotype labels to compare.
#' @param GSEArepLbl Label for the GSEA report.
#' @param GSEAresOutDir Directory to save GSEA results.
#' @param gmtFile Path to the GMT file containing gene sets.
#' @param pathToBat Optional path to the GSEA batch file.
#' @param ChipFile Optional path to the chip annotation file.
#' @param gmtAbbr Optional abbreviation for the GMT file.
#' @param scoring_scheme Scoring scheme to use. Choices: `"weighted"`, `"classic"`, `"weighted_p2"`, 
#' `"weighted_p1.5"`.
#' @param collapse Collapse mode. Choices: `"Collapse"`, `"No_Collapse"`, `"Remap_Only"`.
#' @param metric Metric for ranking genes. Choices: `"Signal2Noise"`, `"tTest"`, `"Cosine"`, 
#' `"Euclidean"`, `"Manhattan"`, `"Pearson"`, `"Spearman"`, `"Ratio_of_Classes"`, `"Diff_of_Classes"`, 
#' `"log2_Ratio_of_Classes"`.
#' @param modeProbe Probe mode. Choices: `"Max_probe"`, `"Median_of_probe"`, `"Mean_of_probe"`, 
#' `"Sum_of_probes"`, `"Abs_max_of_probes"`.
#' @param normMode Normalization mode. Choices: `"meandiv"`, `"none"`.
#' @param nperm Number of permutations to perform.
#' @param permuteOn Permutation type. Choices: `"gene_set"`, `"phenotype"`.
#' @param rnd_type Randomization type. Choices: `"no_balance"`, `"equalize_and_balance"`.
#' @param sortOn Sorting basis. Choices: `"real"`, `"abs"`.
#' @param orderOn Sorting order. Choices: `"descending"`, `"ascending"`.
#' @param create_gcts Logical; whether to create GCT files.
#' @param create_svgs Logical; whether to create SVG plots.
#' @param include_only_symbols Logical; whether to include only gene symbols.
#' @param make_sets Logical; whether to generate gene sets.
#' @param medianOn Logical; whether to use median for collapsing probes.
#' @param save_rnd_lists Logical; whether to save randomized gene lists.
#' @param zip_report Logical; whether to zip the GSEA report.
#' @param numMarkers Number of top markers to include.
#' @param plot_top_x Number of top gene sets to plot.
#' @param rnd_seed Random seed. Choices: `"timestamp"`, `"149"`.
#' @param set_max Maximum gene set size.
#' @param set_min Minimum gene set size.
#' @param showOutputOnConsole Logical; whether to show output on console.
#'
#' @return No return value. Results are saved to the specified output directory.
#' @export
#'
#' @examples
#' runGSEA("data.gct", "labels.cls", "Control_vs_Treatment", "GSEA_Rep1",
#'         "results/", "genesets.gmt")


runGSEA <- function(gctExprFile, clsLabelsFile, clsComparison,
                    GSEArepLbl, GSEAresOutDir, gmtFile, pathToBat = NA,
                    ChipFile = NA, gmtAbbr = NA, scoring_scheme = "weighted",
                    collapse = "Collapse", metric = "Signal2Noise", 
                    modeProbe = "Max_probe", normMode = "meandiv", nperm = 1000,
                    permuteOn = "gene_set", rnd_type = "no_balance",
                    sortOn = "real", orderOn = "descending", create_gcts = FALSE,
                    create_svgs = FALSE, include_only_symbols = TRUE,
                    make_sets = TRUE, medianOn = FALSE, save_rnd_lists = FALSE,
                    zip_report = FALSE, numMarkers = 100, plot_top_x = 25,
                    rnd_seed = "timestamp", set_max = 500, set_min = 15,
                    showOutputOnConsole = TRUE) {
  
  filePaths <- list(gctExprFile, clsLabelsFile, GSEAresOutDir, gmtFile, ChipFile)
  names(filePaths) <- c("gctExprFile", "clsLabelsFile", "GSEAresOutDir", "gmtFile", "ChipFile")
  
  checkFilePathsExist <- lapply(filePaths, file.exists)
  
  if(any(unlist(checkFilePathsExist) == FALSE)){
    whichCheck <- names(filePaths)[unlist(checkFilePathsExist) == FALSE]
    stop(paste0("One or more of the files supplied to runGSEA do not exist. Check: \n", 
                paste(whichCheck, collapse = "\n")))
  }
  
  checkFilePaths <- lapply(filePaths, rtocppath)
  
  logicalsList <- list(create_gcts, create_svgs, include_only_symbols,
                       make_sets, medianOn, save_rnd_lists, zip_report)
  names(logicalsList) <- c("create_gcts", "create_svgs", "include_only_symbols",
                           "make_sets", "medianOn", "save_rnd_lists", "zip_report")
  
  if(all(unlist(logicalsList) %in% c(1L, 0L)) == FALSE){
    whichCheck <- names(logicalsList)[! unlist(logicalsList) %in% c(1L, 0L)]
    stop(paste0("Non-logical value(s) supplied when expecting logical value(s). Check: \n", 
                paste(whichCheck, collapse = "\n")))
  }
  
  checkLogicals <- lapply(logicalsList, rtocpbool)
  
  scoring_scheme <- match.arg(scoring_scheme, 
                              choices = c("weighted", "classic", "weighted_p2", "weighted_p1.5"))
  collapse <- match.arg(collapse, choices = c("Collapse", "No_Collapse", "Remap_Only"))
  metric <- match.arg(metric, 
                      choices = c("Signal2Noise", "tTest", "Cosine", "Euclidean",
                                  "Manhattan", "Pearson", "Spearman", "Ratio_of_Classes",
                                  "Diff_of_Classes", "log2_Ratio_of_Classes"))
  modeProbe <- match.arg(modeProbe, choices = c("Max_probe", "Median_of_probe", 
                                                "Mean_of_probe", "Sum_of_probes", 
                                                "Abs_max_of_probes"))
  normMode <- match.arg(normMode, choices = c("meandiv", "none"))
  permuteOn <- match.arg(permuteOn, choices = c("gene_set", "phenotype"))
  rnd_type <- match.arg(rnd_type, choices = c("no_balance", "equalize_and_balance"))
  sortOn <- match.arg(sortOn, choices = c("real", "abs"))
  orderOn <- match.arg(orderOn, choices = c("descending", "ascending"))
  rnd_seed <- match.arg(rnd_seed, choices = c("timestamp", "149"))
  
  clsString = paste0(checkFilePaths["clsLabelsFile"], "#", clsComparison)
  
  commandRun = paste("gsea-cli.bat GSEA",
                     "-res", checkFilePaths["gctExprFile"],
                     "-cls", clsString,
                     "-gmx", checkFilePaths["gmtFile"],
                     "-collapse", collapse,
                     "-chip", checkFilePaths["ChipFile"], 
                     "-mode", modeProbe,
                     "-norm", normMode,
                     "-nperm", as.character(nperm),
                     "-permute", permuteOn, 
                     "-rnd_type", rnd_type,
                     "-scoring_scheme", scoring_scheme,
                     "-rpt_label", GSEArepLbl,
                     "-metric", metric,
                     "-sort", sortOn,
                     "-order", orderOn, 
                     "-create_gcts", checkLogicals["create_gcts"],
                     "-create_svgs", checkLogicals["create_svgs"],
                     "-include_only_symbols",  checkLogicals["include_only_symbols"],
                     "-make_sets", checkLogicals["make_sets"],
                     "-median", checkLogicals["medianOn"],
                     "-save_rnd_lists", checkLogicals["save_rnd_lists"],
                     "-zip_report", checkLogicals["zip_report"],
                     "-num", as.character(numMarkers),
                     "-plot_top_x", as.character(plot_top_x),
                     "-rnd_seed", rnd_seed,
                     "-set_max", as.character(set_max),
                     "-set_min", as.character(set_min),
                     "-out", checkFilePaths["GSEAresOutDir"], collapse = " ")
  
  setwd(pathToBat)
  
  system(commandRun, intern = TRUE,
         ignore.stdout = FALSE, ignore.stderr = FALSE,
         wait = TRUE, input = NULL, show.output.on.console = showOutputOnConsole,
         minimized = FALSE, invisible = TRUE, timeout = 0)
  
  return(commandRun)
}


#' Run Preranked Gene Set Enrichment Analysis (GSEAPrerank) via Command Prompt
#'
#' This function performs a preranked Gene Set Enrichment Analysis (GSEA) using a ranked gene list.
#' It supports various customization options including scoring schemes, probe handling, normalization,
#' and output formatting. Defaults are set to current GSEA defaults
#'
#' @param rnkFile Path to the RNK file containing the ranked gene list.
#' @param GSEArepLbl Label for the GSEA report.
#' @param GSEAresOutDir Directory to save GSEA results.
#' @param gmtFile Path to the GMT file containing gene sets.
#' @param pathToBat Optional path to the GSEA batch file.
#' @param altDelim Optional alternate delimiter for the RNK file.
#' @param ChipFile Optional path to the chip annotation file.
#' @param gmtAbbr Optional abbreviation for the GMT file.
#' @param scoring_scheme Scoring scheme to use. Choices: `"weighted"`, `"classic"`, `"weighted_p2"`, 
#' `"weighted_p1.5"`.
#' @param collapse Collapse mode. Choices: `"Collapse"`, `"No_Collapse"`, `"Remap_Only"`.
#' @param modeProbe Probe mode. Choices: `"Max_probe"`, `"Median_of_probe"`, `"Mean_of_probe"`, 
#' `"Sum_of_probes"`, `"Abs_max_of_probes"`.
#' @param normMode Normalization mode. Choices: `"meandiv"`, `"none"`.
#' @param nperm Number of permutations to perform.
#' @param create_svgs Logical; whether to create SVG plots.
#' @param include_only_symbols Logical; whether to include only gene symbols.
#' @param make_sets Logical; whether to generate gene sets.
#' @param zip_report Logical; whether to zip the GSEA report.
#' @param numMarkers Number of top markers to include.
#' @param plot_top_x Number of top gene sets to plot.
#' @param rnd_seed Random seed. Choices: `"timestamp"`, `"149"`.
#' @param set_max Maximum gene set size.
#' @param set_min Minimum gene set size.
#' @param showOutputOnConsole Logical; whether to show output on console.
#'
#' @return No return value. Results are saved to the specified output directory.
#' @export
#'
#' @examples
#' runGSEAPrerank("ranked_genes.rnk", "GSEA_Prerank_Rep1", "results/",
#'                "genesets.gmt")


runGSEAPrerank <- function(rnkFile, GSEArepLbl, 
                           GSEAresOutDir, gmtFile, pathToBat = NA, altDelim = NA,
                           ChipFile = NA, gmtAbbr = NA, scoring_scheme = "weighted",
                           collapse = "Collapse", modeProbe = "Max_probe", 
                           normMode = "meandiv", nperm = 1000, create_svgs = FALSE, 
                           include_only_symbols = TRUE, make_sets = TRUE,
                           zip_report = FALSE, numMarkers = 100, plot_top_x = 25,
                           rnd_seed = "timestamp", set_max = 500, set_min = 15,
                           showOutputOnConsole = TRUE) {
  
  filePaths <- list(rnkFile, GSEAresOutDir, gmtFile, ChipFile)
  names(filePaths) <- c("rnkFile", "GSEAresOutDir", "gmtFile", "ChipFile")
  
  checkFilePathsExist <- lapply(filePaths, file.exists)
  
  if(any(unlist(checkFilePathsExist) == FALSE)){
    whichCheck <- names(filePaths)[unlist(checkFilePathsExist) == FALSE]
    stop(paste0("One or more of the files supplied to runGSEAPrerank do not exist. Check: \n", 
                paste(whichCheck, collapse = "\n")))
  }
  
  checkFilePaths <- lapply(filePaths, rtocppath)
  
  logicalsList <- list(create_svgs, include_only_symbols, make_sets, zip_report)
  names(logicalsList) <- c("create_svgs", "include_only_symbols", "make_sets", "zip_report")
  
  if(all(unlist(logicalsList) %in% c(1L, 0L)) == FALSE){
    whichCheck <- names(logicalsList)[! unlist(logicalsList) %in% c(1L, 0L)]
    stop(paste0("Non-logical value(s) supplied when expecting logical value(s). Check: \n", 
                paste(whichCheck, collapse = "\n")))
  }
  
  checkLogicals <- lapply(logicalsList, rtocpbool)
  
  scoring_scheme <- match.arg(scoring_scheme, 
                              choices = c("weighted", "classic", "weighted_p2", "weighted_p1.5"))
  collapse <- match.arg(collapse, choices = c("Collapse", "No_Collapse", "Remap_Only"))
  
  modeProbe <- match.arg(modeProbe, choices = c("Max_probe", "Median_of_probe", "Mean_of_probe", 
                                                "Sum_of_probes", "Abs_max_of_probes"))
  normMode <- match.arg(normMode, choices = c("meandiv", "none"))
  rnd_seed <- match.arg(rnd_seed, choices = c("timestamp", "149"))
  
  commandRun = paste("gsea-cli.bat GSEAPreranked",
                     "-rnk", checkFilePaths["rnkFile"],
                     "-gmx", checkFilePaths["gmtFile"],
                     "-collapse", collapse,
                     "-chip", checkFilePaths["ChipFile"], 
                     "-mode", modeProbe,
                     "-norm", normMode,
                     "-nperm", as.character(nperm),
                     "-scoring_scheme", scoring_scheme,
                     "-rpt_label", GSEArepLbl,
                     "-create_svgs", checkLogicals["create_svgs"],
                     "-include_only_symbols",  checkLogicals["include_only_symbols"],
                     "-make_sets", checkLogicals["make_sets"],
                     "-zip_report", checkLogicals["zip_report"],
                     "-plot_top_x", as.character(plot_top_x),
                     "-rnd_seed", rnd_seed,
                     "-set_max", as.character(set_max),
                     "-set_min", as.character(set_min),
                     "-out", checkFilePaths["GSEAresOutDir"], collapse = " ")
  
  if(!is.na(altDelim)){
    if(nchar(altDelim) > 1){
      stop(paste0("altDelim is larger than one character"))
    }
    commandRun <- paste0(commandRun, " -altDelim ", altDelim) 
  }
  
  setwd(pathToBat)
  
  system(commandRun, intern = TRUE,
         ignore.stdout = FALSE, ignore.stderr = FALSE,
         wait = TRUE, input = NULL, show.output.on.console = showOutputOnConsole,
         minimized = FALSE, invisible = TRUE, timeout = 0)
  
  return(commandRun)
}
