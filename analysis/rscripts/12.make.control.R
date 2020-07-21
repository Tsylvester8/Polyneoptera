make.control <- function(mainType = NULL,
                         dataFile = NULL,
                         treeFile = NULL,
                         outDir = NULL,
                         maxChrNum = NULL,
                         minChrNum = NULL,
                         branchMul = NULL,
                         simulationsNum = NULL,
                         logFile = NULL,
                         optimizePointsNum = NULL, 
                         optimizeIterNum = NULL,
                         pars = list(),
                         control.path = NULL){
  control <- c()
  
  if(is.null(mainType)){
    stop("specify the model")
  }
  if(!is.element(el = mainType,
                 set = c("All_Models",
                         "Run_Fix_Param",
                         "Optimize_Model"))){
    stop("main type should be either 'All_Models', 'Run_Fix_Param' or 'Optimize_Model'")
  }
  if(is.null(dataFile)){
    stop("provide the path to the chromosome counts")
  }
  if(is.null(treeFile)){
    stop("provide the path to the tree file")
  }
  if(is.null(outDir)){
    stop("provide the path to tree file")
  }
  if(mainType == "Optimize_Model"){
    if(length(pars) == 0){
      stop("model parameters are not given")
    }
  }
  if(mainType == "Run_Fix_Param"){
    if(length(pars) == 0){
      stop("model parameters are not given")
    }
  }
  control <- c(paste("_mainType", mainType, sep = " "),
               paste("_dataFile", dataFile, sep = " "),
               paste("_treeFile", treeFile, sep = " "),
               paste("_outDir", outDir, sep = " "))
  
  if(!is.null(maxChrNum)){
    if(!is.numeric(maxChrNum)){
      stop("Chromosome number should be numeric")
    }else{
      control <- c(control,
                   paste("_maxChrNum", maxChrNum))
    }
  }
  if(!is.null(minChrNum)){
    if(!is.numeric(minChrNum)){
      stop("Chromosome number should be numeric")
    }else{
      control <- c(control,
                   paste("_minChrNum", minChrNum))
    }
  }
  if(!is.null(branchMul)){
    control <- c(control,
                 paste("_branchMul", branchMul, sep = " "))
  }
  if(!is.null(simulationsNum)){
    control <- c(control,
                 paste("_simulationsNum", simulationsNum, sep = " "))
  }
  if(!is.null(logFile)){
    control <- c(control,
                 paste("_branchMul", logFile, sep = " "))
  }
  if(!is.null(optimizePointsNum)){
    if(is.null(optimizeIterNum)){
      stop("optimizePointsNum and optimizeIterNum should be given together")
    }
  }
  if(!is.null(optimizeIterNum)){
    if(is.null(optimizePointsNum)){
      stop("optimizePointsNum and optimizeIterNum should be given together")
    }
  }
  if(!is.null(optimizePointsNum) && !is.null(optimizeIterNum)){
    control <- c(control,
                 paste("_optimizePointsNum", optimizePointsNum, sep = " "),
                 paste("_optimizeIterNum", optimizeIterNum, sep = " "))
  }
  if(mainType %in% c("Optimize_Model","Run_Fix_Param")){
    if(length(pars) == 0){
      stop("model parameters are not given")
    }
    if(length(pars) != 0){
      pars.vector <- c()
      for (i in 1:(length(pars))){
        pars.vector[(i)] <- paste("_", names(pars)[i], " ", pars[i][[1]], sep = "")
      }
      control <- c(control, pars.vector)
    }
  }
  if(is.null(control.path)){
    stop("path to save the control file is not given")
  }
  write(control, file= control.path)
}

