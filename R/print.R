#' Print TrainDI
#' @param x trainDI object
#' @param ... other params
#' @export

print.trainDI = function(x, ...){
  cat(paste0("DI of ", nrow(x$train), " observation \n"))
  cat(paste0("Predictors:"), x$variables, "\n\n")

  cat("AOA Threshold: ")
  cat(x$threshold)

}


#' Show TrainDI
#' @param x trainDI object
#' @param ... other params
#' @export

show.trainDI = function(x, ...){
  print.trainDI(x)
}




#' Print AOA
#' @param x aoa object
#' @param ... other params
#' @export


print.aoa = function(x, ...){

  cat("AOA:\n")

  print(x$AOA)

  cat("\n\nPredictor Weights:\n")

  print(x$parameters$weight)

  cat("\n\nAOA Threshold: ")
  cat(x$parameters$threshold)



}


#' Show AOA
#' @param x aoa object
#' @param ... other params
#' @export


show.aoa = function(x, ...){
  print.aoa(x)
}
