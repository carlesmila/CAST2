#' Estimate the area of applicability
#'
#' @description
#' this function estimates the area of inapplicability index (AOII) and the derived
#' area of applicability (AOI) of spatial prediction models by
#' considering the distance of new data (i.e. a Raster Stack of spatial predictors
#' used in the models) in the predictor variable space to the data used for model
#' training. Predictors can be weighted in the ideal case based on the internal
#' variable importance of the machine learning algorithm used for model training.
#'
#' @param train a data.frame containing the data used for model training
#' @param predictors A RasterStack, RasterBrick or data.frame containing the data
#' the model was meant to make predictions for.
#' @param weight A data.frame containing weights for each variable
#' @param model A caret model used to extract weights from (based on variable importance)
#' @param variables character vector of predictor variables. if "all" then all variables
#' of the train dataset are used. Check varImp(model).
#' @param threshold Numeric or character indicating the quantile (e.g. "90%"). Used to
#' define the AOA
#' @param scale logical. If TRUE uncertainty is scaled between 0 and 1. See Details.
#' @param clstr Numeric or character. Spatial cluster affiliation for each data point. Should be used if replicates are present.
#' @param cl Cluster object created with parallel::makeCluster. To run things in parallel.
#' @details The "Area of Inapplicability Index" (AOII) and the corresponding Area of Applicability (AOA) are calculated.
#' Interpretation of results: If a location is very similar to the properties
#' of the training data it will have a low distance in the predictor variable space
#' (high inapplicability (aoii) index) while locations that are very different in their properties
#' will have a high Area of Inapplicability Index (AOII).
#' If scale is FALSE then AOII is returned as distance scaled by the average distance to a nearest training data point as
#' observed in the training data. The further the distance in this predicor space, the larger the AOII gets.
#' To get the AOA, a threshold to the AOII is applied.
#' @return A RasterStack or data.frame with the AOII and AOA
#' @author
#' Hanna Meyer
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(raster)
#' library(caret)
#'
#' # prepare sample data:
#' dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
#' dat <- aggregate(dat[,c("VW","Easting","Northing")],by=list(as.character(dat$SOURCEID)),mean)
#' pts <- st_as_sf(dat,coords=c("Easting","Northing"))
#' pts$ID <- 1:nrow(pts)
#' studyArea <- stack(system.file("extdata","predictors_2012-03-25.grd",package="CAST"))[[1:8]]
#' trainDat <- extract(studyArea,pts,df=TRUE)
#' trainDat <- merge(trainDat,pts,by.x="ID",by.y="ID")
#'
#' # visualize data spatially:
#' spplot(scale(studyArea))
#' plot(studyArea$DEM)
#' plot(pts[,1],add=TRUE,col="black")
#'
#' # first calculate uncertainty based on a set of variables with equal weights:
#' variables <- c("DEM","Easting","Northing")
#' AOA <- aoa(trainDat,studyArea,variables=variables)
#' spplot(AOA$AOAI, col.regions=viridis(100))+
#' spplot(AOA$AOA,col.regions=c("transparent","grey"))
#'
#' # or weight variables based on variable improtance from a (simple) trained model:
#' set.seed(100)
#' model <- train(trainDat[,which(names(trainDat)%in%variables)],
#' trainDat$VW,method="rf",importance=TRUE,tuneLength=1)
#' prediction <- predict(studyArea,model)
#' plot(varImp(model,scale=FALSE))
#' #
#' # note that coordinates are the major predictors here,
#' # so uncertainty becomes higher when moving away from the training data:
#' par(mfrow=c(1,2))
#' plot(prediction,main="predicted VW")
#' AOA <- aoa(trainDat,studyArea,model=model,variables=variables)
#' spplot(AOA$AOAI, col.regions=viridis(100))
#' spplot(AOA$AOAI, col.regions=viridis(100))+
#' spplot(AOA$AOA,col.regions=c("transparent","grey"))
#' }
#' @export aoa
#' @aliases aoa

aoa <- function (train, predictors, weight=NA, model=NA,
                 variables="all", threshold="90%",scale=FALSE,
                 clstr=NULL,cl=NULL){
  ### if not specified take all variables from train dataset as default:
  if(nrow(train)<=1){stop("at least two training points need to be specified")}
  if(length(variables)==1&&variables=="all"){
    variables=names(train)
  }
  #### Prepare output as either as RasterLayer or vector:
  out <- NA
  if (class(predictors)=="RasterStack"|class(predictors)=="RasterBrick"|
      class(predictors)=="RasterLayer"){
    out <- predictors[[1]]
    names(out) <- "uncertainty"
  }
  #### Extract weights from trained model:
  weight <- tryCatch(if(model$modelType=="Classification"){
    as.data.frame(t(apply(caret::varImp(model,scale=F)$importance,1,mean)))
  }else{
    as.data.frame(t(caret::varImp(model,scale=F)$importance[,"Overall"]))
  }, error=function(e) e)
  if(!inherits(weight, "error")){
    names(weight)<- rownames(caret::varImp(model,scale=F)$importance)
  }else{
    message("note: variables were not weighted either because no weights or model were given,
    or no variable importance could be retrieved from the given model.
    Check caret::varImp(model)")
  }
  #### order data:
  if (class(predictors)=="RasterStack"|class(predictors)=="RasterBrick"|
      class(predictors)=="RasterLayer"){
    predictors <- predictors[[na.omit(match(variables, names(predictors)))]]
  }else{
    predictors <- predictors[,na.omit(match(variables, names(predictors)))]
  }
  train <- train[,na.omit(match(variables, names(train)))]
  if(!inherits(weight, "error")){
    weight <- weight[,na.omit(match(variables, names(weight)))]
    if (any(weight<0)){
      weight[weight<0]<-0
      message("negative weights were set to 0")
    }
  }
  #### Scale data and weight predictors if applicable:
  train <- scale(train)
  scaleparam <- attributes(train)
  if(!inherits(weight, "error")){
    train <- sapply(1:ncol(train),function(x){train[,x]*unlist(weight[x])})
  }
  if (class(predictors)=="RasterStack"|class(predictors)=="RasterBrick"|
      class(predictors)=="RasterLayer"){
    predictors <- raster::as.data.frame(predictors)
  }
  predictors <- scale(predictors,center=scaleparam$`scaled:center`,#scaleparam$`scaled:center`
                      scale=scaleparam$`scaled:scale`)

  if(!inherits(weight, "error")){
    predictors <- sapply(1:ncol(predictors),function(x){predictors[,x]*unlist(weight[x])})
  }
  #### For each pixel caclculate distance to each training point and search for
  #### min distance:
  if(!is.null(cl)){ # if parallel then use parapply:
    parallel::clusterExport(cl=cl, varlist=c("train"))
    mindist <- parallel::parApply(cl=cl,predictors,1,FUN=function(x){
      tmp <- NA
      for (i in 1:nrow(train)){
        current_min <- dist(rbind(x,train[i,]))
        current_min <- pmin(current_min,tmp,na.rm=T)
        tmp <- current_min
      }
      return(current_min)
    })
  }else{ # ...if not in parallel loop over train data:
    tmp <- NA
    for (i in 1:nrow(train)){
      mindist <- apply(predictors,1,function(x){dist(rbind(x,train[i,]))})
      mindist <- pmin(mindist,tmp,na.rm=T)
      tmp <- mindist
    }
  }
  trainDist <- as.matrix(dist(train))
  diag(trainDist) <- NA
  # If data are highly clustered (repliates) make sure that distance to data from same
  # cluster are excluded
  if (!is.null(clstr)){
    for (i in 1:nrow(trainDist)){
      trainDist[i,clstr==clstr[i]] <- NA
    }
  }
  #scale the distance to nearest training point by average minimum distance as observed in training data
  trainDist_min <- apply(trainDist,1,FUN=function(x){min(x,na.rm=T)})
  trainDist_mean <- mean(trainDist_min)
  trainDist_quantiles <- quantile((trainDist_mean-trainDist_min)/trainDist_mean)
  mindist <- (mindist-trainDist_mean)/trainDist_mean
  # define threshold for AOA
  thres <- quantile(trainDist_mean-trainDist_min,
                    probs=as.numeric(gsub("%","",threshold))/100)/trainDist_mean

  #### scale distances if appliable:
  if (class(out)=="RasterLayer"){
    if(scale){
      raster::values(out) <- scales::rescale(mindist, to = c(0, 1))
    }else{
      raster::values(out) <- mindist
    }
  } else{
    if(scale){
      out <- scales::rescale(mindist, to = c(0, 1))
    }else{
      out <- mindist
    }
  }
  #### Create Mask for DOA and return statistics
  if (class(out)=="RasterLayer"){
    masked_result <- out
    raster::values(masked_result) <- 1
    masked_result[out>thres] <- 0
    masked_result <- raster::mask(masked_result,out)
    out <- raster::stack(out,masked_result)
    names(out) <- c("AOAI","AOA")
    aoa_stats <- list("mean"=trainDist_mean,
                      "quantiles"=t(data.frame(trainDist_quantiles)),
                      "selected_quantile" =thres)
    attributes(out)$aoa_stats <- NULL
    attributes(out)$aoa_stats <- aoa_stats
  }
  return(out)
}