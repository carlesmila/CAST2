#' Plot euclidean nearest neighbor distances in geographic space or feature space
#'
#' @description Density plot of nearest neighbor distances in geographic space or feature space between training data as well as between training data and prediction locations.
#' Optional, the nearest neighbor distances between training data and test data or between training data and CV iterations is shown.
#' The plot can be used to check the suitability of a chosen CV method to be representative to estimate map accuracy. Alternatively distances can also be calculated in the multivariate feature space.
#'
#' @param x object of class sf, training data locations
#' @param modeldomain raster or sf object defining the prediction area (see Details)
#' @param cvfolds optional. List of row indices of x that are held back in each CV iteration. See e.g. ?createFolds or ?createSpaceTimeFolds
#' @param testdata optional. object of class sf: Data used for independent validation
#' @param samplesize numeric. How many random prediction samples should be used? Only required if modeldomain is a raster (see Details)
#' @param scale logical. Present distances on log scale?
#' @param distance "geo" or "feature". Should the distance be computed in geographic space or in the normalized multivariate predictor space (see Details)
#' @param variables character vector defining the predictor variables used if distance="feature. If not provided all variables included in modeldomain are used.
#' @return A plot and a data.frame containing the distances
#'
#'
#'
#' @details The modeldomain is a sf polygon or a raster that defines the prediction area. The function takes a regular point sample (amount defined by samplesize) from the spatial extent.
#'     If distance = "feature", the argument modeldomain (and if provided then also the testdata) has to include predictors. Predictor values for x are optional if modeldomain is a raster. If not provided they are extracted from the modeldomain rasterStack.
#'
#' @import ggplot2
#'
#' @author Hanna Meyer, Edzer Pebesma, Marvin Ludwig
#' @examples
#' \dontrun{
#' library(sf)
#' library(raster)
#' library(caret)
#'
#' ########### prepare sample data:
#' dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
#' dat <- aggregate(dat[,c("DEM","TWI", "NDRE.M", "Easting", "Northing")],
#' by=list(as.character(dat$SOURCEID)),mean)
#' pts <- dat[,-1]
#' pts <- st_as_sf(pts,coords=c("Easting","Northing"))
#' st_crs(pts) <- 26911
#' pts_train <- pts[1:29,]
#' pts_test <- pts[30:42,]
#' studyArea <- raster::stack(system.file("extdata","predictors_2012-03-25.grd",package="CAST"))
#' studyArea = studyArea[[c("DEM","TWI", "NDRE.M", "NDRE.Sd", "Bt")]]
#'
#' ########### Distance between training data and new data:
#' dist <- plot_geodist(pts_train,studyArea)
#'
#' ########### Distance between training data, new data and test data:
#' mapview(pts_train,col.regions="blue")+mapview(pts_test,col.regions="red")
#' dist <- plot_geodist(pts_train,studyArea,testdata=pts_test)
#'
#' ########### Distance between training data, new data and CV folds:
#' folds <- createFolds(1:nrow(pts_train),k=3,returnTrain=FALSE)
#' plot_geodist(x=pts_train, modeldomain=studyArea, cvfolds=folds)
#'
#' ########### Distances in the feature space:
#' plot_geodist(x=pts_train, modeldomain=studyArea,
#' distance = "feature",variables=c("DEM","TWI", "NDRE.M"))
#'
#' plot_geodist(x=pts_train, modeldomain=studyArea, cvfolds = folds, testdata = pts_test,
#' distance = "feature",variables=c("DEM","TWI", "NDRE.M"))
#' }
#' @export

plot_geodist <- function(x,
                         modeldomain,
                         samplesize=2000,
                         cvfolds=NULL,
                         testdata=NULL,
                         scale=TRUE,
                         distance = "geo",
                         variables=NULL){


  # input formatting ------------
  x <- sf::st_transform(x,4326) # not sure if thus is needed

  if(distance == "feature"){
    if(is.null(variables)){
      variables <- names(modeldomain)
    }
    if(!any(names(x)%in%variables)){
      x <- sf::st_as_sf(raster::extract(modeldomain, x, df = TRUE, sp = TRUE))
    }
    if(!is.null(testdata)){
      if(!any(names(testdata)%in%variables)){
        testdata <- sf::st_as_sf(raster::extract(modeldomain, testdata, df = TRUE, sp = TRUE))
      }
    }
  }

  if(inherits(modeldomain, "Raster")){
    modeldomain <- sampleFromArea(modeldomain, samplesize, distance,variables)
  }


  # required steps ----
  # always do sample-to-sample and sample-to-prediction
  s2s <- sample2sample(x, distance,variables)
  s2p <- sample2prediction(x, modeldomain, distance, samplesize,variables)

  dists <- rbind(s2s, s2p)

  # optional steps ----
  ##### Distance to test data:
  if(!is.null(testdata)){
    s2t <- sample2test(x, testdata, distance,variables)
    dists <- rbind(dists, s2t)
  }

  ##### Distance to CV data:
  if(!is.null(cvfolds)){
    cvd <- cvdistance(x, cvfolds, distance,variables)
    dists <- rbind(dists, cvd)
  }

  # Compile and plot data ----
  xlabs <- "geographic distance (m)"
  if(distance=="feature"){ xlab <- "feature space distance"}
  what <- "" #just to avoid check note
  if (distance=="feature"){unit ="unitless"}
  p <- ggplot2::ggplot(data=dists, aes(x=dist, group=what, fill=what)) +
    ggplot2::geom_density(adjust=1.5, alpha=.4) +
    ggplot2::scale_fill_discrete(name = "distance function") +
    ggplot2::xlab(xlabs) +
    ggplot2::theme(legend.position="bottom")

  if(scale){
    p <- p+scale_x_log10(labels=round)
  }
  print(p)
  return(dists)
}




# Sample to Sample Distance

sample2sample <- function(x, distance,variables){

  if(distance == "geo"){
    sf::sf_use_s2(TRUE)
    d <- sf::st_distance(x)
    diag(d) <- Inf
    min_d <- apply(d, 1, min)
    sampletosample <- data.frame(dist = min_d,
                                 what = factor("sample-to-sample"),
                                 type = "geo")


  }else if(distance == "feature"){
    x <- x[,variables]
    x <- sf::st_drop_geometry(x)
    scaleparam <- attributes(scale(x))
    x <- data.frame(scale(x))

    # sample to sample feature distance
    d <- c()
    for (i in 1:nrow(x)){
      trainDist <- FNN::knnx.dist(x[i,],x,k=1)
      trainDist[i] <- NA
      d <- c(d,min(trainDist,na.rm=T))
    }
    sampletosample <- data.frame(dist = d,
                                 what = factor("sample-to-sample"),
                                 type = "feature")

  }
  return(sampletosample)
}


# Sample to Prediction


sample2prediction = function(x, modeldomain, distance, samplesize,variables){

  if(distance == "geo"){
    modeldomain <- sf::st_transform(modeldomain, sf::st_crs(x))
    sf::sf_use_s2(TRUE)
    d0 <- sf::st_distance(modeldomain, x)
    min_d0 <- apply(d0, 1, min)
    sampletoprediction <- data.frame(dist = min_d0,
                                     what = factor("sample-to-prediction"),
                                     type = "geo")

  }else if(distance == "feature"){
    x <- x[,variables]
    x <- sf::st_drop_geometry(x)
    scaleparam <- attributes(scale(x))
    x <- data.frame(scale(x))

    modeldomain <- modeldomain[,variables]
    modeldomain <- sf::st_drop_geometry(modeldomain)
    modeldomain <- data.frame(scale(modeldomain,center=scaleparam$`scaled:center`,
                                    scale=scaleparam$`scaled:scale`))

    target_dist_feature <- c()
    for (i in 1:nrow(modeldomain)){
      trainDist <- FNN::knnx.dist(modeldomain[i,],x,k=1)
      target_dist_feature <- c(target_dist_feature,min(trainDist,na.rm=T))
    }
    sampletoprediction <- data.frame(dist = target_dist_feature,
                                     what = "sample-to-prediction",
                                     type = "feature")
  }

  return(sampletoprediction)
}


# sample to test


sample2test <- function(x, testdata, distance,variables){

  if(distance == "geo"){
    testdata <- sf::st_transform(testdata,4326)
    d_test <- sf::st_distance(testdata, x)
    min_d_test <- apply(d_test, 1, min)

    dists_test <- data.frame(dist = min_d_test,
                             what = factor("sample-to-test"),
                             type = "geo")


  }else if(distance == "feature"){
    x <- x[,variables]
    x <- sf::st_drop_geometry(x)
    scaleparam <- attributes(scale(x))
    x <- data.frame(scale(x))

    testdata <- testdata[,variables]
    testdata <- sf::st_drop_geometry(testdata)
    testdata <- data.frame(scale(testdata,center=scaleparam$`scaled:center`,
                                 scale=scaleparam$`scaled:scale`))


    test_dist_feature <- c()
    for (i in 1:nrow(testdata)){
      testDist <- FNN::knnx.dist(testdata[i,],x,k=1)
      test_dist_feature <- c(test_dist_feature,min(testDist,na.rm=T))
    }
    dists_test <- data.frame(dist = test_dist_feature,
                             what = "sample-to-test",
                             type = "feature")


  }
  return(dists_test)
}



# between folds

cvdistance <- function(x, cvfolds, distance,variables){

  if(distance == "geo"){
    d_cv <- c()
    for (i in 1:length(cvfolds)){
      d_cv_tmp <- sf::st_distance(x[cvfolds[[i]],], x[-cvfolds[[i]],])
      d_cv <- c(d_cv,apply(d_cv_tmp, 1, min))
    }

    dists_cv <- data.frame(dist = d_cv,
                           what = factor("CV-distances"),
                           type = "geo")


  }else if(distance == "feature"){
    x <- x[,variables]
    x <- sf::st_drop_geometry(x)
    x <- data.frame(scale(x))

    d_cv <- c()
    for(i in 1:length(cvfolds)){
      testdata_i <- x[cvfolds[[i]],]
      traindata_i <- x[-cvfolds[[i]],]

      for (k in 1:nrow(testdata_i)){
        trainDist <- FNN::knnx.dist(testdata_i[k,],traindata_i,k=1)
        trainDist[k] <- NA
        d_cv <- c(d_cv,min(trainDist,na.rm=T))
      }
    }


    dists_cv <- data.frame(dist = d_cv,
                           what = factor("CV-distances"),
                           type = "feature")

  }
  return(dists_cv)
}





sampleFromArea <- function(modeldomain, samplesize, distance,variables){

  ##### Distance to prediction locations:
  # regularly spread points (prediction locations):
  # see https://edzer.github.io/OGH21/
  if(inherits(modeldomain, "Raster")){
    if(samplesize>raster::ncell(modeldomain)){
      samplesize <- raster::ncell(modeldomain)
      message(paste0("samplesize for new data shouldn't be larger than number of pixels.
              Samplesize was reduced to ",raster::ncell(modeldomain)))
    }
    modeldomainextent <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(modeldomain)))
  }else{
    modeldomainextent <- modeldomain
  }

  sf::sf_use_s2(FALSE)
  sf::st_as_sf(modeldomainextent) |>
    sf::st_transform(4326) -> bb
  methods::as(bb, "Spatial") |>
    sp::spsample(n =samplesize, type = "regular")  |>
    sf::st_as_sfc() |>
    sf::st_set_crs(4326) -> predictionloc

  predictionloc <- sf::st_as_sf(predictionloc)

  if(distance == "feature"){

    if(!any(names(x)%in%variables)){
      x <- sf::st_as_sf(raster::extract(modeldomain, x, df = TRUE, sp = TRUE))
    }
    predictionloc <- sf::st_as_sf(raster::extract(modeldomain, predictionloc, df = TRUE, sp = TRUE))
    predictionloc <- na.omit(predictionloc)
  }

  return(predictionloc)

}
