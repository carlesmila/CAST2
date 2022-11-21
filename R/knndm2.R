#' EXPERIMENTAL2 K-fold Nearest Neighbour Distance Matching function
#' @description
#' This function implements the kNNDM algorithm and returns the necessary
#' indices to perform a kNNDM CV for map validation.
#'
#' @author Carles Milà
#' @param tpoints sf or sfc point object. Contains the training points samples.
#' @param modeldomain raster or sf object defining the prediction area (see Details).
#' @param ppoints sf or sfc point object. Contains the target prediction points. Optional. Alternative to modeldomain (see Details).
#' @param k integer. Final number of folds desired for CV. Defaults to 5.
#' @param samplesize numeric. How many points in the modeldomain should be sampled as prediction points?
#' Only required if modeldomain is used instead of ppoints.
#' @param sampling character. How to draw prediction points from the modeldomain? See spsample. Use sampling = "Fibonacci" for global applications.
#' Only required if modeldomain is used instead of ppoints.
#'
#' @return An object of class \emph{knndm} consisting of a list of six elements:
#' indx_train, indx_test (indices of the observations to use as
#' training/test data in each kNNDM CV iteration), Gij (distances for
#' G function construction between prediction and target points), Gj
#' (distances for G function construction during LOO CV), Gjstar (distances
#' for modified G function during kNNDM CV), clusters (list of cluster IDs).
#'
#' @details TBC.
#' @note Experimental cycle. Parallel implementation pending. Missing references.
#' @export
#'
#' @examples
#' ########################################################################
#' # Example 1: Simulated data - Weak clustering
#' ########################################################################
#'
#' library(sf)
#' library(ggplot2)
#'
#' # Simulate 5000 clustered training points in a 100x100 square
#' set.seed(1234)
#' simarea <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
#' simarea <- sf::st_polygon(simarea)
#' train_points <- clustered_sample(simarea, 5000, 200, 5)
#' pred_points <- sf::st_sample(simarea, 1000, type = "regular")
#'
#' # Run kNNDM for the whole domain, here the prediction points are known.
#' knndm_pred <- knndm2(train_points, ppoints = pred_points)
#' plot(knndm_pred)
#' ggplot() +
#'   geom_sf(data = simarea, alpha = 0) +
#'   geom_sf(data = train_points, col = knndm_pred$clusters)
#'
#' ########################################################################
#' # Example 2: Simulated data - Strong clustering
#' ########################################################################
#'
#' library(sf)
#'
#' # Simulate 5000 clustered training points in a 100x100 square
#' set.seed(1234)
#' simarea <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
#' simarea <- sf::st_polygon(simarea)
#' train_points <- clustered_sample(simarea, 5000, 50, 5)
#' pred_points <- sf::st_sample(simarea, 1000, type = "regular")
#'
#' # Run kNNDM for the whole domain, here the prediction points are known.
#' knndm_pred <- knndm2(train_points, ppoints = pred_points)
#' plot(knndm_pred)
#' ggplot() +
#'   geom_sf(data = simarea, alpha = 0) +
#'   geom_sf(data = train_points, col = knndm_pred$clusters)
#'
#' ########################################################################
#' # Example 3: Real- world example; using a modeldomain instead of previously
#' # sampled prediction locations
#' ########################################################################
#' \dontrun{
#' library(sf)
#' library(raster)
#'
#' ### prepare sample data:
#' dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
#' dat <- aggregate(dat[,c("DEM","TWI", "NDRE.M", "Easting", "Northing","VW")],
#'    by=list(as.character(dat$SOURCEID)),mean)
#' pts <- dat[,-1]
#' pts <- st_as_sf(pts,coords=c("Easting","Northing"))
#' st_crs(pts) <- 26911
#' studyArea <- raster::stack(system.file("extdata","predictors_2012-03-25.grd",package="CAST"))
#'
#' knndm_folds <- knndm2(pts, k=3, modeldomain=studyArea)
#' plot(knndm_folds)
#' ggplot() +
#'   geom_sf(data = pts, col = knndm_folds$clusters)
#'
#' #use for cross-validation:
#' library(caret)
#' ctrl <- trainControl(method="cv",
#'    index=knndm_folds$indx_train)
#' model_knndm <- train(dat[,c("DEM","TWI", "NDRE.M")],
#'    dat$VW,
#'    method="rf",
#'    trControl = ctrl)
#' model_knndm
#'}
knndm2 <- function(tpoints, modeldomain = NULL, ppoints = NULL,
                   k = 5, samplesize = 1000, sampling = "regular"){

  # Prior checks
  if (is.null(ppoints) & !is.null(modeldomain)) {

    message(paste0(samplesize, " prediction points are sampled from the modeldomain"))
    ppoints <- sampleFromArea(modeldomain, samplesize, type = "geo",
                              variables = NULL, sampling)
  }
  if (any(class(tpoints) %in% "sfc")) {
    tpoints <- sf::st_sf(geom = tpoints)
  }
  if (any(class(ppoints) %in% "sfc")) {
    ppoints <- sf::st_sf(geom = ppoints)
  }
  if(sf::st_crs(tpoints) != sf::st_crs(ppoints)){
    tpoints <- sf::st_transform(tpoints,sf::st_crs(ppoints))
    message("tpoints and ppoints must have the same CRS. tpoints have been transformed.")
  }

  # Gj: NND function for a cluster per point, i.e. LOO CV
  clust <- 1:nrow(tpoints)
  distmat <- sf::st_distance(tpoints)
  units(distmat) <- NULL
  Gj <- distclust(distmat, clust)

  # Gij: prediction to training NN distances
  Gij <- sf::st_distance(ppoints, tpoints)
  Gij <- apply(Gij, 1, min)
  attr(Gij,"units") <- NULL

  # Hierarchical clustering
  hc <- stats::hclust(d = stats::as.dist(distmat), method="average")

  # Build grid of number of clusters to try - we sample low numbers more intensively
  clustgrid <- data.frame(nk = as.integer(round(exp(c(seq(log(k), log(nrow(tpoints)),
                                                          length.out = 50))))))
  clustgrid$stat <- NA
  clustgrid <- clustgrid[!duplicated(clustgrid$nk),]
  clustgroups <- list()

  # We test each number of clusters
  for(nk in clustgrid$nk){

    # Cut nk clusters
    clust_nk <- stats::cutree(hc, k=nk)
    tabclust <- data.frame(clust_nk = unique(clust_nk))
    tabclust2 <- data.frame(clust_nk = clust_nk)
    tabclust$clust_k <- rep(1:k, ceiling(nk/k))[1:nk]
    tabclust2 <- merge(tabclust2, tabclust, by = "clust_nk")
    clust_k <- tabclust2$clust_k

    # Compute statistic
    Gjstar_i <- distclust(distmat, clust_k)
    clustgrid$stat[clustgrid$nk==nk] <- twosamples::wass_stat(Gjstar_i, Gij)
    clustgroups[[paste0("nk", nk)]] <- clust_k
    rm("clust_k", "clust_nk", "Gjstar_i", "tabclust", "tabclust2")
  }

  # Final configuration
  k_final <- clustgrid$nk[which.min(clustgrid$stat)]
  clust <- clustgroups[[paste0("nk", k_final)]]
  Gjstar <- distclust(distmat, clust)
  cfolds <- CAST::CreateSpacetimeFolds(data.frame(clust=clust), spacevar = "clust", k = k)
  res <- list(clusters = clust,
              indx_train = cfolds$index, indx_test = cfolds$indexOut,
              Gij = Gij, Gj = Gj, Gjstar = Gjstar)
  class(res) <- c("knndm", "list")
  res
}

# Helper function: Compute out-of-fold NN distance
distclust <- function(distm, folds){

  alldist <- c()
  for(f in unique(folds)){
    alldist <- c(alldist, apply(distm[f == folds, f != folds, drop=FALSE], 1, min))
  }
  alldist
}
