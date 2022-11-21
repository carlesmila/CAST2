#' EXPERIMENTAL K-fold Nearest Neighbour Distance Matching function
#' @description
#' This function implements the kNNDM algorithm and returns the necessary
#' indices to perform a kNNDM CV for map validation.
#'
#' @author Carles Milà
#' @param tpoints sf or sfc point object. Contains the training points samples.
#' @param modeldomain raster or sf object defining the prediction area (see Details).
#' @param ppoints sf or sfc point object. Contains the target prediction points. Optional. Alternative to modeldomain (see Details).
#' @param k integer. Final number of folds desired for CV. Defaults to 10.
#' @param kclust integer. Intermediate number of clusters for agglomerative clustering. Lower numbers result in faster results,
#' but the quality of the results may be compromised. Defaults to half of the training points.
#' Consider increasing that number if G*j(r) < Gij(r) for short distances, which may happen in cases with weak clustering.
#' @param maxsize numeric. maximum proportion of training samples to have in a fold. Defaults to 1/k*2
#' @param samplesize numeric. How many points in the modeldomain should be sampled as prediction points?
#' Only required if modeldomain is used instead of ppoints.
#' @param sampling character. How to draw prediction points from the modeldomain? See spsample. Use sampling = "Fibonacci" for global applications.
#'  Only required if modeldomain is used instead of ppoints.
#' @param cl TBC.
#' @param verbose logical. Show progress, TRUE by default.
#'
#' @return An object of class \emph{knndm} consisting of a list of six elements:
#' indx_train, indx_test (indices of the observations to use as
#' training/test data in each kNNDM CV iteration), Gij (distances for
#' G function construction between prediction and target points), Gj
#' (distances for G function construction during LOO CV), Gjstar (distances
#' for modified G function during kNNDM CV), clusters (ordered list of cluster
#' IDs for each training point).
#'
#' @details TBC.
#' @note Experimental cycle. Parallel implementation pending. Missing references.
#' @export
#'
#' @examples
#' ########################################################################
#' # Example 1: Simulated data
#' ########################################################################
#'
#' library(sf)
#' library(ggplot2)
#'
#' # Simulate 10000 clustered training points in a 100x100 square
#' set.seed(1234)
#' simarea <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
#' simarea <- sf::st_polygon(simarea)
#' train_points <- clustered_sample(simarea, 500, 50, 5)
#' pred_points <- sf::st_sample(simarea, 1000, type = "regular")
#'
#' # Run kNNDM for the whole domain, here the prediction points are known.
#' knndm_pred <- knndm(train_points, ppoints = pred_points)
#' plot(knndm_pred)
#' ggplot() +
#'   geom_sf(data = simarea, alpha = 0) +
#'   geom_sf(data = train_points, col = knndm_pred$clusters)
#'
#' ########################################################################
#' # Example 2: Real- world example; using a modeldomain instead of previously
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
#' knndm_folds <- knndm(pts, k=5, modeldomain=studyArea)
#' plot(knndm_folds)
#' ggplot() +
#'   geom_sf(data = knndm_folds, col = knndm_folds$clusters)
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
knndm <- function(tpoints, modeldomain = NULL, ppoints = NULL,
                  k = 10, kclust = NULL, maxsize = NULL,
                  samplesize = 1000, sampling = "regular",
                  cl = NULL, verbose = T){

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
  if (is.null(maxsize)) {
    maxsize <- 1/k*2
  }
  if (is.null(kclust)) {
    kclust <- 50
  }
  if(sf::st_crs(tpoints) != sf::st_crs(ppoints)){
    tpoints <- sf::st_transform(tpoints,sf::st_crs(ppoints))
    message("tpoints and ppoints must have the same CRS. tpoints have been transformed.")
  }
  if (!is.null(cl)) {
    parallel::clusterExport(cl=cl, "distclust")
  }

  # Gj: NND function for a cluster per point, i.e. LOO CV
  clust <- 1:nrow(tpoints)
  distmat <- sf::st_distance(tpoints)
  units(distmat) <- NULL
  # Gj <- distclust(distmat, clust)
  diag(distmat) <- NA
  system.time(Gj <- apply(distmat, 1, min, na.rm = T))

  # Gij: prediction to training NN distances
  Gij <- sf::st_distance(ppoints, tpoints)
  Gij <- apply(Gij, 1, min)
  attr(Gij,"units") <- NULL

  # Check if there's is clustering in the data in the first place
  testks <- stats::ks.test(Gj, Gij, alternative = "great")
  if(testks$p.value >= 0.05){
    clust <- sample(rep(1:k, ceiling(nrow(tpoints)/k)), size = nrow(tpoints), replace=F)
    warning("No evidence of clustering has been found, a random CV assignment is returned")

  }else{

    # We cluster points as a starting point
    hc <- stats::hclust(d = stats::as.dist(distmat), method="average")
    clust <- stats::cutree(hc, k=kclust)

    # Helper function: Compute statistic for each candidate merger
    assess_merge <- function(i, listpos, clust, distmat, Gij){
      clust_i <- clust
      clust_i[clust_i ==listpos[i,2]] <- listpos[i,1]
      Gjstar_i <- distclust(distmat, clust_i)
      listpos$stat[i] <- twosamples::wass_stat(Gjstar_i, Gij)
    }

    # We begin the loop until we have the target number of clusters
    while(length(unique(clust))>k){

      # List possible merging combinations, discard too big cluster merging and try random 100
      listpos <- expand.grid(V1 = unique(clust), V2 = unique(clust))
      listpos <- listpos[listpos$V1 != listpos$V2,]
      listpos$stat <- NA
      listpos$Freq <- sapply(listpos$V1, function(x) sum(clust==x)) +
        sapply(listpos$V2, function(x) sum(clust==x))
      listpos <- listpos[listpos$Freq/length(clust) <= maxsize,]
      if(nrow(listpos)>200){
        listpos <- listpos[sample(1:nrow(listpos), 200),]
      }

      # Compute statistic for each candidate merger
      if (!is.null(cl)){
        listpos$stat <- parallel::parSapply(cl = cl, X = 1:nrow(listpos), FUN = assess_merge,
                                            listpos = listpos, clust = clust,
                                            distmat = distmat, Gij = Gij)
      }else{
        listpos$stat <- sapply(1:nrow(listpos), assess_merge, listpos = listpos,
                               clust = clust, distmat = distmat, Gij = Gij)
      }

      # Choose the one that minimizes the statistic and has the lowest number of points
      listpos <- listpos[order(listpos$Freq),]
      clust[clust == listpos$V2[which.min(listpos$stat)]] <- listpos$V1[which.min(listpos$stat)]
      if(isTRUE(verbose)){
        print(paste0("Current number of clusters: ", length(unique(clust))))
      }
    }
  }

  # Final configuration
  clust <- as.integer(as.factor(clust))
  Gjstar <- distclust(distmat, clust)
  cfolds <- CAST::CreateSpacetimeFolds(data.frame(clust=clust), spacevar = "clust", k = k)
  res <- list(clusters = clust,
              indx_train = cfolds$index, indx_test = cfolds$indexOut,
              Gij = Gij, Gj = Gj, Gjstar = Gjstar)
  class(res) <- c("knndm", "list")
  res
}

# Helper function
distclust <- function(distm, folds){

  alldist <- c()
  for(f in unique(folds)){
    alldist <- c(alldist, apply(distm[f == folds, f != folds, drop=FALSE], 1, min))
  }
  alldist
}

