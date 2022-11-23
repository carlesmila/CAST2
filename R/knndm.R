#' EXPERIMENTAL K-fold Nearest Neighbour Distance Matching
#' @description
#' This function implements the kNNDM algorithm and returns the necessary
#' indices to perform a kNNDM CV for map validation.
#'
#' @author Carles Milà
#' @param tpoints sf or sfc point object. Contains the training points samples.
#' @param modeldomain raster or sf object defining the prediction area (see Details).
#' @param ppoints sf or sfc point object. Contains the target prediction points. Optional. Alternative to modeldomain (see Details).
#' @param k integer. Final number of folds desired for CV. Defaults to 5.
#' @param maxp numeric. Maximum fold size allowed, defaults to 1/k*2. Partitions leading to larger folds than `maxp` are discarded.
#' Increase if a better fit is required at the cost of having less balanced folds.
#' @param linkf character. Link function for agglomerative clustering. Defaults to "ward.D2" for balanced clusters.
#' Check `stats::hclust` for other options.
#' @param samplesize numeric. How many points in the modeldomain should be sampled as prediction points?
#' Only required if modeldomain is used instead of ppoints.
#' @param sampling character. How to draw prediction points from the modeldomain? See spsample. Use "Fibonacci" for global applications.
#' Only required if modeldomain is used instead of ppoints.
#'
#' @return An object of class \emph{knndm} consisting of a list of seven elements:
#' indx_train, indx_test (indices of the observations to use as
#' training/test data in each kNNDM CV iteration), Gij (distances for
#' G function construction between prediction and target points), Gj
#' (distances for G function construction during LOO CV), Gjstar (distances
#' for modified G function during kNNDM CV), clusters (list of cluster IDs), and
#' stat (Wasserstein distance statistic).
#'
#' @details TBC.
#' @note Experimental cycle. Documentation to be completed. Missing references. Large distance matrices?
#' @export
#'
#' @examples
#' ########################################################################
#' # Example 1: Simulated data - Randomly-distributed training points
#' ########################################################################
#'
#' library(sf)
#' library(ggplot2)
#'
#' # Simulate 1000 random training points in a 100x100 square
#' set.seed(1234)
#' simarea <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
#' simarea <- sf::st_polygon(simarea)
#' train_points <- sf::st_sample(simarea, 1000, type = "random")
#' pred_points <- sf::st_sample(simarea, 1000, type = "regular")
#'
#' # Run kNNDM for the whole domain, here the prediction points are known.
#' knndm_folds <- knndm(train_points, ppoints = pred_points)
#' knndm_folds$stat # Integral between the two curves
#' plot(knndm_folds)
#' ggplot() +
#'   geom_sf(data = simarea, alpha = 0) +
#'   geom_sf(data = train_points, col = knndm_folds$clusters)
#'
#' ########################################################################
#' # Example 2: Simulated data - Clustered training points
#' ########################################################################
#'
#' library(sf)
#' library(ggplot2)
#'
#' # Simulate 1000 clustered training points in a 100x100 square
#' set.seed(1234)
#' simarea <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
#' simarea <- sf::st_polygon(simarea)
#' train_points <- clustered_sample(simarea, 1000, 50, 5)
#' pred_points <- sf::st_sample(simarea, 1000, type = "regular")
#'
#' # Run kNNDM for the whole domain, here the prediction points are known.
#' knndm_folds <- knndm(train_points, ppoints = pred_points)
#' knndm_folds$stat # Integral between the two curves
#' plot(knndm_folds)
#' ggplot() +
#'   geom_sf(data = simarea, alpha = 0) +
#'   geom_sf(data = train_points, col = knndm_folds$clusters)
#'
#' ########################################################################
#' # Example 3: Real- world example; using a modeldomain instead of previously
#' # sampled prediction locations
#' ########################################################################
#' \dontrun{
#' library(sf)
#' library(raster)
#' library(ggplot2)
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
#' knndm_folds <- knndm(pts, modeldomain=studyArea)
#' knndm_folds$stat # Integral between the two curves
#' plot(knndm_folds)
#' plot_geodist(pts, studyArea, cvfolds = knndm_folds$indx_test)
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
knndm <- function(tpoints, modeldomain = NULL, ppoints = NULL,
                  k = 5, maxp = NULL, linkf = "ward.D2",
                  samplesize = 1000, sampling = "regular"){

  # Prior checks
  if (is.null(ppoints) & !is.null(modeldomain)) {

    message(paste0(samplesize, " prediction points are sampled from the modeldomain"))
    ppoints <- sampleFromArea(modeldomain, samplesize, type = "geo", variables = NULL, sampling)
  }
  if (any(class(tpoints) %in% "sfc")) {
    tpoints <- sf::st_sf(geom = tpoints)
  }
  if (any(class(ppoints) %in% "sfc")) {
    ppoints <- sf::st_sf(geom = ppoints)
  }
  if(sf::st_crs(tpoints) != sf::st_crs(ppoints)){
    tpoints <- sf::st_transform(tpoints, sf::st_crs(ppoints))
    message("tpoints and ppoints must have the same CRS. tpoints have been transformed.")
  }
  if (is.null(maxp)) {
    maxp <- 1/k*2
  }

  # Gj: NND function for a cluster per point, i.e. LOO CV
  clust <- 1:nrow(tpoints)
  distmat <- sf::st_distance(tpoints)
  units(distmat) <- NULL
  Gj <- distclust(distmat, clust)

  # Gij: prediction to training NN distances
  Gij <- sf::st_distance(ppoints, tpoints)
  Gij <- apply(Gij, 1, min)
  units(Gij) <- NULL

  # Hierarchical clustering
  hc <- stats::hclust(d = stats::as.dist(distmat), method=linkf)

  # Build grid of number of clusters to try - we sample low numbers more intensively
  clustgrid <- data.frame(nk = as.integer(round(exp(seq(log(k), log(nrow(tpoints)),
                                                          length.out = 50)))))
  clustgrid$stat <- NA
  clustgrid <- clustgrid[!duplicated(clustgrid$nk),]
  clustgroups <- list()

  # We test each number of clusters
  for(nk in clustgrid$nk){

    # Cut nk clusters
    clust_nk <- stats::cutree(hc, k=nk)
    tabclust <- as.data.frame(table(clust_nk))
    tabclust <- tabclust[order(tabclust$Freq, decreasing=T),]
    tabclust$clust_k <- NA

    # We don't merge big clusters
    clust_i <- 1
    for(i in 1:nrow(tabclust)){
      if(tabclust$Freq[i] >= nrow(tpoints)/k){
        tabclust$clust_k[i] <- clust_i
        clust_i <- clust_i + 1
      }
    }
    rm("clust_i")
    clust_i <- setdiff(1:k, unique(tabclust$clust_k))
    tabclust$clust_k[is.na(tabclust$clust_k)] <- rep(clust_i, ceiling(nk/length(clust_i)))[1:sum(is.na(tabclust$clust_k))]
    tabclust2 <- data.frame(ID = 1:length(clust_nk), clust_nk = clust_nk)
    tabclust2 <- merge(tabclust2, tabclust, by = "clust_nk")
    tabclust2 <- tabclust2[order(tabclust2$ID),]
    clust_k <- tabclust2$clust_k

    # Compute statistic if not exceeding limit
    if(!any(table(clust_k)/length(clust_k)>maxp)){
      Gjstar_i <- distclust(distmat, clust_k)
      clustgrid$stat[clustgrid$nk==nk] <- twosamples::wass_stat(Gjstar_i, Gij)
      clustgroups[[paste0("nk", nk)]] <- clust_k
    }
  }

  # Final configuration
  k_final <- clustgrid$nk[which.min(clustgrid$stat)]
  stat_final <- min(clustgrid$stat, na.rm=T)
  clust <- clustgroups[[paste0("nk", k_final)]]
  Gjstar <- distclust(distmat, clust)
  cfolds <- CAST::CreateSpacetimeFolds(data.frame(clust=clust), spacevar = "clust", k = k)

  # Output
  res <- list(clusters = clust,
              indx_train = cfolds$index, indx_test = cfolds$indexOut,
              Gij = Gij, Gj = Gj, Gjstar = Gjstar, stat = stat_final)
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
