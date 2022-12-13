#' NN distance based block CV
#' @description
#' This function estimates the optimal cross-validation folds.
#' Therefore the function calculates the nearest neighbor distance
#' between training locations and prediction locations (Gij), as well as
#' between different cross-validation folds  (Gjstar). The Wasserstein-Test is used
#' to select the CV-folds yielding the lowest distance between Gjstar and Gij.
#' @param x A sf object containing the training points.
#' @param modeldomain A raster or terra object containing the predictor variables.
#' @param samplesize Optional. The number of samplepoints generated in the modeldomain to describe the prediction space.
#' If NULL (default), the samplesize equals the number of cells in the predictor raster divided by ten.
#' @param sampling The type in which the samplepoints are generated. The default value is "Fibonacci", which is recommended for global studies.
#' @param n_steps The number of steps until termination.
#' @param maxn The maximal number of clusters to be fitted. By default 10.
#' @param cl A cluster object.
#' @details The nearest neighbor distance distribution between training and prediction locations,
#' as well as between cross-validation folds are calculated. Based on that, the
#' cross-validation split yielding the lowest distance between those distributions is chosen.
#' @details TBC.
#' @note Experimental.
#' @export
#'
#' @examples
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
#' bnndm_folds <- bnndm(pts, modeldomain=studyArea)
#' bnndm_folds$stat # Integral between the two curves
#' plot(bnndm_folds)
#' plot_geodist(pts, studyArea, cvfolds = bnndm_folds$indx_test)
#' ggplot() +
#'   geom_sf(data = pts, col = bnndm_folds$clusters)
#'
#' #use for cross-validation:
#' library(caret)
#' ctrl <- trainControl(method="cv",
#'    index=bnndm_folds$indx_train)
#' model_bnndm <- train(dat[,c("DEM","TWI", "NDRE.M")],
#'    dat$VW,
#'    method="rf",
#'    trControl = ctrl)
#' model_bnndm
#'}


bnndm <- function(x=NULL,
                     modeldomain=NULL,
                     samplesize=1000,
                     sampling = "regular",
                     n_steps = 10,
                     maxn = 10,
                     cl=NULL) {


  # helper functions ---------------------------------------------------

  # Sample prediction location from the study area:
  sampleFromArea <- function(modeldomain, samplesize, sampling){

    ##### Distance to prediction locations:
    # regularly spread points (prediction locations):
    # see https://edzer.github.io/OGH21/
    if(samplesize>terra::ncell(modeldomain)){
      samplesize <- terra::ncell(modeldomain)
      message(paste0("samplesize for new data shouldn't be larger than number of pixels.
              Samplesize was reduced to ",terra::ncell(modeldomain)))
    }


    sf::sf_use_s2(FALSE)

    bb <- terra::as.polygons(modeldomain) |>
      sf::st_as_sf() |>
      sf::st_geometry() |>
      sf::st_union() |>
      sf::st_transform(4326)

    methods::as(bb, "Spatial") |>
      sp::spsample(n = samplesize, type = sampling)  |>
      sf::st_as_sfc() |>
      sf::st_set_crs(4326) -> predictionloc

    predictionloc <- sf::st_as_sf(predictionloc)

    return(predictionloc)

  }

  # Compute out-of-fold NN distance
  distclust <- function(distm, folds){

    alldist <- c()
    for(f in unique(folds)){
      alldist <- c(alldist, apply(distm[f == folds, f != folds, drop=FALSE], 1, min))
    }
    alldist
  }

  # Compute CV distances
  cvdistance <- function(x, cvfolds){
    d_cv <- c()
    d_cv <- lapply(cvfolds, function(y) {
      d_cv_tmp <- sf::st_distance(x[y,], x[-y,])
      c(d_cv,apply(d_cv_tmp, 1, min))
    }) |> unlist()

    d_cv
  }


  # input formatting --------------------------------------------------------

  possible_n <- 2:maxn

  if (!inherits(x, "sf")) {
    x <- sf::st_as_sf(x) |>
      sf::st_geometry() |>
      sf::st_as_sf()
  } else {
    x <- sf::st_geometry(x) |>
      sf::st_as_sf()
  }

  x <- sf::st_transform(x,4326)

  if(inherits(modeldomain,"Raster")) {
    modeldomain <- terra::rast(modeldomain)
  }

  modeldomain <- terra::project(modeldomain, terra::crs(x))


  # calculate NN distance of samples to prediction --------------------------

  # sample prediction points
  modeldomain <- sampleFromArea(modeldomain, samplesize, sampling)

  # Gj: NND function for a cluster per point, i.e. LOO CV
  clust <- 1:nrow(x)
  distmat <- sf::st_distance(x)
  units(distmat) <- NULL
  Gj <- distclust(distmat, clust)

  # Gij: prediction to training NN distances
  Gij <- sf::st_distance(modeldomain, x)
  Gij <- apply(Gij, 1, min)
  units(Gij) <- NULL

  # Check if there's is clustering in the data in the first place
  testks <- stats::ks.test(Gj, Gij, alternative = "great")
  if(testks$p.value >= 0.05){

    # if no evidence of clustering found: random 10 fold CV
    k = 10
    blocks <- sample(rep(1:k, ceiling(nrow(x)/k)), size = nrow(x), replace=F)
    Gjstar <- distclust(distmat, blocks)
    stat_final <- twosamples::wass_stat(Gjstar, Gij)
    cfolds <- CAST::CreateSpacetimeFolds(data.frame(blocks=blocks), spacevar = "blocks", k = k)

    warning("No evidence of clustering has been found, a random CV assignment is returned")


  } else {


    # assign initial groups to the tpoints -----------------------------------

    coords <- sf::st_coordinates(x) |> as.matrix()
    clust <- lapply(possible_n, function(y) {
      c <- stats::kmeans(x=coords, centers = y)
      c$cluster})

    clust <- do.call(cbind.data.frame, clust)
    colnames(clust) <- paste0("groups_", 1:ncol(clust))
    pts_clust <- cbind(x, clust)

    if (n_steps == 1) {

      # CV based on kmeans-clustering
      pts_id_w <- sf::st_drop_geometry(pts_clust)

      spatialfolds <- apply(pts_id_w, 2, function(y) {
        inp <- as.data.frame(y)
        k <- length(unique(inp[[1]]))
        CAST::CreateSpacetimeFolds(inp, spacevar = colnames(inp), k=k)
        })
    } else {


      # assign folds based on probabilities -----------------------------------

      # define probabilities based on n_steps
      pmin <- 1/possible_n
      probs <- lapply(pmin, function(pmin) {
        pr <- seq(pmin, 1, length.out=n_steps) |>
          as.data.frame() |>
          t()
      })

      # sample folds based on those probabilities
      pts_df <- sf::st_drop_geometry(pts_clust)
      pts_f <- lapply(1:length(pmin), function(x) x=pts_df)

      for(i in 1:ncol(pts_df)) {
        for(k in 1:nrow(pts_df)) {

          # use the initial groups as possible folds
          all_folds = unique(pts_df[,i])
          current_fold = pts_df[k,i]
          other_folds = all_folds[all_folds != current_fold]

          for (p in 1:length(probs[[i]])) {

            # use the defined probabilities
            current_prob = probs[[i]][p]
            other_probs =  rep((1-current_prob)/(length(other_folds)),
                               length(other_folds))

            # sample the current fold based on the current probability
            pts_f[[i]][k,paste0("fold_", p)] = sample(c(current_fold, other_folds),
                                                      size = 1,
                                                      prob = c(current_prob, other_probs))
          }
        }
        pts_f[[i]] <- pts_f[[i]][,grepl("fold_", colnames(pts_f[[i]]))]
        colnames(pts_f[[i]]) <- paste0(i+1,"_fold_p_", 1:ncol(pts_f[[i]]))
      }

      # collapse the folds and with colnames: folds _ number of probability step
      pts_id_w <- do.call(cbind.data.frame, pts_f)

      # create caret folds based on those folds
      spatialfolds <- list()
      df_c <- list()
      for (i in 1:ncol(pts_id_w)) {
        df_c[[i]] <- as.data.frame(pts_id_w[,i])
        names(df_c[[i]]) <- paste0("cl_",i)
        pts_f[[i]] <- length(unique(df_c[[i]][!is.na(df_c[[i]])]))
        spatialfolds[[i]] <- CAST::CreateSpacetimeFolds(df_c[[i]],
                                                        spacevar = paste0("cl_",i),
                                                        k=pts_f[[i]])
      }
    }



    # calculate NN distance of CV-folds (Gjstar) ---------------------------------------

    # loop over folds to calculate CV-distance for each block size
    if(!is.null(cl)) {
      parallel::clusterExport(cl=cl, varlist=c("spatialfolds","cvdistance", "x"), envir=environment())
      Gjstar <- parallel::parLapply(cl, spatialfolds, function(y) {
        cvdistance(x=x, cvfolds=y$indexOut)
      })
    } else {
      Gjstar <- lapply(spatialfolds, function(y){
        cvdistance(x=x, cvfolds=y$indexOut) })
    }

    # calculate the WS.distance of Gjstar to Gij for each CV-fold design
    ws_dist <- lapply(Gjstar, function(y) twosamples::wass_stat(Gij, y)) |>
      do.call(what="rbind.data.frame")

    names(ws_dist) <- "D"
    ws_dist$n <- 1:nrow(ws_dist)


    # select the CV-fold design yielding the lowest WS.distance
    n_sel <- ws_dist[ws_dist$D==min(ws_dist$D),][["n"]]

    if(length(n_sel)>1) {
      n_sel <- n_sel[[1]]
    }

    cfolds <- spatialfolds[[n_sel]] # the selected caret folds
    Gjstar <- Gjstar[[n_sel]] # the Gjstar of the selected design
    stat_final <- ws_dist[ws_dist$D==min(ws_dist$D),"D"][[1]] # the WS.distance of the selected design
    blocks <- pts_id_w[,n_sel] # the assignment of tpoints to CV-folds
  }

  # Output
  res <- list(blocks = blocks,
              indx_train = cfolds$index, indx_test = cfolds$indexOut,
              Gij = Gij, Gj = Gj, Gjstar = Gjstar, stat = stat_final)
  class(res) <- c("knndm", "list")
  res

}


