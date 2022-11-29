#' NN distance based CV
#' @description
#' This function estimates the optimal cross-validation folds.
#' Therefore the function calculates the nearest neighbor distance
#' between training locations and prediction locations, as well as
#' between different cross-validation folds.
#' The resulting numerical vectors are then used to calculate
#' density distributions. The magnitude of difference between
#' the prediction to sample nearest neighbor density distribution to
#' that of different cross-validation folds is quantified using
#' the Kolmogorov-Smirnov Test. The cross-validation folds yielding
#' the lowest magnitude of distance towards the sample to prediction
#' nearest neighbor density distribution is then chosen as the final
#' split.
#' The function returns the density distribution plots,
#' the training points with assigned cross-validation folds,
#' the grid that was used to split the data, as well as
#' a caret object inheriting the folds.
#' @param x A sf object containing the training points.
#' @param modeldomain A raster or terra object containing the predictor variables.
#' @param samplesize Optional. The number of samplepoints generated in the modeldomain to describe the prediction space.
#' If NULL (default), the samplesize equals the number of cells in the predictor raster divided by ten.
#' @param sampling The type in which the samplepoints are generated. The default value is "Fibonacci", which is recommended for global studies.
#' @param w The window size. By default 60.
#' @param n_steps The number of steps until termination.
#' @param maxn The maximal number of clusters to be fitted. One of {4,9,16,25,36,49,...}. By default 36.
#' @param n_cores The number of cores used. By default 6.
#' @details The nearest neighbor distance distribution between training and prediction locations,
#' as well as between cross-validation folds are calculated. Based on that, the
#' cross-validation split yielding the lowest distance between those distributions is chosen.
#' @export block_nndm
#' @aliases block_nndm


block_nndm_par <- function(x=samplepoints,
                   modeldomain=predictors,
                   samplesize=1000,
                   sampling = "regular",
                   w = 60,
                   n_steps = 6,
                   maxn = 36,
                   n_cores=6) {

  x=train_points[[3]]
  modeldomain=raster::stack("data/species_stack.grd")
  samplesize=1000;
  sampling = "regular";
  w = 60;
  n_steps = 6;
  maxn = 36;
  n_cores=6

  # input formatting --------------------------------------------------------

  # input formatting
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

  if(is.null(samplesize)) {
    samplesize = terra::ncell(modeldomain) / 10
  } else {
    samplesize = samplesize
  }


  # calculate nn distance of samples to prediction --------------------------

  # define functions

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
    modeldomainextent <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(modeldomain)))

    sf::sf_use_s2(FALSE)
    sf::st_as_sf(modeldomainextent) |>
      sf::st_transform(4326) -> bb
    methods::as(bb, "Spatial") |>
      sp::spsample(n = samplesize, type = sampling)  |>
      sf::st_as_sfc() |>
      sf::st_set_crs(4326) -> predictionloc

    predictionloc <- sf::st_as_sf(predictionloc)

    return(predictionloc)

  }

  # calculate sample2prediction distance
  sample2prediction <- function(x, modeldomain, samplesize){

    modeldomain <- sf::st_transform(modeldomain, sf::st_crs(x))
    sf::sf_use_s2(TRUE)
    d0 <- sf::st_distance(modeldomain, x)
    min_d0 <- apply(d0, 1, min)
    sampletoprediction <- data.frame(dist = min_d0,
                                     what = factor("sample-to-prediction"),
                                     dist_type = "geo")


    return(sampletoprediction)
  }

  suppressWarnings(modeldomain <- sampleFromArea(modeldomain, samplesize, sampling))

  suppressWarnings(s2p <- sample2prediction(x=x, modeldomain=modeldomain,
                           samplesize=samplesize))


  # split into CV folds ----------------------------------------
  bb <- sf::st_bbox(x) |> sf::st_as_sfc()

  ugt <- sf::st_make_grid(bb, n=w) |> sf::st_as_sf()
  ugt$id <- as.numeric(rownames(ugt))

  # assign groups (folds)
  possible_n <- c(2:sqrt(maxn))^2
  pmin = 1/possible_n
  pmax = 1
  probs = lapply(pmin, seq, 1, length.out = n_steps)
  probs = lapply(probs, function(x) {
    x[[length(x)]] <- 1
    return(x)
  })

  group_mat <- function(x) {
    nrow_m2 = w/nrow(x)
    m2 <- matrix(1,nrow=nrow_m2,ncol=nrow_m2)
    kronecker(x,m2)
  }

  m <- list()
  for (i in seq_along(possible_n)) {
    m[[i]] <- matrix(1:possible_n[[i]], nrow=possible_n[[i]]^0.5,byrow=TRUE)
  }

  m_g <- lapply(m, group_mat)

  # sample random numbers based on probabilities
  cl <- parallel::makeCluster(mc <- getOption("cl.cores", n_cores))
  parallel::clusterExport(cl=cl, varlist=c("m_g", "probs"), envir=environment())

  m_f <- snow::parLapply (cl, seq_along(m_g), function(i) {
   return(lapply (seq_along(probs[[i]]), function(e) {
      vapply(m_g[[i]],function(current_fold) {
        other_folds <- unique(as.vector(m_g[[i]]))[-current_fold]
        current_prob <- probs[[i]][[e]]
        if(current_prob >= 1) {
          other_probs <- rep(0, length(other_folds))
        } else {
          other_probs <-  rep((1-current_prob)/(length(other_folds)),
                              length(other_folds))
        }
        sample(c(current_fold, other_folds),
               size = 1,
               prob = c(current_prob, other_probs))
      },numeric(1))
    }))
  })

  parallel::stopCluster(cl)

  # assign sampled numbers to sf objects
  m_fv <- lapply(m_f, function(x) lapply(x, as.vector))
  m_fd <- lapply(m_fv, function(x) do.call("cbind.data.frame", x))
  m_fd <- lapply(m_fd, function(x) {
    colnames(x) <- paste0("p", 1:ncol(x))
    return(x)})

  ugt_folds <- lapply(m_fd, function(x) {
    sf::st_as_sf(x, sf::st_geometry(ugt))
  })

  # join folds to points
  pts_id <- lapply(ugt_folds, sf::st_join,x=x)
  pts_id_w <- lapply(pts_id,sf::st_drop_geometry)
  for (i in seq_along(pts_id_w)) {
    f <- (i+1)^2
    colnames(pts_id_w[[i]])<- paste0(f,"_fold_p_", 1:ncol(pts_id_w[[i]]))
  }

  pts_id_w <- do.call(cbind.data.frame, pts_id_w)

  spatialfolds <- list()
  df_c <- list()
  for (i in 1:ncol(pts_id_w)) {
    df_c[[i]] <- as.data.frame(pts_id_w[,i])
    names(df_c[[i]]) <- paste0("cl_",i)
    f[[i]] <- nrow(unique(df_c[[i]]))
    spatialfolds[[i]] <- CAST::CreateSpacetimeFolds(df_c[[i]],
                                                    spacevar = paste0("cl_",i),
                                                    k=f[[i]])
  }


  # calculate nn distance of cv folds ---------------------------------------

  cvdistance <- function(x, cvfolds){
    d_cv <- c()
    for (i in 1:length(cvfolds)){
      d_cv_tmp <- sf::st_distance(x[cvfolds[[i]],], x[-cvfolds[[i]],])
      d_cv <- c(d_cv,apply(d_cv_tmp, 1, min))
    }

    dists_cv <- data.frame(dist = d_cv,
                           what = factor("CV-distances"),
                           dist_type = "geo")
    return(dists_cv)
  }


  # loop over folds to calculate CV-distance for each block size
  cl <- parallel::makeCluster(mc <- getOption("cl.cores", n_cores))
  parallel::clusterExport(cl=cl, varlist=c("spatialfolds", "cvdistance", "x"), envir=environment())
  cv_dist <- snow::parLapply(cl, spatialfolds, function(y){
    cvdistance(x=x, cvfolds=y$indexOut)
  })
  parallel::stopCluster(cl)


  # calculate the difference to sample2prediction distribution for each chunk -----

  # filter for cv fold sizes with insufficient size
  n_m = list()
  for(i in seq_along(cv_dist)) {
    n_m[[i]] <- ifelse(nrow(cv_dist[[i]]) < 5, i, NaN)
  }

  rem <- n_m[which(!is.na(n_m))] |> unlist()

  if(!is.null(rem)) {
    cv_dist = cv_dist[-rem]
    possible_n = possible_n[-(rem-2)]
  }

  ws_dist <- lapply(cv_dist, function(y) twosamples::wass_stat(s2p$dist, y$dist)) |>
    do.call(what="rbind.data.frame")

  names(ws_dist) <- "D"
  ws_dist$n <- 1:nrow(ws_dist)

  n_sel <- ws_dist[ws_dist$D==min(ws_dist$D),][["n"]]

  if(length(n_sel)>1) {
    n_sel <- n_sel[[1]]
  }


  # check if random CV yields lower WS stat ---------------------------
  pts_id_rand <- x
  pts_id_rand$n <- sample(1:10, nrow(pts_id_rand), replace=TRUE)
  rand_f <- CAST::CreateSpacetimeFolds(pts_id_rand, spacevar="n",k=10)
  cv_dist_rand <- cvdistance(x=x, cvfolds=rand_f$index)
  ws_dist_rand <- twosamples::wass_stat(s2p$dist, cv_dist_rand$dist)


  # create plot and return grid / assigned points ---------------------------

  if(ws_dist_rand < min(ws_dist$D)) {
    print("random CV")
    spatialfold <- rand_f
    pts_cv <- pts_id_rand
    sf::st_geometry(pts_cv) <- "geom"
    cv_dist_sel <- cv_dist_rand
  } else {
    # caret cv folds
    spatialfolds <- spatialfolds[[n_sel]]

    # grid and assigned points
    pts_cv <- sf::st_as_sf(data.frame(cv_fold=pts_id_w[,n_sel]), pts_id[[1]]$geometry)
    sf::st_geometry(pts_cv) <- "geom"

    # nnd plot and distances
    cv_dist_sel <- cv_dist[[n_sel]]}

  cv_s2p <- rbind(cv_dist_sel, s2p)

  plot.nnd_ecdf <- function(x) {
    xlab="geographic distances"
    ggplot2::ggplot(x, aes(x=dist, group=what, fill=what)) +
      ggplot2::stat_ecdf(aes(ymin=0,ymax=..y..), geom="ribbon", alpha=0.6) +
      ggplot2::stat_ecdf(geom="step") +
      ggplot2::scale_fill_discrete(name = "distance function") +
      ggplot2::xlab(xlab) +
      ggplot2::theme(legend.position="bottom",
                     plot.margin = unit(c(0,0.5,0,0),"cm"))
  }

  plot.nnd <- function(x){
    what <- "" #just to avoid check note
    unit ="m"
    xlab="geographic distances"
    ggplot2::ggplot(data=x, aes(x=as.numeric(dist), group=what, fill=what)) +
      ggplot2::geom_density(adjust=1.5, alpha=.4, stat="density") +
      ggplot2::scale_fill_discrete(name = "distance function") +
      ggplot2::xlab(xlab) +
      ggplot2::theme(legend.position="bottom",
                     plot.margin = unit(c(0,0.5,0,0),"cm"))
  }

  nnd_plot <- plot.nnd(cv_s2p)
  nnd_plot_ecdf <- plot.nnd_ecdf(cv_s2p)
  nnd_distances <- cv_dist_sel

  rtrn <- list(indx_train = spatialfolds$index, indx_test = spatialfolds$indexOut, pts_cv,
               nnd_plot, nnd_plot_ecdf, nnd_distances, ws_dist[ws_dist$D==min(ws_dist$D),"D"])
    names(rtrn) <- c("indx_train","indx_test", "pts_with_cv", "nnd_plot", "nnd_plot_ecdf", "nnd_distances", "WS_distance")
  return(rtrn)

}


