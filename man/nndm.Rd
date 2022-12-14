% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nndm.R
\name{nndm}
\alias{nndm}
\title{Nearest Neighbour Distance Matching (NNDM) algorithm}
\usage{
nndm(
  tpoints,
  modeldomain = NULL,
  ppoints = NULL,
  samplesize = 1000,
  sampling = "regular",
  phi = "max",
  min_train = 0
)
}
\arguments{
\item{tpoints}{sf or sfc point object. Contains the training points samples.}

\item{modeldomain}{raster or sf object defining the prediction area (see Details).}

\item{ppoints}{sf or sfc point object. Contains the target prediction points. Optional. Alternative to modeldomain (see Details).}

\item{samplesize}{numeric. How many points in the modeldomain should be sampled as prediction points? Only required if modeldomain is used instead of ppoints.}

\item{sampling}{character. How to draw prediction points from the modeldomain? See spsample. Use sampling = "Fibonacci" for global applications." Only required if modeldomain is used instead of ppoints.}

\item{phi}{Numeric. Estimate of the landscape autocorrelation range in the
same units as the tpoints and ppoints for projected CRS, in meters for geographic CRS.
Per default (phi="max"), the size of the prediction area is used. See Details.}

\item{min_train}{Numeric between 0 and 1. Minimum proportion of training
data that must be used in each CV fold. Defaults to 0 (i.e. no restrictions).}
}
\value{
An object of class \emph{nndm} consisting of a list of six elements:
indx_train, indx_test, and indx_exclude (indices of the observations to use as
training/test/excluded data in each NNDM LOO CV iteration), Gij (distances for
G function construction between prediction and target points), Gj
(distances for G function construction during LOO CV), Gjstar (distances
for modified G function during NNDM LOO CV), phi (landscape autocorrelation range).
indx_train and indx_test can directly be used as "index" and "indexOut" in
caret's \code{\link{trainControl}} function or used to initiate a custom validation strategy in mlr3.
}
\description{
This function implements the \emph{NNDM} algorithm and returns the necessary
indices to perform a NNDM LOO CV for map validation.
}
\details{
NNDM proposes a LOO CV scheme such that the nearest neighbour distance distribution function between the test and training data during the CV process is matched to the nearest neighbour
distance distribution function between the prediction and training points. Details of the method can be found in Mil?? et al. (2022).

Specifying \emph{phi} allows limiting distance matching to the area where this is assumed to be relevant due to spatial autocorrelation.
Distances are only matched up to \emph{phi}. Beyond that range, all data points are used for training, without exclusions.
When \emph{phi} is set to "max", nearest neighbor distance matching is performed for the entire prediction area. Euclidean distances are used for projected
and non-defined CRS, great circle distances are used for geographic CRS (units in meters).

The \emph{modeldomain} is a sf polygon or a raster that defines the prediction area. The function takes a regular point sample (amount defined by \emph{samplesize)} from the spatial extent.
As an alternative use \emph{ppoints} instead of \emph{modeldomain}, if you have already defined the prediction locations.
}
\note{
NNDM is a variation of LOOCV and therefore may take a long time for large training data sets.
You may need to consider alternatives following the ideas of Mil?? et al. (2022) for large data sets.
}
\examples{
library(sf)

# Simulate 100 random training and test points in a 100x100 square
set.seed(123)
poly <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
sample_poly <- sf::st_polygon(poly)
train_points <- sf::st_sample(sample_poly, 100, type = "random")
pred_points <- sf::st_sample(sample_poly, 100, type = "random")

# Run NNDM for the whole domain, here the prediction points are known.
nndm_pred <- nndm(train_points, ppoints=pred_points)
nndm_pred
plot(nndm_pred)

# ...or run NNDM with a known autocorrelation range of 10
# to restrict the matching to distances lower than that.
nndm_pred <- nndm(train_points, ppoints=pred_points, phi = 10)
nndm_pred
plot(nndm_pred)


########################################################################
# Example 2: Real- world example; using a modeldomain instead of previously
# sampled prediction locations
########################################################################
\dontrun{
library(raster)

### prepare sample data:
dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
dat <- aggregate(dat[,c("DEM","TWI", "NDRE.M", "Easting", "Northing","VW")],
   by=list(as.character(dat$SOURCEID)),mean)
pts <- dat[,-1]
pts <- st_as_sf(pts,coords=c("Easting","Northing"))
st_crs(pts) <- 26911
studyArea <- raster::stack(system.file("extdata","predictors_2012-03-25.grd",package="CAST"))

nndm_folds <- nndm(pts, modeldomain= studyArea)
plot(nndm_folds)

#use for cross-validation:
library(caret)
ctrl <- trainControl(method="cv",
   index=nndm_folds$indx_train,
   indexOut=nndm_folds$indx_test,
   savePredictions='final')
model_nndm <- train(dat[,c("DEM","TWI", "NDRE.M")],
   dat$VW,
   method="rf",
   trControl = ctrl)
model_nndm
global_validation(model_nndm)
}

}
\references{
\itemize{
\item Mil??, C., Mateu, J., Pebesma, E., Meyer, H. (2022): Nearest Neighbour Distance Matching Leave-One-Out Cross-Validation for map validation. Methods in Ecology and Evolution 00, 1??? 13.
\item Meyer, H., Pebesma, E. (2022): Machine learning-based global maps of ecological variables and the challenge of assessing them. Nature Communications. 13.
}
}
\seealso{
\code{\link{plot_geodist}}
}
\author{
Carles Mil??
}
