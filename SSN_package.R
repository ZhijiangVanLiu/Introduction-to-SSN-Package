# Author: Zhijiang (Van) Liu
# Date: 4/11/2019

###############################
# Introduction to SSN package #
###############################


# Load SSN library
install.packages('SSN')
vignette('SSN', package = "SSN")
require("SSN")


#################################################
# 4. Manipulating the SpatialStreamNetwork object

#---------------------------------------
# 4.1. Data available in the SSN package

# Save copy of the MiddleFork04 data set provided in the SSN package
# to workding directory
file.copy(system.file("lsndata/MiddleFork04.ssn", package = "SSN"),
          to = getwd(), recursive = TRUE, copy.mode = FALSE)

#--------------------------------------------------------------
# 4.2. Importing and subsetting the SpatialStreamNetwork object

# Import the MiddleFork04 data
mf04p <- importSSN("./MiddleFork04.ssn", predpts = "pred1km")

# Import more prediction sites one set at a time
mf04p <- importPredpts(mf04p, "Knapp", "ssn")
mf04p <- importPredpts(mf04p, "CapeHorn", "ssn")

#-------------------------------------------
# 4.3. Generating an additive function value

# Observe the attributes of the line segments
names(mf04p@data)
head(mf04p@data)
head(mf04p@data[, c("h2oAreaKm2", "afvArea")])
summary(mf04p@data[, c("h2oAreaKm2", "afvArea")])

# Use additive.function to compute the additive function value based on h2oAreaKm2,
# and name it computed.afv
mf04p <- additive.function(mf04p, "h2oAreaKm2", "computed.afv")
names(mf04p@data)
head(mf04p@data[, c("h2oAreaKm2", "afvArea", "computed.afv")])

#-----------------------------------
# 4.4. Calculating distance matrices

# Creates the distance matrices of observed sites
createDistMat(mf04p, o.write = TRUE)
# Creates the relevant parts of the matrix N before fitting models
# or making predictions
createDistMat(mf04p, predpts = "Knapp", o.write = TRUE,
              amongpreds = TRUE)
createDistMat(mf04p, predpts = "CapeHorn", o.write = TRUE,
              amongpreds = TRUE)

# Import distance matrices of observed sites
distObs <- getStreamDistMat(mf04p)
str(distObs)
distObs$dist.net2[5:10,5:10]

# Compute stream distance matrix from asymmetric distance matrix
strDistNet2 <- distObs$dist.net2 + t(distObs$dist.net2)
strDistNet2[5:10,5:10]

# Check flow-connectedness
FCIndNet2 <- 1 - sign(distObs$dist.net2 * t(distObs$dist.net2))
FCIndNet2[5:10,5:10]
FUIndNet2 <- sign(distObs$dist.net2 * t(distObs$dist.net2))
FUIndNet2[5:10,5:10]

# Obtain distance between observed sites and prediction locations
distPred1km <- getStreamDistMat(mf04p, Name = "pred1km")
str(distPred1km)
distPred1km$dist.net1.a[1:5,1:5]
distPred1km$dist.net1.b[1:5,1:5]

# Get distances among just the prediction locations
createDistMat(mf04p, predpts = "CapeHorn", o.write = TRUE,
              amongpreds = TRUE)
distCape <- getStreamDistMat(mf04p, Name = "CapeHorn")
str(distCape)
distCape$dist.net2[1:5,1:5]


#################################################
# 5. Data analysis using SSN

#-------------------------------
# 5.1. Exploratory data analysis
# Obtain the names of the variables in the point.data data.frame for each 
# observed and prediction data set in mf04
names(mf04p)

# The default plotting function for a SpatialStreamNetwork object is a map
plot(mf04p, lwdLineCol = "afvArea", lwdLineEx = 10, lineCol = "green",
     pch = 19, xlab = "x-coordinate (m)", ylab = "y-coordinate (m)",
     asp = 1)

# Summer_mn, is the average summer stream temperature
brks <- plot(mf04p, "Summer_mn", lwdLineCol = "afvArea",
             lwdLineEx = 15, lineCol = "black", xlab = "x-coordinate" ,
             ylab = "y-coordinate", asp = 1, cex = 1.5)

# plot the stream lines
plot(as.SpatialLines(mf04p), col = "blue")
# add the observed locations with size proportional
# to mean summer temperature
plot(as.SpatialPoints(mf04p), pch = 19,
     cex = as.SpatialPointsDataFrame(mf04p)$Summer_mn/9 , add = TRUE)
# add the prediction locations on the 1 km spacing
plot(as.SpatialPoints(mf04p, data = "pred1km"), cex = 1.5, add = TRUE)
# add the dense set of points for block prediction on Knapp segment
plot(as.SpatialPoints(mf04p, data = "Knapp"), pch = 19, cex = 0.3,
     col = "red", add = TRUE)
# add the dense set of points for block prediction on Knapp segment
plot(as.SpatialPoints(mf04p, data = "CapeHorn"), pch = 19, cex = 0.3,
     col = "orange", add = TRUE)

# Torgegram for the mean summer stream temperature
mf04.Torg <- Torgegram(mf04p, "Summer_mn", nlag = 20, maxlag = 50000)
plot(mf04.Torg)

#-------------------
# 5.2. Model fitting

# glmssn0: a non-spatial model
mf04.glmssn0 <- glmssn(Summer_mn ~ ELEV_DEM + SLOPE, mf04p,
                       CorModels = NULL, use.nugget = TRUE)
summary(mf04.glmssn0)
# Compare to OLS model
summary(lm(Summer_mn ~ ELEV_DEM + SLOPE, getSSNdata.frame(mf04p)))

# glmssn1: spatial model including a mixture of tail-up, tail-down, 
#          and Euclidean covariance models
mf04.glmssn1 <- glmssn(Summer_mn ~ ELEV_DEM + SLOPE, mf04p,
                       CorModels = c("Exponential.tailup", "Exponential.taildown",
                                     "Exponential.Euclid"), addfunccol = "afvArea")
summary(mf04.glmssn1)

# The variable MaxOver20 is a binary variable indicating if a stream temperature value was over
# 20 degrees C. Here is an example of fitting a spatial generalized linear model to the
# binary data by using the family = "binomial" argument.
mf04.glmssnBin <- glmssn(MaxOver20 ~ ELEV_DEM + SLOPE, mf04p,
                         CorModels = c("Mariah.tailup", "Spherical.taildown"),
                         family = "binomial", addfunccol = "afvArea")
summary(mf04.glmssnBin)

# The variable C16 represents the number of summer days when the stream temperature was
# greater than 16 degrees C. Here is an example of fitting a spatial generalized linear model
# to these count data by using the family = "poisson" argument
mf04.glmssnPoi <- glmssn(C16 ~ ELEV_DEM + SLOPE, mf04p,
                         CorModels = c("LinearSill.tailup", "LinearSill.taildown"),
                         family = "poisson", addfunccol = "afvArea")
summary(mf04.glmssnPoi)

#-------------------------------
# 5.3. Residuals and diagnostics

# The result of the residuals function is an influenceSSN object, which is
# an exact copy of the glmssn object, except that residual diagnostics are appended as new
# columns to the point.data data.frame containing the observed data. 
mf04.resid1 <- residuals(mf04.glmssn1)
names(getSSNdata.frame(mf04.resid1))

# The default plotting method for an influenceSSN object is a map with color-coded 
# raw residuals
plot(mf04.resid1)

# Histograms can be plotted of the raw response variable and residuals
par(mfrow = c(1, 2))
hist(mf04.resid1)
hist(mf04p, "Summer_mn")

#---------------------
# 5.4. Model selection

AIC(mf04.glmssn0)
AIC(mf04.glmssn1)

#---------------------
# 5.5. Prediction

# Point-wise prediction using kriging
mf04.pred1km <- predict(mf04.glmssn1, "pred1km")
plot(mf04.pred1km, SEcex.max = 1, SEcex.min = .5/3*2,
     breaktype = "user", brks = brks)

# A dense set of points along a single tributary, called Knapp Creek, of the Middle Fork river. 
# We predict all points using the predict function.
plot(mf04p, "Summer_mn", pch = 1, cex = 3,
     xlab = "x-coordinate", ylab = "y-coordinate",
     xlim = c(-1511000,-1500000), ylim = c(2525000,2535000))
mf04.glmssn1.Knapp <- predict(mf04.glmssn1, "Knapp")
plot(mf04.glmssn1.Knapp, "Summer_mn", add = TRUE,
     xlim = c(-1511000,-1500000), ylim = c(2525000,2535000))

# The prediction sites are on a dense evenly-spaced grid because they are used to approximate
# the integrals involved with block prediction
mf04.glmssn1.BPKnapp <- BlockPredict(mf04.glmssn1, "Knapp")
mf04.glmssn1.BPKnapp
mf04.glmssn1.BPCapeHorn <- BlockPredict(mf04.glmssn1, "CapeHorn")
mf04.glmssn1.BPCapeHorn





