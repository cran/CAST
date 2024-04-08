## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)

## ----ffs-parallel-unix, eval = FALSE------------------------------------------
#  data("splotdata")
#  spatial_cv = CreateSpacetimeFolds(splotdata, spacevar = "Biome", k = 5)
#  ctrl <- trainControl(method="cv",index = spatial_cv$index)
#  
#  ffsmodel <- ffs(predictors = splotdata[,6:16],
#                  response = splotdata$Species_richness,
#                  tuneLength = 1,
#                  method = "rf",
#                  trControl = ctrl,
#                  ntree = 20,
#                  seed = 1,
#                  cores = 4)

## ----parallel-cv, eval = FALSE------------------------------------------------
#  library(doParallel)
#  
#  data("splotdata")
#  spatial_cv = CreateSpacetimeFolds(splotdata, spacevar = "Biome", k = 4)
#  ctrl <- trainControl(method="cv",index = spatial_cv$index)
#  
#  cl <- makeCluster(4)
#  registerDoParallel(cl)
#  ffsmodel <- ffs(predictors = splotdata[,6:16],
#                  response = splotdata$Species_richness,
#                  tuneLength = 4,
#                  method = "rf",
#                  trControl = ctrl,
#                  ntree = 20,
#                  seed = 1,
#                  cores = 1)
#  stopCluster(cl)
#  

## ----ffs-ranger, eval = FALSE-------------------------------------------------
#  
#  ffsmodel <- ffs(predictors = splotdata[,6:16],
#                  response = splotdata$Species_richness,
#                  method = "ranger")
#  

## ----message = FALSE, warning=FALSE-------------------------------------------
library(CAST)
library(caret)
library(terra)
library(sf)

data("splotdata")
predictors <- rast(system.file("extdata","predictors_chile.tif",package="CAST"))
splotdata <- st_drop_geometry(splotdata)

## ----message = FALSE, warning=FALSE-------------------------------------------
set.seed(10)
model <- train(splotdata[,names(predictors)],
               splotdata$Species_richness,
               method="rf",
               tuneLength = 1,
               importance=TRUE,
               ntrees = 20,
               trControl = trainControl(method="cv"), number = 5)
prediction <- predict(predictors, model, na.rm=TRUE)

## -----------------------------------------------------------------------------
tdi = trainDI(model, verbose = FALSE)
print(tdi)
class(tdi)
str(tdi)

## ----eval = FALSE-------------------------------------------------------------
#  # you can save the trainDI object for later application
#  saveRDS(tdi, "path/to/file")

## ----fig.show="hold", out.width="30%"-----------------------------------------
r1 = crop(predictors, c(-75.66667, -67, -30, -17.58333))
r2 = crop(predictors, c(-75.66667, -67, -45, -30))
r3 = crop(predictors, c(-75.66667, -67, -55.58333, -45))


plot(r1[[1]],main = "Tile 1")
plot(r2[[1]],main = "Tile 2")
plot(r3[[1]],main = "Tile 3")

## ----fig.show="hold", out.width="30%"-----------------------------------------

aoa_r1 = aoa(newdata = r1, trainDI = tdi)

plot(r1[[1]], main = "Tile 1: Predictors")
plot(aoa_r1$DI, main = "Tile 1: DI")
plot(aoa_r1$AOA, main = "Tile 1: AOA")


## ----eval = FALSE-------------------------------------------------------------
#  
#  library(parallel)
#  
#  tiles_aoa = mclapply(list(r1, r2, r3), function(tile){
#    aoa(newdata = tile, trainDI = tdi)
#  
#  }, mc.cores = 3)
#  

## ----echo = FALSE-------------------------------------------------------------
tiles_aoa = lapply(list(r1, r2, r3), function(tile){
  aoa(newdata = tile, trainDI = tdi)
})

## ----fig.show="hold", out.width="30%"-----------------------------------------
plot(tiles_aoa[[1]]$AOA, main = "Tile 1")
plot(tiles_aoa[[2]]$AOA, main = "Tile 2")
plot(tiles_aoa[[3]]$AOA, main = "Tile 3")

## ----eval = FALSE-------------------------------------------------------------
#  # Simple Example Code for raster tiles on the hard drive
#  
#  tiles = list.files("path/to/tiles", full.names = TRUE)
#  
#  tiles_aoa = mclapply(tiles, function(tile){
#    current = terra::rast(tile)
#    aoa(newdata = current, trainDI = model_random_trainDI)
#  
#  }, mc.cores = 3)

