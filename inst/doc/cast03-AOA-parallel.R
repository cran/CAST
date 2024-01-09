## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)

## ----message = FALSE, warning=FALSE-------------------------------------------
library(CAST)
library(caret)
library(terra)
library(sf)

## ----message = FALSE, warning=FALSE-------------------------------------------
data("splotdata")
predictors <- rast(system.file("extdata","predictors_chile.tif",package="CAST"))
splotdata <- st_drop_geometry(splotdata)

## ----message = FALSE, warning=FALSE-------------------------------------------
set.seed(10)
model_random <- train(splotdata[,names(predictors)],
                      splotdata$Species_richness,
                      method="rf",
                      importance=TRUE,
                      ntrees = 50,
                      trControl = trainControl(method="cv"))
prediction_random <- predict(predictors,model_random,na.rm=TRUE)

## -----------------------------------------------------------------------------
model_random_trainDI = trainDI(model_random)
print(model_random_trainDI)

## ----eval = FALSE-------------------------------------------------------------
#  saveRDS(model_random_trainDI, "path/to/file")

## ----fig.show="hold", out.width="30%"-----------------------------------------
r1 = crop(predictors, c(-75.66667, -67, -30, -17.58333))
r2 = crop(predictors, c(-75.66667, -67, -45, -30))
r3 = crop(predictors, c(-75.66667, -67, -55.58333, -45))


plot(r1[[1]],main = "Tile 1")
plot(r2[[1]],main = "Tile 2")
plot(r3[[1]],main = "Tile 3")

## ----fig.show="hold", out.width="30%"-----------------------------------------

aoa_r1 = aoa(newdata = r1, trainDI = model_random_trainDI)

plot(r1[[1]], main = "Tile 1: Predictors")
plot(aoa_r1$DI, main = "Tile 1: DI")
plot(aoa_r1$AOA, main = "Tile 1: AOA")


## ----eval = FALSE-------------------------------------------------------------
#  
#  library(parallel)
#  
#  tiles_aoa = mclapply(list(r1, r2, r3), function(tile){
#    aoa(newdata = tile, trainDI = model_random_trainDI)
#  
#  }, mc.cores = 3)
#  

## ----echo = FALSE-------------------------------------------------------------
tiles_aoa = lapply(list(r1, r2, r3), function(tile){
  aoa(newdata = tile, trainDI = model_random_trainDI)
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

