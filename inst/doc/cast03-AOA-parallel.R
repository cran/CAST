## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)

## ---- message = FALSE, warning=FALSE------------------------------------------
library(CAST)
library(caret)
library(terra)
#library(sp)
library(sf)
library(viridis)
library(latticeExtra)
library(gridExtra)

## ----message = FALSE, warning=FALSE-------------------------------------------

generate_random_response <- function(raster, predictornames =
names(raster), seed = sample(seq(1000), 1)){
  operands_1 = c("+", "-", "*", "/")
  operands_2 = c("^1","^2")
  
  expression <- paste(as.character(predictornames, sep=""))
  # assign random power to predictors
  set.seed(seed)
  expression <- paste(expression,
                      sample(operands_2, length(predictornames),
replace = TRUE),
                      sep = "")
  
  # assign random math function between predictors (expect after the last one)
  set.seed(seed)
  expression[-length(expression)] <- paste(expression[-
length(expression)],
                                           sample(operands_1,
length(predictornames)-1, replace = TRUE),
                                           sep = " ")
  print(paste0(expression, collapse = " "))
  # collapse
  e = paste0("raster$", expression, collapse = " ")
  
  response = eval(parse(text = e))
  names(response) <- "response"
  return(response)
  
}

## ---- message = FALSE, warning=FALSE------------------------------------------
predictors <- rast(system.file("extdata","bioclim.grd",package="CAST"))
response <- generate_random_response (predictors, seed = 10)


mask <- predictors[[1]]
values(mask)[!is.na(values(mask))] <- 1
mask <- st_as_sf(as.polygons(mask))
mask <- st_make_valid(mask)

set.seed(15)
samplepoints <- clustered_sample(mask,75,15,radius=25000)


trainDat <- extract(predictors,samplepoints,na.rm=TRUE)
trainDat$response <- extract(response,samplepoints,na.rm=FALSE)$response
trainDat <- data.frame(trainDat,samplepoints)
trainDat <- na.omit(trainDat)

## ----message = FALSE, warning=FALSE-------------------------------------------
set.seed(10)
model_random <- train(trainDat[,names(predictors)],
               trainDat$response,
               method="rf",
               importance=TRUE,
               trControl = trainControl(method="cv"))
prediction_random <- predict(predictors,model_random,na.rm=TRUE)

## -----------------------------------------------------------------------------
model_random_trainDI = trainDI(model_random)
print(model_random_trainDI)

## ---- eval = FALSE------------------------------------------------------------
#  saveRDS(model_random_trainDI, "path/to/file")

## ----  fig.show="hold", out.width="30%"---------------------------------------
r1 = crop(predictors, c(3496791,4073906,2143336,3579086))
r2 = crop(predictors, c(4073906,4651021,2143336,3579086))
r3 = crop(predictors, c(4651021,5228136,2143336,3579086))


plot(r1[[1]],main = "Tile 1")
plot(r2[[1]],main = "Tile 2")
plot(r3[[1]],main = "Tile 3")

## ----  fig.show="hold", out.width="30%"---------------------------------------

aoa_r1 = aoa(newdata = r1, trainDI = model_random_trainDI)

plot(r1[[1]], main = "Tile 1: Predictors")
plot(aoa_r1$DI, main = "Tile 1: DI")
plot(aoa_r1$AOA, main = "Tile 1: AOA")


## ---- eval = FALSE------------------------------------------------------------
#  
#  library(parallel)
#  
#  tiles_aoa = mclapply(list(r1, r2, r3), function(tile){
#    aoa(newdata = tile, trainDI = model_random_trainDI)
#  
#  }, mc.cores = 3)
#  

## ---- echo = FALSE------------------------------------------------------------
tiles_aoa = lapply(list(r1, r2, r3), function(tile){
  aoa(newdata = tile, trainDI = model_random_trainDI)
})

## ----  fig.show="hold", out.width="30%"---------------------------------------
plot(tiles_aoa[[1]]$AOA, main = "Tile 1")
plot(tiles_aoa[[2]]$AOA, main = "Tile 2")
plot(tiles_aoa[[3]]$AOA, main = "Tile 3")

## ---- eval = FALSE------------------------------------------------------------
#  # Simple Example Code for raster tiles on the hard drive
#  
#  tiles = list.files("path/to/tiles", full.names = TRUE)
#  
#  tiles_aoa = mclapply(tiles, function(tile){
#    current = terra::rast(tile)
#    aoa(newdata = current, trainDI = model_random_trainDI)
#  
#  }, mc.cores = 3)

