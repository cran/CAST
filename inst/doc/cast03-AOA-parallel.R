## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)

## ---- message = FALSE, warning=FALSE------------------------------------------
library(CAST)
library(virtualspecies)
library(caret)
library(raster)
library(sp)
library(sf)
library(viridis)
library(latticeExtra)
library(gridExtra)

## ---- message = FALSE, warning=FALSE------------------------------------------
predictors <- stack(system.file("extdata","bioclim.grd",package="CAST"))
response <- generateSpFromPCA(predictors,
                              means = c(3,1),sds = c(2,2), plot=F)$suitab.raster


mask <- predictors[[1]]
values(mask)[!is.na(values(mask))] <- 1
mask <- rasterToPolygons(mask)

# Generate Clustered Training Samples
csample <- function(x,n,nclusters,maxdist,seed){
  set.seed(seed)
  cpoints <- sp::spsample(x, n = nclusters, type="random")
  result <- cpoints
  result$clstrID <- 1:length(cpoints)
  for (i in 1:length(cpoints)){
    ext <- rgeos::gBuffer(cpoints[i,], width = maxdist)
    newsamples <- sp::spsample(ext, n = (n-nclusters)/nclusters, 
                               type="random")
    newsamples$clstrID <- rep(i,length(newsamples))
    result <- rbind(result,newsamples)
    
  }
  result$ID <- 1:nrow(result)
  return(result)
}


samplepoints <- csample(mask,75,15,maxdist=0.20,seed=15)


trainDat <- extract(predictors,samplepoints,df=TRUE)
trainDat$response <- extract (response,samplepoints)
trainDat <- merge(trainDat,samplepoints,by.x="ID",by.y="ID")
trainDat <- trainDat[complete.cases(trainDat),]

## ----message = FALSE, warning=FALSE-------------------------------------------
set.seed(10)
model_random <- train(trainDat[,names(predictors)],
               trainDat$response,
               method="rf",
               importance=TRUE,
               trControl = trainControl(method="cv"))
prediction_random <- raster::predict(predictors,model_random)

## -----------------------------------------------------------------------------
model_random_trainDI = trainDI(model_random)
print(model_random_trainDI)

## ---- eval = FALSE------------------------------------------------------------
#  saveRDS(model_random_trainDI, "path/to/file")

## -----------------------------------------------------------------------------
r1 = crop(predictors, c(0,7,42.33333,54.83333))
r2 = crop(predictors, c(7,14,42.33333,54.83333))
r3 = crop(predictors, c(14,21,42.33333,54.83333))

grid.arrange(spplot(r1[[1]], main = "Tile 1"),
             spplot(r2[[1]], main = "Tile 2"),
             spplot(r3[[1]], main = "Tile 3"), ncol = 3)

## -----------------------------------------------------------------------------

aoa_r1 = aoa(newdata = r1, trainDI = model_random_trainDI)

grid.arrange(spplot(r1[[1]], main = "Tile 1: Predictors"),
             spplot(aoa_r1$DI, main = "Tile 1: DI"),
             spplot(aoa_r1$AOA, main = "Tile 1: AOA"), ncol = 3)


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

## -----------------------------------------------------------------------------
grid.arrange(spplot(tiles_aoa[[1]]$AOA, main = "Tile 1"),
             spplot(tiles_aoa[[2]]$AOA, main = "Tile 2"),
             spplot(tiles_aoa[[3]]$AOA, main = "Tile 3"), ncol = 3)

## ---- eval = FALSE------------------------------------------------------------
#  # Simple Example Code for raster tiles on the hard drive
#  
#  tiles = list.files("path/to/tiles", full.names = TRUE)
#  
#  tiles_aoa = mclapply(tiles, function(tile){
#    current = raster::stack(tile)
#    aoa(newdata = current, trainDI = model_random_trainDI)
#  
#  }, mc.cores = 3)

