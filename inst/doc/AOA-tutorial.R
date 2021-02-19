## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(fig.width = 8.83)

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

## ----message = FALSE,include=FALSE, warning=FALSE-----------------------------
RMSE = function(a, b){
    sqrt(mean((a - b)^2,na.rm=T))
}

## ---- message = FALSE, warning=FALSE------------------------------------------
predictors <- stack(system.file("extdata","bioclim.grd",package="CAST"))
spplot(stretch(predictors,0,1),col.regions=viridis(100))

## ----message = FALSE, warning=FALSE-------------------------------------------
response <- generateSpFromPCA(predictors,
                              means = c(3,1),sds = c(2,2), plot=F)$suitab.raster

## ----message = FALSE, warning=FALSE-------------------------------------------
mask <- predictors[[1]]
values(mask)[!is.na(values(mask))] <- 1
mask <- rasterToPolygons(mask)
set.seed(15)
samplepoints <- spsample(mask,20,"random")
spplot(response,col.regions=viridis(100),
            sp.layout=list("sp.points", samplepoints, col = "red", first = FALSE, cex=2))

## ----message = FALSE, warning=FALSE-------------------------------------------
trainDat <- extract(predictors,samplepoints,df=TRUE)
trainDat$response <- extract (response,samplepoints)
trainDat <- trainDat[complete.cases(trainDat),]

## ----message = FALSE, warning=FALSE-------------------------------------------
set.seed(10)
model <- train(trainDat[,names(predictors)],
               trainDat$response,
               method="rf",
               importance=TRUE,
               trControl = trainControl(method="cv"))
print(model)


## ----message = FALSE, warning=FALSE-------------------------------------------
plot(varImp(model,scale = F),col="black")

## ----message = FALSE, warning=FALSE-------------------------------------------
prediction <- predict(predictors,model)
truediff <- abs(prediction-response)
spplot(stack(prediction,response),main=c("prediction","reference"))

## ----message = FALSE, warning=FALSE-------------------------------------------
AOA <- aoa(predictors, model)
attributes(AOA)$aoa_stats

## ----message = FALSE, warning=FALSE-------------------------------------------
grid.arrange(
  spplot(truediff,col.regions=viridis(100),main="true prediction error"),
  spplot(AOA$DI,col.regions=viridis(100),main="DI"),
  spplot(prediction, col.regions=viridis(100),main="prediction for AOA")+ spplot(AOA$AOA,col.regions=c("grey","transparent")), ncol=3)

## ----clusteredpoints,message = FALSE, include=FALSE---------------------------
#For a clustered sesign:
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

## ----message = FALSE, warning=FALSE-------------------------------------------

samplepoints <- csample(mask,75,15,maxdist=0.20,seed=15)
spplot(response,col.regions=viridis(100),
            sp.layout=list("sp.points", samplepoints, col = "red", first = FALSE, cex=2))

## ----message = FALSE, warning=FALSE-------------------------------------------

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
prediction_random <- predict(predictors,model_random)
print(model_random)

## ----message = FALSE, warning=FALSE-------------------------------------------
folds <- CreateSpacetimeFolds(trainDat, spacevar="clstrID",k=10)
set.seed(15)
model <- train(trainDat[,names(predictors)],
                 trainDat$response,
                     method="rf",
                 importance=TRUE,
                 tuneGrid = expand.grid(mtry = c(2:length(names(predictors)))),
                 trControl = trainControl(method="cv",index=folds$index))
  print(model)
  
prediction <- predict(predictors,model)

## ----message = FALSE, warning=FALSE-------------------------------------------
AOA_spatial <- aoa(predictors, model)

AOA_random <- aoa(predictors, model_random)

## ----message = FALSE, warning=FALSE-------------------------------------------
grid.arrange(spplot(AOA_spatial$DI,col.regions=viridis(100),main="DI"),
  spplot(prediction, col.regions=viridis(100),main="prediction for AOA \n(spatial CV error applies)")+
         spplot(AOA_spatial$AOA,col.regions=c("grey","transparent")),
  spplot(prediction_random, col.regions=viridis(100),main="prediction for AOA \n(random CV error applies)")+
         spplot(AOA_random$AOA,col.regions=c("grey","transparent")),
ncol=3)

## ----message = FALSE, warning=FALSE-------------------------------------------
###for the spatial CV:
RMSE(values(prediction)[values(AOA_spatial$AOA)==1],values(response)[values(AOA_spatial$AOA)==1])
RMSE(values(prediction)[values(AOA_spatial$AOA)==0],values(response)[values(AOA_spatial$AOA)==1])
model$results

###and for the random CV:
RMSE(values(prediction_random)[values(AOA_random$AOA)==1],values(response)[values(AOA_random$AOA)==1])
RMSE(values(prediction_random)[values(AOA_random$AOA)==0],values(response)[values(AOA_random$AOA)==1])
model_random$results

## ----message = FALSE, warning=FALSE-------------------------------------------
AOA_calib <- calibrate_aoa(AOA_spatial,model,window.size = 5,length.out = 5, multiCV=TRUE,showPlot=FALSE)
AOA_calib$plot
spplot(AOA_calib$AOA$expected_RMSE,col.regions=viridis(100),main="expected RMSE")+
 spplot(AOA$AOA,col.regions=c("grey","transparent"))

## ---- message = FALSE, warning=FALSE------------------------------------------
dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
# calculate average of VW for each sampling site:
dat <- aggregate(dat[,c("VW","Easting","Northing")],by=list(as.character(dat$SOURCEID)),mean)
# create sf object from the data:
pts <- st_as_sf(dat,coords=c("Easting","Northing"))

##### Extract Predictors for the locations of the sampling points
studyArea <- stack(system.file("extdata","predictors_2012-03-25.grd",package="CAST"))
st_crs(pts) <- crs(studyArea)
trainDat <- extract(studyArea,pts,df=TRUE)
pts$ID <- 1:nrow(pts)
trainDat <- merge(trainDat,pts,by.x="ID",by.y="ID")
# The final training dataset with potential predictors and VW:
head(trainDat)

## ---- message = FALSE, warning=FALSE------------------------------------------
predictors <- c("DEM","NDRE.Sd","TWI","Bt")
response <- "VW"

model <- train(trainDat[,predictors],trainDat[,response],
               method="rf",tuneLength=3,importance=TRUE,
               trControl=trainControl(method="LOOCV"))
model

## ---- message = FALSE, warning=FALSE------------------------------------------
#Predictors:
spplot(stretch(studyArea[[predictors]]))

#prediction:
prediction <- predict(studyArea,model)

## ---- message = FALSE, warning=FALSE------------------------------------------
AOA <- aoa(studyArea,model)

#### Plot results:
grid.arrange(spplot(AOA$DI,col.regions=viridis(100),main="DI with sampling locations (red)")+
               spplot(as_Spatial(pts),zcol="ID",col.regions="red"),
  spplot(prediction, col.regions=viridis(100),main="prediction for AOA \n(LOOCV error applies)")+ spplot(AOA$AOA,col.regions=c("grey","transparent")),ncol=2)


