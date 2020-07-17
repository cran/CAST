## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(fig.width = 8.83)
user_hanna <- Sys.getenv("USER") %in% c("hanna")

## ---- message = FALSE, warning=FALSE------------------------------------------
#install.packages("CAST")
library(CAST)

## ---- message = FALSE, warning=FALSE------------------------------------------
help(CAST)

## ---- message = FALSE, warning=FALSE------------------------------------------
data <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
head(data)

## ---- message = FALSE, warning=FALSE------------------------------------------

library(sf)
data_sp <- unique(data[,c("SOURCEID","Easting","Northing")])
data_sp <- st_as_sf(data_sp,coords=c("Easting","Northing"),crs=26911)
plot(data_sp,axes=T,col="black")

## ---- message = FALSE, warning=FALSE, eval=user_hanna-------------------------
#...or plot the data with mapview (recommended!):
library(mapview)
mapviewOptions(basemaps = c("Esri.WorldImagery"))
mapview(data_sp)

## ---- message = FALSE, warning=FALSE------------------------------------------
library(lubridate)
library(ggplot2)
trainDat <- data[data$altitude==-0.3&
                   year(data$Date)==2012&
                   week(data$Date)%in%c(10:12),]
ggplot(data = trainDat, aes(x=Date, y=VW)) +
  geom_line(aes(colour=SOURCEID))

## ---- message = FALSE, warning=FALSE------------------------------------------
library(caret)
predictors <- c("DEM","TWI","Precip_cum","cday",
                "MaxT_wrcc","Precip_wrcc","BLD",
                "Northing","Easting","NDRE.M")
set.seed(10)
model <- train(trainDat[,predictors],trainDat$VW,
               method="rf",tuneGrid=data.frame("mtry"=2),
               importance=TRUE,ntree=50,
               trControl=trainControl(method="cv",number=3))

## ---- message = FALSE, warning=FALSE------------------------------------------
library(raster)
predictors_sp <- stack(system.file("extdata","predictors_2012-03-25.grd",package="CAST"))
prediction <- predict(predictors_sp,model)
spplot(prediction)

## ---- message = FALSE, warning=FALSE------------------------------------------
model

## ---- message = FALSE, warning=FALSE------------------------------------------
set.seed(10)
indices <- CreateSpacetimeFolds(trainDat,spacevar = "SOURCEID",
                                k=3)
set.seed(10)
model_LLO <- train(trainDat[,predictors],trainDat$VW,
                   method="rf",tuneGrid=data.frame("mtry"=2), importance=TRUE,
                   trControl=trainControl(method="cv",
                                          index = indices$index))
model_LLO

## ---- message = FALSE, warning=FALSE------------------------------------------
plot(varImp(model_LLO))

## ---- message = FALSE, warning=FALSE------------------------------------------
set.seed(10)
ffsmodel_LLO <- ffs(trainDat[,predictors],trainDat$VW,metric="Rsquared",
                    method="rf", tuneGrid=data.frame("mtry"=2),
                    verbose=FALSE,ntree=50,
                    trControl=trainControl(method="cv",
                                           index = indices$index))
ffsmodel_LLO
ffsmodel_LLO$selectedvars

## ---- message = FALSE, warning=FALSE------------------------------------------
plot_ffs(ffsmodel_LLO)

## ---- message = FALSE, warning=FALSE------------------------------------------
prediction_ffs <- predict(predictors_sp,ffsmodel_LLO)
spplot(prediction_ffs)

## ---- message = FALSE, warning=FALSE------------------------------------------
library(latticeExtra)

### AOA for which the spatial CV error applies:
AOA <- aoa(predictors_sp,ffsmodel_LLO)

spplot(prediction_ffs,main="prediction for the AOA \n(spatial CV error applied)")+
spplot(AOA$AOA,col.regions=c("grey","transparent"))

### AOA for which the random CV error applies:
AOA_random <- aoa(predictors_sp,model)
spplot(prediction,main="prediction for the AOA \n(random CV error applied)")+
spplot(AOA_random$AOA,col.regions=c("grey","transparent"))


