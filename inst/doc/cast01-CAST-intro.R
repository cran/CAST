## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(fig.width = 8.83)

## ----message = FALSE, warning=FALSE-------------------------------------------
#install.packages("CAST")
library(CAST)

## ----message = FALSE, warning=FALSE-------------------------------------------
help(CAST)

## ----message = FALSE, warning=FALSE-------------------------------------------
library(geodata)
library(terra)
library(sf)
library(caret)
library(tmap)
library(viridis)

## ----message = FALSE, warning=FALSE-------------------------------------------
data(splotdata)
head(splotdata)

## ----message = FALSE, warning=FALSE-------------------------------------------

wc <- worldclim_global(var="bio",res = 10,path=tempdir())
elev <- elevation_global(res = 10, path=tempdir())
predictors_sp <- crop(c(wc,elev),st_bbox(splotdata))
names(predictors_sp) <- c(paste0("bio_",1:19),"elev")

#note: if you prefer to work on a smaller dataset and not download any data, 
#here is a subset of Chile:
#predictors_sp <- rast(system.file("extdata", "predictors_chile.tif",package="CAST"))

## ----message = FALSE, warning=FALSE-------------------------------------------
plot(predictors_sp$elev)
plot(splotdata[,"Species_richness"],add=T)

## ----message = FALSE, warning=FALSE-------------------------------------------
predictors <- c("bio_1", "bio_4", "bio_5", "bio_6",
                "bio_8", "bio_9", "bio_12", "bio_13",
                "bio_14", "bio_15", "elev")

# note that to use the data for model training we have to get rid of the 
# geometry column of the sf object: st_drop_geometry(splotdata)

set.seed(10) # set seed to reproduce the model
model_default <- train(st_drop_geometry(splotdata)[,predictors],
               st_drop_geometry(splotdata)$Species_richness,
               method="rf",tuneGrid=data.frame("mtry"=2),
               importance=TRUE, ntree=50,
               trControl=trainControl(method="cv",number=3, savePredictions = "final"))

## ----message = FALSE, warning=FALSE-------------------------------------------
prediction <- predict(predictors_sp,model_default,na.rm=TRUE)
plot(prediction)

## ----message = FALSE, warning=FALSE-------------------------------------------
model_default
global_validation(model_default)

## ----message = FALSE, warning=FALSE-------------------------------------------
set.seed(10)
indices_LLO <- CreateSpacetimeFolds(splotdata,spacevar = "Country",
                                k=3)

## ----message = FALSE, warning=FALSE-------------------------------------------
set.seed(10)
indices_knndm <- knndm(splotdata,predictors_sp,k=3)

## ----message = FALSE, warning=FALSE-------------------------------------------
plot(geodist(splotdata,predictors_sp,cvfolds =model_default$control$indexOut))+ 
  scale_x_log10(labels=round)
plot(geodist(splotdata,predictors_sp,cvfolds =indices_LLO$indexOut))+ 
  scale_x_log10(labels=round)
plot(geodist(splotdata,predictors_sp,cvfolds =indices_knndm$indx_test))+ 
  scale_x_log10(labels=round)

## ----message = FALSE, warning=FALSE-------------------------------------------
model <- train(st_drop_geometry(splotdata)[,predictors],
                   st_drop_geometry(splotdata)$Species_richness,
                   method="rf",
                   tuneGrid=data.frame("mtry"=2), 
                   importance=TRUE,
                   trControl=trainControl(method="cv",
                                          index = indices_knndm$indx_train,
                                          savePredictions = "final"))
model
global_validation(model)

## ----message = FALSE, warning=FALSE-------------------------------------------
plot(varImp(model))

## ----message = FALSE, warning=FALSE-------------------------------------------
set.seed(10)
ffsmodel <- ffs(st_drop_geometry(splotdata)[,predictors],
                    st_drop_geometry(splotdata)$Species_richness,
                    method="rf", 
                    tuneGrid=data.frame("mtry"=2),
                    verbose=FALSE,
                    ntree=50,
                    trControl=trainControl(method="cv",
                                           index = indices_knndm$indx_train,
                                           savePredictions = "final"))
ffsmodel
global_validation(ffsmodel)
ffsmodel$selectedvars

## ----message = FALSE, warning=FALSE-------------------------------------------
plot(ffsmodel)

## ----message = FALSE, warning=FALSE-------------------------------------------
prediction_ffs <- predict(predictors_sp, ffsmodel, na.rm=TRUE)
plot(prediction_ffs)

## ----message = FALSE, warning=FALSE-------------------------------------------
### AOA for which the spatial CV error applies:
AOA <- aoa(predictors_sp,ffsmodel,LPD = TRUE,verbose=FALSE)

tm_shape(prediction)+
  tm_raster(title="Species \nrichness",palette=viridis(50),style="cont")+
  tm_shape(AOA$AOA)+
  tm_raster(palette=c("1"=NA,"0"="grey"),style="cat",legend.show = FALSE)+
  tm_layout(frame=FALSE,legend.outside = TRUE)+
  tm_add_legend(type="fill",col="grey",border.lwd=0, labels="Outside \nAOA")
    


## ----message = FALSE, warning=FALSE-------------------------------------------
plot(c(AOA$DI,AOA$LPD))

errormodel_DI <- errorProfiles(model_default,AOA,variable="DI")
errormodel_LPD <- errorProfiles(model_default,AOA,variable="LPD")

plot(errormodel_DI)
plot(errormodel_LPD)


## ----message = FALSE, warning=FALSE-------------------------------------------
expected_error_LPD = terra::predict(AOA$LPD, errormodel_LPD)
plot(expected_error_LPD)

