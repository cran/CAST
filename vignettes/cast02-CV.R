## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE,
                      fig.align = 'center')

## ----load packages------------------------------------------------------------
library("sf")
library("CAST")
library("caret")
library("gridExtra")
library("ggplot2")
library("knitr")

set.seed(1234)

## ----read data, fig.width=5, fig.height=7-------------------------------------
# Read data
temperature <- read_sf("https://github.com/carlesmila/RF-spatial-proxies/raw/main/data/temp/temp_train.gpkg")
pm25 <- read_sf("https://github.com/carlesmila/RF-spatial-proxies/raw/main/data/AP/PM25_train.gpkg")
spain <- read_sf("https://github.com/carlesmila/RF-spatial-proxies/raw/main/data/boundaries/spain.gpkg")

# df versions
temperature_df <- as.data.frame(st_drop_geometry(temperature))
pm25_df <- as.data.frame(st_drop_geometry(pm25))

# Plot them
p1 <- ggplot() +
  geom_sf(data = spain, fill = "grey", alpha = 0.1) +
  geom_sf(data = temperature, aes(col = temp)) +
  scale_colour_distiller(palette = "RdYlBu") +
  theme_bw() +
  labs(col = "") +
  ggtitle("Average 2019 temperature (ºC)")
p2 <- ggplot() +
  geom_sf(data = spain, fill = "grey", alpha = 0.1) +
  geom_sf(data = pm25, aes(col = PM25)) +
  scale_colour_viridis_c(option = "cividis") +
  theme_bw() +
  labs(col = "") +
  ggtitle(expression(Average~2019~PM[2.5]~(mu*g/m^3)))
grid.arrange(p1, p2, nrow=2)

## ----geodist density, fig.width=6, fig.height=7-------------------------------
# Random 5-fold CV
fold5_temp <- createFolds(1:nrow(temperature), k=5, returnTrain=FALSE)
fold5_pm25 <- createFolds(1:nrow(pm25), k=5, returnTrain=FALSE)

# Explore geographic predictive conditions
predcond_temp <- geodist(temperature, modeldomain = spain, cvfolds = fold5_temp)
predcond_pm25 <- geodist(pm25, modeldomain = spain, cvfolds = fold5_pm25)

# Plot density functions
p1 <- plot(predcond_temp) + ggtitle("Temperature")
p2 <- plot(predcond_pm25) + ggtitle(expression(PM[2.5]))
grid.arrange(p1, p2, nrow=2)

## ----geodist ecdf, fig.width=6, fig.height=7----------------------------------
# Plot ECDF functions
p1 <- plot(predcond_temp, stat = "ecdf") + ggtitle("Temperature")
p2 <- plot(predcond_pm25, stat = "ecdf") + ggtitle(expression(PM[2.5]))
grid.arrange(p1, p2, nrow=2)

## ----NNDM LOO temp, fig.width=6, fig.height=4---------------------------------
temp_nndm <- nndm(temperature, modeldomain = spain, samplesize = 1000)
print(temp_nndm)
plot(temp_nndm, type = "simple")

## ----NNDM LOO pm25-1, fig.width=6, fig.height=4-------------------------------
pm25_nndm <- nndm(pm25, modeldomain = spain, samplesize = 1000)
print(pm25_nndm)
plot(pm25_nndm, type = "simple")

## ----NNDM LOO pm25-2, fig.height=4, fig.width=6, echo=FALSE-------------------
# The CV iteration with the most excluded data
id_point <- which.max(sapply(pm25_nndm$indx_exclude, length)) 
pm25_plot <- pm25
pm25_plot$set <- ""
pm25_plot$set[pm25_nndm$indx_train[[id_point]]] <- "train"
pm25_plot$set[pm25_nndm$indx_exclude[[id_point]]] <- "exclude"
pm25_plot$set[pm25_nndm$indx_test[[id_point]]] <- "test"
pm25_plot <- pm25_plot[order(pm25_plot$set),] # highlight test point

# And plot
ggplot() +
  geom_sf(data = spain, fill = "grey", alpha = 0.1) +
  geom_sf(data = pm25_plot, aes(col = set)) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw()

## ----model fitting temp NNDM--------------------------------------------------
# LOO CV
temp_loo_ctrl <- trainControl(method="LOOCV", savePredictions=TRUE)
temp_loo_mod <- train(temperature_df[c("dem", "ndvi", "lst_day",  "lst_night")],
                      temperature_df[,"temp"],
                      method="rf", importance=FALSE,
                      trControl=temp_loo_ctrl, ntree=100, tuneLength=1)
temp_loo_res <- global_validation(temp_loo_mod)

# NNDM LOO CV
temp_nndm_ctrl <- trainControl(method="cv", 
                               index=temp_nndm$indx_train, # Obs to fit the model to
                               indexOut=temp_nndm$indx_test, # Obs to validate
                               savePredictions=TRUE)
temp_nndm_mod <- train(temperature_df[c("dem", "ndvi", "lst_day",  "lst_night")],
                       temperature_df[,"temp"],
                       method="rf", importance=FALSE,
                       trControl=temp_nndm_ctrl, ntree=100, tuneLength=1)
temp_nndm_res <- global_validation(temp_nndm_mod)

## ----model fitting pm25 NNDM--------------------------------------------------
# LOO CV
pm25_loo_ctrl <- trainControl(method="LOOCV", savePredictions=TRUE)
pm25_loo_mod <- train(pm25_df[c("popdens", "primaryroads", "ntl", "imd")],
                      pm25_df[,"PM25"],
                      method="rf", importance=FALSE,
                      trControl=pm25_loo_ctrl, ntree=100, tuneLength=1)
pm25_loo_res <- global_validation(pm25_loo_mod)

# NNDM LOO CV
pm25_nndm_ctrl <- trainControl(method="cv", 
                               index=pm25_nndm$indx_train, # Obs to fit the model to
                               indexOut=pm25_nndm$indx_test, # Obs to validate
                               savePredictions=TRUE)
pm25_nndm_mod <- train(pm25_df[c("popdens", "primaryroads", "ntl", "imd")],
                       pm25_df[,"PM25"],
                       method="rf", importance=FALSE,
                       trControl=pm25_nndm_ctrl, ntree=100, tuneLength=1)
pm25_nndm_res <- global_validation(pm25_nndm_mod)

## ----parse results LOO CV, echo=FALSE-----------------------------------------
rbind(
  data.frame(outcome="Temperature", validation="LOO CV",
             t(as.data.frame(temp_loo_res))),
  data.frame(outcome="Temperature", validation="NNDM LOO CV",
             t(as.data.frame(temp_nndm_res))),
  data.frame(outcome="PM2.5", validation="LOO CV",
             t(as.data.frame(pm25_loo_res))),
  data.frame(outcome="PM2.5", validation="NNDM LOO CV",
             t(as.data.frame(pm25_nndm_res)))
) |> 
  kable(digits=2, row.names = FALSE)

## ----kNNDM temp, fig.width=6, fig.height=4------------------------------------
temp_knndm <- knndm(temperature, k = 5, modeldomain = spain, samplesize = 1000, 
                   clustering = "hierarchical", linkf = "ward.D2")
print(temp_knndm)
plot(temp_nndm, type = "simple")

## ----kNNDM pm25, fig.width=6, fig.height=4------------------------------------
pm25_knndm <- knndm(pm25, k = 5, modeldomain = spain, samplesize = 1000, 
                    clustering = "hierarchical", linkf = "ward.D2")
print(pm25_knndm)
plot(pm25_knndm, type = "simple")

## ----kNNDM pm25 v2, fig.width=6, fig.height=4---------------------------------
pm25_knndm_v2 <- knndm(pm25, k = 5, modeldomain = spain, samplesize = 1000, 
                       clustering = "kmeans")
print(pm25_knndm_v2)

## ----kNNDM viz, fig.width=7, fig.height=3, echo=FALSE-------------------------
p1 <- ggplot() +
  geom_sf(data = spain, fill = "grey", alpha = 0.1) +
  geom_sf(data = temperature, aes(col = as.factor(temp_knndm$clusters))) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("Temperature kNNDM") +
  theme_bw() + theme(legend.position = "none")
p2 <- ggplot() +
  geom_sf(data = spain, fill = "grey", alpha = 0.1) +
  geom_sf(data = pm25, aes(col = as.factor(pm25_knndm$clusters))) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle(expression(PM[2.5]~kNNDM)) +
  theme_bw() + theme(legend.position = "none")
grid.arrange(p1, p2, nrow = 1)

## ----model fitting temp kNNDM-------------------------------------------------
# Random 5-fold CV
temp_rndmk_ctrl <- trainControl(method="cv", number=5, savePredictions=TRUE)
temp_rndmk_mod <- train(temperature_df[c("dem", "ndvi", "lst_day",  "lst_night")],
                      temperature_df[,"temp"],
                      method="rf", importance=FALSE,
                      trControl=temp_rndmk_ctrl, ntree=100, tuneLength=1)
temp_rndmk_res <- global_validation(temp_rndmk_mod)

# kNNDM 5-fold CV
temp_knndm_ctrl <- trainControl(method="cv", 
                                index=temp_knndm$indx_train,
                                savePredictions=TRUE)
temp_knndm_mod <- train(temperature_df[c("dem", "ndvi", "lst_day",  "lst_night")],
                        temperature_df[,"temp"],
                        method="rf", importance=FALSE,
                        trControl=temp_knndm_ctrl, ntree=100, tuneLength=1)
temp_knndm_res <- global_validation(temp_knndm_mod)

## ----model fitting pm25 kNNDM-------------------------------------------------
# Random 5-fold CV
pm25_rndmk_ctrl <- trainControl(method="cv", number=5, savePredictions=TRUE)
pm25_rndmk_mod <- train(pm25_df[c("popdens", "primaryroads", "ntl", "imd")],
                        pm25_df[,"PM25"],
                        method="rf", importance=FALSE,
                        trControl=pm25_rndmk_ctrl, ntree=100, tuneLength=1)
pm25_rndmk_res <- global_validation(pm25_rndmk_mod)

# kNNDM 5-fold CV
pm25_knndm_ctrl <- trainControl(method="cv", 
                                index=pm25_knndm$indx_train,
                                savePredictions=TRUE)
pm25_knndm_mod <- train(pm25_df[c("popdens", "primaryroads", "ntl", "imd")],
                        pm25_df[,"PM25"],
                        method="rf", importance=FALSE,
                        trControl=pm25_knndm_ctrl, ntree=100, tuneLength=1)
pm25_knndm_res <- global_validation(pm25_knndm_mod)

## ----parse results kfold CV, echo=FALSE---------------------------------------
rbind(
  data.frame(outcome="Temperature", validation="Random 5-fold",
             t(as.data.frame(temp_rndmk_res))),
  data.frame(outcome="Temperature", validation="kNNDM 5-fold",
             t(as.data.frame(temp_knndm_res))),
  data.frame(outcome="PM2.5", validation="Random 5-fold",
             t(as.data.frame(pm25_rndmk_res))),
  data.frame(outcome="PM2.5", validation="kNNDM 5-fold",
             t(as.data.frame(pm25_knndm_res)))
) |> 
  kable(digits=2, row.names = FALSE)

