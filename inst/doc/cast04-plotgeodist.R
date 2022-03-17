## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,fig.width=6.1, fig.height=3.5)


## ---- message = FALSE, warning=FALSE------------------------------------------
library(CAST)
library(caret)
library(raster)
library(sf)
library(rnaturalearth)
library(ggplot2)

## ----message = FALSE, warning=FALSE-------------------------------------------
ee <- st_crs("+proj=eqearth")
co <- ne_countries(returnclass = "sf")
co.ee <- st_transform(co, ee)

## ----message = FALSE, warning=FALSE, results='hide'---------------------------
sf_use_s2(FALSE)
set.seed(10)
pts_random <- st_sample(co, 1000)
### See points on the map:
ggplot() + geom_sf(data = co.ee, fill="#00BFC4",col="#00BFC4") +
  geom_sf(data = pts_random, color = "#F8766D",size=0.5, shape=3) +
  guides(fill = FALSE, col = FALSE) +
  labs(x = NULL, y = NULL)


## ----message = FALSE, include=FALSE-------------------------------------------
# adjusted from from https://github.com/carlesmila/NNDMpaper/blob/main/code/sim_utils.R
clustered_sample <- function(sarea, nsamples, nparents, radius){

  # Number of offspring per parent
  nchildren <- round((nsamples-nparents)/nparents, 0)

  # Simulate parents
  parents <- st_sf(geometry=st_sample(sarea, nparents, type="random"))
  res <- parents
  res$parent <- 1:nrow(parents)

  # Simulate offspring
  for(i in 1:nrow(parents)){

    # Generate buffer and cut parts outside of the area of study
    buf <- st_buffer(parents[i,], dist=radius)
    buf <- st_intersection(buf, sarea)

    # Simulate children
    children <- st_sf(geometry=st_sample(buf, nchildren, type="random"))
      children$parent <- i
    res <- rbind(res, children)
  }

  return(res)
}

## ----message = FALSE, warning=FALSE, results='hide'---------------------------
set.seed(10)
sf_use_s2(FALSE)
pts_clustered <- clustered_sample(co, 1000, 20, 8)

ggplot() + geom_sf(data = co.ee, fill="#00BFC4",col="#00BFC4") +
  geom_sf(data = pts_clustered, color = "#F8766D",size=0.5, shape=3) +
  guides(fill = FALSE, col = FALSE) +
  labs(x = NULL, y = NULL)


## ----message = FALSE, warning=FALSE, results='hide'---------------------------
dist_random <- plot_geodist(pts_random,co,
                            sampling="Fibonacci",
                            showPlot = FALSE)
dist_clstr <- plot_geodist(pts_clustered,co,
                           sampling="Fibonacci",
                           showPlot = FALSE)

dist_random$plot+scale_x_log10(labels=round)+ggtitle("Randomly distributed reference data")
dist_clstr$plot+scale_x_log10(labels=round)+ggtitle("Clustered reference data")
             

## ----message = FALSE, warning=FALSE, results='hide'---------------------------
randomfolds <- caret::createFolds(1:nrow(pts_clustered))

## ----message = FALSE, warning=FALSE, results='hide',echo=FALSE----------------
for (i in 1:nrow(pts_clustered)){
  pts_clustered$randomCV[i] <- which(unlist(lapply(randomfolds,function(x){sum(x%in%i)}))==1)
}

ggplot() + geom_sf(data = co.ee, fill="#00BFC4",col="#00BFC4") +
  geom_sf(data = pts_clustered, color = rainbow(max(pts_clustered$randomCV))[pts_clustered$randomCV],size=0.5, shape=3) +
  guides(fill = FALSE, col = FALSE) +
  labs(x = NULL, y = NULL)+ggtitle("random fold membership shown by color")

## ----message = FALSE, warning=FALSE, results='hide'---------------------------
dist_clstr <- plot_geodist(pts_clustered,co,
                           sampling="Fibonacci", 
                           cvfolds= randomfolds, 
                           showPlot=FALSE)
dist_clstr$plot+scale_x_log10(labels=round)


## ----message = FALSE, warning=FALSE, results='hide'---------------------------
spatialfolds <- CreateSpacetimeFolds(pts_clustered,spacevar="parent",k=length(unique(pts_clustered$parent)))

## ----message = FALSE, warning=FALSE, results='hide',echo=FALSE----------------
ggplot() + geom_sf(data = co.ee, fill="#00BFC4",col="#00BFC4") +
  geom_sf(data = pts_clustered, color = rainbow(max(pts_clustered$parent))[pts_clustered$parent],size=0.5, shape=3) +
  guides(fill = FALSE, col = FALSE) +
  labs(x = NULL, y = NULL)+ ggtitle("spatial fold membership by color")

## ----message = FALSE, warning=FALSE, results='hide'---------------------------
dist_clstr <- plot_geodist(pts_clustered,co,
                           sampling="Fibonacci",
                           cvfolds= spatialfolds$indexOut, 
                           showPlot=FALSE)
dist_clstr$plot+scale_x_log10(labels=round)

             

## ----message = FALSE, warning=FALSE, results='hide'---------------------------
predictors_global <- stack(system.file("extdata","bioclim_global.grd",package="CAST"))

plot(predictors_global)

## ----message = FALSE, warning=FALSE, results='hide'---------------------------

# use random CV:
dist_clstr_rCV <- plot_geodist(pts_clustered,predictors_global,
                               type = "feature", 
                               sampling="Fibonacci",
                               cvfolds = randomfolds,
                               showPlot=FALSE)

# use spatial CV:
dist_clstr_sCV <- plot_geodist(pts_clustered,predictors_global,
                               type = "feature", sampling="Fibonacci",
                               cvfolds = spatialfolds$indexOut,
                               showPlot=FALSE)

# Plot results:
dist_clstr_rCV$plot+scale_x_log10(labels=round)+ggtitle("Clustered reference data and random CV")
dist_clstr_sCV$plot+scale_x_log10(labels=round)+ggtitle("Clustered reference data and spatial CV")

