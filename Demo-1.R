libray(terra)
library(terra)
#loading libraries
library(tidyverse)
library(sf)
library(caret)
# setting Working directory
setwd("C:\\Users\\ADMIN\\Desktop\\UAV")

setwd("C:\\Users\\ADMIN\\Desktop\\UAV\\DS4-20260624")
#Read raster image

ms_image <- rast("DS4_UAV_Multispectral_Image.tif")
ms_image

#plot the raster
plot(ms_image)

#to asssign band names (if they are not available)
names(ms_image) <- c("green","red","rededge","nir")

#plot multicolor composite
plotRGB(ms_image,r="red",g="green",b="rededge",stretch="hist",main="True_colour_composite")

#read vector polygon
aoi <- st_read("DS4_AOI.gpkg")
aoi
plot(aoi)

subplots <-sf::st_read("DS4_Subplots.gpkg")
subplots
plot(subplots)
plot(subplots["layer"]) #plots subplots according to data

# plot together Field AOI and subplots
plot(aoi)
plot(subplots["layer"],add=TRUE)

#cutout raster into field boundary
ms_image_crop <- terra::crop(ms_image,aoi)
ms_image_mask <-terra::mask(ms_image_crop,aoi)
plotRGB(ms_image_mask,r="red",g="green",b="rededge",stretch="hist",main="True_colour_composite")

#VI calculation

eps<-1e-10 #value to avoid zero divition
ndvi <-(ms_image_mask$nir-ms_image_mask$red)/(ms_image_mask$nir+ms_image_mask$red)
ndvi
names(ndvi)<-"ndvi"
plot(ndvi)


gndre<-(ms_image_mask$nir-ms_image_mask$green)/(ms_image_mask$nir+ms_image_mask$green)
gndre
plot(gndre)
names(gndre) <- "gndre"

gndvi<-(ms_image_mask$nir-ms_image_mask$rededge)/(ms_image_mask$nir+ms_image_mask$rededge)

names(gndvi)<-"gndvi"


#create final raster stack for extraction

final_stack<-c(ms_image_mask,ndvi,gndre,gndvi)

plot(final_stack)

#  Extract raster statistics for each subplot ------------------------------
# Each field subplot covers many pixels. We summarize those pixels using simple
# statistics such as mean, median, quartiles, and standard deviation.
extract_one_stat <- function(r, polygons, fun, suffix, ...) {
  out <- terra::extract(
    r,
    terra::vect(polygons),
    fun = fun,
    na.rm = TRUE,
    ...
  ) %>%
    tibble::as_tibble() %>%
    dplyr::select(-ID)
  
  names(out) <- paste0(names(out), "_", suffix)
  out
}

q25_fun <- function(x, na.rm = TRUE) {
  stats::quantile(x, probs = 0.25, na.rm = na.rm, names = FALSE)
}

q75_fun <- function(x, na.rm = TRUE) {
  stats::quantile(x, probs = 0.75, na.rm = na.rm, names = FALSE)
}

extract_mean<-extract_one_stat(r=final_stack,polygons= subplots,fun=mean,suffix="mean")
extract_mean





