#'
#' 0. Overview
#'
#' This script calculates the Greenspace Availibility Index (GAVI) based on
#' three different vegetation indices: The Normalized Difference Vegetation
#' Index (NDVI), the Leaf Area Index (LAI), and a binary greenspace raster
#' (1=green; 0=non-green). The data can be retrieved from Zenodo
#' (https://zenodo.org/record/14579522).
#'

library(sf)
library(dplyr)
library(terra)
library(CGEI)


#'
#' # 1. Load data
#'

# AOI
aoi <- st_read("https://zenodo.org/records/14579522/files/01_Nbg_Stadtteile.gpkg") %>%
  summarise()

# NDVI
ndvi <- rast("https://zenodo.org/records/14579522/files/03_ndvi_10m.tif")

# LAI
lai <- rast("https://zenodo.org/records/14579522/files/03_lai_10m.tif")

# Greenspace (binary)
lulc <- rast("https://zenodo.org/records/14579522/files/03_lulc_10m.tif")


#'
#' # 2. GAVI
#'

# Combine rasters
rast_vec <- c(ndvi, lai, lulc)

# First we need to calculate lacunarity for all three rasters
# First I do it for all window sizes for a plot
lac_large <- CGEI::lacunarity(rast_vec,
                              plot = TRUE,
                              plot_path = "01_analysis/0101_data/",
                              cores = 22L,
                              progress = TRUE)

# Now only for the 5 relevant levels: 50m, 100m, 200m, 300m, 400m
lac_small <- CGEI::lacunarity(rast_vec,
                              r_vec = c(50, 100, 200, 300, 400)/10 * 2 + 1,
                              cores = 22L,
                              progress = TRUE)

# Calculate GAVI
gavi_rast <- gavi(x = rast_vec, lac_small, cores = 22, progress = TRUE)

# Save
writeRaster(gavi_rast, "01_analysis/0101_data/01_gavi.tif", overwrite = TRUE)
