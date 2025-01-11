#'
#' 0. Overview
#'
#' This script calculates the Viewshed Greenness Visibilty Index (VGVI) based on
#' a Digital Surface Model (DSM), a Digital Terrain Model (DTM) and a binary
#' greenspace raster (1=green; 0=non-green). The data can be retrieved from
#' Zenodo (https://zenodo.org/record/14633167).
#'

library(sf)
library(dplyr)
library(terra)
library(CGEI)


#'
#' # 1. Load data
#'

# AOI
aoi <- st_read("https://zenodo.org/records/14633167/files/01_Nbg_Stadtteile.gpkg") %>%
  summarise()

# DSM
dsm <- rast("https://zenodo.org/records/14633167/files/03_dsm_1m.tif")

# DTM
dtm <- rast("https://zenodo.org/records/14633167/files/03_dtm_1m.tif")

# Greenspace (binary)
lulc <- rast("https://zenodo.org/records/14633167/files/03_ndvi_01_1m.tif")


#'
#' # 2. VGVI
#'

# First we need all valid cells in a 5x5m grid. These are where the DSM is not
# higher than 1.8m + DTM
dsm5 <- aggregate(dsm, 5)
dtm5 <- aggregate(dtm, 5)

obs <- dsm5 <= (dtm5 + 2.2)
obs <- obs %>%
  crop(aoi) %>%
  mask(aoi)

obs_vals <- values(obs)
obs_vals <- which(as.vector(obs_vals))

# Convert the cell coordinates as an SF
obs_sf <- terra::xyFromCell(obs, obs_vals)
obs_sf <- st_as_sf(as_tibble(obs_sf), coords = c("x", "y"), crs = crs(dsm)) %>%
  mutate(id = 1:n()) %>%
  relocate(id)

# Calculate VGVI
vgvi_sf <- vgvi(observer = obs_sf,
                dsm_rast = dsm, dtm_rast = dtm, greenspace_rast = lulc,
                max_distance = 500, observer_height = 2.2, # I set it to 2.2m so that we can also look above cars on the road
                m = 1, b = 3, mode = "exponential",
                cores = 22, progress = TRUE)

# Now rasterize the VGVI
vgvi_sf <- vgvi_sf[aoi,]
vgvi_rast <- sf_interpolat_IDW(observer = vgvi_sf,
                               v = "VGVI",
                               aoi = aoi,
                               raster_res = 5,
                               n = 10, beta = 2, max_distance = 500,
                               na_only = TRUE,
                               cores = 22, progress = TRUE)

vgvi_rast_10 <- aggregate(vgvi_rast, 2)

vgvi_rast_10 <- vgvi_rast_10 %>%
  crop(aoi) %>%
  mask(aoi)

vgvi_rast_10 <- CGEI:::reclassify_jenks(vgvi_rast_10, 9)
vgvi_rast_10 <- as.int(vgvi_rast_10)

# Save
writeRaster(vgvi_rast_10, "01_analysis/0101_data/01_vgvi.tif", overwrite = TRUE)
