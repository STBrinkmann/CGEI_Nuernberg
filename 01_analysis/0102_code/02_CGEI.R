library(dplyr)
library(sf)
library(terra)


#'
#' # 1. Load data
#'

data_repo <- "https://github.com/STBrinkmann/SpatData_Nbg/raw/refs/heads/master/01_analysis/0102_data/02_processed/"

# AOI
aoi <- read_sf(paste0(data_repo, "01_Nbg_Bezirke.gpkg")) %>%
  filter(code_bez != "97") %>% # This is not available in the GAVI
  summarise()

# 1. VGVI
vgvi <- rast("01_analysis/0101_data/01_vgvi.tif") %>%
  crop(aoi) %>%
  mask(aoi)

# 2. GAVI
gavi <- rast("01_analysis/0101_data/01_gavi.tif") %>%
  crop(aoi) %>%
  mask(aoi)

# 3. GACI
gaci <- rast("01_analysis/0101_data/01_gaci.tif") %>%
  crop(aoi) %>%
  mask(aoi)

#'
#' # 2. Calculate Composite Greenspace Exposure Index (CGEI)
#'

# First resample to same extent
vgvi <- vgvi %>%
  crop(gavi) %>%
  resample(gavi, method = "med")

gaci <- gaci %>%
  crop(gavi) %>%
  resample(gavi, method = "med")

# We will just use equal weights for the three indices. Howerver, you could also
# use different weights. E.g. put more weight on the GACI to account for the
# higher importance of access to greenspace.

cgei <- (vgvi + gavi + gaci) / 3
names(cgei) <- "cgei"

# Save
writeRaster(cgei, "01_analysis/0101_data/02_cgei.tif", overwrite = TRUE)



