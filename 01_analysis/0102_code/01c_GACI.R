library(sf)
library(dplyr)
library(osmdata)
library(osrm)
library(terra)
library(pbmcapply)

# Helper function
st_multipoint_to_point <- function(x) {
  x.points <- x[which(sf::st_geometry_type(x) == "POINT"), ]
  x.multipoints <- x[which(sf::st_geometry_type(x) == "MULTIPOINT"), ]

  for (i in 1:nrow(x.multipoints)) {
    x.points <- x.points %>%
      dplyr::add_row(sf::st_cast(x.multipoints[i, ], "POINT"))
  }

  return(x.points)
}


#'
#' # 1. Load data
#'

data_repo <- "https://github.com/STBrinkmann/SpatData_Nbg/raw/refs/heads/master/01_analysis/0102_data/02_processed/"

# AOI
aoi <- read_sf(paste0(data_repo, "01_Nbg_Stadtteile.gpkg")) %>%
  st_buffer(0.1) %>%
  summarise()

# Water mask
water_mask <- read_sf(paste0(data_repo, "02_Nbg_Water.gpkg"))

# Parks (>1ha)
osm_parks <- read_sf(paste0(data_repo, "02_Nbg_Greenspace.gpkg")) %>%
  mutate(name = stringr::str_pad(1:n(), 3, pad = "0")) %>%
  mutate(size_ha = st_area(geom) %>%
           units::set_units(value = ha) %>%
           as.numeric()) %>%
  filter(size_ha > 1) %>%
  mutate(size_log = log(size_ha)) %>%
  relocate(geom, .after = last_col())


#'
#' # 2. Estimate park access points
#'
#' Parks can be accessed on a road network through access points. These are the
#' intersections of the park boundary with the road network.
#'

# Add a 5km buffer around the AOI to avoid edge effects
aoi_bbox <- aoi %>%
  st_buffer(5000) %>%
  st_transform(4326) %>%
  st_bbox()

# 1. Download OSM road network
osm_roads_raw <- opq(bbox = aoi_bbox) %>%
  add_osm_feature(key = "highway") %>%
  osmdata_sf(quiet = FALSE) %>%
  osm_poly2line()

# Remove features that are not made for walking (e.g. motorway, ...)
osm_roads <- osm_roads_raw$osm_lines %>%
  select(highway) %>%
  filter(!highway %in% c("motorway", "motorway_link", "trunk", "trunk_link", "raceway")) %>%
  st_transform(st_crs(aoi)) %>%
  st_intersection(st_buffer(aoi , 5000))

# Break lines in segments to increase accuracy
cores <- 12
# ---- WINDWOS ----
if (Sys.info()[["sysname"]] == "Windows") {
  cl <- parallel::makeCluster(cores)
  osm_roads <- suppressWarnings(split(osm_roads, seq(from = 1, to = nrow(osm_roads), by = 200)))
  osm_roads <- parallel::parLapply(cl, osm_roads, fun = function(x){
    nngeo::st_segments(sf::st_cast(x, "LINESTRING"), progress = FALSE)
  })
  parallel::stopCluster(cl)
} else { # ---- Linux and macOS ----
  osm_roads <- suppressWarnings(split(osm_roads, seq(from = 1, to = nrow(osm_roads), by = 200))) %>%
    parallel::mclapply(function(x){
      nngeo::st_segments(sf::st_cast(x, "LINESTRING"), progress = FALSE)
    },
    mc.cores = cores, mc.preschedule = TRUE)
}

osm_roads <- st_as_sf(dplyr::as_tibble(data.table::rbindlist(osm_roads)))
osm_roads <- st_set_geometry(osm_roads, "geom")

# Round coordinates to 0 digits.
st_geometry(osm_roads) <- st_geometry(osm_roads) %>%
  lapply(function(x) round(x, 0)) %>%
  st_sfc(crs = st_crs(osm_roads))

# 2. Park access points (park centroids + streets that lead into park)
osm_park_centroids <- osm_parks %>%
  st_centroid() %>%
  select(name, size_ha, size_log)

osm_park_road_interesect <- osm_parks %>%
  st_boundary() %>%
  st_intersection(osm_roads) %>%
  select(name, size_ha, size_log) %>%
  st_multipoint_to_point()

park_accesspoints <- rbind(osm_park_centroids, osm_park_road_interesect)


#'
#' # 3. Routing
#'

# First derive observer locations for a regular point grid each 10m
observer <- st_make_grid(aoi, cellsize = 50, what = "centers") %>%
  st_as_sf() %>%
  st_set_geometry("geom")
observer <- observer[aoi,]
observer <- observer %>%
  mutate(id = 1:n()) %>%
  relocate(id, .after = last_col())

# Observers in parks always have good access
obs_in_parks <- observer[osm_parks, ]

# So remove them from the observer list
observer <- observer[-obs_in_parks$id, ]

# Also remove those inside water
observer <- observer[-(observer[water_mask, ]$id), ]


# I am using the Open Source Routing Machine (OSRM) for routing. The OSRM server
# is running in a Docker container. The following code downloads the road network.
# Make sure that Docker is installed (https://docs.docker.com/engine/install/ubuntu/)
if(!file.exists("01_analysis/0101_data/osrm/mittelfranken-latest.osm.pbf")){
  dir.create("01_analysis/0101_data/osrm/", showWarnings = FALSE)
  # Change timeout
  default_timeout <- getOption("timeout")
  options(timeout = max(600, default_timeout))
  download.file(url = "https://download.geofabrik.de/europe/germany/bayern/mittelfranken-latest.osm.pbf",
                destfile = "01_analysis/0101_data/osrm/mittelfranken-latest.osm.pbf", mode="wb")
  options(timeout = default_timeout)
}

# Once the download is complete, open a terminal/bash and navigate to the
# "01_analysis/0101_data/" folder.
cat("cd", file.path(getwd(), "01_analysis/0101_data/osrm/"))

# Run these lines to start a local Docker OSRM instance:

# sudo docker run -t -v "${PWD}:/data" osrm/osrm-backend osrm-extract -p /opt/foot.lua /data/mittelfranken-latest.osm.pbf
# sudo docker run -t -v "${PWD}:/data" osrm/osrm-backend osrm-partition /data/mittelfranken-latest.osrm
# sudo docker run -t -v "${PWD}:/data" osrm/osrm-backend osrm-customize /data/mittelfranken-latest.osrm
# sudo docker run -t -i -p 5000:5000 -v "${PWD}:/data" osrm/osrm-backend osrm-routed --max-table-size=50000 --algorithm mld /data/mittelfranken-latest.osrm

# Function that calculates the walking distance to the big parks for a distinct observer
observer_park_duration <- function(j){
  # Get knn with k = sqrt of number of park access points
  suppressMessages({
    knn_access <- nngeo::st_nn(observer[j,], park_accesspoints,
                               k = 2*round(sqrt(nrow(park_accesspoints))),
                               progress = FALSE) %>%
      unlist()
  })
  # knn_access <- 1:nrow(park_accesspoints)

  # Calculate duration using local docker instance
  duration_table <- osrmTable(src = observer[j,],
                              dst = park_accesspoints[knn_access,],
                              osrm.server = "http://0.0.0.0:5000/", osrm.profile = "foot")

  # Get the smallest distance of each unique park
  out <- tibble(
    id = observer[j,]$id,
    dist = as.numeric(duration_table$durations),
    size_ha = park_accesspoints[knn_access,]$size_ha,
    size_log = park_accesspoints[knn_access,]$size_log,
    name = park_accesspoints[knn_access,]$name
  ) %>%
    group_by(name) %>%
    filter(dist == min(dist)) %>%
    distinct() %>%
    ungroup() %>%
    mutate(dist_01 = dist / mean(dist) * 100,
           size_01 = size_log / mean(size_log) * 100) %>%
    mutate(dist_01 = 1 - dist_01 / max(dist_01)) %>%
    mutate(size_01 = size_01 / max(size_01)) %>%
    mutate(w = (1.91/0.85*dist_01 + size_01) / 2) %>%
    filter(w == max(w)) %>%
    select(id, name, dist, size_ha, w)

  #invisible(gc())
  return(out)
}

# Split observer in list of 5000 ids
obs_seq <- split(1:nrow(observer), ceiling(seq_along(1:nrow(observer))/250))

# Output table that will be loaded with data
tibble(
  id = as.numeric(),
  name = as.character(),
  dist = as.numeric(),
  size_ha = as.numeric(),
  w = as.numeric()
) %>%
  readr::write_csv("01_analysis/0101_data/osrm/access_table.csv", col_names = TRUE)

# Parallel calculate distance to nearest big park
pb = pbmcapply::progressBar(min = 0, max = length(obs_seq)); stepi = 1
for(this_seq in obs_seq){
  # Get distance to closest parks
  dist_table <- parallel::mclapply(this_seq, observer_park_duration,
                                   mc.cores = 18) %>%
    do.call(rbind, .)

  # Append to table
  readr::write_csv(dist_table, "01_analysis/0101_data/osrm/access_table.csv",
                   append = TRUE)

  setTxtProgressBar(pb, stepi); stepi = stepi + 1
}
close(pb)

# Load all distance data
dist_table <- readr::read_csv("01_analysis/0101_data/osrm/access_table.csv")

# Merge with observer locations
dist_sf <- inner_join(observer, dist_table)

# End the OSRM server (!) and clean up
unlink("01_analysis/0101_data/osrm/", recursive = TRUE)


#'
#' # 4. Convert to Jenks-raster
#'

# First bring it to a 10m raster
gaci_raw <- CGEI::sf_interpolat_IDW(observer = dist_sf,
                                    v = "w",
                                    aoi = aoi,
                                    max_distance = Inf,
                                    n = 20,
                                    raster_res = 10,
                                    na_only = TRUE,
                                    cores = 22,
                                    progress = TRUE)

park_cells <- terra::extract(gaci_raw, osm_parks, cells = TRUE, ID = FALSE)$cell
gaci_raw[park_cells] <- NA

# And classify using Jenks
gaci <- CGEI:::reclassify_jenks(gaci_raw, n_classes = 9)
gaci[park_cells] <- 9

gaci <- gaci %>%
  crop(aoi) %>%
  mask(aoi)

# Save
writeRaster(gaci, "01_analysis/0101_data/01_gaci.tif", overwrite = TRUE)

















