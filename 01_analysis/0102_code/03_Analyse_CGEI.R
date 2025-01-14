library(tidyverse)
library(sf)
library(terra)
source("https://raw.githubusercontent.com/STBrinkmann/CGEI-Vancouver/refs/heads/main/workflow/theme_publication.R")


#'
#' # 1. Load data
#'
#' ## 1.1 Shapefiles
#'

data_repo <- "https://github.com/STBrinkmann/SpatData_Nbg/raw/refs/heads/master/01_analysis/0102_data/02_processed/"

# Bezirke
sf_bezirke <- read_sf(paste0(data_repo, "01_Nbg_Bezirke.gpkg")) %>%
  filter(code_bez != "97") %>% # This is not available in the GAVI
  group_by(code_stadtteil, code_bez, name_bez) %>%
  summarise() %>%
  ungroup()

# Stadtteile
sf_stadtteile <- read_sf(paste0(data_repo, "01_Nbg_Stadtteile.gpkg")) %>%
  st_intersection(st_geometry(sf_bezirke)) %>%
  group_by(code_stadtteil, name_stadtteil) %>%
  summarise() %>%
  ungroup()

# AOI
sf_aoi <- sf_stadtteile %>%
  st_buffer(1) %>%
  summarise()


#'
#' ## 1.2 Census
#'

# Path for the "Zensus DEU 2022"
census_path <- "~/Documents/Data/GER/Census2022/Grid"

# Only take the 100m files
census_files <- list.files(census_path, pattern = "100m", recursive = TRUE)

# Only mean age and Ausländeranteil are of interest for now
census_files <- census_files[c(3, 4)]

# Load in
census_raw <- lapply(census_files, function(x) {
  data.table::fread(file.path(census_path, x), na.strings = "-", dec = ",") %>%
    select(-any_of("werterlaeuternde_Zeichen")) %>%
    mutate(across(where(is.character), ~ gsub(",", ".", .x))) %>%
    mutate(across(-GITTER_ID_100m, as.numeric)) %>%
    mutate(across(-GITTER_ID_100m, ~tidyr::replace_na(., 0)))
})

# Left join the list
census_raw <- purrr::reduce(census_raw, dplyr::left_join, by = c("GITTER_ID_100m", "x_mp_100m", "y_mp_100m"))
census_raw_sf <- st_as_sf(census_raw, coords = c("x_mp_100m", "y_mp_100m"), crs = 3035) %>%
  st_transform(st_crs(sf_aoi))
census_raw_sf <- census_raw_sf[sf_aoi,]

# Make polygons
census_raw_sf <- census_raw_sf %>%
  st_buffer(50, endCapStyle = "SQUARE")


#'
#' ## 1.3 CGEI
#'

cgei <- rast("01_analysis/0101_data/02_cgei.tif")

# Add to the Bezirke
cgei_bezirke <- terra::extract(cgei, sf_bezirke) %>%
  group_by(ID) %>%
  reframe(CGEI = mean(cgei, na.rm = TRUE),
          VGVI = mean(vgvi, na.rm = TRUE),
          GAVI = mean(gavi, na.rm = TRUE),
          GACI = mean(gaci, na.rm = TRUE)) %>%
  mutate(code_bez = sf_bezirke$code_bez) %>%
  select(code_bez, CGEI, VGVI, GAVI, GACI)

# Also for the Stadtteile
cgei_stadtteile <- terra::extract(cgei, sf_stadtteile) %>%
  group_by(ID) %>%
  reframe(CGEI = mean(cgei, na.rm = TRUE),
          VGVI = mean(vgvi, na.rm = TRUE),
          GAVI = mean(gavi, na.rm = TRUE),
          GACI = mean(gaci, na.rm = TRUE)) %>%
  mutate(code_stadtteil = sf_stadtteile$code_stadtteil) %>%
  select(code_stadtteil, CGEI, VGVI, GAVI, GACI)

# Add to the census
cgei_census <- terra::extract(cgei, census_raw_sf) %>%
  group_by(ID) %>%
  reframe(CGEI = mean(cgei, na.rm = TRUE),
          VGVI = mean(vgvi, na.rm = TRUE),
          GAVI = mean(gavi, na.rm = TRUE),
          GACI = mean(gaci, na.rm = TRUE)) %>%
  mutate(GITTER_ID_100m = census_raw_sf$GITTER_ID_100m) %>%
  select(GITTER_ID_100m, CGEI, VGVI, GAVI, GACI)


#'
#' # 2. Analysis
#'
#' ## 2.1 Maps
#'

# Stadtteile
sf_cgei_stadtteile <- sf_stadtteile %>%
  left_join(cgei_stadtteile) %>%
  relocate(CGEI, VGVI, GAVI, GACI, .after = name_stadtteil)

ggplot(sf_cgei_stadtteile) +
  geom_sf(aes(fill = CGEI), color = "white", size = 0.2) +
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(11, "RdYlGn"),
    limits  = c(2.6, 7),
    breaks  = c(3, 5, 7),
    labels  = c("niedrig", "mittel", "hoch")
  ) +
  geom_sf(data = sf_aoi, fill = NA, color = "black", lwd = 0.4) +
  labs(
    title = "Composite Greenspace Exposure Index – Nürnberg",
    fill  = "CGEI"
  ) +
  ggthemes::theme_map() +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = "right"
  )

# Bezirke
sf_cgei_bezirke <- sf_bezirke %>%
  left_join(cgei_bezirke) %>%
  relocate(CGEI, VGVI, GAVI, GACI, .after = name_bez)

ggplot(sf_cgei_bezirke) +
  geom_sf(aes(fill = CGEI), color = "white", size = 0.2) +
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(11, "RdYlGn"),
    limits  = c(1, 9),
    breaks  = c(2, 4.5, 7),
    labels  = c("niedrig", "mittel", "hoch")
  ) +
  geom_sf(data = sf_stadtteile, fill = NA, color = "gray50", lwd = 0.2) +
  geom_sf(data = sf_aoi, fill = NA, color = "black", lwd = 0.4) +
  ggthemes::theme_map() +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = "none"
  )
ggsave("01_analysis/0101_data/03_CGEI_bezirke.png", width = 10, height = 10, dpi = 300)
ggsave("01_analysis/0101_data/03_CGEI_bezirke.svg", width = 10, height = 10)


#'
#' ## 2.2 Census
#'

census_ses <- census_raw_sf %>%
  st_drop_geometry() %>%
  left_join(cgei_census)

#'
#' ### 2.2.1 Mean age
#'

# Coarse classes
census_ses <- census_ses %>%
  mutate(age_class = case_when(
    Durchschnittsalter < 18 ~ "<18",
    Durchschnittsalter >= 18 & Durchschnittsalter < 30 ~ "18-30",
    Durchschnittsalter >= 30 & Durchschnittsalter < 40 ~ "30-40",
    Durchschnittsalter >= 40 & Durchschnittsalter < 50 ~ "40-50",
    Durchschnittsalter >= 50 & Durchschnittsalter < 60 ~ "50-60",
    Durchschnittsalter >= 60 ~ ">60",
  )) %>%
  mutate(age_class = factor(age_class, levels = c("<18", "18-30", "30-40", "40-50", "50-60", ">60")))

# Convert to shares
census_ses %>%
  group_by(age_class) %>%
  mutate(CGEI_class = mean(CGEI)) %>%
  ungroup() %>%
  ggplot(aes(x = age_class, y = CGEI, fill = CGEI_class)) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(11, "RdYlGn"),
    limits  = c(2, 7.5),
    breaks  = c(2.5, 4.75, 7),
    labels  = c("niedrig", "mittel", "hoch")
  ) +
  labs(x = "Durchschnittsalter", y = "CGEI") +
  theme_Publication() +
  theme(legend.position = "none")
ggsave("01_analysis/0101_data/03_CGEI_age.png", width = 10, height = 6, dpi = 300)
ggsave("01_analysis/0101_data/03_CGEI_age.svg", width = 10, height = 6)


#'
#' ### 2.2.2 Ausländeranteil
#'

auslaender_jenks <- classInt::classIntervals(census_ses$AnteilAuslaender, n = 6, style = "fisher", intervalClosure = "right")$brks
census_ses <- census_ses %>%
  mutate(auslaender_class = cut(AnteilAuslaender, breaks = auslaender_jenks, include.lowest = TRUE)) %>%
  mutate(auslaender_class = factor(auslaender_class, levels = levels(auslaender_class),
                                   labels = paste0("bis ", round(auslaender_jenks[-1]), "%")))

census_ses %>%
  group_by(auslaender_class) %>%
  mutate(CGEI_class = mean(CGEI)) %>%
  ungroup() %>%
  ggplot(aes(x = auslaender_class, y = CGEI, fill = CGEI_class)) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(11, "RdYlGn"),
    limits  = c(2, 7.5),
    breaks  = c(2.5, 4.75, 7),
    labels  = c("niedrig", "mittel", "hoch")
  ) +
  labs(x = "Ausländeranteil", y = "CGEI") +
  theme_Publication() +
  theme(legend.position = "none")
ggsave("01_analysis/0101_data/03_CGEI_auslaender.png", width = 10, height = 6, dpi = 300)
ggsave("01_analysis/0101_data/03_CGEI_auslaender.svg", width = 10, height = 6)
