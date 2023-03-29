# *------------------------------------------------------------------
# | PROGRAM NAME:
# | FILE NAME: .R
# | DATE:
# | CREATED BY:  Jim Stagge
# *----------------------------------------------------------------
# | PURPOSE:
# |
# |
# *------------------------------------------------------------------

cat("Starting 05_model_3_regions")

###########################################################################
## Set the Paths
###########################################################################
require(here)

### Path for Data and Output
data_path <- file.path(here(), "data")
#output_path <- "/fs/ess/PAS1921"
output_path <- file.path(here(), "output")

### Set up output folders
write_output_path <- output_path
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

### Set up figure folder
write_figures_path <- file.path(output_path, "figures")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)

###########################################################################
###  Load functions
###########################################################################
require(tidyverse)
require(tictoc)
require(viridis)

### For Dates
require(lubridate)

### For mapping
require(rnaturalearth)
#dyn.load("/usr/local/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
#dyn.load("/usr/local/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
require(sf)
#library(maps)


require(rnoaa)
require(tidyverse)
library(sf)
library(stars)
library(rnaturalearth)
library(rnaturalearthhires)

require(scales)
#require(lfstat)

require(chronosphere)

#require(ggthemes)
require(raster)
select <- dplyr::select

require(mgcv)
require(gratia)

###########################################################################
###  Read in data
###########################################################################
gauge_path <- file.path(output_path, "gauge_data")

djf_df <- readRDS(file.path(gauge_path, "djf_df.rds"))

### Read in station info
stations <- readRDS(file.path(gauge_path, "all_stations.rds"))
gauges_sf <- readRDS(file.path(gauge_path, "gauges_sf.rds"))
gauges_sf_proj <- readRDS(file.path(gauge_path, "gauges_sf_proj.rds"))


### Read in plot grid
plot_grid_sf <- readRDS(file.path(gauge_path, "plot_grid_sf.rds"))

###########################################################################
###  Create folder
###########################################################################
model_folder <- file.path(output_path, "model_3_region")
dir.create(model_folder, recursive=TRUE, showWarnings = FALSE)

write_figures_path <- file.path(write_figures_path, "model_3_region")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)

###########################################################################
###  Read in models 1 and 2
###########################################################################
model_climatol <- readRDS(file.path(output_path, "model_1_climatol/model_climatol.rds"))
model_enso <- readRDS(file.path(output_path, "model_2_enso/model_enso.rds"))

###########################################################################
###  Clean up a potential problem site
###########################################################################
djf_df <- djf_df %>%
	filter(id != "MXN00008118") %>%
	mutate(id = as.factor(id))

###########################################################################
###  Download background data
###########################################################################
### Download data for plot
usa_states <- ne_states(country = 'United States of America', returnclass = 'sf')
ca_states <- ne_states(country = 'Canada', returnclass = 'sf')
mx_states <- ne_states(country = 'Mexico', returnclass = 'sf')

lakes <- ne_download(scale = 50, type = 'lakes', category = 'physical', returnclass = 'sf')
rivers <- ne_download(scale = 10, type = 'rivers_lake_centerlines', category = 'physical', returnclass = 'sf')

#background_orig <- ne_download(scale = 50, type = 'HYP_50M_SR_W', category = 'raster')
#background_orig <- ne_download(scale = 10, type = 'HYP_LR_SR_W_DR', category = 'raster')

### Convert background for easy plotting
#bb = st_bbox(c(xmin = -115, ymin = 25, xmax = -100, ymax = 50), crs = st_crs(background_orig))
#background <- st_as_stars(background_orig, proxy= FALSE) %>% st_crop(bb) %>%st_rgb(3)

##########################################################################
###  Create regions
###########################################################################

### You could raise this further north to cut out the Uintas, which are unique (E-W)
### Define the Mexico (Sierra Madre Occidental) polygon
x_coord <- c(-108.5,  -112,  -106.8, -104.3, -108.5)
y_coord <- c(26, 30.2, 30.5, 26.2, 26)

mx_pol = st_sfc(st_polygon(list(cbind(x_coord, y_coord))))
st_crs(mx_pol) = 4326

### Define the Utah Wyoming (Wasatch, Wyoming, Teton, Wind River Ranges)
#x_coord <- c(-112.3,  -112.3,  -108, -108, -112.3)
#y_coord <- c(39.8, 45.5, 45.5, 39.8, 39.8)

### You could raise this further north to cut out the Uintas, which are unique (E-W)
x_coord <- c(-112.3,  -112.3,  -108, -108, -112.3)
y_coord <- c(41.2, 45.5, 45.5, 41.2, 41.2)

wyo_pol = st_sfc(st_polygon(list(cbind(x_coord, y_coord))))
st_crs(wyo_pol) = 4326

### Create north and south polygons for plotting
n_pol <- st_sfc(st_polygon(list(cbind(c(-116,  -116,  -99, -99, -116), c(37, 50, 50, 37, 37)))))
st_crs(n_pol) = 4326

n_pol <- st_sfc(st_polygon(list(cbind(c(-116,  seq(-116, -99, length.out = 100), seq(-99, -116, length.out = 100)), c(37, rep(50,100), rep(37,100))))))
st_crs(n_pol) = 4326

s_pol <- st_sfc(st_polygon(list(cbind(c(-116,  -116,  -99, -99, -116), c(25, 37, 37, 25, 25)))))
st_crs(s_pol) = 4326

s_pol <- st_sfc(st_polygon(list(cbind(c(-116,  seq(-116, -99, length.out = 100), seq(-99, -116, length.out = 100)), c(25, rep(37,100), rep(25,100))))))
st_crs(s_pol) = 4326

plot_df <- plot_grid_sf

### Test plot
p <- ggplot(data = plot_df) %>%
	#+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	#+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = elevation)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = elevation), colour = NA) %>%
	+ geom_sf(data = n_pol, fill = "white", alpha = 0.2, size = 1, colour = "orange") %>%
	+ geom_sf(data = s_pol, fill = "white", alpha = 0.2, size = 1, colour = "red") %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_pol, fill = "red", alpha = 0.5) %>%
	+ geom_sf(data = wyo_pol, fill = "orange", alpha = 0.5) %>%
	+ scale_fill_viridis(name = "Elevation (m)") %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#p

### Save plot
ggsave(file.path(write_figures_path, "region_plot.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "region_plot.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "region_plot.svg"), p,  width = 4, height = 6.5)

### Extract regions
mx_test <- c(st_contains(mx_pol , gauges_sf%>% st_transform(st_crs(mx_pol))))[[1]]
wyo_test <- c(st_contains(wyo_pol, gauges_sf%>% st_transform(st_crs(mx_pol))))[[1]]

### Test plot
p <- ggplot(data = plot_df) %>%
	#+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	#+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	#+ geom_raster(alpha = 0.75, aes(x = x_proj, y=y_proj, fill = elevation), alpha = 0.2) %>%
	+ geom_sf(alpha = 0.95, aes(fill = elevation), colour = NA) %>%
	+ geom_sf(data = n_pol, fill = "white", alpha = 0.1, size = 1, colour = "orange") %>%
	+ geom_sf(data = s_pol, fill = "white", alpha = 0.1, size = 1, colour = "red") %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_pol, fill = "red", alpha = 0.7) %>%
	+ geom_sf(data = wyo_pol, fill = "red", alpha = 0.7) %>%
	+ geom_sf(data = gauges_sf[mx_test,], aes(colour = elevation), size = 0.6) %>%
	+ geom_sf(data = gauges_sf[wyo_test,], aes(colour = elevation), size = 0.6) %>%
	+ scale_fill_viridis(name = "Elevation (m)", limits = c(0,3400), oob = squish) %>%
	+ scale_colour_viridis(name = "Elevation (m)", limits = c(0,3400), oob = squish) %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#p

### Save plot
ggsave(file.path(write_figures_path, "region_plot_elev.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "region_plot_elev.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "region_plot_elev.svg"), p,  width = 4, height = 6.5)


##########################################################################
###  Create column for regions and save to file
###########################################################################

### Create column for regions
djf_df <- djf_df %>%
	mutate(region = case_when(
		id %in% gauges_sf$id[wyo_test] ~ "wyo",
		id %in% gauges_sf$id[mx_test] ~ "mx",
		latitude >= 37 ~ "north",
		latitude < 37 ~ "south",
		TRUE  ~ NA_character_)
	) %>%
	mutate(region = factor(region))

### Save the djf dataframe with new region column
saveRDS(djf_df, file = file.path(write_output_path, "djf_df.rds"))


### Also run for plot_grid
#plot_grid_sf = st_as_sf(plot_grid, coords = c("lon", "lat"), crs = 4326, agr = "constant")

mx_test <- c(st_contains(mx_pol %>% st_transform(st_crs(plot_grid_sf)) , plot_grid_sf))[[1]]
mx_test <- unique(c(mx_test, c(st_touches(mx_pol %>% st_transform(st_crs(plot_grid_sf)) , plot_grid_sf))[[1]] ))
wyo_test <- c(st_intersects(wyo_pol %>% st_transform(st_crs(plot_grid_sf)) ,  plot_grid_sf))[[1]]
wyo_test <- unique(c(wyo_test, c(st_touches(wyo_pol %>% st_transform(st_crs(plot_grid_sf)) , plot_grid_sf))[[1]] ))

### Create column for regions
plot_grid_sf <- plot_grid_sf %>%
	mutate(region = case_when(
		point %in% plot_grid_sf$point[wyo_test] ~ "wyo",
		point %in% plot_grid_sf$point[mx_test] ~ "mx",
		lat >= 37 ~ "north",
		lat < 37 ~ "south",
		TRUE  ~ NA_character_)
	) %>%
	mutate(region = factor(region))

ggplot(plot_grid_sf, aes(fill = region)) + geom_sf() + 	geom_sf(data = mx_pol, fill = "red", alpha = 0.5) + geom_sf(data = wyo_pol, fill = "orange", alpha = 0.5)
#saveRDS(djf_df, file = file.path(write_output_path, "plot_grid.rds"))

##########################################################################
###  Create a transect plot
###########################################################################
### Wyoming transect
wyo_profile  <- st_sfc(st_linestring(rbind(c(-113,43.5),c(-110,43.5),c(-108, 43.5)))) %>%
	st_set_crs(4326)

### Mexico transect
mx_profile  <- st_sfc(st_linestring(rbind(c(-110,28),c(-108,28),c(-105, 28)))) %>%
	st_set_crs(4326)


### Test plot
p <- ggplot(data = plot_df) %>%
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = elevation)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = elevation), colour = NA) %>%
	#+ geom_sf(data = n_pol, fill = "white", alpha = 0.2, size = 1, colour = "orange") %>%
	#+ geom_sf(data = s_pol, fill = "white", alpha = 0.2, size = 1, colour = "red") %>%
	+ geom_sf(data = lakes, fill = "#fb8072", alpha = 0.4) %>%
	+ geom_sf(data = rivers, colour = "#fb8072", alpha = 0.4) %>%	#"#4a80f5"
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_pol, fill = NA, colour = "white", alpha = 1, size = 1.1, linetype = "twodash") %>%
	+ geom_sf(data = wyo_pol, fill = NA, colour = "white",  alpha = 1, size = 1.1, linetype = "twodash") %>%
	+ geom_sf(data = wyo_profile, colour = "#e41a1c", size = 0.9) %>%
	+ geom_sf(data = mx_profile, colour = "#e41a1c", size = 0.9) %>%
	+ scale_fill_viridis(name = "Elevation (m)") %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#p


### Save plot
ggsave(file.path(write_figures_path, "region_plot_transect.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "region_plot_transect.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "region_plot_transect.svg"), p,  width = 4, height = 6.5)


### Test plot
p <- ggplot(data = plot_df) %>%
	#+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	#+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = elevation)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = elevation), colour = NA) %>%
	#+ geom_sf(data = n_pol, fill = "white", alpha = 0.2, size = 1, colour = "orange") %>%
	#+ geom_sf(data = s_pol, fill = "white", alpha = 0.2, size = 1, colour = "red") %>%
#	+ geom_sf(data = rivers, colour = "#4a80f5", alpha = 0.15) %>%
#	+ geom_sf(data = lakes, fill = "#9bbff4", alpha = 0.15) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_pol, fill = "red", alpha = 0.5) %>%
	+ geom_sf(data = wyo_pol, fill = "orange", alpha = 0.5) %>%
	+ geom_sf(data = wyo_profile, colour = "grey", size = 1.5) %>%
	+ geom_sf(data = mx_profile, colour = "grey", size = 1.5) %>%
	+ geom_sf(data = gauges_sf[mx_test,], aes(colour = elevation), size = 0.6) %>%
	+ geom_sf(data = gauges_sf[wyo_test,], aes(colour = elevation), size = 0.6) %>%	+ scale_fill_viridis(name = "Elevation (m)") %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#p


### Save plot
ggsave(file.path(write_figures_path, "region_plot_transect_elev.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "region_plot_transect_elev.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "region_plot_transect_elev.svg"), p,  width = 4, height = 6.5)


### Test plot
p <- ggplot(data = plot_df) %>%
	#+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	#+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	#+ geom_raster(alpha = 0.9, aes(x = x_proj, y=y_proj, fill = elevation)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = elevation), colour = NA) %>%
	#+ geom_sf(data = n_pol, fill = "white", alpha = 0.2, size = 1, colour = "orange") %>%
	#+ geom_sf(data = s_pol, fill = "white", alpha = 0.2, size = 1, colour = "red") %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = wyo_pol, fill = "white", colour = NA, alpha = 0.4) %>%
	+ geom_sf(data = wyo_pol, colour = "red", fill = NA, size = 0.5) %>%
	#+ geom_sf(data = gauges_sf[mx_test,], aes(colour = elevation), size = 0.6) %>%
	+ geom_sf(data = gauges_sf[wyo_test,], aes(colour = elevation), size = 0.6) %>%
	+ geom_sf(data = wyo_profile, colour = "black", size = 1.5) %>%
	+ geom_sf(data = mx_profile, colour = "black", size = 1.5) %>%
	+ scale_fill_viridis(name = "Elevation (m)", limits = c(0,3400)) %>%
	+ scale_colour_viridis(name = "Elevation (m)", limits = c(0,3400)) %>%
	+ coord_sf(xlim = c(-16e5, -7e5), ylim = c(19e5, 28e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 2)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#p

### Save plot
ggsave(file.path(write_figures_path, "wyo_region_plot.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "wyo_region_plot.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "wyo_region_plot.svg"), p,  width = 4, height = 6.5)


### Test plot
p <- ggplot(data = plot_df) %>%
	#+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	#+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	#+ geom_raster(alpha = 0.9, aes(x = x_proj, y=y_proj, fill = elevation)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = elevation), colour = NA) %>%
	#+ geom_sf(data = n_pol, fill = "white", alpha = 0.2, size = 1, colour = "orange") %>%
	#+ geom_sf(data = s_pol, fill = "white", alpha = 0.2, size = 1, colour = "red") %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_pol, fill = "white", colour = NA, alpha = 0.4) %>%
	+ geom_sf(data = mx_pol, colour = "red", fill = NA, size = 0.5) %>%
	#+ geom_sf(data = gauges_sf[mx_test,], aes(colour = elevation), size = 0.6) %>%
	+ geom_sf(data = gauges_sf[mx_test,], aes(colour = elevation), size = 0.6) %>%
	+ geom_sf(data = wyo_profile, colour = "black", size = 1.5) %>%
	+ geom_sf(data = mx_profile, colour = "black", size = 1.5) %>%
	+ scale_fill_viridis(name = "Elevation (m)", limits = c(0,3400)) %>%
	+ scale_colour_viridis(name = "Elevation (m)", limits = c(0,3400)) %>%
	+ coord_sf(xlim =  c(-17.5e5, -7e5), ylim = c(2e5, 10.5e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 2)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#p

### Save plot
ggsave(file.path(write_figures_path, "mx_region_plot.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "mx_region_plot.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "mx_region_plot.svg"), p,  width = 4, height = 6.5)




##########################################################################
###  Plot Elevation Profiles
###########################################################################
plot_wyo_profile <- plot_grid_sf %>%
	filter(lat >= 43.4 & lat <= 43.5 & lon >= -113 & lon <= -108) %>%
	as.data.frame()

p <- ggplot(plot_wyo_profile, aes(x=lon, y=elevation)) %>%
	+ geom_line() %>%
	+ scale_x_continuous(name = "Longitude") %>%
	+ scale_y_continuous(name = "Elevation (m)") %>%
	+ theme_bw(9)

### Save plot
ggsave(file.path(write_figures_path, "wyo_transect_elevation.png"), p,  width = 5, height = 3, dpi = 600)
ggsave(file.path(write_figures_path, "wyo_transect_elevation.pdf"), p,  width = 5, height = 3)
ggsave(file.path(write_figures_path, "wyo_transect_elevation.svg"), p,  width = 5, height = 3)

plot_mx_profile <- plot_grid_sf %>%
	filter(lat >= 27.9500 & lat <= 28.05001 & lon >= -110 & lon <= -105) %>%
	as.data.frame()

p <- ggplot(plot_mx_profile, aes(x=lon, y=elevation)) %>%
	+ geom_line() %>%
	+ scale_x_continuous(name = "Longitude") %>%
	+ scale_y_continuous(name = "Elevation (m)") %>%
	+ theme_bw(9)

### Save plot
ggsave(file.path(write_figures_path, "mx_transect_elevation.png"), p,  width = 5, height = 3, dpi = 600)
ggsave(file.path(write_figures_path, "mx_transect_elevation.pdf"), p,  width = 5, height = 3)
ggsave(file.path(write_figures_path, "mx_transect_elevation.svg"), p,  width = 5, height = 3)


##########################################################################
###  Fit Wyoming model
###########################################################################
djf_df <- djf_df %>%
	drop_na(djf_mm_daily, x_proj, y_proj, elevation, enso)

### Regional model for WYO
#submodel_wyo <-  bam(djf_mm_daily  ~ te(x_proj, y_proj, enso, d=c(2,1), bs = c("tp", "tp")),
#						  data = djf_df %>% filter(region == "wyo") ,
#						 family =  tw())

#BIC(submodel_wyo)
#AIC(submodel_wyo)

#appraise(submodel_wyo)

submodel_wyo <-  bam(djf_mm_daily  ~ ti(x_proj, y_proj, d=c(2), bs = c("tp")) +
   ti(x_proj, y_proj, enso, d=c(2,1), bs = c("tp", "tp")) +
   ti(enso, bs = "tp") + s(elevation, bs = "tp")  + s(year, bs = "tp", k = 4),
						 						  data = djf_df %>% filter(region == "wyo") ,
						 						 family =  tw())

BIC(submodel_wyo)
AIC(submodel_wyo)
summary(submodel_wyo)


### Save the model objects
saveRDS(submodel_wyo, file = file.path(model_folder, "submodel_wyo.rds"))

### Save some residual information
resid_model <- resid(submodel_wyo)

resid_stats <- data.frame(mean_resid = mean(resid_model, na.rm=TRUE),
	median_resid = median(resid_model, na.rm=TRUE),
	mae = mean(abs(resid_model), na.rm=TRUE),
	rmse = sqrt(mean((resid_model)^2)),
	aic = AIC(submodel_wyo),
	bic = BIC(submodel_wyo),
	r_square = summary(submodel_wyo)$r.sq
	)

write_csv(resid_stats, file.path(model_folder, "submodel_wyo_resid_stats.csv"))


p <- appraise(submodel_wyo)
#p

### Save plot
ggsave(file.path(write_figures_path, "submodel_wyo_appraise.png"), p,  width = 7, height = 7, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_wyo_appraise.pdf"), p,  width = 7, height = 7)
ggsave(file.path(write_figures_path, "submodel_wyo_appraise.svg"), p,  width = 7, height = 7)

p <- draw(submodel_wyo)
p

### Save plot
ggsave(file.path(write_figures_path, "submodel_wyo_draw.png"), p,  width = 7, height = 7, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_wyo_draw.pdf"), p,  width = 7, height = 7)
ggsave(file.path(write_figures_path, "submodel_wyo_draw.svg"), p,  width = 7, height = 7)


##########################################################################
###  Plot Wyoming model
###########################################################################
plot_wyo_profile <- plot_grid_sf %>%
	filter(lat >= 43.4 & lat <= 43.5 & lon >= -113 & lon <= -108) %>%
	as.data.frame() %>%
	mutate(year = 1975)

### Add in ENSO terms
plot_wyo_profile <- expand_grid(plot_wyo_profile, data.frame(enso = seq(-2,3,0.5)))

plot_wyo_profile$pred <- predict(submodel_wyo, newdata = plot_wyo_profile, type = "response")

plot_wyo_profile <- plot_wyo_profile %>%
	mutate(pred_mm =  pred* 90)

p <- ggplot(plot_wyo_profile, aes(x=lon, y=pred, colour = enso, group = enso)) %>%
	+ geom_line(aes(y=pred_mm, colour = enso, group = enso)) %>%
	+ scale_colour_distiller(name = "ENSO", type = "div", palette = "PuOr", direction = -1, limits = c(-3,3)) %>%
	+ scale_x_continuous(name = "Longitude") %>%
	+ scale_y_continuous(name = "Predicted DJF \nPrecipitation (mm)") %>%
	+ theme_bw(9) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


### Save plot
ggsave(file.path(write_figures_path, "submodel_wyo_transect.png"), p,  width = 5, height = 4.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_wyo_transect.pdf"), p,  width = 5, height = 4.5)
ggsave(file.path(write_figures_path, "submodel_wyo_transect.svg"), p,  width = 5, height = 4.5)


### I could put elevation on second axis

p <- ggplot(plot_wyo_profile, aes(x=lon)) %>%
	+ geom_line(aes(y=pred_mm, colour = enso, group = enso)) %>%
	+ geom_line(data = plot_wyo_profile %>% filter(enso == 0), aes(y = elevation/8+1000), colour = "grey50", size = 0.25) %>%
#	+ scale_colour_gradient2(name = "ENSO") %>%
	+ scale_colour_distiller(name = "ENSO", type = "div", palette = "PuOr", direction = -1, limits = c(-3,3)) %>%
	+ scale_x_continuous(name = "Longitude") %>%
	+ scale_y_continuous(name = "Predicted DJF \nPrecipitation (mm)" , sec.axis = sec_axis(~.*8 - 1000, name="Elevation (m)")) %>%
	+ coord_cartesian() %>%
	+ theme_bw(9) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


	plot_wyo_profile <- plot_wyo_profile %>%
		mutate(elev_adj = elevation/5 + 400/5)


	p <- ggplot(plot_wyo_profile, aes(x=lon)) %>%
		+ geom_line(aes(y=pred_mm, colour = enso, group = enso)) %>%
		+ geom_line(data = plot_wyo_profile %>% filter(enso == 0), aes(y = elev_adj), colour = "grey50", size = 0.25) %>%
	#	+ scale_colour_gradient2(name = "ENSO") %>%
		+ scale_colour_distiller(name = "ENSO", type = "div", palette = "PuOr", direction = -1, limits = c(-3,3)) %>%
		+ scale_x_continuous(name = "Longitude") %>%
		+ scale_y_continuous(name = "Predicted DJF \nPrecipitation (mm)" , sec.axis = sec_axis(~.*5-400, name="Elevation (m)", breaks = seq(0,4000, by = 1000))) %>%
		+ coord_cartesian() %>%
		+ theme_bw(9) %>%
		+ theme(legend.position="bottom") %>%
		+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Save plot
ggsave(file.path(write_figures_path, "submodel_wyo_transect_elev.png"), p,  width = 6, height = 4.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_wyo_transect_elev.pdf"), p,  width = 6, height = 4.5)
ggsave(file.path(write_figures_path, "submodel_wyo_transect_elev.svg"), p,  width = 6, height = 4.5)



##########################################################################
###  Plot Wyoming model subset
###########################################################################
plot_wyo_profile <- plot_grid_sf %>%
filter(lat >= 43.4 & lat <= 43.5 & lon >= -113 & lon <= -108) %>%
as.data.frame() %>%
mutate(year = 1975)

### Add in ENSO terms
plot_wyo_profile <- expand_grid(plot_wyo_profile, data.frame(enso = seq(-2,2,0.5)))

plot_wyo_profile$pred <- predict(submodel_wyo, newdata = plot_wyo_profile, type = "response")

plot_wyo_profile <- plot_wyo_profile %>%
	mutate(pred_mm =  pred* 90)

p <- ggplot(plot_wyo_profile, aes(x=lon, y=pred, colour = enso, group = enso)) %>%
	+ geom_line(aes(y=pred_mm, colour = enso, group = enso)) %>%
	+ scale_colour_distiller(name = "ENSO", type = "div", palette = "PuOr", direction = -1, limits = c(-2,2)) %>%
	+ scale_x_continuous(name = "Longitude") %>%
	+ scale_y_continuous(name = "Predicted DJF \nPrecipitation (mm)") %>%
	+ theme_bw(9) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


### Save plot
ggsave(file.path(write_figures_path, "submodel_wyo_transect_subset.png"), p,  width = 5, height = 4.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_wyo_transect_subset.pdf"), p,  width = 5, height = 4.5)
ggsave(file.path(write_figures_path, "submodel_wyo_transect_subset.svg"), p,  width = 5, height = 4.5)


### I could put elevation on second axis
plot_wyo_profile_adj <- plot_wyo_profile %>%
	mutate(elev_adj = elevation/5 + 400/5)

p <- ggplot(plot_wyo_profile_adj, aes(x=lon)) %>%
	+ geom_line(aes(y=pred_mm, colour = enso, group = enso)) %>%
	+ geom_line(data = plot_wyo_profile_adj %>% filter(enso == 0), aes(y = elev_adj), colour = "grey50", size = 0.25) %>%
#	+ scale_colour_gradient2(name = "ENSO") %>%
	+ scale_colour_distiller(name = "ENSO", type = "div", palette = "PuOr", direction = -1, limits = c(-2,2)) %>%
	+ scale_x_continuous(name = "Longitude") %>%
	+ scale_y_continuous(name = "Predicted DJF \nPrecipitation (mm)" , sec.axis = sec_axis(~.*5-400, name="Elevation (m)", breaks = seq(0,4000, by = 1000))) %>%
	+ coord_cartesian() %>%
	+ theme_bw(9) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Save plot
ggsave(file.path(write_figures_path, "submodel_wyo_transect_elev_subset.png"), p,  width = 6, height = 4.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_wyo_transect_elev_subset.pdf"), p,  width = 6, height = 4.5)
ggsave(file.path(write_figures_path, "submodel_wyo_transect_elev_subset.svg"), p,  width = 6, height = 4.5)



### I could put elevation on second axis
plot_wyo_profile_adj <- plot_wyo_profile %>%
	mutate(elev_adj = elevation/5.5 + 600/5.5)

p <- ggplot(plot_wyo_profile_adj, aes(x=lon)) %>%
	+ geom_line(aes(y=pred_mm, colour = enso, group = enso)) %>%
	+ geom_line(data = plot_wyo_profile_adj %>% filter(enso == 0), aes(y = elev_adj), colour = "grey50", size = 0.25) %>%
#	+ scale_colour_gradient2(name = "ENSO") %>%
	+ scale_colour_distiller(name = "ENSO", type = "div", palette = "PuOr", direction = -1, limits = c(-2,2)) %>%
	+ scale_x_continuous(name = "Longitude") %>%
	+ scale_y_continuous(name = "Predicted DJF \nPrecipitation (mm)" , sec.axis = sec_axis(~.*5.5-600, name="Elevation (m)", breaks = seq(0,4000, by = 1000))) %>%
	+ coord_cartesian(ylim=c(0,650)) %>%
	+ theme_bw(11) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) %>%
	+  theme(legend.position="none")

	### Save plot
	ggsave(file.path(write_figures_path, "submodel_wyo_transect_elev_subset2.png"), p,  width = 6, height = 4.5, dpi = 600)
	ggsave(file.path(write_figures_path, "submodel_wyo_transect_elev_subset2.pdf"), p,  width = 6, height = 4.5)
	ggsave(file.path(write_figures_path, "submodel_wyo_transect_elev_subset2.svg"), p,  width = 6, height = 4.5)



##########################################################################
###  Plot spatial and elevation terms
###########################################################################

### Create object to hold results
#plot_df <- expand_grid(x_proj = x_seq, y_proj = y_seq, elevation = 1300, id = "USS0009E11S", enso = 0)
plot_df <- plot_grid_sf %>%
	mutate(elevation = 1300, enso = 0, year = 1975) %>%
	filter(region == "wyo")

### Extract the terms
partial_df<- predict(submodel_wyo, newdata = plot_df, type = "terms")

plot_df$spatial_partial <- partial_df[,1]
plot_df$spatial_partial_int <- (partial_df[,1]+coef(submodel_wyo)[1])

### Convert back to 3 month precip. Remember tweedie uses a log link
plot_df <- plot_df %>%
	mutate(spatial_partial_mm_djf = exp(spatial_partial) * 90) %>%
	mutate(spatial_partial_int_mm_djf = exp(spatial_partial_int) * 90)


### Create plot
p <- ggplot(data = plot_df) %>%
	#+ geom_stars(data = background, alpha = 0.1) %>%
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = spatial_partial)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = spatial_partial), colour = NA) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4", alpha = 0.2) %>%
	+ geom_sf(data = rivers, colour = "#4a80f5", alpha = 0.2) %>%
	+ geom_sf(data = wyo_pol, colour = "grey", fill = NA, alpha = 0.5) %>%
	#+ geom_sf(size = 0.1, alpha = 0.5) %>%
	#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
#	+ scale_fill_distiller(name = "Spatial Effect\n (Log space)", type = "div", palette = "RdBu", direction = 1, limits = c(-1.8,1.8)) %>%
	+ scale_fill_gradientn(name = "Spatial Effect\n (Log space)", colors = ipccPrec(11), limits = c(-1.8,1.8), oob=scales::squish) %>%
	+ coord_sf(xlim = c(-16e5, -7e5), ylim = c(19e5, 28e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#p

### Save plot
ggsave(file.path(write_figures_path, "submodel_wyo_spatial_term_log.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_wyo_spatial_term_log.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "submodel_wyo_spatial_term_log.svg"), p,  width = 4, height = 6.5)


###########################################################################
###  Plot the elevation term
###########################################################################
### Separate terms
djf_wyo <- djf_df %>%
	filter(region == "wyo")

plot_df <- data.frame(x_proj = median(djf_wyo$x_proj), y_proj = median(djf_wyo$y_proj),  enso = 0, year = 1975)
plot_df <- expand_grid(plot_df, data.frame(elevation = seq(min(djf_wyo$elevation), max(djf_wyo$elevation))))
plot_df



### Create an object to plot the rug plot for elevation
rug_df <- djf_df %>%
	filter(region == "wyo") %>%
	group_by(id) %>%
	summarize(elevation = median(elevation, na.rm=TRUE))

#partial_df<- predict(model_climatol, newdata = plot_df, type = "response", se.fit=TRUE)
partial_df<- predict(submodel_wyo, newdata = plot_df, type = "terms", se.fit=TRUE)

plot_df$elev_partial <- partial_df$fit[,4]
plot_df$elev_partial_int <- partial_df$fit[,4]+coef(submodel_wyo)[1]
plot_df$elev_partial_se <- partial_df$se.fit[,4]
plot_df <- plot_df %>%
	mutate(lower = elev_partial_int - 1.96*elev_partial_se) %>%
	mutate(upper = elev_partial_int + 1.96*elev_partial_se)

### Convert back to 3 month precip. Remember tweedie uses a log link
plot_df <- plot_df %>%
	mutate(elev_partial_mm_djf = exp(elev_partial) * 90) %>%
	mutate(elev_partial_int_mm_djf = exp(elev_partial_int) * 90) %>%
	mutate(lower_mm_djf = exp(lower) * 90) %>%
	mutate(upper_mm_djf = exp(upper) * 90)

### Create elevation plot
p <- ggplot(plot_df, aes(x=elevation)) %>%
	+ geom_hline(yintercept = 0, linetype = "dashed") %>%
	+ geom_ribbon(aes(y= elev_partial_int_mm_djf, ymin = lower_mm_djf, ymax = upper_mm_djf), alpha = 0.9, fill = "grey80")%>%
	+ geom_line(aes(y= elev_partial_int_mm_djf)) %>%
	+ geom_rug(data = rug_df, sides="b", alpha = 0.2)  %>%
	+ scale_x_continuous(name = "Elevation (m)") %>%
	+ scale_y_continuous(name = "Elevation Effect \nAdd. DJF Precipitation (mm)", breaks = seq(0,1000,100)) %>%
	+ coord_cartesian(xlim = c(0,3000)) %>%
	+ theme_bw(9)

### Save plot
ggsave(file.path(write_figures_path, "submodel_wyo_elevation_term.png"), p,  height = 3, width = 4, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_wyo_elevation_term.pdf"), p,  height = 3, width = 4)
ggsave(file.path(write_figures_path, "submodel_wyo_elevation_term.svg"), p,  height = 3, width = 4)


### Create elevation plot in log space
p <- ggplot(plot_df, aes(x=elevation)) %>%
	+ geom_hline(yintercept = 0, linetype = "dashed") %>%
	+ geom_ribbon(aes(y= elev_partial_int, ymin = lower, ymax = upper), alpha = 0.9, fill = "grey80")%>%
	+ geom_line(aes(y= elev_partial_int)) %>%
	+ geom_rug(data = rug_df, sides="b")  %>%
	+ scale_x_continuous(name = "Elevation (m)") %>%
	+ scale_y_continuous(name = "Elevation Effect \n(log space)") %>%
	+ coord_cartesian(xlim = c(0,3000)) %>%
	+ theme_bw(9)

### Save plot
ggsave(file.path(write_figures_path, "submodel_wyo_elevation_term_log.png"), p,  height = 3, width = 4, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_wyo_elevation_term_log.pdf"), p,  height = 3, width = 4)
ggsave(file.path(write_figures_path, "submodel_wyo_elevation_term_log.svg"), p,  height = 3, width = 4)


p <- draw(submodel_wyo, select = 4, residuals = TRUE)
#p
ggsave(file.path(write_figures_path, "submodel_wyo_elevation_term_gratia.png"), p,  height = 3, width = 4, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_wyo_elevation_term_gratia.pdf"), p,  height = 3, width = 4)
ggsave(file.path(write_figures_path, "submodel_wyo_elevation_term_gratia.svg"), p,  height = 3, width = 4 )



###########################################################################
###  Plot ENSO effect in precipitation
###########################################################################
### Instead, do each grid with its elevation and lat / lon for ENSO = 0
### Then subtract the other ENSOs from it

### Create estimate for enso = 0
### Create estimate for enso = 0
plot_base <- plot_grid_sf %>%
	filter(region == "wyo") %>%
	mutate(enso = 0, year = 1975)

plot_base <- plot_base %>%
	mutate(baseline = predict(submodel_wyo, newdata = plot_base, type = "response")) %>%
	as.data.frame()

#fpp <- predict(model_enso, newdata = plot_base, type = "terms")

### Create estimates for enso != 0
enso_list <- seq(-2,3,0.5)

for (j in seq(1, length(enso_list))){
	enso_j <- enso_list[[j]]

	plot_temp <- plot_grid_sf %>%
		filter(region == "wyo") %>%
		mutate(year = 1975) %>%
		mutate(enso = enso_j)

	if(j == 1){
		plot_enso <- plot_temp
	} else {
		plot_enso <- plot_enso %>%
			bind_rows(plot_temp)
	}
}

plot_enso <- plot_enso %>%
	filter(enso != 0) %>%
	left_join(plot_base %>% select(point, baseline), by = "point")

plot_enso <- plot_enso %>%
	mutate(pred = predict(submodel_wyo, newdata = plot_enso, type = "response")) %>%
	mutate(resid = 90*(pred - baseline)) %>%
	mutate(enso_label = paste0("MEI = ", enso))%>%
	mutate(enso_label = factor(enso_label, levels = paste0("MEI = ", enso_list)))

plot_lims <- max(abs(plot_enso$resid))
plot_lims <- quantile(abs(plot_enso$resid), 0.99)

### Create plot
p <- ggplot(data = plot_enso) %>%
	#+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	#+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	#+ geom_stars(data = background, alpha = 0.1) %>%
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = resid)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = resid), colour = NA) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5, size = 0.15) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5, size = 0.15) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5, size = 0.15) %>%
	#+ geom_sf(size = 0.1, alpha = 0.5) %>%
	#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
	+ facet_grid(.~enso_label) %>%
#+ scale_fill_distiller(name = "ENSO Effect\n (mm)", type = "div", palette = "RdBu", direction = 1, limits = c(-plot_lims,plot_lims), oob = squish) %>%
	+ scale_fill_gradientn(name = "ENSO Effect\n (mm)", colors = ipccPrec(11), limits = c(-plot_lims,plot_lims), oob=scales::squish) %>%
	+ coord_sf(xlim = c(-16e5, -7e5), ylim = c(19e5, 28e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p

### Save plot
ggsave(file.path(write_figures_path, "submodel_wyo_enso_effect_mm.png"), p,  width = 12, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_wyo_enso_effect_mm.pdf"), p,  width = 12, height = 6.5)
ggsave(file.path(write_figures_path, "submodel_wyo_enso_effect_mm.svg"), p,  width = 12, height = 6.5)



enso_list <- seq(-2,3,1)

### Fewer breaks
plot_sub <- plot_enso %>%
	filter(enso %in% enso_list)%>%
	mutate(enso_label = factor(enso_label, levels = paste0("MEI = ", enso_list)))

p <- ggplot(data = plot_sub) %>%
	#+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	#+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	#+ geom_stars(data = background, alpha = 0.1) %>%
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = resid)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = resid), colour = NA) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	#+ geom_sf(size = 0.1, alpha = 0.5) %>%
	#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
	+ facet_grid(.~enso_label) %>%
#	+ scale_fill_distiller(name = "ENSO Effect\n (mm)", type = "div", palette = "RdBu", direction = 1, limits = c(-plot_lims,plot_lims), oob = squish) %>%
	+ scale_fill_gradientn(name = "ENSO Effect\n (mm)", colors = ipccPrec(11), limits = c(-plot_lims,plot_lims), oob=scales::squish) %>%
	+ coord_sf(xlim = c(-16e5, -7e5), ylim = c(19e5, 28e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Save plot
ggsave(file.path(write_figures_path, "submodel_wyo_enso_effect_mm_subset.png"), p,  width = 12, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_wyo_enso_effect_mm_subset.pdf"), p,  width = 12, height = 6.5)
ggsave(file.path(write_figures_path, "submodel_wyo_enso_effect_mm_subset.svg"), p,  width = 12, height = 6.5)



enso_list <- seq(-2,2,1)

### Fewer breaks
plot_sub <- plot_enso %>%
	filter(enso %in% enso_list)%>%
	mutate(enso_label = factor(enso_label, levels = paste0("MEI = ", enso_list)))

p <- ggplot(data = plot_sub) %>%
	#+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	#+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	#+ geom_stars(data = background, alpha = 0.1) %>%
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = resid)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = resid), colour = NA) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	#+ geom_sf(size = 0.1, alpha = 0.5) %>%
	#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
	+ facet_grid(.~enso_label) %>%
#	+ scale_fill_distiller(name = "ENSO Effect\n (mm)", type = "div", palette = "RdBu", direction = 1, limits = c(-plot_lims,plot_lims), oob = squish) %>%
	+ scale_fill_gradientn(name = "ENSO Effect\n (mm)", colors = ipccPrec(11), limits = c(-plot_lims,plot_lims), oob=scales::squish) %>%
	+ coord_sf(xlim = c(-16e5, -7e5), ylim = c(19e5, 28e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Save plot
ggsave(file.path(write_figures_path, "submodel_wyo_enso_effect_mm_subset_noextreme.png"), p,  width = 12, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_wyo_enso_effect_mm_subset_noextreme.pdf"), p,  width = 12, height = 6.5)
ggsave(file.path(write_figures_path, "submodel_wyo_enso_effect_mm_subset_noextreme.svg"), p,  width = 12, height = 6.5)




##########################################################################
###  Fit Mexico model
###########################################################################
#### Do the same for Mexico
submodel_mx <-  bam(djf_mm_daily  ~  ti(x_proj, y_proj, d=c(2), bs = c("tp")) +
	 ti(x_proj, y_proj, enso, d=c(2,1), bs = c("tp", "tp")) +
	 ti(enso, bs = "tp") + s(elevation, bs = "tp")  + s(year, bs = "tp", k = 5),
						  data = djf_df %>% filter(region == "mx") ,
						 family =  tw())


BIC(submodel_mx)
AIC(submodel_mx)
summary(submodel_mx)



### Save the model objects
saveRDS(submodel_mx, file = file.path(model_folder, "submodel_mx.rds"))
# submodel_mx <- readRDS(file.path(model_folder, "submodel_mx.rds"))

### Save some residual information
resid_model <- resid(submodel_mx)

resid_stats <- data.frame(mean_resid = mean(resid_model, na.rm=TRUE),
	median_resid = median(resid_model, na.rm=TRUE),
	mae = mean(abs(resid_model), na.rm=TRUE),
	rmse = sqrt(mean((resid_model)^2)),
	aic = AIC(submodel_mx),
	bic = BIC(submodel_mx),
	r_square = summary(submodel_mx)$r.sq
	)

write_csv(resid_stats, file.path(model_folder, "submodel_mx_resid_stats.csv"))


p <- appraise(submodel_mx)
p

### Save plot
ggsave(file.path(write_figures_path, "submodel_mx_appraise.png"), p,  width = 7, height = 7, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_mx_appraise.pdf"), p,  width = 7, height = 7)
ggsave(file.path(write_figures_path, "submodel_mx_appraise.svg"), p,  width = 7, height = 7)

p <- draw(submodel_mx)
p

### Save plot
ggsave(file.path(write_figures_path, "submodel_mx_draw.png"), p,  width = 7, height = 7, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_mx_draw.pdf"), p,  width = 7, height = 7)
ggsave(file.path(write_figures_path, "submodel_mx_draw.svg"), p,  width = 7, height = 7)


##########################################################################
###  Plot Mexico model
###########################################################################
plot_mx_profile <- plot_grid_sf %>%
	filter(lat >= 27.9500 & lat <= 28.05001 & lon >= -110 & lon <= -105) %>%
	as.data.frame()


### Add in ENSO terms
plot_mx_profile <- expand_grid(plot_mx_profile, data.frame(enso = seq(-2,3,0.5)), year = 1975)

plot_mx_profile$pred <- predict(submodel_mx, newdata = plot_mx_profile, type = "response")

plot_mx_profile <- plot_mx_profile %>%
	mutate(pred_mm =  pred* 90)

p <- ggplot(plot_mx_profile, aes(x=lon, y=pred, colour = enso, group = enso)) %>%
	+ geom_line(aes(y=pred_mm, colour = enso, group = enso)) %>%
	+ scale_colour_distiller(name = "ENSO", type = "div", palette = "PuOr", direction = -1, limits = c(-3,3)) %>%
	+ scale_x_continuous(name = "Longitude") %>%
	+ scale_y_continuous(name = "Predicted DJF \nPrecipitation (mm)") %>%
	+ theme_bw(9) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


### Save plot
ggsave(file.path(write_figures_path, "submodel_mx_transect.png"), p,  width = 5, height = 4.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_mx_transect.pdf"), p,  width = 5, height = 4.5)
ggsave(file.path(write_figures_path, "submodel_mx_transect.svg"), p,  width = 5, height = 4.5)



	plot_mx_profile <- plot_mx_profile %>%
		mutate(elev_adj = elevation/7 + 1400/7)


	p <- ggplot(plot_mx_profile, aes(x=lon)) %>%
		+ geom_line(aes(y=pred_mm, colour = enso, group = enso)) %>%
		+ geom_line(data = plot_mx_profile %>% filter(enso == 0), aes(y = elev_adj), colour = "grey50", size = 0.25) %>%
	#	+ scale_colour_gradient2(name = "ENSO") %>%
		+ scale_colour_distiller(name = "ENSO", type = "div", palette = "PuOr", direction = -1, limits = c(-3,3)) %>%
		+ scale_x_continuous(name = "Longitude") %>%
		+ scale_y_continuous(name = "Predicted DJF \nPrecipitation (mm)" , sec.axis = sec_axis(~.*7-1400, name="Elevation (m)", breaks = seq(0,4000, by = 1000))) %>%
		+ coord_cartesian() %>%
		+ theme_bw(9) %>%
		+ theme(legend.position="bottom") %>%
		+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


		### Save plot
		ggsave(file.path(write_figures_path, "submodel_mx_transect_elev.png"), p,  width = 6, height = 4.5, dpi = 600)
		ggsave(file.path(write_figures_path, "submodel_mx_transect_elev.pdf"), p,  width = 6, height = 4.5)
		ggsave(file.path(write_figures_path, "submodel_mx_transect_elev.svg"), p,  width = 6, height = 4.5)




##########################################################################
###  Plot Mexico model subset
###########################################################################
plot_mx_profile <- plot_grid_sf %>%
	filter(lat >= 27.9500 & lat <= 28.05001 & lon >= -110 & lon <= -105) %>%
	as.data.frame()

### Add in ENSO terms
plot_mx_profile <- expand_grid(plot_mx_profile, data.frame(enso = seq(-2,2,0.5), year = 1975))

plot_mx_profile$pred <- predict(submodel_mx, newdata = plot_mx_profile, type = "response")

plot_mx_profile <- plot_mx_profile %>%
	mutate(pred_mm =  pred* 90)

p <- ggplot(plot_mx_profile, aes(x=lon, y=pred, colour = enso, group = enso)) %>%
	+ geom_line(aes(y=pred_mm, colour = enso, group = enso)) %>%
	+ scale_colour_distiller(name = "ENSO", type = "div", palette = "PuOr", direction = -1, limits = c(-2,2)) %>%
	+ scale_x_continuous(name = "Longitude") %>%
	+ scale_y_continuous(name = "Predicted DJF \nPrecipitation (mm)") %>%
	+ theme_bw(9) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


### Save plot
ggsave(file.path(write_figures_path, "submodel_mx_transect_subset.png"), p,  width = 5, height = 4.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_mx_transect_subset.pdf"), p,  width = 5, height = 4.5)
ggsave(file.path(write_figures_path, "submodel_mx_transect_subset.svg"), p,  width = 5, height = 4.5)


### I could put elevation on second axis

	plot_mx_profile_adj <- plot_mx_profile %>%
		mutate(elev_adj = elevation/7 + 1400/7)


	p <- ggplot(plot_mx_profile_adj, aes(x=lon)) %>%
		+ geom_line(aes(y=pred_mm, colour = enso, group = enso)) %>%
		+ geom_line(data = plot_mx_profile_adj %>% filter(enso == 0), aes(y = elev_adj), colour = "grey50", size = 0.25) %>%
	#	+ scale_colour_gradient2(name = "ENSO") %>%
		+ scale_colour_distiller(name = "ENSO", type = "div", palette = "PuOr", direction = -1, limits = c(-2,2)) %>%
		+ scale_x_continuous(name = "Longitude") %>%
		+ scale_y_continuous(name = "Predicted DJF \nPrecipitation (mm)" , sec.axis = sec_axis(~.*7-1400, name="Elevation (m)", breaks = seq(0,4000, by = 1000))) %>%
		+ coord_cartesian() %>%
		+ theme_bw(9) %>%
		+ theme(legend.position="bottom") %>%
		+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


### Save plot
ggsave(file.path(write_figures_path, "submodel_mx_transect_elev_subset.png"), p,  width = 6, height = 4.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_mx_transect_elev_subset.pdf"), p,  width = 6, height = 4.5)
ggsave(file.path(write_figures_path, "submodel_mx_transect_elev_subset.svg"), p,  width = 6, height = 4.5)




### I could put elevation on second axis

	plot_mx_profile_adj <- plot_mx_profile %>%
		mutate(elev_adj = elevation/5.5 + 600/5.5)

	p <- ggplot(plot_mx_profile_adj, aes(x=lon)) %>%
		+ geom_line(aes(y=pred_mm, colour = enso, group = enso)) %>%
		+ geom_line(data = plot_mx_profile_adj %>% filter(enso == 0), aes(y = elev_adj), colour = "grey50", size = 0.25) %>%
	#	+ scale_colour_gradient2(name = "ENSO") %>%
		+ scale_colour_distiller(name = "ENSO", type = "div", palette = "PuOr", direction = -1, limits = c(-2,2)) %>%
		+ scale_x_continuous(name = "Longitude") %>%
		+ scale_y_continuous(name = "Predicted DJF \nPrecipitation (mm)" , sec.axis = sec_axis(~.*5.5-600, name="Elevation (m)", breaks = seq(0,4000, by = 1000))) %>%
		+ coord_cartesian(ylim = c(0,650)) %>%
		+ theme_bw(11) %>%
		+ theme(legend.position="bottom") %>%
		+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


### Save plot
ggsave(file.path(write_figures_path, "submodel_mx_transect_elev_subset_legend.png"), p,  width = 6, height = 4.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_mx_transect_elev_subset_legend.pdf"), p,  width = 6, height = 4.5)
ggsave(file.path(write_figures_path, "submodel_mx_transect_elev_subset_legend.svg"), p,  width = 6, height = 4.5)


p <- p +  theme(legend.position="none")

### Save plot
ggsave(file.path(write_figures_path, "submodel_mx_transect_elev_subset2.png"), p,  width = 6, height = 4.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_mx_transect_elev_subset2.pdf"), p,  width = 6, height = 4.5)
ggsave(file.path(write_figures_path, "submodel_mx_transect_elev_subset2.svg"), p,  width = 6, height = 4.5)



##########################################################################
###  Plot spatial and elevation terms
###########################################################################
### Create object to hold results
#plot_df <- expand_grid(x_proj = x_seq, y_proj = y_seq, elevation = 1300, id = "USS0009E11S", enso = 0)
plot_df <- plot_grid_sf %>%
	mutate(elevation = 1300, enso = 0, year = 1975) %>%
	filter(region == "mx")

### Extract the terms
partial_df<- predict(submodel_mx, newdata = plot_df, type = "terms")

plot_df$spatial_partial <- partial_df[,1]
plot_df$spatial_partial_int <- (partial_df[,1]+coef(submodel_mx)[1])

### Convert back to 3 month precip. Remember tweedie uses a log link
plot_df <- plot_df %>%
	mutate(spatial_partial_mm_djf = exp(spatial_partial) * 90) %>%
	mutate(spatial_partial_int_mm_djf = exp(spatial_partial_int) * 90)


### Create plot
p <- ggplot(data = plot_df) %>%
	#+ geom_stars(data = background, alpha = 0.1) %>%
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = spatial_partial)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = spatial_partial), colour = NA) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4", alpha = 0.2) %>%
	+ geom_sf(data = rivers, colour = "#4a80f5", alpha = 0.2) %>%
	+ geom_sf(data = mx_pol, colour = "orange", fill = NA, alpha = 0.5) %>%
	#+ geom_sf(size = 0.1, alpha = 0.5) %>%
	#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
	#+ scale_fill_distiller(name = "Spatial Effect\n (Log space)", type = "div", palette = "RdBu", direction = 1, limits = c(-1.8,1.8)) %>%
	+ scale_fill_gradientn(name = "Spatial Effect\n (Log space)", colors = ipccPrec(11), limits = c(-1.8,1.8), oob=scales::squish) %>%

	+ coord_sf(xlim =  c(-17.5e5, -7e5), ylim = c(2e5, 10.5e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#p

### Save plot
ggsave(file.path(write_figures_path, "submodel_mx_spatial_term_log.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_mx_spatial_term_log.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "submodel_mx_spatial_term_log.svg"), p,  width = 4, height = 6.5)



###########################################################################
###  Plot the elevation term
###########################################################################
### Separate terms
djf_mx <- djf_df %>%
	filter(region == "mx")

plot_df <- data.frame(x_proj = median(djf_mx$x_proj), y_proj = median(djf_mx$y_proj),  enso = 0, year = 1975)
plot_df <- expand_grid(plot_df, data.frame(elevation = seq(min(djf_mx$elevation), max(djf_mx$elevation))))
plot_df

### Create an object to plot the rug plot for elevation
rug_df <- djf_df %>%
	filter(region == "mx") %>%
	group_by(id) %>%
	summarize(elevation = median(elevation, na.rm=TRUE))

#partial_df<- predict(model_climatol, newdata = plot_df, type = "response", se.fit=TRUE)
partial_df<- predict(submodel_mx, newdata = plot_df, type = "terms", se.fit=TRUE)

plot_df$elev_partial <- partial_df$fit[,4]
plot_df$elev_partial_int <- partial_df$fit[,4]+coef(submodel_mx)[1]
plot_df$elev_partial_se <- partial_df$se.fit[,4]
plot_df <- plot_df %>%
	mutate(lower = elev_partial_int - 1.96*elev_partial_se) %>%
	mutate(upper = elev_partial_int + 1.96*elev_partial_se)

### Convert back to 3 month precip. Remember tweedie uses a log link
plot_df <- plot_df %>%
	mutate(elev_partial_mm_djf = exp(elev_partial) * 90) %>%
	mutate(elev_partial_int_mm_djf = exp(elev_partial_int) * 90) %>%
	mutate(lower_mm_djf = exp(lower) * 90) %>%
	mutate(upper_mm_djf = exp(upper) * 90)

### Create elevation plot
p <- ggplot(plot_df, aes(x=elevation)) %>%
	+ geom_hline(yintercept = 0, linetype = "dashed") %>%
	+ geom_ribbon(aes(y= elev_partial_int_mm_djf, ymin = lower_mm_djf, ymax = upper_mm_djf), alpha = 0.9, fill = "grey80")%>%
	+ geom_line(aes(y= elev_partial_int_mm_djf)) %>%
	+ geom_rug(data = rug_df, sides="b", alpha = 0.2)  %>%
	+ scale_x_continuous(name = "Elevation (m)") %>%
	+ scale_y_continuous(name = "Elevation Effect \nAdd. DJF Precipitation (mm)", breaks = seq(0,1000,100)) %>%
	+ coord_cartesian(xlim = c(0,3000)) %>%
	+ theme_bw(9)

### Save plot
ggsave(file.path(write_figures_path, "submodel_mx_elevation_term.png"), p,  height = 3, width = 4, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_mx_elevation_term.pdf"), p,  height = 3, width = 4)
ggsave(file.path(write_figures_path, "submodel_mx_elevation_term.svg"), p,  height = 3, width = 4)


### Create elevation plot in log space
p <- ggplot(plot_df, aes(x=elevation)) %>%
	+ geom_hline(yintercept = 0, linetype = "dashed") %>%
	+ geom_ribbon(aes(y= elev_partial_int, ymin = lower, ymax = upper), alpha = 0.9, fill = "grey80")%>%
	+ geom_line(aes(y= elev_partial_int)) %>%
	+ geom_rug(data = rug_df, sides="b")  %>%
	+ scale_x_continuous(name = "Elevation (m)") %>%
	+ scale_y_continuous(name = "Elevation Effect \n(log space)") %>%
	+ coord_cartesian(xlim = c(0,3000)) %>%
	+ theme_bw(9)

### Save plot
ggsave(file.path(write_figures_path, "submodel_mx_elevation_term_log.png"), p,  height = 3, width = 4, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_mx_elevation_term_log.pdf"), p,  height = 3, width = 4)
ggsave(file.path(write_figures_path, "submodel_mx_elevation_term_log.svg"), p,  height = 3, width = 4)


p <- draw(submodel_mx, select = 4 , residuals = TRUE)
#p
ggsave(file.path(write_figures_path, "submodel_mx_elevation_term_gratia.png"), p,  height = 3, width = 4, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_mx_elevation_term_gratia.pdf"), p,  height = 3, width = 4)
ggsave(file.path(write_figures_path, "submodel_mx_elevation_term_gratia.svg"), p,  height = 3, width = 4 )



###########################################################################
###  Plot ENSO effect in precipitation
###########################################################################
### Instead, do each grid with its elevation and lat / lon for ENSO = 0
### Then subtract the other ENSOs from it

### Create estimate for enso = 0
	plot_base <- plot_grid_sf %>%
		filter(region == "mx") %>%
		mutate(enso = 0, year = 1975)

	plot_base <- plot_base %>%
		mutate(baseline = predict(submodel_mx, newdata = plot_base, type = "response")) %>%
		as.data.frame()

### Create estimates for enso != 0
enso_list <- seq(-2,3,0.5)

for (j in seq(1, length(enso_list))){
	plot_temp <- plot_grid_sf %>%
		filter(region == "mx") %>%
		#mutate(elevation = 1300, id = "USS0009E11S") %>%
		mutate(enso = enso_list[[j]], year = 1975)

	if(j == 1){
		plot_enso <- plot_temp
	} else {
		plot_enso <- plot_enso %>%
			bind_rows(plot_temp)
	}
}

plot_enso <- plot_enso %>%
	filter(enso != 0) %>%
	left_join(plot_base %>% select(point, baseline), by = "point")

plot_enso <- plot_enso %>%
	mutate(pred = predict(submodel_mx, newdata = plot_enso, type = "response")) %>%
	mutate(resid = 90*(pred - baseline)) %>%
	mutate(enso_label = paste0("MEI = ", enso))%>%
	mutate(enso_label = factor(enso_label, levels = paste0("MEI = ", enso_list)))

plot_lims <- max(abs(plot_enso$resid))
plot_lims <- quantile(abs(plot_enso$resid), 0.99)




### Create plot
p <- ggplot(data = plot_enso) %>%
	#+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	#+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	#+ geom_stars(data = background, alpha = 0.1) %>%
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = resid)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = resid), colour = NA) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5, size = 0.15) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5, size = 0.15) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5, size = 0.15) %>%
	#+ geom_sf(size = 0.1, alpha = 0.5) %>%
	#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
	+ facet_grid(.~enso_label) %>%
#	+ scale_fill_distiller(name = "ENSO Effect\n (mm)", type = "div", palette = "RdBu", direction = 1, limits = c(-plot_lims,plot_lims), oob = squish) %>%
	+ scale_fill_gradientn(name = "ENSO Effect\n (mm)", colors = ipccPrec(11), limits = c(-plot_lims,plot_lims), oob=scales::squish) %>%

	+ coord_sf(xlim =  c(-17.5e5, -7e5), ylim = c(2e5, 10.5e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#p


### Save plot
ggsave(file.path(write_figures_path, "submodel_mx_enso_effect_mm.png"), p,  width = 12, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_mx_enso_effect_mm.pdf"), p,  width = 12, height = 6.5)
ggsave(file.path(write_figures_path, "submodel_mx_enso_effect_mm.svg"), p,  width = 12, height = 6.5)



enso_list <- seq(-2,3,1)

### Fewer breaks
plot_sub <- plot_enso %>%
	filter(enso %in% enso_list)%>%
	mutate(enso_label = factor(enso_label, levels = paste0("MEI = ", enso_list)))

p <- ggplot(data = plot_sub) %>%
	#+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	#+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	#+ geom_stars(data = background, alpha = 0.1) %>%
#	+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = resid)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = resid), colour = NA) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	#+ geom_sf(size = 0.1, alpha = 0.5) %>%
	#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
	+ facet_grid(.~enso_label) %>%
#	+ scale_fill_distiller(name = "ENSO Effect\n (mm)", type = "div", palette = "RdBu", direction = 1, limits = c(-plot_lims,plot_lims), oob = squish) %>%
	+ scale_fill_gradientn(name = "ENSO Effect\n (mm)", colors = ipccPrec(11), limits = c(-plot_lims,plot_lims), oob=scales::squish) %>%

	+ coord_sf(xlim =  c(-17.5e5, -7e5), ylim = c(2e5, 10.5e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Save plot
ggsave(file.path(write_figures_path, "submodel_mx_enso_effect_mm_subset.png"), p,  width = 12, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "submodel_mx_enso_effect_mm_subset.pdf"), p,  width = 12, height = 6.5)
ggsave(file.path(write_figures_path, "submodel_mx_enso_effect_mm_subset.svg"), p,  width = 12, height = 6.5)






###########################################################################
###  Plot Annual effect in precipitation
###########################################################################
p <- draw(submodel_mx, select = 5) + theme_bw(9)

				ggsave(file.path(write_figures_path, "submodel_mx_annual.png"), p,  height = 4, width = 6.5, dpi = 600)
				ggsave(file.path(write_figures_path, "submodel_mx_annual.svg"), p,  height = 4, width = 6.5)

p <- draw(submodel_wyo, select = 5) + theme_bw(9)

								ggsave(file.path(write_figures_path, "submodel_wyo_annual.png"), p,  height = 4, width = 6.5, dpi = 600)
								ggsave(file.path(write_figures_path, "submodel_wyo_annual.svg"), p,  height = 4, width = 6.5)
