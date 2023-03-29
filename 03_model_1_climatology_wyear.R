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

cat("Starting 03_model_1_climatology")

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

require(scales)

require(rnoaa)
require(tidyverse)
library(sf)
library(stars)
library(rnaturalearth)
library(rnaturalearthhires)

#require(lfstat)

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
###  Create simple climatology
###########################################################################
model_folder <- file.path(output_path, "model_1_climatol_trend")
dir.create(model_folder, recursive=TRUE, showWarnings = FALSE)

### Calculate median DJF precip
djf_summary <- djf_df %>%
	group_by(id) %>%
	summarize(djf_mm_daily_median = median(djf_mm_daily, na.rm=TRUE), djf_mm_daily_mean = mean(djf_mm_daily, na.rm=TRUE), n = n(), elevation = median(elevation), latitude = median(latitude), longitude = median(longitude), x_proj = median(x_proj), y_proj = median(y_proj)) #%>%
	#filter(djf_mm_daily_median > 0)

### Row 126 seems very strange outlier
djf_summary <- djf_summary %>%
	filter(id != "MXN00008118") %>%
	mutate(id = as.factor(id))


###########################################################################
###  Create model 1 - Climatology
###########################################################################, bs = "ds" , k = c(15,15)

### Create the model
model_climatol  <-  bam(djf_mm_daily  ~ te(x_proj,y_proj, bs = "tp", k = c(15,15))+
	 					 s(elevation, bs = "ad", k = 15) +
						  te(x_proj, y_proj, year, d=c(2,1), k = c(8,3), bs = c("tp", "tp")),
						  data = djf_df,
						 family =  tw(),
						 samfrac = 0.1,
						 # discrete=TRUE,
					  nthreads=8)


cat("Finished Climatology Model Fitting")
#model_climatol <- readRDS(file.path(model_folder, "model_climatol_trend.rds"))
#model_climatol_df <- readRDS(file.path(write_output_path, "model_1_climatol_trend/model_climatol_trend_df.rds"))

### Check the AIC / BIC
AIC(model_climatol)
BIC(model_climatol)

### Save the fit
model_climatol_df <- djf_df
model_climatol_df$pred <- fitted(model_climatol)
model_climatol_df$resid <- resid(model_climatol)

### Save the model objects
saveRDS(model_climatol, file = file.path(model_folder, "model_climatol_trend.rds"))
saveRDS(model_climatol_df, file = file.path(model_folder, "model_climatol_trend_df.rds"))


### Save some residual information
resid_model <- resid(model_climatol)

resid_stats <- data.frame(mean_resid = mean(resid_model, na.rm=TRUE),
	median_resid = median(resid_model, na.rm=TRUE),
	mae = mean(abs(resid_model), na.rm=TRUE),
	rmse = sqrt(mean((resid_model)^2)),
	aic = AIC(model_climatol),
	bic = BIC(model_climatol),
	r_square = summary(model_climatol)$r.sq
	)

write_csv(resid_stats, file.path(model_folder, "model_climatol_resid_stats.csv"))


###########################################################################
###  Plot the results
###########################################################################
write_figures_path <- file.path(write_figures_path, "model_1_climatol_trend")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)

p <- draw(model_climatol) + theme_bw(9)
ggsave(file.path(write_figures_path, "draw_model.png"), p,  width = 10, height = 6.5, dpi = 600)

p <- appraise(model_climatol) + theme_bw(9)
ggsave(file.path(write_figures_path, "appraise_model.png"), p,  width = 10, height = 6.5, dpi = 600)

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


###########################################################################
###  Plot Model Comparison
###########################################################################
#plot_df <- gauges_sf %>%
#	inner_join(model_climatol_df %>% select(id, djf_mm_daily_median, pred, resid), by = "id")

plot_df <- djf_summary

my_breaks <- c(2, 5, 10, 20, 50, 100, 200, 500)

p <- ggplot(data = plot_df) %>%
	#+ geom_stars(data = background, alpha = 0.2) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
#	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_point(alpha = 0.95, aes(x = x_proj, y=y_proj, colour = djf_mm_daily_median * 90), size = 0.7) %>%
	#+ geom_sf(aes(colour = djf_mm_daily_median * 90), size = 0.7) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	#+ scale_fill_identity() %>%
	+ scale_colour_viridis(name = "Median\nWinter (DJF)\nPrecip (mm)", trans = "log", option =  "inferno" , direction = 1, limits = c(4, 750), breaks = my_breaks) %>%
	#+ scale_colour_distiller(name = "Median\nWinter (DJF)\nPrecip (m)", trans = "log", palette = "RdYlBu", direction = -1,  limits = c(4, 750), breaks = my_breaks) %>%
	#+ scale_colour_distiller(name = "Median\nWinter (DJF)\nPrecip (m)", trans = "log", palette = "RdPu", direction = 1, limits = c(4, 750), breaks = my_breaks) %>%
#	+ coord_sf(xlim = c(-115, -100), ylim = c(29, 50), expand = FALSE) %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p


### Save plot
ggsave(file.path(write_figures_path, "true_precip_clim.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "true_precip_clim.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "true_precip_clim.svg"), p,  width = 4, height = 6.5)

p <- p +  scale_colour_distiller(name = "Median\nWinter (DJF)\nPrecip (mm)", trans = "log", palette = "RdYlBu", direction = -1,  limits = c(4, 750), breaks = my_breaks)


### Save plot
ggsave(file.path(write_figures_path, "true_precip_clim_alt.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "true_precip_clim_alt.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "true_precip_clim_alt.svg"), p,  width = 4, height = 6.5)


plot_df <- model_climatol_df %>%
	filter(year == 1975)

p <- ggplot(data = plot_df) %>%
	#+ geom_stars(data = background, alpha = 0.2) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
#	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	#+ geom_sf(aes(colour = pred * 90), size = 0.7) %>%
	+ geom_point(alpha = 0.95, aes(x = x_proj, y=y_proj, colour = pred * 90), size = 0.7) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ scale_fill_identity() %>%
	+ scale_colour_viridis(name = "Predicted\nWinter (DJF)\nPrecip (mm)", trans = "log", option =  "inferno" , direction = 1, limits = c(4, 750), breaks = my_breaks) %>%
	#+ scale_colour_distiller(name = "Median\nWinter (DJF)\nPrecip (m)", trans = "log", palette = "RdYlBu", direction = -1,  limits = c(4, 750), breaks = my_breaks) %>%
	#+ scale_colour_distiller(name = "Median\nWinter (DJF)\nPrecip (m)", trans = "log", palette = "RdPu", direction = 1, limits = c(4, 750), breaks = my_breaks) %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#p

### Save plot
ggsave(file.path(write_figures_path, "pred_precip_clim_1975.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "pred_precip_clim_1975.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "pred_precip_clim_1975.svg"), p,  width = 4, height = 6.5)

p <- p +  scale_colour_distiller(name = "Predicted\nWinter (DJF)\nPrecip (mm)", trans = "log", palette = "RdYlBu", direction = -1,  limits = c(4, 750), breaks = my_breaks)


### Save plot
ggsave(file.path(write_figures_path, "pred_precip_clim_alt.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "pred_precip_clim_alt.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "pred_precip_clim_alt.svg"), p,  width = 4, height = 6.5)




##### Started here

### Plot predicted versus observed
plot_df <- model_climatol_df

p <- ggplot(plot_df, aes(x= djf_mm_daily * 90, y = pred * 90)) %>%
	+ geom_abline(intercept = 0, slope = 1, colour = "blue") %>%
	+ geom_point(alpha = 0.7, size = 0.1) %>%
	+ scale_y_continuous(name = "Predicted Winter (DJF) Precip (mm)", breaks = seq(0,2000,200)) %>%
	+ scale_x_continuous(name = "Observed Winter (DJF) Precip (mm)", breaks = seq(0,2000,200)) %>%
	+ coord_fixed(ratio = 1) %>%
	+ theme_bw(9) %>%
	  + theme(panel.grid.minor = element_blank())
#	p

### Save plot
ggsave(file.path(write_figures_path, "clim_pred_vs_obs.png"), p,  width = 4, height = 4, dpi = 600)
ggsave(file.path(write_figures_path, "clim_pred_vs_obs.pdf"), p,  width = 4, height = 4)
ggsave(file.path(write_figures_path, "clim_pred_vs_obs.svg"), p,  width = 4, height = 4)


log_breaks <- c(1, 10, 100,  1000, 10000)

p <- ggplot(plot_df, aes(x= djf_mm_daily * 90, y = pred * 90)) %>%
	+ geom_hex() %>%
	+ geom_abline(intercept = 0, slope = 1, colour = "#d95f02") %>%
	+ scale_y_continuous(name = "Predicted Winter (DJF) Precip (mm)", breaks = seq(0,2000,200)) %>%
	+ scale_x_continuous(name = "Observed Winter (DJF) Precip (mm)", breaks = seq(0,2000,200)) %>%
	+ scale_fill_viridis(name = "Count", trans = "log", breaks = log_breaks, option = "cividis") %>%
	+ coord_fixed(ratio = 1) %>%
	+ theme_bw(9) %>%
	  + theme(panel.grid.minor = element_blank())
	p


	### Save plot
	ggsave(file.path(write_figures_path, "clim_pred_vs_obs_hex.png"), p,  width = 5, height = 4.5, dpi = 600)
	ggsave(file.path(write_figures_path, "clim_pred_vs_obs_hex.pdf"), p,  width = 5, height = 4.5)
	ggsave(file.path(write_figures_path, "clim_pred_vs_obs_hex.svg"), p,  width = 5, height = 4.5)




log_breaks <- c(0.1,0.5, 1, 5, 10, 50, 100, 500, 1000, 5000)

p <- ggplot(plot_df, aes(x= djf_mm_daily * 90, y = pred * 90)) %>%
	+ geom_abline(intercept = 0, slope = 1, colour = "blue") %>%
	+ geom_point() %>%
	+  scale_y_log10(name = "Predicted Winter (DJF) Precip (mm)", breaks = log_breaks ) %>%
	+ scale_x_log10(name = "Observed Winter (DJF) Precip (mm)", breaks = log_breaks) %>%
  + theme_bw(9) %>%
 	+ coord_fixed(ratio = 1, ylim = c(8, 2000), xlim = c(2, 2000)) %>%
  + annotation_logticks(sides = 'lrbt') %>%
  + theme(panel.grid.minor = element_blank())


### Save plot
ggsave(file.path(write_figures_path, "clim_pred_vs_obs_log.png"), p,  width = 4, height = 4, dpi = 600)
ggsave(file.path(write_figures_path, "clim_pred_vs_obs_log.pdf"), p,  width = 4, height = 4)
ggsave(file.path(write_figures_path, "clim_pred_vs_obs_log.svg"), p,  width = 4, height = 4)


###########################################################################
###  Plot prediction gridded
###########################################################################
### Create
plot_df <- plot_grid_sf %>%
	mutate(id = "USS0009E11S") %>%
	mutate(year = 1975)

est_response <- predict(model_climatol, newdata = plot_df, type = "response")

### Extract the terms
plot_df <- plot_df %>%
	mutate(pred = est_response * 90)

### Create plot
p <- ggplot(data = plot_df ) %>%
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = pred)) %>%
	#+ geom_stars(data = background, alpha = 0.1) %>%
	+ geom_sf(alpha = 0.95, aes(fill = pred), colour = NA) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4", alpha = 0.5) %>%
	+ geom_sf(data = rivers, colour = "#4a80f5", alpha = 0.5) %>%
	#+ geom_stars(data = background, alpha = 0.1) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ scale_fill_viridis(name = "Predicted\nWinter (DJF)\nPrecip (mm)", trans = "log", option =  "inferno" , direction = 1, limits = c(4, 750), breaks = my_breaks) %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	#+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


### Save plot
ggsave(file.path(write_figures_path, "pred_precip_clim_grid.png"), p,  width = 7, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "pred_precip_clim_grid.pdf"), p,  width = 7, height = 6.5)
ggsave(file.path(write_figures_path, "pred_precip_clim_grid.svg"), p,  width = 7, height = 6.5)



###########################################################################
###  Plot Spatial terms
###########################################################################
require(raster)

### Create object to hold results
plot_df <- plot_grid_sf %>%
	mutate(elevation = 1300, id = "USS0009E11S", year = 1975)

### Extract the terms
partial_df<- predict(model_climatol, newdata = plot_df, type = "terms")

plot_df$spatial_partial <- partial_df[,1]
plot_df$spatial_partial_int <- (partial_df[,1]+coef(model_climatol)[1])

### Convert back to 3 month precip. Remember tweedie uses a log link
plot_df <- plot_df %>%
	mutate(spatial_partial_mm_djf = exp(spatial_partial) * 90) %>%
	mutate(spatial_partial_int_mm_djf = exp(spatial_partial_int) * 90)


### Create plot
p <- ggplot(data = plot_df) %>%
	#+ geom_stars(data = background, alpha = 0.1) %>%
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = spatial_partial)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = spatial_partial), colour = NA) %>%
	+ geom_sf(data = rivers, colour = "#4a80f5", alpha = 0.15) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4", alpha = 0.15) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	#+ geom_sf(size = 0.1, alpha = 0.5) %>%
	#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
	+ scale_fill_distiller(name = "Spatial Effect\n (Log space)", type = "div", palette = "RdBu", direction = 1, limits = c(-3.25, 3.25), oob = squish) %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#p

### Save plot
ggsave(file.path(write_figures_path, "map_spatial_term_log.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "map_spatial_term_log.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "map_spatial_term_log.svg"), p,  width = 4, height = 6.5)

### Create plot
p <- ggplot(data = plot_df) %>%
	#+ geom_stars(data = background, alpha = 0.1) %>%
#	+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = spatial_partial_int_mm_djf)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = spatial_partial_int_mm_djf), colour = NA) %>%
	+ geom_sf(data = rivers, colour = "#4a80f5", alpha = 0.15) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4", alpha = 0.15) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	#+ geom_sf(size = 0.1, alpha = 0.5) %>%
	#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
	#+ scale_fill_distiller(name = "Spatial Effect\n Add. DJF \nPrecipitation (mm)", type = "div", palette = "RdBu", direction = 1, limits = c(-160,160)) %>%
	+ scale_fill_viridis(name = "Spatial Effect\n Add. DJF \nPrecipitation (mm)", trans = "log", option =  "inferno" , direction = -1, limits = c(4, 750), breaks = my_breaks, oob = squish) %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#p

### Save plot
ggsave(file.path(write_figures_path, "map_spatial_term.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "map_spatial_term.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "map_spatial_term.svg"), p,  width = 4, height = 6.5)




###########################################################################
###  Plot Elevation term
###########################################################################

### Separate terms
plot_df <- expand_grid(x_proj = median(djf_summary$x_proj), y_proj = median(djf_summary$y_proj), longitude = -108, latitude = 40, elevation = seq(0,3500), id = "USS0009E11S", year = 1975)

#partial_df<- predict(model_climatol, newdata = plot_df, type = "response", se.fit=TRUE)
partial_df<- predict(model_climatol, newdata = plot_df, type = "terms", se.fit=TRUE)

plot_df$elev_partial <- partial_df$fit[,2]
plot_df$elev_partial_int <- partial_df$fit[,2]+coef(model_climatol)[1]
plot_df$elev_partial_se <- partial_df$se.fit[,2]
plot_df <- plot_df %>%
	mutate(lower = elev_partial_int - 1.96*elev_partial_se) %>%
	mutate(upper = elev_partial_int + 1.96*elev_partial_se)

### Convert back to 3 month precip. Remember tweedie uses a log link
plot_df <- plot_df %>%
	mutate(elev_partial_mm_djf = exp(elev_partial) * 90) %>%
	mutate(elev_partial_int_mm_djf = exp(elev_partial_int) * 90) %>%
	mutate(lower_mm_djf = exp(lower) * 90) %>%
	mutate(upper_mm_djf = exp(upper) * 90)

### Create an object to plot the rug plot for elevation
rug_df <- djf_df %>%
	group_by(id) %>%
	summarize(elevation = median(elevation, na.rm=TRUE))

p <- ggplot(plot_df, aes(x=elevation)) %>%
	+ geom_hline(yintercept = 0, linetype = "dashed") %>%
	+ geom_ribbon(aes(y= elev_partial_int_mm_djf, ymin = lower_mm_djf, ymax = upper_mm_djf), alpha = 0.9, fill = "grey80")%>%
	+ geom_line(aes(y= elev_partial_int_mm_djf)) %>%
	+ geom_rug(data = rug_df, sides="b", alpha = 0.2)  %>%
	+ scale_x_continuous(name = "Elevation (m)") %>%
	+ scale_y_continuous(name = "Elevation Effect \nAdd. DJF Precipitation (mm)", breaks = seq(0,1000,100)) %>%
	+ theme_bw(9)

#p

### Save plot
ggsave(file.path(write_figures_path, "elevation_term.png"), p,  height = 3, width = 4, dpi = 600)
ggsave(file.path(write_figures_path, "elevation_term.pdf"), p,  height = 3, width = 4)
ggsave(file.path(write_figures_path, "elevation_term.svg"), p,  height = 3, width = 4)


p <- ggplot(plot_df, aes(x=elevation)) %>%
	+ geom_hline(yintercept = 0, linetype = "dotted", colour = "grey30") %>%
	+ geom_ribbon(aes(y= elev_partial_int, ymin = lower, ymax = upper), alpha = 0.9, fill = "grey80")%>%
	+ geom_line(aes(y= elev_partial_int)) %>%
	+ geom_rug(data = rug_df, sides="b", alpha = 0.2)  %>%
	+ scale_x_continuous(name = "Elevation (m)") %>%
	+ scale_y_continuous(name = "Elevation Effect \n(log space)") %>%
	+ theme_bw(9)
#p

### Save plot
ggsave(file.path(write_figures_path, "elevation_term_log.png"), p,  height = 3, width = 4, dpi = 600)
ggsave(file.path(write_figures_path, "elevation_term_log.pdf"), p,  height = 3, width = 4)
ggsave(file.path(write_figures_path, "elevation_term_log.svg"), p,  height = 3, width = 4)


p <- draw(model_climatol, select = 2, residuals = TRUE)
#p
ggsave(file.path(write_figures_path, "elevation_term_gratia.png"), p,  height = 4, width = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "elevation_term_gratia.pdf"), p,  height = 4, width = 6.5)
ggsave(file.path(write_figures_path, "elevation_term_gratia.svg"), p,  height = 4, width = 6.5)


cat("Finished Model 1 Climatology")


####################
### Plot long term trend term
####################
wyo_ids <- c("USC00104588", "USC00108937")

### Extract the terms
plot_wyo_ts <- djf_df %>%
	filter(id %in% wyo_ids) %>%
	mutate(class = "obs") %>%
	mutate(group = paste0(class, "_", id)) %>%
	as.data.frame()

plot_wyo_pred <-  plot_wyo_ts %>%
		select(id, x_proj, y_proj, elevation) %>%
		unique() %>%
		expand_grid(year = seq(1900, 2020, by =1))

partial_df<- predict(model_climatol, newdata = plot_wyo_pred, type = "link", se.fit=TRUE)

plot_wyo_pred$pred_terms <- exp(partial_df$fit)
plot_wyo_pred$lower <- exp(partial_df$fit - 1.96*partial_df$se.fit)
plot_wyo_pred$upper <- exp(partial_df$fit + 1.96*partial_df$se.fit)


	plot_wyo_pred <- plot_wyo_pred %>%
		mutate(pred_terms = exp(partial_df$fit)*90) %>%
		mutate(lower = exp(partial_df$fit - 1.96*partial_df$se.fit)*90) %>%
		mutate(upper = exp(partial_df$fit + 1.96*partial_df$se.fit)*90) %>%
		mutate(class = "pred")  %>%
		mutate(group = paste0(class, "_", id)) %>%
		as.data.frame() %>%
		filter(id == "USC00104588")

	p <- ggplot(plot_wyo_ts, aes(x = year, colour = class, group = group)) %>%
		+ geom_ribbon(data = plot_wyo_pred , aes(ymin=lower, ymax = upper, colour = class, fill = class),  alpha =0.05) %>%
		+ geom_line(aes(y = djf_mm_daily * 90)) %>%
		+ geom_line(data = plot_wyo_pred, aes(y=pred_terms)) %>%
		+ scale_colour_manual(breaks = c("obs", "pred"), values = c("black", "red")) %>%
		+ scale_fill_manual(breaks = c("obs", "pred"), values = c("black", "red")) %>%
		+ scale_y_continuous(name = "Precip (mm)") %>%
		+ scale_x_continuous(name = "Year", breaks = seq(1900, 2050, by = 25)) %>%
		+ theme_classic(9) %>%
			+  theme(legend.position="none")

	ggsave(file.path(write_figures_path,"wyo_ts.png"), p,  height = 4, width = 6.5, dpi = 600)
	ggsave(file.path(write_figures_path,"wyo_ts.svg"), p,  height = 4, width = 6.5)



#mx_ids <- c("MXN00026109", "MXN00026075")
#	mx_ids <-"MXN00026109"
	mx_ids <-"MXN00026075"


	### Extract the terms
	plot_mx_ts <- djf_df %>%
		filter(id %in% mx_ids) %>%
		mutate(class = "obs") %>%
		mutate(group = paste0(class, "_", id)) %>%
		as.data.frame()

	plot_mx_pred <-  plot_mx_ts %>%
		select(id, x_proj, y_proj, elevation) %>%
		unique() %>%
		expand_grid(year = seq(1900, 2020, by =1))


	partial_df<- predict(model_climatol, newdata = plot_mx_pred, type = "link", se.fit=TRUE)


plot_mx_pred$pred_terms <- exp(partial_df$fit)
plot_mx_pred$lower <- exp(partial_df$fit - 1.96*partial_df$se.fit)
plot_mx_pred$upper <- exp(partial_df$fit + 1.96*partial_df$se.fit)


	plot_mx_pred <- plot_mx_pred %>%
		mutate(pred_terms = exp(partial_df$fit)*90) %>%
		mutate(lower = exp(partial_df$fit - 1.96*partial_df$se.fit)*90) %>%
		mutate(upper = exp(partial_df$fit + 1.96*partial_df$se.fit)*90) %>%
		mutate(class = "pred")  %>%
		mutate(group = paste0(class, "_", id)) %>%
		as.data.frame()

	p <- ggplot(plot_mx_ts, aes(x = year, colour = class, group = group)) %>%
		+ geom_ribbon(data = plot_mx_pred, aes(ymin=lower, ymax = upper, colour = class, fill = class),  alpha =0.05) %>%
		+ geom_line(aes(y = djf_mm_daily * 90)) %>%
		+ geom_line(data = plot_mx_pred, aes(y=pred_terms)) %>%
		+ scale_colour_manual(breaks = c("obs", "pred"), values = c("black", "red")) %>%
		+ scale_fill_manual(breaks = c("obs", "pred"), values = c("black", "red")) %>%
		+ scale_y_continuous(name = "Precip (mm)") %>%
		+ scale_x_continuous(name = "Year", breaks = seq(1900, 2050, by = 25)) %>%
		+ theme_classic(9) %>%
			+  theme(legend.position="none")

		ggsave(file.path(write_figures_path,"mx_ts.png"), p,  height = 4, width = 6.5, dpi = 600)
		ggsave(file.path(write_figures_path,"mx_ts.svg"), p,  height = 4, width = 6.5)


####################
###  Plot long term trend spatially
####################

		new_df <- plot_grid_sf

		wmo_base <- new_df %>%
			mutate(year = 1975)

		wmo_base <- wmo_base %>%
				mutate(base = predict(model_climatol, newdata = wmo_base, type = "response")*90) %>%
				as.data.frame()


		for (j in seq(1900, 2020, by = 10)){

		new_df <- new_df %>%
			mutate(year = j)

		### Extract the terms
		plot_df <- new_df %>%
			mutate(pred = predict(model_climatol, newdata = new_df, type = "response")*90)

		plot_df <- plot_df %>%
			left_join(wmo_base %>% select(point, base) , by = "point") %>%
			mutate(diff = pred - base)


		p <- ggplot(plot_df ) %>%
			#+ geom_stars(data = background, alpha = 0.1) %>%
		#	+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = spatial_partial_int_mm_djf)) %>%
			+ geom_sf( alpha = 0.95, aes(fill = pred), colour = NA) %>%
			+ geom_sf(data = rivers, colour = "#4a80f5", alpha = 0.15) %>%
			+ geom_sf(data = lakes, fill = "#9bbff4", alpha = 0.15) %>%
			+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
			+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
			+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
			#+ geom_sf(size = 0.1, alpha = 0.5) %>%
			#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
			#+ scale_fill_distiller(name = "Spatial Effect\n Add. DJF \nPrecipitation (mm)", type = "div", palette = "RdBu", direction = 1, limits = c(-160,160)) %>%
			+ scale_fill_viridis(name = "Spatial Effect\n Add. DJF \nPrecipitation (mm)", trans = "log", option =  "inferno" , direction = -1, limits = c(4, 750), breaks = my_breaks) %>%
			+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
			+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
			+ xlab("Longitude") %>%
			+ ylab("Latitude") %>%
			+ theme_bw(8) %>%
			+ theme(legend.position="bottom") %>%
			+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

		ggsave(file.path(write_figures_path, paste0("byyear_",j,".png")), p,  height = 4, width = 6.5, dpi = 600)



		p <- ggplot(plot_df ) %>%
			#+ geom_stars(data = background, alpha = 0.1) %>%
		#	+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = spatial_partial_int_mm_djf)) %>%
			+ geom_sf( alpha = 0.95, aes(fill = diff), colour = NA) %>%
			+ geom_sf(data = rivers, colour = "#4a80f5", alpha = 0.15) %>%
			+ geom_sf(data = lakes, fill = "#9bbff4", alpha = 0.15) %>%
			+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
			+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
			+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
			#+ geom_sf(size = 0.1, alpha = 0.5) %>%
			#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
			+ scale_fill_distiller(name = "Diff \nPrecipitation (mm)", type = "div", palette = "RdBu", direction = 1, limits = c(-120,120), , oob = squish) %>%
			#+ scale_fill_viridis(name = "Spatial Effect\n Add. DJF \nPrecipitation (mm)", trans = "log", option =  "inferno" , direction = -1, limits = c(4, 750), breaks = my_breaks) %>%
			+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
			+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
			+ xlab("Longitude") %>%
			+ ylab("Latitude") %>%
			+ theme_bw(8) %>%
			+ theme(legend.position="bottom") %>%
			+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

		ggsave(file.path(write_figures_path, paste0("byyear_",j,"_diff.png")), p,  height = 4, width = 6.5, dpi = 600)


		}


		####################
		###  Plot long term trend north to south
		####################
		plot_transect <- data.frame(lon = -110, lat = seq(46, 25, by = -1))
		plot_transect <- plot_transect %>% mutate(id = seq(1, dim(plot_transect)[1]))
		### Create an sf object of gauge locations
		plot_transect_sf = st_as_sf(plot_transect, coords = c("lon", "lat"), crs = 4326 )

		### Reproject into USA CONUS Albert conic
		plot_transect_sf_proj <- plot_transect_sf %>%
			st_transform(5070)

		proj_coord <- st_coordinates(plot_transect_sf_proj)
		proj_coord <- data.frame(id = plot_transect_sf_proj$id, x_proj = proj_coord[,1], y_proj = proj_coord[,2])

		plot_transect <- plot_transect %>%
			left_join(proj_coord, by = "id")

		plot_df <- expand_grid(plot_transect, data.frame(year = seq(1900,2020, by = 1)))
		plot_df$elevation <- 1000

	partial_df<- predict(model_climatol, newdata = plot_df, type = "iterms")
	intercept <- coef(model_climatol)[[1]][[1]]

	plot_df <- plot_df %>%
		mutate(year_pred_logspace = partial_df[,3]) %>%
		mutate(year_pred = exp(year_pred_logspace + intercept))

		p <- ggplot(plot_df, aes(x = year, y = year_pred_logspace, colour = lat, group = id)) %>%
			+ geom_line() %>%
			+ scale_colour_viridis() %>%
				+ scale_y_continuous(name = "Precip (mm)") %>%
			+ scale_x_continuous(name = "Year", breaks = seq(1900, 2050, by = 25)) %>%
			+ theme_classic(9)

			p <- ggplot(plot_df, aes(x = year, y = year_pred, colour = lat, group = id)) %>%
				+ geom_line() %>%
				+ scale_colour_viridis() %>%
					+ scale_y_continuous(name = "Precip (mm)") %>%
				+ scale_x_continuous(name = "Year", breaks = seq(1900, 2050, by = 25)) %>%
				+ theme_classic(9)

				new_df <- plot_grid_sf %>%
				 filter(lon >= -110.05 & lon <= -109.95)

				 plot_df <- expand_grid(st_drop_geometry(new_df), data.frame(year = seq(1850,2020, by = 1)))
		 	 baseline_df <- st_drop_geometry(new_df) %>% mutate(year = 1975)

				response_df <- predict(model_climatol, newdata = plot_df, type = "response")
					response_baseline <- baseline_df %>%
					 		mutate(baseline = predict(model_climatol, newdata = baseline_df, type = "response"))
		 	partial_df<- predict(model_climatol, newdata = plot_df, type = "terms")
		 	intercept <- coef(model_climatol)[[1]][[1]]

		 	plot_df <- plot_df %>%
		 		mutate(year_pred_logspace = partial_df[,3]) %>%
		 		mutate(year_pred = exp(year_pred_logspace + intercept)) %>%
				mutate(response = response_df)

				plot_df <- plot_df %>%
					left_join(response_baseline %>% select(point, baseline), by = "point") 	%>%
					mutate(diff_baseline = response - baseline)

		 		p <- ggplot(plot_df, aes(x = year, y = year_pred_logspace, colour = lat, group = point)) %>%
		 			+ geom_line() %>%
		 			+ scale_colour_viridis() %>%
		 				+ scale_y_continuous(name = "Precip (mm)") %>%
		 			+ scale_x_continuous(name = "Year", breaks = seq(1900, 2050, by = 25)) %>%
		 			+ theme_classic(9)

			p
			dev.new()
		 			p <- ggplot(plot_df, aes(x = year, y = year_pred, colour = lat, group = point)) %>%
		 				+ geom_line() %>%
		 				+ scale_colour_viridis() %>%
		 					+ scale_y_continuous(name = "Precip (mm)") %>%
		 				+ scale_x_continuous(name = "Year", breaks = seq(1900, 2050, by = 25)) %>%
		 				+ theme_classic(9)

						p
						dev.new()
			 			p <- ggplot(plot_df, aes(x = year, y = response, colour = lat, group = point)) %>%
			 				+ geom_line() %>%
			 				+ scale_colour_viridis() %>%
			 					+ scale_y_continuous(name = "Precip (mm)") %>%
			 				+ scale_x_continuous(name = "Year", breaks = seq(1900, 2050, by = 25)) %>%
			 				+ theme_classic(9)
							p
							dev.new()
							p <- ggplot(plot_df, aes(x = year, y = diff_baseline * 90, colour = lat, group = point)) %>%
				 				+ geom_line() %>%
				 				+ scale_colour_viridis() %>%
				 					+ scale_y_continuous(name = "Precip difference (mm)") %>%
				 				+ scale_x_continuous(name = "Year", breaks = seq(1850, 2050, by = 25)) %>%
				 				+ theme_classic(9)
								p
								dev.new()

								p <- ggplot(plot_df %>% filter(lat %in% unique(plot_df$lat)[seq(1, 225, by = 10)]), aes(x = year, y = diff_baseline * 90, colour = lat, group = point)) %>%
									+ geom_line() %>%
									+ scale_colour_viridis(option = "plasma", name = "Latitude") %>%
										+ scale_y_continuous(name = "Precip difference (mm)") %>%
									+ scale_x_continuous(name = "Year", breaks = seq(1850, 2050, by = 25)) %>%
									+ theme_classic(9)
									p

						p <- 	p + coord_cartesian(xlim = c(1895, 2020), ylim = c(-100,100))

						ggsave(file.path(write_figures_path, "year_diff.png"), p,  height = 4, width = 6.5, dpi = 600)
						ggsave(file.path(write_figures_path,"year_diff.svg"), p,  height = 4, width = 6.5)



####################
### The END
####################
