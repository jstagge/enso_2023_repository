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

cat("Starting 04_model_2_enso")

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

#require(lfstat)
require(scales)

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
model_folder <- file.path(output_path, "model_2_enso_trend")
dir.create(model_folder, recursive=TRUE, showWarnings = FALSE)

### Set up figure folder
write_figures_path <- file.path(write_figures_path, "model_2_enso_trend")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)

###########################################################################
###  Read in model 1 - Climatology
###########################################################################
model_climatol <- readRDS(file.path(output_path, "model_1_climatol/model_climatol.rds"))


###########################################################################
###  Extract Climatology Model coefficients for use as initial values
###########################################################################
### Extract coefficients to use as initial values
#coef_climatol <- coef(model_climatol)
#head(coef_climatol)
#tail(coef_climatol)

### Extract distribution information
#power_climatol <- as.numeric(substr(summary(model_climatol)$family$family, 11, 14))
#power_climatol

#phi_climatol <- summary(model_climatol)$dispersion
#phi_climatol


###########################################################################
###  Predict values from climatology to use for a potential anomaly model
###########################################################################
### I don't love this approach, it will use a normal distribution, when realistically it doesn't make sense (there is a lower limit in some cases, but it's not zero)
#djf_df$climatol_pred <- predict(model_climatol, newdata = djf_df, type = "response")

#djf_df <- djf_df %>%
#	mutate(anom = djf_mm_daily - climatol_pred)

###########################################################################
###  Clean up a potential problem site
###########################################################################
djf_df <- djf_df %>%
	filter(id != "MXN00008118") %>%
	mutate(id = as.factor(id))


###########################################################################
###  Loop through and test best ENSO model
###########################################################################


###########################################################################
###  Create model 2 - ENSO
###########################################################################
### Clear off NAs
djf_df <- djf_df %>%
	drop_na(djf_mm_daily, x_proj, y_proj, elevation, enso)

#+ s(id, bs = 're')
#### Create ENSO model
### Fit the model
model_enso <-   bam(djf_mm_daily  ~ te(x_proj,y_proj, bs = "tp", k = c(15,15))+
 						  s(elevation, bs = "ad", k = 15) +
						  te(x_proj, y_proj, enso, d=c(2,1), k = c(15,2), bs = c("tp", "cr")) +
						 te(x_proj, y_proj, year, d=c(2,1), k = c(8,3), bs = c("tp", "tp")),
						  data = djf_df ,
						 family =  tw(),
						 samfrac = 0.1,
						 # discrete=TRUE,
						  nthreads=8)

#model_enso <- readRDS(file.path(model_folder, "model_enso.rds"))
#model_enso_df <- readRDS(file.path(model_folder, "model_enso_df.rds"))

### Save the model
model_enso_df <- djf_df
model_enso_df$pred <- fitted(model_enso, new_data = model_enso_df)
model_enso_df$resid <- resid(model_enso, newdata = model_enso_df)

### Save the model objects
saveRDS(model_enso, file = file.path(model_folder, "model_enso.rds"))
saveRDS(model_enso_df, file = file.path(model_folder, "model_enso_df.rds"))

### Save some residual information
resid_model <- resid(model_enso)

resid_stats <- data.frame(mean_resid = mean(resid_model, na.rm=TRUE),
	median_resid = median(resid_model, na.rm=TRUE),
	mae = mean(abs(resid_model), na.rm=TRUE),
	rmse = sqrt(mean((resid_model)^2)),
	aic = AIC(model_enso),
	bic = BIC(model_enso),
	r_square = summary(model_enso)$r.sq
	)

write_csv(resid_stats, file.path(model_folder, "model2_resid_stats.csv"))

#quantile(resid_model, seq(0.1, 0.9, 0.1))

AIC(model_enso)
BIC(model_enso)
###########################################################################
###  Plot the results
###########################################################################

p <- draw(model_enso) + theme_bw(9)
ggsave(file.path(write_figures_path, "draw_model.png"), p,  width = 10, height = 6.5, dpi = 600)

p <- appraise(model_enso) + theme_bw(9)
ggsave(file.path(write_figures_path, "appraise_model.png"), p,  width = 10, height = 6.5, dpi = 600)



###########################################################################
###  Plot Model Skill
###########################################################################
### Plot predicted versus observed
plot_df <- model_enso_df

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
###  Plot the spatial only term
###########################################################################
require(raster)

### Create object to hold results
#plot_df <- expand_grid(x_proj = x_seq, y_proj = y_seq, elevation = 1300, id = "USS0009E11S", enso = 0)
plot_df <- plot_grid_sf %>%
	mutate(elevation = 1300, id = "USS0009E11S", enso = 0, year = 1975)

### Extract the terms
partial_df<- predict(model_enso, newdata = plot_df, type = "terms")

plot_df$spatial_partial <- partial_df[,1]
plot_df$spatial_partial_int <- (partial_df[,1]+coef(model_enso)[1])

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
	+ scale_fill_distiller(name = "Spatial Effect\n (Log space)", type = "div", palette = "RdBu", direction = 1, limits = c(-2,2), oob = squish) %>%
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
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = spatial_partial)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = spatial_partial), colour = NA) %>%
	#+ geom_sf(data = rivers, colour = "#4a80f5", alpha = 0.15) %>%
	#+ geom_sf(data = lakes, fill = "#9bbff4", alpha = 0.15) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	#+ geom_sf(size = 0.1, alpha = 0.5) %>%
	#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
	+ scale_fill_distiller(name = "Spatial Effect\n (Log space)", type = "div", palette = "RdBu", direction = 1, limits = c(-2,2), oob=scales::squish) %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#p

### Save plot
ggsave(file.path(write_figures_path, "map_spatial_term_log_noriver.png"), p,  width = 4, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "map_spatial_term_log_noriver.pdf"), p,  width = 4, height = 6.5)
ggsave(file.path(write_figures_path, "map_spatial_term_log_noriver.svg"), p,  width = 4, height = 6.5)


### Plot break
my_breaks <- c(2, 5, 10, 20, 50, 100, 200, 500)

### Create plot
p <- ggplot(data = plot_df) %>%
	#+ geom_stars(data = background, alpha = 0.1) %>%
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = spatial_partial_int_mm_djf)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = spatial_partial_int_mm_djf), colour = NA) %>%
	#+ geom_sf(data = rivers, colour = "#4a80f5", alpha = 0.15) %>%
	#+ geom_sf(data = lakes, fill = "#9bbff4", alpha = 0.15) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	#+ geom_sf(size = 0.1, alpha = 0.5) %>%
	#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
	#+ scale_fill_distiller(name = "Spatial Effect\n Add. DJF \nPrecipitation (mm)", type = "div", palette = "RdBu", direction = 1, limits = c(-160,160)) %>%
	+ scale_fill_viridis(name = "Spatial Effect\n Add. DJF \nPrecipitation (mm)", trans = "log", option =  "inferno" , direction = -1, limits = c(4, 750), breaks = my_breaks, oob=scales::squish) %>%
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


p <- draw(model_enso, select = 1, residuals = FALSE) + theme_bw(9)
#p
ggsave(file.path(write_figures_path, "spatial_term_gratia.png"), p,  height = 4, width = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "spatial_term_gratia.pdf"), p,  height = 4, width = 6.5)
ggsave(file.path(write_figures_path, "spatial_term_gratia.svg"), p,  height = 4, width = 6.5)



###########################################################################
###  Plot the elevation term
###########################################################################
### Separate terms
plot_df <- expand_grid(x_proj = median(djf_df$x_proj), y_proj = median(djf_df$y_proj), longitude = -108, latitude = 40, elevation = seq(0,3500), id = "USS0009E11S", enso = 0, year = 1975)

### Create an object to plot the rug plot for elevation
rug_df <- djf_df %>%
	group_by(id) %>%
	summarize(elevation = median(elevation, na.rm=TRUE))

#partial_df<- predict(model_climatol, newdata = plot_df, type = "response", se.fit=TRUE)
partial_df<- predict(model_enso, newdata = plot_df, type = "terms", se.fit=TRUE)

plot_df$elev_partial <- partial_df$fit[,2]
plot_df$elev_partial_int <- partial_df$fit[,2]+coef(model_enso)[1]
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

### Create elevation plot
p <- ggplot(plot_df, aes(x=elevation)) %>%
	+ geom_hline(yintercept = 0, linetype = "dashed") %>%
	+ geom_ribbon(aes(y= elev_partial_int_mm_djf, ymin = lower_mm_djf, ymax = upper_mm_djf), alpha = 0.9, fill = "grey80")%>%
	+ geom_line(aes(y= elev_partial_int_mm_djf)) %>%
	+ geom_rug(data = rug_df, sides="b", alpha = 0.2)  %>%
	+ scale_x_continuous(name = "Elevation (m)") %>%
	+ scale_y_continuous(name = "Elevation Effect \nAdd. DJF Precipitation (mm)", breaks = seq(0,1000,100)) %>%
	+ theme_bw(9)

### Save plot
ggsave(file.path(write_figures_path, "elevation_term.png"), p,  height = 3, width = 4, dpi = 600)
ggsave(file.path(write_figures_path, "elevation_term.pdf"), p,  height = 3, width = 4)
ggsave(file.path(write_figures_path, "elevation_term.svg"), p,  height = 3, width = 4)


### Create elevation plot in log space
p <- ggplot(plot_df, aes(x=elevation)) %>%
	+ geom_hline(yintercept = 0, linetype = "dashed") %>%
	+ geom_ribbon(aes(y= elev_partial_int, ymin = lower, ymax = upper), alpha = 0.9, fill = "grey80")%>%
	+ geom_line(aes(y= elev_partial_int)) %>%
	+ geom_rug(data = djf_df, sides="b")  %>%
	+ scale_x_continuous(name = "Elevation (m)") %>%
	+ scale_y_continuous(name = "Elevation Effect \n(log space)") %>%
	+ theme_bw(9)

### Save plot
ggsave(file.path(write_figures_path, "elevation_term_log.png"), p,  height = 3, width = 4, dpi = 600)
ggsave(file.path(write_figures_path, "elevation_term_log.pdf"), p,  height = 3, width = 4)
ggsave(file.path(write_figures_path, "elevation_term_log.svg"), p,  height = 3, width = 4)


p <- draw(model_enso, select = 2, residuals = TRUE)
#p
ggsave(file.path(write_figures_path, "elevation_term_gratia.png"), p,  height = 3, width = 4, dpi = 600)
ggsave(file.path(write_figures_path, "elevation_term_gratia.pdf"), p,  height = 3, width = 4)
ggsave(file.path(write_figures_path, "elevation_term_gratia.svg"), p,  height = 3, width = 4 )


###########################################################################
###  Plot ENSO effect
###########################################################################

enso_list <- seq(-2,3,0.5)

for (j in seq(1, length(enso_list))){
	plot_temp <- plot_grid_sf %>%
		mutate(elevation = 1300, id = "USS0009E11S", year = 1975) %>%
		mutate(enso = enso_list[[j]])

	if(j == 1){
		plot_df <- plot_temp
	} else {
		plot_df <- plot_df %>%
			bind_rows(plot_temp)
	}
}

partial_df<- predict(model_enso, newdata = plot_df, type = "terms")
plot_df$enso_partial <- partial_df[,3]
plot_df$enso_partial_int <- (partial_df[,3]+coef(model_enso)[1])

plot_df <- plot_df %>%
	mutate(enso_label = paste0("MEI = ", enso))%>%
	mutate(enso_label = factor(enso_label, levels = paste0("MEI = ", enso_list)))

### Convert back to 3 month precip. Remember tweedie uses a log link
plot_df <- plot_df %>%
	mutate(enso_partial_mm_djf = exp(enso_partial) * 90) %>%
	mutate(enso_partial_int_mm_djf = exp(enso_partial_int) * 90)

### Create plot
p <- ggplot(data = plot_df) %>%
	#+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	#+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	#+ geom_stars(data = background, alpha = 0.1) %>%
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = enso_partial)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = enso_partial), colour = NA) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5, size = 0.15) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5, size = 0.15) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5, size = 0.15) %>%
	#+ geom_sf(size = 0.1, alpha = 0.5) %>%
	#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
	+ facet_grid(.~enso_label) %>%
#	+ scale_fill_distiller(name = "ENSO Effect\n (Log space)", type = "div", palette = "RdBu", direction = 1, limits = c(-1.8,1.8), oob=scales::squish) %>%
	+ scale_fill_gradientn(name = "ENSO Effect\n (Log space)", colors = ipccPrec(11), limits = c(-1.8,1.8), oob=scales::squish) %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(7) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Save plot
ggsave(file.path(write_figures_path, "enso_term_log.png"), p,  width = 8, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "enso_term_log.pdf"), p,  width = 8, height = 6.5)
ggsave(file.path(write_figures_path, "enso_term_log.svg"), p,  width = 8, height = 6.5)


### Fewer breaks
enso_list <- seq(-2,3,1)

plot_sub <- plot_df %>%
	filter(enso %in% enso_list)%>%
	mutate(enso_label = factor(enso_label, levels = paste0("MEI = ", enso_list)))

p <- ggplot(data = plot_sub) %>%
	#+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	#+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	#+ geom_stars(data = background, alpha = 0.1) %>%
	#+ geom_raster(alpha = 0.95, aes(x = x_proj, y=y_proj, fill = enso_partial)) %>%
	+ geom_sf(alpha = 0.95, aes(fill = enso_partial), colour = NA) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5, size = 0.15) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5, size = 0.15) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5, size = 0.15) %>%
	#+ geom_sf(size = 0.1, alpha = 0.5) %>%
	#+ scale_fill_gradient2(name = "ENSO Partial\nLa Nina\nENSO = -2") %>%
	+ facet_grid(.~enso_label) %>%
	#+ scale_fill_distiller(name = "ENSO Effect\n (Log space)", type = "div", palette = "RdBu", direction = 1, limits = c(-1.8,1.8), oob=scales::squish) %>%
	+ scale_fill_gradientn(name = "ENSO Effect\n (Log space)", colors = ipccPrec(11), limits = c(-1.8,1.8), oob=scales::squish) %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(7) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Save plot
ggsave(file.path(write_figures_path, "enso_term_log_subset.png"), p,  width = 8, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "enso_term_log_subset.pdf"), p,  width = 8, height = 6.5)
ggsave(file.path(write_figures_path, "enso_term_log_subset.svg"), p,  width = 8, height = 6.5)





###########################################################################
###  Plot ENSO effect in precipitation
###########################################################################
### Instead, do each grid with its elevation and lat / lon for ENSO = 0
### Then subtract the other ENSOs from it

### Create estimate for enso = 0
plot_base <- plot_grid_sf %>%
	mutate(enso = 0, year = 1975)

plot_base <- plot_base %>%
	mutate(baseline = predict(model_enso, newdata = plot_base, type = "response")) %>%
	as.data.frame()

#fpp <- predict(model_enso, newdata = plot_base, type = "terms")

### Create estimates for enso != 0
enso_list <- seq(-2,3,0.5)

for (j in seq(1, length(enso_list))){
	enso_j <- enso_list[[j]]

	plot_temp <- plot_grid_sf %>%
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
	mutate(pred = predict(model_enso, newdata = plot_enso, type = "response")) %>%
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
	+ scale_fill_gradientn(name = "ENSO Effect\n (mm)", colors = ipccPrec(11), limits = c(-plot_lims,plot_lims), oob = squish) %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p


### Save plot
ggsave(file.path(write_figures_path, "enso_effect_mm.png"), p,  width = 12, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "enso_effect_mm.pdf"), p,  width = 12, height = 6.5)
ggsave(file.path(write_figures_path, "enso_effect_mm.svg"), p,  width = 12, height = 6.5)



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
	#+ scale_fill_distiller(name = "ENSO Effect\n (mm)", type = "div", palette = "RdBu", direction = 1, limits = c(-plot_lims,plot_lims), oob = squish) %>%
	+ scale_fill_gradientn(name = "ENSO Effect\n (mm)", colors = ipccPrec(11), limits = c(-plot_lims,plot_lims), oob = squish) %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Save plot
ggsave(file.path(write_figures_path, "enso_effect_mm_subset.png"), p,  width = 12, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "enso_effect_mm_subset.pdf"), p,  width = 12, height = 6.5)
ggsave(file.path(write_figures_path, "enso_effect_mm_subset.svg"), p,  width = 12, height = 6.5)




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
	select(id, x_proj, y_proj, elevation, enso, year) %>%
	mutate(enso =0)


partial_df<- predict(model_enso, newdata = plot_wyo_pred, type = "link", se.fit=TRUE)

plot_wyo_pred$pred_terms <- exp(partial_df$fit)
plot_wyo_pred$lower <- exp(partial_df$fit - 1.96*partial_df$se.fit)
plot_wyo_pred$upper <- exp(partial_df$fit + 1.96*partial_df$se.fit)


	plot_wyo_pred <- plot_wyo_pred %>%
		mutate(pred_terms = exp(partial_df$fit)*90) %>%
		mutate(lower = exp(partial_df$fit - 1.96*partial_df$se.fit)*90) %>%
		mutate(upper = exp(partial_df$fit + 1.96*partial_df$se.fit)*90) %>%
		mutate(class = "pred")  %>%
		mutate(group = paste0(class, "_", id)) %>%
		as.data.frame() #%>%
	#	filter(id == "USC00104588")

	p <- ggplot(plot_wyo_ts, aes(x = year, colour = class, group = group)) %>%
		+ geom_ribbon(data = plot_wyo_pred , aes(ymin=lower, ymax = upper, colour = class, fill = class),  alpha =0.15, colour = NA) %>%
		+ geom_line(aes(y = djf_mm_daily * 90)) %>%
		+ geom_line(data = plot_wyo_pred, aes(y=pred_terms)) %>%
		+ scale_colour_manual(breaks = c("obs", "pred"), values = c("black", "red")) %>%
		+ scale_fill_manual(breaks = c("obs", "pred"), values = c("black", "red")) %>%
		+ scale_y_continuous(name = "DJF Precip (mm)") %>%
		+ scale_x_continuous(name = "Year", breaks = seq(1900, 2050, by = 25)) %>%
		+ theme_classic(9) %>%
		+ theme(legend.position = "none")


	ggsave(file.path(write_figures_path,"wyo_ts.png"), p,  height = 1.7, width = 3.6, dpi = 600)
	ggsave(file.path(write_figures_path,"wyo_ts.svg"), p,  height = 1.7, width = 3.6)



	### Extract the terms
	plot_wyo_ts <- djf_df %>%
		filter(id %in% wyo_ids) %>%
		as.data.frame()

	plot_wyo_pred <-  plot_wyo_ts %>%
			select(id, x_proj, y_proj, elevation) %>%
			unique() %>%
			filter(id == "USC00104588") %>%
			expand_grid(year = seq(1900, 2020, by = 1), enso = 0)

	partial_df<- predict(model_enso, newdata = plot_wyo_pred, type = "link", se.fit=TRUE)

	plot_wyo_pred$pred_terms <- exp(partial_df$fit)
	plot_wyo_pred$lower <- exp(partial_df$fit - 1.96*partial_df$se.fit)
	plot_wyo_pred$upper <- exp(partial_df$fit + 1.96*partial_df$se.fit)

		p <- ggplot(plot_wyo_pred, aes(x = year)) %>%
			+ geom_ribbon(data = plot_wyo_pred , aes(ymin=lower, ymax = upper),  alpha =0.05, fill = "grey50", colour = "grey80") %>%
			+ geom_line(aes(y = pred_terms), colour = "grey20") %>%
			+ scale_y_continuous(name = "Precip (mm)") %>%
			+ scale_x_continuous(name = "Year", breaks = seq(1900, 2050, by = 25)) %>%
			+ theme_classic(9)

		ggsave(file.path(write_figures_path,"wyo_time_term.png"), p,  height = 4, width = 6.5, dpi = 600)
		ggsave(file.path(write_figures_path,"wyo_time_term.svg"), p,  height = 4, width = 6.5)


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
			select(id, x_proj, y_proj, elevation, enso, year)


	plot_mx_pred <-  plot_mx_ts %>%
		select(id, x_proj, y_proj, elevation) %>%
		unique() %>%
		expand_grid(year = seq(1900, 2020, by = 1)) %>%
		mutate(enso = 0)

	partial_df<- predict(model_enso, newdata = plot_mx_pred, type = "link", se.fit=TRUE)


plot_mx_pred$pred_terms <- exp(partial_df$fit) *90
plot_mx_pred$lower <- exp(partial_df$fit - 1.96*partial_df$se.fit) *90
plot_mx_pred$upper <- exp(partial_df$fit + 1.96*partial_df$se.fit) *90


	plot_mx_pred <- plot_mx_pred %>%
		mutate(class = "pred")  %>%
		mutate(group = paste0(class, "_", id)) %>%
		as.data.frame()

	p <- ggplot(plot_mx_ts, aes(x = year, colour = class, group = group, linetype = id)) %>%
		+ geom_ribbon(data = plot_mx_pred, aes(ymin=lower, ymax = upper, colour = class, fill = class),  alpha =0.15, colour = NA) %>%
		+ geom_line(aes(y = djf_mm_daily * 90)) %>%
		+ geom_line(data = plot_mx_pred, aes(y=pred_terms)) %>%
		+ scale_colour_manual(breaks = c("obs", "pred"), values = c("black", "red")) %>%
		+ scale_fill_manual(breaks = c("obs", "pred"), values = c("black", "red")) %>%
		+ scale_y_continuous(name = "Precip (mm)") %>%
		+ scale_x_continuous(name = "Year", breaks = seq(1900, 2050, by = 25)) %>%
		+ theme_classic(9) %>%
		+ theme(legend.position = "none")


	ggsave(file.path(write_figures_path,"mx_ts.png"), p,  height = 1.7, width = 3.6, dpi = 600)
	ggsave(file.path(write_figures_path,"mx_ts.svg"), p,  height = 1.7, width = 3.6)


		### Extract the terms
		plot_mx_ts <- djf_df %>%
			filter(id %in% mx_ids) %>%
			as.data.frame()

		plot_mx_pred <-  plot_mx_ts %>%
				select(id, x_proj, y_proj, elevation) %>%
				unique() %>%
				filter(id == "MXN00026075") %>%
				expand_grid(year = seq(1900, 2020, by = 1), enso = 0)

		partial_df<- predict(model_enso, newdata = plot_mx_pred, type = "link", se.fit=TRUE)

		plot_mx_pred$pred_terms <- exp(partial_df$fit) *90
		plot_mx_pred$lower <- exp(partial_df$fit - 1.96*partial_df$se.fit) *90
		plot_mx_pred$upper <- exp(partial_df$fit + 1.96*partial_df$se.fit) *90

			p <- ggplot(plot_mx_pred, aes(x = year)) %>%
				+ geom_ribbon(data = plot_mx_pred , aes(ymin=lower, ymax = upper),  alpha =0.05, fill = "grey50", colour = "grey80") %>%
				+ geom_line(aes(y = pred_terms), colour = "grey20") %>%
				+ scale_y_continuous(name = "Precip (mm)") %>%
				+ scale_x_continuous(name = "Year", breaks = seq(1900, 2050, by = 25)) %>%
				+ theme_classic(9)

			ggsave(file.path(write_figures_path,"mx_time_term.png"), p,  height = 4, width = 6.5, dpi = 600)
			ggsave(file.path(write_figures_path,"mx_time_term.svg"), p,  height = 4, width = 6.5)



####################
###  Plot long term trend spatially
####################

		new_df <- plot_grid_sf

		wmo_base <- new_df %>%
			mutate(year = 1975, enso = 0)

		wmo_base <- wmo_base %>%
				mutate(base = predict(model_enso, newdata = wmo_base, type = "response")*90) %>%
				as.data.frame()


		for (j in seq(1900, 2020, by = 10)){

		new_df <- new_df %>%
			mutate(year = j, enso = 0)

		### Extract the terms
		plot_df <- new_df %>%
			mutate(pred = predict(model_enso, newdata = new_df, type = "response")*90)

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

						new_df <- plot_grid_sf %>%
						 filter(lon >= -110.05 & lon <= -109.95)

						 plot_df <- expand_grid(st_drop_geometry(new_df), data.frame(year = seq(1850,2020, by = 1))) %>% mutate(enso = 0)
				 	 baseline_df <- st_drop_geometry(new_df) %>% mutate(year = 1975) %>% mutate(enso = 0)

						response_df <- predict(model_enso, newdata = plot_df, type = "response")
							response_baseline <- baseline_df %>%
							 		mutate(baseline = predict(model_enso, newdata = baseline_df, type = "response"))
				 	partial_df<- predict(model_enso, newdata = plot_df, type = "terms")
				 	intercept <- coef(model_enso)[[1]][[1]]

				 	plot_df <- plot_df %>%
				 		mutate(year_pred_logspace = partial_df[,3]) %>%
				 		mutate(year_pred = exp(year_pred_logspace + intercept)) %>%
						mutate(response = response_df)

						plot_df <- plot_df %>%
							left_join(response_baseline %>% select(point, baseline), by = "point") 	%>%
							mutate(diff_baseline = response - baseline)

										p <- ggplot(plot_df %>% filter(lat %in% unique(plot_df$lat)[seq(1, 225, by = 10)]), aes(x = year, y = diff_baseline * 90, colour = lat, group = point)) %>%
											+ geom_line() %>%
											+ scale_colour_viridis(option = "plasma", name = "Latitude") %>%
												+ scale_y_continuous(name = "Precip difference (mm)") %>%
											+ scale_x_continuous(name = "Year", breaks = seq(1850, 2050, by = 25)) %>%
											+ theme_classic(9)
											p

								p <- 	p + coord_cartesian(xlim = c(1895, 2020), ylim = c(-100,100))

								ggsave(file.path(write_figures_path, "year_diff.png"), p,  height = 3, width = 6.5, dpi = 600)
								ggsave(file.path(write_figures_path,"year_diff.svg"), p,  height = 3, width = 6.5)



####################
### The END
####################
