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
#require(rnaturalearth)
#dyn.load("/usr/local/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
#dyn.load("/usr/local/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
#require(sf)
#library(maps)


#require(rnoaa)
require(tidyverse)
#library(sf)
#library(stars)
#library(rnaturalearth)
#library(rnaturalearthhires)
require(viridis)
require(raster)
#require(lfstat)
require(mgcv)
#require(ggthemes)
#require(raster)
select <- dplyr::select

#require(mgcv)
require(gratia)


###########################################################################
###  Set output folder
###########################################################################
write_output_path <- file.path(output_path, "gauge_data")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

write_figures_path <- file.path(write_figures_path, "gauge_data")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)



###########################################################################
###  Read in data
###########################################################################
### Read in all the DJF time series
djf_all <- readRDS(file.path(write_output_path, "djf_all.rds"))

### Read in station info
stations <- readRDS(file.path(write_output_path, "all_stations.rds"))

###########################################################################
###  Read in ENSO
###########################################################################
mei_ext <- read_tsv(file.path(data_path, "mei_ext/mei_ext.txt"))
mei_orig <- read_tsv(file.path(data_path, "mei_orig/mei_orig.txt"))
mei_v2 <- read_tsv(file.path(data_path, "mei_v2/mei_v2.txt"))

enso_df <- mei_ext %>%
	filter(YEAR <= 1949) %>%
	bind_rows(mei_orig)

head(enso_df)

enso_dj <- enso_df %>%
	select(YEAR, DECJAN) %>%
	rename(year = YEAR) %>%
	rename(enso = DECJAN)

head(enso_dj)



##########################################################################
###  Calculate differences
###########################################################################
check_diff <- mei_ext %>%
	select(YEAR, DECJAN) %>%
	inner_join(mei_orig %>% select(YEAR, DECJAN), by = "YEAR")

hydroGOF::gof(check_diff$DECJAN.y,check_diff$DECJAN.x)

##########################################################################
###  Describe the ENSO values
###########################################################################
plot_df <- mei_ext %>%
		mutate(data = "mei_ext") %>%
		bind_rows(mei_orig %>% mutate(data = "mei_v1")) %>%
		bind_rows(mei_v2 %>% mutate(data = "mei_v2"))

### Plot Dec Jan MEI time series
p <- ggplot(plot_df, aes(x=YEAR, y= DECJAN, colour = data, group = data)) %>%
	+ geom_hline(yintercept = 0, colour = "grey50", alpha = 0.5) %>%
	+ geom_line() %>%
	+ scale_colour_brewer(name = "Source", palette = "Set2", type = "qual") %>%
	+ scale_x_continuous(name = "Year", breaks = seq(1800, 2040, by = 20)) %>%
	+ scale_y_continuous(name = "MEI") %>%
	+ theme_classic(9) %>%
	+ theme(legend.position="bottom")

### Save plot
ggsave(file.path(write_figures_path, "mei_ts_dj_alldata.png"), p,  width = 5, height = 2.5, dpi = 600)
ggsave(file.path(write_figures_path, "mei_ts_dj_alldata.pdf"), p,  width = 5, height = 2.5)
ggsave(file.path(write_figures_path, "mei_ts_dj_alldata.svg"), p,  width = 5, height = 2.5)


### Plot Dec Jan MEI time series
p <- ggplot(plot_df %>% filter(data != "mei_v2"), aes(x=YEAR, y= DECJAN, colour = data, group = data)) %>%
	+ geom_hline(yintercept = 0, colour = "grey50", alpha = 0.5) %>%
	+ geom_line() %>%
	+ scale_colour_brewer(name = "Source", palette = "Set2", type = "qual") %>%
	+ scale_x_continuous(name = "Year", breaks = seq(1800, 2040, by = 20)) %>%
	+ scale_y_continuous(name = "MEI") %>%
	+ theme_classic(9) %>%
	+ theme(legend.position="bottom")

### Save plot
ggsave(file.path(write_figures_path, "mei_ts_dj_no_meiv2.png"), p,  width = 5, height = 2.5, dpi = 600)
ggsave(file.path(write_figures_path, "mei_ts_dj_no_meiv2.pdf"), p,  width = 5, height = 2.5)
ggsave(file.path(write_figures_path, "mei_ts_dj_no_meiv2.svg"), p,  width = 5, height = 2.5)


### Merge the MEI Ext and MEI v1
enso_df <- mei_ext %>%
	filter(YEAR <= 1949) %>%
	bind_rows(mei_orig)

head(enso_df)

enso_dj <- enso_df %>%
	select(YEAR, DECJAN) %>%
	rename(year = YEAR) %>%
	rename(enso = DECJAN)

head(enso_dj)

### Plot Dec Jan MEI time series
p <- ggplot(plot_df %>% filter(data != "mei_v2"), aes(x=YEAR, y= DECJAN)) %>%
	+ geom_hline(yintercept = 0, colour = "grey50", alpha = 0.5) %>%
	+ geom_line(aes(colour = data, group = data)) %>%
	+ geom_line(data = enso_dj, aes(x=year, y=enso), colour = "black") %>%
	+ scale_colour_brewer(name = "Source", palette = "Set2", type = "qual") %>%
	+ scale_x_continuous(name = "Year", breaks = seq(1800, 2040, by = 20)) %>%
	+ scale_y_continuous(name = "MEI") %>%
	+ theme_classic(9) %>%
	+ theme(legend.position="bottom")

	p <- ggplot(enso_dj, aes(x=year, y=enso)) %>%
		+ geom_hline(yintercept = 0, colour = "grey50", alpha = 0.5, linetype = "dashed") %>%
		+ geom_line( colour = "black", size = 0.25) %>%
		+ scale_x_continuous(name = "Year", breaks = seq(1800, 2040, by = 20)) %>%
		+ scale_y_continuous(name = "MEI") %>%
		+ theme_classic(9) %>%
		+ theme(legend.position="bottom")

### Save plot
ggsave(file.path(write_figures_path, "mei_ts_dj_comb.png"), p,  width = 5.5, height = 2, dpi = 600)
ggsave(file.path(write_figures_path, "mei_ts_dj_comb.pdf"), p,  width = 5.5, height = 2)
ggsave(file.path(write_figures_path, "mei_ts_dj_comb.svg"), p,  width = 5.5, height = 2)



#### plot with shaded regions
df <- enso_dj
df$grp <- "orig"

new_df <- do.call(what = "rbind",
                  args = sapply(X = 1:(nrow(x = df) -1),
                                FUN = function(i)
                                {
                                  f <- lm(formula = (year ~ enso),
                                          data = df[i:(i + 1),])
                                  if (f$qr$rank < 2)
                                  {
                                    return(NULL)
                                  }
                                  r <- predict(object = f,
                                               newdata = data.frame(enso = 0))
                                  if(df[i,]$year < r & r < df[i+1,]$year)
                                  {
                                    return(data.frame(year = r,
                                                      enso = 0))
                                  } else
                                  {
                                    return(NULL)
                                  }
                                }))

new_df$grp <- "new"

df_mod <- rbind(df, new_df)

p <- ggplot(data = df_mod, aes(x = year, y = enso))  %>%
	+   geom_area(data = subset(x = df_mod, subset = (enso <= 0)), fill = "#542788", alpha = 0.7) %>%
	+  geom_area(data = subset(x = df_mod, subset = (enso >= 0)),  fill = "#b35806", alpha = 0.7) %>%
	+ geom_hline(yintercept = 0, colour = "grey50", alpha = 0.5, linetype = "dashed") %>%
	+ geom_line( data = df, colour = "black", size = 0.25) %>%
	+ scale_x_continuous(name = "Year", breaks = seq(1800, 2040, by = 20)) %>%
	+ scale_y_continuous(name = "MEI") %>%
	+ theme_classic(9) %>%
	+ theme(legend.position="bottom")

p

### Save plot
ggsave(file.path(write_figures_path, "mei_ts_dj_comb_shaded.png"), p,  width = 5.5, height = 2, dpi = 600)
ggsave(file.path(write_figures_path, "mei_ts_dj_comb_shaded.pdf"), p,  width = 5.5, height = 2)
ggsave(file.path(write_figures_path, "mei_ts_dj_comb_shaded.svg"), p,  width = 5.5, height = 2)





### Plot distribution
p <- ggplot(enso_dj, aes(x=enso)) %>%
	+ geom_histogram(binwidth = 0.25, colour = "black", fill = "grey80") %>%
	+ scale_x_continuous(name = "MEI Value") %>%
	+ scale_y_continuous("Count (years)") %>%
	+ theme_classic(9)
p


	### Save plot
	ggsave(file.path(write_figures_path, "mei_dist.png"), p,  width = 3.6, height = 2.55, dpi = 600)
	ggsave(file.path(write_figures_path, "mei_dist.pdf"), p,  width = 3.6, height = 2.55)
	ggsave(file.path(write_figures_path, "mei_dist.svg"), p,  width = 3.6, height = 2.55)


	p <- p +  geom_density(aes(y=0.25 * ..count..), colour = "blue")

		### Save plot
		ggsave(file.path(write_figures_path, "mei_dist_line.png"), p,  width = 3.6, height = 2.55, dpi = 600)
		ggsave(file.path(write_figures_path, "mei_dist_line.pdf"), p,  width = 3.6, height = 2.55)
		ggsave(file.path(write_figures_path, "mei_dist_line.svg"), p,  width = 3.6, height = 2.55)


	#hist(unique(djf_df$enso))

	### Save some MEI information
	mei_stats <- data.frame(mean = mean(enso_dj$enso, na.rm=TRUE),
		median = median(enso_dj$enso, na.rm=TRUE),
		min = min(enso_dj$enso, na.rm=TRUE),
		max = max(enso_dj$enso, na.rm=TRUE),
		perc_lt_neg1 = mean(enso_dj$enso <= -1, na.rm=TRUE),
		perc_gt_pos1 = mean(enso_dj$enso >= 1, na.rm=TRUE))

	write_csv(mei_stats, file.path(write_output_path, "mei_stats.csv"))

	### 17% >= 1
	### 17% <= -1
	### 66% between them



plot_df <- mei_ext %>%
		mutate(data = "mei_ext") %>%
		bind_rows(mei_orig %>% mutate(data = "mei_v1")) %>%
		bind_rows(mei_v2 %>% mutate(data = "mei_v2"))

month_name <- data.frame(mid_month = seq(1, 12, 1),
	month_char = c("DECJAN", "JANFEB", "FEBMAR", "MARAPR", "APRMAY", "MAYJUN", "JUNJULY", "JULYAUG", "AUGSEP", "SEPOCT", "OCTNOV", "NOVDEC"))

enso_long <- plot_df %>%
	pivot_longer(c(-YEAR,-data),  values_to = "values", names_to = "month") %>%
	rename(month_char = month) %>%
	rename(enso = values) %>%
	left_join(month_name, by = "month_char")

enso_long <- enso_long %>%
		mutate(date = as.Date(paste0(YEAR, "-", mid_month, "-01")))

		p <- ggplot(enso_long, aes(x=date, y= enso, colour = data, group = data)) %>%
			+ geom_hline(yintercept = 0, colour = "grey50", alpha = 0.5) %>%
			+ geom_line() %>%
			+ scale_colour_brewer(name = "Source", palette = "Set2", type = "qual") %>%
			+ scale_x_date(name = "Date", breaks = seq(as.Date("1800-01-01"), as.Date("2040-01-01"), by = "10 years")) %>%
			+ scale_y_continuous(name = "MEI") %>%
			+ theme_classic(9) %>%
			+ theme(legend.position="bottom")


### Save plot
ggsave(file.path(write_figures_path, "mei_ts_all_data.png"), p,  width = 5.5, height = 2, dpi = 600)
ggsave(file.path(write_figures_path, "mei_ts_all_data.pdf"), p,  width = 5.5, height = 2)
ggsave(file.path(write_figures_path, "mei_ts_all_data.svg"), p,  width = 5.5, height = 2)

#head(enso_long)

###########################################################################
###  Clean up the data and join
###########################################################################
### Convert water_year to numerical
djf_all <- djf_all %>%
	rename(year = water_year) %>%
	mutate(year = as.numeric(as.character(year)))
### Because ENSO DJ refers to the January year, we are matched up


p <- ggplot(djf_all, aes(x=year, y= djf_mm_daily*30, group = id)) + geom_line() + scale_y_continuous(name = "Precip (mm per month)") + facet_grid(source ~ .)
ggsave(file.path(write_figures_path, "all_time_series_vertical.png"), p,  width = 10, height = 15, dpi = 600)


p <- ggplot(djf_all, aes(x=year, y= djf_mm_daily*30, group = id)) + geom_line() + scale_y_continuous(name = "Precip (mm per month)") + facet_grid(. ~ source )
ggsave(file.path(write_figures_path, "all_time_series_horizontal.png"), p,  width = 25, height = 7, dpi = 600)

### Filter them out in addition to requiring at most one missing week
djf_all <- djf_all %>%
	filter(n >=83) %>%
	filter(djf_mm_daily < 500 & djf_mm_daily >= 0)

### 4851 gauges
length(unique(djf_all$id))


###########################################################################
###  Clean up the stations and join
###########################################################################
### Filter to only those gauges
stations <- stations %>%
	filter(element == "PRCP") %>%
	filter(id %in% unique(djf_all$id))

head(stations)

###
djf_df <- djf_all %>%
	left_join(stations %>% select(id, latitude, longitude, elevation, state, name), by = "id") %>%
	mutate(id = as.factor(id))

head(djf_df)

djf_df <- djf_df %>%
	drop_na() %>%
	left_join(enso_dj, by = "year")

head(djf_df)
dim(djf_df)

### Remove negative elevations
djf_df <- djf_df %>%
	filter(elevation != -999.9)

###########################################################################
###  Make sure at least 20 years of data
###########################################################################
### Calculate median DJF precip
djf_summary <- djf_df %>%
	group_by(id) %>%
	summarize(djf_mm_daily_median = median(djf_mm_daily, na.rm=TRUE), djf_mm_daily_mean = mean(djf_mm_daily, na.rm=TRUE), n = n())

### Make sure there is at least 20 years of precip
djf_summary <- djf_summary %>%
	filter(n >= 20)

### Cut to only those with more than 20 years
djf_df <- djf_summary %>%
	select(id) %>%
	left_join(djf_df, by = "id")

### Number of gauges
dim(djf_summary)
length(unique(djf_df$id))

### 3864 gauges

###########################################################################
###  Save stations and djf_df
###########################################################################
### For mapping
#dyn.load("/usr/local/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
#dyn.load("/usr/local/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
require(sf)
require(rnaturalearth)
require(rnaturalearthhires)
#require(stars)
#library(maps)

### Convert to sf object
### Create an sf object of gauge locations
gauges_sf = st_as_sf(stations, coords = c("longitude", "latitude"), crs = 4326 )



### Reproject into USA CONUS Albert conic
gauges_sf_proj <- gauges_sf %>%
	st_transform(5070)

proj_coord <- st_coordinates(gauges_sf_proj)
proj_coord <- data.frame(id = gauges_sf_proj$id, x_proj = proj_coord[,1], y_proj = proj_coord[,2])

djf_df <- djf_df %>%
	left_join(proj_coord, by = "id")

### Save objects
saveRDS(djf_df, file = file.path(write_output_path, "djf_df.rds"))
saveRDS(gauges_sf, file = file.path(write_output_path, "gauges_sf.rds"))
saveRDS(gauges_sf_proj, file = file.path(write_output_path, "gauges_sf_proj.rds"))


###########################################################################
###  Download background data
###########################################################################
### Download data for plot
usa_states <- ne_states(country = "United States of America", returnclass='sf') %>%
	filter(name != "Hawaii" & name != "Alaska")
ca_states <- ne_states(country = "Canada", returnclass='sf')
mx_states <- ne_states(country = "Mexico", returnclass='sf')

all_states <- usa_states %>%
	select(name) %>%
	bind_rows(ca_states %>% select(name)) %>%
	bind_rows(mx_states %>% select(name))

lakes <- ne_download(scale = 50, type = 'lakes', category = 'physical', returnclass = 'sf')
rivers <- ne_download(scale = 10, type = 'rivers_lake_centerlines', category = 'physical', returnclass = 'sf')
#land <- ne_download(scale = 50, type = 'land', category = 'physical', returnclass = 'sf')
#background_orig <- ne_download(scale = 50, type = 'NE1_50M_SR_W', category = 'raster')
#background_orig <- ne_download(scale = 10, type = 'NE1_LR_LC_SR_W', category = 'raster')
### This was the one - no longer working
#background_orig <- ne_download(scale = 50, type = 'HYP_50M_SR_W', category = 'raster')
#background_orig <- ne_download(scale = 10, type = 'HYP_LR_SR_W_DR', category = 'raster')

### Convert background for easy plotting
#bb = st_bbox(c(xmin = -115, ymin = 34, xmax = -100, ymax = 50), crs = st_crs(background_orig))
#background <- st_as_stars(background_orig, proxy= FALSE) %>% st_crop(bb) %>%st_rgb(3)

###########################################################################
###  Create a continent boundary
###########################################################################

n_amer_bound <- all_states %>%
	st_transform(5070) %>%
	mutate(continent = "n_amer") %>%
	st_buffer(10000) %>%
    group_by(continent) %>%
    summarise(geometry = sf::st_union(geometry)) %>%
    ungroup() %>%
	st_simplify()

#plot(na_bound[1])

n_amer_bound_lat_long <- n_amer_bound %>%
	st_transform(4326)

#dev.new()
#plot(na_bound_lat_long[1])


### Save objects
saveRDS(all_states, file = file.path(write_output_path, "all_states.rds"))
saveRDS(n_amer_bound, file = file.path(write_output_path, "n_amer_bound.rds"))
saveRDS(n_amer_bound_lat_long, file = file.path(write_output_path, "n_amer_bound_lat_long.rds"))


###########################################################################
###  Create Plots
###########################################################################
### Create plot
p <- ggplot(data = gauges_sf) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = rivers, colour = "#4a80f5", alpha = 0.15) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4", alpha = 0.15) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = mx_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(alpha = 0.7, size = 0.3) %>%
	+ scale_fill_identity() %>%
	+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
	+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ theme_bw(8) %>%
	+ theme(legend.position="bottom") %>%
	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p

### Save plot
ggsave(file.path(write_figures_path, "gauge_overview.png"), p,  width = 6, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "gauge_overview.pdf"), p,  width = 6, height = 6.5)

### Reproject
p <- p + coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)
ggsave(file.path(write_figures_path, "gauge_overview_projected.png"), p,  width = 6, height = 6.5, dpi = 600)


p <- ggplot(data = gauges_sf) %>%
#+ geom_stars( data = foo, alpha = 0.95, aes(fill = elev)) %>%
#+ geom_tile(alpha = 0.95, aes(x = lon, y=lat, fill = elevation)) %>%
+ geom_sf(aes(colour = elevation), size = 0.7, alpha = 0.9) %>%
+ geom_sf(data = canada_states, fill = NA, alpha = 0.5) %>%
+ geom_sf(data = mexico_states, fill = NA, alpha = 0.5) %>%
+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
#+ geom_sf(data = lakes, fill = "#9bbff4", alpha = 0.2) %>%
#+ geom_sf(data = rivers, colour = "#fb8072", alpha = 0.6) %>%	#"#4a80f5"
+ scale_colour_viridis(name = "Elevation (m)", limits = c(-75, 3750)) %>%
+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
+ xlab("Longitude") %>%
+ ylab("Latitude") %>%
+ theme_bw(8) %>%
+ theme(legend.position="bottom") %>%
+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p

### Save plot
ggsave(file.path(write_figures_path, "gauge_elevation.png"), p,  width = 6, height = 6.5, dpi = 600)
#ggsave(file.path(write_figures_path, "gauge_elevation.pdf"), p,  width = 6, height = 6.5)


### Density of Elevations
elev_df <- djf_df %>%
	group_by(id) %>%
	summarize(elevation = median(elevation))

p <- ggplot(elev_df, aes(x = elevation)) + geom_density() + theme_classic()

### Save plot
ggsave(file.path(write_figures_path, "elevation_dist.png"), p,  width = 6, height = 4.5, dpi = 600)
ggsave(file.path(write_figures_path, "elevation_dist.pdf"), p,  width = 6, height = 4.5)

### Summarize elevation
mean(elev_df$elevation, na.rm=TRUE) ##1327
median(elev_df$elevation, na.rm=TRUE) ##1260
min(elev_df$elevation, na.rm=TRUE) #3
max(elev_df$elevation, na.rm=TRUE)  ##3536


### Create plot of precipitation
plot_df <- gauges_sf %>%
	inner_join(djf_summary, by = "id")

p <- ggplot(data = plot_df) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(colour = djf_mm_daily_mean * 30)) %>%
	+ scale_fill_identity() %>%
	+ coord_sf(xlim = c(-116, -99), ylim = c(28, 52), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ scale_colour_viridis(name = "Mean DJF\nPrecip\n(mm/month)")

p

### Save plot
ggsave(file.path(write_figures_path, "gauge_mean_precip.png"), p,  width = 6, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "gauge_mean_precip.pdf"), p,  width = 6, height = 6.5)


### Reproject
p <- p + coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)
ggsave(file.path(write_figures_path, "gauge_mean_precip_projected.png"), p,  width = 6, height = 6.5, dpi = 600)


p <- ggplot(data = plot_df) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(colour = djf_mm_daily_median * 30)) %>%
	+ scale_fill_identity() %>%
	+ coord_sf(xlim = c(-116, -99), ylim = c(28, 52), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ scale_colour_viridis(name = "Median DJF\nPrecip\n(mm/month)")

p

### Save plot
ggsave(file.path(write_figures_path, "gauge_median_precip.png"), p,  width = 6, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "gauge_median_precip.pdf"), p,  width = 6, height = 6.5)

### Reproject
p <- p + coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)
ggsave(file.path(write_figures_path, "gauge_median_precip_projected.png"), p,  width = 6, height = 6.5, dpi = 600)



p <- ggplot(data = plot_df) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(colour = log10(djf_mm_daily_median * 30))) %>%
	+ scale_fill_identity() %>%
	+ coord_sf(xlim = c(-116, -99), ylim = c(28, 52), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude") %>%
	+ scale_colour_viridis(name = "Median DJF\nPrecip\n(mm/month)")

p

### Save plot
ggsave(file.path(write_figures_path, "gauge_median_log10_precip.png"), p,  width = 6, height = 6.5, dpi = 600)
ggsave(file.path(write_figures_path, "gauge_median_log10_precip.pdf"), p,  width = 6, height = 6.5)


### Reproject
p <- p + coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)
ggsave(file.path(write_figures_path, "gauge_median_log10_precip_projected.png"), p,  width = 6, height = 6.5, dpi = 600)



	### Plot distribution
	p <- ggplot(plot_df, aes(x=elevation)) %>%
		+ geom_histogram( 	binwidth = 150, colour = "black", fill = "grey80") %>%
		+ scale_x_continuous(name = "Elevation (m)") %>%
		+ scale_y_continuous(name = "Stations") %>%
		+ theme_classic(9)
	p


		### Save plot
		ggsave(file.path(write_figures_path, "elev_dist.png"), p,  width = 4.5, height = 4.5, dpi = 600)
		ggsave(file.path(write_figures_path, "elev_dist.pdf"), p,  width = 4.5, height = 4.5)
		ggsave(file.path(write_figures_path, "elev_dist.svg"), p,  width = 4.5, height = 4.5)


		p <- p +  geom_density(aes(y=150 * ..count..), colour = "blue")

			### Save plot
			ggsave(file.path(write_figures_path, "elev_dist_line.png"), p, width = 4.5, height = 4.5, dpi = 600)
			ggsave(file.path(write_figures_path, "elev_dist_line.pdf"), p,   width = 4.5, height = 4.5)
			ggsave(file.path(write_figures_path, "elev_dist_line.svg"), p,   width = 4.5, height = 4.5)

			### Save plot
			ggsave(file.path(write_figures_path, "elev_dist_line_wide.png"), p, width = 7.5, height = 1.5, dpi = 600)
			ggsave(file.path(write_figures_path, "elev_dist_line_wide.pdf"), p,   width = 7.5, height = 1.5)
			ggsave(file.path(write_figures_path, "elev_dist_line_wide.svg"), p,   width = 7.5, height = 1.5)



###########################################################################
###  Create object for plotting
###########################################################################
require(raster)
require(gstat)
require(rgdal)
require(sf)
require(stars)

### Load in north american bounds
n_amer_bound <- readRDS(file.path(output_path, "gauge_data/n_amer_bound.rds"))

### Create a bounding box
lon_lim <- c(-116, -99)
lat_lim <- c(25, 50)


### Download elevations for continental US
elev_us <- raster::getData('alt', country='USA', mask=TRUE)
elev_mx <- raster::getData('alt', country='MEX', mask=TRUE)
elev_can <- raster::getData('alt', country='CAN', mask=TRUE)
elev_comb <- raster::merge(elev_us[[1]], elev_us[[2]], elev_mx, elev_can)

### Convert to stars
elev_stars <- st_as_stars(elev_comb)

### Make coarse
elev_coarse = st_as_stars(st_bbox(elev_stars), dx = 0.1)
elev_coarse = st_warp(elev_stars, elev_coarse, method = "bilinear", use_gdal = TRUE, no_data_value = NA_real_)

elev_coarse <- elev_coarse %>%	st_transform(5070)
elev_coarse <- elev_coarse[n_amer_bound]

plot_grid_sf <- st_as_sf(elev_coarse )
plot_grid_sf <- plot_grid_sf %>%
	rename(elevation = 1) %>%
	mutate(point = seq(1, dim(plot_grid_sf)[1]))


plot_grid_centroid <- plot_grid_sf %>%
	 st_centroid()

plot_grid_centroid <- plot_grid_centroid %>%
	mutate(x_proj = unlist(map(plot_grid_centroid$geometry,1)),
		 y_proj = unlist(map(plot_grid_centroid$geometry,2)))

plot_grid_centroid <- plot_grid_centroid	 %>%
	 st_transform(4326)

plot_grid_centroid <- as.data.frame(plot_grid_centroid, xy=TRUE) %>%
	mutate(lon = unlist(map(plot_grid_centroid$geometry,1)),
	 	lat = unlist(map(plot_grid_centroid$geometry,2))) %>%
		select(point, lon, lat, x_proj, y_proj)

plot_grid_sf <- plot_grid_sf %>%
	left_join(plot_grid_centroid, by = "point")

plot_grid_sf <- plot_grid_sf %>%
	filter(lat >= lat_lim[[1]] & lat <= lat_lim[[2]] & lon >= lon_lim[[1]] & lon <= lon_lim[[2]])


	### Create plot
	p <- ggplot(plot_grid_sf) %>%
		#+ geom_stars( data = foo, alpha = 0.95, aes(fill = elev)) %>%
		#+ geom_tile(alpha = 0.95, aes(x = lon, y=lat, fill = elevation)) %>%
		+ geom_sf(alpha = 0.95, aes(fill = elevation), colour = NA) %>%
		+ geom_sf(data = canada_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = mexico_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4", alpha = 0.2) %>%
		+ geom_sf(data = rivers, colour = "#fb8072", alpha = 0.6) %>%	#"#4a80f5"
		+ scale_fill_viridis(name = "Elevation (m)", limits = c(-75, 3750)) %>%
		+ coord_sf(xlim = c(-19.7e5, -1.7e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
		+ scale_x_continuous(breaks = seq(-120, -90, by = 5)) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude") %>%
		+ theme_bw(8) %>%
		+ theme(legend.position="bottom") %>%
		+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

	p

#ggsave(file.path(write_figures_path, "interp_elev.pdf"), p,  width = 6, height = 6.5)
ggsave(file.path(write_figures_path, "elev_srtm.png"), p,  width = 6, height = 6.5, dpi = 600)

### Save object
#saveRDS(plot_grid, file = file.path(write_output_path, "plot_grid.rds"))
saveRDS(plot_grid_sf, file = file.path(write_output_path, "plot_grid_sf.rds"))
