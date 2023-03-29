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
#dyn.load("/usr/local/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
#dyn.load("/usr/local/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
require(sf)
require(rnaturalearth)
require(rnaturalearthhires)
#library(maps)


require(rnoaa)
require(tidyverse)
library(sf)
library(stars)
library(rnaturalearth)
library(rnaturalearthhires)
require(viridis)
require(lfstat)
require(mgcv)

select <- dplyr::select


###########################################################################
###  Set output folder
###########################################################################
write_output_path <- file.path(output_path, "gauge_data")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

write_figures_path <- file.path(write_figures_path, "gauge_data")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)

###########################################################################
###  Download gauge data
###########################################################################

### Get stations, ghcnd-stations and ghcnd-inventory merged
stations <- ghcnd_stations()

### Save for later
saveRDS(stations, file = file.path(write_output_path, "all_stations.rds"))

### Read in station info
#stations <- readRDS(file.path(write_output_path, "all_stations.rds"))

###########################################################################
###  Filter to region
###########################################################################
### Filter to region
region_df <- stations %>%
	filter(latitude >= 25 & latitude <= 50 & longitude >= -116 & longitude <= -99) %>%
	filter(element == "PRCP") %>%
	filter(elevation > 0)

dim(region_df)

region_df <- region_df %>%
	mutate(max_length = last_year - first_year + 1)


### Create an sf object of gauge locations
gauges_sf = st_as_sf(region_df, coords = c("longitude", "latitude"), crs = 4326 )


##########################################################
### Download data
##########################################################
region_df <- region_df %>%
	filter(max_length >= 20)

all_ids <- region_df$id
length(all_ids)

save_points <- which(seq(1, length(all_ids)) %% 20 == 0)

tic()

#for (j in seq(1, 250)){
for (j in seq(1, length(all_ids))){

	cat(j)
	cat("\n")

	id <- all_ids[j]

	### Download data
	gauge_data <- meteo_tidy_ghcnd( id,  var = "prcp", keep_flags = TRUE)

	### Process to get djf_summary
	### Combine to monthly data
	### Use USGS water year (beginning in October, use the year of the later part
	### This means the ENSO values are off by a year if they preceed

	### Drop nas and cut to only DJF
	### Use USGS water year (beginning in October, use the year of the later part
	### This means the ENSO values are off by a year if they preceed
	djf_temp <- gauge_data %>%
		drop_na(prcp) %>%
		mutate(month = month(date), year = year(date), water_year = water_year(date, origin = 10), prcp_mm_daily = as.numeric(as.character(prcp))/10) %>%
		filter(month == 12 | month == 1 | month == 2)  %>%
		mutate(water_year = as.numeric(as.character(water_year)))

	### Drop Quality control flags
	djf_temp <- djf_temp %>%
		filter(qflag_prcp == " ")

	### Calculate DJF mean and the count of non na values, drop when too few values
	### 89 Days possible, Allow for a full week of missing values
	djf_temp <- djf_temp %>%
		group_by(id, water_year) %>%
		summarize(djf_mm_daily = mean(prcp_mm_daily, na.rm=TRUE), source = names(which.max(table(sflag_prcp))), n = n(), .groups="drop") %>%
		filter(n >= 82)

	if (j == 1) {
		djf_summary <- djf_temp
	} else{
		djf_summary <- djf_summary %>%
			bind_rows(djf_temp)
	}

	if (j %in% save_points){
		saveRDS(djf_summary, file = file.path(write_output_path, "djf_all.rds"))
	}

}

toc()

saveRDS(djf_summary, file = file.path(write_output_path, "djf_all.rds"))
