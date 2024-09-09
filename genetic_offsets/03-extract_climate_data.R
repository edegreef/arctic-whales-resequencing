# Script to extract climate data

# For quick plot, can follow: https://www.bio-oracle.org/code.php
# Used a lot of help from Tom Jenkin's tutorial for plotting and extracting points for multiple layers: https://tomjenkins.netlify.app/tutorials/r-extract-marine-data/ 

### Step 1) Test run with sdmpredictors, and also set up crop boundary
### Step 2) Download climate data rasters/layers 
### Step 3) Make nicer plots and also one for difference between present/future
### Step 4) Extract data from rasters for each sampling location 

setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current-June2024/climate_data")

# Load packages
library(sdmpredictors)
library(raster)
library(tidyverse)
library(sp)

### Step 1) Start in sdmpredictors and do a test run
# Can take a quick look at layers. set marine=TRUE to look at only marine data
list_layers(marine=TRUE)

# It's a very long list, but can save as a csv
list <- list_layers(marine=TRUE)
#write.csv(list, "list_layers_marine.csv") 

# Use variables with "BO22" for the most updated version

# example for some variables:
# BO22_icethickltmin_ss = ice thickness mean
# BO22_icecovermean_ss = ice cover mean
# BO22_tempmean_ss = sea surface temp mean

# Quick test on one of the climate variables.
# load raster for sea surface temp.
temp <- load_layers("BO22_tempmean_ss")

# Plot for fun
# Crop raster to fit region
crop.extent <- extent(-110,-45, 44, 80)
temp.crop <- crop(temp, crop.extent)

# Generate a nice color ramp and plot the map
my.colors <- colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))
plot(temp.crop,col=my.colors(1000),axes=FALSE, box=FALSE)
title("sst mean - present")


## Skip to line 70 if choosing specific variables
# Can create a list of variables based on key words
variables <- c("temp", "ice", "salinity", "dissox", "ph", "phosphate", "chlo", "curvel", "iron", "nitrate", "pp")

# Extract present-day data sets
present <- list_layers(list) %>%
  dplyr::select(dataset_code, layer_code, name, units, description, contains("cellsize"), version) %>%
  dplyr::filter(grepl(paste(variables, collapse = "|"), layer_code))

# Export present-day data sets to csv
#write_csv(present, path="bio-oracle-present-datasets.csv")

# Future Representative Concentration Pathway (RCP) scenarios of interest
# RCP85 == RCP 8.5, the "worst case" scenario (re greenhouse gas emissions)
rcp <- c("RCP26","RCP45","RCP60","RCP85")

# Extract future data sets (available for 2050 and 2100)
future <- list_layers_future(list) %>%
  dplyr::filter(grepl(paste(rcp, collapse = "|"), scenario)) %>% 
  dplyr::filter(year == 2050 | year == 2100) %>% 
  dplyr::filter(grepl(paste(variables, collapse = "|"), layer_code))

# Export future data sets to csv file
#write_csv(future, path="bio-oracle-future-datasets.csv")
## 

### Step 2) Download climate data rasters/layers 

# Select variables
# For narwhal and bowhead whale, choosing 4 variables (sst, ice, salinity, chloro): 
present.vars <- c("BO22_tempmean_ss", #sea surface temp
                  "BO22_icethickmean_ss", #ice thickness
                  "BO22_salinitymean_ss", #salinity
                  "BO22_chlomean_ss" #chlorophyll concentration mean
)

# For these future climate variables, don't necessarily have to use these values yet- will be uploading full rasters into the gradient forest analysis. But extracting some here to take a look at the spread within our sampling sites.

# There is probably a faster way to slap on the RCP60_2050 and 2100
# Also it looks like there are not future data for all variables - like no iron, phosphate, etc. so need to check list before poppin them in.

# Also for some reason the future projections for chlorophyll have different extent (ymin -79.5 instead of -90. not relevant for this project but annoying for data upload)

future.vars.60 <- c("BO22_RCP60_2050_tempmean_ss","BO22_RCP60_2100_tempmean_ss",
                   "BO22_RCP60_2050_icethickmean_ss","BO22_RCP60_2100_icethickmean_ss",
                   "BO22_RCP60_2050_salinitymean_ss","BO22_RCP60_2100_salinitymean_ss")
future.vars.85 <- c("BO22_RCP85_2050_tempmean_ss","BO22_RCP85_2100_tempmean_ss",
                 "BO22_RCP85_2050_icethickmean_ss","BO22_RCP85_2100_icethickmean_ss",
                 "BO22_RCP85_2050_salinitymean_ss","BO22_RCP85_2100_salinitymean_ss")
future.chlo <- c("BO22_RCP60_2050_chlomean_ss", "BO22_RCP60_2100_chlomean_ss", "BO22_RCP85_2050_chlomean_ss", "BO22_RCP85_2100_chlomean_ss")

# Combine present-day and future vectors
data <- c(present.vars)#, future.vars.60, future.vars.85)#, future.chlo)
data.chlo <- c(future.chlo)

# Note that sometimes url::downloads might timeout. a way around this:
#getOption('timeout') # right now my limit is 60 seconds
#options(timeout=300) # increase this to 5 minutes?

# Download rasters!
#data.rasters <- load_layers(data)

# Tip: to avoid re-loading these layers in the future, can save in directory, example:
#data.rasters <- load_layers(data, datadir="C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current-June2024/climate_rasters/")

# To load in the downloaded rasters:
options(sdmpredictors_datadir="C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current-June2024/climate_rasters")

data.rasters <- load_layers(data)
data.rasters.chlo <- load_layers(data.chlo)

# Crop to specified extent
data.rasters.crop <- crop(data.rasters, crop.extent)
data.rasters.crop.chlo <- crop(data.rasters.chlo, crop.extent)

# Merge chlorophyll ones with the rest
data.rasters.all <- stack(data.rasters.crop, data.rasters.crop.chlo)

# Quick plot
plot(data.rasters.all, col=my.colors(100), fill="gray80")

# Plot sea surface mean temp
raster::subset(data.rasters.all, grep("tempmean_ss", names(data.rasters.all), value = TRUE)) %>%
  plot(col = my.colors(100))#, zlim = c(9,17), axes = FALSE, box = FALSE)

# plot ice
raster::subset(data.rasters.all, grep("icethickmean_ss", names(data.rasters.all), value = TRUE)) %>%
  plot(col = my.colors(100))#, zlim = c(9,17), axes = FALSE, box = FALSE)

# plot salinity
raster::subset(data.rasters.all, grep("salinitymean_ss", names(data.rasters.all), value = TRUE)) %>%
  plot(col = my.colors(100))#, zlim = c(9,17), axes = FALSE, box = FALSE)

# plot chloro
raster::subset(data.rasters.all, grep("chlomean_ss", names(data.rasters.all), value = TRUE)) %>%
  plot(col = my.colors(100))#, zlim = c(9,17), axes = FALSE, box = FALSE)


### Step 3) Make nicer plots and also one for difference between present/future
temp <- c("BO22_tempmean_ss")
ice <- c("BO22_icethickmean_ss")
sal <- c("BO22_salinitymean_ss")
chlo <- c("BO22_chlomean_ss")

data.temp <- load_layers(temp)
data.ice <- load_layers(ice)
data.sal <- load_layers(sal)
data.chlo <- load_layers(chlo)

data.temp <- crop(data.temp, crop.extent)
data.ice <- crop(data.ice, crop.extent)
data.sal <- crop(data.sal, crop.extent)
data.chlo <- crop(data.chlo, crop.extent)

# Add in future data
fut.temp <- c("BO22_RCP85_2100_tempmean_ss")
fut.ice <- c("BO22_RCP85_2100_icethickmean_ss")
fut.sal <- c("BO22_RCP85_2100_salinitymean_ss")
fut.chlo <- c("BO22_RCP85_2100_chlomean_ss")

data.fut.temp <- load_layers(fut.temp)
data.fut.ice <- load_layers(fut.ice)
data.fut.sal <- load_layers(fut.sal)
data.fut.chlo <- load_layers(fut.chlo)

data.fut.temp <- crop(data.fut.temp, crop.extent)
data.fut.ice <- crop(data.fut.ice, crop.extent)
data.fut.sal <- crop(data.fut.sal, crop.extent)
data.fut.chlo <- crop(data.fut.chlo, crop.extent)

# Subtract future and present
temp.diff <- (data.fut.temp$`BO22_RCP85_2100_tempmean_ss`-data.temp$`BO22_tempmean_ss`)
ice.diff <- (data.fut.ice$`BO22_RCP85_2100_icethickmean_ss`-data.ice$`BO22_icethickmean_ss`)
sal.diff <- (data.fut.sal$`BO22_RCP85_2100_salinitymean_ss`-data.sal$`BO22_salinitymean_ss`)
chlo.diff <- (data.fut.chlo$`BO22_RCP85_2100_chlomean_ss`-data.chlo$`BO22_chlomean_ss`)

# Quick plot
my.colors <- colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))
plot(temp.diff, col=my.colors(1000))
plot(ice.diff, col=my.colors(1000))
plot(sal.diff, col=my.colors(1000))
plot(chlo.diff, col=my.colors(1000))

# Make nicer plot
library(rnaturalearth)
library(rworldmap)
library(raster)
library(terra)
library(tidyterra)
library(ggspatial)
library(svglite)
library(patchwork)
library(colorspace)

map_outline <- getMap(resolution="high")
crop.extent <- extent(-110,-45, 44, 80)
# Crop outline to boundary and convert to dataframe
map_outline <- crop(map_outline, y=crop.extent) %>% fortify()

data.temp.spat <- rast(data.temp)
data.ice.spat <- rast(data.ice)
data.sal.spat <- rast(data.sal)
data.chlo.spat <- rast(data.chlo)


temp.diff.spat <- rast(temp.diff)
ice.diff.spat <- rast(ice.diff)
sal.diff.spat <- rast(sal.diff)
chlo.diff.spat <- rast(chlo.diff)


# This map shows sum of the 3 species
#svglite("MAP_SST.svg", width=6, height=6)
SST <- ggplot()+
  # To plot present data only
 #geom_spatraster(data=data.temp.spat)+ #present data
 #scale_fill_gradient2(low="#5E85B8",mid="#EDF0C0",high="#C13127", midpoint=7.5)+ 
 #ggtitle("Sea surface temperature")+

  # To plot difference with future data
  geom_spatraster(data=temp.diff.spat)+ 
  scale_fill_gradient2(low="#5E85B8",mid="#EDF0C0",high="#C13127", midpoint=3.5)+ 
  ggtitle("Change in sea surface temperature")+ 

  # Below same for both
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  scale_x_continuous(expand=c(0,0), breaks=c(-100,-80,-60)) +
  scale_y_continuous(expand=c(0,0), breaks=c(50,60,70,80)) +
  xlab("Longitude")+
  ylab("Latitude")+
  labs(fill="°C")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.title = element_text(hjust = 0.5))+
  annotation_scale(height=unit(0.15, "cm"), location="bl", aes(width_hint=0.3), text_cex=0.8)+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
    location = "bl", which_north = "true",
    pad_x = unit(0.1, "cm"), pad_y = unit(0.53, "cm"),
    style = ggspatial::north_arrow_fancy_orienteering())+
  theme(text=element_text(size=15)) 
#dev.off()

ICE <- ggplot()+
  # To plot present data only
 #geom_spatraster(data=data.ice.spat)+ 
 #scale_fill_gradient2(low="#5E85B8",mid="#EDF0C0",high="#C13127", midpoint=3.5)+ 
 #ggtitle("Sea ice thickness")+
  
  # To plot difference with future data
  geom_spatraster(data=ice.diff.spat)+
  scale_fill_gradient2(low="#C13127",mid="#EDF0C0",high="#5E85B8", midpoint=-3)+
  ggtitle("Change in sea ice thickness")+
  
  # Below same for both
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  scale_x_continuous(expand=c(0,0), breaks=c(-100,-80,-60)) +
  scale_y_continuous(expand=c(0,0), breaks=c(50,60,70,80)) +
  xlab("Longitude")+
  ylab("Latitude")+
  labs(fill="m")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.title = element_text(hjust = 0.5))+
  annotation_scale(height=unit(0.15, "cm"), location="bl", aes(width_hint=0.3), text_cex=0.8)+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
                         location = "bl", which_north = "true",
                         pad_x = unit(0.1, "cm"), pad_y = unit(0.53, "cm"),
                         style = ggspatial::north_arrow_fancy_orienteering())+
  theme(text=element_text(size=15))

ICE

SAL <- ggplot()+
  # To plot present data only
  #geom_spatraster(data=data.sal.spat)+
  #scale_fill_gradient2(low="#5E85B8",mid="#EDF0C0",high="#C13127", midpoint=21)+
  #ggtitle("Sea surface salinity")+
  
  # To plot difference with future data
  geom_spatraster(data=sal.diff.spat)+
  scale_fill_gradient2(low="#C13127",mid="#EDF0C0",high="#5E85B8", midpoint=-1.25)+
  ggtitle("Changes in sea surface salinity")+  
  
  # Below same for both
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  scale_x_continuous(expand=c(0,0), breaks=c(-100,-80,-60)) +
  scale_y_continuous(expand=c(0,0), breaks=c(50,60,70,80)) +
  xlab("Longitude")+
  ylab("Latitude")+
  labs(fill="PSS")+

  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.title = element_text(hjust = 0.5))+
  annotation_scale(height=unit(0.15, "cm"), location="bl", aes(width_hint=0.3), text_cex=0.8)+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
                         location = "bl", which_north = "true",
                         pad_x = unit(0.1, "cm"), pad_y = unit(0.53, "cm"),
                         style = ggspatial::north_arrow_fancy_orienteering())+
  theme(text=element_text(size=15))

SAL
CHLO <- ggplot()+
  # To plot present data only
  #geom_spatraster(data=data.chlo.spat)+
  #scale_fill_gradient2(low="#5E85B8",mid="#EDF0C0",high="#C13127", midpoint=1.0)+
  #ggtitle("Chlorophyll concentration")+
  
  # To plot difference with future data
  geom_spatraster(data=chlo.diff.spat)+ #difference with future data
  scale_fill_gradient2(low="#5E85B8",mid="#EDF0C0",high="#C13127", midpoint=-0.1)+
  ggtitle("Change in chlorophyll concentration")+
  
  # Below same for both
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  scale_x_continuous(expand=c(0,0), breaks=c(-100,-80,-60)) +
  scale_y_continuous(expand=c(0,0), breaks=c(50,60,70,80)) +
  xlab("Longitude")+
  ylab("Latitude")+
  labs(fill="mg/m³")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.title = element_text(hjust = 0.5))+
  annotation_scale(height=unit(0.15, "cm"), location="bl", aes(width_hint=0.3), text_cex=0.8)+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
                         location = "bl", which_north = "true",
                         pad_x = unit(0.1, "cm"), pad_y = unit(0.53, "cm"),
                         style = ggspatial::north_arrow_fancy_orienteering())+
  theme(text=element_text(size=15))
CHLO

# combine and deal with axes things
patchwork1 = SST+ICE 
patchwork2 = SAL+CHLO

# probably not the most elegant way, but it works (to reduce the redundancy of axes titles)
svglite("MAPS_CLIMATEVAR_differences_RCP8.5_2100.svg", width=10, height=11)
((patchwork1[[1]] + theme(axis.title.x = element_blank())) +
(patchwork1[[2]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank()))) / 
  ((patchwork2[[1]]) +
     (patchwork2[[2]] + theme(axis.title.y = element_blank())))
dev.off()


### Step 4) Extract data from rasters for each sampling location 
# Load site locations (needs coordinates for each sampling/target site)
sites <- read.csv("site_locations_allsp.csv")

# Set x and y for lat and long
x <- sites$est_longitude
y <- sites$est_latitude

# Set up coordinates as a tibble
pts <- cbind(x, y) %>% as_tibble()

# Set pts to match spatial points with the rasters
pts <- SpatialPoints(pts, proj4string = CRS("+proj=longlat +datum=WGS84"))

# Double check that the project matches (should say "TRUE")
projection(pts)==projection(data.rasters.all)

# Plot points on map with each of the rasters #this one is a bit rough
sp::spplot(data.rasters.all,
           scales = list(draw=TRUE), col.regions = my.colors(100),
           sp.layout = c("sp.points", pts, pch = 21, cex = 1, fill = "yellow", col = "black")
)

# Create dataframe to store extracted data
climate.data <- tibble(ID = 1:nrow(pts@coords), Lon = pts$x, Lat = pts$y)

# Extract data for each point
# Combine rasters if multiple
rasters = raster::stack(data.rasters.all)
nlayers(rasters)

# Extract data for each layer. here adding a buffer (radius in metres) to select mean value for each site location. can play around with the value more later, but for now setting buffer to 100000 (for 100km buffer)
store_data = list()
for (i in 1:nlayers(rasters)){
  store_data[[i]] = raster::extract(rasters[[i]], pts, buffer=100000, fun=mean)
}

# Set up data file with extracted values
# Rename variables in the list and then combine data

# First make variable names a bit better
# Let's check what it is now
names(rasters)

# Rename just remember the order has to match, or just leave as is
names(rasters) <- c("sst", "icethick", "salinity", "chloro", "sst_RCP60_2050", "sst_RCP60_2100", "icethick_RCP60_2050","icethick_RCP60_2100", "salinity_RCP60_2050", "salinity_RCP60_2100", "sst_RCP85_2050", "sst_RCP85_2100", "icethick_RCP85_2050","icethick_RCP85_2100", "salinity_RCP85_2050", "salinity_RCP85_2100","chloro_RCP60_2050", "chloro_RCP60_2100", "chloro_RCP85_2050", "chloro_RCP85_2100")

names(store_data) <- names(rasters)
climate.data <- bind_cols(climate.data, as_tibble(store_data))
climate.data

# Want to put my site names back in
climate.data$ID <- sites$Location_abbrev

# Save climate.data
write_csv(climate.data, "climatedata.temp.ice.sal.chlo.100KBUFFER.allsp.csv")

# Visualize spread of extracted points

# Import if starting from here
climate.data <- read.csv("climatedata.temp.ice.sal.chlo.100KBUFFER.allsp.csv")

# edited to remove year 2050

SST <- climate.data %>% 
  # select colums with Lat, sst, and future sst
  dplyr::select(3,4,9,15) %>% 
  # transform data to long format for plotting
  pivot_longer(names_to = "Variable", values_to = "Values", cols = c("sst", "sst_RCP60_2100", "sst_RCP85_2100")) %>% 
  # plot data
  ggplot()+#aes(x=factor(level=level_order)))+
#  geom_violin(aes(x = Variable, y = Values, fill = Variable), alpha=0.5, show.legend = FALSE)+
  geom_quasirandom(aes(x = Variable, y = Values, color=Lat), pch=19,size=2.5, width=0.4, show.legend = TRUE, alpha =0.8)+
  #scale_y_continuous(expand = c(0,0), limits = c(0,16), breaks = c(seq(0,16,2)))+
  scale_x_discrete(limits = c('sst', 'sst_RCP60_2100', 'sst_RCP85_2100'), labels = c("Present",  "2100 (RCP 6.0)", "2100 (RCP 8.5)"))+
  scale_color_gradientn(colours = colorspace::heat_hcl(7))+
  xlab("")+
  ylab(expression("Celcius ("^o*"C)"))+
  ggtitle("Sea surface temperature")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

SST


ICE <- climate.data %>% 
  dplyr::select(3,5,11,17) %>% 
  pivot_longer(names_to = "Variable", values_to = "Values", cols = c("icethick", "icethick_RCP60_2100", "icethick_RCP85_2100")) %>% 
  ggplot()+#aes(x=factor(level=level_order)))+
 # geom_violin(aes(x = Variable, y = Values, fill = Variable), alpha=0.5, show.legend = FALSE)+
  geom_quasirandom(aes(x = Variable, y = Values, color=Lat), pch=19,size=2.5, width=0.4, show.legend = TRUE, alpha =0.8)+
  scale_x_discrete(limits = c('icethick', 'icethick_RCP60_2100', 'icethick_RCP85_2100'), labels = c("Present",  "2100 (RCP 6.0)", "2100 (RCP 8.5)"))+
  scale_color_gradientn(colours = colorspace::heat_hcl(7))+
  xlab("")+
  ylab("Meters (m)")+
  ggtitle("Sea ice thickness")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ICE

SAL <- climate.data %>% 
  dplyr::select(3,6,13,19) %>% 
  pivot_longer(names_to = "Variable", values_to = "Values", cols = c("salinity", "salinity_RCP60_2100", "salinity_RCP85_2100")) %>% 
  ggplot()+#aes(x=factor(level=level_order)))+
 # geom_violin(aes(x = Variable, y = Values, fill = Variable), alpha=0.5, show.legend = FALSE)+
  geom_quasirandom(aes(x = Variable, y = Values, color=Lat), pch=19,size=2.5, width=0.4, show.legend = TRUE, alpha =0.8)+
  scale_x_discrete(limits = c('salinity', 'salinity_RCP60_2100', 'salinity_RCP85_2100'), labels = c("Present",  "2100 (RCP 6.0)", "2100 (RCP 8.5)"))+
  scale_color_gradientn(colours = colorspace::heat_hcl(7))+
  xlab("")+
  ylab("Practical salinity scale (PSS)")+
  ggtitle("Sea surface salinity")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
SAL

CHLO <- climate.data %>% 
  dplyr::select(3,7,21,23) %>% 
  pivot_longer(names_to = "Variable", values_to = "Values", cols = c("chloro", "chloro_RCP60_2100", "chloro_RCP85_2100")) %>% 
  ggplot()+#aes(x=factor(level=level_order)))+
  #geom_violin(aes(x = Variable, y = Values, fill = Variable), alpha=0.5, show.legend = FALSE)+
 # geom_quasirandom(aes(x = Variable, y = Values), pch=21,size=1, colour="black", fill="black",width=0.2, show.legend = FALSE, alpha =1)+
  geom_quasirandom(aes(x = Variable, y = Values, color=Lat), pch=19,size=2.5, width=0.4, show.legend = TRUE, alpha =0.8)+
  scale_x_discrete(limits = c('chloro', 'chloro_RCP60_2100', 'chloro_RCP85_2100'), labels = c("Present",  "2100 (RCP 6.0)", "2100 (RCP 8.5)"))+
  scale_color_gradientn(colours = colorspace::heat_hcl(7))+
  xlab("")+
  ylab("Concentration (mg/m³)")+
  ggtitle("Chlorophyll concentration")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
CHLO


(SST + ICE)/(SAL+CHLO) + plot_layout(guides="collect")

ggsave("beeswarm_plots_100KMBUFFER_allsp_size3_latcolor.png", width=7, height=6, dpi=500)
