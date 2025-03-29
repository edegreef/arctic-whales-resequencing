# Stacking genetic offsets from multiple species

# Load libraries
library(rnaturalearth)
library(rworldmap)
library(ggspatial)
library(tidyterra)
library(terra)
library(svglite) #saving files
library(raster)
library(tidyverse)
library(climateStability) # for the raster scale 0 to 1 function
library(RColorBrewer)

setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current")

#### Step 1) Scale and combine
# Load offset rasters
beluga_offset <- raster("beluga/results_within_range/raster.offset.RCP85.2100.beluga.easternCanada.RANGECLIP.tif")
narwhal_offset <- raster("narwhal/results_within_range/raster.offset.RCP85.2100.narwhal.easternCanada.RANGECLIP.tif")
bowhead_offset <- raster("bowhead/results_within_range/raster.offset.RCP85.2100.bowhead.easternCanada.RANGECLIP.tif")

# Set color scale
offset.colors = colorRampPalette(c("#4575B4", "#abd9e9","#ffffbf","#fdae61","#f46d43","#a50026"))

# Test that rasters loaded in OK
plot(beluga_offset, col=offset.colors(200))
plot(narwhal_offset, col=offset.colors(200))
plot(bowhead_offset, col=offset.colors(200))

# Convert to spatraster
beluga_offset.spat <- rast(beluga_offset)
narwhal_offset.spat <- rast(narwhal_offset)
bowhead_offset.spat <- rast(bowhead_offset)

# scale offset values to 0-1 range --- use this for the 20% chunks
#beluga_offset_scaled <- rescale0to1(beluga_offset.spat)
#narwhal_offset_scaled <- rescale0to1(narwhal_offset.spat)
#bowhead_offset_scaled <- rescale0to1(bowhead_offset.spat)

# Use the scale() function for SDs
beluga_offset_scaled <- scale(beluga_offset.spat)
narwhal_offset_scaled <- scale(narwhal_offset.spat)
bowhead_offset_scaled <- scale(bowhead_offset.spat)

# Look at plots to make sure it worked
plot(beluga_offset_scaled, col=offset.colors(200))
plot(narwhal_offset_scaled, col=offset.colors(200))
plot(bowhead_offset_scaled, col=offset.colors(200))

# Add offset values together with mosaic merge in terra
temp_combined_offset_scaled <- mosaic(beluga_offset_scaled,narwhal_offset_scaled, fun="sum")
combined_offset_scaled <- mosaic(temp_combined_offset_scaled,bowhead_offset_scaled, fun="sum")

# Quick plot
plot(combined_offset_scaled, col=offset.colors(200))

map_outline <- getMap(resolution="high")
crop.extent <- extent(-110,-45, 44, 80)

# Crop outline to boundary and convert to dataframe
map_outline <- crop(map_outline, y=crop.extent) %>% fortify()

# Nice map showing sum of the 3 species
svglite("combined_all3sp_GF_GeneticOffset_RCP85_2100_landmap_stepsn_SCALED_updatedBelK6_scalefunctSDs_RANGECLIP_noscalebar.svg", width=6, height=6)
ggplot()+
  geom_spatraster(data=combined_offset_scaled)+
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  scale_fill_stepsn(n.breaks = 20, colours = offset.colors(12), na.value="white")+#, limits=grad_lims)+
  scale_x_continuous(expand=c(0,0), breaks=c(-100,-80,-60)) +
  scale_y_continuous(expand=c(0,0), breaks=c(50,60,70,80)) +
  xlab("Longitude")+
  ylab("Latitude")+
  labs(fill="Genetic offset")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))+
  # annotation_scale(height=unit(0.15, "cm"), location="tr", aes(width_hint=0.15), text_cex=0.8)+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
                         location = "tr", which_north = "true",
                         pad_x = unit(0.12, "cm"), pad_y = unit(0.3, "cm"), #pad_y unit 0.6 if including scale bar
                         style = ggspatial::north_arrow_fancy_orienteering())+
  theme(text=element_text(size=15))+
  theme(legend.title = element_text( size=12), legend.text=element_text(size=12))
dev.off()

### Next, make map for 20% chunks
# Test -one way with clamp, but later decide to just set contour range
beluga_offset_scaled.clamp.15 <- clamp(beluga_offset_scaled, lower = 0.15, value=FALSE)
plot(beluga_offset_scaled.clamp.15, col="#FFFFB2")

# site points
narbow_sites <- read.csv("climate_data/site_locations.csv", header=T)
nar <- subset(narbow_sites, Narwhal>0)
bow <- subset(narbow_sites, Bowhead>0)
bel <- readxl::read_xlsx("climate_data/Sites_beluga.xlsx")


# plot
bel_clamp <-ggplot()+
  geom_spatraster_contour_filled(data=beluga_offset_scaled, breaks=seq(0.2,1, 0.2))+
  geom_spatraster_contour(data=beluga_offset_scaled, breaks=seq(0.2,1, 0.2))+
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  geom_point(data=bel, aes(x=Lon_adj, y=Lat_adj),size=1.5,colour="black")+
  scale_fill_brewer(palette = "YlOrRd")+
  scale_x_continuous(expand=c(0,0), breaks=c(-100,-80,-60)) +
  scale_y_continuous(expand=c(0,0), breaks=c(50,60,70,80)) +
  xlab("Longitude")+
  ylab("Latitude")+
  labs(fill="Genetic offset")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), panel.background = element_blank())+
  annotation_scale(height=unit(0.15, "cm"), location="tr", aes(width_hint=0.15), text_cex=0.8)+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
                         location = "tr", which_north = "true",
                         pad_x = unit(0.12, "cm"), pad_y = unit(0.6, "cm"),
                         style = ggspatial::north_arrow_fancy_orienteering())


nar_clamp <- ggplot()+
  geom_spatraster_contour_filled(data=narwhal_offset_scaled, breaks=seq(0.2,1, 0.2))+
  geom_spatraster_contour(data=narwhal_offset_scaled, breaks=seq(0.2,1, 0.2))+
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  geom_point(data=nar, aes(x=est_longitude, y=est_latitude),size=1.5,colour="black")+
  scale_fill_brewer(palette = "YlOrRd")+
  scale_x_continuous(expand=c(0,0), breaks=c(-100,-80,-60)) +
  scale_y_continuous(expand=c(0,0), breaks=c(50,60,70,80)) +
  xlab("Longitude")+
  ylab("Latitude")+
  labs(fill="Genetic offset")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), panel.background = element_blank())+
  annotation_scale(height=unit(0.15, "cm"), location="tr", aes(width_hint=0.15), text_cex=0.8)+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
                         location = "tr", which_north = "true",
                         pad_x = unit(0.12, "cm"), pad_y = unit(0.6, "cm"),
                         style = ggspatial::north_arrow_fancy_orienteering())

bow_clamp <- ggplot()+
  geom_spatraster_contour_filled(data=bowhead_offset_scaled, breaks=seq(0.2,1, 0.2))+
  geom_spatraster_contour(data=bowhead_offset_scaled, breaks=seq(0.2,1, 0.2))+
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  geom_point(data=bow, aes(x=est_longitude, y=est_latitude),size=1.5,colour="black")+
  scale_fill_brewer(palette = "YlOrRd")+
  scale_x_continuous(expand=c(0,0), breaks=c(-100,-80,-60)) +
  scale_y_continuous(expand=c(0,0), breaks=c(50,60,70,80)) +
  xlab("Longitude")+
  ylab("Latitude")+
  labs(fill="Genetic offset")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), panel.background = element_blank())+
 # annotation_scale(height=unit(0.15, "cm"), location="tr", aes(width_hint=0.15), text_cex=0.8)+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
                         location = "tr", which_north = "true",
                         pad_x = unit(0.12, "cm"), pad_y = unit(0.3, "cm"),
                         style = ggspatial::north_arrow_fancy_orienteering())

# Plot together
library(patchwork)
bel_clamp + nar_clamp + bow_clamp

# To save separately
svglite("beluga/beluga_K6_GF_GeneticOffset_RCP85_2100_highlight2.svg", width=6, height=6)
bel_clamp
dev.off()

svglite("narwhal/narwhal_GF_GeneticOffset_RCP85_2100_highlight2.svg", width=6, height=6)
nar_clamp
dev.off()

svglite("bowhead/bowhead_GF_GeneticOffset_RCP85_2100_highlight2.svg", width=6, height=6)
bow_clamp
dev.off()

