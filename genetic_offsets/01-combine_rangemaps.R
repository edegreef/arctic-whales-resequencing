# Combining range maps for the 3 species

library(tidyverse)
library(raster)
library(rworldmap)
library(mapproj)
library(ggspatial)

# Load range shape files
narwhal_range <- shapefile('C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/maps/range_maps/redlist_species_data_narwhal/data_0.shp')
bowhead_range <- shapefile('C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/maps/range_maps/redlist_species_data_bowhead/data_0.shp')
beluga_range <- shapefile('C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/maps/range_maps/redlist_species_data_beluga/data_0.shp')


# Prepare map
map_outline <- getMap(resolution="high")
crop.extent <- extent(-110,-45, 44, 80)

# Crop outline to boundary and convert to dataframe
map_outline <- crop(map_outline, y=crop.extent) %>% fortify()

# Load in data for sites points

# narwhal and bowhead:
narbow_sites <- read.csv("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current-June2024/climate_data/site_locations.csv", header=T)
nar <- subset(narbow_sites, Narwhal>0)
bow <- subset(narbow_sites, Bowhead>0)

# beluga
bel <- readxl::read_xlsx("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current-June2024/climate_data/Sites_beluga.xlsx")

# Plot
ggplot()+
  geom_polygon(data=beluga_range, aes(x=long, y=lat, group=group), fill="#440d55c4", alpha=0.5)+
  geom_polygon(data=bowhead_range, aes(x=long, y=lat, group=group), fill="#0489eaff", alpha=0.5)+
  geom_polygon(data=narwhal_range, aes(x=long, y=lat, group=group), fill="#fce824c4", alpha=0.5)+
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  geom_point(data=nar, aes(x=est_longitude, y=est_latitude), pch=16, size=2, alpha=1,colour="#f2ff00ff")+
 geom_point(data=bow, aes(x=est_longitude, y=est_latitude), pch=16, size=2, alpha=1,colour="#0489eaff")+
  geom_point(data=bel, aes(x=Lon, y=Lat), pch=16, size=2, alpha=1,colour="#440d55ff")+
  xlab("Longitude")+
  ylab("Latitude")+
  coord_fixed(ratio=2.1,xlim=c(-110,-45),ylim=c(44,80))+
  #coord_fixed(ratio=2.3,xlim=c(-110,-45),ylim=c(44,80))+
  #coord_quickmap(expand=F)+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))+
  annotation_scale(height=unit(0.15, "cm"), location="tr", aes(width_hint=0.15), text_cex=0.8)+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
                         location = "tr", which_north = "true",
                         pad_x = unit(0.12, "cm"), pad_y = unit(0.6, "cm"),
                         style = ggspatial::north_arrow_fancy_orienteering())+
  theme(text=element_text(size=15))+
  theme(legend.title = element_text( size=12), legend.text=element_text(size=12))

ggsave("species_range_overlap_shapefiles_crop2_TEMP.png", width=6, height=7, dpi=800)  
# Save as SVG to touch up in inkscape