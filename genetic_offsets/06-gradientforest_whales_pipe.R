# Gradient Forest analyses on whale data

# Lots of help for gradient forest analyses from pgugger's github (https://github.com/pgugger/LandscapeGenomics/blob/master/2019/Exercise4.md)

# For converting vcf to gradient forest file I used beginning of tutorial here (came across this one before the one above, but from same person: https://github.com/pgugger/LandscapeGenomics/blob/master/2018_China/Exercise3.md)

# And help from gradient forest vignette: (https://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf)

# need to make sure rtools is downloaded (https://cran.rstudio.com/bin/windows/Rtools/rtools42/rtools.html)
#install.packages("extendedForest", repos="http://R-Forge.R-project.org")
#install.packages("gradientForest", repos="http://R-Forge.R-project.org")

# load libraries
library(gradientForest)
library(dplyr)
library(vegan)
library(data.table)
library(marmap) # for measuring distances between sites within water

# for climate data stuff
library(raster)
library(tidyverse)
library(sdmpredictors)
library(sp)
library(sf)

# some more to help with editing map plots
library(RStoolbox)
library(rnaturalearth)
library(rworldmap)
library(ggspatial)
library(tidyterra)
library(terra)
library(svglite) #saving files

# more packages for making box plots out of offset raster
library(plyr)
library(reshape2)

# Step 1: Initial run with gradient forest on present data
# Step 2: Create plots for biological space and geographical space
# Step 3: Measuring genetic offset with future data
# Step 4: Create box plots

setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current")

##########  
##### Step 1: Initial run with gradient forest on present data

## 1.1 prepare data files

# Load snp data - note that data file must not have missing data for GF
#snp <- read.table("narwhal_snps_nodup.filteredJul2024_shortfilename_K3_top0.01_scan_4env.forR", header=T, row.names=1)

# use fread() for faster upload
snp_temp <- fread("bowhead/bowhead_RWmap_snps_rmfix_noNEW.PASS.minQ50.miss.biallel.min100kb.autosomes.n19.maf.miss01.imputed.thin1000.forR", header=T)
snp_temp <- as.data.frame(snp_temp)
snp <- snp_temp[,-1]
rownames(snp) <- snp_temp[,1]


# Load environmental data for each site --update this file for each individual
site.clim.dat <- read.csv("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current/climate_data/sample_info_bowhead_with_climdat_100KM.csv")
colnames(site.clim.dat)

# mMake sure ID orders match between sample info and individuals in the snp file!
# Extract sample_ID and location_ID
sample_info <- site.clim.dat[,c("sample_ID", "location_ID")]

# Left_join or something to the climate data
sample_clim <- left_join(sample_info, site.clim.dat, by="sample_ID")

# Pull out present data (I also had future data in site.clim.dat)
pres.clim.points <- sample_clim[c("Lon", "Lat", "sst","icethick", "salinity", "chloro")]

# can rename here if needed
#colnames(pres.clim.points) <- c("Lon", "Lat", "sst","icethick", "salinity", "chloro")

## 1.2 Generate PCNM spatial variables 

# Because sites are within water, using marmap distances instead of linear

# Coordinates
coord <- pres.clim.points[,c('Lon','Lat')]
#pcnm <- pcnm(dist(coord))  # if using linear

### Adding in marmap distances

# Download ocean depth map. Lower resolution number is more fine-scale.
ocean_map <- getNOAA.bathy(lon1 = -105, lon2 = -50, lat1 = 55, lat2 = 85, resolution = 5)
colnames(coord) <- c("x", "y")
trans <- trans.mat(ocean_map, min.depth = 0)
dist <- lc.dist(trans,coord,res="dist")
dist.mat <- as.matrix(dist)
pcnm <- pcnm(dist.mat)

###

# Keep half of positive pcnms
keep <- round(length(which(pcnm$value > 0))/2) 
pcnm.keep <- scores(pcnm)[,1:keep]
pcnm.keep

# Plot the PCNM data

# add coord and pcnm.keep together
pcnm.for.plot <- cbind(coord, pcnm.keep)
map_outline <- getMap(resolution="high")

crop.extent <- extent(-110,-45, 44, 80) #eastern Canadian Arctic for all three sp

# Crop outline to boundary and convert to dataframe
map_outline <- crop(map_outline, y=crop.extent) %>% fortify()

# Plot PCNMs on map- have to go one at a time here (adjust the "PCNM#")
ggplot()+
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  geom_point(data=pcnm.for.plot, aes(x=x, y=y, fill=PCNM4), shape=21, stroke=0.7,color="black",size=abs(20*pcnm.for.plot$PCNM4))+
  scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027")+
  scale_x_continuous(expand=c(0,0), breaks=c(-100,-80,-60)) +
  scale_y_continuous(expand=c(0,0), breaks=c(50,60,70,80)) +
  xlab("Longitude")+
  ylab("Latitude")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.title = element_text(hjust = 0.5))+
  annotation_scale(height=unit(0.15, "cm"), location="tr", aes(width_hint=0.15), text_cex=0.8)+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
                         location = "tr", which_north = "true",
                         pad_x = unit(0.12, "cm"), pad_y = unit(0.6, "cm"),
                         style = ggspatial::north_arrow_fancy_orienteering())+
  theme_classic()+
  theme(text=element_text(size=15),legend.title = element_text( size=12), legend.text=element_text(size=12), panel.border=element_rect(colour="black", fill=NA))
ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current/narwhal/pcnms/narwhal_pcnm4.png", width=4.5, height=4, dpi=400)


# Create file with climate and PCNM variables (no lat/long)
env.gf <- cbind(pres.clim.points[,c("sst", "icethick", "salinity", "chloro")], pcnm.keep)

# Define max number of splits
lev <- floor(log2(0.368*nrow(env.gf)/2))
lev

# Run GF using the climate, spatial, and snp data inputs
gf <- gradientForest(cbind(env.gf, snp), 
                     predictor.vars=colnames(env.gf),
                     response.vars=colnames(snp), 
                     ntree=100, #500 for narwhal, 100 for bowhead 
                     maxLevel=lev, 
                     trace=T, 
                     corr.threshold=0.50,
                     nbin=601) #301 for narwhal, 601 for bowhead
# Experimenting nbin and run time:
# 200 bins for narwhal took 1hr 22 min
# 300 bins for narwhal took 45 min (ntree 500)
# this means, more bins -> less time (to an extent probably).
# For bowhead, ended up using ntree=100, and nbin=600 (otherwise crashes)

# save here too.
#save.image("narwhal_K3_GF_snps_scan_4env_ntree500_nbin301.RData")


## 1.3 Initial gradient forest plots
# tip: use help(par) to see all the plotting options for color, size, etc

# Bar graphs of predictor overall variable importance
# Chose colors to highlight PCNMs and known environmental variables separately (orange for PCNM, blue for variables)
# for narwhal
plot(gf, plot.type = "O", col=c("#80ace1","#80ace1","#80ace1","#80ace1","#f6a565", "#f6a565","#f6a565","#f6a565"), lwd=1.5, cex.axis=0.9)

# for bowhead
plot(gf, plot.type = "O", col=c("#80ace1","#80ace1","#f6a565","#f6a565","#f6a565", "#80ace1","#80ace1","#f6a565","#f6a565"), lwd=1.5, cex.axis=0.9)

# save these as svg

# Organize variables by importance, for other plots
by.importance <- names(importance(gf))
by.importance

# list the values for R2 weighted importance
as.data.frame(importance(gf, type = "Weighted"))

# Splits density plot (showing binned split importance and location on each gradient)
#plot(gf, plot.type="S", imp.vars=by.importance,
#     leg.posn="topright", cex.legend=0.7, cex.axis=0.8,cex.lab = 0.7, line.ylab=0.9, 
#     par.args=list(mgp=c(1.5,0.5,0), mar=c(3.1,1.5,0.1,1)))

# Plot turnover functions, showing how allelic composition changes along spatial or environmental gradients. this is a predictor cumulative plot
plot(gf, plot.type="C", imp.vars=by.importance, 
     show.species=F, common.scale=T, cex.axis=1, cex.lab=1.2, line.ylab=1, lwd=2,col="#7393B3",
     par.args=list(mgp=c(1.5,0.5,0), mar=c(2.5,2,2,2), omi=c(0.2,0.3,0.2,0.4)))
# can't get it to plot 4x2 instead of 2x4 but can adjust in inkscape

# Plot turnover functions for individual loci. each line represent allelic change at a single SNP
# (can take a long time to load. maybe skip for now)
#plot(gf, plot.type="C", imp.vars=by.importance, 
#     show.overall=F, legend=T, leg.posn="topleft", leg.nspecies=5, 
#     cex.lab=1, cex.legend=0.5, cex.axis=1, line.ylab=0.9, 
#     par.args=list(mgp=c(1.5,0.5,0), mar=c(2.5,1,0.1,0.5), omi=c(0,0.3,0,0)))

# Plot R2 measure of the fit of the random forest model 
plot(gf, plot.type="P", show.names=T, horizontal=F, cex.axis=1, cex.labels=0.7, line=2.5)

##########  

##########  
##### Step 2: Create plots for biological space and geographical space

## 2.1 Load climate data rasters
# Note: instead of re-loading the layers from online with sdmpredictors every time I run this during load_layers(), I  downloaded and saved the rasters into a local directory so can load faster from there

# example:
#options(timeout=600) # set timeout limit to 10 min
#data.rasters <- load_layers(present.vars, datadir="C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current/climate_rasters")

# Then set the datadir for sdmpredictors to find the downloaded rasters
options(sdmpredictors_datadir="C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current/climate_rasters")

# Names of present variables of interest
present.vars <- c("BO22_tempmean_ss",
                  "BO22_icethickmean_ss",
                  "BO22_salinitymean_ss", 
                  "BO22_chlomean_ss")

# Can load future vars now too though won't use it until later
future.vars <- c("BO22_RCP85_2100_tempmean_ss",
                 "BO22_RCP85_2100_icethickmean_ss",
                 "BO22_RCP85_2100_salinitymean_ss") 
future.vars.chlo <- c("BO22_RCP85_2100_chlomean_ss") #duno why future chlorophyll data extent is different, so have to load separately

# Options: RCP26, RCP45, RCP60, RCP85: Going to do for RCP45, RCP60, RCP85

# Load rasters
present.rasters <- load_layers(present.vars)
future.rasters.nochlo <- load_layers(future.vars)
future.rasters.chlo <- load_layers(future.vars.chlo)

# Define boundary/extent/crop
crop.extent <- extent(-110,-45, 44, 80) #eastern Canadian Arctic for all three sp

present.rasters.crop <- crop(present.rasters, crop.extent)
future.rasters.nochlo.crop <- crop(future.rasters.nochlo, crop.extent)
future.rasters.chlo.crop <- crop(future.rasters.chlo, crop.extent)

# add the future chlorophyll with other future variables
future.rasters.crop <- stack(future.rasters.nochlo.crop, future.rasters.chlo.crop)

# Rename variables (need to match climate data points)
names(present.rasters.crop) <- c("sst", "icethick", "salinity", "chloro")
names(future.rasters.crop) <- c("sst", "icethick", "salinity", "chloro")

## 2.2 Crop to desired area, using species' range here
# Load range shapefile (obtained from IUCN red list database)
species <- read_sf('C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/maps/range_maps/redlist_species_data_narwhal/data_0.shp')

# Plot for fun, maybe use for a supplemental figure
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))
plot(present.rasters.crop,col=my.colors(200))#,axes=FALSE, box=FALSE) 
plot(future.rasters.crop,col=my.colors(200)) 

# Crop to species range shape file
present.rasters.range <- mask(present.rasters.crop, species)
future.rasters.range <- mask(future.rasters.crop, species)

crop.map2 <- crop(future.rasters.crop,crop.extent)
future.rasters.range <- mask(crop.map2, species)

#plot(present.rasters.range,col=my.colors(1000),axes=FALSE, box=FALSE)
#plot(future.rasters.range,col=my.colors(1000),axes=FALSE, box=FALSE)

# can convert rasters to points with rasterToPoints() if want to keep with vignette, but skipping this b/c want to use other format
#r.points <- rasterToPoints(data.rasters.crop, spatial=FALSE)

# column names
#colnames(r.points) <- c("Lon", "Lat", "sst","icethick", "salinity", "chloro")

## 2.3 Extract climate data from rasters and transform environmental predictors
# start with present data (present.raster.crop for non-range specific, and present.rasters.range for range-specific)
clim.land <- raster::extract(present.rasters.range, 1:ncell(present.rasters.range), df = TRUE)
clim.land <- na.omit(clim.land)

# Use predict function to transform environmental predictors
pred <- predict(gf, clim.land[,-1])  #note the removal of the cell ID column with [,-1])

# Convert predictions to color scale for map
# use pca on predictions
pca <- prcomp(pred, center=T, scale.=F)

# Assign pcs to colors
r <- pca$x[, 1]
g <- pca$x[, 2]
b <- pca$x[, 3]

# Scale colors
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255

# Define raster properties with an existing one
mask<-present.rasters.crop$salinity 
#mask[]<-as.numeric(mask[]>0)

# Assign color to raster
rastR <- rastG <- rastB <- mask
rastR[clim.land$ID] <- r
rastG[clim.land$ID] <- g
rastB[clim.land$ID] <- b

# Stack color raster
rgb.rast <- stack(rastR, rastG, rastB)

######
# Read in beluga coords
coord <- read.csv("beluga/beluga_coords_n140.csv")
pres.clim.points <- coord

# also site by location for later too
site <- readxl::read_xlsx("beluga/Sites_beluga_popinfo.xlsx")

## redo site load for narwhal/bowhead
sites <- read.csv("climate_data/site_locations.csv")

# narwhal
sites <- subset(sites, Narwhal > 0)
names(sites)[names(sites) == 'narwhal_pop'] <- 'Population'

# bowhead
sites <- subset(sites, Bowhead > 0)
sites$Population <- "ECA"

######
#####

#pdf("bowhead_RWmap_GF_Map.pdf")
svglite("beluga/results_within_range/beluga_GF_Map_updated_rangeclip.svg", width=6, height=6)
# initial map/plot
plotRGB(rgb.rast, colNA="gray80")
points(pres.clim.points$Lon, pres.clim.points$Lat, pch=19, col="black")
dev.off()

# PCA legend of biological space
nvs <- dim(pca$rotation)[1]
vec <- c("sst", "icethick", "salinity", "chloro")
lv <- length(vec)
vind <- rownames(pca$rotation) %in% vec
scal <- 40
xrng <- range(pca$x[, 1], pca$rotation[, 1]/scal) * 1.1
yrng <- range(pca$x[, 2], pca$rotation[, 2]/scal) * 1.1

svglite("beluga/results_within_range/beluga_GF_Map_updated_rangeclip_PCA_legend.svg", width=6, height=5)
plot((pca$x[, 1:2]), xlim = xrng, ylim = yrng, 
     pch = ".", cex = 4, col = rgb(r, g, b, max = 255),
     asp = 1)

arrows(rep(0, lv), rep(0, lv), pca$rotation[vec,1]/scal, pca$rotation[vec, 2]/scal, length = 0.0625)
jit <- 0.0015

text(pca$rotation[vec, 1]/scal + jit * sign(pca$rotation[vec, 1]), pca$rotation[vec, 2]/scal + jit * sign(pca$rotation[vec, 2]), labels = vec)
dev.off()

#From vignette:
#"The colors represent genetic variation (allelic composition) as predicted based on the modeled relationships with environmental and spatial variables. Similar colors are predicted to be more similar genetically."


##########  

##########  
##### Step 3: Estimate genetic offsets using future environmental data
# replaced "future.rasters.crop" with "future.rasters.range" to make range specific

## 3.1 prepare data files
clim.land.future <- raster::extract(future.rasters.range, 1:ncell(future.rasters.range), df = TRUE)
clim.land.future <- na.omit(clim.land.future)

# Transform environmental variables
pred.future <- predict(gf, clim.land.future[,-1])

# Noting here that changed pred.adaptive (tutorial) to pred (here).
# Genetic offset:
genetic.offset.adaptive <- sqrt((pred.future[,1]-pred[,1])^2 + 
                                  (pred.future[,2]-pred[,2])^2 + 
                                  (pred.future[,3]-pred[,3])^2 + 
                                  (pred.future[,4]-pred[,4])^2)

# Define raster properties --- the variable doesn't matter here i think
# "future.rasters.range" for species range, "future.rasters.crop" for full ECA
rast.offset <- future.rasters.range$salinity 

# Assign genetic offset values (difference between future and present predictions) to raster 
rast.offset[clim.land.future$ID] <- genetic.offset.adaptive

# Make color scale
offset.colors = colorRampPalette(c("#4575B4", "#abd9e9","#ffffbf","#fdae61","#f46d43","#a50026"))

# Quick plot!
pdf("bowhead/results_within_range/bowhead_RWmap_GF_GeneticOffset_RCP85_2100_rangeclip.pdf")
plot(rast.offset, col=offset.colors(200))
points(pres.clim.points$Lon, pres.clim.points$Lat)
dev.off()

# Add land to map
# Get map outlines from rworldmap package. For "high" resolution, need to also install package rworldxtra
map_outline <- getMap(resolution="high")

# Crop outline to boundary and convert to dataframe
map_outline <- crop(map_outline, y=crop.extent) %>% fortify()

# create SpatRaster to work with tidyterra
rast.offset.spat <- rast(rast.offset)

##### To add crop to species range at this stage
# without rgdal, use "sf::read_sf" to load shape file
#species <- read_sf('C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/maps/range_maps/redlist_species_data_beluga/data_0.shp')

#rast.offset.clip <- mask(rast.offset, species)
#rast.offset.clip.spat <- rast(rast.offset.clip)

#####

# plot and save as svg
svglite("narwhal/results_within_range/narwhal_GF_GeneticOffset_RCP85_2100_landmap_sites_step_updated2_shape_rangeclip.svg", width=6, height=6)
ggplot()+
  geom_spatraster(data=rast.offset.spat)+
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  
  # beluga points
  #geom_point(data=site, aes(x=Lon, y=Lat, shape=Population),size=2,colour="black")+ 
  #scale_shape_manual(values=c(15, 8, 17, 18, 12, 16))+ 
  
  # narwhal points
  geom_point(data=sites, aes(x=est_longitude, y=est_latitude, shape=Population),size=2,colour="black")+
  scale_shape_manual(values=c(17, 15, 18))+
  
  # bowhead points
  #geom_point(data=sites, aes(x=est_longitude, y=est_latitude, shape=Population),size=2,colour="black")+
  #scale_shape_manual(values=c(16))+
  
  #scale_fill_gradientn(colors=offset.colors(100))+
  scale_fill_stepsn(n.breaks = 20, colours = offset.colors(12),  na.value="white")+ #limits=c(0,0.11))+
  scale_x_continuous(expand=c(0,0), breaks=c(-100,-80,-60)) +
  scale_y_continuous(expand=c(0,0), breaks=c(50,60,70,80)) +
  #scale_x_continuous(expand=c(0,0), breaks=c(-105,-50)) +
  #scale_y_continuous(expand=c(0,0), breaks=c(46,78)) +
  xlab("Longitude")+
  ylab("Latitude")+
  labs(fill="Genetic offset")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.title = element_text(hjust = 0.5))+
 # annotation_scale(height=unit(0.15, "cm"), location="tr", aes(width_hint=0.15), text_cex=0.8)+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
                         location = "tr", which_north = "true",
                         pad_x = unit(0.12, "cm"), pad_y = unit(0.3, "cm"), #pad_y unit 0.6 if including scale bar
                         style = ggspatial::north_arrow_fancy_orienteering())+
  theme(text=element_text(size=15))+
  theme(legend.title = element_text( size=12), legend.text=element_text(size=12))
dev.off()



# Save the rast.offset for later too
raster::writeRaster(rast.offset, "beluga/results_within_range/raster.offset.RCP45.2100.beluga.easternCanada.RANGECLIP.tif")

# Can import raster back in like this
#rast.offset.test <- raster("raster.offset.RCP85.2100.narwhal.easternCanada.tif")


##########  

##########  
##### Step 4: Creating box plots - extract genetic offset values

# Treat rast.offset as a raster to extract values by site.
# load site locations
sites <- read.csv("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current/climate_data/site_locations.csv") %>% subset(Narwhal > 0) #%>% subset(Location_abbrev != "GH")

x <- sites$est_longitude
y <- sites$est_latitude

#beluga
x <- site$Lon
y <- site$Lat

pts <- cbind(x, y) %>% as_tibble()
pts <- SpatialPoints(pts, proj4string = CRS("+proj=longlat +datum=WGS84"))
projection(pts)==projection(rast.offset) # make sure projection matches

# Create dataframe to store extracted data
offset.data <- tibble(ID = 1:nrow(pts@coords), Lon = pts$x, Lat = pts$y)

# Extract for each layer (buffer is in meters) 
# try fun=NULL to avoid summarizing
store_data = list()
for (i in 1:nlayers(rast.offset)){
  store_data[[i]] = raster::extract(rast.offset[[i]], pts, buffer=100000, fun=NULL)
}
# Now have genetic offset values for 100km for each site.

# Code below for narwhal, bowhead, belugs separately 

##### Narwhal
# 1=AB, 2=BI, 3=CR, 4=GF, 5=IG, 6=PG, 7=PB, 8=PI, 9=RB, 10=RE, 11=SB.
# make colnames CHA, BI, NHB for subgroups

AB <- as.data.frame(store_data[[1]][[1]])
colnames(AB) <- "AB"
BI <- as.data.frame(store_data[[1]][[2]])
colnames(BI) <- "BI"
CR <- as.data.frame(store_data[[1]][[3]])
colnames(CR) <- "CR"
GF <- as.data.frame(store_data[[1]][[4]])
colnames(GF) <- "GF"
IG <- as.data.frame(store_data[[1]][[5]])
colnames(IG) <- "IG"
PG <- as.data.frame(store_data[[1]][[6]])
colnames(PG) <- "PG"
PB <- as.data.frame(store_data[[1]][[7]])
colnames(PB) <- "PB"
PI <- as.data.frame(store_data[[1]][[8]])
colnames(PI) <- "PI"
RB <- as.data.frame(store_data[[1]][[9]])
colnames(RB) <- "RB"
RE <- as.data.frame(store_data[[1]][[10]])
colnames(RE) <- "RE"
SB <- as.data.frame(store_data[[1]][[11]])
colnames(SB) <- "SB"

# combining data frames with diff. some sites might have NAs, but because the 100km buffer might extend over land 

CHApop <- rbind.fill(GF, RE, SB)
CHApop$pop <- "CHA"

BIpop <- rbind.fill(AB, BI, CR, IG, PG, PB, PI)
BIpop$pop <- "BI"

NHBpop <- rbind.fill(RB)
NHBpop$pop <- "NHB"

combined <- rbind.fill(CHApop, BIpop, NHBpop)
meltData <- melt(combined)

boxplot(data=meltData, value~pop)

narwhal_offset <- meltData %>%
  mutate(pop = fct_reorder(pop, value, .fun="mean")) %>%
  ggplot(aes(x=pop, y=value, fill=pop)) + 
  geom_boxplot(width=0.5, lwd=0.5, color="gray30", fill="#fdf391ff",alpha=1, outlier.size = 0.5) +
  theme(legend.position="none")+
  ylab("Genetic offset")+
  theme_bw()+
  theme(legend.position="none", panel.grid.major.x=element_blank())+
  xlab("Population")

narwhal_offset
#write.csv(meltData, "narwhal/narwhal_meltdata_for_boxplot_RCP85_2100.csv")
ggsave("narwhal/results_within_range/narwhal_boxplot_RCP85_2100_bypop_rangeclip.svg", width=3, height=4, dpi=1000)

# if summary doesn't work use CHApop %>% mutate(new=mean(c(GF, RE, SB), na.rm=T))
#ggsave("narwhal/narwhal_offset_boxplot_RCP85_2100_bypop_yellow.svg", width=2.5, height=4, dpi=900)

save.image("narwhal_K3_GF_snps_scan_4env_ntree500_nbin301_offset_fullpipeRCP85_2100_rangeclip.RData")


##### Bowhead
# 1=AB, 2=CD, 3=CR, 4=CH, 5=GH, 6=HB, 7=IQ, 8=PG, 9=PB, 10=PI, 11=RI, 12=RB, 13=SB, 14=WB, 15=DIS 

# make colnames MID, NOR, RB for subgroups

AB <- as.data.frame(store_data[[1]][[1]])
colnames(AB) <- "AB"
CD <- as.data.frame(store_data[[1]][[2]])
colnames(CD) <- "CD"
CR <- as.data.frame(store_data[[1]][[3]])
colnames(CR) <- "CR"
CH <- as.data.frame(store_data[[1]][[4]])
colnames(CH) <- "CH"
GH <- as.data.frame(store_data[[1]][[5]])
colnames(GH) <- "GH"
HB <- as.data.frame(store_data[[1]][[6]])
colnames(HB) <- "HB"
IQ <- as.data.frame(store_data[[1]][[7]])
colnames(IQ) <- "IQ"
PG <- as.data.frame(store_data[[1]][[8]])
colnames(PG) <- "PG"
PB <- as.data.frame(store_data[[1]][[9]])
colnames(PB) <- "PB"
PI <- as.data.frame(store_data[[1]][[10]])
colnames(PI) <- "PI"
RI <- as.data.frame(store_data[[1]][[11]])
colnames(RI) <- "RI"
RB <- as.data.frame(store_data[[1]][[12]])
colnames(RB) <- "RB"
SB <- as.data.frame(store_data[[1]][[13]])
colnames(SB) <- "SB"
WB <- as.data.frame(store_data[[1]][[14]])
colnames(WB) <- "WB"
DIS <- as.data.frame(store_data[[1]][[15]])
colnames(DIS) <- "DIS"

# combining data frames with diff. lengths will add a lot of NAs but is fine since we  want raw values 

# 1=AB, 2=CD, 3=CR, 4=CH, 5=GH, 6=HB, 7=IQ, 8=PG, 9=PB, 10=PI, 11=RI, 12=RB, 13=SB, 14=WB, 15=DIS 

# combine
combined <- rbind.fill(AB, CD, CR, CH, GH, HB, IQ, PG, PB, PI, RI, RB, SB, WB, DIS)
combined$pop <- "ECA"
meltData <- melt(combined)
boxplot(data=meltData, value~variable)

# melt data by each site # by longitude order for bowhead
#meltData$variable <- factor(meltData$variable , levels=c("GH", "SB", "RI", "PB", "RB", "AB", "CH", "HB", "CD", "PI", "WB", "CR", "PG", "IQ", "DIS")) 

bowhead_offset <- meltData %>%
  ggplot(aes(x=pop, y=value, fill=pop)) + 
  geom_boxplot(width=0.5, lwd=0.5, color="gray30", fill="#81c3f4ff",alpha=1, outlier.size = 0.5) +
  theme(legend.position="none")+
  ylab("Genetic offset")+
  theme_bw()+
  theme(legend.position="none", panel.grid.major.x=element_blank())+
  xlab("Population")

bowhead_offset
#write.csv(meltData, "bowhead/bowhead_meltdata_for_boxplot_RCP85_2100.csv")
ggsave("bowhead/bowhead_boxplot_RCP45_2100_bypop.svg", width=3, height=4, dpi=1000)
# combined %>% mutate(new=mean(c(AB, CD, CR, CH, GH, HB, IQ, PG, PB, PI, RI, RB, SB, WB, DIS), na.rm=T))
#ggsave("bowhead/bowhead_offset_boxplot_RCP85_2100_bypop_lighterblue.svg", width=1.7, height=4, dpi=900)

save.image("bowhead_RWmap_GF_snps_ntree100_nbin601_offset_fullpipeRCP45_2100.RData")


##### Beluga:
####################### part from Claudio (lines 632-803):
#####
# 1=CS (1), 2=EHA (2), 3=EHB (3, 6-7), 4=LWR (4-5), 5=HB (8), 6=Iqaluit (9), 7=JB (10-12),
#8=LanseAuxLoup (13), 9=SL (14-31), 10=UB (32-33), 11=WHB (34)

CS <- as.data.frame(store_data[[1]][[1]])
CS$location <- "CS"
CS$site <- "CS"
colnames(CS)[1] <- "data"

EHA <- as.data.frame(store_data[[1]][[2]])
EHA$location <- "EHA"
EHA$site <- "EHA"
colnames(EHA)[1] <- "data"

EHB1 <- as.data.frame(store_data[[1]][[3]])
EHB1$location <- "EHB1"
EHB1$site <- "EHB"
colnames(EHB1)[1] <- "data"

LWR2 <- as.data.frame(store_data[[1]][[4]])
LWR2$location <- "LWR2"
LWR2$site <- "LWR"
colnames(LWR2)[1] <- "data"

LWR <- as.data.frame(store_data[[1]][[5]])
LWR$location <- "LWR"
LWR$site <- "LWR"
colnames(LWR)[1] <- "data"

EHB4 <- as.data.frame(store_data[[1]][[6]])
EHB4$location <- "EHB4"
EHB4$site <- "EHB"
colnames(EHB4)[1] <- "data"

EHB5 <- as.data.frame(store_data[[1]][[7]])
EHB5$location <- "EHB5"
EHB5$site <- "EHB"
colnames(EHB5)[1] <- "data"

HB <- as.data.frame(store_data[[1]][[8]])
HB$location <- "HB"
HB$site <- "HB"
colnames(HB)[1] <- "data"

Iqaluit <- as.data.frame(store_data[[1]][[9]])
Iqaluit$location <- "Iqaluit"
Iqaluit$site <- "Iqaluit"
colnames(Iqaluit)[1] <- "data"

JB <- as.data.frame(store_data[[1]][[10]])
JB$location <- "JB"
JB$site <- "JB"
colnames(JB)[1] <- "data"

JB2 <- as.data.frame(store_data[[1]][[11]])
JB2$location <- "JB2"
JB2$site <- "JB"
colnames(JB2)[1] <- "data"

JB3 <- as.data.frame(store_data[[1]][[12]])
JB3$location <- "JB3"
JB3$site <- "JB"
colnames(JB3)[1] <- "data"

LaL <- as.data.frame(store_data[[1]][[13]])
LaL$location <- "LaL"
LaL$site <- "LaL"
colnames(LaL)[1] <- "data"

SL1 <- as.data.frame(store_data[[1]][[14]])
SL1$location <- "SL1"
SL1$site <- "SL"
colnames(SL1)[1] <- "data"

SL2 <- as.data.frame(store_data[[1]][[15]])
SL2$location <- "SL2"
SL2$site <- "SL"
colnames(SL2)[1] <- "data"

SL3 <- as.data.frame(store_data[[1]][[16]])
SL3$location <- "SL3"
SL3$site <- "SL"
colnames(SL3)[1] <- "data"

SL4 <- as.data.frame(store_data[[1]][[17]])
SL4$location <- "SL4"
SL4$site <- "SL"
colnames(SL4)[1] <- "data"

SL5 <- as.data.frame(store_data[[1]][[18]])
SL5$location <- "SL5"
SL5$site <- "SL"
colnames(SL5)[1] <- "data"

SL6 <- as.data.frame(store_data[[1]][[19]])
SL6$location <- "SL6"
SL6$site <- "SL"
colnames(SL6)[1] <- "data"

SL7 <- as.data.frame(store_data[[1]][[20]])
SL7$location <- "SL7"
SL7$site <- "SL"
colnames(SL7)[1] <- "data"

SL8 <- as.data.frame(store_data[[1]][[21]])
SL8$location <- "SL8"
SL8$site <- "SL"
colnames(SL8)[1] <- "data"

SL9 <- as.data.frame(store_data[[1]][[22]])
SL9$location <- "SL9"
SL9$site <- "SL"
colnames(SL9)[1] <- "data"

SL10 <- as.data.frame(store_data[[1]][[23]])
SL10$location <- "SL10"
SL10$site <- "SL"
colnames(SL10)[1] <- "data"

SL11 <- as.data.frame(store_data[[1]][[24]])
SL11$location <- "SL11"
SL11$site <- "SL"
colnames(SL11)[1] <- "data"

SL12 <- as.data.frame(store_data[[1]][[25]])
SL12$location <- "SL12"
SL12$site <- "SL"
colnames(SL12)[1] <- "data"

SL13 <- as.data.frame(store_data[[1]][[26]])
SL13$location <- "SL13"
SL13$site <- "SL"
colnames(SL13)[1] <- "data"

SL14 <- as.data.frame(store_data[[1]][[27]])
SL14$location <- "SL14"
SL14$site <- "SL"
colnames(SL14)[1] <- "data"

SL15 <- as.data.frame(store_data[[1]][[28]])
SL15$location <- "SL15"
SL15$site <- "SL"
colnames(SL15)[1] <- "data"

SL16 <- as.data.frame(store_data[[1]][[29]])
SL16$location <- "SL16"
SL16$site <- "SL"
colnames(SL16)[1] <- "data"

SL17 <- as.data.frame(store_data[[1]][[30]])
SL17$location <- "SL17"
SL17$site <- "SL"
colnames(SL17)[1] <- "data"

SL18 <- as.data.frame(store_data[[1]][[31]])
SL18$location <- "SL18"
SL18$site <- "SL"
colnames(SL18)[1] <- "data"

UB1 <- as.data.frame(store_data[[1]][[32]])
UB1$location <- "UB1"
UB1$site <- "UB"
colnames(UB1)[1] <- "data"

UB2 <- as.data.frame(store_data[[1]][[33]])
UB2$location <- "UB2"
UB2$site <- "UB"
colnames(UB2)[1] <- "data"

WHB <- as.data.frame(store_data[[1]][[34]])
WHB$location <- "WHB"
WHB$site <- "WHB"
colnames(WHB)[1] <- "data"

#######################################
# Combine so it's by population instead of sampling site
EHApop <- rbind.fill(EHA)
EHApop$pop <- "EHA"

HBpop <- rbind.fill(Iqaluit,EHB1, EHB4, EHB5, HB, WHB, UB1, UB2, LaL)
HBpop$pop <- "HB"

LWRpop <- rbind.fill(LWR, LWR2)
LWRpop$pop <- "LWR"

CSpop <- rbind.fill(CS)
CSpop$pop <- "CS"

JBpop <- rbind.fill(JB, JB2, JB3)
JBpop$pop <- "JB"

SLpop <- rbind.fill(SL1,SL10, SL11, SL12, SL13, SL14, SL15, SL16, SL17, SL18, SL2, SL3, SL4, SL5, SL6, SL7, SL8, SL9)
SLpop$pop <- "SL"

combined <- rbind.fill(EHApop, HBpop, LWRpop, CSpop, JBpop, SLpop)
meltData <- melt(combined)
boxplot(data=meltData, value~pop)

# box plot
beluga_offset <- meltData %>%
  mutate(pop = fct_reorder(pop, value, .fun="mean")) %>%
  ggplot(aes(x=pop, y=value, fill=pop)) + 
  geom_boxplot(width=0.5, lwd=0.5, color="gray30", fill="#a185a9ff", alpha=1, outlier.size = 0.5) +
  theme(legend.position="none")+
  ylab("Genetic offset")+
  theme_bw()+
  theme(legend.position="none", panel.grid.major.x=element_blank())+
  xlab("Population")

beluga_offset

ggsave("beluga/results_within_range/beluga_K6_boxplot_RCP45_2100_bypop_rangeclip.svg", width=3, height=4, dpi=1000)

#write.csv(meltData, "beluga/beluga_K4_5pop_meltdata_for_boxplot_RCP85_2100.csv")
#ggsave("beluga/beluga_K4_offset_boxplot_RCP85_2100_bypop_purple.svg", width=3.5, height=4, dpi=900)


######
## all boxplots together in one -- did not end up using this in the end, but leaving here in case useful later. tried two different ways (not facet and facet)

beluga_boxdata <- read.csv("beluga/beluga_K4_5pop_meltdata_for_boxplot_RCP85_2100.csv")
beluga_boxdata <- beluga_boxdata %>% select(pop, variable, value)
beluga_boxdata$sp <- "Beluga"
beluga_boxdata$group <- paste(beluga_boxdata$sp, beluga_boxdata$pop, sep="_")

narwhal_boxdata <- read.csv("narwhal/narwhal_meltdata_for_boxplot_RCP85_2100.csv")
narwhal_boxdata <- narwhal_boxdata %>% select(pop, variable, value)
narwhal_boxdata$sp <- "Narwhal"
narwhal_boxdata$group <- paste(narwhal_boxdata$sp, narwhal_boxdata$pop, sep="_")

bowhead_boxdata <- read.csv("bowhead/bowhead_meltdata_for_boxplot_RCP85_2100.csv")
bowhead_boxdata <- bowhead_boxdata %>% select(pop, variable, value)
bowhead_boxdata$sp <- "Bowhead"
bowhead_boxdata$group <- paste(bowhead_boxdata$sp, bowhead_boxdata$pop, sep="_")

meltData_combined <- rbind(beluga_boxdata, narwhal_boxdata, bowhead_boxdata)

ggplot(meltData_combined %>% #arrange(sp, desc(value)) %>%
        mutate(group=factor(group, levels=(unique(group)))), 
      aes(x=group,y=value, fill = sp))+
#  mutate(pop = fct_reorder(pop, value, .fun="mean")) %>%
 # ggplot(aes(x=pop, y=value, fill=pop)) + 
  scale_x_discrete(limits=c("Beluga_CS", "Beluga_SL", "Beluga_HB", "Beluga_JB", "Beluga_EHA", "Narwhal_BI", "Narwhal_CHA", "Narwhal_NHB", "Bowhead_ECA"))+
  geom_boxplot(width=0.5, lwd=0.5, color="gray30",alpha=1, outlier.size = 0.5) +
  theme(legend.position="none")+
  ylab("Genetic offset")+
  theme_bw()+
  theme(legend.position="none", panel.grid.major.x=element_blank())+
  #scale_fill_manual(values=c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61","#D7191C"))+
    scale_fill_manual(values=c("#a185a9ff", "#81c3f4ff", "#fdf391ff"))+
  xlab("Population")

#ggsave("all3sp_box_boxplot_RCP85_2100_bypop_.svg", width=6, height=3, dpi=1000)

library(tidyr)

ggplot(meltData_combined %>% #arrange(pop, value) %>%
           mutate(pop = fct_reorder(pop, value, .fun="mean")),
            #mutate(pop=factor(pop, levels=(unique(pop)))), 
       aes(x=pop,y=value, order=pop,fill = sp))+
  geom_boxplot(width=0.5, lwd=0.5, color="gray30",alpha=1, outlier.size = 0.5) +
  theme(legend.position="none")+
  facet_grid(. ~ sp, scales="free_x", space="free_x")+
  ylab("Genetic offset")+
  theme_bw()+
  theme(legend.position="none", panel.grid.major.x=element_blank(), panel.spacing=unit(2, "lines"))+
  #scale_fill_manual(values=c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61","#D7191C"))+
  scale_fill_manual(values=c("#a185a9ff", "#81c3f4ff", "#fdf391ff"))+
  xlab("Population")+
  scale_y_continuous(limits = c(0,0.04), expand = c(0, 0)) +
ggsave("all3sp_box_boxplot_RCP85_2100_bypop_FACET.svg", width=8, height=3, dpi=1000)

