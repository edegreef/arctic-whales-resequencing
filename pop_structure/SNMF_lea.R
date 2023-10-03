# SNMF analysis using code from LEA manual, and plotting pie chart admixture map using help from Tom Jenkins' tutorial on github (https://github.com/Tom-Jenkins/admixture_pie_chart_map_tutorial)

#library(devtools)
#devtools::install_github("bcm-uga/LEA")
library(LEA)
library(reshape2)

# for making pie chart map:
library(raster)
library(rworldmap)
library(rworldxtra)
library(ggplot2)
library(ggpubr)

# Note: the snmf function doesn't do well with files on dropbox 

setwd("C:/Users/eveli/Desktop/LEA")

# Convert ped file to geno file (automatically outputs in working directory)
ped2geno("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.ped")

# Run snmf for K 1-5
project=NULL
project=snmf("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.geno", K=1:5, entropy=TRUE, repetitions=10, project="new")

# Can also load snmf project if saved
project = load.snmfProject("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.snmfProject")

# Reading in sample info
sample_info <- read.csv("C:/Users/eveli/Dropbox/Whales with Garroway/narwhal/sample_info_narwhal.csv", header=T)

# Narwhal - remove these (kin pair/duplicate)
sample_info <- subset(sample_info, temp_remove!="X")

# Bowhead - remove these 2 individuals (kin pair)
sample_info <- subset(sample_info, sample_ID!="99_01")
sample_info <- subset(sample_info, sample_ID!="ARBMGH_2002_001")

# Plot cross-entropy
plot(project, col="blue", pch=19, cex=1.2)

# Extract the cross-entropy of runs K1:5
ce_k1 = cross.entropy(project, K = 1)
ce_k2 = cross.entropy(project, K = 2)
ce_k3 = cross.entropy(project, K = 3)
ce_k4 = cross.entropy(project, K = 4)
ce_k5 = cross.entropy(project, K = 5)

# Save values in csv
ce_merged = as.data.frame(cbind(ce_k1, ce_k2, ce_k3, ce_k4, ce_k5))
write.csv(ce_merged, "C:/Users/eveli/Dropbox/Whales with Garroway/narwhal/LEA/cross-entropy_K1-5.csv")

# Can move to dropbox here if want to
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/narwhal/LEA")

# Basic barcharts, going to do for runs K=2, 3, 4, 5
# K2
best = which.min(cross.entropy(project, K = 2))
best
my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold", "gray")
barchart(project, K = 2, run = best,
         border = NA, space = 0, col = my.colors,
         xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .6)

# K3
best = which.min(cross.entropy(project, K = 3))
best
barchart(project, K = 3, run = best,
         border = NA, space = 0, col = my.colors,
         xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .6)

# K4
best = which.min(cross.entropy(project, K = 4))
best
barchart(project, K = 4, run = best,
         border = NA, space = 0, col = my.colors,
         xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .6)

# K5
best = which.min(cross.entropy(project, K = 5))
best
barchart(project, K = 5, run = best,
         border = NA, space = 0, col = my.colors,
         xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)

# Moving on to individual ancestry and site pie map
# Use run with the lowest cross-entropy (K=1 had best fit for me, but K=2 CE was very close so I'm going to try K=2)
best = which.min(cross.entropy(project, K = 3))
best

# Extract Q-matrix for the best run
qmatrix = as.data.frame(Q(project, K = 3, run = best))
head(qmatrix)

# Label column names of qmatrix
ncol(qmatrix)
cluster_names = c()
for (i in 1:ncol(qmatrix)){
  cluster_names[i] = paste("Cluster", i)
}
cluster_names
colnames(qmatrix) = cluster_names
head(qmatrix)

# Add individual IDs
qmatrix$Ind = sample_info$sample

# Add site IDs
qmatrix$Site = sample_info$location_ID
head(qmatrix)

write.csv(qmatrix, "qmatrix_K3.csv", row.names=FALSE)

# Convert dataframe to long format
qlong = melt(qmatrix, id.vars=c("Ind","Site"))
head(qlong)

write.csv(qlong, "qmatrix_K3_qlong.csv", row.names=FALSE)

# change order of sites by using the factor function
#site.order = c("Iceland","Arctic","Labrador","Newfoundland","Scotian_shelf")
#qlong$Site_ord = factor(qlong$Site, levels = site.order)

# adjust facet labels
#levels(qlong$Site)
#facet.labs = c("Iceland","Arctic","Labrador","Newfoundland","Scotian_shelf")
#levels(qlong$Site) = facet.labs
#levels(qlong$Site)


# Load in qlong if starting from this step
setwd("LEA/pop structure")
qlong <- read.csv("qmatrix_K2_qlong.csv")
qmatrix <- read.csv("qmatrix_K2_poporder.csv")


# Define colour palette
pal = colorRampPalette(c("#F64F2E", "#1EB4C4", "#807DBA"))
cols = pal(length(unique(qlong$variable)))

# Plot admixture barplot 
admix.bar = ggplot(data=qlong, aes(x=Ind, y=value, fill=variable))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = c(0,0))+
  facet_wrap(~Site, scales = "free", ncol = 2)+
  scale_fill_manual(values = cols)+
  ylab("Admixture proportion")+
  xlab("Individual")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(colour="black", size=12),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
admix.bar
ggsave("1.admixture_barplot_sites_K2.png", width=10, height=9, dpi=300)


# Next, making pie chart maps

# Calculate mean admixture proportions for each site
head(qmatrix)
clusters = grep("Cluster", names(qmatrix)) # indexes of cluster columns
avg_admix = aggregate(qmatrix[, clusters], list(qmatrix$Site), mean)

# Order alphabetically by site
avg_admix = avg_admix[order(as.character(avg_admix$Group.1)), ]
avg_admix

# Convert dataframe from wide to long format
avg_admix = melt(avg_admix, id.vars = "Group.1")
head(avg_admix)

# Define a function to plot pie charts using ggplot for each site
pie_charts = function(admix_df, site, cols){
  ggplot(data = subset(admix_df, Group.1 == site),
         aes(x = "", y = value, fill = variable))+
    geom_bar(width = 1, stat = "identity", colour = "black", show.legend = FALSE)+
    coord_polar(theta = "y")+
    scale_fill_manual(values = cols)+
    theme_void()
}

# Run function on one site
pie_charts(avg_admix, site = "AB", cols = cols)

# Apply function to all sites using for loop
# narwhal
subsites = sort(c("AB", "BI", "CR", "GF", "IG", "PB", "PG", "PI", "RB", "RE", "SB"))

# bowhead
subsites = sort(c("AB", "CD", "CH", "CR", "GR", "HB", "IQ", "MD", "NF", "PB", "PG", "PI", "RB", "RI", "SB", "SH", "WB", "WG"))


pies = list()
for (i in subsites){
  pies[[i]] = pie_charts(admix_df = avg_admix, site = i, cols = cols) 
}
pies


# Import csv file containing coordinates
#coords = read.csv("F:/NBW-me/reseq/lea_rerun/site_mean_coords.csv", header=T)
# or if using from sample_info df:
coords <- sample_info[,c("location_ID", "est_latitude", "est_longitude")]

# Remove duplicates
coords <- coords %>% distinct()

# Order alphabetically by site
coords = coords[order(coords$location_ID), ] 
coords


# Check order matches coords order
as.character(avg_admix$Group.1) == as.character(coords$location_ID)

# Set map boundary (xmin, xmax, ymin, ymax)
boundary = extent(-170, -10, 45, 85)# bowhead
boundary = extent(-105,-50, 55, 85)
boundary

# Get map outlines from rworldmap package
map.outline = getMap(resolution = "high")

# Crop to boundary and convert to dataframe
map.outline = crop(map.outline, y = boundary) %>% fortify()

# Plot basemap
basemap = ggplot()+
  geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="grey88",
               colour="black", size=0.5)+
  coord_quickmap(expand=F)+
  ggsn::north(map.outline, symbol = 10, scale = 0.06, location = "topleft")+
  #ggsn::scalebar(data = map.outline, dist = 1000, dist_unit = "km", height = 0.01,
  #               transform = TRUE, model = "WGS84", 
  #               location = "bottomleft", anchor = c(x = -70, y = 55),
  #               st.bottom = FALSE, st.size = 4, st.dist = 0.015)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(
    axis.text = element_text(colour="black", size=12),
    axis.title = element_text(colour="black", size=14),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
    panel.grid.minor = element_line(colour="grey90", size=0.5),
    panel.grid.major = element_line(colour="grey90", size=0.5),
    legend.text = element_text(size=12),
    legend.title = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    legend.position = "top"
  )
basemap


# Then add pie charts to basemap

# Extract coordinates for each site
coord.list = list()
for (i in subsites){
  coord.list[[i]] = c(subset(coords, location_ID == i)$est_longitude, subset(coords, location_ID == i)$est_latitude)
}
coord.list

# Define pie chart sizes
radius = 2.5

# Convert ggplot pie charts to annotation_custom layers
pies.ac = list()
for (i in 1:length(subsites)){
  pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies[[i]]),
                                   xmin = coord.list[[i]][[1]] - radius,
                                   xmax = coord.list[[i]][[1]] + radius,
                                   ymin = coord.list[[i]][[2]] - radius,
                                   ymax = coord.list[[i]][[2]] + radius)
}

# Add layers to basemap
pie.map = basemap + pies.ac
pie.map
ggsave("2.piecharts_map_K2_colours.png", width = 8, height = 10, dpi = 1000)

# Combine ggplots
ggarrange(admix.bar + theme(legend.position = "right") + labs(title = "Individual admixture proportions", tag = "A"),
          pie.map + labs(title = "Mean admixture proportions for each site", tag = "B"))
ggsave("3.admixture_bar_map_K2_colours.png", width = 15, height = 6, dpi = 700)
