# Jump to line 203 for isolation-by-distance part

# 1) Make narwhal map with distances

# modules
library(raster)
library(rworldmap)
library(rworldxtra)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggsn)
library(rgeos)
library(marmap)

setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/narwhal/fst/marmap")

# Set map boundary (xmin, xmax, ymin, ymax)
boundary = extent(-105,-50, 55, 85)
boundary

# Get map outlines from rworldmap package
map.outline = getMap(resolution = "high")

# Crop to boundary and convert to dataframe
map.outline = crop(map.outline, y = boundary) %>% fortify()

# Download ocean depth map. Lower resolution number is more fine-scale.
ocean_map_0 <- getNOAA.bathy(lon1 = -105, lon2 = -50, lat1 = 55, lat2 = 85, resolution = 5)

ocean_map <- ocean_map_0
dim(ocean_map)

# Add in the line in the Northwestern passages b/c narwhal range doesn't connect here
ocean_map[99900:100050] <-1000
ocean_map[99250:99400] <-1000
ocean_map[98600:98750] <-1000
ocean_map[97950:98100] <-1000
ocean_map[97300:97450] <-1000

autoplot(ocean_map, geom=c("raster")) + 
  scale_fill_stepsn(n.breaks=5,colors=c("#125ca1", "#70AED3","#BEDAEC","white", "gray"))


###
# Load site info
site_info <-  read.csv("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/narwhal/site_coords_only_narwhal.csv", header=T)

# Set site point colors
#site_manual_fill <- c("#4575B4", "#ABD9E9", "#FEE090", "#F46D43", "#A50026")

# Plot with ggplot2
nar_ocean_map <- autoplot(ocean_map, geom=c("raster")) + 
  scale_fill_stepsn(n.breaks=5,colors=c("#125ca1", "#70AED3","#BEDAEC","white", "gray"))+
  #geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), 
   #            fill="gray75",colour="gray40", size=0.5)+
  geom_point(data=site_info, aes(x=est_longitude, y=est_latitude),size=2, colour="green")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

nar_ocean_map

# Can stop here, but want to calculate distances between sites outside continental shelf
# Calculate distance with -500 depth minimum, making path impossible in waters shallower than -500 meters depth


site_coords <- site_info[,c("est_longitude", "est_latitude")]
colnames(site_coords) <- c("x", "y")

trans <- trans.mat(ocean_map,min.depth=-5)
out <- lc.dist(trans,site_coords,res="path")

# 'out' has lists for all the points in distance lines
# Extract the output for adding to ggplot map
A <- as.data.frame(out[1])
B <- as.data.frame(out[2])
C <- as.data.frame(out[3])
D <- as.data.frame(out[4])
E <- as.data.frame(out[5])
F <- as.data.frame(out[6])
G <- as.data.frame(out[7])
H <- as.data.frame(out[8])
I <- as.data.frame(out[9])
J <- as.data.frame(out[10])
K <- as.data.frame(out[11])
L <- as.data.frame(out[12])
M <- as.data.frame(out[13])
N <- as.data.frame(out[14])
O <- as.data.frame(out[15])
P <- as.data.frame(out[16])
Q <- as.data.frame(out[17])
R <- as.data.frame(out[18])
S <- as.data.frame(out[19])
T <- as.data.frame(out[20])
U <- as.data.frame(out[21])
V <- as.data.frame(out[22])
W <- as.data.frame(out[23])
X <- as.data.frame(out[24])
Y <- as.data.frame(out[25])
Z <- as.data.frame(out[26])
AA <- as.data.frame(out[27])
AB <- as.data.frame(out[28])

AC <- as.data.frame(out[29])
AD <- as.data.frame(out[30])
AE <- as.data.frame(out[31])
AF <- as.data.frame(out[32])
AG <- as.data.frame(out[33])
AH <- as.data.frame(out[34])
AI <- as.data.frame(out[35])
AJ <- as.data.frame(out[36])
AK <- as.data.frame(out[37])
AL <- as.data.frame(out[38])
AM <- as.data.frame(out[39])
AN <- as.data.frame(out[40])
AO <- as.data.frame(out[41])
AP <- as.data.frame(out[42])
AQ <- as.data.frame(out[43])
AR <- as.data.frame(out[44])
AS <- as.data.frame(out[45])


# Add the distance lines to map (adding geom_point again at the end to be on top of the lines)

nar_ocean_map2 <- autoplot(ocean_map_0, geom=c("raster")) + 
  scale_fill_stepsn(n.breaks=5,colors=c("#125ca1", "#70AED3","#BEDAEC","white", "gray"))+
  #geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), 
  #            fill="gray75",colour="gray40", size=0.5)+
  geom_point(data=site_info, aes(x=est_longitude, y=est_latitude),size=2, colour="green")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

nar_ocean_map2 +
  geom_path(data=A, aes(x=x, y=y), color="black", lwd=0.7,linetype="solid")+
  geom_path(data=B, aes(x=x, y=y), color="black", lwd=0.7,linetype="solid")+
  geom_path(data=C, aes(x=x, y=y), color="black", lwd=0.7,linetype="solid")+
  geom_path(data=D, aes(x=x, y=y), color="black", lwd=0.7,linetype="solid")+
  geom_path(data=E, aes(x=x, y=y), color="black", lwd=0.7,linetype="solid")+
  geom_path(data=F, aes(x=x, y=y), color="black", lwd=0.7,linetype="solid")+
  geom_path(data=H, aes(x=x, y=y), color="black", lwd=0.7,linetype="solid")+
  geom_path(data=J, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=K, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=L, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=M, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=N, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=O, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=P, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=Q, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=R, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=S, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=T, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=U, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=V, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=W, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=X, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=Y, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=Z, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AA, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AB, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AC, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AD, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AE, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AF, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AG, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AH, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AI, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AJ, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AK, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AL, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AM, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AN, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AO, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AP, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AQ, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
 geom_path(data=AR, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+
  geom_path(data=AS, aes(x=x, y=y), color="black", lwd=0.7, linetype="solid")+

  geom_point(data=site_info, aes(x=est_longitude, y=est_latitude),size=2, colour="yellow")
  

#ggsave("8sites/site_map_lines_NAR_base_8sites.png", width=8, height=8, dpi=1000)

# Actual distance to use in IBD analyses
dist <- lc.dist(trans,site_coords,res="dist")
dist

# Convert to matrix
mat <- as.matrix(dist)

# units in km

# 1=AB, 2=BI, 3=CR, 4=GF, 5=PI, 6=RB, 7=RE, 8=SB --OLD
# 1=PG, 2=RB, 3=BI, 4=PB, 5=SB, 6=IG, 7=CR, 8=PI, 9=AB, 10=RE, 12=GF
# rename numbers to locations

colnames(mat) <- c("PG", "RB", "BI", "PB", "SB", "IG", "CR", "PI", "AB", "RE", "GF")
rownames(mat) <- c("PG", "RB", "BI", "PB", "SB", "IG", "CR", "PI", "AB", "RE", "GF")

# Save distance matrix to csv#
write.csv(mat, "distance_matrix_NAR_11sites_noNPWbarrier.csv")


# see if we have IBD

# 2) Isolation-by-distance
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(RColorBrewer)
library(ggrepel)

# Load in data (in data, each pop pair is a separate row, with FST values and distance as columns)
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/narwhal/fst/reich")

data <- read.csv("Reichs_fst_output_narwhal_10sites_locationcorrected_wdist.csv")#, header=T)
data$pair


# since we have some sites with very low sample size, put all put

# IBD without PG and with NWPbarrier

# Also make highlight set for sites with less than 3 samples
data_highlight_main <- data %>% filter(less3!="X")
data_highlight_smol <- data %>% filter(less3=="X")


# Plot IBD with distance
IBD <- ggplot(data, aes(x=dist_km_NWPbarrier,y=Reich_fst))+
    geom_point(colour="gray50",size=2,alpha=0.8) +
geom_smooth(method = "lm", se=TRUE, fill="#9ec2fb", formula=y~x, linetype=1, color="#2953d1")+
 # geom_pointrange(aes(ymin=ll, ymax=ul), size=1, color="black")+
  theme_classic()+
  #ggtitle("Isolation by distance - Narwhal")+
  xlab("Distance (km)")+
  ylab(expression(paste(italic("F")~""[ST])))+
  theme(text=element_text(family="serif"))
 #geom_point(data=data_highlight_smol, aes(x=dist_km_NWPbarrier,y=Reich_fst), color="red",alpha=0.8)+
 # geom_point(data=data_highlight_main, aes(x=dist_km_NWPbarrier,y=Reich_fst), color="royalblue",alpha=0.8)+
 # geom_smooth(data=data_highlight_smol, aes(x=dist_km_NWPbarrier,y=Reich_fst),method = "lm", se=TRUE, fill="pink", formula=y~x, linetype=1, color="darkred")+
 # geom_smooth(data=data_highlight_main, aes(x=dist_km_NWPbarrier,y=Reich_fst),method = "lm", se=TRUE, fill="lightblue", formula=y~x, linetype=1, color="darkblue")
  
IBD
ggsave("IBD_plot_reich_10sites_loccorrected_NWPbarrier.png", width=5, height=4, dpi=900)
#ggsave("8sites/IBD_plot_reich_hwefilter_forppt.png", width=5, height=4, dpi=1000)
ggsave("IBD_plot_reich_10sites_loccorrected_NWPbarrier_size3_color_serif.png", width=3.5, height=2.7, dpi=1200)


# 3) Look at mantel correlation
library(ade4)

# Load data (distance matrix, and fst matrix)
distance_data <- read.csv("distance_matrix_nobarrier_8sites_mantel_prep.csv")
fst_data <- read.csv("fst_matrix_8sites_mantel_prep_locationcorrected.csv")

distance_dist <- dist(distance_data)
fst_dist <- dist(fst_data)

# Fst and distance
distance_mantel <- mantel.rtest(distance_dist, fst_dist, nrepet=9999)
distance_mantel


