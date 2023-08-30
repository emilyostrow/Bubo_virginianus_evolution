library(rEEMSplots)
library(rgdal)
library(broom)
library(rworldmap) 
library(rworldxtra)
library("ggplot2")
library("dplyr")
require(scatterpie)
library(reemsplots2)
library(misc.wrappers)
library(rnaturalearthhires)
library(rnaturalearthdata)
library(rnaturalearth)
library(sp)
library(tibble)
library(sf)
library(pophelper)
library(grid)


####getting the data####
mcmcpath = "./sim3"
plotpath = "./plot/sim3rEEMSplots"

##from https://github.com/OKersten/PuffPopGen/blob/master/Nuclear/11.EEMS.md
##exporting files with coordinates the eems way
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"

coord <- read.table("datapath.coord")
coords_merc <- sp::spTransform(
  SpatialPoints(coord, CRS(projection_none)),
  CRS(projection_mercator)
)

# `coords_merc` is a SpatialPoints structure
# but we only need the coordinates themselves
coords_merc <- coords_merc@coords

# So to summarize start with points in one projection
coord
# And project them to the target projection
coords_merc

##the one that works##
eems.plots(mcmcpath, 
           plotpath, 
           longlat = T, 
           out.png = F, res = 800,
           plot.height = 8, plot.width = 7,
           add.outline = TRUE, col.outline = "blue", lwd.outline = 1.5,  
#           add.demes = TRUE, col.demes = "red", min.cex.demes = 1.5, max.cex.demes = 1.8, 
           projection.in = projection_none,
           projection.out = projection_mercator,
           add.map = TRUE, col.map = "black", lwd.map = 1,
           add.abline = TRUE, add.r.squared = TRUE,
          m.plot.xy = {points(coords_merc, col = "purple")},
          q.plot.xy = {points(coords_merc, col = "purple")}
)



####StructureCode####
points<- read.csv("../InitialSmmaryAnalyses/StructureOrderPoints.csv" )
world <- map_data("world")
states <- map_data("state")
sfiles <- "structureOutk3v1.txt"
single <- readQ(files=sfiles)
single <- single$structureOutk3v1
together <- data.frame(long=points$decimallongitude, lat=points$decimallatitude)
together$A <- single$Cluster1
together$B <- single$Cluster2
together$C <- single$Cluster3

##creating a LGM Layer based on https://r-graph-gallery.com/168-load-a-shape-file-into-r.html

my_spdf <- readOGR( 
  dsn= paste0(getwd(),"/digital_maps_02_all_other_files/") , 
  layer="lgm_global",
  verbose=FALSE
)
par(mar=c(0,0,0,0))
plot(my_spdf, col="#f2f2f2", bg="skyblue", lwd=0.25, border=0 )
# 'fortify' the data to get a dataframe format required by ggplot2

spdf_fortified <- tidy(my_spdf)

# Plot it

iceExtent <- ggplot() +
  geom_polygon(data = spdf_fortified, aes( x = long, y = lat, group = group), fill="#69b3a2", color="white", alpha=.3) +
  theme_void() +geom_polygon(data = states, aes(x=long, y = lat, group = group), 
                             fill = NA, color="black") +
  coord_fixed(xlim = c(-175,-65),  ylim = c(15, 70), ratio = 1.2)
iceExtent

##creating a transparent structure map
NAmap <- ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill = NA, color="black") +
  coord_fixed(xlim = c(-175,-65),  ylim = c(15, 70), ratio = 1.2)+
  geom_scatterpie(aes(x=long, y=lat, r=1),
                  data=together, cols=c("A", "B", "C"), alpha=.8)+
  theme_void()+
  theme(legend.position="bottom")
NAmap

##creating a transparent eems map
plots <- make_eems_plots(mcmcpath, longlat = T)

eemsMap <- plots$qrates02 +
  geom_polygon(data = world, aes(x=long, y = lat, group = group),fill=NA, color="black") +
  coord_fixed(xlim = c(-175,-65),  ylim = c(15, 70), ratio = 1.2)+
  theme_void()

eemsMap

##putting one on top of the other
#following https://stackoverflow.com/questions/11508902/plotting-discrete-and-continuous-scales-in-same-ggplot
##data from: Quaternary Glaciations – Extent and chronology, a closer look
##data from: https://booksite.elsevier.com/9780444534477/digital_maps.php
grid.newpage()
pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) ) 
print( eemsMap + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
print( iceExtent + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
print( NAmap + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )


