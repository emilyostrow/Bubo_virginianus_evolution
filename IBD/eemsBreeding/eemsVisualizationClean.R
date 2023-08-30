#from https://github.com/dipetkov/eems/blob/master/README.md
#points outer from http://www.birdtheme.org/useful/v3tool.html
#devtools::install_github("dipetkov/eems/plotting/rEEMSplots")
#install.packages("rworldxtra")
library(rEEMSplots)
library(rgdal)
library(rworldmap) 
library(rworldxtra)
library("ggplot2")
library("dplyr")
require(scatterpie)
#install.packages("remotes")
#remotes::install_github("dipetkov/reemsplots2")
#remotes::install_github("JeffWeinell/misc.wrappers")ropensci/rnaturalearthhires
#remotes::install_github("ropensci/rnaturalearthhires")
library(reemsplots2)
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

coord <- read.table("data.coord")
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
together$breeding <- points$BreedingByGonads
together <- together[points$BreedingByGonads==1,]

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

#eemsMap

##pussing one on top of the other
#following https://stackoverflow.com/questions/11508902/plotting-discrete-and-continuous-scales-in-same-ggplot
grid.newpage()
pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) ) 
print( eemsMap + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
print( NAmap + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
#dev.off()
