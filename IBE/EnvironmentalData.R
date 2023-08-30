##Raster data from merraclim on dryad 2.5 arcmin

library(velox)
library(raster)
library(sp)
library(rgeos)
library(ggbiplot)
library(FactoMineR)
library(vegan)
library(ggvegan)

rasterList <- list.files("./2_5m_mean_00s/")

points <- read.csv("NASamplinglocs.csv", stringsAsFactors = F)


bio1 <- raster(paste0("2_5m_mean_00s/",rasterList[1]), format="+proj=longlat +datum=WGS84")
bio2 <- raster(paste0("2_5m_mean_00s/",rasterList[2]), format="+proj=longlat +datum=WGS84")
bio3 <- raster(paste0("2_5m_mean_00s/",rasterList[3]), format="+proj=longlat +datum=WGS84")
bio4 <- raster(paste0("2_5m_mean_00s/",rasterList[4]), format="+proj=longlat +datum=WGS84")
bio5 <- raster(paste0("2_5m_mean_00s/",rasterList[5]), format="+proj=longlat +datum=WGS84")
bio6 <- raster(paste0("2_5m_mean_00s/",rasterList[6]), format="+proj=longlat +datum=WGS84")
bio7 <- raster(paste0("2_5m_mean_00s/",rasterList[7]), format="+proj=longlat +datum=WGS84")
bio8 <- raster(paste0("2_5m_mean_00s/",rasterList[8]), format="+proj=longlat +datum=WGS84")
bio9 <- raster(paste0("2_5m_mean_00s/",rasterList[9]), format="+proj=longlat +datum=WGS84")
bio10 <- raster(paste0("2_5m_mean_00s/",rasterList[10]), format="+proj=longlat +datum=WGS84")
bio11 <- raster(paste0("2_5m_mean_00s/",rasterList[11]), format="+proj=longlat +datum=WGS84")
bio12 <- raster(paste0("2_5m_mean_00s/",rasterList[12]), format="+proj=longlat +datum=WGS84")
bio13 <- raster(paste0("2_5m_mean_00s/",rasterList[13]), format="+proj=longlat +datum=WGS84")
bio14 <- raster(paste0("2_5m_mean_00s/",rasterList[14]), format="+proj=longlat +datum=WGS84")
bio15 <- raster(paste0("2_5m_mean_00s/",rasterList[15]), format="+proj=longlat +datum=WGS84")
bio16 <- raster(paste0("2_5m_mean_00s/",rasterList[16]), format="+proj=longlat +datum=WGS84")
bio17 <- raster(paste0("2_5m_mean_00s/",rasterList[17]), format="+proj=longlat +datum=WGS84")
bio18 <- raster(paste0("2_5m_mean_00s/",rasterList[18]), format="+proj=longlat +datum=WGS84")
bio19 <- raster(paste0("2_5m_mean_00s/",rasterList[19]), format="+proj=longlat +datum=WGS84")

climLayers <- list(bio1, bio2, bio3, bio4,bio5,bio6,bio7,bio8,bio9,bio10,bio11,bio12,bio13,bio14,bio15,bio16,bio17,bio18,bio19)


veloxClim <- list()
for (i in (1:19)){
  veloxClim[[i]] <- velox(climLayers[[i]])
  i <- i+1
}



coords <- na.omit( data.frame(points$decimallongitude, points$decimallatitude))
sp <- SpatialPoints(coords)

Valuesdist <- list()

climdistPoints <- data.frame(rnorm(114, 1,1))

##extract using polygon


##extract using points

for (i in (1:19)){
  climdistPoints[,i] <- veloxClim[[i]]$extract_points(sp)
}


head(climdistPoints)
colnames(climdistPoints) <- c("bio1", "bio2", "bio3", "bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")

rownames(climdistPoints) <- points$Filename

head(climdistPoints)

corEnvPoints <- cor(climdistPoints)
corEnvPoints



climPCA <- rda(climdistPoints)
plot(climPCA$Ybar[,1],climPCA$Ybar[,2])
contrib=climPCA$CA$v

#Isolate most important PC's
contrib2=contrib[,1:3]

xx=rowSums(abs(contrib2))
print(xx[order(xx,decreasing = T)])


##not keeping variables within .8 corr
keepClim <- c("bio1","bio4","bio7","bio12","bio13","bio14")

exportclim <- climdistPoints[,keepClim]
write.csv(exportclim, "climVar.csv")



