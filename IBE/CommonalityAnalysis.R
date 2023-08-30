##following https://github.com/lukemusher/Southern_Amazon_cophylogeography/blob/main/Commonality%20Analysis/Galbula/CA.Galbula.logistic.least_cost.dxy.R

#=============================================================================================================#
# Script created by Luke Musher, edited by EN Ostrow
# Script created in version R 3.3.2 
# Adapted from scripts from Glenn Seeholzer (Seeholzer and Brumfield 2017) as well as
#     from Prunier et al. 2014 
# This script: runs multivariate logistic regressions testing for isolation-by-distance, 
#			isolation-by-environment, and isolation-by-barrier in Amazonian birds Where 
#			dispersal distance (least cost-path) is used in place of geographic (euclidiean) distance.  
#			Runs a binomial GLM and commonality analysis that subsamples and removes duplicate localities
#     assesing the relative important of these processes in shaping 
#			genetic structure of Galbula
#
# Bootstrap procedure modified from code in supplementary of Prunier et al. (2015), which was based on methods in Peterman et al. (2014)
#
#	Prunier JG, Colyn M, Legendre X, Nimon KF, Flamand MC (2015) Multicollinearity in spatial genetics: Separating the wheat from the chaff using commonality analyses. Molecular Ecology, 24, 263–283.
#	Peterman WE, Connette GM, Semlitsch RD, Eggert LS (2014) Ecological resistance surfaces predict fine-scale genetic differentiation in a terrestrial woodland salamander. Molecular Ecology, 23, 2402–2413.
#
# Usage notes: run line by line
#=============================================================================================================#
rm(list=ls())
#set working directory
library(descr) # addon for Commonality analysis
library(plyr)
library(StAMPP)
library(gdistance)
library(ecodist)
library(maptools)
library(hierfstat)
library(pcadapt)
library(vcfR)
library(adegenet)
source("supporting.functions/CA_logistic.R")
source("supporting.functions/CA_semiStdCoef.R")
source('supporting.functions/FUN.add.alpha.R', chdir = TRUE)
source('supporting.functions/MMRR.R', chdir = F)


buboVCF <- read.vcfR("m3MAC1Unlinkgq30.vcf.gz")
data <- read.csv("StructureOrderPoints.csv", stringsAsFactors = F)
buboGen <- vcfR2genlight(buboVCF)
pop(buboGen) <- rep(1, times=(length(colnames(buboVCF@gt))-1))
sample.div <- stamppNeisD(buboGen, pop = FALSE)

####	ENV
####	import bioclim and enviro rasters
files <- list.files(path="../EnvironmentalData/2_5m_mean_00s/", pattern='tif', full.names=TRUE)
merra = stack(files) #create raster stack of BIOCLIM rasters
envValues <- raster::extract(merra, data[,c('decimallongitude','decimallatitude')])
envdata = data.frame(pop=data[,'LongName'],as.data.frame(envValues)) #extract BIOCLIM data for each individual
tmp = aggregate(. ~ pop, data= envdata,mean,na.rm=T) #take mean for each locality
my.prc = prcomp(tmp[,-1], center=T, scale = T, na.omit=T) #do DXY
pc.values = predict(my.prc)#[,1:ncol(tmp[,-1])]
pop.env.raw = data.frame(pc.values)
rownames(pop.env.raw) = tmp[,1]



# ####	Generate geographic DISTANCE MATRIX
tmp = aggregate(. ~ LongName, data= data[,c('LongName','decimallongitude','decimallatitude')],mean,na.rm=T)
pop.loc = as.matrix(tmp[,c('decimallongitude','decimallatitude')])
rownames(pop.loc) = tmp[,'LongName']
DIST.raw = spDists(pop.loc,longlat=TRUE) ; colnames(DIST.raw) = rownames(DIST.raw) = rownames(pop.loc)

#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------ Create distance matrices of population data for predictors ------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#

# Libraries
library(yhat) # Commonality analysis

standardize = function(x){tmp=(x-mean(x,na.rm=T))/sd(x,na.rm=T);diag(tmp)=0;return(tmp)}

prep.standardize <- function(x){x=data.matrix(x)
x=x[lower.tri(x, diag = FALSE)]
x=(x-mean(x))/sqrt(var(x)) 
x}

prep <- function(x){x=data.matrix(x)
x=x[lower.tri(x, diag = FALSE)]
x}

################################################################################################
####	Raw Distance Matrices
################################################################################################
##not dxy but not changing all the variables
mDXY = as.matrix(dist(sample.div))
mDXY[upper.tri(mDXY)] = mDXY[lower.tri(mDXY)]
mDIST = as.matrix(DIST.raw)
mEnv = as.matrix(dist(pop.env.raw))

################################################################################################
####	Distance matrices formatted for Commonality Analysis
################################################################################################
####	z-transformation for commonality analysis
#DXY = as.numeric(prep(mDXY))

# identify end of first "valley" in the kernal density
# this is the threshold above which we define successful divergence
# below this threshold we consider values equivalent to failure to diverge

hist(DXY, prob=T, col="skyblue1",breaks=75); rug(DXY)
lines(density(DXY, adj=.25),lwd=3)
abline(v=-1.85,col="red",lwd=2,lty=2)
DXY[DXY> -1.85]<-1
DXY[DXY< -1.85]<-0

### z-transform
DIST = as.numeric(prep(mDIST))
ENV = as.numeric(prep(mEnv))

####new analysis####
################################################################################################
####	Lists of distance matrices formatted for Multiple Matrix Regression on Distance Matrices (MMRR, Wang 2013)
################################################################################################
#MMRR to get regression coefficients and significance values for multivariate model
#data for MMRR

#Do Commonality Analysis to get Unique, Common, and Total contribution of each predictor variable
#for univariate and regressions
predictors = list(DIST=mDIST,ENV=mEnv)
predictors = lapply(predictors,standardize)
x = predictors
resBootstrap.list = list()
fit = MMRR(Y = mDXY, X = x, nperm=1000) ###This is a multiple matrix regression without commonality analysis

################################################################################################
####	Distance matrices formatted for Commonality Analysis
################################################################################################
####	for commonality analysis
DXY = as.numeric(prep.standardize(mDXY))
DIST = as.numeric(prep.standardize(mDIST))
ENV = as.numeric(prep.standardize(mEnv))
#BAR = as.numeric(prep(mbar))

data=data.frame(DXY = DXY, DIST = DIST, ENV = ENV)
ca = regr(lm(DXY ~ DIST + ENV,data=data)) #this is the commonality analysis
