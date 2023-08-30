#BiocManager::install("LEA")
library(LEA)
library(stringr)
library(vcfR)
library(ggplot2)
library(adegenet)
library(dartR)
library(radiator)
library(SNPfiltR)

##creating a breeding only vcf
setwd("~/Dropbox/My Mac (Emilys-MacBook-Pro-2.local)/Documents/Research/Bubo_virginianus/Chapter2/FinalAllindsAnalyses/EnvironmentalData")
allInds <- read.vcfR("m3MAC1Unlinkgq30.vcf")
points <- read.csv("StructureOrderPoints.csv", stringsAsFactors = F, header = T)
#breedingInds <- points[points$BreedingByGonads==1,]
nameList <- c("FORMAT",points$FileName)
for (i in 2:length(nameList)){
  nameList[i] <- paste0("B_virginianus_",nameList[i])
}
nameList[16] <- paste0(nameList[16],"_Comb")
nameList[44] <- paste0(nameList[44],"_Comb")
nameList[61] <- paste0(nameList[61],"_Comb")

allInds@gt <- allInds@gt[,colnames(allInds@gt) %in% nameList]
allInds@gt<-allInds@gt[!duplicated(allInds@fix[,1]),]
allInds@fix<-allInds@fix[!duplicated(allInds@fix[,1]),]
allInds

missing_by_snp(allInds)
filter_biallelic(allInds)





#polySNPs <- is.polymorphic(allInds)
#polySNPs <- na.omit(polySNPs[polySNPs==TRUE])
#breedinggen <- vcfR2genlight(allInds)
#filteredBreedingGen <- gl.keep.loc(breedinggen,loc.list=names(polySNPs))
#write_genlight(filteredBreedingGen)
#genomic_converter(filteredBreedingGen, output = "vcf")
#gl2vcf(filteredBreedingGen, outfile="output.vcf", outpath = "./", verbose = 5)
#write.vcf(allInds, "BreedingFilteredBubo.vcf.gz")


##Importing and converting data
vcf2lfmm("m3MAC1Unlinkgq30.vcf")
vcf2geno("m3MAC1Unlinkgq30.vcf")
envData <- read.csv("climVar.csv", stringsAsFactors = F, header = T)
row.names(envData) <- envData[,1]
envData <- envData[,2:7]
#breedingenvData <- envData[points$BreedingByGonads==1,]
write.env(envData, "climVar.env")
#write.env(breedingenvData, "breedingclimVar.env")



pc = pca("m3MAC1Unlinkgq30.lfmm")

## ----results='hide'-----------------------------------------------------------
# Perfom Tracy-Widom tests on all eigenvalues.
# create file: tuto.tracyWidom - tracy-widom test information.  
tw = tracy.widom(pc)

## ----results='asis'-----------------------------------------------------------
# display p-values for the Tracy-Widom tests (first 5 pcs). 
tw$pvalues[1:5]

## ----fig.width=4, fig.height=4, echo=TRUE-------------------------------------
# plot the percentage of variance explained by each component
plot(tw$percentage)

## ----results='hide'-----------------------------------------------------------
# main options
# K = number of ancestral populations
# entropy = TRUE: computes the cross-entropy criterion, 
# CPU = 4 the number of CPUs.
project = NULL
project = snmf("m3MAC1Unlinkgq30.geno",
               K = 1:10, 
#                K=3,
               entropy = TRUE, 
               repetitions = 10,
               project = "new")

## ----fig.width=4, fig.height=4, echo=TRUE-------------------------------------
# plot cross-entropy criterion for all runs in the snmf project
plot(project, col = "blue", pch = 19, cex = 1.2)

## ----fig.width=10, fig.height=4, echo=TRUE------------------------------------
# select the best run for K = 4
best = which.min(cross.entropy(project, K = 10))
my.colors <- c("#40004b", "#00441b", 
               "#e7d4e8", "#d9f0d3", 
               "#9970ab", "#5aae61", 
               "#762a83", "#1b7837", 
               "#c2a5cf", "#a6dba0")
barchart(project, K = 10, run = best,
         border = NA, space = 0, 
         col = my.colors, 
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order), 
     labels = bp$order, las=1, 
     cex.axis = .4)

## ----fig.width=8, fig.height=6, echo=TRUE, results='hide'---------------------
# Population differentiation tests
p = snmf.pvalues(project, 
                 entropy = TRUE, 
                 ploidy = 2, 
                 K = 3)
pvalues = p$pvalues
par(mfrow = c(2,1))
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)

## -----------------------------------------------------------------------------
# creation of a genotypic matrix  with missing genotypes
#dat = as.numeric(tutorial.R)
#dat[sample(1:length(dat), 100)] <-  9
#dat <- matrix(dat, nrow = 50, ncol = 400)
#write.lfmm(dat, "genoM.lfmm")

## ----results='hide'-----------------------------------------------------------
project.missing = snmf("m3MAC1Unlinkgq30.lfmm", K = 3, entropy = TRUE, repetitions = 10,project = "new")

## -----------------------------------------------------------------------------
# select the run with the lowest cross-entropy value
best = which.min(cross.entropy(project.missing, K = 3))

# Impute the missing genotypes
impute(project.missing, "m3MAC1Unlinkgq30.lfmm", method = 'mode', K = 3, run = best)


# Proportion of correct imputation results
dat.imp = read.lfmm("m3MAC1Unlinkgq30.lfmm_imputed.lfmm")
#mean( tutorial.R[dat == 9] == dat.imp[dat == 9] )

## ----results='hide'-----------------------------------------------------------
# main options: 
# K the number of latent factors
# Runs with K = 6 and 5 repetitions.
project = NULL
project = lfmm("m3MAC1Unlinkgq30.lfmm", 
               "climVar.env", 
               K = 3, 
               repetitions = 5, 
               project = "new")

#The project is saved into :
#  m5unlinkedbreeding2_climVar.lfmmProject 

#To load the project, use:
#  project = load.lfmmProject("m5unlinkedbreeding2_climVar.lfmmProject")

#To remove the project, use:
 # remove.lfmmProject("m5unlinkedbreeding2_climVar.lfmmProject")

## -----------------------------------------------------------------------------
# compute adjusted p-values
p = lfmm.pvalues(project, K = 3, d=6)
pvalues = p$pvalues

## ----fig.width=8, fig.height=6, echo=TRUE-------------------------------------
# GWAS significance test
par(mfrow = c(2,1))
hist(pvalues, col = "lightblue")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .7)

## -----------------------------------------------------------------------------
for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("Expected FDR:", alpha))
  L = length(pvalues)
  
  # return a list of candidates with expected FDR alpha.
  # Benjamini-Hochberg's algorithm:
  w = which(sort(pvalues) < alpha * (1:L) / L)
  candidates = order(pvalues)[w]
  
  # estimated FDR and True Positive Rate
  Lc = length(candidates)
  estimated.FDR = sum(candidates <= 350)/Lc
  print(paste("Observed FDR:", 
              round(estimated.FDR, digits = 2)))    
  estimated.TPR = sum(candidates > 350)/50
  print(paste("Estimated TPR:", 
              round(estimated.TPR, digits = 2)))  
}

## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
# Fitting an LFMM with K = 3 factors
mod <- lfmm2("m3MAC1Unlinkgq30.lfmm_imputed.lfmm", "climVar.env", K = 3)

## ----fig.width=8, fig.height=6, echo=TRUE, results='hide'---------------------
# Computing P-values and plotting their minus log10 values 
pv <- lfmm2.test(object = mod, "m3MAC1Unlinkgq30.lfmm_imputed.lfmm", "climVar.env", linear = TRUE)


####plotting and identifying outlier SNPs####
alpha <- .05
L = length(pv$pvalues)
snippetys <- vector()

for (i in 1:9){
w = which(sort(pv$pvalues[i,]) < alpha * (1:L) / L)
SNPnames <- names(w)
candidates <- as.numeric(str_remove_all(SNPnames, "Response V"))

plot(-log10(pv$pvalues[i,]), col = "grey", cex = .6, pch = 19)
points(candidates, -log10(pv$pvalues[i,candidates]), col = "red")
title(main = paste0("Variable ", i))
snippetys <- c(snippetys, candidates)
}

snippetysReal <- sort(unique(snippetys))

####filtering by Clim SNPs####
m3 <- read.vcfR("m3MAC1Unlinkgq30.vcf")
SNPlist <- as.character(snippetysReal)

m3Clim <- m3
m3Clim@gt<-m3@gt[snippetysReal,]
m3Clim@fix<-m3@fix[snippetysReal,]
vcfR::write.vcf(m3Clim, "clim.vcf.gz")
envSNPSgen <- vcfR2genlight(m3Clim)

####PCA with environmental Variables###
pca <- glPca(envSNPSgen, nf=2)

BuboPCA.scores<-as.data.frame(pca$scores)

ggplot(BuboPCA.scores, aes(x=PC1, y=PC2)) +
  geom_point(cex = 4, alpha=.75)+
  theme_classic()+
  theme(legend.position = "bottom", legend.key.size = unit(0, 'cm'),
        legend.title = element_text(size=8), legend.text = element_text(size=6.5))

ggplot(BuboPCA.scores, aes(x=PC1, y=PC3)) +
  geom_point(cex = 4, alpha=.75)+
  theme_classic()+
  theme(legend.position = "bottom", legend.key.size = unit(0, 'cm'),
        legend.title = element_text(size=8), legend.text = element_text(size=6.5))

ggplot(BuboPCA.scores, aes(x=PC2, y=PC3)) +
  geom_point(cex = 4, alpha=.75)+
  theme_classic()+
  theme(legend.position = "bottom", legend.key.size = unit(0, 'cm'),
        legend.title = element_text(size=8), legend.text = element_text(size=6.5))


barplot(100*RADBuboPCA2$eig/sum(RADBuboPCA2$eig), col = heat.colors(56), main="PCA Eigenvalues", ylab="Percent of variance explained", xlab="Eigenvalues")
