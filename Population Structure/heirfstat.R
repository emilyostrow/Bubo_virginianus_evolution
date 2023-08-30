library(pophelper)
library(diveRsity)
library(adegenet)
library(vcfR)
library(hierfstat)

##following https://popgen.nescent.org/DifferentiationSNP.html
####setting up the data####
##read in structure populations
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sfiles <- "structureOutk3v1.txt"
slist <- readQ(files=sfiles)

##read in points and vcf
points <- read.csv("StructureOrderPoints.csv", stringsAsFactors = F, header = T)
allindsvcf <- read.vcfR("m3MAC2Unlinkgq30.vcf.gz")
allindsgenind <- vcfR2genind(allindsvcf)

###identify which individuals are 90% or greater in a population
structurePops <- cbind(slist$structureOutk3v1$Cluster1,slist$structureOutk3v1$Cluster2,slist$structureOutk3v1$Cluster3,rep(4, times=114))
row.names(structurePops) <- points$FileName
clust1Vec90 <- structurePops[,1]>=.9
clust2Vec90 <- structurePops[,2]>=.9
clust3Vec90 <- structurePops[,3]>=.9

clust1Names <- row.names(structurePops[clust1Vec90,])
clust2Names <- row.names(structurePops[clust2Vec90,])
clust3Names <- row.names(structurePops[clust3Vec90,])

###identifying which structure population each individual most closely aligns with
for (i in 1:length(structurePops[,1])){
  if (structurePops[i,1]>structurePops[i,2]&&structurePops[i,1]>structurePops[i,3]){
    structurePops[i,4] <- 1
  }
  else if (structurePops[i,2]>structurePops[i,1]&&structurePops[i,2]>structurePops[i,3]){
    structurePops[i,4] <- 2
  }
  else{
    structurePops[i,4] <- 3
  }
}
structurePopsOriginal <- structurePops
#write.table(cbind(locs$LongName,structurePops[,4]), "structureMajPops.txt",quote = F, sep="\t")

allindsgenind@pop <- as.factor(structurePops[,4])
allindsHierf <- genind2hierfstat(allindsgenind)
allindsFst <- basic.stats(allindsHierf)

pairwiseFst <- genet.dist(allindsHierf, method = "WC84")
pairwiseFstmat <- as.matrix(pairwiseFst)
colnames(pairwiseFstmat) <- c("east","northwest","southwest")
row.names(pairwiseFstmat) <- c("east","northwest","southwest")

pops90comb <- c(clust1Names,clust2Names,clust3Names)
for (i in 1:length(pops90comb)){
  pops90comb[i] <-  paste0("B_virginianus_",pops90comb[i])
}
pops90comb[33] <- "B_virginianus_UAMX6198_Comb"
pops90comb[4] <- "B_virginianus_475289_Comb"
pops90comb[14] <-  "B_virginianus_ENO310_Comb"

allindsvcf@gt <- allindsvcf@gt[,colnames(allindsvcf@gt) %in% c("FORMAT",pops90comb)]
ninetypgenind <- vcfR2genind(allindsvcf)

ninetypgenind@pop <- as.factor(c(rep(2,times=23),rep(1,times=16),3,3,1,3,1))
ninetypHierf <- genind2hierfstat(ninetypgenind)
ninetypindsFst <- basic.stats(ninetypHierf)


ninetyppairwiseFst <- genet.dist(ninetypHierf, method = "WC84")
ninetyppairwiseFstmat <- as.matrix(ninetyppairwiseFst)
colnames(ninetyppairwiseFstmat) <- c("east","northwest","southwest")
row.names(ninetyppairwiseFstmat) <- c("east","northwest","southwest")

pairwiseFstmat
ninetyppairwiseFstmat

distsList <- list()
####randomize populations for testing significance####
for (i in 1:100){
  structurePopsOriginal[,4] <- sample(structurePopsOriginal[,4])
  allindsgenind@pop <- as.factor(structurePopsOriginal[,4])
  allindsHierf <- genind2hierfstat(allindsgenind)
  distsList[[i]] <- genet.dist(allindsHierf, method = "WC84")
  
}

##for each population comparison, put the numbers into a vector
clust12vec <- vector()
clust13vec <- vector()
clust23vec <- vector()

for (i in 1:length(distsList)){
  if (distsList[[i]][1]<0){
    distsList[[i]][1]=0
  }
  clust12vec[i] <- distsList[[i]][1]
}
for (i in 1:length(distsList)){
  if (distsList[[i]][2]<0){
    distsList[[i]][2]=0
  }
  clust13vec[i] <- distsList[[i]][2]
}
for (i in 1:length(distsList)){
  if (distsList[[i]][3]<0){
    distsList[[i]][3]=0
  }
  clust23vec[i] <- distsList[[i]][3]
}

##do a t-test to test the significance of the population differentiation to the 'background' differentiation in GHOW
t.test(clust12vec,mu=pairwiseFst[1], alternative = "less")
t.test(clust13vec,mu=pairwiseFst[2], alternative = "less")
t.test(clust23vec,mu=pairwiseFst[3], alternative = "less")

##average background differentiation
mean(c(clust12vec,clust13vec,clust23vec))
 

#### identify fixed differences for each population ####
ninetypFst <- basic.stats(ninetypHierf)
fixedvec <- ninetypFst$perloc$Fst
names(fixedvec) <- rownames(ninetypFst$perloc)
overallFixed <- names(na.omit(fixedvec[fixedvec==1]))

##comparing 1 and 2
pops90comb <- c(clust1Names,clust2Names)
for (i in 1:length(pops90comb)){
  pops90comb[i] <-  paste0("B_virginianus_",pops90comb[i])
}
pops90comb[33] <- "B_virginianus_UAMX6198_Comb"
pops90comb[4] <- "B_virginianus_475289_Comb"
pops90comb[14] <-  "B_virginianus_ENO310_Comb"

allindsvcf <- read.vcfR("m3MAC1Unlinkgq30.vcf.gz")
allindsvcf@gt <- allindsvcf@gt[,colnames(allindsvcf@gt) %in% c("FORMAT",pops90comb)]
ninetypgenind <- vcfR2genind(allindsvcf)

ninetypgenind@pop <- as.factor(c(rep(1,times=18),rep(2,times=23)))
ninetypHierf <- genind2hierfstat(ninetypgenind)
ninetypFst <- basic.stats(ninetypHierf)

fixedvec <- ninetypFst$perloc$Fst
names(fixedvec) <- rownames(ninetypFst$perloc)
eightyp12 <- names(na.omit(fixedvec[fixedvec>=.8]))
sort(fixedvec, decreasing = T)

##comparing 1 and 3
pops90comb <- c(clust1Names,clust3Names)
for (i in 1:length(pops90comb)){
  pops90comb[i] <-  paste0("B_virginianus_",pops90comb[i])
}
pops90comb[4] <- "B_virginianus_475289_Comb"
pops90comb[14] <-  "B_virginianus_ENO310_Comb"

allindsvcf <- read.vcfR("m3MAC1Unlinkgq30.vcf.gz")
allindsvcf@gt <- allindsvcf@gt[,colnames(allindsvcf@gt) %in% c("FORMAT",pops90comb)]
ninetypgenind <- vcfR2genind(allindsvcf)

ninetypgenind@pop <- as.factor(c(rep(1,times=18),3,3,3))
ninetypHierf <- genind2hierfstat(ninetypgenind)
ninetypFst <- basic.stats(ninetypHierf)

fixedvec <- ninetypFst$perloc$Fst
names(fixedvec) <- rownames(ninetypFst$perloc)
Fixed13 <- names(na.omit(fixedvec[fixedvec==1]))
eightyp13 <- names(na.omit(fixedvec[fixedvec>=.8]))

##comparing 2 and 3
pops90comb <- c(clust2Names,clust3Names)
for (i in 1:length(pops90comb)){
  pops90comb[i] <-  paste0("B_virginianus_",pops90comb[i])
}
pops90comb[15] <- "B_virginianus_UAMX6198_Comb"

allindsvcf <- read.vcfR("m3MAC1Unlinkgq30.vcf.gz")
allindsvcf@gt <- allindsvcf@gt[,colnames(allindsvcf@gt) %in% c("FORMAT",pops90comb)]
ninetypgenind <- vcfR2genind(allindsvcf)

ninetypgenind@pop <- as.factor(c(rep(2,times=23),3,3,3))
ninetypHierf <- genind2hierfstat(ninetypgenind)
ninetypFst <- basic.stats(ninetypHierf)

fixedvec <- ninetypFst$perloc$Fst
names(fixedvec) <- rownames(ninetypFst$perloc)
Fixed23 <- names(na.omit(fixedvec[fixedvec==1]))
eightyp23 <- names(na.omit(fixedvec[fixedvec>=.8]))
write.csv(Fixed23, "fixedNorthSouth.csv")
write.csv(Fixed13, "fixedNorthEast.csv")
write.csv(eightyp12, "80pSouthEast.csv")
write.csv(eightyp13, "80pNorthEast.csv")
write.csv(eightyp23, "80pNorthSouth.csv")

