#devtools::install_github("DevonDeRaad/SNPfiltR")
library(SNPfiltR)
library(vcfR)
library(adegenet)
library(tsne)
library(ggplot2)
library(PopGenome)

#### Assessing filtering based on completeness and MAC filtering ####
n <- read.vcfR("m3M2N2.vcf")


inds <- colnames(n@gt)[2:119]
reps <- rep(1, times=118)
pops <- data.frame(inds, reps)
names(pops)[1] <- "id"
names(pops)[2] <- "pop"


#retry
y<-missing_by_sample(n, pops, cutoff=0.9)
z<-missing_by_snp(n)


#filter by SNP
filtered95<-missing_by_snp(y, cutoff = .95)
filtered90<-missing_by_snp(y, cutoff = .9)
filtered85<-missing_by_snp(y, cutoff = .85)
filtered80<-missing_by_snp(y, cutoff = .8)
filtered75<-missing_by_snp(y, cutoff = .75)


filteredgq95 <- hard_filter(filtered95,gq=30)
filteredgq90 <- hard_filter(filtered90,gq=30)
filteredgq85 <- hard_filter(filtered85,gq=30)
filteredgq80 <- hard_filter(filtered80,gq=30)
filteredgq75 <- hard_filter(filtered75,gq=30)

filteredgqsingletons95 <- min_mac(filteredgq95,2)
filteredgqsingletons90 <- min_mac(filteredgq90,2)
filteredgqsingletons85 <- min_mac(filteredgq85,2)
filteredgqsingletons80 <- min_mac(filteredgq80,2)
filteredgqsingletons75 <- min_mac(filteredgq75,2)

##unlinked SNPs
#filter to one SNP per locus
filteredgq95@gt<-filteredgq95@gt[!duplicated(filteredgq95@fix[,1]),]
filteredgq95@fix<-filteredgq95@fix[!duplicated(filteredgq95@fix[,1]),]

filteredgq90@gt<-filteredgq90@gt[!duplicated(filteredgq90@fix[,1]),]
filteredgq90@fix<-filteredgq90@fix[!duplicated(filteredgq90@fix[,1]),]

filteredgq85@gt<-filteredgq85@gt[!duplicated(filteredgq85@fix[,1]),]
filteredgq85@fix<-filteredgq85@fix[!duplicated(filteredgq85@fix[,1]),]

filteredgq80@gt<-filteredgq80@gt[!duplicated(filteredgq80@fix[,1]),]
filteredgq80@fix<-filteredgq80@fix[!duplicated(filteredgq80@fix[,1]),]

filteredgq75@gt<-filteredgq75@gt[!duplicated(filteredgq75@fix[,1]),]
filteredgq75@fix<-filteredgq75@fix[!duplicated(filteredgq75@fix[,1]),]

filteredgqsingletons95@gt<-filteredgqsingletons95@gt[!duplicated(filteredgqsingletons95@fix[,1]),]
filteredgqsingletons95@fix<-filteredgqsingletons95@fix[!duplicated(filteredgqsingletons95@fix[,1]),]

filteredgqsingletons90@gt<-filteredgqsingletons90@gt[!duplicated(filteredgqsingletons90@fix[,1]),]
filteredgqsingletons90@fix<-filteredgqsingletons90@fix[!duplicated(filteredgqsingletons90@fix[,1]),]

filteredgqsingletons85@gt<-filteredgqsingletons85@gt[!duplicated(filteredgqsingletons85@fix[,1]),]
filteredgqsingletons85@fix<-filteredgqsingletons85@fix[!duplicated(filteredgqsingletons85@fix[,1]),]

filteredgqsingletons80@gt<-filteredgqsingletons80@gt[!duplicated(filteredgqsingletons80@fix[,1]),]
filteredgqsingletons80@fix<-filteredgqsingletons80@fix[!duplicated(filteredgqsingletons80@fix[,1]),]

filteredgqsingletons75@gt<-filteredgqsingletons75@gt[!duplicated(filteredgqsingletons75@fix[,1]),]
filteredgqsingletons75@fix<-filteredgqsingletons75@fix[!duplicated(filteredgqsingletons75@fix[,1]),]


finalSNPs <- data.frame(matrix(NA, ncol=5, nrow=2))
colnames(finalSNPs) <- c("95 Percent","90 Percent","85 Percent","80 Percent","75 Percent")
row.names(finalSNPs) <- c("No MAC Filter", "MAC2 Filter")

finalSNPs[1,1] <- nrow(filteredgq95@gt)
finalSNPs[1,2] <- nrow(filteredgq90@gt)
finalSNPs[1,3] <- nrow(filteredgq85@gt)
finalSNPs[1,4] <- nrow(filteredgq80@gt)
finalSNPs[1,5] <- nrow(filteredgq75@gt)
finalSNPs[2,1] <- nrow(filteredgqsingletons95@gt)
finalSNPs[2,2] <- nrow(filteredgqsingletons90@gt)
finalSNPs[2,3] <- nrow(filteredgqsingletons85@gt)
finalSNPs[2,4] <- nrow(filteredgqsingletons80@gt)
finalSNPs[2,5] <- nrow(filteredgqsingletons75@gt)

finalSNPs

write.vcf(filteredgq95,"filteredgq95.vcf.gz")
write.vcf(filteredgq90,"filteredgq90.vcf.gz")
write.vcf(filteredgq85,"filteredgq85.vcf.gz")
write.vcf(filteredgq80,"filteredgq80.vcf.gz")
write.vcf(filteredgq75,"filteredgq75.vcf.gz")
write.vcf(filteredgqsingletons95,"filteredgqsingletons95.vcf.gz")
write.vcf(filteredgqsingletons90,"filteredgqsingletons90.vcf.gz")
write.vcf(filteredgqsingletons85,"filteredgqsingletons85.vcf.gz")
write.vcf(filteredgqsingletons80,"filteredgqsingletons80.vcf.gz")
write.vcf(filteredgqsingletons75,"filteredgqsingletons75.vcf.gz")

#### Assessing Popgen Stats Using Filtered Datasets ####

##unzip exported vcfs and move to appropriate folders
coords<- read.csv("StructureOrderPoints.csv" )
pops <- read.csv("structurePopAssignments.csv")

filt95 <- readData("./95/", format="VCF",big.data = T, include.unknown = T)
filt90 <- readData("./90/", format="VCF",big.data = T, include.unknown = T)
filt85 <- readData("./85/", format="VCF",big.data = T, include.unknown = T)
filt80 <- readData("./80/", format="VCF",big.data = T, include.unknown = T)
filt75 <- readData("./75/", format="VCF",big.data = T, include.unknown = T)
filt95MAC2 <- readData("./95MAC2/", format="VCF",big.data = T, include.unknown = T)
filt90MAC2 <- readData("./90MAC2/", format="VCF",big.data = T, include.unknown = T)
filt85MAC2 <- readData("./85MAC2/", format="VCF",big.data = T, include.unknown = T)
filt80MAC2 <- readData("./80MAC2/", format="VCF",big.data = T, include.unknown = T)
filt75MAC2 <- readData("./75MAC2/", format="VCF",big.data = T, include.unknown = T)

SW <- pops[pops$popID==3,3]
NW <- pops[pops$popID==2,3]
east <- pops[pops$popID==1,3]

filt95  <- set.populations(filt95, list(SW,NW,east))
filt90  <- set.populations(filt90, list(SW,NW,east))
filt85  <- set.populations(filt85, list(SW,NW,east))
filt80  <- set.populations(filt80, list(SW,NW,east))
filt75  <- set.populations(filt75, list(SW,NW,east))
filt95MAC2  <- set.populations(filt95MAC2, list(SW,NW,east))
filt90MAC2  <- set.populations(filt90MAC2, list(SW,NW,east))
filt85MAC2  <- set.populations(filt85MAC2, list(SW,NW,east))
filt80MAC2  <- set.populations(filt80MAC2, list(SW,NW,east))
filt75MAC2  <- set.populations(filt75MAC2, list(SW,NW,east))



filt95 <- F_ST.stats(filt95) 
filt95 <- diversity.stats(filt95)
filt95 <- linkage.stats(filt95)
filt95 <- neutrality.stats(filt95)
filt90 <- F_ST.stats(filt90) 
filt90 <- diversity.stats(filt90)
filt90 <- linkage.stats(filt90)
filt90 <- neutrality.stats(filt90)
filt85 <- F_ST.stats(filt85)
filt85 <- diversity.stats(filt85)
filt85 <- linkage.stats(filt85)
filt85 <- neutrality.stats(filt85)
filt80 <- F_ST.stats(filt80)  
filt80 <- diversity.stats(filt80)
filt80 <- linkage.stats(filt80)
filt80 <- neutrality.stats(filt80)
filt75 <- F_ST.stats(filt75) 
filt75 <- diversity.stats(filt75)
filt75 <- linkage.stats(filt75)
filt75 <- neutrality.stats(filt75)
filt95MAC2 <- F_ST.stats(filt95MAC2) 
filt95MAC2 <- diversity.stats(filt95MAC2)
filt95MAC2 <- linkage.stats(filt95MAC2)
filt95MAC2 <- neutrality.stats(filt95MAC2)
filt90MAC2 <- F_ST.stats(filt90MAC2) 
filt90MAC2 <- diversity.stats(filt90MAC2)
filt90MAC2 <- linkage.stats(filt90MAC2)
filt90MAC2 <- neutrality.stats(filt90MAC2)
filt85MAC2 <- F_ST.stats(filt85MAC2) 
filt85MAC2 <- diversity.stats(filt85MAC2)
filt85MAC2 <- linkage.stats(filt85MAC2)
filt85MAC2 <- neutrality.stats(filt85MAC2)
filt80MAC2 <- F_ST.stats(filt80MAC2) 
filt80MAC2 <- diversity.stats(filt80MAC2)
filt80MAC2 <- linkage.stats(filt80MAC2)
filt80MAC2 <- neutrality.stats(filt80MAC2)
filt75MAC2 <- F_ST.stats(filt75MAC2) 
filt75MAC2 <- diversity.stats(filt75MAC2)
filt75MAC2 <- linkage.stats(filt75MAC2)
filt75MAC2 <- neutrality.stats(filt75MAC2)


#Pi
statTable <- data.frame(matrix(NA, ncol=13, nrow=10))
colnames(statTable) <- c("SNPs","PiEast","PiNW","PiSW","WTEast","WTNW","WTSW","WallsBEast","WallsBNW","WallsBSW","TajDEast","TajDNW","TajDSW")
row.names(statTable) <- c("95 Percent","90 Percent","85 Percent","80 Percent","75 Percent","95 PercentMAC2","90 PercentMAC2","85 PercentMAC2","80 PercentMAC2","75 PercentMAC2")

statTable[1,1] <- nrow(filteredgq95@gt)
statTable[2,1] <- nrow(filteredgq90@gt)
statTable[3,1] <- nrow(filteredgq85@gt)
statTable[4,1] <- nrow(filteredgq80@gt)
statTable[5,1] <- nrow(filteredgq75@gt)
statTable[6,1] <- nrow(filteredgqsingletons95@gt)
statTable[7,1] <- nrow(filteredgqsingletons90@gt)
statTable[8,1] <- nrow(filteredgqsingletons85@gt)
statTable[9,1] <- nrow(filteredgqsingletons80@gt)
statTable[10,1] <- nrow(filteredgqsingletons75@gt)

statTable[1,2] <- filt95@nuc.diversity.within[1]/statTable[1,1]
statTable[2,2] <- filt90@nuc.diversity.within[1]/statTable[2,1]
statTable[3,2] <- filt85@nuc.diversity.within[1]/statTable[3,1]
statTable[4,2] <- filt80@nuc.diversity.within[1]/statTable[4,1]
statTable[5,2] <- filt75@nuc.diversity.within[1]/statTable[5,1]
statTable[6,2] <- filt95MAC2@nuc.diversity.within[1]/statTable[6,1]
statTable[7,2] <- filt90MAC2@nuc.diversity.within[1]/statTable[7,1]
statTable[8,2] <- filt80MAC2@nuc.diversity.within[1]/statTable[8,1]
statTable[9,2] <- filt85MAC2@nuc.diversity.within[1]/statTable[9,1]
statTable[10,2] <- filt75MAC2@nuc.diversity.within[1]/statTable[10,1]

statTable[1,3] <- filt95@nuc.diversity.within[2]/statTable[1,1]
statTable[2,3] <- filt90@nuc.diversity.within[2]/statTable[2,1]
statTable[3,3] <- filt85@nuc.diversity.within[2]/statTable[3,1]
statTable[4,3] <- filt80@nuc.diversity.within[2]/statTable[4,1]
statTable[5,3] <- filt75@nuc.diversity.within[2]/statTable[5,1]
statTable[6,3] <- filt95MAC2@nuc.diversity.within[2]/statTable[6,1]
statTable[7,3] <- filt90MAC2@nuc.diversity.within[2]/statTable[7,1]
statTable[8,3] <- filt80MAC2@nuc.diversity.within[2]/statTable[8,1]
statTable[9,3] <- filt85MAC2@nuc.diversity.within[2]/statTable[9,1]
statTable[10,3] <- filt75MAC2@nuc.diversity.within[2]/statTable[10,1]

statTable[1,4] <- filt95@nuc.diversity.within[3]/statTable[1,1]
statTable[2,4] <- filt90@nuc.diversity.within[3]/statTable[2,1]
statTable[3,4] <- filt85@nuc.diversity.within[3]/statTable[3,1]
statTable[4,4] <- filt80@nuc.diversity.within[3]/statTable[4,1]
statTable[5,4] <- filt75@nuc.diversity.within[3]/statTable[5,1]
statTable[6,4] <- filt95MAC2@nuc.diversity.within[3]/statTable[6,1]
statTable[7,4] <- filt90MAC2@nuc.diversity.within[3]/statTable[7,1]
statTable[8,4] <- filt80MAC2@nuc.diversity.within[3]/statTable[8,1]
statTable[9,4] <- filt85MAC2@nuc.diversity.within[3]/statTable[9,1]
statTable[10,4] <- filt75MAC2@nuc.diversity.within[3]/statTable[10,1]

#Wattersons Theta
statTable[1,5] <- filt95@theta_Watterson[1]/statTable[1,1]
statTable[2,5] <- filt90@theta_Watterson[1]/statTable[2,1]
statTable[3,5] <- filt85@theta_Watterson[1]/statTable[3,1]
statTable[4,5] <- filt80@theta_Watterson[1]/statTable[4,1]
statTable[5,5] <- filt75@theta_Watterson[1]/statTable[5,1]
statTable[6,5] <- filt95MAC2@theta_Watterson[1]/statTable[6,1]
statTable[7,5] <- filt90MAC2@theta_Watterson[1]/statTable[7,1]
statTable[8,5] <- filt80MAC2@theta_Watterson[1]/statTable[8,1]
statTable[9,5] <- filt85MAC2@theta_Watterson[1]/statTable[9,1]
statTable[10,5] <- filt75MAC2@theta_Watterson[1]/statTable[10,1]

statTable[1,6] <- filt95@theta_Watterson[2]/statTable[1,1]
statTable[2,6] <- filt90@theta_Watterson[2]/statTable[2,1]
statTable[3,6] <- filt85@theta_Watterson[2]/statTable[3,1]
statTable[4,6] <- filt80@theta_Watterson[2]/statTable[4,1]
statTable[5,6] <- filt75@theta_Watterson[2]/statTable[5,1]
statTable[6,6] <- filt95MAC2@theta_Watterson[2]/statTable[6,1]
statTable[7,6] <- filt90MAC2@theta_Watterson[2]/statTable[7,1]
statTable[8,6] <- filt80MAC2@theta_Watterson[2]/statTable[8,1]
statTable[9,6] <- filt85MAC2@theta_Watterson[2]/statTable[9,1]
statTable[10,6] <- filt75MAC2@theta_Watterson[2]/statTable[10,1]

statTable[1,7] <- filt95@theta_Watterson[3]/statTable[1,1]
statTable[2,7] <- filt90@theta_Watterson[3]/statTable[2,1]
statTable[3,7] <- filt85@theta_Watterson[3]/statTable[3,1]
statTable[4,7] <- filt80@theta_Watterson[3]/statTable[4,1]
statTable[5,7] <- filt75@theta_Watterson[3]/statTable[5,1]
statTable[6,7] <- filt95MAC2@theta_Watterson[3]/statTable[6,1]
statTable[7,7] <- filt90MAC2@theta_Watterson[3]/statTable[7,1]
statTable[8,7] <- filt80MAC2@theta_Watterson[3]/statTable[8,1]
statTable[9,7] <- filt85MAC2@theta_Watterson[3]/statTable[9,1]
statTable[10,7] <- filt75MAC2@theta_Watterson[3]/statTable[10,1]

#WallsB
statTable[1,8] <- filt95@Wall.B[1]
statTable[2,8] <- filt90@Wall.B[1]
statTable[3,8] <- filt85@Wall.B[1]
statTable[4,8] <- filt80@Wall.B[1]
statTable[5,8] <- filt75@Wall.B[1]
statTable[6,8] <- filt95MAC2@Wall.B[1]
statTable[7,8] <- filt90MAC2@Wall.B[1]
statTable[8,8] <- filt80MAC2@Wall.B[1]
statTable[9,8] <- filt85MAC2@Wall.B[1]
statTable[10,8] <- filt75MAC2@Wall.B[1]

statTable[1,9] <- filt95@Wall.B[2]
statTable[2,9] <- filt90@Wall.B[2]
statTable[3,9] <- filt85@Wall.B[2]
statTable[4,9] <- filt80@Wall.B[2]
statTable[5,9] <- filt75@Wall.B[2]
statTable[6,9] <- filt95MAC2@Wall.B[2]
statTable[7,9] <- filt90MAC2@Wall.B[2]
statTable[8,9] <- filt80MAC2@Wall.B[2]
statTable[9,9] <- filt85MAC2@Wall.B[2]
statTable[10,9] <- filt75MAC2@Wall.B[2]

statTable[1,10] <- filt95@Wall.B[3]
statTable[2,10] <- filt90@Wall.B[3]
statTable[3,10] <- filt85@Wall.B[3]
statTable[4,10] <- filt80@Wall.B[3]
statTable[5,10] <- filt75@Wall.B[3]
statTable[6,10] <- filt95MAC2@Wall.B[3]
statTable[7,10] <- filt90MAC2@Wall.B[3]
statTable[8,10] <- filt80MAC2@Wall.B[3]
statTable[9,10] <- filt85MAC2@Wall.B[3]
statTable[10,10] <- filt75MAC2@Wall.B[3]

#TajimasD
statTable[1,11] <- filt95@Tajima.D[1]
statTable[2,11] <- filt90@Tajima.D[1]
statTable[3,11] <- filt85@Tajima.D[1]
statTable[4,11] <- filt80@Tajima.D[1]
statTable[5,11] <- filt75@Tajima.D[1]
statTable[6,11] <- filt95MAC2@Tajima.D[1]
statTable[7,11] <- filt90MAC2@Tajima.D[1]
statTable[8,11] <- filt80MAC2@Tajima.D[1]
statTable[9,11] <- filt85MAC2@Tajima.D[1]
statTable[10,11] <- filt75MAC2@Tajima.D[1]

statTable[1,12] <- filt95@Tajima.D[2]
statTable[2,12] <- filt90@Tajima.D[2]
statTable[3,12] <- filt85@Tajima.D[2]
statTable[4,12] <- filt80@Tajima.D[2]
statTable[5,12] <- filt75@Tajima.D[2]
statTable[6,12] <- filt95MAC2@Tajima.D[2]
statTable[7,12] <- filt90MAC2@Tajima.D[2]
statTable[8,12] <- filt80MAC2@Tajima.D[2]
statTable[9,12] <- filt85MAC2@Tajima.D[2]
statTable[10,12] <- filt75MAC2@Tajima.D[2]

statTable[1,13] <- filt95@Tajima.D[3]
statTable[2,13] <- filt90@Tajima.D[3]
statTable[3,13] <- filt85@Tajima.D[3]
statTable[4,13] <- filt80@Tajima.D[3]
statTable[5,13] <- filt75@Tajima.D[3]
statTable[6,13] <- filt95MAC2@Tajima.D[3]
statTable[7,13] <- filt90MAC2@Tajima.D[3]
statTable[8,13] <- filt80MAC2@Tajima.D[3]
statTable[9,13] <- filt85MAC2@Tajima.D[3]
statTable[10,13] <- filt75MAC2@Tajima.D[3]

write.csv(statTable, "FilteringPopgenStats.csv")
