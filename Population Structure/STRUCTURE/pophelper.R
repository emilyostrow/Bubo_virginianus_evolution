# install dependencies
#install.packages(c("ggplot2","gridExtra","label.switching","tidyr","remotes"),repos="https://cloud.r-project.org")
library(ggplot2)
library(gridExtra)
# install pophelper package from GitHub
#remotes::install_github('royfrancis/pophelper')
# load library for use
library(pophelper)
##following http://www.royfrancis.com/pophelper/articles/index.html

# STRUCTURE files (do not use this command to read local files)
sfiles <- c("structureOutk1v1.txt", "structureOutk1v2.txt", "structureOutk1v3.txt", "structureOutk2v1.txt", "structureOutk2v2.txt", "structureOutk2v3.txt", "structureOutk3v1.txt", "structureOutk3v2.txt", "structureOutk3v3.txt", "structureOutk4v1.txt", "structureOutk4v2.txt", "structureOutk4v3.txt", "structureOutk5v1.txt","structureOutk5v2.txt","structureOutk5v3.txt")

# basic usage
slist <- readQ(files=sfiles)
readQ(files=sfiles,filetype="structure")
class(slist)
# view head of first converted file
head(slist[[1]])
# qlist attributes
attributes(slist)
# dataframe attributes
attributes(slist[[1]])

# basic usage
tr1 <- tabulateQ(qlist=slist)
tabulateQ(tlist)
tabulateQ(alist, writetable=TRUE)

# choose files
# files=choose.files(multi=TRUE)
tabulateQ(qlist=readQ(files))

head(tabulateQ(slist))
sr1 <- summariseQ(tr1)
summariseQ(tr1, writetable=TRUE, exportpath=getwd())


####calculating best k####
evannoMethodStructure(data=sr1, exportplot=T, exportpath=getwd())
em <- evannoMethodStructure(summariseQ(tabulateQ(slist)))

####graphs####
slist1 <- alignK(slist[c(4,5,6,7,8,9, 10, 11, 12)])
p1 <- plotQ(slist1,imgoutput="join",returnplot=T,exportplot=F,basesize=11)
grid.arrange(p1$plot[[1]])

