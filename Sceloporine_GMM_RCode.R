### Geometric Morphometrics of mAXilla ###

# Read in landmark tps file
library(curl)
library(geomorph)
f <- curl("https://raw.githubusercontent.com/rangerjohn99/Sceloporine-GMM-Repo/main/Maxilla%20Lateral/TPS%20Files/Sceloporine%20Maxilla%20Landmarks%20No%20Missing%20LM.tps")
raw_data <- readland.tps(f, specID = c("imageID")) # the function "readland.tps" reads the landmark data in tps format and returns a 3D array of the coordinate data
plot(raw_data)
head(raw_data)

# Read in csv file of specimens
f2 <- curl("https://raw.githubusercontent.com/rangerjohn99/Sceloporine-GMM-Repo/main/Maxilla%20Lateral/TPS%20Files/Sceloporine%20Maxilla%20Specimens%20No%20Missing%20LM.CSV")
specimenList <- read.csv(f2, header = FALSE, sep = ",", stringsAsFactors = TRUE) # this is a matrix of each specimen  
head(specimenList)
names(specimenList)[2] <- "Specimen"
specimen <- gsub("\\\\", "", specimenList$Specimen) # backslashes are terrible to work with
specimen <- gsub("'C:UsersKempLabBoxKemp LabGeometric Morphometrics ProjectProcessed ImagesMaxillaMaxilla Lateral", "", specimen)
library(dplyr)
genus <-gsub("_.*","",specimen) # make a separate genus vector
speciesV1 <-gsub("_M-.*","",specimen) # make a separate species vector part 1
species <-gsub("_CJB.*","",speciesV1) # make a separate species vector part 2



### GENERALIZED PROCRUSTES ANALYSIS ### alligns all the landmarks of all specimens

GPA_landmarks <- gpagen(raw_data) # performs Generalized Procrustes analysis of landmarks and creates aligned Procrustes coordinates
plot(GPA_landmarks)


## CREATE GMM DATAFRAMES

GMM_data <-geomorph.data.frame(coords=GPA_landmarks$coords,
                               size=GPA_landmarks$Csize, species=species, genus=genus, specimen = specimen)


## PRINCIPAL COMPONENT ANALYSIS ##

GMM_data$coords <- two.d.array(GMM_data$coords) #get the data in XY format for PCA

Sceloporine_PCA <- prcomp(GMM_data$coords) #PC analysis



# PLOT PCA #

PC_scores <- as.data.frame(Sceloporine_PCA$x)
PC_scores <- cbind(PC_scores, genus= GMM_data$genus)
percentage <- round(Sceloporine_PCA$sdev / sum(Sceloporine_PCA$sdev) * 100, 2) # find percentage variance explained by PC's
percentage <- paste( colnames(PC_scores), "(", paste( as.character(percentage), "%", ")", sep="") )

library(ggplot2)
library(ggforce)
p<-ggplot(PC_scores,aes(x=PC1,y=PC2,color=genus)) + 
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=species), size = 1) +
  geom_point(size =3)+ xlab(percentage[1]) + ylab(percentage[2]) +
  theme_classic()
p

## Discriminant Function Analysis ##

library(Morpho)

DFA<-CVA(GMM_data$coords, GMM_data$genus, cv = FALSE) # performs CVA (canonical variation analysis) aka Discriminant function analysis

barplot(DFA$Var[,2]) # Variance explained by the canonical roots

# Plot first two DF axes #

DFA_cva <- data.frame(DFA$CVscores, genus = DFA$groups)

ggplot(DFA_cva, aes(CV.1, CV.2)) +
  geom_point(aes(color = genus)) + theme_classic()

# alternative plot species

plot(DFA$CVscores, col=as.numeric(GMM_data$speices), pch=as.numeric(GMM_data$species), typ="n",asp=1, 
     xlab=paste("1st canonical axis", paste(round(DFA$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(DFA$Var[2,2],1),"%"))) 
text(DFA$CVscores, as.character(GMM_data$species), col=as.numeric(GMM_data$species), cex=.7)

# alternative plot genus

plot(DFA$CVscores, col=as.numeric(GMM_data$genus), pch=as.numeric(GMM_data$genus), typ="n",asp=1, 
     xlab=paste("1st canonical axis", paste(round(DFA$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(DFA$Var[2,2],1),"%"))) 
text(DFA$CVscores, as.character(GMM_data$genus), col=as.numeric(GMM_data$genus), cex=.7)

# Plot Mahalanobis distances as dendrogram #

dendroS=hclust(DFA$Dist$GroupdistMaha)
dendroS$labels=levels(GMM_data$genus)
par(mar=c(6.5,4.5,1,1))
dendroS=as.dendrogram(dendroS)
plot(dendroS, main='',sub='', xlab="",
     ylab='Mahalanobis distance')

## Missing landmarks dataset ##

# Read in landmark tps file

library(curl)
library(geomorph)
f3 <- curl("https://raw.githubusercontent.com/rangerjohn99/Sceloporine-GMM-Repo/main/Maxilla%20Lateral/TPS%20Files/Sceloporine%20Maxilla%20Landmarks.TPS")
raw_data2 <- readland.tps(f3, specID = c("imageID"), negNA = TRUE) # the function "readland.tps" reads the landmark data in tps format and returns a 3D array of the coordinate data
plot(raw_data2)
head(raw_data2)

# Read in csv file of specimens

f4 <- curl("https://raw.githubusercontent.com/rangerjohn99/Sceloporine-GMM-Repo/main/Maxilla%20Lateral/TPS%20Files/Sceloporine%20Maxilla%20Specimens%20Correct%20N%20Updated.CSV")
specimenList2 <- read.csv(f4, header = FALSE, sep = ",", stringsAsFactors = TRUE) # this is a matrix of each specimen  
head(specimenList2)
names(specimenList2)[2] <- "Specimen2"
specimen2 <- gsub("\\\\", "", specimenList2$Specimen2) # backslashes are terrible to work with
specimen2 <- gsub("'C:UsersKempLabBoxKemp LabGeometric Morphometrics ProjectProcessed ImagesMaxillaMaxilla Lateral", "", specimen2)
library(dplyr)
genus2 <-gsub("_.*","", specimen2) # make a separate genus vector
speciesV2 <-gsub("_M-.*","",specimen2) # make a separate species vector part 1
species2 <-gsub("_CJB.*","",speciesV2) # make a separate species vector part 2

### GENERALIZED PROCRUSTES ANALYSIS ### alligns all the landmarks of all specimens

estimate.missing(raw_data2, method = c("TPS", "Reg")) # need to estimate missing values before GPA       
GPA_landmarks2 <- gpagen(raw_data2) # performs Generalized Procrustes analysis of landmarks and creates aligned Procrustes coordinates
plot(GPA_landmarks2)

## CREATE GMM DATAFRAMES

GMM_data2 <-geomorph.data.frame(coords=GPA_landmarks2$coords,
                               size=GPA_landmarks2$Csize, species=species2, genus=genus2, specimen = specimen2)

## PRINCIPAL COMPONENT ANALYSIS ##

GMM_data2$coords <- two.d.array(GMM_data2$coords) #get the data in XY format for PCA

Sceloporine_PCA2 <- prcomp(GMM_data2$coords) #PC analysis

# PLOT PCA #

PC_scores2 <- as.data.frame(Sceloporine_PCA2$x)
PC_scores2 <- cbind(PC_scores2, genus= GMM_data2$genus)
percentage2 <- round(Sceloporine_PCA2$sdev / sum(Sceloporine_PCA2$sdev) * 100, 2) # find percentage variance explained by PC's
percentage2 <- paste(colnames(PC_scores2), "(", paste( as.character(percentage2), "%", ")", sep="") )

library(ggplot2)
library(ggforce)
p2<-ggplot(PC_scores2,aes(x=PC1,y=PC2,color=genus)) + 
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=species), size = 1) +
  geom_point(size =3)+ xlab(percentage[1]) + ylab(percentage[2]) +
  theme_classic()
p2

#+ find which landmarks are important in discriminating between groups (genera) (David will get more code for this)

#+ estimate location of missing landmarks (David will get more code for this)

#+ predict group assignment based on DFA