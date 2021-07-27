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







