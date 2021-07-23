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
specimen <- gsub("\\\\", "", specimenList$Specimen)
specimen <- gsub("'C:UsersKempLabBoxKemp LabGeometric Morphometrics ProjectProcessed ImagesMaxillaMaxilla Lateral", "", specimen)
library(dplyr)
genus <-gsub("_.*","",specimen) #make a separate genus vector
speciesV1 <-gsub("_M-.*","",specimen) #make a separate species vector part 1
species <-gsub("_CJB.*","",speciesV1) #make a separate species vector part 2

