### Geometric Morphometrics of Maxilla ###

# Read in landmark tps file__________________________________________________________________________ 
library(curl)
library(geomorph)
f <- curl("https://raw.githubusercontent.com/rangerjohn99/Sceloporine-GMM-Repo/main/Maxilla%20Lateral/TPS%20Files/Sceloporine%20Maxilla%20Landmarks%20No%20Missing%20LM.tps")
raw_data <- readland.tps(f, specID = c("imageID")) # the function "readland.tps" reads the landmark data in tps format and returns a 3D array of the coordinate data
plot(raw_data)
head(raw_data)

# Read in csv file of specimens__________________________________________________________________________ 
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



### GENERALIZED PROCRUSTES ANALYSIS ### alligns all the landmarks of all specimens__________________________________________________________________________ 

GPA_landmarks <- gpagen(raw_data) # performs Generalized Procrustes analysis of landmarks and creates aligned Procrustes coordinates
plot(GPA_landmarks)



## CREATE GMM DATAFRAMES__________________________________________________________________________ 

GMM_data <-geomorph.data.frame(coords=GPA_landmarks$coords,
                               size=GPA_landmarks$Csize, species=species, genus=genus, specimen = specimen)


## PRINCIPAL COMPONENT ANALYSIS ##__________________________________________________________________________ 

GMM_data$coords <- two.d.array(GMM_data$coords) #get the data in XY format for PCA

Sceloporine_PCA <- prcomp(GMM_data$coords) #PC analysis

Sceloporine_PCA2 <- gm.prcomp(GPA_landmarks$coords)
summary(Sceloporine_PCA2)
plot(Sceloporine_PCA2, main = "PCA", col=as.numeric(as.factor(genus)), cex = 1.5, cex.lab = 1.5, font.lab = 2)
text(Sceloporine_PCA2$x, as.character(as.factor(genus)),col=as.numeric(as.factor(genus)), cex=.7)

## VIEW SHAPES ON PRINCIPAL COMPONENTS ##__________________________________________________________________________ 
scallinks <- matrix(c(1,rep(2:11, each=2),1), nrow=11, byrow=TRUE)
plotRefToTarget(Sceloporine_PCA2$shapes$shapes.comp1$min, Sceloporine_PCA2$shapes$shapes.comp1$max, method = "points", mag = 1, links = scallinks)
plotRefToTarget(Sceloporine_PCA2$shapes$shapes.comp2$min, Sceloporine_PCA2$shapes$shapes.comp2$max, method = "vector", mag = 1, links = scallinks)



# PLOT PCA #__________________________________________________________________________ 

PC_scores <- as.data.frame(Sceloporine_PCA$x)
PC_scores <- cbind(PC_scores, genus= GMM_data$genus, est = as.factor("FALSE"))
percentage <- round(Sceloporine_PCA$sdev^2 / sum(Sceloporine_PCA$sdev^2) * 100, 2) # find percentage variance explained by PC's
percentage <- paste( colnames(PC_scores), "(", paste( as.character(percentage), "%", ")", sep="") )

library(ggplot2)
library(ggforce)
p<-ggplot(PC_scores,aes(x=PC1,y=PC2,color=genus)) + 
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=genus), size = 1) +
  geom_point(size =3)+ xlab(percentage[1]) + ylab(percentage[2]) +
  theme_classic()
p

plot3d(PC_scores[,1:3])#interactive 3D plots
text3d(PC_scores[,1:3],texts=(PC_scores$genus),pos=4,cex=.5)

# Calculating Mean Shapes for Genera  
Sceloporus_mshape <- mshape(GPA_landmarks$coords[,,1:17])
Urosaurus_mshape <- mshape(GPA_landmarks$coords[,,18:22])
Uta_mshape <- mshape(GPA_landmarks$coords[,,23:34])


#### Canonical Variate Analysis ####__________________________________________________________________________ 

library(Morpho)

cva.1<-CVA(PC_scores[,1:14], groups = genus, rounds = 10000)
DFA_cva <- data.frame(cva.1$CVscores, genus = genus)
ggplot(DFA_cva, aes(CV.1, CV.2)) +
  geom_point(aes(color = genus)) + xlab(paste("1st Canonical Axis", paste(round(cva.1$Var[1,2],1),"%"))) + ylab(paste("2nd Canonical Axis", paste(round(cva.1$Var[2,2],1),"%"))) + 
  theme_classic()

# Mahalanobis distance probabilities
cva.1$Dist$probsMaha
?showPC
# Visualize a shape change from score -5 to 5 for CV1/CV2______________________________________________________________
cvall <- cva.1
proc<-procSym(GPA_landmarks$coords)
cvvis5 <- 5*cvall$CVvis[,1]+cvall$Grandm
cvvisNeg5 <- -5*cvall$CVvis[,1]+cvall$Grandm
cvvis5 <- showPC(cvvis5,proc$PCs[,1:14],proc$mshape)
cvvisNeg5 <- showPC(cvvisNeg5,proc$PCs[,1:14],proc$mshape)
deformGrid2d(cvvis5,cvvisNeg5,ngrid = 0, wireframe = 1:11)
CVA1p<-cvvis5
CVA1n<-cvvisNeg5
cvvis5_2 <- 5*cvall$CVvis[,2]+cvall$Grandm
cvvisNeg5_2 <- -2.5*cvall$CVvis[,2]+cvall$Grandm
cvvis5_2 <- showPC(cvvis5_2,proc$PCs[,1:14],proc$mshape)
cvvisNeg5_2 <- showPC(cvvisNeg5_2,proc$PCs[,1:14],proc$mshape)
deformGrid2d(cvvis5_2,cvvisNeg5_2,ngrid = 0, wireframe = 1:11,col1 = 1, col2 = 2)
CVA2p<-cvvis5_2
CVA2n<-cvvisNeg5_2

### RANDOM FOREST CLASSIFICATION ###:Non-parametric
PC_scores_rf <- data.frame(Sceloporine_PCA$x, genus= as.factor(genus))
library(randomForest)
set.seed(123)
Scelop.rf <- randomForest(genus ~., data=PC_scores_rf)
print(Scelop.rf)
rf_acc <- Scelop.rf$confusion
rf_acc <- 1-rf_acc[,4] # percent correct classification
rf_acc


## Missing landmarks dataset ##______________________________________________________________

# Read in landmark tps file

library(curl)
library(geomorph)
f3 <- curl("https://raw.githubusercontent.com/rangerjohn99/Sceloporine-GMM-Repo/main/Maxilla%20Lateral/TPS%20Files/Sceloporine%20Maxilla%20Landmarks.TPS")
raw_data2 <- readland.tps(f3, specID = c("imageID"), negNA = TRUE) # the function "readland.tps" reads the landmark data in tps format and returns a 3D array of the coordinate data
plot(raw_data2)
head(raw_data2)


missing_landmarks <- apply(is.na(raw_data2), 3, which) #find which rows have missing landmarks
## Get names of elements with length > 0
specimens_missing_landmarks <- names(missing_landmarks)[lapply(missing_landmarks, length) > 0]

# Read in csv file of specimens

f4 <- curl("https://raw.githubusercontent.com/rangerjohn99/Sceloporine-GMM-Repo/main/Maxilla%20Lateral/TPS%20Files/Sceloporine%20Maxilla%20Specimens%20Correct%20N%20Updated.CSV")
specimenList2 <- read.csv(f4, header = FALSE, sep = ",", stringsAsFactors = TRUE) # this is a matrix of each specimen  
head(specimenList2)
names(specimenList2)[2] <- "Specimen2"
specimen2 <- gsub("\\\\", "", specimenList2$Specimen2) # backslashes are terrible to work with
specimen2 <- gsub("'C:UsersKempLabBoxKemp LabGeometric Morphometrics ProjectProcessed ImagesMaxillaMaxilla Lateral", "", specimen2)
specimen2 <-gsub("_maxilla.*","", specimen2) # make a separate genus vector
specimen2 <-gsub("_maxiila.*","", specimen2) # make a separate genus vector
library(dplyr)
genus2 <-gsub("_.*","", specimen2) # make a separate genus vector
speciesV2 <-gsub("_M-.*","",specimen2) # make a separate species vector part 1
species2 <-gsub("_CJB.*","",speciesV2) # make a separate species vector part 2



### ESTIMATE LOCATION OF MISSING LANDMARKS______________________________________________________________

estimated_landmarks <- estimate.missing(raw_data2, method = c("TPS")) # need to estimate missing values before GPA       


### GENERALIZED PROCRUSTES ANALYSIS ### alligns all the landmarks of all specimens______________________________________________________________

GPA_landmarks2 <- gpagen(estimated_landmarks) # performs Generalized Procrustes analysis of landmarks and creates aligned Procrustes coordinates

plot(GPA_landmarks2) # a bit messy but that's expected


### SUBSET TPS TO GENERATE DATASET INCLUDING ONLY SPECIMENS WITH ESTIMATED LANDMARKS______________________________________________________________

library(tidyverse)
new <- list(land = estimated_landmarks, size = GPA_landmarks2$Csize,  species=species2, genus = genus2, specimen = specimen2)

# I first subset the original landmark data.
# Data Processing

GPA_sub <- new

# replace the "land" array (16 x 2 x 361 dimensions) in rawData with a list of 361 16 x 2 arrays
splitLand <- list(GPA_sub)
for (i in 1:62){
  
  splitLand[[i]] <-GPA_sub$land[1:11,1:2,i]
  
}
GPA_sub[["land"]] <- splitLand

# Specify a vector of the specimen names you want to extract

missing_landmarks <- apply(is.na(raw_data2), 3, which) #find which rows have missing landmarks
  ## Get names of elements with length > 0
specimens_missing_landmarks <- names(missing_landmarks)[lapply(missing_landmarks, length) > 0]
specimens_missing_landmarks <- gsub("\\\\", "", specimens_missing_landmarks) # backslashes are terrible to work with
specimens_missing_landmarks <- gsub("C:UsersKempLabBoxKemp LabGeometric Morphometrics ProjectProcessed ImagesMaxillaMaxilla Lateral", "", specimens_missing_landmarks)
specimens_missing_landmarks <-gsub("_maxilla.*","", specimens_missing_landmarks) # make a separate genus vector
specimens_missing_landmarks <-gsub("_maxiila.*","", specimens_missing_landmarks) # make a separate genus vector

SpecimensToExtract <- c(specimens_missing_landmarks)

# Create a tibble of just specimen names from your dataset and add a row ID column

specimens <- tibble(specimenName=as.character(new$specimen))

specimens_id <- rowid_to_column(specimens, "ID")

# Filter that tibble to include only the ones you want to extract

specimens_sub <- filter(specimens_id,specimenName %in% SpecimensToExtract)

# Set up vectors of variables to pull out for each specimen

land <- vector("list",nrow(specimens_sub))

species_s <- vector()

specimen_s <- vector()

genus_s <- vector()

size_s <- vector()


# Extract from rawData just the specimens for which you want data

j <- 0
for (i in specimens_sub$ID){
  
  j <- j + 1
  
  land[[j]] <- GPA_sub[["land"]][[i]]
  
  species_s <- c(species_s,as.character(GPA_sub[["species"]][[i]]))
  
  specimen_s <- c(specimen_s,as.character(GPA_sub[["specimen"]][[i]]))
  
  size_s <- c(size_s,as.character(GPA_sub[["size"]][[i]]))
  
  genus_s <- c(genus_s,as.character(GPA_sub[["genus"]][[i]]))
  
}

# Assemble the extracted data into a format like that of your original dataset

GMM_data_s <- list("land"=unlist(land),"species"=as.factor(species_s),"specimen"=as.factor(specimen_s), "genus"=as.factor(genus_s), "size" =size_s) #raw data of fossil specimens


# the lines below recast the list in "land" to the original array format and attributes

dim(GMM_data_s$land) <- c(11,2,nrow(specimens_sub))

attributes(GMM_data_s$land)$dimnames[[3]] <- GMM_data_s$specimenName


## CREATE GMM DATAFRAMES______________________________________________________________

GMM_data2 <-geomorph.data.frame(coords=GMM_data_s$land,
                               size=GMM_data_s$size, species=GMM_data_s$species, genus=GMM_data_s$genus, specimen = GMM_data_s$specimen)

GMM_data2$coords <- two.d.array(GMM_data2$coords) #get the data in XY format for PCA


## PRINCIPAL COMPONENT ANALYSIS ##______________________________________________________________

# PROJECT ESTIMATED DATA #
dimnames(GMM_data$coords) <- NULL
Sceloporine_PCA <- prcomp(GMM_data$coords) #PC analysis

Est_PCA <- predict(Sceloporine_PCA, GMM_data2$coords)
Est_PC_scores <- as.data.frame(Est_PCA)
Est_PC_scores <- cbind(Est_PC_scores, genus=GMM_data_s$genus, est = as.factor("TRUE"))

All_PC_scores <- rbind(PC_scores, Est_PC_scores) # create a new dataframe with the original PC scores and the PC scores the estimated specimens

# PLOT PCA #
speciescolors <- c("#666600", "#C9D42D" ,"#42CEBC" ,"#F2AB1F" ,"#864ED0" ,"#261ACE", "#086AFD", "#08FD6A", "#0C8B3F", "#E50CF5", "#FF5E00","#FF0000", "#FF6A6A", "#D5930F", "#9E1616")
speciesshapes <- c(rep(16,16), rep(20,30))

library(ggplot2)
library(ggforce)
p2<-ggplot(All_PC_scores,aes(x=PC1,y=PC2,color=genus, shape = est)) + 
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=species), size = 1) +
  geom_point(size =3)+ xlab(percentage[1]) + ylab(percentage[2]) +
  #scale_color_manual(name = "Species", breaks=levels(All_PC_scores$genus),  values=c(speciescolors, "black", "black", "black", "black", "black")) +
  #scale_shape_manual(values = c(speciesshapes), guide = 'none') + 
  theme_classic()
p2


#### Linear Discriminant Analysis (LDA) ####________________________________________________________________________

# Run LDA Analysis
Scelop.lda<-lda(PC_scores[,1:14], grouping = genus,CV=F) #Performs LDA

plot(Scelop.lda, col = as.integer(genus))

# See LDA Results  
# Classifying Unknown Specimens in LDA
Test.pred<-predict(Scelop.lda, Est_PC_scores[,1:14]) #Performs Prediction
Est.posteriors<-Test.pred$posterior
Est.scores<-Test.pred$x

table(Est_PC_scores$genus, Test.pred$class)


# Classifying Unknown Specimens in CVA ####________________________________________________________________________
 
CVAProb<-typprobClass(Est_PC_scores[,1:14], PC_scores[,1:14], groups = as.factor(PC_scores$genus),small = T, method = "wilson", cova = T, outlier = 0.001, cv =T)
CVAProb$probsCV
CVAProb$groupaffin
CVAProb$probs

table(Est_PC_scores$genus, CVAProb$groupaffin)


#+ find which landmarks are important in discriminating between groups (genera) (David will get more code for this)


#+ predict group assignment based on KNN & RF?






