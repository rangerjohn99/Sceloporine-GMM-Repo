### Geometric Morphometrics of Dentary medial ###---------------------------------

# Read in landmark tps file---------------------------------
library(curl)
library(geomorph)
f <- curl("https://raw.githubusercontent.com/rangerjohn99/Sceloporine-GMM-Repo/main/Dentary%20Medial/TPS%20files/Sceloporine%20Dentary%20Medial%20TPS%20Images%20Semilandmarks%20With%20Landmarked%20Curves%20Overlapped%20Deleted.TPS")
raw_data <- readland.tps(f, specID = c("imageID")) # the function "readland.tps" reads the landmark data in tps format and returns a 3D array of the coordinate data
plot(raw_data)
head(raw_data)


ff <- curl("https://raw.githubusercontent.com/rangerjohn99/Sceloporine-GMM-Repo/main/Dentary%20Medial/TPS%20files/Sceloporine%20Dentary%20Medial%20TPS%20Images%20Semilandmarks%20With%20Landmarked%20Curves%20Overlapped%20Deleted%20Sliders%20File.NTS")
sliders <- readland.nts(ff)
sliders <- as.data.frame(sliders)

data.super <- gpagen(raw_data, curves = sliders, ProcD = FALSE, PrinAxes = FALSE) 
plot(data.super)

# Read in csv file of specimens---------------------------------

specimen<- dimnames(data.super$coords)[[3]]
specimen <- gsub("\\\\", "", specimen) # backslashes are terrible to work with
specimen <- gsub("C:UsersKempLabBoxKemp LabGeometric Morphometrics ProjectProcessed ImagesDentary Medial", "", specimen)
library(dplyr)
genus <-gsub("_.*","",specimen) # make a separate genus vector
speciesV1 <-gsub("_M-.*","",specimen) # make a separate species vector part 1
species <-gsub("_CJB.*","",speciesV1) # make a separate species vector part 2


## CREATE GMM DATAFRAMES---------------------------------
GMM_data <-geomorph.data.frame(coords=data.super$coords,
                               size=data.super$Csize, species=species, genus=genus, specimen = specimen)

## PRINCIPAL COMPONENT ANALYSIS ##---------------------------------
GMM_data$coords <- two.d.array(GMM_data$coords) #get the data in XY format for PCA

Sceloporine_PCA <- prcomp(GMM_data$coords) #PC analysis

Sceloporine_PCA2 <- gm.prcomp(data.super$coords)
summary(Sceloporine_PCA2)
plot(Sceloporine_PCA2, main = "PCA", col=as.numeric(as.factor(genus)), cex = 1.5, cex.lab = 1.5, font.lab = 2)
text(Sceloporine_PCA2$x, as.character(as.factor(genus)),col=as.numeric(as.factor(genus)), cex=.7)

## VIEW SHAPES ON PRINCIPAL COMPONENTS ##_---------------------------------
scallinks <- matrix(c(1,rep(2:79, each=2),1), nrow=79, byrow=TRUE) #link landmarks
plotRefToTarget(Sceloporine_PCA2$shapes$shapes.comp1$min, Sceloporine_PCA2$shapes$shapes.comp1$max, method = "points", mag = 1) #shapes corresponding to PC1 min and max; grey=min 
plotRefToTarget(Sceloporine_PCA2$shapes$shapes.comp2$min, Sceloporine_PCA2$shapes$shapes.comp2$max, method = "points", mag = 1) #shapes corresponding to PC2 min and max; grey=min

# PLOT PCA #---------------------------------

PC_scores <- as.data.frame(Sceloporine_PCA$x)
PC_scores <- cbind(PC_scores, genus= GMM_data$genus, est = as.factor("FALSE"))
percentage <- round(Sceloporine_PCA$sdev^2 / sum(Sceloporine_PCA$sdev^2) * 100, 2) # find percentage variance explained by PC's
percentage <- paste( colnames(PC_scores), "(",paste(as.character(percentage), "%", ")", sep=""))

library(ggplot2)
library(ggforce)
p<-ggplot(PC_scores,aes(x=PC1,y=PC2,color=genus)) + 
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=genus), size = 1) +
  geom_point(size =3)+ xlab(percentage[1]) + ylab(percentage[2]) +
  theme_classic()
p


# Calculating Mean Shapes for Genera  ---------------------------------
Sceloporus_mshape <- mshape(data.super$coords[,,2:32])
Urosaurus_mshape <- mshape(data.super$coords[,,33:45])
Uta_mshape <- mshape(data.super$coords[,,46:65])

plotRefToTarget(Urosaurus_mshape, Uta_mshape, method = "points", mag = 1) #shapes corresponding to PC2 min and max; grey=Urosaurus
plotRefToTarget(Sceloporus_mshape, Urosaurus_mshape, method = "points", mag = 1) #shapes corresponding to PC2 min and max; grey=Sceloporus
plotRefToTarget(Sceloporus_mshape, Uta_mshape, method = "points", mag = 1) #shapes corresponding to PC2 min and max; grey=Sceloporus

#### Canonical Variate Analysis ####---------------------------------

library(Morpho)

cva.1<-CVA(PC_scores[,1:14], groups = genus, cv = TRUE)
DFA_cva <- data.frame(cva.1$CVscores, genus = genus)
ggplot(DFA_cva, aes(CV.1, CV.2)) +
  geom_point(aes(color = genus)) + xlab(paste("1st Canonical Axis", paste(round(cva.1$Var[1,2],1),"%"))) + ylab(paste("2nd Canonical Axis", paste(round(cva.1$Var[2,2],1),"%"))) + 
  theme_classic()

#typprobs <- typprobClass(cva.1$CVscores,groups= as.factor( genus)) #classification accuracy from jackknife Crossvalidation, < 0.01 posterior prob assigned to no class
#print(typprobs)

### RANDOM FOREST CLASSIFICATION ###:Non-parametric---------------------------------
PC_scores_rf <- data.frame(Sceloporine_PCA$x[,1:5], genus= as.factor(genus))
library(randomForest)
set.seed(123)
Scelop.rf <- randomForest(genus ~., data=PC_scores_rf)
print(Scelop.rf)
rf_acc <- Scelop.rf$confusion
rf_acc <- 1-rf_acc[,5] # percent correct classification
rf_acc



