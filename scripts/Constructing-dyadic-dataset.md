# Making dyadic dataset for inter species comparisson

Based on the workflow described here  https://github.com/nuorenarra/Analysing-dyadic-data-with-brms 



## Table of Contents

[TOC]

------



## 1. Read in the data

```R
library(phyloseq)
library(tidyverse)
library(janitor)
library(microbiome)
library(sp)
library(rgdal)


#Read in microbiome data and associated sample data in phyloseq format. These data sets contain either 70 samples of 70 individuals
micdata<-readRDS( "ps_css.rds")
metadata<- sample_data(micdata)
summarize_phyloseq(ps_css)
```



## 2. Construct matrices to create the dayadic dataset



### 2.1 Microbiome dissimilarity/distance matrices (BC, WU)

```R
#make a key for the order of sample names and their associated individual IDs.
key <- data.frame(ID=sample_data(micdata)$identifier)

#Make Bray curtis matrix from microbiome data using vegdist function embedded in phyloseq::distance
BCM<- as.matrix(phyloseq::distance(micdata, method = "bray", type = "samples"))
#Make unweighted unifrac matrix
WUM <- as.matrix(phyloseq::distance(micdata, method = "wunifrac", type = "samples"))

#Save matrix to ready matrices folder
saveRDS(BCM,"BCM.rds")
saveRDS(WUM,"WUM.rds")
```



### 2.2 Sex combination  matrix

```R
#Create data frame with each Individual name (character) and their Age (Character)
sex_frame<-metadata[,c("identifier","sex")]
sex_frame$identifier<-as.character(sex_frame$identifier)
sex_frame$sex<-as.character(sex_frame$sex)

#Create an empty character matrix to fill with characters
SEXM<-array(as.character(NA),c(nrow(sex_frame),nrow(sex_frame)))

for(i in 1:nrow(sex_frame)){
  for(j in 1:nrow(sex_frame)){ 
    if(sex_frame$sex[i]=="F" & sex_frame$sex[i]==sex_frame$sex[j]){
      SEXM[i,j]= "FF"}
    if(sex_frame$sex[i]=="M" & sex_frame$sex[i]==sex_frame$sex[j]){
      SEXM[i,j]= "MM"}
    if( sex_frame$sex[i]!=sex_frame$sex[j]){
      SEXM[i,j]= "FM"}}}

#Name rown amd colnames with individual names 
rownames(SEXM)<-key$ID
colnames(SEXM)<-key$ID

#Save matrix to ready matrices folder
saveRDS(SEXM,"SEXM-combo.rds") 
```


### 2.3 Species combination matrix

```R
#Create data frame with each Individual name (character) and their species_codeent (Character)
species_frame<-metadata[,c("identifier","species_code")]
species_frame$identifier<-as.character(species_frame$identifier)
species_frame$species_code<-as.character(species_frame$species_code)

# Change species_codeent levels from water,control,mal to W,C,M

species_frame$species_code <- with(species_frame, factor(species_code, levels = c('KP', 'WP', 'MP'), labels = c("CP", "CM", "CT"))) 

# B: before species_codeent; A:after species_codeent; C: control

#Create an empty character matrix to fill with characters
SPECIESM<-array(as.character(NA),c(nrow(species_frame),nrow(species_frame)))

for(i in 1:nrow(species_frame)){
  for(j in 1:nrow(species_frame)){
    if(species_frame$species_code[i]=="CP" & species_frame$species_code[i]==species_frame$species_code[j]){
      SPECIESM[i,j]= "CPCP"}
    if(species_frame$species_code[i]=="CM" & species_frame$species_code[i]==species_frame$species_code[j]){
      SPECIESM[i,j]= "CMCM"}
    if(species_frame$species_code[i]=="CT" & species_frame$species_code[i]==species_frame$species_code[j]){
      SPECIESM[i,j]= "CTCT"}
    if( species_frame$species_code[i]=="CP" & species_frame$species_code[j]=="CM"){
      SPECIESM[i,j]= "CPCM"}
    if( species_frame$species_code[i]=="CM" & species_frame$species_code[j]=="CP"){
      SPECIESM[i,j]= "CPCM"}
    if( species_frame$species_code[i]=="CP" & species_frame$species_code[j]=="CT"){
      SPECIESM[i,j]= "CPCT"}
    if( species_frame$species_code[i]=="CT" & species_frame$species_code[j]=="CP"){
      SPECIESM[i,j]= "CPCT"}
    if( species_frame$species_code[i]=="CM" & species_frame$species_code[j]=="CT"){
      SPECIESM[i,j]= "CMCT"}
    if( species_frame$species_code[i]=="CT" & species_frame$species_code[j]=="CM"){
      SPECIESM[i,j]= "CMCT"}}}

#Name rown and colnames with individual names 
rownames(SPECIESM)<-key$ID
colnames(SPECIESM)<-key$ID

#Save matrix to ready matrices folder
saveRDS(SPECIESM,"SPECIES-combo.rds") 
```



### 2.4 Age combination matrix

```R
#Create data frame with each Individual name (character) and their Age (Character)
age_frame<-metadata[,c("identifier","age")]
age_frame$identifier<-as.character(age_frame$identifier)
age_frame$age<-as.character(age_frame$age)

#Create an empty character matrix to fill with characters
AGEM<-array(as.character(NA),c(nrow(age_frame),nrow(age_frame)))

for(i in 1:nrow(age_frame)){
  for(j in 1:nrow(age_frame)){ 
    if(age_frame$age[i]=="A" & age_frame$age[i]==age_frame$age[j]){
      AGEM[i,j]= "AA"}
    if(age_frame$age[i]=="J" & age_frame$age[i]==age_frame$age[j]){
      AGEM[i,j]= "JJ"}
    if(age_frame$age[i]!=age_frame$age[j]){
      AGEM[i,j]= "AJ"}}}

#Name rown and colnames with individual names 
rownames(AGEM)<-key$ID
colnames(AGEM)<-key$ID

#Save matrix to ready matrices folder
saveRDS(AGEM,"AGE-combo.rds") 
```




### 2.5 Nest sharing  matrix (nest similarity) 

```R
#Create data frame with each Individual name (character) and their nest ID (Character)
nest_frame<- metadata[,c("identifier","nest")]
nest_frame$identifier<-as.character(nest_frame$identifier)
nest_frame$nest<-as.character(nest_frame$nest)

#Create an empty numeric matrix to fill with distances
nestM<-array(0,c(nrow(nest_frame),nrow(nest_frame)))

#Derive matrix with binary Age similarity between each sample
for(i in 1:nrow(nest_frame)){
  for(j in 1:nrow(nest_frame)){ 
    if(nest_frame$nest[i]==nest_frame$nest[j]){
      nestM[i,j]= 1
    } else{
      nestM[i,j]= 0}}} 

#Name rown and colnames with individual names 
all(rownames(nestM)==key$ID)
rownames(nestM)<-key$ID
colnames(nestM)<-key$ID

nestM

#Save matrix to ready matrices folder
saveRDS(nestM,"nestM.rds")
```



### 2.6 Year similarity matrix 

```R
#Create data frame with each Individual name (character) and their nest ID (Character)
year_frame<- metadata[,c("identifier","year")]
year_frame$identifier<-as.character(year_frame$identifier)
year_frame$year<-as.character(year_frame$year)

#Create an empty numeric matrix to fill with distances
yearM<-array(0,c(nrow(year_frame),nrow(year_frame)))

#Derive matrix with binary Age similarity between each sample
for(i in 1:nrow(year_frame)){
  for(j in 1:nrow(year_frame)){ 
    if(year_frame$year[i]==year_frame$year[j]){
      yearM[i,j]= 1
    } else{
      yearM[i,j]= 0}}} 

#Name rown and colnames with individual names 
all(rownames(yearM)==key$ID)
rownames(yearM)<-key$ID
colnames(yearM)<-key$ID


#Save matrix to ready matrices folder
saveRDS(yearM,"yeartM.rds")
```



### 2.7 Distance (meters) matrix between individuals  
```R
# Create a dataframe with your data
distance_df <- metadata[, c("identifier", "utm_1", "utm_2")]
distance_df <- na.omit(distance_df) #remove NA
#individuals_to_prune <- c("FJ02235", "FJ02223", "FH80666", "FH69257")#no gps points

# Create a SpatialPointsDataFrame

coordinates <- distance_df[, c("utm_1", "utm_2")]
proj4string <- CRS("+proj=utm +zone=38K +datum=WGS84") # coordinate reference system
spdf <- SpatialPointsDataFrame(coordinates, data = distance_df, proj4string = proj4string)

# Calculate distances
distance_matrix <- spDists(spdf)

# Assign proper sample names to the matrix
key_distance<-data.frame(identifier=distance_df$identifier)
rownames(distance_matrix)<-key$identifier # assign ring_number to rows
colnames(distance_matrix)<-key$identifier#assign rin_number to columns

#Save matrix to ready matrices folder
saveRDS(distance_matrix,"DistM.rds")
```





## 3. Unravel matrices into one dyadic data frame



### 3.1 Build dyadic dataset

```R
#Read in microbial distance matrices if not in already
BCM <- readRDS("BCM.rds") # bray-curtis
WUM <- readRDS("WUM.rds") # weighted unifrac

# Ready in matrices of other variables
speciesM_combo <- readRDS("SPECIES-combo.rds")
speciesM_sim <- readRDS("SPECIES-sim.rds")
ageM_combo<-readRDS("AGE-combo.rds") # age
ageM_sim <-readRDS("AGE-sim.rds") # age
sexM_combo <- readRDS("SEXM-combo.rds")
sexM_sim <- readRDS("SEXM-sim.rds")
nestM <- readRDS("nestM.rds")
yearM <- readRDS("yeartM.rds")
distM <- readRDS("DistM.rds")

#First unravel the matrices into vectors matching the lower quantile of each matrix. 
#From numeric matrices, this can be done by making a list (c()) of the distance object (dist()) derived from the matrix. 
#as.dist() by default includes only the lower quantile of the matrix and excludes the diagonal.
#From categorical matrices, this can be done by making a list (c()) of the lower quantile of the matrix with lower.tri() -function.

BC <- c(as.dist(BCM))
WU <- c(as.dist(WUM))
species_combo <- c(speciesM_combo[lower.tri(speciesM_combo)])
species_sim <- c(speciesM_sim[lower.tri(speciesM_sim)])
sex_combo <- c(sexM_combo[lower.tri(sexM_combo)])
sex_sim <- c(sexM_sim[lower.tri(sexM_sim)])
age_combo <- c(ageM_combo[lower.tri(ageM_combo)])
age_sim <- c(ageM_sim[lower.tri(ageM_sim)])
nest <- c(nestM[lower.tri(nestM)])
year <- c(yearM[lower.tri(yearM)])
distance <-c(as.dist(distM))


#Combine these vectors into a data frame
data.dyad<-data.frame(BC=BC , WU=WU, species_combo=species_combo, species_sim=species_sim, age_combo=age_combo, age_sim=age_sim,  
                      sex_combo=sex_combo, sex_sim=sex_sim, nest=nest, year=year)

#Add the identities of both individuals in each dyad as separate columns into the data frame and exclude 
#self-comparisons (as these are not meaningful). 
```



### 3.2 Add  Sample ID combinations to the data set

```R
# extracting Individual-combinations present in the matrices
list<-expand.grid(key$ID,key$ID) 

# This created individual-to-same-individual pairs as well. Get rid of these
list<-list[which(list$Var1!=list$Var2),] 

# this still has both quantiles in--> add 'unique' key 
list$key <- apply(list, 1, function(x)paste(sort(x), collapse='')) 
list<-subset(list, !duplicated(list$key)) 

# sanity check that the Individual name combinations are in the same exact order as the lower quantile value vector of the matrices
i=34
BCM[which(rownames(BCM)==list$Var1[i]),which(colnames(BCM)==list$Var2[i])]==BC[i]
WUM[which(rownames(WUM)==list$Var1[i]),which(colnames(WUM)==list$Var2[i])]==WU[i]
ageM_combo[which(rownames(ageM_combo)==list$Var1[i]),which(colnames(ageM_combo)==list$Var2[i])]==age_combo[i]
sexM_combo[which(rownames(sexM_combo)==list$Var1[i]),which(colnames(sexM_combo)==list$Var2[i])]==sex_combo[i]

# add the names of both individuals participating in each dyad into the data frame
data.dyad$IDA<-list$Var2
data.dyad$IDB<-list$Var1


data.dyad$sex_combo <- as.factor(data.dyad$sex_combo)
data.dyad$species_combo <- as.factor(data.dyad$species_combo)
data.dyad$age_combo <- as.factor(data.dyad$age_combo)
data.dyad$age_sim <- as.factor(data.dyad$age_sim)
data.dyad$nest <- as.factor(data.dyad$nest)
data.dyad$year <- as.factor(data.dyad$year)
str(data.dyad)
                  
saveRDS(data.dyad, "data.dyad.rds")
```


