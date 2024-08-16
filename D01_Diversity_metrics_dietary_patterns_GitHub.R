################################################
# Compute Functional Diversity and Phylogentic Species Metrics
# by Simonetta Selva
#
# Created: May 2nd, 2023
# Project: Bumblebee 2022
################################################



# Here, all Functional Diversity Metrics and the Phylogenetic Species Variability are calculated 
# per individual bumblebee or if wanted per site. They are needed in all other script. 
##### preparation ----
# clear environment
rm(list=ls())
#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(caret)
library(vegan)
library(reshape2)
library(FD)
library(GGally)
# Create directory for results
#dir.create("ANALYSES/OUTPUT/Intermediate_files/diversity.metrics/")

# Load data
Data.dietary.patterns.leg <- read_csv("BB22_pollen_composition_full_plant_traits.csv") %>% 
  mutate(ID = as.character(ID)) %>% dplyr::filter(bborgan=="L")
Data.dietary.patterns.leg$ID.short = as.factor(substring(Data.dietary.patterns.leg$ID,1, nchar(Data.dietary.patterns.leg$ID)-1)) #remove information on leg and body from ID
Data.dietary.patterns.leg %>% dplyr::group_by(bbspecies) %>% dplyr::summarise(n=n_distinct(species))

##### SPECIES RICHNESS --------
### Calculations
species.richness.leg=Data.dietary.patterns.leg %>% group_by(ID, location, landscape, replicate, bbspecies) %>% dplyr::summarise(species.richness=n_distinct(plant.species))
species.richness.leg$loc.id=paste(species.richness.leg$location, species.richness.leg$replicate, sep="")
### Export
write.csv(species.richness.leg, "species.richness.csv")

# FUNCTIONAL DIVERSITY METRICS ----
# choose only to numeric data
Data.dietary.patterns.leg.numeric <- Data.dietary.patterns.leg %>% 
  dplyr::reframe(ID.short = ID.short, #remove information on leg and body from ID
            ID = as.factor(ID),
            location = as.factor(location),
            landscape = as.factor(landscape),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            site = as.factor(site),
            plant.species = plant.species,
            binom.abund = binom.abund,
            abundance = Abundance,
            Flowering_duration = Flowering_months_duration,
            Flowering_start = start_flowering,
            growth_form_numeric = growth_form_numeric,
            structural_blossom_numeric = structural_blossom_numeric,
            sugar.concentration = sugar.concentration,
            symmetry_numeric = symmetry_numeric,
            plant_height_m = plant_height_m)
Data.dietary.patterns.leg.numeric$Flowering_duration <- as.numeric(Data.dietary.patterns.leg.numeric$Flowering_duration)

## Correlation analysis ----
# correlation analysis of the variables
M <- cor(Data.dietary.patterns.leg.numeric[,c(11:17)], 
         use = "complete.obs")

### FIGURE S3 ----
corrplot::corrplot(M, type="upper", order="hclust", 
                   col = COL2('RdBu', 10),mar = c(1, 1, 1, 1)) # plot correlation with p-values

# remove plant height, flowering start  and symmetry look at correlation
M <- cor(Data.dietary.patterns.leg.numeric[,c(11, 13:15)], use = "complete.obs")
corrplot::corrplot(M, type="upper", order="hclust", 
                   col = COL2('RdBu', 10),mar = c(1, 1, 1, 1)) # plot correlation with p-values

# remove plant height, flowering start from data set
Data.dietary.patterns.leg.numeric1 <- Data.dietary.patterns.leg.numeric%>% 
  dplyr::select(-Flowering_start,
                -plant_height_m,
                -symmetry_numeric) 

## compute FDs ----


# the following code has to be run 4 times with each time different input variables i and j 
# to calculate FD for each bumblebee species and on a site or ID level
# i = "B.lapidarius" or "B.pascuorum"
# j = "site" or "ID.short"
# It can also be performed in a loop, but due to computational time it is easier to them individually. 

FD.bb=list()
for(i in c("B.pascuorum", "B.lapidarius")){
  cat(i)
  # filter for bumblebee species, remove plant species entries that do not have traits and prepare data
  BB22_full.loop <- Data.dietary.patterns.leg.numeric1[Data.dietary.patterns.leg.numeric1$bbspecies == i,]%>%
    filter(plant.species!= "Fabaceae sp.", 
           plant.species!= "Cyclamen sp.",
           plant.species!= "Petunia sp.",
           plant.species!= "Mandevilla atroviolacea" )%>%
    mutate(Flowering_duration = as.numeric(Flowering_duration))%>% 
    droplevels() # remove duplicates
  
  # select columns used in computing FDs
  BB22_full.loop.species <- BB22_full.loop %>% 
    dplyr::select(plant.species, Flowering_duration, structural_blossom_numeric, sugar.concentration) %>% 
  distinct() # remove duplicates
  
  # impute missing data
  trt.mis.pred <- preProcess(as.data.frame(BB22_full.loop.species[,-c(1)]), "knnImpute")
  traits <- predict(trt.mis.pred, BB22_full.loop.species[,-c(1)]); head(trt.mis.pred)
  traits <- as.data.frame(traits)
  rownames(traits)  <- BB22_full.loop.species$plant.species
  traits <- traits[order(row.names(traits)), ] # reorder traits into alphabetical order (needed for FD package)
  
  # relative abundance data
  BB22_full.ab <- BB22_full.loop%>%
    group_by(ID, plant.species, location, landscape, replicate, bbspecies)%>%
    summarise(abundance = sum(abundance))%>% 
    distinct() # remove duplicates
  
  # bring into wide format
  wide <- dcast(BB22_full.ab, ID ~ plant.species, value.var="abundance")[,-1]
  rownames(wide)  <- levels(BB22_full.ab$ID)
  wide[is.na(wide)] <- 0
  
  # create matrix (abundance and presence/absence)
  # abundance
  sp.ab <- as.matrix(wide)
  write.table(sp.ab, file = paste("community_matrix_", i, ".txt", sep=""), sep = ",")
  
  # presence/absence
  sp.pa <- decostand(wide, "pa")
  sp.pa <- as.matrix(sp.pa) #turn into matrix
  
  # compute FD
  fd.weig <- FD::dbFD(x = traits , a = sp.ab, w.abun = T) # weighted 
  fd.bino <- FD::dbFD(x = traits , a = sp.pa) # not weighted 
  species.richness.leg.bb=species.richness.leg[species.richness.leg$bbspecies==i,]
  df.FD <- data.frame(
                      FRic.w = fd.weig$FRic, 
                      FRic = fd.bino$FRic, 
                      FEve.w = fd.weig$FEve, 
                      FEve = fd.bino$FEve, 
                      FDiv.w = fd.weig$FDiv, 
                      FDiv = fd.bino$FDiv)
  df.FD$bbspecies = paste(i)
  df.FD$ID=rownames(df.FD)
  df.FD=merge(df.FD, BB22_full.loop[, c(2,3,4,5)])
  df.FD=merge(df.FD, species.richness.leg.bb, by=c("ID", "location", "landscape", "replicate", "bbspecies"))
  FD.bb[[i]] = df.FD
} 

FD.bb.unlist=do.call(what = rbind, FD.bb)
## Export
write.csv(FD.bb.unlist, "Functional.metrics.csv")


# PHYLOGENETIC SPECIES METRICS ----

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(V.PhyloMaker2)
library(picante)

# import phylo tree from (file "Diet Pattern")
phylo <- data.frame(species = Data.dietary.patterns.leg$species, genus = Data.dietary.patterns.leg$genus, family = Data.dietary.patterns.leg$family)

# 2. phylogentic tree
# phylogenetic hypotheses under three scenarios based on a backbone phylogeny
tree <- phylo.maker(phylo) # new version of package

 
## compute PD ----
# loop for two species and two resolution levels (site and ID)
species <- c("B.lapidarius", "B.pascuorum")

PD.bb=list()
for (i in species) {
    cat(i)
  BB22_full.loop <- Data.dietary.patterns.leg.numeric1[Data.dietary.patterns.leg.numeric1$bbspecies == i,]
    sp.ab <- read.delim(paste("community_matrix_", i, ".txt", sep=""), sep = ",", 
                        check.names = F)
    colnames(sp.ab) <- sub(" ", "_", colnames(sp.ab)) # match plant species names
    
    PD_var = psv(sp.ab,tree,compute.var=TRUE,scale.vcv=TRUE)
    PD_ric = psr(sp.ab,tree,compute.var=TRUE,scale.vcv=TRUE)
    PD_eve = pse(sp.ab,tree,scale.vcv=TRUE)
    PD_clu = psc(sp.ab,tree,scale.vcv=TRUE)
    PD.bb.temp = data.frame(ID = rownames(PD_var),
               nbsp = PD_var$SR,
               PVar = PD_var$PSVs,
               PRic = PD_ric$PSR,
               PEve = PD_eve$PSEs,
               PClu = PD_clu$PSCs,
               bbspecies=paste(i))
    PD.bb.temp=merge(PD.bb.temp, BB22_full.loop[, c(2,3,4,5)])
    PD.bb[[i]]=PD.bb.temp
}

PD.bb.unlisted=do.call(rbind, PD.bb)

write.csv(PD.bb.unlisted, "Phylogenetic.metrics.csv")
