library(rgdal)
library(raster)
###DIVERSITY MAPS

    ###First stack rasters with points
### function to rasterize points
rasterize_species= function (x,mask=selected_mask,r=raster_base) {
 # r<-raster(ncol=300,nrow=400,resolution=resolution,ext=extent(selected_mask))
 # res(r)<-resolution 
  r<-extend(r,selected_mask)
  r<-crop(r,extent(selected_mask))
  r<-mask(r,selected_mask)
  values(r)<-0
  map<-readOGR(dsn=maps_folder,layer=x)
  r<-rasterize(map,r,1,update=T,background=0)
  r<-mask(r,selected_mask)
  valor<-unique(getValues(r))
  
  if(length(valor)==1&&is.na(valor)==TRUE){
    
    
  }
  else if (length(valor)==2&&valor[2]==0){
    
  }
  else {
    writeRaster(r,paste(x,".asc",sep=""))
    return (raster(paste(x,".asc",sep="")))
  }
}
###load mask to use (AF) ##change path
selected_mask<-readOGR("/Users/andreapaz/Dropbox/Guano's_files/AF_shapefile_from_Guano","atlantic_forest") 
selected_mask = aggregate(selected_mask, dissolve = TRUE)
###list species
setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/SDM/Points_only/") ##to folder with points


resolution<-0.008333333 
maps_folder<-getwd()
distribution_files<-list.files(path=getwd(), pattern= "*.shp$")
species_names_points<-sub(".shp","",distribution_files) ##really file names
raster_base<-raster(models[139]) ### this needs to be created with a different path
#Rasterize distribution files and keep only those of interest
layers<-lapply(species_names_points,rasterize_species)
##get species names from with grep removing the Points part. 
species_names_points<-sub("Points_","",species_names_points)
#assign species names to layers with species distributions
layers[sapply(layers,    is.null)]<-NULL
Stack_points<-stack(layers)
###load SDMs and make correct size
setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/SDM/OR_AUC/OR_AUC_T10/") ##to folder with SDMs
models<-list.files()
for(i in 1:length(models)){
  r<-raster(models[i])
  #par(mfrow=c(1,4))
  #plot(r)
  r<-extend(r,selected_mask)
 # plot(r)
  r<-crop(r,extent(selected_mask))
  #plot(r)
  r<-mask(r,selected_mask)
 # plot(r)
  writeRaster(r,paste("AF_",models[i],sep=""),format="ascii")
}
##when everything is created just:
setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/SDM/OR_AUC/OR_AUC_T10/")
models<-list.files(pattern="AF_")
Stack_models<-stack(models)
species_names_models<-sub("_OR_AUC_T10","",models) 
species_names_models<-sub("AF_","",species_names_models)
species_names_models<-sub(".asc","",species_names_models,fixed=T)
setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/SDM/Points_only/") 
points<-list.files(pattern=".asc")
Stack_points<-stack(points)
species_names_points<-sub(".asc","",points,fixed=T) 
species_names_points<-sub("Points_","",species_names_points)
species_names_full<-c(species_names_points,species_names_models)

##To get richness and community composition
Stack_full<-stack(Stack_points,Stack_models)
TD<-calc(Stack_full,sum,na.rm=T)
TD<-crop(TD,selected_mask)
TD<-mask(TD,selected_mask)
plot(TD)
writeRaster(TD,"richness_hylids_complete_0_0083.asc")

#### Phylogenetic Diversity

####Load Phylogeny, generate list of species
##To read a phylogeny in newick format use read.tree instead of read.nexus
#In working directory always a file with the same name: phylogeny.nex if not then change name here
setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/PD_data/")
##read phylogeny
user_phylogeny<-ape::read.tree("Jetz_Pyron_amph_shl_new_Consensus_7238.tre")

####MMake sure the names match between maps and phylogeny
setdiff(user_phylogeny$tip.label, species_names_full)
 ##if it doesnt use this code
###Trim phylogeny to match distribution data (remove non hylids)
pruned.tree<-ape::drop.tip(user_phylogeny, setdiff(user_phylogeny$tip.label, species_names_full))
##check that all disrtribution data is in tree
test<-as.data.frame(species_names_full)
rownames(test)<-species_names_full
check_names<-geiger::name.check(pruned.tree, test, data.names=NULL) #this must be OK
### Trim the distribution data to match phylogeny
#species_names_full<-species_names_full[species_names_full %in% user_phylogeny$tip,label]
#distribution_files<-sub("*$","_OR_AUC_T10.asc",species_names)

#########################
####START ANALYSES#######
#########################

tabla<-as.data.frame(species_names_full)
colnames(tabla)<-"Grilla"

#Create empty raster of the desired area and resolution to assign pixel numbers
r<-Stack_full[[1]]
res(r)<-res(Stack_full) #resolution
#r[is.na(r[])]<-0
#r<-crop(r,selected_mask)
#r<-mask(r,selected_mask)
#r<-aggregate(r,fact=10,fun=max)
grilla=r
names(grilla)="grilla"
grilla[1:ncell(grilla)]<-1:ncell(grilla)
lista_grilla<-list(grilla)
names(lista_grilla)<-"Grilla"
#Change names of models tu species names and add the empty raster in the beggining (pixel number)
names(Stack_full)<-as.vector(tabla$Grilla)
Stack<-stack(lista_grilla$Grilla,Stack_full) ###Full stack of all maps
#
#Turn maps into dataframe for computation of PD
community_data<-as.data.frame(Stack)
##remove rows that are NA 
community_data1<-community_data[-which(rowSums(!is.na(community_data[,2:159]))==0),]
species_names<-colnames(community_data1)[2:length(community_data1)] #Store species names


#Store 
setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/")
write.table(community_data1, file = "communities_hylids.txt", append = FALSE,row.names=F,quote=F,sep="\t") #tratar de agregar nombre de mascara

community_data1[is.na(community_data1)] = 0
#In the community data frame NA must be eliminated done before is this if you load community data and not maps?
#community_data=na.omit(community_data)
#head(community_data)

#III-Phylogenetic diversity computation 
#computes only Faith�s PD others may be added

pd.result <-picante::pd(community_data1[,2:ncol(community_data1)],pruned.tree,include.root = F) 
##if it fails switch to
##PhyloMeasures::pd.query(tree, matrix,)

#Add the pixel PD value to data frame
community_data1$pd<-pd.result[[1]]

#Write the new matrix to a file to avoid rerunning the script for potential further analyses
write.table(community_data, file = "communities_hylids_and_pd.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#Generate a raster containing PD information per pixel

#1-First generate an empty raster using a base model for resolution and area

values(r)<-0
pd_ras<-r
values(pd_ras)<-NA #Eliminate every value they will be replaced by PD values further down


#2- Assign PD values to raster
pd_ras[community_data1$grilla]<-community_data1$pd

#3- Save raster to file 

writeRaster(pd_ras,"hylids_PD.asc",format="ascii")

#4-Optional plotting map in R 
plot(pd_ras)


###FUNCTIONAL DIVERSITY


# Function to calculate the four Functional diversity indices (So far we only use the FRic index since no abundance data is available from distribution maps)
library(geometry)
library(ape)
library(ade4)

##Load maps exactly as for PD if not yet loaded
setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/SDM/OR_AUC/OR_AUC_T10/")
models<-list.files(pattern="AF_")
Stack_models<-stack(models)
species_names_models<-sub("_OR_AUC_T10","",models) 
species_names_models<-sub("AF_","",species_names_models)
species_names_models<-sub(".asc","",species_names_models,fixed=T)
setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/SDM/Points_only/") 
points<-list.files(pattern=".asc")
Stack_points<-stack(points)
species_names_points<-sub(".asc","",points,fixed=T) 
species_names_points<-sub("Points_","",species_names_points)
species_names_full<-c(species_names_points,species_names_models)

##To get richness and community composition
Stack_full<-stack(Stack_points,Stack_models)

##polygon with mask is called selected_maks from above


####Communities####



#Trait data: select file containing trait data
trait<-read.csv("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/Species_traits_imputed_lilian.csv",h=T)
###Trim trait data to match distribution data #2 don t match the ocean ones
trait<-trait[trait$Species %in% species_names_full,]
#trait<-na.omit(trait)
### Trim the distribution datta to match trait data , This is not needed in this case
#species_names_full<-species_names_full[species_names_full %in% trait1$Species]

tabla<-as.data.frame(species_names_full)
colnames(tabla)<-"Grilla"
#Load all distribution maps and make them same extent as will be used
r<-Stack_full[[1]]
#r<-extend(r,AF,value=0)
#r[is.na(r[])]<-0
#r<-crop(r,AF)
#r<-mask(r,AF)
#r<-aggregate(r,fact=10,fun=max)
grilla=r
names(grilla)="grilla"
grilla[1:ncell(grilla)]<-1:ncell(grilla)
lista_grilla<-list(grilla)
names(lista_grilla)<-"Grilla"
#layers<-lapply(distribution_files,load_species)
#names(layers)<-as.vector(tabla$Grilla)
#layers[sapply(layers,is.null)]<-NULL
#Change names of models tu species names and add the empty raster in the beggining (pixel number)
names(Stack_full)<-as.vector(tabla$Grilla)
Stack<-stack(lista_grilla$Grilla,Stack_full) ###Full stack of all maps

#Turn maps into dataframe for computation of FD
marco<-as.data.frame(Stack)
marco<-marco[-which(rowSums(!is.na(marco[,2:159]))==0),]
marco[is.na(marco)] = 0
species_names<-colnames(marco)[2:length(marco)] #Store species names
marco1<-marco
marco<-marco[-which(rowSums(marco[,2:159])==0),]




###compute Functional dispersion including categorical data
#Trait data: select file containing trait data
trait1<-read.csv("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/Species_traits_imputed_lilian.csv",h=T)
#FD_dis<-FD::gowdis(trait1[,c(3,4,7,8)])
#FDDis<-FD::fdisp(FD_dis,marco)
##Trim trait data to match distribution data #2 don t match the ocean ones
#trait1<-trait1[trait1$Species %in% names(marco)[2:159],]
rownames(trait1)<-trait1[,1]
trait1<-trait1[,c(4,7,8)]
for (i in 1:3)
{
  trait1[,i]<-log(trait1[,i]) ## This could be updated to another mathematic transformation if needed
}

FDindices<-FD::dbFD(trait1,marco[,2:159],calc.FDiv=F)
map<-raster(Stack_full[[1]])
fddis_ras<-r
fdric_ras<-r
values(fddis_ras)<-NA # Erase alll values from the distribution map
values(fdric_ras)<-NA
#Assign to the empty raster the values of FD that correspond to each pixel
marco$fdDis<-FDindices$FDis
marco$fdRic<-FDindices$FRic
fddis_ras[marco$grilla]<-marco$fdDis
fdric_ras[marco$grilla]<-marco$fdRic
writeRaster(fddis_ras,"hylid_functional_dispersion_cont.asc", format="ascii")
writeRaster(fdric_ras,"hylid_functional_richness_cont.asc", format="ascii")

###try with categorical

trait2<-read.csv("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/Species_traits_imputed_lilian.csv",h=T)
rownames(trait2)<-trait2[,1]
trait2<-trait2[,c(3,4,7,8)] #rep mode, size, head w tibiaL
for (i in 2:4)
{
  trait2[,i]<-log(trait2[,i]) ## This could be updated to another mathematic transformation if needed
}
##check if cat is cat
FDindices_cat<-FD::dbFD(trait2,marco[,2:159],calc.FDiv=F,corr="cailliez")
marco$fdDis<-FDDis[,2] ##check
###prob remove

# Select traits 
##remove non trait info like typo of model , for richness remove the categorical rep mode
##remove dpulicted size column using data from Celio , keep size, HW and TL
trait<-trait[,c(1,4,7,8)]
str(trait)
head(trait)
summary(trait)
traitnb=dim(trait)[2]-1 #Number of traits, minus 1 because the first column contains species names
spnb=dim(trait)[1] #Number of species (one species per row thus number of species=number of rows)
### log transform traits

for (i in 2:traitnb)
{
  trait[,i]<-log(trait[,i]) ## This could be updated to another mathematic transformation if needed
}


traitc=trait[,2:traitnb] #Select only columns containign trait info and not species names
traitcn=apply(traitc,2,as.numeric) #Ensure all trait values are numeric


###Confirm number of species with map=number of species with traits
####Confirm names of species with map are the same as names of species with traits
trait[,1]<-as.character(trait[,1])
differences<-setdiff(trait[,1],species_names)

# Calculate FRic, FEve, FDiv and FSpe
if (length(differences)==0){
  FD=as.numeric() #Create an empty vector to store functional diversity data
  # Community=unique(marco$grilla)
  FD.i=FDind(traitcn,marco[,1:spnb])
  for (i in 1:length(Community)) {
    
    abund.i=marco[marco$grilla==Community[i],1:spnb]
    if ((sum(abund.i)-abund.i[1])>=5){
      abundm.i=as.matrix(abund.i)
      
      FD.i=FDind(traitcn,abundm.i)
      FD=rbind(FD,FD.i)
    }
    else 
    {    FD.i=data.frame(Nbsp=NA, FRic=NA, FEve=NA, FDiv=NA, FSpe=NA) ; row.names(FD.i)=row.names(abund.i)
    FD=rbind(FD,FD.i)
    }
  }
  
  marco$fd<-FD[,2]
  map<-raster(Stack_full[[1]])
  fd_ras<-r
  values(fd_ras)<-NA # Erase alll values from the distribution map
  
  #Assign to the empty raster the values of FD that correspond to each pixel
  fd_ras[marco$grilla]<-marco$fd
  
  
  #Compute taxonomic diversity as number of species done before
  #  richness<-stack(layers)
  # TD<-calc(richness,sum,na.rm=T)
  
  
  #Plot both maps in R to the left functional diversity and to the right taxonomic diversity
  # par(mfrow=c(1,2))
  
  #plot(fd_ras, main="Functional Diversity")
  #  plot(TD, main="Species Richness")
  
  #Write rasters (FD and TD) to file
  setwd("~/")
  writeRaster(fd_ras,"hylid_functional_richness.asc", format="ascii")
  
  
} else print("There are differences in the species names in the trait and distriubtion databases:");differences

