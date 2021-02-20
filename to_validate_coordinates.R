####Validate Hylidae points for Ch2
library(dismo)
library(sp)
library(rgdal)
library(rgeos)
setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/data_curation/")
    ###No points in oceans (validated by eye). Manually fixed accents, ñ,  and info on loc description
    
###Load points and make into shapefile
      ##origi
points<-read.table("merged_all_hylids_from_leyla.txt",h=T)

      ##extras
points<-read.table("all_species_to_complete.txt",h=T)

points_coord<-subset(points,!is.na(Latitude))
points_export<-points_coord
coordinates(points_coord) <- ~Longitude+Latitude
###Load admin boundaries world
setwd("~/Desktop/PC")
world_adm<-readOGR(".","world_adm0")
crs(points_coord) <- crs(world_adm)
plot(world_adm)
points(points_coord)

  ###Validate political boundaries
      
ovr <- over(points_coord, world_adm)
cntr <- ovr$NAME
cntr<-as.character(cntr)

    ###Ocean?
    i <- which(is.na(cntr))
    i
      ####- Country level 0 division
    j <- which(cntr != as.character(points_coord$Pais))
    cbind(cntr, as.character(points_coord$Pais))[j,]
      #### Load Brazil divisions for Department/canton level 1 division
    setwd("~/Desktop/PC/temp/")
    level1<-readOGR(".", "adm_1_cut") #this must change
    crs(points_coord) <- crs(level1)
    ovr1 <- over(points_coord, level1)
    cntr1 <- ovr1$NameR 
    cntr1<-as.character(cntr1)
    h<-which(cntr1 != as.character(points_coord$X1st_divR))
   
      #### Province + assign missing (in Col this wil be municipio) in any case lavel 2 division
 ### Validate species IUCN distribution
    
    ###Visualize problems
    plot(points_coord)
    plot(world_adm, add=T, border='blue', lwd=2)
    points(points_coord[j, ], col='red', pch=20, cex=2) ##wrong country
    points(points_coord[i,],col="blue",pch=20,cex=2) ##ocean points
    points(points_coord[h,],col="pink",pch=20,cex=2)
    
    ###ADD Flags of problems!
    points_export$Country_Flag<-rep(NA,length(points_export$Genero))
    points_export$Country_Flag[j]<-cntr[j]
    points_export$Ocean_Flag<-rep(NA,length(points_export$Genero))
    points_export$Ocean_Flag[i]<-cntr[i]
#    points_export$Dept_Flag<-rep(NA,length(points_export$FID))
 #   points_export$Dept_Flag[h]<-cntr1[h]
  #  points_export$Div1_map<-rep(NA,length(points_export$FID))
  #  points_export$Div1_map<-ovr1$NameR
   # bla<-subset(points_export,!is.na(Dept_Flag))
    #bla1<-subset(bla,X2nd_division_municipality_canton_province!="Guapi")
    #bla[,c("FID","X1st_divR","Dept_Flag","Country")]
    
    ###Export flags
   # setwd("~/Dropbox/mass extinctions/Archivos Andrea")
    write.csv(points_export,"Base_curada_para_pegar.csv")
    
    
    #####Validate species distributions
    
      ###Load all distribution shapefiles from folder
    setwd("/Users/andreapaz/Dropbox/**Tesis_PHD/Ch2_traits/data_curation/all_hylids/")
    maps<-list.files(pattern="*.shp")
    IUCN_poly<-lapply(maps,readOGR)
      ###test
    points_export$Species_Flag<-rep(NA,length(points_export$Genero))
    points_export$Species_GbufferFlag<-rep(NA,length(points_export$Genero))
    tabla_puntos<-data.frame(rep(NA,length(unique(points_export$Genus_species))))
    tabla_puntos$Especie<-unique(points_export$Genus_species)
    tabla_puntos$Numero_puntos<-rep(NA, length(tabla_puntos$Especie))
    tabla_puntos<-tabla_puntos[,2:3]
for(i in 1:length(IUCN_poly)){
  species<-IUCN_poly[[i]]
  name<- gsub(" ", "_",IUCN_poly[[i]]$BINOMIAL)
  name<-name[1]
  db<-subset(points_export,Genus_species==name)
  if(length(db$Genero)>0)
     {
    coordinates(db) <- ~Longitude+Latitude
  crs(db) <- crs(species)
  ovr <- over(db,species)
  ovr1<-over(db,gBuffer(species,width=0.1,byid=T))
  cntr <- ovr$BINOMIAL
  cntr<-as.character(cntr)
  cntr1<-ovr1$BINOMIAL
  cntr1<-as.character(cntr1)
  k <- which(is.na(cntr))
  k1 <- which(is.na(cntr1))
  position<-db@data$Especie[k]
  position1<-db@data$Especie[k1]
  index<-which(!is.na(match(points_export$Especie,position)))
  index1<-which(!is.na(match(points_export$Especie,position1)))
  points_export$Species_Flag[index]<-"Outside_of_range"
  points_export$Species_GbufferFlag[index1]<-"Outside_of_range"
  tabla_puntos$Numero_puntos[tabla_puntos$Especie== unique(db@data$Genus_species)]<-length(db@data$Genus_species)
 # if(length(index1)>0){
  pdf(paste("~/Dropbox/**Tesis_PHD/Ch2_traits/data_curation/New_species_dist/", name, ".pdf", sep=""))
  buffered<-gBuffer(species,width=4,byid=T)
  plot(buffered,lty=0)
  plot(world_adm,add=T)
  plot(species,col="red",add=T)
  points(db,col="blue",pch=16,cex=0.5)
  dev.off()
#    }
  }
  else 
    print(paste(name,"not in db"))
}
    ###Export species flags
    setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/data_curation/")
    write.csv(points_export,"db_adding_species_flags.csv")
    write.csv(tabla_puntos,"Points_per_species_new.csv")
    
  
    }