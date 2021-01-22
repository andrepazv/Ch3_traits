### Get areas for all species and summ stats per group
setwd("/Users/andreapaz/Dropbox/**Tesis_PHD/Mata_atlantica_Diversity_correlates/Alphahulls/Alpha_maps/")
#list all folders (one per group)
groups<-list.files()
    #Create data frame to store results
data_areas<- data.frame(NA,NA,NA)
colnames(data_areas)<-c("group","species","area")
for (i in groups){
  setwd("/Users/andreapaz/Dropbox/**Tesis_PHD/Mata_atlantica_Diversity_correlates/Alphahulls/Alpha_maps/")
setwd(i)
models<-list.files(pattern="map_.*\\.asc",full.names=T)
for(j in models){
  x<-raster(j)
  crs(x)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  x[x==0]<-NA
  y<-raster::area(x,na.rm=T)  
  y<-sum(getValues(y),na.rm=T) 
  ## 
 new<- c(i,j,y)
  data_areas[nrow(data_areas) + 1, ] <- new
}
}
data_areas<-data_areas[2:length(data_areas$group),]
data_areas$area<-as.numeric(data_areas$area)

write.csv(data_areas,"all_areas.csv")


###GET TOTAL area of distribution
all_maps<-list.files(pattern="complete",full.names = TRUE)
for (i in all_maps){
  x<-raster::raster(i)
  raster::crs(x)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  x[x==0]<-NA
  y<-raster::area(x,na.rm=T)  
  y<-sum(raster::getValues(y),na.rm=T) 
  print(paste(i," area is: ",y))
}
 
    