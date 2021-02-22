#setwd and libraries
library(raster)

##Load diversity rasters (1km T10)
setwd("Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/")
pd_raster<-raster("hylids_PD.asc")
sr_raster<-raster("richness_hylids_complete_0_0083.asc")
fdRic_raster<-raster("hylid_functional_richness_cont.asc")
fdDis_raster<-raster("hylid_functional_dispersion_cat1k.asc")
stack_SR_PD<-stack(pd_raster,sr_raster)
stack_SR_FD<-stack(fdRic_raster,sr_raster)
stack_PD_FD<-stack(fdRic_raster,pd_raster)
stack_PD_FDdis<-stack(fdDis_raster,pd_raster)
srpd_df<-as.data.frame(stack_SR_PD)
srfd_df<-as.data.frame(stack_SR_FD)
pdfd_df<-as.data.frame(stack_PD_FD)
pdfddis_df<-as.data.frame(stack_PD_FDdis)
#Generate regression models
#LINEAR
linear_model_sr_pd<-lm(srpd_df)
linear_model_sr_fd<-lm(srfd_df)
linear_model_pd_fd<-lm(pdfd_df)
linear_model_pd_fddis<-lm(pdfddis_df)
#Compute and save residuals from both regression models
linear_residuals_sr_pd<-resid(linear_model_sr_pd)
linear_residuals_sr_fd<-resid(linear_model_sr_fd)
linear_residuals_pd_fd<-resid(linear_model_pd_fd)
linear_residuals_pd_fddis<-resid(linear_model_pd_fddis)

#Add residuals to dataframe
srpd_df$PD_linear_residuals<-rep(NA,length(srpd_df[,1]))
srpd_df$PD_linear_residuals[as.numeric(names(linear_residuals_sr_pd))]<-linear_residuals_sr_pd[names(linear_residuals_sr_pd)]
#add residuals of fd to data frame
srfd_df$FD_linear_residuals<-rep(NA,length(srfd_df[,1]))
srfd_df$FD_linear_residuals[as.numeric(names(linear_residuals_sr_fd))]<-linear_residuals_sr_fd[names(linear_residuals_sr_fd)]
#add residuals of pd/fd to data frame
pdfd_df$PDFD_linear_residuals<-rep(NA,length(pdfd_df[,1]))
pdfd_df$PDFD_linear_residuals[as.numeric(names(linear_residuals_pd_fd))]<-linear_residuals_pd_fd[names(linear_residuals_pd_fd)]
#add residuals of pd/fddis to data frame ##check
pdfddis_df$PDFDdis_linear_residuals<-rep(NA,length(pdfddis_df[,1]))
pdfddis_df$PDFDdis_linear_residuals[as.numeric(names(linear_residuals_pd_fddis))]<-linear_residuals_pd_fddis[names(linear_residuals_pd_fddis)]

###plot residuals of PD~TD against Frich
stack_SR_PD_FD<-stack(pd_raster,sr_raster,fdRic_raster)
srpdfd_df<-as.data.frame(stack_SR_PD_FD)
linear_model_sr_pd<-lm(srpdfd_df$hylids_PD~srpdfd_df$richness_hylids_complete_0_0083)
linear_residuals_sr_pd<-resid(linear_model_sr_pd)
srpdfd_df$PD_linear_residuals<-rep(NA,length(srpdfd_df[,1]))
srpdfd_df$PD_linear_residuals[as.numeric(names(linear_residuals_sr_pd))]<-linear_residuals_sr_pd[names(linear_residuals_sr_pd)]
plot(srpdfd_df$PD_linear_residuals,srpdfd_df$hylid_functional_richness_cont,xlab="Phylogenetic diversity res",ylab="Functional richness")
abline(linear_model_pd_fddis,col="red")
#Generate rasters with residual values for PD
linear_residuals_sr_pd_raster<-raster(pd_raster)
values(linear_residuals_sr_pd_raster)<-NA
values(linear_residuals_sr_pd_raster)<- srpd_df$PD_linear_residuals
name_group<-"residuals_SR_PD_1km_T10"
writeRaster(linear_residuals_sr_pd_raster,name_group,"ascii")
##FD
linear_residuals_sr_fd_raster<-fdRic_raster
values(linear_residuals_sr_fd_raster)<-NA
values(linear_residuals_sr_fd_raster)<- srfd_df$FD_linear_residuals
name_group<-"residuals_SR_FD_1km_T10"
writeRaster(linear_residuals_sr_fd_raster,name_group,"ascii")

##PD FD
linear_residuals_pd_fd_raster<-fdRic_raster
values(linear_residuals_pd_fd_raster)<-NA
values(linear_residuals_pd_fd_raster)<- pdfd_df$PDFD_linear_residuals
name_group<-"residuals_PD_FD_1km_T10"
writeRaster(linear_residuals_pd_fd_raster,name_group,"ascii")

##PD FDdis
linear_residuals_pd_fddis_raster<-fdDis_raster
values(linear_residuals_pd_fddis_raster)<-NA
values(linear_residuals_pd_fddis_raster)<- pdfddis_df$PDFDdis_linear_residuals
name_group<-"residuals_PD_FDdis_1km_T10"
writeRaster(linear_residuals_pd_fddis_raster,name_group,"ascii")
##use the stack to plot elev with residuals
plot(pdfddis_df$hylids_PD,pdfddis_df$hylid_functional_dispersion_cat1k,xlab="Phylogenetic diversity",ylab="Functional dispersion")
abline(linear_model_pd_fddis,col="red")
# par(mar=(c(4,4,1,1)))
#  plot(Tana_dataframe$dataframe_SR.Tanagers_RICHNESS_ALL,Tana_dataframe$dataframe_PD.Tanagers_PD,xlab="Tanagers species richness",ylab="Tanagers phylogenetc diversity")
#  abline(linear_model,col="red")


###Explore relationship with elevation
##load elev and regress all rasters to elev first show wether simple relationship is or not the same


setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/")
library(raster)
elev<-raster("wc2.1_30s_elev.tif")
lat<-raster::raster("AF_latitude.asc")
# mask<-rgdal::readOGR("/Users/andreapaz/Dropbox/Guano's_files/AF_shapefile_from_Guano","atlantic_forest") 

##load residuals rasters
pdSr_res<-raster("residuals_SR_PD_1km_T10.asc")
fdricSr_res<-raster("residuals_SR_FD_1km_T10.asc")
fdricPd_res<-raster("residuals_PD_FD_1km_T10.asc")
fddisPd_res<-raster("residuals_PD_FDdis_1km_T10.asc")

elev<-crop(elev,sr_raster)
extent(elev)<-extent(sr_raster)
stack_elev_SR<-stack(sr_raster,elev,lat)
stack_elev_PD<-stack(pd_raster,elev,lat)
stack_elev_FD<-stack(fdRic_raster,elev,lat)
stack_elev_FDdis<-stack(fdRic_raster,elev,lat)
elevSR_df<-as.data.frame(stack_elev_SR)
elevPD_df<-as.data.frame(stack_elev_PD)
elevFD_df<-as.data.frame(stack_elev_FD)
elevFDdis_df<-as.data.frame(stack_elev_FDdis)
###residuals
stack_elev_SRPD<-stack(pdSr_res,elev,lat)
stack_elev_SRFD<-stack(fdricSr_res,elev,lat)
stack_elev_PDFD<-stack(fdricPd_res,elev,lat)
stack_elev_PDFDdis<-stack(fddisPd_res,elev,lat)
elevSRPD_df<-as.data.frame(stack_elev_SRPD)
elevSRFD_df<-as.data.frame(stack_elev_SRFD)
elevPDFD_df<-as.data.frame(stack_elev_PDFD)
elevPDFDdis_df<-as.data.frame(stack_elev_PDFDdis)
#LINEAR regressions just to plot trend
#loess_model_srpd_elev<-loess(elevSRPD_df)
linear_model_sr_fd<-lm(srfd_df)
linear_model_pd_fd<-lm(pdfd_df)
##use the stack to plot elev with residuals
plot(elevSRPD_df$wc2.1_30s_elev,elevSRPD_df$residuals_SR_PD_1km_T10,xlab="Elevation (masl)",ylab="Residuals of the PD~SR regression")
#abline(loess_model_srpd_elev,col="grey")
plot(elevSRFD_df$wc2.1_30s_elev,elevSRFD_df$residuals_SR_FD_1km_T10,xlab="Elevation (masl)",ylab="Residuals of the FDrich~SR regression")
plot(elevPDFD_df$wc2.1_30s_elev,elevPDFD_df$residuals_PD_FD_1km_T10,xlab="Elevation (masl)",ylab="Residuals of the FDrich~PD regression")
plot(elevPDFDdis_df$wc2.1_30s_elev,elevPDFDdis_df$residuals_PD_FDdis_1km_T10,xlab="Elevation (masl)",ylab="Residuals of the FDdis~PD regression")

library(ggplot2)
ggplot()+
  geom_point(data = elevSRFD_df, aes(x = wc2.1_30s_elev, y = residuals_SR_FD_1km_T10, color =AF_latitude)) +
  theme_bw() +
  labs(y="Residuals of the PD~TD regression", x = "Elevation (masl)",color="Latitude") +
  scale_color_gradient(low="blue", high="red")

###same with residuals are those different?
srpd_df$PD_linear_residuals
srfd_df$FD_linear_residuals
pdfd_df$PDFD_linear_residuals


####explore elvation and null model result

significant_pixelxFDdis<-raster("Null_model/Significant_pixels_FDdis1k.asc")
temp_min<-raster("~/Dropbox/Andrea_Lab/trasnfer_819_PC/Functional_chap/bio_6.tif") #bio6
precip_min<- raster("~/Dropbox/Andrea_Lab/trasnfer_819_PC/Functional_chap/bio_14.tif") #bio14
isotherm<- raster("~/Dropbox/Andrea_Lab/trasnfer_819_PC/Functional_chap/bio_2.tif") #daily tem range bio2
precipCV<- raster("~/Dropbox/Andrea_Lab/trasnfer_819_PC/Functional_chap/bio_15.tif") #precip CV bio15
temp_min<-crop(temp_min,sr_raster)
extent(temp_min)<-extent(sr_raster)
precip_min<-crop(precip_min,sr_raster)
extent(precip_min)<-extent(sr_raster)
isotherm<-crop(isotherm,sr_raster)
extent(isotherm)<-extent(sr_raster)
precipCV<-crop(precipCV,sr_raster)
extent(precipCV)<-extent(sr_raster)
precipCV<-mask(precipCV,sr_raster)
stack_tempmin_FDdissign<-stack(significant_pixelxFDdis,temp_min)
temminFDdissig_df<-as.data.frame(stack_tempmin_FDdissign)
stack_precipmin_FDdissign<-stack(significant_pixelxFDdis,precip_min)
precipminFDdissig_df<-as.data.frame(stack_precipmin_FDdissign)

stack_isotherm_FDdissign<-stack(significant_pixelxFDdis,isotherm)
isothermFDdissig_df<-as.data.frame(stack_isotherm_FDdissign)

stack_precipCV_FDdissign<-stack(significant_pixelxFDdis,precipCV)
precipCVFDdissig_df<-as.data.frame(stack_precipCV_FDdissign)

plot(elevFDdissig_df$wc2.1_30s_elev,elevFDdissig_df$Significant_pixels_FDdis1k,xlab="Elevation (masl)",ylab="Significant pixels")
plot(temminFDdissig_df$bio_6,temminFDdissig_df$Significant_pixels_FDdis1k,xlab="Minimum temperature of the coldest month",ylab="Significant pixels")
plot(precipminFDdissig_df$bio_14,precipminFDdissig_df$Significant_pixels_FDdis1k,xlab="Minimum precipitation of the driest month",ylab="Significant pixels")
plot(isothermFDdissig_df$bio_2,isothermFDdissig_df$Significant_pixels_FDdis1k,xlab="Isothermality",ylab="Significant pixels")
plot(precipCVFDdissig_df$bio_15,precipCVFDdissig_df$Significant_pixels_FDdis1k,xlab="Precipitation seasonality (CV)",ylab="Significant pixels")


####Do analysis by elevation band

##!- Create a raster of latitudes
library(raster)
setwd("Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/")
sr_raster<-raster("richness_hylids_complete_0_0083.asc")
r <- sr_raster
lat <- r
xy <- coordinates(r)
lat[] <- xy[, 2]
writeRaster(lat,"AF_latitude",format="ascii")




sr_raster<-raster("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/richness_hylids_complete_0_0083.asc")
selected_mask<-rgdal::readOGR("/Users/andreapaz/Dropbox/Guano's_files/AF_shapefile_from_Guano","atlantic_forest") 
selected_mask = aggregate(selected_mask, dissolve = TRUE)

bioclim<-crop(bioclims,sr_raster)
bioclim<-mask(bioclim,selected_mask)
for(i in 1:19){
  writeRaster(bioclim[[i]],paste(names(bioclim)[i],".asc",sep=""))
}
setwd("/Users/andreapaz/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/Bioclim_AF/")
bioclim<-list.files()
bioclim<-raster::stack(bioclim)
setwd("/Users/andreapaz/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/")
lat<-raster::raster("AF_latitude.asc")
elev<-raster::raster("wc2.1_30s_elev.tif")
elev<-raster::crop(elev,sr_raster)
elev<-mask(elev,selected_mask)
extent(elev)<-extent(sr_raster)
Stack_intro<-raster::stack(bioclim,lat,elev)
Stack_intro_df<-as.data.frame(Stack_intro)

library(ggplot2)
#ggplot() + 
 # geom_point(data = Stack_intro_df, aes(x = Stack_intro_df$wc2.1_30s_elev, y = Stack_intro_df$bio_1, color =Stack_intro_df$AF_latitude)) +
  #theme_bw() +
  #labs(y="Annual Temperature ()", x = "Elevation (masl)",color="Latitude")+
  #scale_color_gradient(low="blue", high="red")

ggplot()+
  geom_point(data = Stack_intro_df, aes(x = Stack_intro_df$wc2.1_30s_elev, y = Stack_intro_df$bio_6, color =Stack_intro_df$AF_latitude)) +
  theme_bw() +
  labs(y="Minimum temperature of the coldest month", x = "Elevation (masl)",color="Latitude")+
 # scale_color_gradient(low="blue", high="red")
scale_color_gradient2(low="yellow",mid="blue", high="red",midpoint=-20)
ggplot()+
geom_point(data = Stack_intro_df, aes(x = Stack_intro_df$wc2.1_30s_elev, y = Stack_intro_df$bio_17, color =Stack_intro_df$AF_latitude)) +
  theme_bw() +
  labs(y="Precipitation of Driest Quarter", x = "Elevation (masl)",color="Latitude") +
#  scale_color_gradient(low="blue", high="red")
scale_color_gradient2(low="yellow",mid="blue", high="red",midpoint=-20)
#ggplot()+
 # geom_point(data = Stack_intro_df, aes(x = Stack_intro_df$wc2.1_30s_elev, y = Stack_intro_df$bio_12, color =Stack_intro_df$AF_latitude)) +
  #theme_bw() +
#  labs(y="Annual precipitation", x = "Elevation (masl)",color="Latitude") +
 # scale_color_gradient(low="blue", high="red")


#ggplot()+
 # geom_point(data = Stack_intro_df, aes(x = Stack_intro_df$wc2.1_30s_elev, y = Stack_intro_df$bio_11, color =Stack_intro_df$AF_latitude)) +
  #theme_bw() +
  #labs(y="Minimum temperature of the coldest quarter", x = "Elevation (masl)",color="Latitude") +
  #scale_color_gradient(low="blue", high="red")

ggplot()+
  geom_point(data = Stack_intro_df, aes(x = Stack_intro_df$wc2.1_30s_elev, y = Stack_intro_df$bio_5, color =Stack_intro_df$AF_latitude)) +
  theme_bw() +
  labs(y="Maximum temperature of the warmest month", x = "Elevation (masl)",color="Latitude") +
  scale_color_gradient2(low="yellow",mid="blue", high="red", midpoint=-20)

ggplot()+
  geom_point(data = Stack_intro_df, aes(x = Stack_intro_df$wc2.1_30s_elev, y = Stack_intro_df$bio_16, color =Stack_intro_df$AF_latitude)) +
  theme_bw() +
  labs(y="Precipitation of wettest quarter", x = "Elevation (masl)",color="Latitude") +
  scale_color_gradient2(low="yellow",mid="blue", high="red", midpoint=-20)

####to create graphs per altitudinal band 
##divide in deciles
library(dplyr)
all_lat<-getValues(lat)
all_lat<-as.data.frame(all_lat)
all_lat$lat_decile <- ntile(all_lat$all_lat, 10)  
for(i in 1:10){
  temp<-subset(all_lat,lat_decile==i)
  temp_extent<-extent(lat)
  temp_extent[3]<-min(temp$all_lat)
  temp_extent[4]<-max(temp$all_lat)
  raster_temp<-crop(lat,temp_extent)
  writeRaster(raster_temp,paste("lat_subset",i,".asc",sep=""),format="ascii")
}
##first load the shapefiles of latitude 1-10 (divided in deciles but that can be modified)
lats<-list.files(pattern="lat")[2:11]
##get and process elevation for extent
elev1<-raster("wc2.1_30s_elev.tif")
elev1<-crop(elev1,sr_raster)
extent(elev1)<-extent(sr_raster)
elev1<-mask(elev1,sr_raster)

##load residuals
pdSr_res1<-raster("residuals_SR_PD_1km_T10.asc")
fdricSr_res1<-raster("residuals_SR_FD_1km_T10.asc")
fdricPd_res1<-raster("residuals_PD_FD_1km_T10.asc")
fddisPd_res1<-raster("residuals_PD_FDdis_1km_T10.asc")
##then within a loop crop the raster and plot for each quantile the residuals
for (i in 1:length(lats)){
  decile<-raster::raster(lats[i])
  pdSr_res<-crop(pdSr_res1,decile)
  fdricSr_res<-crop(fdricSr_res1,decile)
  fdricPd_res<-crop(fdricPd_res1,decile)
  fddisPd_res<-crop(fddisPd_res1,decile)
  elev<-crop(elev1,decile)
 
###residuals
stack_elev_SRPD<-stack(pdSr_res,elev)
stack_elev_SRFD<-stack(fdricSr_res,elev)
stack_elev_PDFD<-stack(fdricPd_res,elev)
stack_elev_PDFDdis<-stack(fddisPd_res,elev)
elevSRPD_df<-as.data.frame(stack_elev_SRPD)
elevSRFD_df<-as.data.frame(stack_elev_SRFD)
elevPDFD_df<-as.data.frame(stack_elev_PDFD)
elevPDFDdis_df<-as.data.frame(stack_elev_PDFDdis)
##use the stack to plot elev with residuals
par(mfrow=c(2,2))
plot(ylim=c(-400,200),xlim=c(0,2500),pch=16,cex=0.5,elevSRPD_df$wc2.1_30s_elev,elevSRPD_df$residuals_SR_PD_1km_T10,xlab="Elevation (masl)",ylab="Residuals of the PD~SR regression")
plot(ylim=c(-4,4),xlim=c(0,2500),pch=16,cex=0.5,elevPDFD_df$wc2.1_30s_elev,elevPDFD_df$residuals_PD_FD_1km_T10,xlab="Elevation (masl)",ylab="Residuals of the FDrich~PD regression")
plot(ylim=c(-4,4),xlim=c(0,2500),pch=16,cex=0.5,elevSRFD_df$wc2.1_30s_elev,elevSRFD_df$residuals_SR_FD_1km_T10,xlab="Elevation (masl)",ylab="Residuals of the FDrich~SR regression")
plot(ylim=c(-0.2,0.2),xlim=c(0,2500),pch=16,cex=0.5,elevPDFDdis_df$wc2.1_30s_elev,elevPDFDdis_df$residuals_PD_FDdis_1km_T10,xlab="Elevation (masl)",ylab="Residuals of the FDdis~PD regression")

}
      ###this means 10 plots for each residual category : 40 plots total this will be supplemental 
          ###maybe look at extremes and show deciles 1 and 10?


