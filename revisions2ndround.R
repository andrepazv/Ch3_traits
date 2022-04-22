########################################################
##########Analysis reviews ecography April 2022#########
########################################################
library(rgdal)
library(raster)
library(ggplot2)

###compare residuals per region in the AF#####

##load vegeation layer from Figure 1 IBGE
vegetation_AF<-readOGR("/Users/andreapaz/Dropbox/SIG&Modelamiento/datos_sig/Por_paises/Brasil/Cobertura vegetal do Brasil (Escala_1-5000000)/veg_AF_only.shp")

#load asci with residuals
setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/")
res_PD_SR<-raster("gam_residuals_SR_PD_1km_T10.asc")
res_FD_SR<-raster("gam_residuals_SR_FD_1km_T10.asc")
res_FDis_PD<-raster("gam_residuals_PD_FDdis_1km_T10.asc")

##cut per region and plot

##keep only the attribute of interest "Tipo2" and its unique values to cut#
regions<-unique(vegetation_AF$TIPO2)
##only a few are of concern
regions<-regions[c(4,5,8,12)]
res_PD_SR_reg<-list()
res_FD_SR_reg<-list()
res_FDis_PD_reg<-list()
for (i in 1:length(regions)){
  reg_veg<-subset(vegetation_AF,TIPO2==regions[i])
  res_PD_SR_reg[[i]]<-crop(res_PD_SR,reg_veg)
  res_FD_SR_reg[[i]]<-crop(res_FD_SR,reg_veg)
  res_FDis_PD_reg[[i]]<-crop(res_FDis_PD,reg_veg)
}
res_PD_SR_reg<-lapply(res_PD_SR_reg,values)
df<-data.frame(List=unlist(res_PD_SR_reg),Region=rep(c("Contact","Evergreen Forest","Deciduous Forest","Araucaria Forest"),times=c(length(res_PD_SR_reg[[1]]),length(res_PD_SR_reg[[2]]),length(res_PD_SR_reg[[3]]),length(res_PD_SR_reg[[4]]))))
ggplot(df,aes(Region,List))+geom_boxplot(outlier.shape=NA)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(y="PD~TD residuals") + ylim(-150,150)

res_FD_SR_reg<-lapply(res_FD_SR_reg,values)
df<-data.frame(List=unlist(res_FD_SR_reg),Region=rep(c("Contact","Evergreen Forest","Deciduous Forest","Araucaria Forest"),times=c(length(res_FD_SR_reg[[1]]),length(res_FD_SR_reg[[2]]),length(res_FD_SR_reg[[3]]),length(res_FD_SR_reg[[4]]))))
ggplot(df,aes(Region,List))+geom_boxplot(outlier.shape=NA)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(y="FDrich~TD residuals") + ylim(-3,3)

res_FDis_PD_reg<-lapply(res_FDis_PD_reg,values)
df<-data.frame(List=unlist(res_FDis_PD_reg),Region=rep(c("Contact","Evergreen Forest","Deciduous Forest","Araucaria Forest"),times=c(length(res_FDis_SR_reg[[1]]),length(res_FDis_PD_reg[[2]]),length(res_FDis_PD_reg[[3]]),length(res_FDis_PD_reg[[4]]))))
ggplot(df,aes(Region,List))+geom_boxplot(outlier.shape=NA)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
labs(y="FDis~PD residuals")+ ylim(-0.1,0.1)


#######################################
######Generate all GAMM models#########
#######################################

#### Number 1 would be just to model the non linear relationship between diversity metrics

library(raster)
library(mgcv)
###Load maps of different diversity metrics

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
srpd_df$grid<-1:length(srpd_df$hylids_PD)
srpd_df1<-na.omit(srpd_df)
srfd_df<-as.data.frame(stack_SR_FD)
srfd_df$grid<-1:length(srfd_df$richness_hylids_complete_0_0083)
srfd_df1<-na.omit(srfd_df)
pdfd_df<-as.data.frame(stack_PD_FD)
pdfd_df$grid<-1:length(pdfd_df$hylids_PD)
pdfd_df1<-na.omit(pdfd_df)
pdfddis_df<-as.data.frame(stack_PD_FDdis)
pdfddis_df$grid<-1:length(pdfddis_df$hylid_functional_dispersion_cat1k)
pdfddis_df1<-na.omit(pdfddis_df)
#Generate regression models
#LINEAR
linear_model_sr_pd<-lm(srpd_df)
linear_model_sr_fd<-lm(srfd_df)
linear_model_pd_fd<-lm(pdfd_df)
linear_model_pd_fddis<-lm(pdfddis_df)

gam_model_sr_pd<-gam(hylids_PD~te(richness_hylids_complete_0_0083),data=srpd_df1)
gam_model_sr_fd<-gam(hylid_functional_richness_cont~te(richness_hylids_complete_0_0083),data=srfd_df1)
gam_model_pd_fd<-gam(hylid_functional_richness_cont~te(hylids_PD),data=pdfd_df1)
gam_model_pd_fddis<-gam(hylid_functional_dispersion_cat1k ~te(hylids_PD),data=pdfddis_df1)

library(data.table)

preds <- predict(gam_model_sr_pd,se.fit=TRUE)
ggplot()+
  geom_point(data=srpd_df1,aes(x=richness_hylids_complete_0_0083,y=hylids_PD))+
  geom_line( aes(x=srpd_df1$richness_hylids_complete_0_0083, y=preds$fit), size=1, col="blue")+
#  geom_line(aes(x=linear_model_sr_pd$model$richness_hylids_complete_0_0083,y=linear_model_sr_pd$fitted.values),col="red")+
  theme_bw()+labs(x="Taxonomic diversity",y="Phylogenetic diversity")+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16))

preds <- predict(gam_model_sr_fd,se.fit=TRUE)
ggplot()+
  geom_point(data=srfd_df1,aes(x=richness_hylids_complete_0_0083,y=hylid_functional_richness_cont))+
  geom_line(aes(x=srfd_df1$richness_hylids_complete_0_0083, y=preds$fit), size=1, col="blue")+
#  geom_line(aes(x=linear_model_sr_fd$model$richness_hylids_complete_0_0083,y=linear_model_sr_fd$fitted.values),col="red")+
  theme_bw()+labs(x="Taxonomic diversity",y="Functional richness")+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16))

preds <- predict(gam_model_pd_fd,se.fit=TRUE)
ggplot()+
  geom_point(data=pdfd_df1,aes(x=hylids_PD,y=hylid_functional_richness_cont))+
  geom_line(aes(x=pdfd_df1$hylids_PD, y=preds$fit), size=1, col="blue")+
#  geom_line(aes(x=linear_model_pd_fd$model$hylids_PD,y=linear_model_pd_fd$fitted.values),col="red")+
  theme_bw()+labs(x="Phylogenetic diversity",y="Functional richness")+
theme(axis.text=element_text(size=16),axis.title=element_text(size=16))


preds <- predict(gam_model_pd_fddis,se.fit=TRUE)
ggplot()+
  geom_point(data=pdfddis_df1,aes(x=hylids_PD,y=hylid_functional_dispersion_cat1k))+
  geom_line( aes(x=pdfddis_df1$hylids_PD, y=preds$fit), size=1, col="blue")+

#  geom_line(aes(x=linear_model_pd_fddis$model$hylids_PD,y=linear_model_pd_fddis$fitted.values),col="red")+
  theme_bw()+labs(x="Phylogenetic diversity",y="Functional dispersion")+
theme(axis.text=element_text(size=16),axis.title=element_text(size=16))


#Compute and save residuals from linear models
linear_residuals_sr_pd<-resid(linear_model_sr_pd)
linear_residuals_sr_fd<-resid(linear_model_sr_fd)
linear_residuals_pd_fd<-resid(linear_model_pd_fd)
linear_residuals_pd_fddis<-resid(linear_model_pd_fddis)

#Compute and save residuals from gam models
gam_residuals_sr_pd<-resid(gam_model_sr_pd)
gam_residuals_sr_fd<-resid(gam_model_sr_fd)
gam_residuals_pd_fd<-resid(gam_model_pd_fd)
gam_residuals_pd_fddis<-resid(gam_model_pd_fddis)

###ADD linear residuals to dataframe 
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

###ADD GAM residuals to data frame 

#Add residuals to dataframe
srpd_df$GAM_residuals<-rep(NA,length(srpd_df[,1]))
srpd_df$GAM_residuals[srpd_df1$grid]<-gam_residuals_sr_pd
#add residuals of fd to data frame
srfd_df$GAM_residuals<-rep(NA,length(srfd_df[,1]))
srfd_df$GAM_residuals[srfd_df1$grid]<-gam_residuals_sr_fd
#add residuals of pd/fd to data frame
pdfd_df$GAM_residuals<-rep(NA,length(pdfd_df[,1]))
pdfd_df$GAM_residuals[pdfd_df1$grid]<-gam_residuals_pd_fd
#add residuals of pd/fddis to data frame ##check
pdfddis_df$GAM_residuals<-rep(NA,length(pdfddis_df[,1]))
pdfddis_df$GAM_residuals[pdfddis_df1$grid]<-gam_residuals_pd_fddis


#Generate rasters with residual values for PD
gam_residuals_sr_pd_raster<-raster(pd_raster)
values(gam_residuals_sr_pd_raster)<-NA
values(gam_residuals_sr_pd_raster)<- srpd_df$GAM_residuals
name_group<-"gam_residuals_SR_PD_1km_T10"
writeRaster(gam_residuals_sr_pd_raster,name_group,"ascii")

gam_residuals_sr_fd_raster<-raster(pd_raster)
values(gam_residuals_sr_fd_raster)<-NA
values(gam_residuals_sr_fd_raster)<- srfd_df$GAM_residuals
name_group<-"gam_residuals_SR_FD_1km_T10"
writeRaster(gam_residuals_sr_fd_raster,name_group,"ascii")

gam_residuals_pd_fd_raster<-raster(pd_raster)
values(gam_residuals_pd_fd_raster)<-NA
values(gam_residuals_pd_fd_raster)<- pdfd_df$GAM_residuals
name_group<-"gam_residuals_PD_FD_1km_T10"
writeRaster(gam_residuals_pd_fd_raster,name_group,"ascii")


gam_residuals_pd_fddis_raster<-raster(pd_raster)
values(gam_residuals_pd_fddis_raster)<-NA
values(gam_residuals_pd_fddis_raster)<- pdfddis_df$GAM_residuals
name_group<-"gam_residuals_PD_FDdis_1km_T10"
writeRaster(gam_residuals_pd_fddis_raster,name_group,"ascii")

plot(gam_residuals_pd_fddis_raster)


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

gam_residuals_sr_pd_raster
gam_residuals_sr_fd_raster
gam_residuals_pd_fd_raster
gam_residuals_pd_fddis_raster


elev<-crop(elev,sr_raster)
extent(elev)<-extent(sr_raster)
###MAYBE
stack_elev_SR<-stack(sr_raster,elev,lat)
stack_elev_PD<-stack(pd_raster,elev,lat)
stack_elev_FD<-stack(fdRic_raster,elev,lat)
stack_elev_FDdis<-stack(fdRic_raster,elev,lat)
elevSR_df<-as.data.frame(stack_elev_SR)
elevPD_df<-as.data.frame(stack_elev_PD)
elevFD_df<-as.data.frame(stack_elev_FD)
elevFDdis_df<-as.data.frame(stack_elev_FDdis)
###residuals
stack_elev_SRPD<-stack(gam_residuals_sr_pd_raster,elev,lat)
stack_elev_SRFD<-stack(gam_residuals_sr_fd_raster,elev,lat)
stack_elev_PDFD<-stack(gam_residuals_pd_fd_raster,elev,lat)
stack_elev_PDFDdis<-stack(gam_residuals_pd_fddis_raster,elev,lat)
elevSRPD_df<-as.data.frame(stack_elev_SRPD)
elevSRFD_df<-as.data.frame(stack_elev_SRFD)
elevPDFD_df<-as.data.frame(stack_elev_PDFD)
elevPDFDdis_df<-as.data.frame(stack_elev_PDFDdis)
###SET UP GAM for elev * lat
  gam_model_srpdres_elevLat<-gam(layer~ti(wc2.1_30s_elev)+ti(AF_latitude)+ti(wc2.1_30s_elev,AF_latitude),data= elevSRPD_df)
  gam_model_srfdres_elevLat<-gam(layer~ti(wc2.1_30s_elev)+ti(AF_latitude)+ti(wc2.1_30s_elev,AF_latitude),data= elevSRFD_df)
   gam_model_pdfdres_elevLat<-gam(layer~ti(wc2.1_30s_elev)+ti(AF_latitude)+ti(wc2.1_30s_elev,AF_latitude),data= elevPDFD_df)
 gam_model_pdfddisres_elevLat<-gam(layer~ti(wc2.1_30s_elev)+ti(AF_latitude)+ti(wc2.1_30s_elev,AF_latitude),data= elevPDFDdis_df)
  


####ADD here to build Figure 6 , residuals * elevation per altitudinal band
 
 ##then within a loop crop the raster and plot for each quantile the residuals
 elev1<-elev
 lats<-list.files(pattern="lat")[2:11]
 
 
 for (i in 1:length(lats)){
   decile<-raster::raster(lats[i])
   pdSr_res<-crop(gam_residuals_sr_pd_raster,decile)
   fdricSr_res<-crop(gam_residuals_sr_fd_raster,decile)
   fdricPd_res<-crop(gam_residuals_pd_fd_raster,decile)
   fddisPd_res<-crop(gam_residuals_pd_fddis_raster,decile)
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
   par(mfrow=c(1,3))
   plot(main=names(lats[i]),ylim=c(-4,4),xlim=c(0,2500),pch=16,cex=0.5,elevSRFD_df$wc2.1_30s_elev,elevSRFD_df$layer,xlab="Elevation (masl)",ylab="Residuals of the FDrich~SR gam model")
   # plot(ylim=c(-4,4),xlim=c(0,2500),pch=16,cex=0.5,elevPDFD_df$wc2.1_30s_elev,elevPDFD_df$layer,xlab="Elevation (masl)",ylab="Residuals of the FDrich~PD gam model")
   plot(ylim=c(-400,200),xlim=c(0,2500),pch=16,cex=0.5,elevSRPD_df$wc2.1_30s_elev,elevSRPD_df$layer,xlab="Elevation (masl)",ylab="Residuals of the PD~SR gam model")
   plot(ylim=c(-0.2,0.2),xlim=c(0,2500),pch=16,cex=0.5,elevPDFDdis_df$wc2.1_30s_elev,elevPDFDdis_df$layer,xlab="Elevation (masl)",ylab="Residuals of the FDdis~PD gam model")
   
 }
