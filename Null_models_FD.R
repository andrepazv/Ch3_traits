####Modified by Andrea Paz in 2021 from PD randomizations Written by Andrea Paz and David Urbina #####

rm(list=ls())
##PD randomizations for each community
library(FD) #to compute pd
library(parallel) #this code uses parallel computing
library(raster) #to generate output rasters
#=======Functions==================

#In this function m refers to number of randomizations, n to number of species in a community, sn to the names of species in the global pool 
#and traits to the trait database with only traits and rownmaes with species names
give_complete_community_better= function(m=10000, n, sn,trait) {
  #Generate community matrix
  v.names <- paste("random_community",1:m,sep="") #row names
  random_communities<-matrix(0,nrow=m,ncol=length(sn),dimnames=list(v.names,sn))
  for (k in 1:m) {
      sampling_names <- sample(sn, n)
      random_communities[k, sampling_names] <- 1
  }
  FD.result<-FD::dbFD(trait,random_communities,calc.FDiv=F)
#  pd.result<-pd(random_communities,phylo,include.root=TRUE)
  return(FD.result)
}


#===========================================================
##load and process trait database, this must march community data
#trait1<-read.csv("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/Species_traits_imputed_lilian.csv",h=T)
setwd("C:\\Users\\AnaCarnaval\\Dropbox\\Andrea_Lab\\trasnfer_819_PC\\Functional_chap")
setwd("Dropbox/Andrea_Lab/trasnfer_819_PC/Functional_chap/")
trait1<-read.csv("Species_traits_imputed_lilian.csv")
rownames(trait1)<-trait1[,1]
trait1<-trait1[,c(4,7,8)]
for (i in 1:3)
{
  trait1[,i]<-log(trait1[,i]) ## This could be updated to another mathematic transformation if needed
}

##Load community compostition matrix (each row represents one pixel and each column a species, the last column is the observed PD)
##The first column must contain the pixel numbers, these will become the rownames
communities<-read.table("communities_hylids_and_fd10k_cont.txt",h=T,row.names=1) #10k
communities<-read.table("communities_hylids_and_fd1k.txt",h=T,row.names=1)#1k
communities<-communities[,c(1:158,160)] ##leaving only species and Fdric
communities<-communities[,1:159] ##for 1 k matrix is different
names(communities)[159]<-"fdRic"##for 1k only
species_per_pixel <- unique(apply(communities[,1:length(communities[1,]) -1],1,sum))
species_names<-names(communities)[1:(length(communities[1,])-1)]

#MAC
b<-mclapply(species_per_pixel,function(x) { give_complete_community_better(n=x, sn= species_names,trait=trait1) }, mc.cores=7 )
#windows doesn t work
#b<-parLapply(species_per_pixel,function(x) { give_complete_community_better(n=x, sn= species_names,trait=traits1) }, mc.cores=16 )

save(b,file="random_hylid1k_FDrich")
#Create randomization maps
  #First compare observed PD with expected PD
communities1<-communities
communities1$sp_number<-apply(communities[,1:length(communities[1,]) -1],1,sum)
communities1$p_values<-NA
communities1$p_values_lower<-NA
communities1$p_values_higher<-NA

for (i in 1:length(species_per_pixel))
{
  subset_per_count <- subset(communities1, communities1$sp_number == species_per_pixel[i])
  match_rows<-(rownames(subset_per_count))
  observed_values<-communities1[match_rows,"fdRic"] 
  number_match<-length(observed_values)
  p_values_higher<-vector()
  p_values_lower<-vector()
  for (j in 1:number_match){
    cond<-lapply(b,function(x) x[1][[1]][1]==species_per_pixel[i])
    p_values_higher[j]<-((sum(b[unlist(cond)][[1]][[3]]>=observed_values[j]))/(10000+1))*2 
    p_values_lower[j]<-((sum(b[unlist(cond)][[1]][[3]]<=observed_values[j]))/(10000+1))*2
  }

  communities1[match_rows,"p_values_lower"]<-p_values_lower
  communities1[match_rows,"p_values_higher"]<-p_values_higher
  
}

#Chose P-value between the inferior and superior value.
a<-ifelse(communities1$p_values_lower < communities1$p_values_higher, communities1$p_values <- communities1$p_values_lower,  communities1$p_values <- communities1$p_values_higher)
communities1$p_values<-a

#Create a column indicating whether the observed FD is Higher or Lower than expected.
communities1$desviacion<-NA
communities1$desviacion[communities1$p_values==communities1$p_values_higher]<-"Superior"
communities1$desviacion[communities1$p_values==communities1$p_values_lower]<-"Inferior"

#Subsetting the data frame to obtain one per category (one higher one lower)

comunidades_sup<-subset(communities1,desviacion=="Superior")
comunidades_inf<-subset(communities1,desviacion=="Inferior")

#Create rasters and plots for all significantly different than expected pixels, all signifantly higher pixels and all signifcantly lower pixels 

#1-Create an empty raster using a base to ensure correct extent, size and resolution
r<-raster("OR_AUC_T10/AF_Scinax_alter_OR_AUC_T10.asc")
values(r)<-NA #eresa all values from base map
p_value_ras<-r #1
p_value_sup_ras<-r #2
p_value_inf_ras<-r #3




#2- Assign P-values to pixels
p_value_ras[as.integer(rownames(communities1))]<-communities1$p_values#1
p_value_sup_ras[as.integer(rownames(comunidades_sup))]<-comunidades_sup$p_values #2
p_value_inf_ras[as.integer(rownames(comunidades_inf))]<-comunidades_inf$p_values #3

#3- Save rasters to file
writeRaster(p_value_ras,"P_value_randomization_fdrich_hylids.asc")
writeRaster(p_value_sup_ras,"P_value_sup_randomization_fdrich_hylids.asc")
writeRaster(p_value_inf_ras,"P_value_inf_randomization_fdrich_hylids.asc")

#4-Plot maps in R
par(mfrow=c(1,3))
plot(p_value_ras)
plot(p_value_sup_ras)
plot(p_value_inf_ras)

#5 Reclassify rasters to only show significant pixels
  ##assuming alfa of 0.05
  ##Create reclassifying matrix. Pixels significantly higher will have a value of 1 and pixels significantly lower a value of -1
  m<-c(0,0.05,1,0.05,1,NA)
  rclmat<-matrix(m,ncol=3,byrow=TRUE)
  p_value_sup<-reclassify(p_value_sup_ras, rclmat, include.lowest=T)
   m<-c(0,0.05,-1,0.05,1,NA)
  rclmat<-matrix(m,ncol=3,byrow=TRUE)
  p_value_inf<-reclassify(p_value_inf_ras, rclmat,include.lowest=T)
  
 significant_pixels<-mosaic(p_value_sup,p_value_inf, fun=mean)
writeRaster(significant_pixels,"Significant_pixels.asc")
  
  
###if using Fdispersion

####Modified by Andrea Paz in 2021 from PD randomizations Written by Andrea Paz and David Urbina #####

rm(list=ls())
##PD randomizations for each community
library(FD) #to compute pd
library(parallel) #this code uses parallel computing
library(raster) #to generate output rasters
#=======Functions==================

#In this function m refers to number of randomizations, n to number of species in a community, sn to the names of species in the global pool 
#and traits to the trait database with only traits and rownmaes with species names
give_complete_community_better= function(m=10000, n, sn,trait) {
  #Generate community matrix
  v.names <- paste("random_community",1:m,sep="") #row names
  random_communities<-matrix(0,nrow=m,ncol=length(sn),dimnames=list(v.names,sn))
  for (k in 1:m) {
    sampling_names <- sample(sn, n)
    random_communities[k, sampling_names] <- 1
  }
  FD.result<-FD::dbFD(trait,random_communities,calc.FDiv=F,calc.FRic = F,corr="cailliez")
  #  pd.result<-pd(random_communities,phylo,include.root=TRUE)
  return(FD.result)
}


#===========================================================
##load and process trait database, this must march community data
#trait1<-read.csv("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/Species_traits_imputed_lilian.csv",h=T)
setwd("C:\\Users\\AnaCarnaval\\Dropbox\\Andrea_Lab\\trasnfer_819_PC\\Functional_chap")
setwd("Dropbox/Andrea_Lab/trasnfer_819_PC/Functional_chap/")
trait1<-read.csv("Species_traits_imputed_lilian.csv")
rownames(trait1)<-trait1[,1]
trait1<-trait1[,c(3,4,7,8)]
for (i in 2:4)
{
  trait1[,i]<-log(trait1[,i]) ## This could be updated to another mathematic transformation if needed
}

##Load community compostition matrix (each row represents one pixel and each column a species, the last column is the observed PD)
##The first column must contain the pixel numbers, these will become the rownames
communities<-read.table("communities_hylids_and_fd1k_cat.txt",h=T,row.names=1)
communities<-communities[,c(1:158,159)] ##leaving only species and FdDisC Always check names first
species_per_pixel <- unique(apply(communities[,1:length(communities[1,]) -1],1,sum))
species_names<-names(communities)[1:(length(communities[1,])-1)]

#MAC
b<-mclapply(species_per_pixel,function(x) { give_complete_community_better(n=x, sn= species_names,trait=trait1) }, mc.cores=7 )
#windows
b<-parLapply(species_per_pixel,function(x) { give_complete_community_better(n=x, sn= species_names,trait=traits1) }, mc.cores=16 )

save(b,file="random_hylid1k_cat")
#Create randomization maps
#First compare observed PD with expected PD
communities1<-communities
communities1$sp_number<-apply(communities[,1:length(communities[1,]) -1],1,sum)
communities1$p_values<-NA
communities1$p_values_lower<-NA
communities1$p_values_higher<-NA

for (i in 1:length(species_per_pixel))
{
  subset_per_count <- subset(communities1, communities1$sp_number == species_per_pixel[i])
  match_rows<-(rownames(subset_per_count))
  observed_values<-communities1[match_rows,"fdDisC"]
  number_match<-length(observed_values)
  p_values_higher<-vector()
  p_values_lower<-vector()
  for (j in 1:number_match){
    cond<-lapply(b,function(x) x[1][[1]][1]==species_per_pixel[i])
    p_values_higher[j]<-((sum(b[unlist(cond)][[1]][[4]]>=observed_values[j]))/(10000+1))*2 
    p_values_lower[j]<-((sum(b[unlist(cond)][[1]][[4]]<=observed_values[j]))/(10000+1))*2
  }
  
  communities1[match_rows,"p_values_lower"]<-p_values_lower
  communities1[match_rows,"p_values_higher"]<-p_values_higher
  
}

#Chose P-value between the inferior and superior value.
a<-ifelse(communities1$p_values_lower < communities1$p_values_higher, communities1$p_values <- communities1$p_values_lower,  communities1$p_values <- communities1$p_values_higher)
communities1$p_values<-a

#Create a column indicating whether the observed FD is Higher or Lower than expected.
communities1$desviacion<-NA
communities1$desviacion[communities1$p_values==communities1$p_values_higher]<-"Superior"
communities1$desviacion[communities1$p_values==communities1$p_values_lower]<-"Inferior"

#Subsetting the data frame to obtain one per category (one higher one lower)

comunidades_sup<-subset(communities1,desviacion=="Superior")
comunidades_inf<-subset(communities1,desviacion=="Inferior")

#Create rasters and plots for all significantly different than expected pixels, all signifantly higher pixels and all signifcantly lower pixels 

#1-Create an empty raster using a base to ensure correct extent, size and resolution
r<-raster("OR_AUC_T10/AF_Hypsiboas_leptolineatus_OR_AUC_T10.asc")
values(r)<-NA #eresa all values from base map
p_value_ras<-r #1
p_value_sup_ras<-r #2
p_value_inf_ras<-r #3




#2- Assign P-values to pixels
p_value_ras[as.integer(rownames(communities1))]<-communities1$p_values#1
p_value_sup_ras[as.integer(rownames(comunidades_sup))]<-comunidades_sup$p_values #2
p_value_inf_ras[as.integer(rownames(comunidades_inf))]<-comunidades_inf$p_values #3

#3- Save rasters to file
writeRaster(p_value_ras,"P_value_randomization_fdis_hylids1k.asc")
writeRaster(p_value_sup_ras,"P_value_sup_randomization_fdis_hylids1k.asc")
writeRaster(p_value_inf_ras,"P_value_inf_randomization_fddis_hylids1k.asc")

#4-Plot maps in R
par(mfrow=c(1,3))
plot(p_value_ras)
plot(p_value_sup_ras)
plot(p_value_inf_ras)

#5 Reclassify rasters to only show significant pixels
##assuming alfa of 0.05
##Create reclassifying matrix. Pixels significantly higher will have a value of 1 and pixels significantly lower a value of -1
m<-c(0,0.05,1,0.05,2,NA)
rclmat<-matrix(m,ncol=3,byrow=TRUE)
p_value_sup<-reclassify(p_value_sup_ras, rclmat, include.lowest=T)
m<-c(0,0.05,-1,0.05,2,NA)
rclmat<-matrix(m,ncol=3,byrow=TRUE)
p_value_inf<-reclassify(p_value_inf_ras, rclmat,include.lowest=T)

significant_pixels<-mosaic(p_value_sup,p_value_inf, fun=mean)
writeRaster(significant_pixels,"Significant_pixels_FDdis1k.asc")


plot(significant_pixels)








