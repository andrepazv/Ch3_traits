###obtain the bioclimatic predictors
###WARNING this code needs numbers of predictor and repsonse veriables at several instances, they must be manually changed

#library(reticulate)
##use_python("/usr/local/bin/python3") this might be needed in case the python path is not recognized automatically by computer
setwd("~/")
###MODIS data

#source_python("~/Dropbox/**Tesis_PHD/Near_real_time_model/Python_GEE_scripts/GetMODIS.py")
source_python("~/Andrea/Near_real_time_model/GetMODIS.py")

###Chirps data
#source_python("~/Dropbox/**Tesis_PHD/Near_real_time_model/Python_GEE_scripts/GetChirps.py")
source_python("~/Andrea/Near_real_time_model/GetChirps.py")

#######################################################################################
###WARNING: after sourcing python scripts must stop until all layers are generated! 
###If not then not all layers will be dowloaded prior to analyses
#######################################################################################


##Get diversity data (response variables)
library(raster)
#setwd("~/Dropbox/**Tesis_PHD/Near_real_time_model/Layers_yearly/")
setwd("~/Andrea/Near_real_time_model/Layers_yearly/")
div.files <- list.files("Diversity_maps/",
                        pattern = "asc$", full.names = T) ##try to set a CRS for all these diversity layers!
test_ras<-raster(div.files[1])

##obtain all layers that are in google dirve after python GEE processing. 
#This can fail randomly because of an unresolved bug in the googledrive package (for now on computer)
#library(googledrive)
#ayerList<-drive_ls("~/Layers")
#LayerList<-as.vector(LayerList$name)
#All_layers<-lapply(LayerList,drive_download,overwrite=T)

##if already in computer then:
LayerList<-list.files(pattern=".tif")

All_layers<-lapply(LayerList,raster)
All_layers<-lapply(All_layers,resample,y=test_ras)
##get layer names to select 
All_names<-unlist(lapply(All_layers,names))
Mean_layers<-stack(All_layers)[[which(grepl("mean",All_names))]]
Max_layers<-stack(All_layers)[[which(grepl("max",All_names))]]
##here add a replacement for NA
Min_layers<-stack(All_layers)[[which(grepl("min",All_names))]]
Monthly_precip<-stack(All_layers)[[which(grepl("total",All_names))]]
###These are already one per year !
Precip_max<-stack(All_layers)[[which(grepl("Max_Precip",All_names))]] 
Precip_min<-stack(All_layers)[[which(grepl("Min_Precip",All_names))]]
###Select years to use as training 2002-2014 (as trial run) and so predict for 2015,2016,2017 and 2018 if available
##If changing years then just change numbers in train/test vectors
trainV <- c(2000:2014)
testV<-c(2015:2019)
##Separate each stack of variables in train vs. test years
##do check if this is catching only precip or other things..

Monthly_precip_train<-Monthly_precip[[which(grepl(paste(trainV,collapse="|"),names(Monthly_precip)))]] 
Monthly_precip_test<-Monthly_precip[[which(grepl(paste(testV,collapse="|"),names(Monthly_precip)))]] 

Max_precip_train<-Precip_max[[which(grepl(paste(trainV,collapse="|"),names(Precip_max)))]] 
Max_precip_test<-Precip_max[[which(grepl(paste(testV,collapse="|"),names(Precip_max)))]] 

Min_precip_train<-Precip_min[[which(grepl(paste(trainV,collapse="|"),names(Precip_min)))]] 
Min_precip_test<-Precip_min[[which(grepl(paste(testV,collapse="|"),names(Precip_min)))]] 


Mean_layers_train<- Mean_layers[[which(grepl(paste(trainV,collapse="|"),names(Mean_layers)))]]
Mean_layers_test<- Mean_layers[[which(grepl(paste(testV,collapse="|"),names(Mean_layers)))]]

Max_layers_train<- Max_layers[[which(grepl(paste(trainV,collapse="|"),names(Max_layers)))]]
Max_layers_test<- Max_layers[[which(grepl(paste(testV,collapse="|"),names(Max_layers)))]]

Min_layers_train<- Min_layers[[which(grepl(paste(trainV,collapse="|"),names(Min_layers)))]]
Min_layers_test<- Min_layers[[which(grepl(paste(testV,collapse="|"),names(Min_layers)))]]

##For train years get the mean value of each variable (to get only one variable per month per category e.g one total precip for january, one for february etc.)

for (i in 1:12){
  LayerTemp<-mean(Mean_layers_train[[which(grepl(paste("_",i,"$",sep=""),names(Mean_layers_train)))]],na.rm=T)
  names(LayerTemp)<-paste("mean",i,sep="_")
  assign(paste("mean_train",i,sep="_"),LayerTemp)
  LayerTemp1<-mean(Max_layers_train[[which(grepl(paste("_",i,"$",sep=""),names(Max_layers_train)))]],na.rm=T)
  names(LayerTemp1)<-paste("max",i,sep="_")
  assign(paste("max_train",i,sep="_"),LayerTemp1)
  LayerTemp2<-mean(Min_layers_train[[which(grepl(paste("_",i,"$",sep=""),names(Min_layers_train)))]],na.rm=T)
  names(LayerTemp2)<-paste("min",i,sep="_")
  assign(paste("min_train",i,sep="_"),LayerTemp2)
  ##precip layers
  LayerTemp3<-mean(Monthly_precip_train[[which(grepl(paste("_",i,"$",sep=""),names(Monthly_precip_train)))]],na.rm=T)
  names(LayerTemp3)<-paste("total",i,sep="_")
  assign(paste("total_train",i,sep="_"),LayerTemp3 )
}
LayerTemp4<-mean(Max_precip_train,na.rm=T)
names(LayerTemp4)<-"Max_Precip"
max_precip_train<-LayerTemp4

LayerTemp5<-mean(Min_precip_train,na.rm=T)
names(LayerTemp5)<-"Min_Precip"
min_precip_train<-LayerTemp5 
##Put all train variables in a single stack
All_train<-stack(lapply(ls(pattern="mean_train|min_train|max_train|total_train|max_precip_train|min_precip_train"),get))

###Get test variables and separate per test year
for (i in testV){
  LayersTemp<-stack(Mean_layers_test[[which(grepl(i,names(Mean_layers_test)))]])
  LayersTemp1<-stack(Max_layers_test[[which(grepl(i,names(Max_layers_test)))]])
  LayersTemp2<-stack(Min_layers_test[[which(grepl(i,names(Min_layers_test)))]])
  LayersTemp3<-stack(Monthly_precip_test[[which(grepl(i,names(Monthly_precip_test)))]])
  LayersTemp4<-stack(Max_precip_test[[which(grepl(i,names(Max_precip_test)))]])
  LayersTemp5<-stack(Min_precip_test[[which(grepl(i,names(Min_precip_test)))]])
  if(nlayers(LayersTemp)>0){
    names(LayersTemp)<-gsub(paste(i,"_",sep=""),"",names(LayersTemp))
    assign(paste("mean_test",i,sep="_"),LayersTemp)
  }
  
  if(nlayers(LayersTemp1)>0){
    LayersTemp1<-subset(LayersTemp1,sort(names(LayersTemp1)))
    names(LayersTemp1)<-gsub(paste(i,"_",sep=""),"",names(LayersTemp1))
    assign(paste("max_test",i,sep="_"),LayersTemp1)
  }
  if(nlayers(LayersTemp2)>0){
    LayersTemp2<-subset(LayersTemp2,sort(names(LayersTemp2)))
    names(LayersTemp2)<-gsub(paste(i,"_",sep=""),"",names(LayersTemp2))
    assign(paste("min_test",i,sep="_"),LayersTemp2)
  }
  if(nlayers(LayersTemp3)>0){
    LayersTemp3<-subset(LayersTemp3,sort(names(LayersTemp3)))
    names(LayersTemp3)<-gsub(paste(i,"_",sep=""),"",names(LayersTemp3))
    assign(paste("total_test",i,sep="_"),LayersTemp3)
  }
  if(nlayers(LayersTemp4)>0){
    names(LayersTemp4)<-gsub(paste("_",i,sep=""),"",names(LayersTemp4))
    assign(paste("max_precip_test",i,sep="_"),LayersTemp4)
  }
  if(nlayers(LayersTemp5)>0){
    names(LayersTemp5)<-gsub(paste("_",i,sep=""),"",names(LayersTemp5))
    assign(paste("min_precip_test",i,sep="_"),LayersTemp5)
  }
  
  
  
  
}


###Calculate biovars to make predictions with those. 
library(dismo) #this package to calculate bioclimatic variables

###At this point check that numbers correspond to the right variables
precip<-stack(All_train[[39:50]])
precip<-stack(precip[[c(1,5:12,2:4)]]) ##to organize by month jan-Dec
tmax<-stack(All_train[[2:13]])
tmax<-stack(tmax[[c(1,5:12,2:4)]])
tmin<-stack(All_train[[27:38]])
tmin<-stack(tmin[[c(1,5:12,2:4)]])
bioclimatic<-biovars(precip,tmin,tmax)
library(usdm)
viRes<-vifstep(bioclimatic)
bioclimatic<-subset(bioclimatic,as.character(viRes@results$Variables))

# Loading libraries
library(raster)
library(caret)
library(tools)
library(tidyr)
library(corrplot)
library(rasterVis)
library(nnet)
library(caTools)
library(caretEnsemble)
library(CAST)


# Loading input data

# Change this path to where repsonse variables are
#setwd("")


####################################
#### 1) Custom functions for later use

# rs.as.df: Takes a list of raster filenames and reads their values
#           straight into dataframe columns

rs.as.df <- function(raslist){
  if (length(raslist) < 1){
    print("Oops! File list is empty!")
    break
  }
  nomes <- file_path_sans_ext(basename(raslist)) # get file names sans extensions and path
  rs.list <- lapply(raslist,raster)
  rs.ncell <- lapply(rs.list,ncell)
  if (length(unique(rs.ncell)) != 1){ # if raster number of cells are different, stop
    print("Ooops! Rasters must have the same number of cells!")
  } else { # All rasters have the same number of cells, go forward
    rs.vals <- lapply(rs.list,values) # read values for each raster in rlist as list item
    names(rs.vals) <- nomes # assigns clean file names to list items
    rs.df <- as.data.frame(rs.vals) # convert list to dataframe
    ###Add this set if doing spatial blocks
    coords<-coordinates(rs.list[[1]]) ##this is only needed if spatial blocks
    coords<-as.data.frame(coords)##this is only needed if spatial blocks
    rs.df$long<-coords$x##this is only needed if spatial blocks
    rs.df$lat<-coords$y##this is only needed if spatial blocks
    return(rs.df)
  }
}

###Function borrowed and modified from ENMeval
get.block <- function(occs){
  rownames(occs) <- 1:nrow(occs)
  # SPLIT occs POINTS INTO FOUR SPATIAL GROUPS
  noccs <- nrow(occs)
  n1 <- ceiling(nrow(occs)/2)
  n2 <- floor(nrow(occs)/2)
  n3 <- ceiling(n1/2)
  n4 <- ceiling(n2/2)
  grpA <- occs[order(occs[, 2]),][1:n1,]
  grpB <- occs[rev(order(occs[, 2])),][1:n2,]
  grp1 <- grpA[order(grpA[, 1]),][1:(n3),]
  grp2 <- grpA[!rownames(grpA)%in%rownames(grp1),]
  grp3 <- grpB[order(grpB[, 1]),][1:(n4),]
  grp4 <- grpB[!rownames(grpB)%in%rownames(grp3),]
  r <- data.frame()
  if (nrow(grp1) > 0) grp1$grp <- 1; r <- rbind(r, grp1)
  if (nrow(grp2) > 0) grp2$grp <- 2; r <- rbind(r, grp2)
  if (nrow(grp3) > 0) grp3$grp <- 3; r <- rbind(r, grp3)
  if (nrow(grp4) > 0) grp4$grp <- 4; r <- rbind(r, grp4)
  occ.grp <- r[order(as.numeric(rownames(r))),]$grp
  return(occ.grp)
}
##################################################
# 2) Reading in dependent variables and predictors

# Get list of diversity and predictor filenames in specified path,
# ending in 'asc', with full path names
div.files <- list.files("Diversity_atlantic/",
                        pattern = "asc$", full.names = T)


div.df <- rs.as.df(div.files)  ## if doing spatial blocks this will have 2 extra cols at the end long and lat
#make the stack of predictors a data frame
#pred.df <- as.data.frame(All_train) ##with variables original
##with bioclimatic variables
pred.df <- as.data.frame(bioclimatic)

# Lets save the names of responses and predictors, we'll need it later
pred.names <- names(pred.df)
div.names <- names(div.df)
### to add if spatial blocking
div.names<-div.names[-c(length(div.names),length(div.names)-1)]
#pred.names <- pred.names[-c(length(pred.names),length(pred.names)-1)]

# Bind responses and predictors together
input.df <- cbind(div.df,pred.df)

##############################
# 3) Model fitting

# Are there NA values? YES.
#summary(input.df)

# Lets keep only the full cases per group 
for (i in 1:length(div.names)){
  ##create data frame with just the predicted value in this iteration and remove all NA in both predictors and response
  input.df.na<-input.df%>%drop_na(i)
  input.df.na<-input.df.na%>%drop_na(4:12) ##this assumes 1 response (plus lat long)  and 9 predictors else must change numbers
  
  ###for spatial blocking lets just keep the coordinates in separate df. They are there twice once for predictors once for responses
  occs<-input.df.na[,which(names(input.df.na)=="long"|names(input.df.na)=="lat")]
  ###keep a single set
  occs<-occs[,c(1:2)]
  ###But before must remove lat long from the df. 
  input.df.na<-input.df.na[,-which(names(input.df.na)=="long"|names(input.df.na)=="lat")]
  # Lets standardize all predictor variables, as ML methods wil tend to behave better
  # We'll use the preProcess function of the 'caret'package
  # Here is where having a list of variable names is useful
  pre.pred <- preProcess(input.df.na[,pred.names], methods=c('scale','center'))
  input.df.na[,pred.names] <- predict(pre.pred, input.df.na[,pred.names])
  ## Now we can build training and test samples. 
  #train.df <- input.df.na[train.part,]
  #test.df <- input.df.na[-train.part,]
  ####An alternative to this is to use block partitioning and then test with those, 
  ###this means creating several partitions and running the model several times
  ####for this test run lets do 4 partitions thus running 4 times with 3 parts for training. 
  ####will do for only some groups to get fast results and then run the whole thing.
  ###Will don insde loop of 4.
  
  
  ## Always good to start with exploratory graphical analysis.
  ## Lets evaluate big_pd
  # We can use the response and predictor names to build plotting and training formulas
  # and passing them to other functions
  
  # formula for 'bignonias_pd' ##repeated down maybe erase there
  form <- as.formula(paste(div.names[i], # name of response
                           "~", # as function of
                           paste(pred.names, collapse=" + "))) # all predictors
  
  plot.df <- gather(input.df.na[,c(i,2:length(input.df.na))], key = predictor, value = value, -1)
  
  ggplot(plot.df,aes(value, get(names(input.df.na)[i]))) + geom_point() + facet_wrap(~predictor)
  # ALl correlations for all variables
  #corrplot(cor(input.df.na),tl.cex = 0.7)
  
  # A good strategy for optimizing ML models is repeated k-fold
  # With this approach, we can split the data into subsets (folds),
  # and then take turns using the folds for training and validation
  # Repeated folds means we do many random splits and use all of them in turn
  # The parameterization that gives the best model parameters is then chosen
  # This gives us more robust models, but takes more time to run
  
  blocks<- get.block(occs)
  index<-list()
  indexOut<-list()
  for (h in 1:4){
    indexOut[[h]]<-which(blocks==h)
    index[[h]]<-which(blocks!=h)
  }
  indices<-list(index=index,indexOut=indexOut)
  # Lets use 3 repetitions of 10 folds
  # fit.control <- trainControl(method='repeatedcv',number=10,repeats=3,savePredictions="final")
  ##Now using blocks
  set.seed(1979) # Lets make results reproducible
  #indices <- CreateSpacetimeFolds(trainDat,spacevar = "SOURCEID", k=4)
  fit.control <- trainControl(method='cv',index=indices$index,indexOut=indices$indexOut,savePredictions="final")
  
  
  # 'caret' is just a wrapper for the many ML algorithms implemented in R,
  # with the advantage that we can use the same syntax for any model
  
  # Example: fitting a Random Forests model to predict 
  # We can use the response and predictor names to build training formulas
  # and passing them to caret, so we can automate model building later
  
  
  # We build the formula programatically as it would help to automate over all
  # diversity variables later
  
  form <- as.formula(paste(div.names[i], # name of response
                           "~", # as function of
                           paste(pred.names, collapse=" + "))) # all predictors
  
  ###Run only 'rf' to do the spatial crossvalidation
  # m.bigRF <- train(form, data=input.df.na, trControl=fit.control, method="rf",importance=T)
  # m.bigRF
  
  ###using CAST: https://cran.r-project.org/web/packages/CAST/vignettes/CAST-intro.html
  #CAST’s forward feature selection (ffs) selects variables that make sense in view to the selected CV method and excludes those which are counterproductive in view to the selected CV method. When we use LLO as CV method, ffs selects variables that lead in combination to the highest LLO performance (i.e. the best spatial model). All variables that have no spatial meaning or are even counterproductive won’t improve or even reduce the LLO performance and are therefore excluded from the model by the ffs.
  
  # m.bigRF<- ffs(input.df.na[,pred.names],input.df.na[,i],metric="Rsquared",
  #  method="svmRadial",tuneLength=1, verbose=FALSE,
  # trControl=fit.control)
  #m.bigRF
  ###although it might be good to run ffs can t do for stack... so pre filtering vars gives a somewhat sim result
  #def more standard 
  # How did we do?
  
  # m.bigRF$results
  #Testing the stack method with the new partitions
  tuneList<-list(
    rf=caretModelSpec(method="rf",importance=T),
    nnet=caretModelSpec(method="nnet", linout = TRUE,importance=T),
    svmRadial=caretModelSpec(method="svmRadial", linout = TRUE,importance=T),
    # gam=caretModelSpec(method="gam", linout = TRUE,importance=T),
    glm=caretModelSpec(method="glm")
  )
  models <- caretList(form, data=input.df.na, trControl=fit.control, tuneList=tuneList)
  
  results <- resamples(models)
  
  
  # stack using glm
  stackControl <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions="final")
  
  
  
  
  stack.glm <- caretEnsemble(
    models, 
    metric="RMSE",
    trControl= stackControl)
  # summary(greedy_ensemble)
  
  # Calling the model object gives us the tuning results
  stack.glm
  
  #write.csv(sum_RMSE_stack,paste("~/Dropbox/**Tesis_PHD/Near_real_time_model/results/","result_stack",names(input.df.na[i]),".csv",sep=""))
  write.csv(stack.glm$error,paste("~/Andrea/Near_real_time_model/results/","result_Ensemble_",names(input.df.na[i]),".csv",sep=""))
  
  #write.csv(x,paste("~/Andrea/Near_real_time_model/results/","result_RF_evalm",names(input.df.na[i]),".csv",sep=""))
  
  
  ### We can also generate a raster prediction and compare with the actual raster.
  
  ## Stack all predictor bands:
  #pred.stack <- stack(All_train)
  pred.stack<-stack(bioclimatic)
  pred.stack<-scale(pred.stack,center=TRUE,scale=TRUE)
  pred.raster.stack <-  raster::predict(pred.stack,stack.glm, na.rm=T)
  eval.stack.stack <- stack(raster(div.files[i]),pred.raster.stack)
  names(eval.stack.stack) <- c("Observed","Predicted")
  observed.stack<-unstack(eval.stack.stack)[[1]]
  predicted.stack<-unstack(eval.stack.stack)[[2]]
  writeRaster(observed.stack,paste(" ~/Andrea/Near_real_time_model/results/","result_stack_observed",names(input.df.na[i]),".asc",sep=""),format="ascii",overwrite=T)
  writeRaster(predicted.stack,paste(" ~/Andrea/Near_real_time_model/results/","result_stack_predicted_original",names(input.df.na[i]),".asc",sep=""),format="ascii")
  jpeg(paste("~/Andrea/Near_real_time_model/results/","result_stack_original",names(input.df.na[i]),".jpg",sep="_"))
  
  
  myPlot.stack<-levelplot(eval.stack.stack)
  print(myPlot.stack)
  dev.off()
  
  ##Stack predictors for new years (Trying with one first)
  
  test.files<-ls(pattern="test")
  ##create empty container to store all yearly predictions
  all_yearly_project<-stack()
  ##the next assumes the test_2014 onwards variables (stacks) contain all predictors and not just means
  
  
  ### Make test with biovars #### To put in different script
  
  ###At this point check that numbers correspond to the right variables
  for (j in testV){
    if( length(which(grepl(paste("test",j,sep="_"),test.files)))>0){
      temp<-mget(test.files[which(grepl(paste("test",j,sep="_"),test.files))])
      
      precip<-temp[[which(names(temp)==paste("total_test",j,sep="_"))]]
      precip<-stack(precip[[c(1,5:12,2:4)]]) ##to organize by month jan-Dec
      tmax<-temp[[which(names(temp)==paste("max_test",j,sep="_"))]]
      tmax<-stack(tmax[[c(1,5:12,2:4)]])
      tmin<-temp[[which(names(temp)==paste("min_test",j,sep="_"))]]
      tmin<-stack(tmin[[c(1,5:12,2:4)]])
      bioclimatic_test<-biovars(precip,tmin,tmax)
      plot(bioclimatic_test[[1]]*0.02-273.15,main=paste("bio1",j))
      plot(bioclimatic_test[[12]],main=paste("bio12",j))
      temp<-stack(bioclimatic_test)
      temp<-subset(bioclimatic_test,as.character(viRes@results$Variables))
      pred.stack <- scale(temp, center=TRUE, scale=TRUE)
      pred.raster.stack <-  raster::predict(pred.stack,stack.glm, na.rm=T)
      names(pred.raster.stack)<-paste("test",j,sep="_")
      # Stack obs and pred
      eval.stack.stack <- stack(raster(div.files[i]),pred.raster.stack)
      #Stack the yearly projection
      all_yearly_project<-stack(all_yearly_project,pred.raster.stack)
      names(eval.stack.stack) <- c("Observed","Predicted")
      
      observed.stack<-unstack(eval.stack.stack)[[1]]
      predicted.stack<-unstack(eval.stack.stack)[[2]]
      # writeRaster(observed.stack,paste("~/Dropbox/**Tesis_PHD/Near_real_time_model/results/","result_stack_observed",names(input.df.na[i]),".asc",sep=""),format="ascii",overwrite=T)
      #  writeRaster(predicted.stack,paste("~/Dropbox/**Tesis_PHD/Near_real_time_model/results/",paste("result_stack_predicted",j,sep="_"),names(input.df.na[i]),".asc",sep=""),format="ascii")
      writeRaster(observed.stack,paste(" ~/Andrea/Near_real_time_model/results/","result_stack_observed",names(input.df.na[i]),".asc",sep=""),format="ascii",overwrite=T)
      writeRaster(predicted.stack,paste(" ~/Andrea/Near_real_time_model/results/",paste("result_stack_predicted",j,sep="_"),names(input.df.na[i]),".asc",sep=""),format="ascii")
      
      # Plot using rastervis
      # jpeg(paste("~/Dropbox/**Tesis_PHD/Near_real_time_model/results/","result_stack",names(input.df.na[i]),j,".jpg",sep="_"))
      jpeg(paste("~/Andrea/Near_real_time_model/results/","result_stack",names(input.df.na[i]),j,".jpg",sep="_"))
      
      
      myPlot.stack<-levelplot(eval.stack.stack)
      print(myPlot.stack)
      dev.off()
    }
    else 
      print(paste("No test layers for year",j,sep=" "))
  }
  
  ##create one figure with all projections
  #names(all_yearly_project)<-testV 
  eval.stack.stack <- stack(raster(div.files[i]),all_yearly_project)
  
  # jpeg(paste("~/Dropbox/**Tesis_PHD/Near_real_time_model/results/","result_stack",names(input.df.na[i]),"all",".jpg",sep="_"))
  jpeg(paste("~/Andrea/Near_real_time_model/results/","result_stack",names(input.df.na[i]),"all",".jpg",sep="_"))
  
  
  myPlot.stack<-levelplot(eval.stack.stack)
  print(myPlot.stack)
  dev.off()
  
  
  
  
  # Get variable importance for model
  stack.varimp <- varImp(stack.glm)
  # stack.varimp<-stack.varimp$importance
  
  
  # write.csv(stack.varimp,paste("~/Dropbox/**Tesis_PHD/Near_real_time_model/",names(input.df.na[i]),"stack_variable_contribution.csv",sep=""))
  
  #  save.image(paste("~/Dropbox/**Tesis_PHD/Near_real_time_model/",names(input.df.na[i]),"_workspace.RData",sep=""))
  write.csv(stack.varimp,paste("~/Andrea/Near_real_time_model/results/",names(input.df.na[i]),"stack_variable_contribution",".csv",sep=""))
  
  save.image(paste("~/Andrea/Near_real_time_model/results/",names(input.df.na[i]),"_workspace.RData",sep=""))
  
}


