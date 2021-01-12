#------------------------------------------------------------------#
#-------------- Mapping floristic patterns using the --------------#
#----------------- GIS-layers as predictors -----------------------#
#-------(Amazon Landsat TM/ETM+ imangery bands, CEC, HAND) --------#
#----------------- in Jurua, Brazilian Amazonia -------------------#
#------------------------------------------------------------------#
# by Gabriela Zuquim 
#We thank Pablo Perez Chaves for support in code development

#please cite doi:

Packages <- c("raster", "caret", "mapview", "sf", "CAST", "vegan", 
              "cluster", "mvtnorm","dendextend","dplyr",
              "randomForest","rasterVis")
lapply(Packages, library, character.only = T)


#------------------------------------------------------------------#
##### Ordinations ####
# define community data
spdata <-  read.csv("Jurua_Palm_data_v2.csv", header=TRUE, stringsAsFactors=FALSE)

#------------------------------------------------------------------#

#####   PART 1B:  Ordination plots with classification results #####
# checked 2020-02-04

# run ordination and classification 
#dist matrix all species
  # run distances
  dist.pa <- vegdist(spdata[,-1],binary=T,diag=T,upper=F) # presence-absence floristic matrices
  dist.pa <- stepacross(dist.pa) # calculate extended dissimilarities
  
  # run ordinations
  mds.pa <- monoMDS(dist.pa, y = cmdscale(dist.pa, k=3), k = 3, model = "global", threshold = 0.8, maxit = 200, weakties = TRUE, stress = 1,
                    scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7, sratmax=0.99999) 
  
  #Merging the floristic ordinations with the environmental table
  #the ordination method was non-metric multidimensional scaling (NMDS) 
  
  nmds_pa<-as.data.frame(mds.pa$points)
  names(nmds_pa)<-c("nmds1pa","nmds2pa","nmds3pa")
  
  # Make the ordination plot
  plot(nmds_pa$nmds1pa, nmds_pa$nmds2pa, cex=1, pch=1,
       lwd=2, xlab="NMDS1", ylab="NMDS 2")
  
#End of ordinations#

#Get data of predictors obtained from GIS-layers#
  envi <-  read.csv("plot_predictors.csv")
  envi<-envi[,-1]

  #--------------------------#
  #Random Forest regression
  #--------------------------#
  set.seed(100)
  
  #Create a 10-fold crossvalidation set
  ctrl<-trainControl(method="cv",number=10,savePredictions=TRUE)
  
  #model each ordination axis using all variables (the 30m median Landsat, bioclim, logCEC, Hand)
  nmds1_all<-train(envi[,-c(1:3)],nmds_pa[,"nmds1pa"],method="rf",metric="Rsquared",
                   trControl=ctrl,importance=TRUE,ntree=1500)
  nmds2_all<-train(envi[,-c(1:3)],nmds_pa[,"nmds2pa"],method="rf",metric="Rsquared",
                   trControl=ctrl,importance=TRUE,ntree=1500)
  nmds3_all<-train(envi[,-c(1:3)],nmds_pa[,"nmds3pa"],method="rf",metric="Rsquared",
                   trControl=ctrl,importance=TRUE,ntree=1500)
  
  
  # Forward feature selection (ffs)
  # models all possible combination of pairs of variables until it find the pair with better performance
  # performance can be Rsquared or RSME
  # after the best pair of variables it adds a third variable until find the best
  # it repeats this until there are no more additional variables that contribute
  
  nmds1_ffs <-  ffs(envi[,-c(1:3)],nmds_pa[,"nmds1pa"],method="rf",metric="Rsquared",
                    tuneGrid=data.frame("mtry"=2),
                    trControl=ctrl,ntree=1500)
  nmds2_ffs <-  ffs(envi[,-c(1:3)],nmds_pa[,"nmds2pa"],method="rf",metric="Rsquared",
                    tuneGrid=data.frame("mtry"=2),
                    trControl=ctrl,ntree=1500)
  nmds3_ffs <-  ffs(envi[,-c(1:3)],nmds_pa[,"nmds3pa"],method="rf",metric="Rsquared",
                    tuneGrid=data.frame("mtry"=2),
                    trControl=ctrl,ntree=1500)
  
  #----------------------------------------------#
# Do spatial predictions
  predictors<-brick("jurua_predictors_selected.tif")
  names(predictors)<-c("Band1_median15", "Band2_median15", "Band3_median15", "Band4_median15", "Band5_median15", "Band7_median15",
                       "HAND_50","logCEC")

  r.nmds1.ffs <- predict(predictors,nmds1_ffs)
  r.nmds2.ffs <- predict(predictors,nmds2_ffs)
  r.nmds3.ffs <- predict(predictors,nmds3_ffs)
  
  #Visually assess how each ordination axis looks like
  library(rasterVis)
  library(RColorBrewer)
  plot(r.nmds1.ffs,col=brewer.pal(n = 10, name = "Spectral"),legend=T, cex.axis=1.5,legend.width=1.5)
  plot(r.nmds2.ffs,col=brewer.pal(n = 10, name = "Spectral"),legend=T, cex.axis=1.5,legend.width=1.5)
  plot(r.nmds3.ffs,col=brewer.pal(n = 10, name = "Spectral"),legend=T, cex.axis=1.5,legend.width=1.5)
  
  # stack floristic predictions
  r.nmds.ffs<-stack(r.nmds1.ffs,r.nmds2.ffs,r.nmds3.ffs)
  names(r.nmds.ffs)<-c(paste0("nmds1",plant[i]), paste0("nmds2",plant[i]), paste0("nmds3",plant[i]))
  
  #Export the floristic ordinations
  writeRaster(r.nmds.ffs,filename="floristic_predictions_nmds_ffs.tif",format="GTiff",overwrite=T)

  ###End of predictions using Random Forest###

#######Area of Applicability maps#########
##selected variables
  
  r.nmds1.aoa <- aoa(predictors,nmds1_ffs)
  r.nmds2.aoa <- aoa(predictors,nmds2_ffs)
  r.nmds3.aoa <- aoa(predictors,nmds3_ffs)
  
  #stack and save results
  r.nmds.aoa<-stack(r.nmds1.aoa,r.nmds2.aoa,r.nmds3.aoa)
  names(r.nmds.aoa)<-c("DI_nmds1","aoa_nmds1","DI_nmds2","aoa_nmds2","DI_nmds3","aoa_nmds3")
 #visualize DI and AOA maps
  plot(r.nmds.aoa)
  writeRaster(r.nmds.aoa,filename="floristic_predictions_nmds_aoa.tif",format="GTiff",overwrite=T)
