# # Import guildes
options("scipen"=100, "digits"=4)
require(sf)
require(data.table)
require(tidyverse)
require(ggplot2)
require(reshape2)
require(qgraph)
require(vegan)
require(raster)

# # # # Import obs
source("C://Dossier_Ervan/R/functions_ER.R")
system.time(obs <- FILTER.OBS(GUILD=2))
st_geometry(obs) <- "geometry"
st_crs(obs)  # should be EPSG: 21781 

# # Subset obs in river buffer
buf <- LOAD("data/buffer_river_200m_filtered_BV_polygons_sf.Rdata")

obs <- st_transform(obs,st_crs(buf))
system.time(test <- st_join(obs, buf[,c("in.riv")], join = st_intersects)) # 22 sec
table(is.na(test$in.riv)) # 111663 obs falling outside 
table(is.na(test$grid.id),test$group)

# Subset obs in buffer & keep only 
obs_sf <- subset(test,!is.na(in.riv))  # 212'424 obs

save(obs_sf,file="data/guild2/obs_guild2_sf_05.02.20.Rdata")
rm(list=c("obs2","obs","test"))

### Compute quality in river buffer
load("data/guild2/obs_guild2_sf_05.02.20.Rdata")

div <- data.table::setDT(obs_sf)[,.(N=length(unique(taxonid)),NG=length(unique(group)),Qobs=ifelse(sum(w[!duplicated(taxonid)])>=1,1,0)),by=.(grid.id)]

## Reassign spatial information
load("C:/Dossier_Ervan/R/Grid100/grid100_sf.Rdata")
div <- dplyr::right_join(grid_sf[,c("grid.id","BV04","canton","lac","geometry")],div,by="grid.id")
sf::st_geometry(div) <- "geometry"
grid <- raster::raster("C:/Dossier_Ervan/R/Grid100/grid100.tif")
raster::crs(grid) <- sp::CRS('+init=EPSG:21781')
Q_ra <- fasterize::fasterize(div[div$Qobs==1,],grid,field="NG") # quality
mapview::mapview(Q_ra,legend = T,maxpixels = 8384956,method="ngb")

raster::writeRaster(Q_ra,"data/G2_IST_04.02.20.tif",format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE) 

#########################
## Modélisation du SOLL   ds zones ouvertes
#########################
## Assign environmental variables & PPS
load("data/guild2/obs_guild2_sf_05.02.20.Rdata")
div <- data.table::setDT(obs_sf)[,.(N=length(unique(taxonid)),NG=length(unique(group)),Qobs=ifelse(sum(w[!duplicated(taxonid)])>=1,1,0)),by=.(grid.id)]
table(div$Qobs)

## Create a set of observation and pseudo-absence
VAR <- c("bio12_p", "bio6_tminc","bio4_ts", "gdd3Y_8110_ngb5_mwconic_", "topos","ndviMEAN","alt","sinuosity", "slope", "debit")

## Reassign spatial information
library(data.table)
load("C:/Dossier_Ervan/R/Grid100/grid100_sf_with_enviro.Rdata")
sf::st_geometry(grid_sf) <- "geometry"


# Merge envir variable to 'buffer' zone
buf2 <- merge(setDT(buf)[,c("grid.id", "alt", "riv.id","sinu","slope","debit","in.riv")],setDT(grid_sf)[,c("grid.id","BV04","canton","geometry","centro",..VAR[1:6])],by="grid.id",all.x=T)

summary(buf2)

## Les grand cours d'eaux ont un débit > 50 m3/sec (classe 4)
buf2 <- within(buf2,debit[is.na(debit)] <- 4)
summary(buf2[buf2$in.riv==1,])
sub <- setDT(buf2)[is.na(canton)]
st_geometry(sub) <- "geometry"
mapview::mapview(sub)

## Pour variables enviro manquantes = réassign valeurs par plus proche voisin
# st_geometry(buf2) <- "centro"
# repNA <- function(X) {
#   ID <- is.na(buf2[,X])
#   Var <- names(buf2)[X]
#   system.time(temp <- st_join(buf2[ID,"grid.id"],buf2[,X], join = nngeo::st_nn,k = 8, maxdist = 150))  # 1 sec
#   names(temp) <- c("grid.id","V1","centro")
#   TT <- setDT(temp)[,mean(V1,na.rm=T),by=grid.id]
#   buf2 <-  within(buf2,..Var[ID] <- TT$V1[match(buf2[ID,]$grid.id,TT$grid.id)])
# }


div2 <- dplyr::left_join(setDT(grid_sf)[grid.id%in%buf$grid.id,c("grid.id","BV04","canton","geometry",..VAR)],div,by="grid.id")
sf::st_geometry(div2) <- "geometry"

# Data for calibration
Q1 <- na.omit(setDT(div2)[Qobs==1,c("grid.id","Qobs","BV04",..VAR)])

# Pseudo-absence
Q0 <- dplyr::sample_n(na.omit(data.table::setDT(div2)[Qobs==0,c("grid.id","Qobs","BV04",..VAR)]),nrow(Q1),replace=F)
OBS <- rbind(Q1,Q0)

## Ici, le buffer represent l'espace guilde
env.pot <- na.omit(setDT(div2)[,c("grid.id","BV04",..VAR)]) ## 1.1 mio ha
env.pot <- merge(env.pot,div2[,c("grid.id","Qobs")],by="grid.id",all.x=T)  # reassing Qobs 

rm(list=c("obs","div","div2","Q0","Q1","PPS","grid_sf","buf","grid","obs_sf","Q_ra"))
.rs.restartR() 
# load("~/R/G14_PPS/G14_env_SOLL.RData")  # if crash

## Define receiving df and formula
fmla <- as.formula(paste("Qobs ~ ", paste(VAR,collapse = "+ ")))
PRED <- env.pot[,c("grid.id","BV04","Qobs")]

### Random forest
system.time(rf1 <- ranger::ranger(fmla, data=OBS,num.trees=1000)) # 1.5 min
system.time(PRED$rf <- predict(rf1,env.pot,type="response")$predictions) # 1,5 min
# rm(list=c("rf1"))

## Boosted regression trees
require(xgboost)
require(caret)
require(data.table)
set.seed(666)  # Pour la 'reproductibilité'
# Split data into training and test
inTrain <- caret::createDataPartition(y = OBS$grid.id, p = 0.85, list = FALSE)  # 85% des données 
training <- OBS[inTrain,]
testing <- OBS[-inTrain,]

# # ## Select best parameters
# xgb_trcontrol = trainControl(method = "cv", number = 5, allowParallel = TRUE, verboseIter = FALSE, returnData = FALSE)
# 
# xgbGrid <- expand.grid(nrounds = c(100,200),
#                        max_depth = c(5, 10),
#                        colsample_bytree = seq(0.5, 0.9, length.out = 3),
#                        ## valeurs par défaut :
#                        eta = c(0.001, 0.3, 0.5),
#                        gamma=c(0, 0.2, 0.5,2),
#                        min_child_weight = 1,
#                        subsample = 1
# )
# 
# system.time(xgb_model <- train(training[,3:13], as.factor(training$Qobs), trControl = xgb_trcontrol, tuneGrid = xgbGrid,method = "xgbTree"))  # 40 min
# B <- xgb_model$bestTune

## Parametres pour 02 (06.02.20)
B <- data.frame(nrounds = 200, max_depth = 10, eta = 0.3, gamma = 0.2,colsample_bytree = 0.7, min_child_weight = 1, subsample = 1)
params <- list(booster = "gbtree", objective = "binary:logistic", eta=B$eta, gamma=B$gamma, max_depth=B$max_depth, min_child_weight=B$min_child_weight, subsample=B$subsample, colsample_bytree=B$colsample_bytree)

# Using the inbuilt xgb.cv function, let's calculate the best nround for this model. In addition, this function also returns CV error, which is an estimate of test error.
dtrain <- xgb.DMatrix(as.matrix(training[,..VAR]), label=training$Qobs)
dtest <- xgb.DMatrix(as.matrix(testing[,..VAR]), label=testing$Qobs)

xgbcv <- xgb.cv( params = params, data = dtrain, nrounds = 500, nfold = 5, showsd = T, stratified = T, print_every_n = 10, early_stopping_rounds = 20, maximize = F)
MIN <- which.min(xgbcv$evaluation_log$test_error_mean)

#first default - model training
xgb1 <- xgb.train (params = params, data = dtrain, nrounds = MIN, watchlist = list(val=dtest,train=dtrain),  print_every_n= 10, early_stopping_rounds = 20, maximize = F , eval_metric = "error")

TEST <-  xgb.DMatrix(as.matrix(env.pot[,..VAR]),label=env.pot$Qobs)
PRED$brt <- predict (xgb1,TEST)

# GAM model
fmlaG <- as.formula(paste("Qobs ~ ", paste(VAR,collapse = "+ ")))
BAM <- mgcv::bam(fmlaG,data=OBS,family=binomial)
PRED$gam <- as.numeric(predict(BAM,env.pot,type="response"))
# rm(list=c("BAM"))

### Average model
PRED$mean <- as.numeric(rowMeans(PRED[,c("rf","brt","gam")],na.rm=T))

## Find best threshold for predictions
library(PresenceAbsence)
DD <- PRED[PRED$grid.id%in%OBS$grid.id,c("grid.id","Qobs","mean")]
OPTI <- PresenceAbsence::optimal.thresholds(DD[!is.na(DD$Qobs),])
OPTI
M <- mean(OPTI$mean)

## Set predictions
PRED$Qpred <- 0
PRED <- within(PRED,Qpred[mean>=M]<-1)

confusionMatrix (as.factor(PRED[!is.na(PRED$Qobs),]$Qpred),as.factor(PRED[!is.na(PRED$Qobs),]$Qobs))  ## Accuracy 0.87 

## Importance variable
IMPVAR <- function(x) {
  M1 <- predict(rf1,x,type="response")$predictions
  M2 <- predict(xgb1,as.matrix(setDT(x)[,..VAR]))
  M3 <- as.numeric(predict(BAM,x,type="response"))
  COR <- matrix(NA,length(VAR),3)
  for (a in 1:length(VAR)){
    temp <- setDT(x)
    IDX <- which(names(temp)==VAR[a])
    temp[,IDX] <- sample(setDT(x)[,get(VAR[a])])
    COR[a,1] <- cor.test(M1, predict(rf1,temp,type="response")$predictions)[[4]]
    COR[a,2] <- cor.test(M2,predict(xgb1,as.matrix(setDT(temp)[,..VAR])))[[4]]
    COR[a,3] <- cor.test(M3,as.numeric(predict(BAM,temp,type="response")))[[4]]
  }
  COR}

library(foreach)
library(doParallel)
registerDoParallel(makeCluster(3))
system.time(test <- foreach(i=1:10, .packages=c('dplyr','mgcv','dismo','ranger','xgboost','data.table')) %dopar% IMPVAR(dplyr::sample_n(env.pot,5000))) # 4.5 min
save(test,file="data/G2-importance_variable.Rdata")

load("data/G14-importance_variable.Rdata")
bb <- apply(matrix(unlist(lapply(test,function(x){
  b <- x[,1]
  b
})),ncol=10,byrow=F),1,mean)



LAB.VAR <- LABEL.VAR(VAR)

BB <- data.frame(var=LAB.VAR,cor=1-bb)
BB$var <- factor(BB$var,levels=LAB.VAR[order(bb)])
ggplot(BB,aes(x=var,y=cor)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=65,size=10,vjust=0.9,hjust=1)) + labs(x="",y="Importance")


# ## Performance du mod?le moyen
# load("data/G24-25_quality_predicted_20.01.20.Rdata")
# SUB <- dplyr::sample_n(PRED,5000)
# (SUB[SUB$obs==1,c("qual","mean")])
# resi = PRED[PRED$obs==1,]$mean - PRED[PRED$obs==1,]$qual
# RMSE = sqrt(mean((PRED[PRED$obs==1,]$mean - PRED[PRED$obs==1,]$qual)^2))
# 1 - var(resi) / var(PRED[PRED$obs==1,]$qual)  # R2


# ### Add spatial information and define potential habitat
load("C:/Dossier_Ervan/R/Grid100/grid100_sf.Rdata")
PRED <- merge(PRED,grid_sf[,c("grid.id","canton","geostat","lac","geometry")],by="grid.id",all.x=T)
sf::st_geometry(PRED) <- "geometry"
rm(grid_sf)

# Force IST & PPS to be in SOLL
PRED <- within(PRED,Qpred[Qobs==1]<-1)

## Ajoute les zones alluviales référencées
za <- st_read("D:/SIG/2-Eaux.dynamiques/data/Auen/Auen_fusionne_sans_rives_lacustres.shp")
za$id <- seq(1,nrow(za),1)
sf::st_geometry(PRED) <- "geometry"
za <- sf::st_transform(za,sf::st_crs(PRED))
system.time(PRED <- sf::st_join(PRED,za[,"id"],join=sf::st_intersects,left=T)) # 36 sec
dim(setDT(PRED)[duplicated(grid.id)])
PRED <- setDT(PRED)[!duplicated(grid.id)]
PRED <- within(PRED,Qpred[!is.na(id)]<-1)

save(PRED,file=paste0(getwd(),"/data/G2_quality_predicted_06.02.20.Rdata"))

## Rasterization and leaflet
load("data/G2_quality_predicted_06.02.20.Rdata")
grid <- raster::raster("D://SIG/data/grid100/grid100.tif") 
raster::crs(grid) <-  sp::CRS('+init=EPSG:21781')
sf::st_geometry(PRED) <- "geometry"
identical(sf::st_crs(grid),sf::st_crs(PRED))
SOLL <- fasterize::fasterize(PRED,grid,field="Qpred") # predicted quality
raster::writeRaster(SOLL,paste0(getwd(),"/data/G2_SOLL_06.02.20.tif"),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)


################################################
# Comparaison par cluster de bassins versants
################################################
source("C://Dossier_Ervan/R/functions_ER.R",echo=F)
library(sf)
library("tidyverse")
library(data.table)

load("D://SIG/data/Bassins_versants/BV_unifies_01.20.Rdata") # data = TT
# load("C:/Dossier_Ervan/R/Grid100/grid100_sf.Rdata")
# AREA <- data.table::setDT(grid_sf)[lac==0,.N,by=BV04]
# TT <- merge(TT,AREA,by="BV04")
BV <- LOAD("data/G2_clusterBV_06.02.20.Rdata")

BV <- dplyr::left_join(BV,TT[,c("BV04","area")],by="BV04",all.x=T)  ## area of BV without lacs
st_geometry(BV) <- "geometry"
BV$CLUST <- as.factor(paste(BV$bioreg,BV$clust,sep="_"))

load(paste0(getwd(),"/data/G2_quality_predicted_06.02.20.Rdata"))
PRED <- merge(PRED,setDT(BV)[,c("BV04","CLUST","area")],by="BV04",all.x=T)
head(data.table::setDT(PRED)[,.N,by=CLUST][order(N)])
table(is.na(PRED$CLUST))
PRED <- data.table::setDT(PRED)[!is.na(CLUST)]

## Classement par comparaison du IST
res <- data.frame("BV04"=NA,"Qobs"=NA,"EG"=NA,"BV_area"=NA,"Qprop"=NA,"bench"=NA,"deficit"=NA,"ha_to_add"=NA,"ha_available"=NA,"ha_to_add_cor"=NA,"Qpred"=NA,"reali"=NA,"rank_clust"=NA,"clust"=NA)
for (cl in levels(factor(PRED$CLUST))){
  SUB <- data.table::setDT(PRED)[CLUST==cl,]
  TAB <- data.table::setDT(SUB)[,.("Qobs"=length(Qobs[Qobs==1 &!is.na(Qobs)]),"EG"=length(Qobs),"BV_area"=round(unique(area)),"Qpred"=length(Qpred[Qpred==1])),by="BV04"]
  TAB$ha_available <- TAB$Qpred-TAB$Qobs
  TAB$Qprop = round(TAB$Qobs*100/TAB$EG,1)
  TAB <- TAB[order(TAB$Qprop,decreasing=T),]
  TAB$bench = round(quantile(TAB$Qprop,.90),1)
  deficit <- round(TAB$bench - TAB$Qprop,2)
  TAB$deficit <- round((TAB$bench - TAB$Qprop)*100/TAB$bench,0)
  TAB$ha_to_add <- ifelse(deficit/100*TAB$EG<0,0,round(deficit/100*TAB$EG))
  TAB$ha_to_add_cor <- ifelse(TAB$ha_to_add>TAB$ha_available,TAB$ha_available,TAB$ha_to_add)
  TAB$reali <- round(TAB$Qobs*100/TAB$Qpred,1)
  TAB$rank_clust<- 1:nrow(TAB)
  TAB$clust=cl
  TAB <- TAB[,c("BV04","Qobs","EG","BV_area","Qprop","bench","deficit","ha_to_add","ha_available","ha_to_add_cor","Qpred","reali","rank_clust","clust")]
  res <- rbind(res,TAB)
}
res <- res[-1,]

## Proportion IST
res$propIST <- round(res$Qobs*100/res$BV_area,1)
res <- within(res,propIST[Qobs==0] <- NA)

## Proportion SOLL
res$propSOLL <- round(res$Qpred*100/res$BV_area,1)
res <- within(res,propSOLL[Qpred==0] <- NA)
write.csv2(res,file="data/G2_tab_BV.csv",row.names = F)

# Add geometry
res <- merge(res,BV[,c("BV04","geometry")],by="BV04",al.x=T)
sf::st_geometry(res) <- "geometry"

# ## Map of clusters
# library(mapview)
# res2 <- aggregate(res,by=list(res$clust),FUN=unique,do_union=F)
# sf::st_geometry(res2) <- "geometry"
# Ncol <- as.numeric(tapply(BV$clust,factor(BV$bioreg),function(x) length(unique(x))))
# library(RColorBrewer)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# COL1=sample(col_vector, sum(Ncol))
# CLUSTERS <- res2[,1]
# names(CLUSTERS) <- c("Cluster ID","geometry")
# CLUST <- mapview::mapview(CLUSTERS,zcol="Cluster ID",col.region=COL1,legend=F,alpha.region=0.2,map.types = "Esri.WorldImagery",popup=leafpop::popupTable(res,feature.id=F),homebutton = FALSE)
# 


## IST
IST <- raster::raster(paste0(getwd(),"/data/G2_IST_04.02.20.tif"))  # done on line 56

# Couche HA to add
COL <- colorRampPalette(c("darkgreen","yellow","red"))
Ncl=5
TT <- c(round(classInt::classIntervals(res$ha_to_add,Ncl)$brks,0))
HA <- mapview::mapview(res,zcol="ha_to_add",at=TT,label="BV",col.region=COL(Ncl*2)[seq(1,Ncl*2,2)],layer.name ="ha to add",alpha.region=0.6,map.types = "Esri.WorldImagery",popup=leafpop::popupTable(res,feature.id=F),homebutton = FALSE)

# Couche HA to add_cor
# TT2 <- round(classInt::classIntervals(res$ha_to_add_cor,Ncl)$brks,0)
HA2 <- mapview::mapview(res,zcol="ha_to_add_cor",at=TT,label="BV",col.region=COL(Ncl*2)[seq(1,Ncl*2,2)],layer.name ="ha to add_cor",alpha.region=0.6,map.types = "Esri.WorldImagery",popup=leafpop::popupTable(res,feature.id=F),homebutton = FALSE)

# Couche realisation
REA <- mapview::mapview(res,zcol="reali",at=c(0,1,25,50,75,100),label="BV",col.region=rev(COL(Ncl*2)[seq(1,Ncl*2,2)]),layer.name ="Taux realisation(%)",alpha.region=0.6,map.types = "Esri.WorldImagery",popup=leafpop::popupTable(res,feature.id=F),homebutton = FALSE)

# Couche proportion PPS ds BV
PROP1 <- mapview::mapview(res,zcol="propIST",at=c(0,1,2,5,10,25),label="BV",col.region=rev(COL(Ncl*2)[seq(1,Ncl*2,2)]),layer.name ="Prop. IST (%)",alpha.region=0.6,na.color="grey40",map.types = "Esri.WorldImagery",popup=leafpop::popupTable(res,feature.id=F),homebutton = FALSE)

# Couche proportion GUILD potentiels ds BV
Ncl=6
MAX <- round(max(res$GUILDprop,na.rm=T))
PROP2 <- mapview::mapview(res,zcol="propSOLL",at=c(0,1,2,5,10,25,MAX),label="BV",col.region=rev(COL(Ncl*2)[seq(1,Ncl*2,2)]),layer.name ="Prop. SOLL (%)",alpha.region=0.6,na.color="grey40",map.types = "Esri.WorldImagery",popup=leafpop::popupTable(res,feature.id=F),homebutton = FALSE)

## SOLL
st_geometry(PRED) <- "geometry"
grid <- raster::raster("D://SIG/data/grid100/grid100.tif")
SOLL <- fasterize::fasterize(PRED,grid,field="Qpred")

# # ## Add geostat
# load("C://Dossier_Ervan/R/Grid100/grid100_sf.Rdata")
# grid <- raster::raster("D://SIG/data/grid100/grid100.tif")
# st_geometry(grid_sf) <- "geometry"
# geo_ra <- fasterize::fasterize(grid_sf,grid,field="geostat")
# mapview::mapviewOptions(raster.size=5522502)
# GEO <- mapview::mapview(geo_ra,legend = F,maxpixels = 7694880,na.color="transparent",layer.name="Geostat",col.regions=grey.colors(72),method="ngb")
# 
# geo <- raster::raster("D://SIG/data/geostat_04_09/Geo_18-19_72.tif")
# GEO <- mapview::mapview(geo,legend = F,maxpixels = 8384956,layer.name ="Geostat",method="ngb",na.color="transparent")

# m <-  mapview::mapview(grid,legend = F,maxpixels = 8384956,na.color="transparent",layer.name="GRID ID",method = "ngb",alpha.region=0) +   mapview::mapview(raster::ratify(SOLL),legend = TRUE,col.regions=c("grey40","springgreen2"),maxpixels =  8384956,na.color="transparent",layer.name="SOLL",method="ngb") + mapview::mapview(IST,legend = T,maxpixels = 8384956 ,map.types = "Esri.WorldImagery",na.color="transparent",layer.name="IST",method="ngb") + HA + HA2 + REA + PPS + PPS2 + CLUST
# m
# mapview::mapshot(m, url ="D://SIG/14-PSS/G14_IST_SOLL_04.02.20_with_clusters.html")

# Without Clusters
m2 <-  mapview::mapview(grid,legend = F,maxpixels = 8384956,na.color="transparent",layer.name="GRID ID",method = "ngb",alpha.region=0,homebutton = FALSE) +   mapview::mapview(raster::ratify(SOLL),legend = TRUE,col.regions=c("grey40","springgreen2"),maxpixels =  8384956,na.color="transparent",layer.name="SOLL",method="ngb",homebutton = FALSE) + mapview::mapview(IST,legend = T,maxpixels = 8384956 ,map.types = "Esri.WorldImagery",na.color="transparent",layer.name="IST",method="ngb",homebutton = FALSE) + HA + HA2 + REA + PROP1 + PROP2
m2
mapview::mapshot(m2, url ="D://SIG/2-Eaux.dynamiques/G2_IST_SOLL_06.02.20.html")


#### Présentation en metres linéaires to add   #####
#### Discusssion avec Stefan 07.02.20 : lisières, rivières, haies sont à exprimer en metres lineaires
source("C://Dossier_Ervan/R/functions_ER.R",echo=F)
library(sf)
library("tidyverse")
library(data.table)

load("D://SIG/data/Bassins_versants/BV_unifies_01.20.Rdata") # data = TT
BV <- LOAD("data/G2_clusterBV_06.02.20.Rdata")
BV <- dplyr::left_join(BV,TT[,c("BV04","area")],by="BV04",all.x=T)  ## area of BV without lacs
st_geometry(BV) <- "geometry"
BV$CLUST <- as.factor(paste(BV$bioreg,BV$clust,sep="_"))

load(paste0(getwd(),"/data/G2_quality_predicted_06.02.20.Rdata"))
PRED <- merge(PRED,setDT(BV)[,c("BV04","CLUST","area")],by="BV04",all.x=T)
head(data.table::setDT(PRED)[,.N,by=CLUST][order(N)])



## Classement par comparaison du IST
res <- data.frame("BV04"=NA,"Qobs"=NA,"EG"=NA,"BV_area"=NA,"Qprop"=NA,"bench"=NA,"deficit"=NA,"ha_to_add"=NA,"ha_available"=NA,"ha_to_add_cor"=NA,"Qpred"=NA,"reali"=NA,"rank_clust"=NA,"clust"=NA)
for (cl in levels(factor(PRED$CLUST))){
  SUB <- data.table::setDT(PRED)[CLUST==cl,]
  TAB <- data.table::setDT(SUB)[,.("Qobs"=length(Qobs[Qobs==1 &!is.na(Qobs)]),"EG"=length(Qobs),"BV_area"=round(unique(area)),"Qpred"=length(Qpred[Qpred==1])),by="BV04"]
  TAB$ha_available <- TAB$Qpred-TAB$Qobs
  TAB$Qprop = round(TAB$Qobs*100/TAB$EG,1)
  TAB <- TAB[order(TAB$Qprop,decreasing=T),]
  TAB$bench = round(quantile(TAB$Qprop,.90),1)
  deficit <- round(TAB$bench - TAB$Qprop,2)
  TAB$deficit <- round((TAB$bench - TAB$Qprop)*100/TAB$bench,0)
  TAB$ha_to_add <- ifelse(deficit/100*TAB$EG<0,0,round(deficit/100*TAB$EG))
  TAB$ha_to_add_cor <- ifelse(TAB$ha_to_add>TAB$ha_available,TAB$ha_available,TAB$ha_to_add)
  TAB$reali <- round(TAB$Qobs*100/TAB$Qpred,1)
  TAB$rank_clust<- 1:nrow(TAB)
  TAB$clust=cl
  TAB <- TAB[,c("BV04","Qobs","EG","BV_area","Qprop","bench","deficit","ha_to_add","ha_available","ha_to_add_cor","Qpred","reali","rank_clust","clust")]
  res <- rbind(res,TAB)
}
res <- res[-1,]

## Proportion IST
res$propIST <- round(res$Qobs*100/res$BV_area,1)
res <- within(res,propIST[Qobs==0] <- NA)

## Proportion SOLL
res$propSOLL <- round(res$Qpred*100/res$BV_area,1)
res <- within(res,propSOLL[Qpred==0] <- NA)
write.csv2(res,file="data/G2_tab_BV.csv",row.names = F)

# Add geometry
res <- merge(res,BV[,c("BV04","geometry")],by="BV04",al.x=T)
sf::st_geometry(res) <- "geometry"

# ## Map of clusters
# library(mapview)
# res2 <- aggregate(res,by=list(res$clust),FUN=unique,do_union=F)
# sf::st_geometry(res2) <- "geometry"
# Ncol <- as.numeric(tapply(BV$clust,factor(BV$bioreg),function(x) length(unique(x))))
# library(RColorBrewer)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# COL1=sample(col_vector, sum(Ncol))
# CLUSTERS <- res2[,1]
# names(CLUSTERS) <- c("Cluster ID","geometry")
# CLUST <- mapview::mapview(CLUSTERS,zcol="Cluster ID",col.region=COL1,legend=F,alpha.region=0.2,map.types = "Esri.WorldImagery",popup=leafpop::popupTable(res,feature.id=F),homebutton = FALSE)
# 


## IST
IST <- raster::raster(paste0(getwd(),"/data/G2_IST_04.02.20.tif"))  # done on line 56

# Couche HA to add
COL <- colorRampPalette(c("darkgreen","yellow","red"))
Ncl=5
TT <- c(round(classInt::classIntervals(res$ha_to_add,Ncl)$brks,0))
HA <- mapview::mapview(res,zcol="ha_to_add",at=TT,label="BV",col.region=COL(Ncl*2)[seq(1,Ncl*2,2)],layer.name ="ha to add",alpha.region=0.6,map.types = "Esri.WorldImagery",popup=leafpop::popupTable(res,feature.id=F),homebutton = FALSE)

# Couche HA to add_cor
# TT2 <- round(classInt::classIntervals(res$ha_to_add_cor,Ncl)$brks,0)
HA2 <- mapview::mapview(res,zcol="ha_to_add_cor",at=TT,label="BV",col.region=COL(Ncl*2)[seq(1,Ncl*2,2)],layer.name ="ha to add_cor",alpha.region=0.6,map.types = "Esri.WorldImagery",popup=leafpop::popupTable(res,feature.id=F),homebutton = FALSE)

# Couche realisation
REA <- mapview::mapview(res,zcol="reali",at=c(0,1,25,50,75,100),label="BV",col.region=rev(COL(Ncl*2)[seq(1,Ncl*2,2)]),layer.name ="Taux realisation(%)",alpha.region=0.6,map.types = "Esri.WorldImagery",popup=leafpop::popupTable(res,feature.id=F),homebutton = FALSE)

# Couche proportion PPS ds BV
PROP1 <- mapview::mapview(res,zcol="propIST",at=c(0,1,2,5,10,25),label="BV",col.region=rev(COL(Ncl*2)[seq(1,Ncl*2,2)]),layer.name ="Prop. IST (%)",alpha.region=0.6,na.color="grey40",map.types = "Esri.WorldImagery",popup=leafpop::popupTable(res,feature.id=F),homebutton = FALSE)

# Couche proportion GUILD potentiels ds BV
Ncl=6
MAX <- round(max(res$GUILDprop,na.rm=T))
PROP2 <- mapview::mapview(res,zcol="propSOLL",at=c(0,1,2,5,10,25,MAX),label="BV",col.region=rev(COL(Ncl*2)[seq(1,Ncl*2,2)]),layer.name ="Prop. SOLL (%)",alpha.region=0.6,na.color="grey40",map.types = "Esri.WorldImagery",popup=leafpop::popupTable(res,feature.id=F),homebutton = FALSE)

## SOLL
st_geometry(PRED) <- "geometry"
grid <- raster::raster("D://SIG/data/grid100/grid100.tif")
SOLL <- fasterize::fasterize(PRED,grid,field="Qpred")

# # ## Add geostat
# load("C://Dossier_Ervan/R/Grid100/grid100_sf.Rdata")
# grid <- raster::raster("D://SIG/data/grid100/grid100.tif")
# st_geometry(grid_sf) <- "geometry"
# geo_ra <- fasterize::fasterize(grid_sf,grid,field="geostat")
# mapview::mapviewOptions(raster.size=5522502)
# GEO <- mapview::mapview(geo_ra,legend = F,maxpixels = 7694880,na.color="transparent",layer.name="Geostat",col.regions=grey.colors(72),method="ngb")
# 
# geo <- raster::raster("D://SIG/data/geostat_04_09/Geo_18-19_72.tif")
# GEO <- mapview::mapview(geo,legend = F,maxpixels = 8384956,layer.name ="Geostat",method="ngb",na.color="transparent")

# m <-  mapview::mapview(grid,legend = F,maxpixels = 8384956,na.color="transparent",layer.name="GRID ID",method = "ngb",alpha.region=0) +   mapview::mapview(raster::ratify(SOLL),legend = TRUE,col.regions=c("grey40","springgreen2"),maxpixels =  8384956,na.color="transparent",layer.name="SOLL",method="ngb") + mapview::mapview(IST,legend = T,maxpixels = 8384956 ,map.types = "Esri.WorldImagery",na.color="transparent",layer.name="IST",method="ngb") + HA + HA2 + REA + PPS + PPS2 + CLUST
# m
# mapview::mapshot(m, url ="D://SIG/14-PSS/G14_IST_SOLL_04.02.20_with_clusters.html")

# Without Clusters
m2 <-  mapview::mapview(grid,legend = F,maxpixels = 8384956,na.color="transparent",layer.name="GRID ID",method = "ngb",alpha.region=0,homebutton = FALSE) +   mapview::mapview(raster::ratify(SOLL),legend = TRUE,col.regions=c("grey40","springgreen2"),maxpixels =  8384956,na.color="transparent",layer.name="SOLL",method="ngb",homebutton = FALSE) + mapview::mapview(IST,legend = T,maxpixels = 8384956 ,map.types = "Esri.WorldImagery",na.color="transparent",layer.name="IST",method="ngb",homebutton = FALSE) + HA + HA2 + REA + PROP1 + PROP2
m2
mapview::mapshot(m2, url ="D://SIG/2-Eaux.dynamiques/G2_IST_SOLL_06.02.20.html")
