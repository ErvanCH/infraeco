library(rgdal)
library(raster)
library(ade4)
library(adehabitatHR)
library(adehabitatMA)
library(ecospat)
library(sf)
library(data.table)
library(tidyverse)
source("C://Dossier_Ervan/R/functions_ER.R",echo=F)

## load
env.pot <- LOAD("envpot","data")
env.pot$bioregion <- env.pot$subreg
env.pot$Qobs <- NULL
BIOR <- split(env.pot,env.pot$bioregion)

AXES <-data.frame(a=c(1,1,2),b=c(2,3,3))
R = 100
mask <- ascgen(SpatialPoints(cbind((0:(R))/R, (0:(R)/R))),
               nrcol = R - 2, count = FALSE)

# source('~/R/functions_ER.R')
# system.time(L <- lapply(BIOR, function(x) pca.rasterize(x,axes=AXES,mask=MASK,R=100,z.th=0.01,save_raster=T))) # 8 min

dest.folder=paste0(getwd(),"/rasterBV/")

library(parallel)
ncores<-detectCores()-1
cl <- makeCluster(ncores,outfiles = 'dest.folder' )

clusterExport(cl,varlist=c('pca.rasterize',
                           'rasterizeBV',
                           'env.overlap',
                           'AXES',
                           'mask',
                           'R'),
              envir=.GlobalEnv)

clusterEvalQ(cl,{
  library(rgdal)
  library(raster)
  library(ade4)
  library(adehabitatHR)
  library(adehabitatMA)
  library(ecospat)
  library(sf)
  library(data.table)
  library(tidyverse)
  library(CircStats)
})

system.time(L <- parLapply(cl,BIOR, pca.rasterize,VAR=VAR,axes=AXES,mask=mask,R=R,z.th=0.01,save_raster=T)) #4 min
save(L,file=paste0("rasterBV/",GUILD,"_comp_bv.RData"))

## Clustering par bioregion
# Clustering
L <- LOAD("comp_bv","rasterBV")
CENTERS <- min(sapply(L,function(x) length(unique(x$BV1))))
test2 <- lapply(L,CLUST,CENTERS)
test <- rbindlist(test2)
test$BASIS_NR <- as.numeric(substr(test$BV04,2,8))
test$BV04 <- as.numeric(substr(test$BV04,2,8))
test$bioreg <- rep(1:12,as.numeric(sapply(test2,nrow)))

## Ensure that clusters don't have single BV (N>5)
load("D://SIG/data/Bassins_versants/BV_unifies_07.20.Rdata")  ## data = TT
BV <- dplyr::right_join(test,TT[,c("BV04","subreg")])  ## add BV without clusters
st_geometry(BV) <- "geometry"
BV$area <- as.numeric(sf::st_area(BV))/10000
BV$CLUST <- as.factor(paste(BV$subreg,BV$clust,sep="_"))

BV <- within(BV,CLUST[is.na(clust)]<-NA)
table(BV$CLUST)

## Nbr de BV par cluster
COUNT <- setDT(BV)[,.N[],by=CLUST]
summary(COUNT$V1)
table(COUNT$V1<10)  # 10 cluster avec moins de 5 BV 

BV$col <- 1
BV <- within(BV,col[CLUST%in%COUNT$CLUST[COUNT$V1<5] | is.na(clust)]<-2)
table(BV$col)
st_geometry(BV) <- "geometry"
BV$centro <- st_geometry(st_centroid(BV))
st_geometry(BV) <- "centro"
SINGLEBV <- BV[BV$col==2,]$BV04
if(length(SINGLEBV)>0) {
for (i in 1:length(SINGLEBV)) {
  SUBREG <- unique(BV[BV$BV04==SINGLEBV[i],]$subreg)
  T <- nngeo::st_nn(BV[BV$BV04==SINGLEBV[i],],BV[BV$col==1 & BV$subreg==SUBREG,],k=1,progress = FALSE)
  BV <- within(BV,CLUST[BV$BV04==SINGLEBV[i]] <- BV[BV$col==1 & BV$subreg==SUBREG,]$CLUST[as.numeric(T)])
}
table(is.na(BV$CLUST))
BV$CLUST <- factor(BV$CLUST)
}

# ## Check result
library(RColorBrewer)
st_geometry(BV) <- "geometry"
BV2 <- aggregate(BV,by=list(BV$CLUST),FUN=unique,do_union=F)
sf::st_geometry(BV2) <- "geometry"
Ncol <- as.numeric(tapply(BV$clust,factor(BV$bioreg),function(x) length(unique(x))))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
COL1=sample(col_vector, sum(Ncol))
CLUSTERS <- BV2[,1]
names(CLUSTERS) <- c("Cluster ID","geometry")
CLUST <- mapview::mapview(CLUSTERS,zcol="Cluster ID",col.region=COL1,legend=F,alpha.region=0.8,map.types = "Esri.WorldImagery",popup=leafpop::popupTable(CLUSTERS,feature.id=F),homebutton = FALSE)
CLUST +BV

## Check isolated BV (not surrounded by any other member of the same cluster)
st_geometry(BV) <- "geometry"
BV$outlier=0
    for (i in 1:1029){
      idx <- unlist(st_intersects(BV[i,], BV))
      idx <- idx[-which(idx==i)]
      test <-  match(BV$CLUST[i],BV$CLUST[idx])
     if (is.na(test)) {
     SUBREG <- unique(BV[i,]$subreg)
     B <- BV[idx,] 
     if (length(B[B$subreg==SUBREG,]$CLUST)>0) {
     BV <- within(BV,CLUST[i] <- names(which.max(table(B[B$subreg==SUBREG,]$CLUST))))
     BV$outlier[i] <- 1
       }}
 }
# table(BV$outlier)
# tapply(BV$subreg,BV$CLUST,function(x) length(unique(x)))

## Check isolated BV (based on distance to centroids of cluster)
st_geometry(BV) <- "centro"
BV$outlier=0
BV$col <- 1
for ( a in unique(BV$CLUST)) {
  centro_clust <- st_geometry(st_centroid(st_union(BV[BV$CLUST==a,])))
  dist2centro <- unlist(nngeo::st_nn(BV[BV$CLUST==a,],centro_clust,returnDist=TRUE,progress=FALSE)[[2]])
  d_mean <- mean(dist2centro)
  outliers <- BV[BV$CLUST==a,]$BV04[which(dist2centro>(mean(dist2centro)+2*sd(dist2centro)))]
  BV <- within(BV,col[BV04%in%outliers] <- 2)
  if (length(outliers)>0){
    for (j in outliers){
      SUBREG <- unique(BV[BV$BV04==j,]$subreg)
      NN <- nngeo::st_nn(BV[BV$BV04==j,],BV[BV$col==1 & BV$subreg==SUBREG,],k=1,progress=FALSE)
      BV <- within(BV,CLUST[BV$BV04==j] <- BV[BV$col==1 & BV$subreg==SUBREG,]$CLUST[as.numeric(NN)])
      BV <- within(BV,outlier[BV$BV04==j] <- 1)
    }}
}
# table(BV$outlier)
# summary(tapply(BV$subreg,BV$CLUST,function(x) length(unique(x))))


# ## Check result
# library(RColorBrewer)
# st_geometry(BV) <- "geometry"
# BV2 <- aggregate(BV,by=list(BV$CLUST),FUN=unique,do_union=F)
# sf::st_geometry(BV2) <- "geometry"
# Ncol <- as.numeric(tapply(BV$clust,factor(BV$bioreg),function(x) length(unique(x))))
# 
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# COL1=sample(col_vector, sum(Ncol))
# CLUSTERS <- BV2[,1]
# names(CLUSTERS) <- c("Cluster ID","geometry")
# CLUST <- mapview::mapview(CLUSTERS,zcol="Cluster ID",col.region=COL1,legend=F,alpha.region=0.8,map.types = "Esri.WorldImagery",popup=leafpop::popupTable(CLUSTERS,feature.id=F),homebutton = FALSE)
# CLUST + BV

save(BV,file=paste0("data/",GUILD,"_clusterBV_",DATE,".RData"))
