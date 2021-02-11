library(sf)
library(raster)
library(rgdal)
library(rgeos)
library(tidyverse)

riv <- st_read("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/Typisierung_LV03/FGT.shp")

sinuosity <- function(polyline,mesh) {
  if (unique(polyline$SHAPE_Leng)<60){
    P2 <- st_sf(data.frame(id=1,sinu=1,polyline[c(1,3,4,5,7,14,18,21,24)],geometry=st_geometry(st_centroid(polyline$geometry))))  
  } else if (unique(polyline$SHAPE_Leng)<80) {
    npoints <-  unique(polyline$SHAPE_Leng)/mesh
    POINTS <- st_sample(polyline[[25]],npoints,type = "regular")
    PP <- st_cast(POINTS,"POINT")  # split into single points
    SINU <- as.numeric(st_distance(PP[1],PP[3])/(mesh*2))
    P2 <- st_sf(data.frame(id=1,sinu=SINU,polyline[c(1,3,4,5,7,14,18,21,24)],geometry=st_geometry(st_centroid(polyline$geometry)))) 
  } else if(unique(polyline$SHAPE_Leng)>=80) {
    npoints <-  unique(polyline$SHAPE_Leng)/mesh
    POINTS <- st_sample(polyline[[25]],npoints,type = "regular")
    PP <- st_cast(POINTS,"POINT")  # split into single points
    
    SINU <- data.frame(id=1:length(PP),sinu=rep(NA,length(PP)))
    for (a in 1:length(PP)){
      if (a <= 2){
        SINU[a,2] <- as.numeric(st_distance(PP[a],PP[a+2])/(mesh*2))
      } else if (a >= length(PP)-2){
        SINU[a,2] <- as.numeric(st_distance(PP[a],PP[a-2])/(mesh*2))
      } else {
        SINU[a,2] <- as.numeric(st_distance(PP[a-2],PP[a+2])/(mesh*4))
      }}
    
    P1 <- merge(st_sf(id=1:length(PP),PP),SINU,by="id")
    P2 <- merge(P1,polyline[c(1,3,4,5,7,14,18,21,24)])
  }
  P2
}
# Debug
# polyline <- as.list(riv[riv$OBJECTID_G==14059367,])
# polyline <- as.list(riv[riv$OBJECTID_G==14078378,])
# sinuosity(polyline,20)
# # options(error = browser)
# res <- apply(riv,1,function(x) 
#   { res <- tryCatch(sinuosity(x,20), error=function(e) { 
#   cat("Failed on x = ", x["OBJECTID_G"], "\n", sep=" /") 
#   stop(e)})
# })
# # options(error=NULL)  
# 
# # Debug
# mesh=20
# for ( i in 1:nrow(riv)) {
#   polyline <- as.list(riv[i,])
#   if (unique(polyline$SHAPE_Leng)<60){
#     P2 <- st_sf(data.frame(id=1,sinu=1,polyline[c(3,4,5,7,14,18,21,24)],geometry=st_centroid(polyline$geometry[1])))  
#   } else if (unique(polyline$SHAPE_Leng)<80) {
#     
#     npoints <-  unique(polyline$SHAPE_Leng)/mesh
#     POINTS <- st_sample(polyline[[25]],npoints,type = "regular")
#     PP <- st_cast(POINTS,"POINT")  # split into single points
#     
#     SINU <- data.frame(id=1:length(PP),sinu=rep(NA,length(PP)))
#     for (a in 1:length(PP)){
#       if (a = 1){
#         SINU[a,2] <- as.numeric(st_distance(PP[a],PP[a+2])/(mesh*2))
#       } else if (a >= length(PP)-2){
#         SINU[a,2] <- as.numeric(st_distance(PP[a],PP[a-2])/(mesh*2))
#       } else {
#         SINU[a,2] <- as.numeric(st_distance(PP[a-2],PP[a+2])/(mesh*2))
#       }}
#   } else {
#   npoints <-  unique(polyline$SHAPE_Leng)/mesh
#   POINTS <- st_sample(polyline[[25]],npoints,type = "regular")
#   PP <- st_cast(POINTS,"POINT")  # split into single points
#   
#   SINU <- data.frame(id=1:length(PP),sinu=rep(NA,length(PP)))
#   for (a in 1:length(PP)){
#     if (a <= 4){
#       SINU[a,2] <- as.numeric(st_distance(PP[a],PP[a+4])/(mesh*4))
#     } else if (a >= length(PP)-4){
#       SINU[a,2] <- as.numeric(st_distance(PP[a],PP[a-4])/(mesh*4))
#     } else {
#       SINU[a,2] <- as.numeric(st_distance(PP[a-2],PP[a+2])/(mesh*4))
#     }}
#   
#   P1 <- merge(st_sf(id=1:length(PP),PP),SINU,by="id")
#   P2 <- merge(P1,polyline[c(3,4,5,7,14,18,21,24)])
#   }
#   if(i==1) {res <- P2}
#   if(i%in%seq(0,180000,100)) {cat(paste(i,"..."))}
#   res <- rbind(res,P2)
#   }

library(parallel)
ncores<-detectCores()-1
if (ncores == 0) {ncores = 1}
cl <- makeCluster(ncores)
clusterExport(cl,varlist=ls(),envir=.GlobalEnv)
clusterEvalQ(cl,{
  library(sf)
  library(raster)
  library(rgdal)
  library(rgeos)
})
system.time(res <- parApply(cl,riv,1,function(x) sinuosity(x,20))) # 11 minutes
res[[105225]]  # check the result

save(res,file=paste0(getwd(),"/data/sinuostiy.Rdata"))

     