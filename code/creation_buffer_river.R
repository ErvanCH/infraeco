### Workflow
# 1. create buffer of 500m around rivers (st_buffer)
# 2. buffer to pixel (fasterize)
# 3. centroid of each buffer cell 
# 4. promote river to polypoints  cast(riv,"POINTS")
# 5. extract altitude at river points & buffer centroides    extract(MNT,rivpts)
# 6. associate buffer cell to nearest riv alt  st_join(buffer.centro,rivpts,join=st_nearest_points)
# 7. substract riv_alt and buffer_alt




# Create buffer
# library(RQGIS3)
# set_env(root = "C:/Program Files/QGIS 3.4",new=T)
# open_app()
# find_algorithms(search_term = "buffer", name_only = TRUE)
# get_usage(alg = "native:buffer")
# buff <- run_qgis(alg = "native:buffer",
#                 INPUT = "C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/Typisierung_LV03/FGT.shp",
#                 DISTANCE=500,
#                 OUTPUT="C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/Export_shp/buffer_rivers/buffer_500m.shp")

# 2. Extract altitude for each buffer cell (fast)
require(data.table)
require(tidyverse)
require(ggplot2)
require(reshape2)
require(qgraph)
require(vegan)
require(raster)
require(sf)
library(rgdal)


# #Making SpatialPolygons from list of polygons
# poly <- SpatialPolygonsDataFrame(SpatialPolygons(buf@polygons),data.frame(id=buf$OBJECTID_G),match.ID = F)
# poly@data$layer <- as.numeric(as.character(poly@data$id ))
# 
# buf2 <- st_read("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/Export_shp/Buffer_river/buffer_riv_200m.shp")
# rast_buf <- function(x) {
#   g <- crop(grid100,extent(st_coordinates(x$geometry)[,1:2]))
#   x_sf <- st_sf(id=x$OBJECTID_G,geometry=st_geometry(x$geometry))
#   s <- fasterize::fasterize(x_sf,g,field="id")
#   s
# }
# 
# system.time(test <- apply(buf2, 1,function(x) rast_buf(x)))
# 
# library('parallel')
# cl <- makeCluster(7)
# z <- clusterEvalQ(cl, library("raster","sf"))
# clusterExport(cl, c("buf2", "grid100","rast_buf"))
# system.time(p <- parApply(cl = cl, buf2, 1, function(x) rast_buf(x)))
# stopCluster(cl)

buf <- st_read("D:/SIG/2-Eaux.dynamiques/Export_shp/Buffer_river/buffer_riv_200m_union.shp")
grid100 = raster("C://Dossier_Ervan/R/Grid100/grid100.tif")
buf_ra <- fasterize::fasterize(buf,grid100)
buf_df <- as(buf_ra, "SpatialPixelsDataFrame")
mnt <- raster("D:/SIG/2-Eaux.dynamiques/data/mnt/mnt100tot1.tif")
# buf_df@data$alt <- extract(grid100,buf_df)
buf_df@data$alt <- extract(mnt,buf_df)
## Assign official grid id
buf_df@data$grid.id <- extract(grid100,buf_df)
# Commute to sf object
temp <- data.frame(buf_df@data,x=buf_df@coords[,1],y=buf_df@coords[,2])
buf_sf <- st_as_sf(temp,coords=c("x","y"),crs=21781)
rm(list=c("buf", "buf_ra", "buf_df","grid100"))# Remove temporary objects

# # # 3. Extract altitude of rivers
# riv <- st_read("D:/SIG/2-Eaux.dynamiques/Typisierung_LV03/FGT.shp")
# system.time(riv_pts <- st_line_sample(riv,density=0.02,type="regular")) ## 300 sec # 1pt every 50m
# 
# riv_sf <- st_sf(riv$OBJECTID_G, riv$SHAPE_Leng,riv_pts)
# riv_pts_sf <- st_cast(riv_sf,"POINT")
# names(riv_pts_sf) <- c("riv.id","length","geometry")
# st_geometry(riv_pts_sf) <- "geometry"
# st_coordinates(riv_pts_sf %>% filter(riv.id==1667031))
# rm(list=c("riv_sf","riv","riv_pts"))
# 
# library(stars)
# mnt <- raster("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/data/mnt/mnt100tot1.tif")
# mnt_sf <- st_as_sf(st_as_stars(mnt), as_points = F, merge = FALSE)
# 
# setDT(riv_pts_sf)[,"pt.id":=seq_len(.N),by=riv.id]
# riv_pts_sf$riv_pt <- paste(riv_pts_sf$riv.id,riv_pts_sf$pt.id,sep="-")
# riv_pts_sf <- st_as_sf(riv_pts_sf)
# mnt_sf <- st_transform(mnt_sf,st_crs(riv_pts_sf))
# 
# system.time(riv_pts_sf <- st_join(riv_pts_sf,mnt_sf,join=st_intersects))  ## 118 sec
# rm(mnt_sf)
# 
# names(riv_pts_sf) <- c("riv.id","length","pt.id","riv_pt","alt.riv","geometry")
# riv <- st_transform(riv_pts_sf,st_crs(buf_sf))
# 
# # ## remove (erroneous) double sampling point for river < 100m
# # table(riv[riv$length<100,]$pt.id)
# # IDX <- which(riv$length<100 & riv$pt.id==2)
# # riv <- riv[-IDX,]   
# save(riv,file="data/guild2/riv_pts_50m.Rdata")


# Find nearest river point
load("C://Dossier_Ervan/R/G2-Eaux_dynamiques/data/guild2/riv_pts_50m.Rdata")
system.time(buf_sf <- st_join(buf_sf, riv, nngeo::st_nn, k = 1, maxdist = 300))  # 45 sec
buf_sf <- buf_sf[,-1]

# # Estimate distance to river for each obs within 500m 
# x <- as.list(test[1,])
# dist2river <-  function(x) {
#   GEOM <- st_geometry(x$geometry)
#   RIV <- test[test$layer%in%as.numeric(x["layer"]),] 
#   st_crs(GEOM) <- st_crs(RIV)
#   dist <- RIV$alt[which.min(st_distance(GEOM,RIV))]
#   dist
# }
# 
# library(parallel)
# ncores<-detectCores()-1
# if (ncores == 0) {ncores = 1}
# cl <- makeCluster(ncores)
# clusterExport(cl,varlist=c("buf_sf","dist2river","test"),envir=.GlobalEnv)
# clusterEvalQ(cl,{
#   library(sf)
#   library(raster)
#   library(rgdal)
#   library(rgeos)
# })
# system.time(buf_sf$alt_riv<- parApply(cl,buf_sf, 1, function(x) dist2river(x)))  # 8.4 hrs buffer 500m et 2.8 hrs buffer 200m
# stopCluster(cl)

# narrow buffer to area with delta alt < 10m 
buf_sf$dZ <- buf_sf$alt-buf_sf$alt.riv
setwd("C://Dossier_Ervan/R/G2-Eaux_dynamiques/")
save(buf_sf,file="data/guild2/buffer_river_with.alt_02.20.Rdata")

## Elimination des points seuls

# Subset elevation > 10 m
buf_sf$keep <- 1
buf_sf <- within(buf_sf,keep[dZ>=10] <- 0)
buf_sf2 <- setDT(buf_sf)[keep==1,]

(nrow(buf_sf)-nrow(buf_sf2))/nrow(buf_sf) # 36% discarded 
rm(list=c("mnt","grid_sf","riv","riv_ra","riv_ra_sf"))

# Avoid discarding pixels falling on rivers

## Extract cells crossed by rivers
## rasterize river
riv <- st_read("D:/SIG/2-Eaux.dynamiques/Typisierung_LV03/FGT.shp")
grid100 = raster("C://Dossier_Ervan/R/Grid100/grid100.tif")
values(grid100) <- NA
library(stars)
riv_ra <- st_rasterize(riv,st_as_stars(grid100,values=NA_real_))
# write_stars(riv_ra, "data/guild2/shp/riverFGT_raster.tif")

# Take cell ID falling on river pixel
riv_ra_sf <- st_as_sf(riv_ra, as_points = F, merge = FALSE)
load("C:/Dossier_Ervan/R/Grid100/grid100_sf.Rdata")
st_geometry(grid_sf) <- "centro"
system.time(riv_ra_sf <- st_join(riv_ra_sf,grid_sf[,"grid.id"],join=st_intersects) ) # 70 sec
table(is.na(riv_ra_sf$grid.id))
mapview::mapview(riv_ra_sf[is.na(riv_ra_sf$grid.id),])  # pixels hors frontières

ID <- riv_ra_sf$grid.id[is.na(match(riv_ra_sf$grid.id,buf_sf2$grid.id))] # compare grid ID that are not in buf_sf2
mapview::mapview(dplyr::sample_n(riv_ra_sf[riv_ra_sf$grid.id%in%ID,],1000))

cell2add <- buf_sf[buf_sf$grid.id%in%ID,]
cell2add$keep <- 1
# buf_sf2$alt.riv <- NULL
# buf_sf2$dZ<- NULL
# buf_sf2$keep <- NA
buf_sf2 <- rbind(buf_sf2,cell2add)
rm(list=c("BB", "cell2add", "grid100", "ID", "LOF", "riv", "riv_ra"))

# Add grid.id falling in zones alluviales fusionnées
za <- st_read("D:/SIG/2-Eaux.dynamiques/data/Auen/Auen_fusionne_sans_rives_lacustres.shp")
st_geometry(grid_sf) <- "centro"
za <- st_transform(za,st_crs(grid_sf))
test <- st_join(za,grid_sf[,c("grid.id")],join=st_intersects)
cell2add <- buf_sf[buf_sf$grid.id%in%test$grid.id,]
cell2add$keep <- 1
buf_sf2 <- rbind(buf_sf2,cell2add)

rm(list=c("buf_sf","cell2add","grid_sf","temp","test","za"))

# Remove isolated points (noize)
temp <- data.table(layer=buf_sf2$riv.id,grid.id=buf_sf2$grid.id,keep=buf_sf2$keep,st_coordinates(buf_sf2$geometry))
temp <- temp[!is.na(X),]
temp <- setorder(temp,layer,X,Y)
temp[,"lof":=dbscan::lof(.SD[,c("X","Y")],3)]
temp$out <- NA
temp <- within(temp,out[is.na(keep) & lof>1.5] <- 1)

# discard_outliers <- function(x) {
#   A <- rep(1,nrow(x))
#   if (nrow(x)>3){
#     LOF <- dbscan::lof(x[,c("X","Y")],3)
#     A[2:(nrow(x)-1)] <- ifelse(LOF[2:(nrow(x)-1)]>1.5,2,1)
#   }
#   A
# }
# temp$out <- as.numeric(NA)
# temp[,"out":=discard_outliers(.SD)]
# table(temp$out)


# riv <- st_read("D:/SIG/2-Eaux.dynamiques/Typisierung_LV03/FGT.shp")
# riv <- st_transform(riv,st_crs(buf_sf2))
# 
# ID <- unique(temp[temp$out==1,]$layer)
# BB <- st_bbox(buf_sf2 %>% filter(riv.id==ID[2])) + c(-400,-400,400,400)
# PT <- st_crop(buf_sf2,BB)
# RIV <- st_crop(riv,BB)
# OUT <- st_as_sf(temp %>% filter(layer==ID[2]),coords=c("X","Y"),crs=st_crs(buf_sf2))
# 
# ggplot() +
#   geom_sf(data=PT,col="grey50") +
#   geom_sf(data=RIV,col="cyan") +
#   geom_sf(data=OUT,aes(col=as.factor(out))) +
#   theme(legend.position = "none")

## rajoute un identifiant indiquant si rivière ou non
buf_sf2$in.riv <- 0
buf_sf2 <- within(buf_sf2,in.riv[grid.id%in%riv_ra_sf$grid.id] <- 1) 
table(buf_sf2$in.riv)
table(duplicated(buf_sf2$grid.id))
buf_sf2 <- buf_sf2[!is.na(buf_sf2$grid.id),]
buf_sf2 <- buf_sf2[!duplicated(buf_sf2$grid.id),]

buf_sf2 <- merge(buf_sf2,temp[,c("grid.id","out")],by="grid.id",all.x=T)
buf_sf2 <- buf_sf2[is.na(buf_sf2$out),]  # remove outliers


buf_sf2$pt.id <- NULL
buf_sf2$keep <- NULL
buf_sf2$out <- NULL
rm(list=c("temp","cell2add","riv","ID"))
save(buf_sf2,file="data/guild2/buffer_river_200m_filtered_02.20.Rdata")

# grid100 = raster("data/guild2/shp/altitude100.tif")
# buf_ra <- st_rasterize(buf_sf2,st_as_stars(grid100,values=NA_real_))
# write_stars(buf_ra, "data/guild2/shp/riverFGT_raster.tif")



# # Illustration
# A <- buf_sf[buf_sf$riv.id%in%c(17848509),"dZ"]
# A$zclass <- cut(A$dZ,c(-500,-10,10,600))
# riv <- st_read("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/Typisierung_LV03/FGT.shp")
# B <- riv[riv$OBJECTID_G%in%c(17848509),][1]
# mnt <- raster("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/data/mnt/mnt100tot1.tif")
# MNT <- crop(mnt,extent(A))
# 
# # Before filtering
# G <- RStoolbox::ggR(MNT,geom_raster = TRUE)  + 
#   scale_fill_gradient(low = "white", high = "black",name = "elevation") +
#   ggnewscale::new_scale("fill") +
#   geom_sf(data=A,aes(col=zclass)) +  
#   scale_colour_manual("Delta altitude",values=c("red","green","orange")) +
#   geom_sf(data=B,col="blue") 
# ggsave(G,file="graphs/buffer_before_filtering.jpg")
# 
# # After filtering
# AA <- buf_sf2[buf_sf2$riv.id%in%c(17848509),"out"]
# G2 <- RStoolbox::ggR(MNT,geom_raster = TRUE)  + 
#   scale_fill_gradient(low = "white", high = "black",name = "elevation") +
#   geom_sf(data=AA,col="green") +  
#   geom_sf(data=AA[AA$out==2,],col="red") 
# G2
# ggsave(G2,file="graphs/buffer_after_filtering.jpg")


# ## Verification à l'Allandon
# bbox=c(xmin=490840,xmax=491454,ymin=114871,ymax=115334)
# buf_crop <- st_crop(buf_sf,extent(bbox))
# buf_crop$zclass <- cut(buf_crop$dZ,c(-500,-10,10,600))
# mnt <- raster("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/data/mnt/mnt100tot1.tif")
# MNT <- crop(mnt,extent(bbox))
# riv <- st_read("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/Typisierung_LV03/FGT.shp")
# 
# buf_crop$out<- NA
# buf_crop[buf_crop$dZ<=10,]$out <- discard_outliers(st_coordinates(buf_crop[buf_crop$dZ<=10,]))
# plot(buf_crop[7])
# 
# RStoolbox::ggR(MNT,geom_raster = TRUE)  +
#   scale_fill_gradient(low = "white", high = "black",name = "elevation") +
#   geom_sf(data=st_crop(riv,extent(bbox)),col="cyan",size=2) +
#   geom_sf(data=buf_crop,aes(col=zclass)) +
#   scale_colour_manual("Delta altitude",values=c("red","green","orange")) +
#   geom_sf(data=buf_crop[buf_crop$dZ<=10,],shape=1,size=2,stroke=2,col="green")

# 2. Add sinuostiy by cell
load("data/guild2/buffer_river_200m_filtered_02.20.Rdata")
table(duplicated(buf_sf2$grid.id))
head(buf_sf2[duplicated(buf_sf2$grid.id),])
buf_sf2 <- buf_sf2[!duplicated(grid.id)]

load("data/guild2/sinuostiy.Rdata") # data = res
res1 <- lapply(res,function(x) st_set_crs(x,st_crs(buf_sf2)))
coord_riv <- lapply(res1,function(x) data.frame(sinu=x$sinu,slope=x$Slope,debit=x$ABFLUSS,length=x$SHAPE_Leng,st_coordinates(x)))
sinu <-  bind_rows(coord_riv)
sinu2 <- st_as_sf(sinu,coords = c("X", "Y"), crs = 21781)
sinu2$debit <- as.numeric(factor(sinu2$debit,levels=c("klein","mittel","gross","fluss")))

rm(list=c("coord_riv","res","res1","sinu"))

st_geometry(buf_sf2) <- "geometry"
system.time(temp <- st_join(buf_sf2[,"grid.id"],sinu2,join = nngeo::st_nn,k = 1, maxdist = 200))  # 33 sec

## 2.1 a) Stat zonale (moyenne) dans les pixels incluant une riviere, nn pour les autres
stat_by_cell <- setDT(temp)[,.(sinu=mean(sinu,na.rm=T),slope=mean(slope,na.rm=T),debit=mean(debit),length=mean(length)),by="grid.id"]
buf_sf3 <- right_join(buf_sf2, stat_by_cell,by="grid.id")
buf_sf3$length.x <- NULL
table(is.na(buf_sf3$sinu))

# ## 2.1 b) Assign value of sinuosity from nn to grid
temp <- setDT(buf_sf3)[is.na(sinu),]
temp <- within(temp,rm(list=c("sinu", "slope", "debit", "length.y")))
st_geometry(temp) <- "geometry"
system.time(temp <- st_join(temp,sinu2, join = nngeo::st_nn,k = 1, maxdist = 300))  # 13 sec
table(is.na(temp$sinu))
names(buf_sf3) <- names(temp)
BB <- setDT(buf_sf3)[!is.na(sinu)]
st_geometry(BB) <- "geometry"
buf_sf3 <- rbind(BB,temp,deparse.level = 1)
table(is.na(buf_sf3$sinu))

# Assign grid information and save (05.02.20)
# buf_sf3 <- left_join(buf_sf3,riv %>% st_drop_geometry(),by=c("riv.id"="OBJECTID_G"))
# table(duplicated(buf_sf3$grid.id))
# buf_sf3$BASIS_NR <- NULL
# save(buf_sf3,file="data/guild2/river_buffer200m_BV_sf.Rdata")


# # Checking
# table(is.na(buf_sf2$sinu))
# nrow(unique(buf_sf2[is.na(buf_sf2$sinu),"riv.id"]))  # 141 section without information
# mean(buf_sf2[is.na(buf_sf2$sinu),]$length,na.rm=T)  # mean length = 35 m
# plot(buf_sf2[is.na(buf_sf2$sinu),][1])
# head(buf_sf2[is.na(buf_sf2$sinu),])
# MISS.RIV <- unique(buf_sf2[is.na(buf_sf2$sinu),]$riv.id)
# summary(riv[riv$OBJECTID_G%in%MISS.RIV,"SHAPE_Leng"])
# 
# ## Toutes les rivières sont ds riv
# riv <- st_read("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/Typisierung_LV03/FGT.shp")
# which(riv$OBJECTID_G%in%MISS.RIV)
# 
# ## Toutes les rivières st ds riv.pts
# load("data/riv_pts.Rdata") # data = riv
# which(riv$riv.id%in%MISS.RIV)
# riv[riv$riv.id==1620128,]
# 
# ## Toutes les rivières st ds sinu2
# which(sinu2$riv.id%in%MISS.RIV)
# 
# ## Les rivières ne st pas ds test2 (rayon de recherche =200m)
# which(test2$riv.id%in%MISS.RIV)
# test2[test2$riv.id==1522124,]
# riv[riv$riv.id==1522124,]
# sinu2[sinu2$riv.id==1522124,]
# riv[riv$riv.id==16754593,]
# st_distance(sinu2[sinu2$riv.id==1522124,],riv[riv$riv.id%in%c(1522124),])
# 
# ggplot() + 
#   geom_sf(data=riv[riv$riv.id%in%c(16754593,1522124),],col=as.numeric(as.factor(riv[riv$riv.id%in%c(16754593,1522124),]$riv.id))) +
#   geom_sf(data=sinu2[sinu2$riv.id==1522124,],shape=3,col=3)

# # # Visual checking
# RIVpt <- riv_pts_sf[riv_pts_sf$riv.id==14081588,]
# riv <- st_read("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/Typisierung_LV03/FGT.shp")
# RIV <- riv[riv$OBJECTID_G==14081588,]
# riv_pts <- st_line_sample(RIV,density=0.01,type="regular")
# BUF <- st_crop(buf_sf2,extent(RIVpt )+c(-200,200,-200,200))
# riv_around <- st_crop(riv,extent(RIVpt )+c(-200,200,-200,200))
# SIN <- st_crop(sinu2,extent(RIV))
# 
# ggplot() +
#   geom_sf(data=st_geometry(BUF),col=as.numeric(as.factor(BUF$riv_pt))) +
#   geom_text(data=st_geometry(BUF),aes(x=st_coordinates(BUF)[,1],y=st_coordinates(BUF)[,2],label=BUF$riv_pt),vjust=-1,size=3) +
#   geom_sf(data=riv_around,col="lightblue")   +
#   geom_sf(data=RIV,col="blue")   +
#   geom_sf(data=RIVpt,shape=4,size=5,col="green")   +
#   geom_sf(data=SIN,col="blue")
rm(list=c("temp","temp2","sinu2","buf_sf2","stat_by_cell"))

# Add bassin versant:
# Rivers are first assign to BV and then all pixels associated get that BV
riv <- st_read("D:/SIG/2-Eaux.dynamiques/Typisierung_LV03/FGT.shp")
hydro <- st_read("Z:/PROJET Ecological Infrastructure/data/sig_data/decoupage_hydrographiquedelasuisse/Hydrografische+Gliederung/Hydrografische Gliederung_LV03/basis04.shp")
hydro <- st_transform(hydro,st_crs(riv))
riv <- st_join(riv[,"OBJECTID_G"] ,hydro[,"BASIS_NR"],join=st_intersects)
table(duplicated(riv$OBJECTID_G))
head(riv[duplicated(riv$OBJECTID_G),])

# Find the longest section of river into BV and assign this BV to the whole river
RIV_ID <- data.frame(riv=unique(riv[duplicated(riv$OBJECTID_G),]$OBJECTID_G),bv=NA)

for (i in 1:nrow(RIV_ID )) {
  RIV <- riv[riv$OBJECTID_G==RIV_ID[i,1],]
  BV <- hydro %>% filter(BASIS_NR%in%RIV$BASIS_NR)
  int <- st_intersection(RIV, BV)
  RIV_ID[i,2] <- int$BASIS_NR[which.max(st_length(int))]
}

riv <- riv[!duplicated(riv$OBJECTID_G),]
riv$BV04 <- riv$BASIS_NR
riv <- within(riv,BV04[match(RIV_ID$riv,OBJECTID_G)] <- RIV_ID$bv)
head(subset(riv,dif!=0))

# # Add bassin versant 2015 
# bv <- st_read("D:/SIG/data/bioregion_ch/région_biogeo_suisse_12classes.shp")  # BV de Luna 
# bv <- st_transform(bv,st_crs(riv))
# riv <- st_join(riv,bv[,"OBJECTID"],join=st_intersects)

# Add cantons
KANT <- st_read("D:/SIG/data/Limites_cantonales/swissBOUNDARIES3D_1_3_TLM_KANTONSGEBIET.shp")
KANT <- st_transform(KANT,st_crs(riv))
temp <- st_join(hydro, KANT[,c("NAME")], join = st_intersects)
table(duplicated(temp$BASIS_NR))

DUP <- data.frame(bv=temp[duplicated(temp$BASIS_NR),]$BASIS_NR,kant=NA)
for (i in 1:nrow(DUP)) {
  BV <- temp[temp$BASIS_NR==DUP[i,1],]
  KAN <- KANT %>% filter(NAME%in%BV$NAME)
  int <- st_intersection(BV, KAN)
  DUP[i,2] <- as.character(int$NAME[which.max(st_area(int))])
}

temp <- temp[!duplicated(temp$BASIS_NR),]
temp$canton <- as.character(temp$NAME)
temp <- within(temp,canton[match(DUP$bv,BASIS_NR)] <- DUP$kant)
riv <- right_join(riv,temp[,c("BASIS_NR","canton")] %>% st_drop_geometry(),by=c("BV04"="BASIS_NR"),all.x=T)
riv <- riv[!is.na(riv$OBJECTID_G),]

rm(list=c("KANT","hydro","bv","KAN","RIV","RIV_ID","temp","DUP"))

buf_sf3 <- left_join(buf_sf4,riv %>% st_drop_geometry(),by=c("riv.id"="OBJECTID_G"))
table(duplicated(buf_sf3$grid.id))
buf_sf3$BASIS_NR <- NULL
save(buf_sf3,file="data/guild2/river_buffer200m_BV_sf_02.20.Rdata")

# Commut buffer into polygons
buf_sf3$centro <- buf_sf3$geometry
buf <- inner_join(grid_sf[,c("grid.id","geometry")],buf_sf3 %>% st_drop_geometry(),by="grid.id")
save(buf,file="data/buffer_river_200m_filtered_BV_polygons_sf.Rdata")  


## Save as raster
grid <- raster::raster("C:/Dossier_Ervan/R/Grid100/grid100.tif")
st_geometry(buf) <- "geometry"
buf_ra <- fasterize::fasterize(buf,grid,field="grid.id")
raster::writeRaster(buf_ra,"data/G2_buffer_river_07.02.20.tif",format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE) 



# Commute df to raster
library(raster)
load("data/guild2/buffer_river_200m_filtered_with_BV_polygons.Rdata")
env <-stack("C:/Dossier_Ervan/R/infraeco_git/data/guild2/shp/env100.tif")
sinu <- env[[1]]
values(sinu) <- NA
values(sinu)[buf$grid.id] <- buf$sinu

slope <- env[[1]]
values(slope) <- NA
values(slope)[buf$grid.id] <- buf$slope

debit <- env[[1]]
values(debit) <- NA
values(debit)[buf$grid.id] <-  buf$debit

riv.length <- env[[1]]
values(riv.length) <- NA
values(riv.length)[buf$grid.id] <-  buf$length

crs(sinu) <- crs(slope) <-crs(debit) <- crs(riv.length) <- crs(env[[1]])

env2 <-  stack(env,sinu,slope,debit,riv.length,full.names = TRUE)
NAME.VAR <-c('bio19_pcoldq', #precipitation of the coldest quarter
             'bio18_pwarmq', #precipitation of the warmest quarter
             'bio17_pdryq', #precipitation of the dryest quarter
             'bio16_pwetq', #precipitation of the wettest quarter
             'bio15_ps', #precipitation seasonality
             'bio14_pdry', #precipitation of the dryest month
             'bio13_pwet', #precipitation of the dryest month
             'bio12_p', # annual precipitation
             'bio11_tcoldq', #temperature of the coldest quarter
             'bio10_twarmq', #temperature of the warmest quarter
             'bio9_tdryq', #temperature of the dryest quarter
             'bio8_twetq', #temperature of the wettest quarter 
             'bio7_tar', #temperature annual range
             'bio6_tminc', # minimal temperature of the coldest month
             'bio5_tmaxw', # minimal temperature of the warmest month
             'bio4_ts', #temperature seasonality
             'bio3_iso', #isothermality
             'bio2_dr', #mean diurnal range
             'bio1_tmean', #average temperature
             'arridity', # precipitation - potential evapotranspiration
             'gdd3Y_8110_ngb5_mwconic_', #growing degree days
             'sdiryy', # radiations
             'slp25',  # slope
             'topos', # topography (curvature)
             'soilPH', #Soil pH
             'ndviSD', #NDVI variation
             'ndviQ80', # 80th percentile of the yearly NDVI's
             'ndviQ50',# 50th percentile of the yearly NDVI's
             'ndviMIN',# minimum of the yearly NDVI's
             'ndviMAX',# maximum of the yearly NDVI's
             'ndviMEAN',# mean of the yearly NDVI's
             'ForestQ95',# 95th percentile of the tree canopy heights
             'ForestQ25',# 25th percentile of the tree canopy heights
             'alt', # altitude
             'bias',# sampling effort
             'mask.tif', # raster mask
             'sinuosity', # mean sinuosity of river portion in cell (1=straight)
             'slope', # mean slope of river portion in cell (units=?)
             'debit', # mean slope of river portion in cell (units=?)
             'river_length') # mean of river section length in cell (not sure  it is useful)
names(env2) <- NAME.VAR
nlayers(env2)
rm(list=c("sinu","slope","debit","riv.length","env"))
writeRaster(env2,file="data/guild2/shp/env100.with.river.tif",overwrite=TRUE,format="GTiff")


# ## Vérification de la continuité du buffer
# ##########################################
# load("data/guild2/buffer_river_with.alt.Rdata")
# hydro <- st_read("Z:/PROJET Ecological Infrastructure/data/sig_data/decoupage_hydrographiquedelasuisse/Hydrografische+Gliederung/Hydrografische Gliederung_LV03/basis04.shp")
# hydro <- st_transform(hydro,st_crs(buf_sf))
# buf_sf <- st_join(buf_sf,hydro[,"BASIS_NR"],join=st_intersects)
# buf_sf$dZclass <-cut(buf_sf$dZ,c(-500,0,10,15,20,600))
# 
# riv2 <- st_read("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/Typisierung_LV03/FGT.shp")
# riv2 <- st_transform(riv2,st_crs(buf_sf2))
# riv2 <- st_join(riv2,riv[,"BV"],join=st_intersects)
# RIV <- riv2 %>% filter(BV%in%50456)
# 
# ID <- buf_sf[buf_sf$BASIS_NR==50456,]$grid.id
# # load("data/guild2/buffer_river_200m_filtered.Rdata")
# ggplot(buf_sf2 %>% filter(BASIS_NR==50456)) + 
#   geom_sf() + 
#   geom_sf(data=buf_sf2 %>% filter(BASIS_NR==50456 & dZ<=10),col=2) 
# 
# mnt <- raster("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/data/mnt/mnt100tot1.tif")
# MNT <- crop(mnt,extent(buf_sf2 %>% filter(BASIS_NR==50456)))
# 
# # Before filtering
# RStoolbox::ggR(MNT,geom_raster = TRUE)  +
#   scale_fill_gradient(low = "white", high = "black",name = "elevation") +
#   ggnewscale::new_scale("fill") +
#   geom_sf(data=buf_sf2 %>% filter(BASIS_NR==50456),aes(col=dZclass,fill=dZclass)) +
#   guides(fill = guide_legend(override.aes = list(shape = 21))) +
#   geom_sf(data=buf_sf2 %>% filter(grid.id%in%c(6767597,6763921)),col=2,shape=21,stroke=1.5) +
#   geom_sf(data=RIV,col="blue")
# 
# ## After filtering
# RStoolbox::ggR(MNT,geom_raster = TRUE)  +
#   scale_fill_gradient(low = "white", high = "black",name = "elevation") +
#   ggnewscale::new_scale("fill") +
#   geom_sf(data=buf_sf %>% filter(BASIS_NR==50456)) +
#   guides(fill = guide_legend(override.aes = list(shape = 21))) +
#   geom_sf(data=buf_sf %>% filter(BASIS_NR==50456 & dZ>10),col=3) +
#   geom_sf(data=buf_sf %>% filter(grid.id%in%c(6767597,6763921)),col=2,shape=21,stroke=1.5) +
#   geom_sf(data=RIV,col="blue")
# 
#   RStoolbox::ggR(MNT,geom_raster = TRUE)  +
#     scale_fill_gradient(low = "white", high = "black",name = "elevation") +
#     ggnewscale::new_scale("fill") +
#     guides(fill = guide_legend(override.aes = list(shape = 21))) +
#     geom_sf(data=buf_sf %>% filter(BASIS_NR==50456 & dZ<=10)) +
#     geom_sf(data=RIV,col="blue")
#   
# library(tmap)
# tmap_mode("view") 
# tm_shape(buf_sf %>% filter(BASIS_NR==50456)) +
#   tm_dots(
#     col = "dZ",
#     palette = "Reds",
#     breaks=c(-500,0,10,15,20,600),
#     id = "grid.id",
#     showNA = T,
#     colorNA = "gold2",
#     alpha = 0.7)  +
#   tm_shape(RIV) +
#   tm_lines(col="blue") 
# 
# write.csv2(buf_sf %>% filter(grid.id%in%c(6767597,6763921)) %>% st_drop_geometry(),"clipboard")

##########################
### Update plausi 09.02.21
##########################
## Narrow initial buffer to buffer of 50m and remove "ruisselet"
## 1. Create buffer and remove "ruisselets"
library(sf)
riv <- st_read("D:/SIG/2-Eaux.dynamiques/Typisierung_LV03/FGT.shp")
table(riv$ABFLUSS)
riv <- riv[riv$ABFLUSS!="klein",]

buf50 <- st_buffer(riv,50)

## 2. Subset part of buffer overlapping buffer 50m
buf <- LOAD("data/buffer_river_200m_filtered_BV_polygons_sf.Rdata")
st_geometry(buf) <- "centro"
system.time(buf <- st_join(buf,buf50[,"OBJECTID_G"],join=st_intersects) ) # 40 sec
table(is.na(buf$OBJECTID_G))

## 3. Add auen, flachmoor und hochmoor
load("C:/Dossier_Ervan/R/Grid100/grid100_sf_with_enviro.Rdata")
st_geometry(grid_sf) <- "geometry"
system.time(buf <- st_join(buf,grid_sf[,c("au","flachm","hochm")],join=st_intersects) ) # 411 sec
buf <- data.table::setDT(buf)[!duplicated(grid.id),]

## 4. Final selection: include in buffer 50 (OBJECTID_G), in auen, out flach un hochmooren
buf_final <-data.table::setDT(buf)[!is.na(OBJECTID_G)|au==1 & flachm==0 & hochm==0,]
head(buf_final)
st_geometry(buf_final) <- "geometry"

grid <- raster::raster("D://SIG/data/grid100/grid100.tif")
EG_ra <- fasterize::fasterize(buf_final,grid)
raster::writeRaster(EG_ra,paste0("shp/G2_EG_09.02.21.tif"),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)





