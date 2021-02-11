library(rgdal)
library(raster)
library(ade4)
library(adehabitatHR)
library(adehabitatMA)
library(ecospat)
library(sf)
library(data.table)
library(tidyverse)

ie.dir<-"C:/Dossier_Ervan/R/G2-Eaux_dynamiques"
load("C:/Dossier_Ervan/R/Grid100/grid100_sf_with_enviro.Rdata")
source('~/R/functions_ER.R', echo=F)

## Definit l'espace guilde
buf <- LOAD("data/river_buffer200m_BV_sf_02.20.Rdata")

# Selection des variables et espace-guilde
VAR <- c("bio12_p", "bio6_tminc","bio4_ts", "gdd3Y_8110_ngb5_mwconic_", "topos","ndviMEAN","alt","sinuosity", "slope", "debit")

test <- merge(setDT(buf[,c("grid.id","sinu","slope","debit")]),setDT(grid_sf)[,c("grid.id","subreg","BV04",..VAR[1:6])],by="grid.id",all.x=T)

test$bioregion <- test$subreg
test$subreg <- NULL
test$geometry <- NULL
test <- na.omit(test)
rm(list=c("AXES", "bav", "bio", "PPS","obs","div", "ID","BIOR", "grid_sf", "ie.dir", "mask",  "R"))
.rs.restartR() 


BIOR <- split(test,test$bioregion)
AXES <-data.frame(a=c(1,1,2),b=c(2,3,3))
R = 100
mask <- ascgen(SpatialPoints(cbind((0:(R))/R, (0:(R)/R))),
               nrcol = R - 2, count = FALSE)

library(rgdal)
library(raster)
library(ade4)
library(adehabitatHR)
library(adehabitatMA)
library(ecospat)
library(sf)
library(data.table)
library(tidyverse)

source('~/R/functions_ER.R')
system.time(L <- lapply(BIOR, function(x) pca.rasterize(x,axes=AXES,mask=MASK,R=100,z.th=0.01,save_raster=T))) # 17 min
save(L,file="data/G2_comp_bv_06.02.20.Rdata")

# library(parallel)
# ncores<-detectCores()-1
# if (ncores == 0) {ncores = 1}
# cl <- makeCluster(ncores)
# clusterExport(cl,varlist=c("BIOREG","AXES","mask","rasterizeBV","env.overlap","pca.rasterize"),envir=.GlobalEnv)
# clusterEvalQ(cl,{
#   library(rgdal)
#   library(raster)
#   library(ade4)
#   library(adehabitatHR)
#   library(adehabitatMA)
#   library(ecospat)
# })
# system.time(test<- parLapply(cl,BIOREG, function(x) pca.rasterize(x)))
# stopCluster(cl)  


## Clustering par bioregion

# Clustering
CLUST <- function(x) {
  set.seed(10)
  L <- setDT(x)[,c("BV1","BV2","D")]
  LL <- unique(L$BV1)
  MT <- matrix(0, length(LL),length(LL))
  colnames(MT) <- rownames(MT) <- LL
  for (i in 1:length(LL)) {
    b <- match(unique(L$BV1)[i],LL)
    IDX <- which(unique(L$BV1)[i]==L$BV1)
    MT[b,match(L$BV2[IDX],LL)] <- L$D[which(unique(L$BV1)[i]==L$BV1)]
  }
  MT[lower.tri(MT)] <- MT[upper.tri(MT)]
  diag(MT) <- 1
  
  CLUST <- kmeans(MT,centers=length(LL)/20)
  CLU <- data.frame(BV04=LL,clust=CLUST$cluster)
  CLU
}  

test2 <- lapply(L,CLUST)
test <- rbindlist(test2)
test$BASIS_NR <- as.numeric(substr(test$BV04,2,8))
test$BV04 <- as.numeric(substr(test$BV04,2,8))
test$bioreg <- rep(1:10,as.numeric(sapply(test2,nrow)))
save(test,file="data/G2_clusterBV_06.02.20.Rdata")

# ## Assign BV to a single bioregion (BV are duplicated on boundaries)
# bav <- st_read("Z:/PROJET Ecological Infrastructure/data/sig_data/decoupage_hydrographiquedelasuisse/Hydrografische+Gliederung/Hydrografische Gliederung_LV03/basis04.shp")  # 1063 unique BV
# st_geometry(bav) <- "geometry"
# 
# 
# bio <- sf::st_read("D:/SIG/data/bioregion_ch/région_biogeo_suisse_6classes.shp")
# st_geometry(bio) <- "geometry"
# bio <- st_transform(bio,st_crs(bav))
# 
# bv <- st_join(bav,bio[,"REGBIOG_C6"],st_intersects)
# table(duplicated(bv$BASIS_NR))
# ID <- bv[duplicated(bv$BASIS_NR),]$BASIS_NR
# test <- st_intersection(bv, bio[,"REGBIOG_C6"])
# test$area <- st_area(test$geometry)
# 
# REG <- setDT(test)[,.("REGBIOG_C6"=REGBIOG_C6.1[which.max(area)]),by=BASIS_NR]  # assing region with max overlap 
# BV <- merge(test4,REG,by="BASIS_NR",all.x=T)
# 
# library(dplyr)
# head(BV[BV$BV04%in%BV[duplicated(BV$BV04),]$BV04,][order(BV$BV04),])
# BV$dif <- BV$bioreg - BV$REGBIOG_C6
# BV <- BV[BV$dif==0,]
# BV <- unique(right_join(bav[,c("BASIS_NR","geometry")],BV,by="BASIS_NR"))
# 
# DUP <- BV[BV$BV04%in%BV[duplicated(BV$BV04),]$BV04 & !is.na(BV$REGBIOG_C6),]
# st_geometry(DUP) <- "geometry"
# ggplot() + geom_sf(data=DUP,col=2) + geom_sf(data=bio,fill=NA)
# 
# B <- DUP[!duplicated(DUP$BV04),]
# st_geometry(B) <- st_geometry(DUP %>% 
#   group_by(BV04) %>%
#   summarise(geometry = sf::st_combine(geometry)) )
# # BB <- lapply(st_geometry(B), function(x) st_cast(x, "POLYGON"))
# # st_geometry(B) <- st_sfc(BB)
# 
# st_geometry(BV) <- "geometry"
# 
# BV2 <- rbind(BV[!BV$BV04%in%BV[duplicated(BV$BV04),]$BV04  & !is.na(BV$REGBIOG_C6),],B)
# BV2$dif <- NULL
# BV2$bioreg <- NULL
# BV2$color <- as.factor(paste(BV2$REGBIOG_C6,BV2$clust,sep="_"))
# BV2 <- within(BV2,color[is.na(BV2$clust)] <- NA)
# 
# tapply(BV2$clust,factor(BV2$REGBIOG_C6),function(x) length(unique(x)))
# 
# mypal <- colorRampPalette(RColorBrewer::brewer.pal(5, "PuBu"))
# mypal2 <- colorRampPalette(RColorBrewer::brewer.pal(5, "OrRd"))
# mypal3 <- colorRampPalette(RColorBrewer::brewer.pal(9, "Purples"))
# mypal4 <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlGn"))
# YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
# mypal5 <- colorRampPalette(YlOrBr)
# mypal6 <- colorRampPalette(RColorBrewer::brewer.pal(9, "Greys"))
# 
# G <- ggplot() + 
#   geom_sf(data=bav,col=NA,fill="cyan") +
#   geom_sf(data=BV2,aes(fill=color)) + 
#   scale_fill_manual(values=c(mypal(7)[3:7],mypal2(17)[4:17],mypal3(14),mypal4(5),mypal5(11)[3:11],mypal6(5))) +  
#   theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_rect(fill = "lightgray")) +
#   geom_sf(data=bio,aes(fill=NA,color=as.factor(REGBIOG_C6)),size=1) +
#   scale_color_manual(values=1:6) 
# plot(G)
# 
# ggsave(G,file="graphs/G14_BV_04.10.2019.png",width=14,height=10,dpi=600)

# library(ggplot2)
# ggplot() + 
#  geom_sf(data=bv,aes(fill=REGBIOG_C6),col=NA) +
#   geom_sf(data=bv[is.na(bv$clust),],col="cyan",fill=NA)+
#   geom_sf(data=bio,col="yellow",fill=NA)
#    
# library(tmap)
# tmap_mode("view")
# tm_shape(bav)+
#   tm_fill(col="grey") +
#   tm_borders(col = "gray", lwd = 0.2) +
#   tm_shape(BV2) +
#   tm_fill(col="color")


