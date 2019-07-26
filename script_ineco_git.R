require(sf)
require(data.table)
require(tidyverse)
require(ggplot2)
require(reshape2)
require(qgraph)
require(vegan)
require(raster)

# # Import guildes
# # Import obs (pour guilde de Stefan - 18/07/2019)
obs <- setDT(read.csv("C:/Dossier_Ervan/R/Ineco_Rproject/data/extraction_Leitarten_17072019.csv"))

guild <- setDT(read.csv2("C:/Dossier_Ervan/R/Ineco_Rproject/data/Leitarten_3Probegilden.csv"))

# Select guild and criterion
select <- guild[G2!="",id_Checklist2017]  # select G2
length(select) ### nbr d'espèces

# Subset par milieu
obs2 <- subset(obs,project_taxon_id%in%select) # project_taxon_id = regroupement des sous-espèces sous 1 taxon id.


# Add milieu
obs2 <- merge(obs2,guild[,c(2,3,6)],by.x="project_taxon_id",by.y="id_Checklist2017",all.x=T,allow.cartesian=TRUE)
obs2$G2 <- factor(obs2$G2)

# Add rivers
library(sp)
riv <- st_read("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/Typisierung_LV03/FGT.shp")

bern <- st_read("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/Export_shp/FGT_20_Ge/Cantons_CH1903.shp")
bern <- st_read("S:/09_Hotes/CRSF/TRANSFERTS/ERR/pourMax/Limites_cantonales/swissBOUNDARIES3D_1_3_TLM_KANTONSGEBIET.shp")
bern <- st_transform(bern,21781)

# Clip data from Bern 
obs3 <- obs2 %>% filter(!is.na(x_swiss))
obs_sf <- obs3 %>%  st_as_sf(coords = c("x_swiss","y_swiss"))
st_crs(obs_sf) <- st_crs(bern)
obs_sf <- st_transform(obs_sf,21781)
compareCRS(obs_sf,bern)

# Clip river in Bern
# riv <- st_transform(riv,21781)
# compareCRS(riv,bern)
# system.time(test <- st_join(riv, bern[bern$NAME=="Bern","NAME"], join = st_intersects)) # 60 sec
# table(is.na(test$NAME.y))
# riv <- test[!is.na(test$NAME.y),]

# ## Buffer of rivers
# # 1. Merge zones alluviales et buffer rivières
# library(RQGIS3)
# set_env(root = "C:/Program Files/QGIS 3.4",new=T)
# open_app()
# qgis_session_info()
# find_algorithms(search_term = "union",
#                 name_only = TRUE)
# get_usage(alg = "native:union")
# params <- get_args_man(alg = "native:union")
# 
# za <- st_read("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/data/au.shp")
# 
# out <- run_qgis(alg = "native:union",
#                 INPUT = "C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/data/au.shp",
#                 OVERLAY="C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/Export_shp/buffer_rivières_500m.shp",
#                 OUTPUT="C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/data/buffer_g2.shp")
buf <- st_read("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/data/buffer_g2.shp")
buf <- st_transform(buf,21781)
system.time(test <- st_join(buf, bern[bern$NAME=="Bern","NAME"], join = st_intersects)) # 60 sec
table(is.na(test$NAME))
buf <- test[!is.na(test$NAME),]

system.time(test <- st_join(obs_sf[obs_sf$v_xy_radius<=100,], buf[,"SHAPE_Le_2"], join = st_intersects)) # 6 sec
table(is.na(test$SHAPE_Le_2))
obs_sf <- test[!is.na(test$SHAPE_Le_2),]

# Sinuosity of neareast river portion

load("C:/Dossier_Ervan/R/Ineco_Rproject/data/sinuostiy.Rdata") # data =res

# 2. Crop buffer where elevation is higher than 10 m
env <-stack("C:/Dossier_Ervan/R/Ineco_Rproject/data/shp/env100.tif")
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
              'mask.tif') # raster mask
names(env) <- NAME.VAR
alt_ra <- subset(env,"alt")
# writeRaster(alt_ra,file="C:/Dossier_Ervan/R/Ineco_Rproject/data/shp/altitude100.tif")

# Transform raster into points (extract centroïd) 
centro_ra <- rasterToPoints(alt_ra, spatial = TRUE)

# Transform raster into sf 
library(stars)  ## see demo here : https://r-spatial.github.io/stars/articles/stars5.html
x = read_stars("C:/Dossier_Ervan/R/Ineco_Rproject/data/shp/altitude100.tif")
# env_sf <- st_as_sf(x, as_points = TRUE, merge = FALSE)
# names(env_sf) <- c(NAME.VAR,"geometry")
alt_sf_pts <- st_as_sf(st_as_stars(alt_ra), as_points = T, merge = FALSE)

# Crop raster to Bern
alt_sf_pts <- st_transform(alt_sf_pts,21781) # time consuming
system.time(test <- st_join(alt_sf_pts, bern[bern$NAME=="Bern","NAME"], join = st_intersects)) # 6 sec
table(is.na(test$NAME))
alt_sf <- test[!is.na(test$NAME),]

# Compute min distance to river and assign altitud of nearest river
alt_sf$id <- 1:nrow(alt_sf)

# riv_sp <- as(riv, 'Spatial')
# dist <- rgeos::gProject(riv_sp,centro_ra[1,])
centro_sf <- st_as_sf(centro_ra)

# Crop raster topo to buffer
buf <- st_transform(buf,st_crs(centro_sf))
system.time(test <- st_join(centro_sf,buf[,"OBJECTID_G"],join = st_within,prepared=F))  # 15 minutes

# Remove duplicates
centro_crop <- test[!duplicated(test[c(1,3)]),]
centro_sf_crop <- centro_crop[!is.na(centro_crop$OBJECTID_G),] # crop to area in buffer


assign.alt <-  function(x) {
  GEOM <- st_geometry(x$geometry)
  #riv <-  get("riv", envir = .GlobalEnv)
  st_crs(GEOM) <- st_crs(riv)
  nn.river <- which.min(st_distance(GEOM,riv))
  
  dist2river <- as.numeric(st_distance(st_cast(riv[nn.river,], "POINT"), GEOM))
  if(any(dist2river<500)){
  pts <- st_cast(riv[nn.river,], "POINT")
  pts <- st_transform(pts,st_crs(centro_sf))
  alt_river <- extract(alt_ra,pts)
  
  alt_point <- alt_river[which.min(dist2river)] 
  
  } else {
  alt_point <- NA
  }
  alt_point
}

assign.alt2 <-  function(x) {
    GEOM <- st_geometry(x$geometry)
    RIV <- riv[riv$OBJECTID_G==x$OBJECTID_G,]
    st_crs(GEOM) <- st_crs(RIV)
    #nn.river <- which.min(st_distance(GEOM,RIV))
    npoints <-  unique(RIV$SHAPE_Leng)/40
    POINTS <- st_sample(RIV,npoints,type = "regular")
    PP <- st_cast(POINTS,"POINT")  # split into single points
    dist2river <- as.numeric(st_distance(PP, GEOM))
    pts <- as(PP[which.min(dist2river)],"Spatial")
    crs(pts) <- crs(alt_ra)
    alt_ra_crop <- crop(alt_ra,extent(RIV))
    alt.riv <- raster::extract(alt_ra_crop,pts)  
    alt.riv
}

centro_crop$alt <- centro_sf$alt
DF <- centro_crop[500:1000,]
DF <- DF[!is.na(DF$OBJECTID_G),]
DF$alt.riv <- NA

f4 <- function(DF) {
  for (i in 1:nrow(DF))
  GEOM <- st_geometry(DF[i,]$geometry)
  RIV <- riv[riv$OBJECTID_G==DF[i,]$OBJECTID_G,]
  GEOM <- st_transform(GEOM,st_crs(RIV))
  #nn.river <- which.min(st_distance(GEOM,RIV))
  npoints <-  unique(RIV$SHAPE_Leng)/40
  POINTS <- st_sample(RIV,npoints,type = "regular")
  PP <- st_cast(POINTS,"POINT")  # split into single points
  
  plot(RIV[1])
  P <- st_cast(RIV[1],"POINT")
  st_crs(P) <- st_crs(RIV[1])
  plot(P,add=T,col=2)
  plot(PP,add=T,col=1)
  
  dist2river <- as.numeric(st_distance(PP, GEOM))
  pts <- as(PP[which.min(dist2river)],"Spatial")
  crs(pts) <- crs(alt_ra)
  alt_ra_crop <- crop(alt_ra,extent(RIV))
  DF[i,"alt.riv"] <- raster::extract(alt_ra_crop,pts)
  }
  
system.time(A <- f4(DF))
system.time(DF[1:100,"alt.riv"] <- apply(DF[1:100,],1, function(x) assign.alt2(x)))

library(microbenchmark)


library(parallel)
ncores<-detectCores()-1
if (ncores == 0) {ncores = 1}
cl <- makeCluster(ncores)
clusterExport(cl,c("alt_ra","alt_sf","assign.alt","f4_par","centro_sf_crop","riv"),envir=.GlobalEnv)
clusterEvalQ(cl,{
  library(sf)
  library(raster)
})

centro_crop$alt <- centro_sf$alt
DF <- centro_crop[500:1000,]
DF <- DF[!is.na(DF$OBJECTID_G),]
DF$alt.riv <- NA

system.time(dist <- parApply(cl,DF, 1, function(x) assign.alt2(x)))

system.time(dist <- parApply(cl,DF, 1, function(x) assign.alt2(x)))
parallel::stopCluster(cl)


# Shannon index by grid cell
in_grid <- st_join(obs_sf, grid_sf, join = st_within)
shannon <- in_grid %>% group_by(id) %>% mutate(div=diversity(as.numeric(table(project_taxon_id))),N=length(unique(project_taxon_id)),indic=ifelse(any(G2=="Z!"),1,0))

# AA <- as.data.table(shannon)
# summary(AA[!is.na(id),div])
# A <- tapply(AA$code,AA$id,length)

subset_grid <- grid_sf[grid_sf$id%in%unique(shannon$id),]
div_in_grid <- unique(left_join(subset_grid, as.data.frame(shannon[,c("id","div","N","indic"),drop=T]),by="id"))
st_write(div_in_grid, "data/shp/guild2_leitarten_in.river.buffer.shp")


library(tmap)
tmap_mode("map")
  tm_shape(riv[,1]) + 
  tm_lines() +
  tm_shape(div_in_grid) +
    tm_fill(
    col = "div",
    palette = "Greens",
    style = "cont",
    contrast = c(0.1, 1),
    title = "Shannon index",
    id = "id",
    showNA = FALSE,
    alpha = 1,
    popup.vars = c(
      "Total obs" = "div"),
    popup.format = list(
      div = list(format = "f", digits = 0)
    )) +
  tm_borders(col = "gray", lwd = 0.2) +
  tm_layout(title= "Eaux dynamiques") +
  tmap_options(basemaps.alpha=0.5) 
  


# View species on boarder of lakes
temp <- read.csv("C:/Dossier_Ervan/SIG/bord_lac_NE.csv")  
list_sp <- as.character(in_grid[in_grid$id%in%temp$id,]$taxon_name)
tab <- data.frame(name=names(table(list_sp)),n=as.numeric(table(list_sp)))
tab <- tab[order(tab$n,decreasing = T),]
head(tab,50)

par(mar=c(12,3,0.2,0))
barplot(tab[1:18,]$n,names.arg=tab[1:18,]$name,las=2,cex.names=0.8,xpd=T)




# Différence (différence d'abondance avant/après 2000)
obs1900 <- obs[!is.na(x_swiss) & milieu_fr=="prairies humides eutrophes" & obs_year<2000]
obs2000 <- obs[!is.na(x_swiss) & milieu_fr=="prairies humides eutrophes" & obs_year>=2000 ]

vars <- c("x_swiss", "y_swiss")

obs1900_sf <- obs1900 %>% mutate_at(vars(x_swiss, y_swiss), as.numeric) %>% 
  st_as_sf(coords = c("x_swiss", "y_swiss")  )
st_crs(obs1900_sf) <- st_crs(obs_sf)
in1900 <- st_join(obs1900_sf, grid_sf, join = st_within)
count1900 <- count(as_tibble(in1900), id) %>%
  print()

obs2000_sf <- obs2000 %>% mutate_at(vars(x_swiss, y_swiss), as.numeric) %>% 
  st_as_sf(coords = c("x_swiss", "y_swiss")  )
st_crs(obs2000_sf) <- st_crs(obs_sf)
in2000 <- st_join(obs2000_sf, grid_sf, join = st_within)
count2000 <- count(as_tibble(in2000), id) %>%
  print()

AA <- left_join(count1900,count2000,by = c("id")) %>% mutate(diff=n.y-n.x) %>% filter(!is.na(id))
AA %>% summarise(min=min(diff,na.rm=T),max=max(diff,na.rm=T),mean=mean(diff,na.rm=T))

tmap_mode("plot")
A <- tm_shape(plant2) +
  tm_fill(
    col = "n",
    palette = "Greens",
    style = "cont",
    contrast = c(0.1, 1),
    title = "Densité de plantes indicatrices (obs/ha) ",
    id = "id",
    showNA = FALSE,
    alpha = 0.8,
    popup.vars = c(
      "Total obs" = "n"
    ),
    popup.format = list(
      n = list(format = "f", digits = 0)
    )) +
  tm_borders(col = "gray", lwd = 0.2) +
  tm_layout(title= "Prairies humides eutrophes")


# # Comparaison des guildes
# BB <- gui3[,name,by=milieu_fr]
# BBB <- dcast(BB,milieu_fr~name)
# dis <- as.matrix(vegdist(decostand(BBB[,-1],"norm"), "euclid",binary=T,diag=T, upper=T))
# 
# #dis[lower.tri(dis, diag = FALSE)] <- NA
# GROUP <- data.frame(id=letters[1:18],name=as.character(as.factor(BBB[,1])))
# GROUP$leg <- paste(GROUP$id,GROUP$name,sep=". ")
# dist_mi <- 1/dis # one over, as qgraph takes similarity matrices as input
# rownames(dist_mi) <- GROUP[,1]
# colnames(dist_mi) <- GROUP[,1]
# 
# # Flux graph 
# png('overlap_guilds.png', width=50, height=20, res=500,unit='cm')
# qgraph(dist_mi, layout='spring', vsize=3,esize=3,groups=GROUP[,3],legend = TRUE, borders = FALSE)
# dev.off()

## Add bassin versant
# hydro <- st_read("Z:/PROJET Ecological Infrastructure/data/sig_data/decoupage_hydrographiquedelasuisse/Hydrografische+Gliederung/Hydrografische Gliederung_LV03/basis04.shp")
# names(hydro)
# print(hydro[1:5], n = 3)
# #plot(hydro[5])
# bassin <- st_geometry(hydro)
# plot(bassin)
# 
# df <- st_sf(id = 1:length(st_geometry(hydro)), area_km2 = st_area(hydro)/10^6 , geometry = st_geometry(hydro))
# class(df)
# 
# Convert obs into sf object
library(data.table)
summary(obs)
obs[,carac := as.numeric(0)]
obs <- within(obs,carac[characteristic==1] <- 1)
obs <- within(obs,carac[characteristic==0 & dominant==1] <- 1)

# obs2 <- obs %>% filter(!is.na(x_swiss) & milieu_fr=="prairies humides eutrophes" & carac==1)
# vars <- c("x_swiss", "y_swiss")
# 
# obs_sf <- obs2 %>% mutate_at(vars(x_swiss, y_swiss), as.numeric) %>% 
#   st_as_sf(coords = c("x_swiss", "y_swiss")  )
# 
# st_crs(obs_sf) <- st_crs(hydro)
# 
# # find points within polygons
# plants_in_bassin <- st_join(obs_sf, df, join = st_within)
# 
# # count plant per bassin
# plant_count <- count(as_tibble(plants_in_bassin), id) %>%
#   print()
# 
# plant_density <- left_join(df, plant_count) %>%
#   mutate(density = as.numeric(n / area_km2))  %>%
#   print()
# 
# library(tmap)
# tmap_mode("view")
# A <- tm_shape(plant_density) +
#   tm_fill(
#     col = "density",
#     palette = "Greens",
#     style = "cont",
#     contrast = c(0.1, 1),
#     title = "Densité (obs./km2) de plantes indicatrices",
#     id = "id",
#     showNA = FALSE,
#     alpha = 0.8,
#     popup.vars = c(
#       "Total obs" = "n",
#       "Obs./km2" = "density"
#     ),
#     popup.format = list(
#       n = list(format = "f", digits = 0),
#       density = list(format = "f", digits = 1)
#     )) +
#   tm_borders(col = "gray", lwd = 0.2) +
#   tm_layout(title= "Prairies humides eutrophes")
# 
# tmap_save(A,filename="maps/Prairies humides eutrophes.png",width=20,height=18,unit="cm")
# 

