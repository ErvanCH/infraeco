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
obs <- setDT(read.csv("C:/Dossier_Ervan/R/infraeco_git/data/guild2/extraction_Leitarten_17072019.csv"))

guild <- setDT(read.csv2("C:/Dossier_Ervan/R/infraeco_git/data/guild2/Leitarten_3Probegilden.csv"))

# Select guild and criterion
select <- guild[G2!="",id_Checklist2017]  # select G2
length(select) ### nbr d'espèces

# Subset par milieu
obs2 <- subset(obs,project_taxon_id%in%select) # project_taxon_id = regroupement des sous-espèces sous 1 taxon id.


# Add milieu
obs2 <- merge(obs2,guild[,c(2,3,6)],by.x="project_taxon_id",by.y="id_Checklist2017",all.x=T,allow.cartesian=TRUE)
obs2$G2 <- factor(obs2$G2)
obs3 <- obs2 %>% filter(!is.na(x_swiss))
obs_sf <- obs3 %>%  st_as_sf(coords = c("x_swiss","y_swiss"))
rm(list=c("obs2","obs3","obs"))

# Subset obs in river buffer
buf <- raster("C:/Dossier_Ervan/R/infraeco_git/data/buffer_river_500m_filtered_za.tif")
library(stars)
buf_sf <- st_as_sf(st_as_stars(buf), as_points = F, merge = FALSE)
names(buf_sf)[1] <- "alt"
st_crs(obs_sf) <- st_crs(buf_sf)

system.time(test <- st_join(obs_sf[obs_sf$v_xy_radius<=100,], buf_sf, join = st_intersects)) # 22 sec
table(is.na(test$alt))
obs_sf <- test[!is.na(test$alt),]
rm(list=c("buf","test"))

# Load environmental variables
## Add river characteristics and crop to river buffer
load("C:/Dossier_Ervan/R/infraeco_git/data/guild2/sinuostiy.Rdata") # data =res
system.time(res2 <- lapply(res, function(x) st_sf(sinu=x$sinu,slope=x$Slope,lenght=x$SHAPE_Leng,geom=st_geometry(x)))) # 2min
rm(list=c("res","riv"))

# Solve issue for small river (<50m) with geometry = LINESTRING
test <- lapply(res2,function(x) dim(st_coordinates(x$geom))[2])
res2[test==3] <- lapply(res2[test==3], function(x){
  x$geom <-  st_centroid(st_geometry(x$geom))
  x
})

# rasterize sinuosity & slope
grid <- buf
values(grid) <- 1:ncell(grid)
system.time(res2 <- mapply(cbind,res2,lapply(res2, function(x) extract(grid,st_coordinates(x$geom))),SIMPLIFY=F))
# 145 sec

summary_list <- function(x)  data.frame(cell_id=aggregate(x[1][[1]],list(x[4][[1]]),mean,na.rm=T)[,1],do.call(cbind,lapply(x[1:3],function(y) aggregate(y,list(x[4][[1]]),mean,na.rm=T)[,2])))[,1:4]
system.time(test <- lapply(res2,summary_list))  # 30 min
rm(res2)
system.time(test2 <- do.call(rbind,test))  # 8 sec  ## commute list to df

# Commute df to raster
sinu <- buf
values(sinu) <- NA
values(sinu)[test2$cell_id] <- test2$sinu

slope <- buf
values(slope) <- NA
values(slope)[test2$cell_id] <- test2$slope

riv.length <- buf
values(riv.length) <- NA
values(riv.length)[test2$cell_id] <- test2$lenght

crs(env[[1]])
crs(sinu) <- crs(env[[1]])
crs(slope) <- crs(env[[1]])
crs(riv.length) <- crs(env[[1]])
env2 <-  addLayer(env,sinu,slope,riv.length)
save(env2,file="C:/Dossier_Ervan/R/infraeco_git/data/guild2/shp/env100_river.tif")

env <-stack("C:/Dossier_Ervan/R/infraeco_git/data/guild2/shp/env100.tif")
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
             'river_length') # mean of river section length in cell (not sure  it is useful)
names(env) <- NAME.VAR

# Statistics by cells:
# 1. Count obs per cell
obs_in_cells <- st_join(obs_sf, buf_sf[,c("id")], join = st_within) # find obs in polygons
plant_count <- count(as_tibble(obs_in_cells), id) 

div_in_cells <- obs_in_cells %>% group_by(id) %>% summarise(div=diversity(as.numeric(table(v_accepted_taxon_id))),N=length(unique(v_accepted_taxon_id)),indic=ifelse(any(G2=="Z!"),1,0))
summary(div_in_cells$div)
centro <- st_centroid(div_in_cells)

# Subset by canton
KANT <- st_read("C:/Dossier_Ervan/SIG/data/Limites_cantonales/swissBOUNDARIES3D_1_3_TLM_KANTONSGEBIET.shp")
st_crs(KANT) <- st_crs(KANT,epsg=21781)
obs_sf <- st_transform(obs_sf,st_crs(KANT))
obs_sf <- st_join(obs_sf, KANT[,c("NAME")], join = st_intersects)
table(is.na(obs_sf$NAME))

library(SSDM)
obs <- data.frame(taxon_id=obs_sf[obs_sf$NAME=="Genève",]$v_accepted_taxon_id,st_coordinates(obs_sf[obs_sf$NAME=="Genève",]))

## Need to mask environment with buffer and add sinuosity..

SDM <- stack_modelling(c('CTA', 'SVM'), obs, env, rep = 1, ensemble.thresh = 0,
                       Xcol = 'X', Ycol = 'Y',
                       Spcol = 'taxon_id', method = "pSSDM", verbose = FALSE)

plot(SDM@diversity.map)

# 2. Sinuosity by cell
load("C:/Dossier_Ervan/R/infraeco_git/data/guild2/sinuostiy.Rdata") # data = res
res1 <- lapply(res,function(x) st_set_crs(x,st_crs(alt_sf)))
coord_riv <- lapply(res1,function(x) data.frame(sinu=x$sinu,slope=x$Slope,debit=x$ABFLUSS,length=x$SHAPE_Leng,st_coordinates(x)))
test <-  bind_rows(coord_riv, .id = "column_label")
test2 <- st_as_sf(test,coords = c("X", "Y"), crs = 21781)

riv_in_cells <- st_join(test2, alt_sf[,c("id")], join = st_within)
riv_in_cells <- riv_in_cells[!is.na(riv_in_cells$id),] # subset to bern area
riv_in_cells$debit_num <- as.numeric(as.factor(riv_in_cells$debit))

stat_by_cell <- setDT(riv_in_cells)[,.(sinu=mean(sinu,na.rm=T),slope=mean(slope,na.rm=T),debit=mean(debit_num),long=mean(length)),by=id]
                                                          
riv_ra <- left_join(alt_sf, stat_by_cell,by="id")
riv_ra <- riv_ra[!is.na(riv_ra$sinu),]

#create grid of 1 ha
temp <- plant_density %>%  filter(density==max(density,na.rm=T))
grid <- temp %>% st_make_grid(cellsize = 100, what = "polygons") %>% st_intersection(temp)

grid_sf <- st_sf(id = 1:length(st_geometry(grid)),geometry = st_geometry(grid))
class(grid_sf)
print(grid_sf, n = 3)
st_crs(grid_sf)

# find points within polygons
in_grid <- st_join(obs_sf, grid_sf, join = st_within)
setDT(in_grid)[id==3712,table(taxon_id)]

# Shannon index by grid cell
in_grid <- in_grid %>% group_by(id) %>% mutate(div=diversity(as.numeric(table(taxon_id))))
AA <- as.data.table(in_grid)
head(AA[!is.na(id),])

plant2 <- left_join(grid_sf, as.data.frame(in_grid %>% select(id,div))[1:2],by="id") %>%
  print()
# # count plant per bassin
# plant_count <- count(as_tibble(plants_in_grid), id) %>%
#   print()
# 
abondance <- data.frame(a=c(200,0,0,0,0,0,151,0,0,0,0,0),b=c(60,40,20,5,10,2,30,5,15,2,12,150))
diversity(t(abondance))

tmap_mode("view")
tm_shape(plant2) +
  tm_fill(
    col = "div",
    palette = "Greens",
    style = "cont",
    contrast = c(0.1, 1),
    title = "Shannon index",
    id = "id",
    showNA = FALSE,
    alpha = 0.8,
    popup.vars = c(
      "Total obs" = "div"
    ),
    popup.format = list(
      div = list(format = "f", digits = 0)
    )) +
  tm_borders(col = "gray", lwd = 0.2) +
  tm_layout(title= "Prairies humides eutrophes")


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

