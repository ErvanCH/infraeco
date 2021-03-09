options("scipen"=100, "digits"=4)
# install.packages("pacman")
pacman::p_load(sf, data.table,tidyverse,ggplot2,reshape2,qgraph,vegan,raster,hablar)
source("code/functions_ER.R")
# LOCATION <- "G://R/"# location of all folders (PC Infoflora)
LOCATION <- "C://Dossier_Ervan/R/"# location of all folders (IF office)

GUILD <- 3

### Select guild
G <- c(1:4,10:13,17,19,22)


for (GUILD in G) {
  dir <- paste0(LOCATION,list.files(LOCATION)[grep(paste0("G",GUILD,"-"),list.files(LOCATION))])
  guild <- readxl::read_xlsx("data/Guildes_revues_ER.xlsx",sheet=1,col_names = T,col_types = "guess")
  guild.info <- guild[which(guild$No==GUILD),]
  guild.info <- guild.info %>% hablar::convert(lgl(RF,GAM,GBM,GBM_xg,ME))
  DATE <- format(Sys.time(), '%d-%m-%y')
  
  # Assign spatial info & enviro variables
  TT <- readxl::read_xlsx("data/var_enviro_by_guild.xlsx")
  VAR <- TT$label[!is.na(TT[,stringr::str_which(names(TT),paste0("G",GUILD,"$"))])]
  
  ## Load IST Infofauna & attach grid info
  load("data/grid100_sf_with_enviro.Rdata")
  FF <- list.files("data/infofauna")[grep(paste0("_",GUILD,"_bio-idx"),list.files("data/infofauna"))]
  ist <- read.csv(paste("data/infofauna", FF,sep="/"))
  IST <- merge(setDT(ist),setDT(grid_sf)[,c("CNHA","grid.id","geometry")],by="CNHA",all.x=T)
  st_geometry(IST) <- "geometry" 


  ### D?sagr?gation des donn?es
  if (is.na(guild.info$IST_defrag)){
    cat("!! Pas de paramètres de défragmentation !!")
  }
  system.time(IST2 <- IST.dsg(IST,guild.info$IST_defrag))


  ## Data for calibration
    ## A. Observations
  div2 <- setDT(grid_sf)[grid.id%in%IST2$grid.id,c("grid.id","BV04","canton","subreg","geometry","centro",..VAR)]
  sf::st_geometry(div2) <- "centro"
  
  Q1 <- na.omit(data.table::setDT(div2)) 
  nrow(IST2)-nrow(Q1) # 18 ha vires
  Q1$Qobs=1

  ## B. Pseudo-absence dans tout l'espace guilde
  eg <- LOAD("EG","data",dir)
  EG <-setDT(grid_sf)[grid.id%in%c(IST$grid.id,eg$grid.id),c("grid.id","BV04","canton","subreg","geometry","centro",..VAR)]
  EG <- EG[!duplicated(EG$grid.id),]
  EG1 <- setDT(EG)[,c("grid.id","CNHA","geometry")]
  save(EG1,file=paste0("/data/",GUILD,"_EG.RData"))
  
  NArep(EG)
  summary(EG)
  Q0 <- dplyr::sample_n(na.omit(data.table::setDT(EG)[!grid.id%in%Q1$grid.id]),nrow(Q1),replace=F)
  Q0$Qobs=0
  OBS <- rbind(Q1,Q0)

  # ## Selection des variables
  # library(rgdal)
  # library(rgeos)
  # library(parallel)
  # library(raster)
  # library(perm)
  # 
  # ncores<-detectCores()-1
  # if (ncores == 0) {ncores = 1}
  # 
  # env.pres<-as.data.frame(OBS[OBS$Qobs==1,..VAR])#env data for presences
  # env.bck<-as.data.frame(OBS[OBS$Qobs==0,..VAR])#env data for background
  # to.do<-1:ncol(env.pres) #index for parallelization
  # 
  # cl <- makeCluster(ncores)
  # clusterExport(cl,varlist=c('to.do','env.bck','env.pres','var.test'),envir=.GlobalEnv)
  # clusterEvalQ(cl,{
  #   library(perm)
  # })
  # 
  # #apply non-parametric t-test between background and presences for each variable
  # system.time(VAR.test<-parLapply(cl,X=to.do,fun=var.test,Nrep=5,mcmc=1000,env.pres,env.bck))  # 10 sec
  # stopCluster(cl)
  # 
  # #results formatting
  # select.var <- data.frame(var=VAR,do.call(rbind,VAR.test))
  # write.table(select.var, file = "data/selection_var.txt",sep='\t',quote=F,row.names = F)
  # 
  # ## Generate report for variable selection
  # rmarkdown::render(input = paste0(getwd(),"/code/var_select.Rmd"), 
  #                   output_format = "html_document",
  #                   output_file = paste0(GUILD,"_var_select2.html"),
  #                   output_dir = paste0(getwd(),"/report"))
  

  ## Espace de projection
  TT2 <- readxl::read_xlsx("data/var_enviro_selected.xlsx")
  VAR <- TT2$label[!is.na(TT2[,stringr::str_which(names(TT2),paste0("G",GUILD,"$"))])]
  
  ## Projection dans l'espace-guilde (to lower computational time)
  env.pot <-setDT(grid_sf)[grid.id%in%c(IST$grid.id,eg$grid.id),c("grid.id","BV04","subreg",..VAR)]
  env.pot <- na.omit(setDT(env.pot)[!duplicated(grid.id),]) 
  
  # nrow(EG) - nrow(env.pot)  # 355 ha perdus sur la fronti?re
  # st_geometry(EG) <- "geometry"
  # mapview::mapview(EG[!EG$grid.id%in%env.pot$grid.id,],legend=F)
  
  ## Projection dans l'espace-guilde (to lower computational time)
  env.pot <- merge(env.pot,setDT(IST)[,c("grid.id","BIOIDX_TXG")],by="grid.id",all.x=T)
  env.pot$ist <- env.pot$BIOIDX_TXG
  env.pot$BIOIDX_TXG <- NULL
  
  env.pot <- merge(env.pot,setDT(OBS)[,c("grid.id","Qobs")],by="grid.id",all.x=T)
  env.pot$cal <- env.pot$Qobs
  env.pot$Qobs <- NULL
  table(duplicated(env.pot$grid.id))
  env.pot <- env.pot[!duplicated(env.pot$grid.id),]
  save(env.pot,file=paste0(dir,"/data/",GUILD,"_envpot.RData"))
  
  rm(list=setdiff(ls(), c("env.pot","OBS","VAR","GUILD","guild.info","dir","DATE")))
  save.image(paste0(dir,"/data/",GUILD,"_enviro4model_",DATE,".RData"))


  ################################################
  # Part 2: Modelling
  ################################################
  # dir <- paste0(LOCATION,list.files(LOCATION)[grep(paste0("G",GUILD,"-"),list.files(LOCATION))])
  # FF <- list.files(paste(dir,"data",sep="/"))[grep("enviro4model",as.vector(list.files(paste(dir,"data",sep="/"))))]
  # IF <- file.info(paste(dir,"data", FF,sep="/"))
  # load(paste(dir,"data", FF[which.max(IF$mtime)],sep="/"))
  rstudioapi::restartSession(command='source("code/Modelling_part.R")')

  
  ################################################
  # Comparaison par cluster de bassins versants
  ################################################
  source("C://Dossier_Ervan/R/functions_ER.R",echo=F)
  library(sf)
  library("tidyverse")
  library(data.table)
  DATE <- format(Sys.time(), '%d-%m-%y')
  GUILD <- 1
  guild <- readxl::read_xlsx("C://Dossier_Ervan/Guildes&co/Guildes_revues_ER.xlsx",sheet=1,col_names = T)
  guild.name <- guild[which(guild$No==GUILD),]$name_guild_de
  
  
  # source("code/comp_BV.R")
  BV<- LOAD("clusterBV","data")
  
  ## Associe les clusters aux PRED
  PRED <- LOAD("quality_predicted","data")
  PRED$subreg <- NULL
  PRED <- merge(setDT(PRED),setDT(BV)[,c("BV04","CLUST","subreg","area")],by="BV04",all.x=TRUE)
  
  OBS <- FILTER.OBS(GUILD)
  
  ##### BENCHMARKING
  guild <- readxl::read_xlsx("C://Dossier_Ervan/Guildes&co/Guildes_revues_ER.xlsx",sheet=1,col_names = TRUE)
  guild.name <- guild[which(guild$No==GUILD),]$name_guild_de
  th.pred <-  guild[which(guild$No==GUILD),]$th.pred
  th.bench <- guild[which(guild$No==GUILD),]$th.bench
  th.qprop <- guild[which(guild$No==GUILD),]$th.qprop
  defrag <- guild[which(guild$No==GUILD),]$defrag
  bench <- guild[which(guild$No==GUILD),]$bench
  
  RES <- BENCH(PRED,OBS,th.pred,th.bench,th.qprop,defrag,sd.min=TRUE)
  PRED <- RES[[2]]
  res <- RES[[1]]
  save(PRED,file=paste0(getwd(),"/data/",GUILD,"_qual4leaf_",format(Sys.time(), '%d-%m-%y'),".Rdata"))
  
  ### Reassign geometry
  setDT(res)[is.na(CLUST),c("Nsp_BV","bench_Nsp","Qobs", "Qpred","EB_sp_defrag", "EB_sp_defrag_weighted"):=0]
  st_geometry(res) <- "geometry"
  save(res,file=paste0("data/",GUILD,"_ha2add_",format(Sys.time(), '%d-%m-%y'),".RData"))
  
  
  ### Generate factsheet
  res <- LOAD("ha2add","data")
  res <- col.bin(res,bench,min.size=5)
  res <- res[[1]]
  
  ## Specie contribution
  # RES <- unique(data.table::data.table("group"=OBS$group,"species"=OBS$name,"bioregion"=OBS$bioregrion,"weight"=OBS$w,prop_area_ha=as.numeric(NA),prop_area_percent=as.numeric(NA)))
  # NN <- data.table::setDT(OBS)[,.(Qobs=ifelse(sum(w[!duplicated(taxonid)])>=1,1,0)),by=.(grid.id)]
  # system.time(for(i in 1:nrow(RES)){
  #   NN2 <-  data.table::setDT(OBS)[name!=RES[i,]$species,.(Qobs=ifelse(sum(w[!duplicated(taxonid)])>=1,1,0)),by=.(grid.id)]
  #   RES[i,]$prop_area_ha <- nrow(NN[Qobs==1]) - nrow(NN2[Qobs==1])
  #   RES[i,]$prop_area_percent <- (nrow(NN[Qobs==1]) - nrow(NN2[Qobs==1]))*100/nrow(NN[Qobs==1])
  # }) # 70 sec
  # names(RES) <- c("group", "speciesCODE", "weight", "prop_area_ha", "Contribution")
  # write.csv(RES,file=paste0("C://Dossier_Ervan/R/INFOFAUNA/guilde_",GUILD,"_spContribution.csv"))
  
  rmarkdown::render(input = paste0(getwd(),"/code/factsheet.Rmd"), 
                    output_format = "html_document",
                    output_file = paste0(GUILD,"_factsheet3.html"),
                    output_dir = paste0(getwd(),"/report"))
  
  GUILD <- 1
  
  # Create leaflets for plausibilisation
  library(leafgl)
  library(leaflet)
  library(leafem)
  library(leafpop)
  library(data.table)
  library(raster)
  library(tidyverse)
  library(sf)
  source("C://Dossier_Ervan/R/functions_ER.R")
  guild <- readxl::read_xlsx("C://Dossier_Ervan/Guildes&co/Guildes_revues_ER.xlsx",sheet=1,col_names = TRUE)
  guild.name <- guild[which(guild$No==GUILD),]$name_guild_de
  bench <- guild[which(guild$No==GUILD),]$bench
  
  res <- LOAD("ha2add","data")
  st_geometry(res) <- "geometry"
  
  res2 <- st_transform(res,2056)  ## goes back to Swiss CH1903+ / LV95
  res2 <- col.bin(res2,bench,min.size=5)
  st_write(res2[[1]][,c("BV.id","CLUST","sp.in.BV","sp.in.bench","observed_qual","potential_qual","Erganzungsbedarf")], paste0("shp/",GUILD,"_EB_2056.shp"),delete_layer = T)
  
  res2 <- st_transform(res,4326)
  V1 <- col.bin(res2,bench,min.size=5)
  
  summary(V1[[1]]$potential_qual-V1[[1]]$Erganzungsbedarf)
  
  m = leaflet() %>%
    addProviderTiles(provider = providers$CartoDB.Positron,group="Positron",layerId="Positron") %>%
    addProviderTiles(provider = providers$Esri.WorldImagery,group="Esri",layerId="Esri") %>%
    addPolygons(data = V1[[1]], group = "Erganzungsbedarf",weight=1,color="grey",opacity=0.8,fillColor = ~col,fillOpacity = 0.8,highlightOptions = highlightOptions(color = "black", weight = 2),popup = popupTable(V1[[1]][,c("BV.id","sp.in.BV","sp.in.bench","observed_qual","potential_qual","Erganzungsbedarf")],feature.id = FALSE,row.numbers = FALSE)) %>%
    addLegend(colors=V1[[2]]$col,labels=V1[[2]]$lab,group = "Erganzungsbedarf", position = "topright",opacity=0.8,title=paste0("Erg?nzungsbedarf [ha] (max:",max(V1[[1]]$Erganzungsbedarf,na.rm=T),")")) %>% 
    addPolygons(data = V1[[1]], group = "BV",weight=1,color="grey",opacity=0.8,fillColor = "transparent",fillOpacity = 0,highlightOptions = highlightOptions(color = "black", weight = 2),popup = popupTable(V1[[1]][,c("BV.id","sp.in.BV","sp.in.bench","observed_qual","potential_qual","Erganzungsbedarf")]))
  
  
  m$dependencies = c(m$dependencies,
                     mapview:::popupLayoutDependencies())
  
  ## ADD IST
  PRED <- LOAD("qual4leaf","data")
  st_geometry(PRED) <- "geometry"
  IST <- st_transform(PRED[!is.na(PRED$ist),],3857)
  grid2 <- raster::raster("D://SIG/data/grid100/grid100_3857.tif")
  I_ra <- fasterize::fasterize(IST,grid2,field="ist")
  col1 <- colorNumeric(palette = "viridis",IST$ist,na.color="#00000000")
  raster::writeRaster(I_ra,paste0("shp/",GUILD,"_IST.tif"),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)
  
  m1 <- m %>% 
    addRasterImage(I_ra,colors=col1,method="ngb",group = "Observed qual") %>% 
    addLegend(pal =  col1, values=IST$ist,group = "Observed qual", position = "topright",opacity=0.8,title="Observed quality")
  
  
  ## ADD SOLL
  P <- setDT(PRED)[quality=="pred",]
  st_geometry(P) <- "geometry"
  P <- st_transform(P,3857)
  grid2 <- raster::raster("D://SIG/data/grid100/grid100_3857.tif")
  P_ra <- fasterize::fasterize(P,grid2,field="Qp")  ## replace by "consensus" when priorisation is done
  raster::writeRaster(P_ra,paste0("shp/",GUILD,"_SOLL2.tif"),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)
  
  # col3 = colourvalues::colour_values_rgb(P$prio, palette = "inferno",include_alpha = FALSE)
  col3 <- colorNumeric(palette = "inferno",P$Qp,na.color="#00000000",reverse=T)
  
  m2 <- m1 %>% 
    addRasterImage(P_ra,colors=col3,method="ngb",group = "Potential qual",maxBytes=4454211) %>% 
    addLegend(color =  col3(1),label="" ,group= "Potential qual", position = "topright",opacity=0.8,title="Potential quality")
  
  
  ### Finalize leaflet
  m3 <- m2 %>% 
    addLayersControl(baseGroups= c("Positron","Esri"),overlayGroups = c("Erganzungsbedarf","BV","Observed qual","Potential qual"),position="topleft") %>% 
    hideGroup(c("Observed qual","Potential qual"))
  m3
  mapview::mapshot(m3, url =paste0("report/",GUILD,"_plausi.html"))
  
  
  ##########################################################
  # Create final leaflets
  rm(list=setdiff(ls(),"GUILD"))
  .rs.restartR()
  
  library(leafgl)
  library(leaflet)
  library(leafem)
  library(leafpop)
  library(data.table)
  library(raster)
  library(tidyverse)
  library(sf)
  source("C://Dossier_Ervan/R/functions_ER.R")
  guild <- readxl::read_xlsx("C://Dossier_Ervan/Guildes&co/Guildes_revues_ER.xlsx",sheet=1,col_names = TRUE)
  guild.name <- guild[which(guild$No==GUILD),]$name_guild_de
  bench <- guild[which(guild$No==GUILD),]$bench
  
  res <- LOAD("ha2add","data")
  st_geometry(res) <- "geometry"
  
  res2 <- st_transform(res,2056)  ## goes back to Swiss CH1903+ / LV95
  res2 <- col.bin(res2,bench,min.size=5)
  st_write(res2[[1]][,c("BV.id","CLUST","sp.in.BV","sp.in.bench","observed_qual","potential_qual","Erganzungsbedarf")], paste0("shp/",GUILD,"_EB_2056.shp"),delete_layer = T)
  
  res2 <- st_transform(res,4326)
  V1 <- col.bin(res2,bench,min.size=5)
  
  summary(V1[[1]]$potential_qual-V1[[1]]$Erganzungsbedarf)
  
  m = leaflet() %>%
    addProviderTiles(provider = providers$CartoDB.Positron,group="Positron",layerId="Positron") %>%
    addProviderTiles(provider = providers$Esri.WorldImagery,group="Esri",layerId="Esri") %>%
    addPolygons(data = V1[[1]], group = "EB",weight=1,color="grey",opacity=0.8,fillColor = ~col,fillOpacity = 0.8,highlightOptions = highlightOptions(color = "black", weight = 2),popup = popupTable(V1[[1]][,c("BV.id","sp.in.BV","sp.in.bench","observed_qual","potential_qual","Erganzungsbedarf")],feature.id = FALSE,row.numbers = FALSE)) %>%
    addLegend(colors=V1[[2]]$col,labels=V1[[2]]$lab,group = "EB", position = "topright",opacity=0.8,title=paste0("Erg?nzungsbedarf [ha] (max:",max(V1[[1]]$Erganzungsbedarf,na.rm=T),")"))
  
  
  m$dependencies = c(m$dependencies,
                     mapview:::popupLayoutDependencies())
  
  ## ADD IST
  PRED <- LOAD("qual4leaf","data")
  st_geometry(PRED) <- "geometry"
  IST <- st_transform(PRED[!is.na(PRED$ist),],3857)
  grid2 <- raster::raster("D://SIG/data/grid100/grid100_3857.tif")
  I_ra <- fasterize::fasterize(IST,grid2,field="ist")
  col1 <- colorNumeric(palette = "viridis",IST$ist,na.color="#00000000")
  raster::writeRaster(I_ra,paste0("shp/",GUILD,"_IST.tif"),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)
  
  m1 <- m %>% 
    addRasterImage(I_ra,colors=col1,method="ngb",group = "Observed qual") %>% 
    addLegend(pal =  col1, values=IST$ist,group = "Observed qual", position = "topright",opacity=0.8,title="Observed quality")
  
  
  ## ADD SOLL
  P <- setDT(PRED)[quality=="pred",]
  # ## Priorisation du SOLL (cf code Priorisation)
  PRIO <- LOAD("prio","data")
  P <- merge(setDT(P),setDT(PRIO)[,c("CNHA","consensus")],by="CNHA",all.x=T)
  st_geometry(P) <- "geometry"
  P <- st_transform(P,3857)
  grid2 <- raster::raster("D://SIG/data/grid100/grid100_3857.tif")
  P_ra <- fasterize::fasterize(P,grid2,field="Qp")  ## replace by "consensus" when priorisation is done
  raster::writeRaster(P_ra,paste0("shp/",GUILD,"_SOLL2.tif"),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)
  
  # col3 = colourvalues::colour_values_rgb(P$prio, palette = "inferno",include_alpha = FALSE)
  col3 <- colorNumeric(palette = "inferno",P$Qp,na.color="#00000000",reverse=T)
  
  m2 <- m1 %>% 
    addRasterImage(P_ra,colors=col3,method="ngb",group = "Potential qual",maxBytes=4454211) %>% 
    addLegend(pal =  col3,values=P$Qp,group = "Potential qual", position = "topright",opacity=0.8,title="Potential quality")
  
  ### Add polygon INFOFAUNA
  dir <- "C://Dossier_Ervan/R/INFOFAUNA"
  FF <- list.files(dir)[grep(paste0("_",GUILD,"_cHull"),list.files(dir))]
  pol <-st_read(paste(dir, FF,sep="/"))
  pol <- st_transform(pol,4326)
  pol$PRIO <- ifelse(pol$proxy_IFed==0,"Regional","National")
  col4 = colourvalues::colour_values_rgb(pol$PRIO,palette="rdbu",include_alpha = FALSE)
  col5 <- colorFactor(palette = "RdBu",pol$PRIO)
  
  m3 <- m2  %>%
    addGlPolygons(data = pol, group = "IST prio",popup="PRIO",fillColor=col4,fillOpacity=0.8) %>%
    addLegend(pal =  col5, values=pol$PRIO,group = "IST prio", position = "topright",opacity=2,title="IST priority")
  
  ### Finalize leaflet
  m3 <- m2 %>% 
    addLayersControl(baseGroups= c("Positron","Esri"),overlayGroups = c("EB","Observed qual","Potential qual","IST prio"),position="topleft") %>% 
    hideGroup(c("Observed qual","Potential qual","IST prio"))
  m3
  mapview::mapshot(m3, url =paste0("report/",GUILD,"_final3.html"))
  
  
  
  ### Write raster for BAFU
  ##########################
  res <- LOAD("ha2add","data")
  res <- col.bin(res,bench,min.size=5)
  res2 <- res[[1]]
  st_geometry(res2) <- "geometry"
  res2 <- st_transform(res2,2056)  ## goes back to Swiss CH1903+ / LV95
  st_write(res2[,c("BV.id","CLUST","sp.in.BV","sp.in.bench","observed_qual","potential_qual","Erganzungsbedarf")], paste0("shp/",GUILD,"_EB.shp"),delete_layer = T)

