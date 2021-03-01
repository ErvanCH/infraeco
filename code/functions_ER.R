# LOAD <- function(fileName){
#   #loads an RData file, and returns it
#   load(fileName)
#   get(ls()[ls() != "fileName"])
# }

# Fonction d?sagr?gation IST (Blaise Petitpierre 10/07/2020)
IST.dsg<-function(IST,min.dist=300,IDX=NULL){
  require(sp)
  require(data.table)
  
  my.obse<-setDT(IST)
  if (length(IDX)>0){
    my.obse$idx <-  round(my.obse[,get(IDX)])
  }  else {
  my.obse$idx <- round(my.obse$BIOIDX_TXG)
  }
  coordinates(my.obse)<-my.obse[,c("CX","CY")]
  to.keep<-c()
  k<-1
  run <- sort(unique(my.obse@data$idx),decreasing = T)
  
  for (i in 1:length(run)){
    IST2run <- sort(unique(my.obse@data$idx),decreasing = T)
    if(length(my.obse@data$CNHA)==0) {break}
    ref.i<-my.obse[which(my.obse@data$idx==max(IST2run)),]
    ref.e<-remove.duplicates(ref.i,zero=min.dist)
    to.keep<-c(to.keep,ref.e$CNHA)
    
    my.obsi<-my.obse[-which(my.obse@data$idx==max(IST2run)),]
    to.rm<-unique(zerodist2(ref.e,my.obsi,zero=min.dist)[,2])
    if(length(to.rm)>0){my.obse<-my.obsi[-to.rm,]
    cat(paste0(length(IST2run),"-"))}
    else{my.obse<-my.obsi}
  }
  
  setDT(IST)[CNHA%in%to.keep]
}

### Load data
LOAD <- function(file,in.folder=NULL,dir=NULL){   #load most recent RData file, and returns it by setting its name
  if (length(dir)==0) {  dir = getwd()  }
  if (length(in.folder)>0) {
    FF <- list.files(paste(dir,in.folder,sep="/"))[grep(file,as.vector(list.files(paste(dir,in.folder,sep="/"))))]
    IF <- file.info(paste(dir,in.folder, FF,sep="/"))
    fileName <- paste(dir,in.folder, FF[which.max(IF$mtime)],sep="/")
  } else {
    fileName = file
  }
  load(fileName)
  get(setdiff(ls(),c("IF","FF","dir","fileName","in.folder","file")))
  }
  
  

LABEL.VAR <- function(X) {
  lab <- read.csv2("D://SIG/data/raster_enviro/label_var_enviro.csv")
  return(lab$label[match(X,lab$var)])
}

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
             'growing_deg', #growing degree days
             'rad', # radiations
             'slope',  # slope
             'topo', # topography (curvature)
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
             'mask.tif')
##  Model diversity
model.div <- function(div.obs,env.pot,VAR) {
  LR <- as.numeric(scale(nrow(div.obs),center=F))/10  ## learning rates by group
  # # Debug
  # div.obs <- div_list[[4]]
  # env.pot <- buf2
  
  env.obs <- na.omit(setDT(div.obs)[W==1,c("grid.id","N","W",..VAR)])
  if (nrow(env.obs)<100) {
    cat("!! Not enough samples for",unique(div.obs$group),"!!")
    PRED <- NA
  } else {
    env.pot <- dplyr::left_join(na.omit(env.pot),env.obs[,1:3],by="grid.id",all.x=TRUE)
    env.pot$obs <- 0
    env.pot <- within(env.pot,obs[grid.id%in%env.obs$grid.id] <- 1)
    fmla <- as.formula(paste("N ~ ", paste(VAR,collapse = "+ ")))
    PRED <- env.pot[,c("grid.id","BV04","N")]
    
    ## Random forest unweighted
    rf1 <- ranger::ranger(fmla, data=env.obs)
    pr <- predict(rf1,env.pot,type="response")
    pr$predictions[env.pot$obs==0] <- pr$predictions[env.pot$obs==0]/max(pr$predictions[env.pot$obs==0],na.rm=TRUE)
    PRED$rf <- pr$predictions
    rm(list=c("rf1","pr"))
    
    # ## Boosted regression tree (BRT)
    ## shell.exec("https://cran.r-project.org/web/packages/dismo/vignettes/brt.pdf")
    testError <- tryCatch(gbm <- dismo::gbm.step(env.obs,gbm.x=VAR,gbm.y="N",family="poisson",learning.rate = LR, bag.fraction = 0.5,max.trees=1000,silent=F), error=function(e) e)
    if(inherits(testError, "error")|is.null(gbm)) {
      gbm <- dismo::gbm.step(env.obs,gbm.x=VAR,gbm.y="N",family="poisson",tree.complexity = 5,learning.rate = LR/100,  site.weights= sqrt(env.obs$W+1),bag.fraction = 0.8,max.trees=1000,silent=F)
    }
    prw <- predict(gbm,env.pot,n.trees=500,,type="response")
    prw[env.pot$obs==0] <- prw[env.pot$obs==0]/max(prw[env.pot$obs==0],na.rm=TRUE)
    PRED$brt <- prw
    rm(list=c("gbm","prw"))
    
    # GAM model
    fmlaG <- as.formula(paste("N ~ ", paste(VAR,collapse = "+ ")))
    BAM <- mgcv::bam(fmlaG,data=env.obs,family=poisson)
    prw <- as.numeric(predict(BAM,env.pot,type="response"))
    prw[env.pot$obs==0] <- prw[env.pot$obs==0]/max(prw[env.pot$obs==0],na.rm=TRUE)
    PRED$gam <- prw
    rm(list=c("BAM","prw"))
    
    ## STACK prediction by group
    PRED
  }
}

RGB_IE <- function(x, y){
  R <- x
  G <- y
  B <- (1-x)/2
  A <- 1- 0.2*exp(-(x^2+y^2)/0.2)
  
  rgb(R, G, B,A)
}

# Compare bassin versants:
rasterizeBV <- function(PCAl,glob,axes,z.th=0.01) {
  xmin <- min(glob[,axes[1]])
  xmax <- max(glob[,axes[1]])
  ymin <- min(glob[,axes[2]])
  ymax <- max(glob[,axes[2]])
  sp<-as.matrix(PCAl)[,axes]  # BVs of a group (=bioreg)
  spr <- data.frame(cbind((sp[,1] - xmin)/abs(xmax - 
                                                xmin), (sp[,2] - ymin)/abs(ymax - ymin)))
  sp.dens <- kernelUD(SpatialPoints(spr[, 1:2]), h = "href", 
                      grid = mask, kern = "bivnorm")
  sp.dens <- raster(xmn = xmin, xmx = xmax, ymn = ymin, 
                    ymx = ymax, matrix(sp.dens$ud, nrow = R))
  z <- sp.dens * nrow(sp)/cellStats(sp.dens, "sum")
  z.uncor <- z/cellStats(z, "max")
  z.thi<-quantile(raster::extract(z.uncor,sp[,1:2]),z.th,na.rm=TRUE)
  z.uncor[z.uncor<z.thi]<-0
  return(z.uncor)
}

env.overlap<-function(rasterselect){
  
  combi <- expand.grid(names(rasterselect),names(rasterselect))
  combi$D <- combi$I <- combi$COR <- NA
  
  for (i in 1:nrow(combi)) {
    z1<-rasterselect[[combi[i,1]]]
    z2<-rasterselect[[combi[i,2]]]
    
    p1 <- as.matrix(z1)/sum(as.matrix(z1))
    p2 <- as.matrix(z2)/sum(as.matrix(z2))
    
    combi[i,3] = 1 - (0.5 * (sum(abs(p1 - p2)))) # Shoener's D
    
    combi[i,4] = 1 - ((sqrt(sum((sqrt(p1) - sqrt(p2))^2)))^2)/2 # Hellinger's I
    
    combi[i,5] = cor(values(z1),values(z2)) # Pearson's correlation
  }
  return(combi)
}

pca.rasterize<-function(X,VAR,axes,mask,R,z.th=0.01,save_raster=F){
  REG <- unique(X$bioregion)
  PCA <- ade4::dudi.pca(setDT(X)[,..VAR],scannf=F,nf=10)  # group = bioreg
  EIGEN <- PCA$eig/sum(PCA$eig)
  
  rasterlist <- list()
  
  for (i in 1:nrow(AXES)) { 
    PCAlist <- split(PCA$li[1:3],X$BV04)
    BB <- sapply(PCAlist,nrow)
    PCAlist[which(names(PCAlist) %in% names(BB)[which(as.numeric(BB)<6)])] <- NULL  # remove BV with less than 6 ha
    # du to the fact that kernelUD() needs at least 5 relocations to fit a home range
    # if (all(names(PCAlist) %in% names(BB))) { ## if not enough ha, drop this step
    #   COMBI <- NA
    #   }  else {
    STACK <- stack(sapply(PCAlist,function(x) rasterizeBV(x,glob=as.matrix(PCA$li),axes=as.numeric(AXES[i,]))))
    if (save_raster==TRUE) {
      if (!any(grepl("rasterBV",list.files(getwd())))){
        dir.create(paste0(getwd(), "/rasterBV"))
      }
      writeRaster(STACK, paste0(getwd(),'/rasterBV/rasterBV_BR_',REG,'_axes',paste(AXES[i,],collapse = ''),'.tif'),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE) }
    names(STACK) <- names(PCAlist)
    rasterlist[[i]] <- STACK
      }
 
  
  OVER <- lapply(rasterlist,env.overlap)  # compute ENVIRONMENTAL OVERLAP INDICES
  COMBI <- OVER[[1]][,1:3]
  WEIG <- c(EIGEN[1]+EIGEN[2],EIGEN[1]+EIGEN[3],EIGEN[2]+EIGEN[3])
  COR <- as.data.frame(do.call(cbind,lapply(OVER,function(x) x[,3:5])))
  COR2 <- as.data.frame(t(t(COR)*WEIG)) ## Weigh indices by eigen values
  COMBI$cor_m <- rowMeans(COR2)  ## mean correlation among indices
  
  names(COMBI) <- c("BV1","BV2","COR","cor_m")
  COMBI <- COMBI[COMBI$COR!=1,]  # remove lines of auto-comparison 
  COMBI$COR <- NULL
  
  COMBI 
  # } # end of loop
}

CLUST <- function(x,CENTERS) {
  set.seed(10)
  B <- setDT(x)[,c("BV1","BV2","cor_m")]
  BB <- unique(B$BV1)
  MT <- matrix(0, length(BB),length(BB))
  colnames(MT) <- rownames(MT) <- BB
  for (i in 1:length(BB)) {
    b <- match(unique(B$BV1)[i],BB)
    IDX <- which(unique(B$BV1)[i]==B$BV1)
    MT[b,match(B$BV2[IDX],BB)] <- B$cor_m[which(unique(B$BV1)[i]==B$BV1)]
  }
  MT[lower.tri(MT)] <- MT[upper.tri(MT)]
  diag(MT) <- 1
  
  CLUST <- kmeans(MT,centers=length(BB)/CENTERS)
  CLU <- data.frame(BV04=BB,clust=CLUST$cluster)
  CLU
}  


FILTER.OBS <- function(GUILD,EXCLUDE=NULL) {
  require(data.table)
  require(sf)
  # PATH <-"C:/Dossier_Ervan/Guildes&co/Notes/"
  # FILES <- list.files(PATH)
  # library(openxlsx)
  # for (i in 1:length(FILES)) {
  #   temp <- openxlsx::read.xlsx(paste0(PATH,FILES[i]),1)
  #   temp$REGION <- NULL
  #   if (i==1) {
  #     TEMP <- temp
  #   } else {
  #     TEMP <- rbind(TEMP,temp)
  #   }
  # }
  # summary(TEMP)
  ## IMport table with notes
  REF <- data.table::setDT(read.csv2("C:/Dossier_Ervan/Guildes&co/Notes/IE_MATRICE_PRIO_20200715.csv"))
  REF[,"IMS_TAXONIDCH":=TAXON_ID_CH]
  ID <- as.character(GUILD)
  REF2 <- subset(REF,GUILD%in%ID & CIBLE==1)
  REF2$name <-  stringr::word(REF2$TAXON_NAME, 1, 2)
  
  # load obs
  # OBS2 <- setDT(read.table("C:/Dossier_Ervan/Guildes&co/G01-23_20200414.tsv",sep="\t",h=TRUE))
  OBS <- setDT(read.csv("C:/Dossier_Ervan/Guildes&co/Notes/G01-23_20200520.csv",h=TRUE))
  
  # Filter by group if needed
  obs2 <- merge(OBS[,-1],unique(REF2[,c("IMS_TAXONIDCH","GROUP")]),by=c("IMS_TAXONIDCH"),all.x=TRUE)
  obs2 <- obs2[!is.na(GROUP) & !is.na(IMS_SWISSCOORDINATE_X),-c("IMS_PRIORITYCH","CIBLE")] # remove obs that are not in the guild
  table(obs2$GROUP)
  obs2[obs2$IMT_GENUS=="Coleosporium",]
  
  # add hectare ID to pictis data frame
  i <- which(obs2$IMS_SWISSCOORDINATE_Y>=100000);length(i)
  obs2[i,'CNHA'] <- paste(floor(obs2[i,IMS_SWISSCOORDINATE_X]/100),floor(obs2[i,IMS_SWISSCOORDINATE_Y]/100),sep='')
  i <- which(obs2$IMS_SWISSCOORDINATE_Y<100000);length(i)
  obs2[i,'CNHA'] <- paste(floor(obs2[i,IMS_SWISSCOORDINATE_X]/100),'0',floor(obs2[i,IMS_SWISSCOORDINATE_Y]/100),sep='')
  obs2$CNHA <- as.numeric(obs2$CNHA)
  
  if (!is.null(EXCLUDE)) {
    obs2 <- obs2[!GROUP%in%EXCLUDE,]
  }
  
  # Add bioregion to do the merge by region
  # library(sf)
  load("C:/Dossier_Ervan/R/Grid100/grid100_light.Rdata")
  obs2 <- merge(obs2,grid_sf[,c("grid.id","CNHA","subreg")],by="CNHA",all.x=TRUE)
  obs2 <- setDT(obs2)[!is.na(grid.id),]  # remove obs outside boundaries
  obs2 <- sf::st_as_sf(obs2,coords=c("IMS_SWISSCOORDINATE_X","IMS_SWISSCOORDINATE_Y"),crs=21781)
  
 if (nrow(REF2[duplicated(REF2$IMS_TAXONIDCH),])>0) {
    print(cat("Duplicated taxon id in reference list"))
    REF2 <- REF2[!duplicated(REF2$IMS_TAXONIDCH),]
  } 
  obs <- merge(obs2,REF2[,c("YEAR_min","YEAR_max","ALT_min","ALT_max","PREC_min","GRP_WEIGHT","IMS_TAXONIDCH")],by="IMS_TAXONIDCH",all.x=TRUE)
  rm(obs2)
  table(obs$GRP_WEIGHT,obs$GROUP) 
  
  names(obs) <- c("IMS_TAXONIDCH", "CNHA", "ISPCHID", "IMS_COORDMAXDEVIATION", 
                  "IMS_ELEVATION", "IMS_DEPTH", "IMS_YEARCOLLECTED", "IMS_MONTHCOLLECTED", 
                  "IMS_DAYCOLLECTED", "IMT_GENUS", "IMT_SUBGENUS", "IMT_SPECIES", 
                  "IMT_INFRARANK", "IMT_SUBSPECIES", "IMS_ORIGINCH", "GROUP", "grid.id", 
                  "subreg", "YEAR_min", "YEAR_max", "ALT_min", "ALT_max", "PREC_min", 
                  "GRP_WEIGHT","geometry")
 
   names(obs) <- c("taxonid","CNHA", "obsid","coord_dev","alt_coll", "depth_coll", "year_coll", "month_coll", "day_coll", "genus", "subgenus","species", "rank", "subsp", "origin", "group","grid.id","subreg","YEAR_min", "YEAR_max", "ALT_min", "ALT_max", "PREC_min", "w","geometry")
  obs <- setDT(obs)[,c("day_coll","subgenus","rank","subsp","origin"):=NULL]
  
  
    # Sans information sur les groupes, on considère que chaque espèce est à 0.2
  obs <- subset(obs,!is.na(group))
  obs$taxonid <- as.factor(obs$taxonid)
  if(length(obs[is.na(obs$w),]$taxonid)>0){
    obs <- within(obs,w[is.na(w)] <- 0.2)}
  obs[,c("MONTH_min","MONTH_max"):=list(1,12)]
  obs$alt_coll <- as.numeric(obs$alt_coll)
  obs$w <- as.numeric(as.character(obs$w))
    # Filtre par années et précision d'obs
  obs[,"keep" :=as.numeric(1)]
  obs[,"YEAR_out":=ifelse(year_coll<YEAR_min,1,0),by=.SD]
  obs[,"COORD_out":=ifelse(coord_dev>PREC_min,1,0),by=.SD]
  obs[,"MONTH_out":=ifelse(month_coll<MONTH_min|month_coll>MONTH_max,1,0),by=.SD]
  obs[,"ALT_out":=ifelse(alt_coll<ALT_min|alt_coll>ALT_max,1,0),by=.SD]
  obs <- within(obs,keep[YEAR_out==1 | COORD_out==1 | MONTH_out==1 | ALT_out==1]<- 0)
  # obs[obs$CNHA==52581687,]
  obs <- subset(obs,keep==1)
  
  # regoupement des EPHE + PLEC + TRIC 
  obs$group <- as.character(obs$group)
  obs <- within(obs, group[group%in%c("EPHE","PLEC","TRIC")] <- "E-P-T")
  obs[,"name":=paste(genus,species,sep=" ")]
  
  # Promote to sf
  # obs <- sf::st_as_sf(obs[keep==1,c("taxonid", "genus", "species", "ssp", "x", "y","id.centre", "group", "w")],coords=c("x","y"),crs=21781)
  obs}




IMPORT.OBS.GM <- function(GUILD) {
  require(data.table)
  require(sf)
  # PATH <-"C:/Dossier_Ervan/Guildes&co/Notes/"
  # FILES <- list.files(PATH)
  # library(openxlsx)
  # for (i in 1:length(FILES)) {
  #   temp <- openxlsx::read.xlsx(paste0(PATH,FILES[i]),1)
  #   temp$REGION <- NULL
  #   if (i==1) {
  #     TEMP <- temp
  #   } else {
  #     TEMP <- rbind(TEMP,temp)
  #   }
  # }
  # summary(TEMP)

  ## Load observations (Mamm/Chiro)
  REF <- readxl::read_xlsx("C:/Dossier_Ervan/Guildes&co/Guildes_mobiles/Liste_sp_guildes_mobiles_update_2020-04-15.xlsx")
  REF$IMS_TAXONIDCH = as.factor(paste(REF$CENTRE,REF$TAXON_ID,sep=":"))
  REF$ALT_min <- as.numeric(REF$ALT_min)
  REF$ALT_max <- as.numeric(REF$ALT_max)
  setDT(REF)[is.na(ALT_max),"ALT_max":=as.numeric(9999)]
  REF[is.na(ALT_min),"ALT_min":=as.numeric(0)]
  REF$PREC_min <- as.numeric(REF$PREC_min)
  REF[is.na(PREC_min),"PREC_min":=as.numeric(9999)]
  ID <- as.character(GUILD)
  REF2 <- subset(REF,GUILD%in%ID)
  
  # Filtrage des obs 
  # # A. Oiseaux sensibles (pas ds cette guilde)
  AVES <- readxl::read_xlsx("C:/Dossier_Ervan/Guildes&co/Guildes_mobiles/ha-Daten_sensible_Vogelarten_2000-2019.xlsx")
  AVES$IMS_TAXONIDCH = as.factor(paste("vogelwarte",AVES$Artid,sep=":"))
  AVES$IMS_COORDMAXDEVIATION <- NA
  AVES$IMS_MONTHCOLLECTED <- NA
  AVES$REPRO <- NA
  AVES <- merge(AVES,REF,by="IMS_TAXONIDCH",all.x=TRUE)
  
  summary(tapply(AVES$Artid,AVES$`Nom fr`,function(x) length(unique(x))))  # v?rifie que les id soient unique
  
  AVES <- setDT(AVES)[GUILD%in%ID,c("IMS_TAXONIDCH","GENUS", "SP", "NAME","koordx_genau","koordy_genau","IMS_COORDMAXDEVIATION","Hohe","jahr","IMS_MONTHCOLLECTED","REPRO", "GROUP","YEAR_min", "YEAR_max", "MONTH_min","MONTH_max","ALT_min", "ALT_max","PREC_min","GRP_WEIGHT","GUILD")]
  names(AVES)  <- c("taxonid", "genus", "species","name","x","y","coord_dev", "alt_coll","year_coll","month_coll","repro","group", "YEAR_min", "YEAR_max", "MONTH_min", "MONTH_max","ALT_min","ALT_max","PREC_min", "w","guild")
  
  # Transform coordinates
  COO <- st_geometry(st_as_sf(AVES,coords=c("x","y")))
  st_crs(COO) <- 2056
  COO2 <- st_transform(COO,21781)
  summary(st_coordinates(COO2))
  AVES$x <- st_coordinates(COO2)[,1]
  AVES$y <- st_coordinates(COO2)[,2]
  
  # B. Autres groupes
  # obs <- read.table(file ="C://Dossier_Ervan/Guildes&co/Guildes_mobiles/202002240938_CL_Gm.tsv", sep = '\t', header = TRUE)
  # ob <- read.table(file ="C://Dossier_Ervan/Guildes&co/Guildes_mobiles/pelophylax-esculentus-aggr.tsv", sep = '\t', header = TRUE)
  # names(ob) <- names(obs)
  # obs <- rbind(obs,ob)
  # save(obs,file="C://Dossier_Ervan/Guildes&co/Guildes_mobiles/obs_GM_06.20.Rdata")
  # ### UPDATE OBS OF CHIRO (24/07)
  # REF <- readxl::read_xlsx("C:/Dossier_Ervan/Guildes&co/Guildes_mobiles/Liste_sp_guildes_mobiles_update_2020-04-15.xlsx")
  # REF$IMS_TAXONIDCH = as.factor(paste(REF$CENTRE,REF$TAXON_ID,sep=":"))
  # REF$ALT_min <- as.numeric(REF$ALT_min)
  # REF$ALT_max <- as.numeric(REF$ALT_max)
  # setDT(REF)[is.na(ALT_max),"ALT_max":=as.numeric(9999)]
  # REF[is.na(ALT_min),"ALT_min":=as.numeric(0)]
  # REF$PREC_min <- as.numeric(REF$PREC_min)
  # REF[is.na(PREC_min),"PREC_min":=as.numeric(9999)]
  # # load("C://Dossier_Ervan/Guildes&co/Guildes_mobiles/obs_GM_06.20.Rdata")
  # obs2 <- merge(obs,REF,by="IMS_TAXONIDCH",all.x=TRUE)
  # 
  # table(obs2$GROUP,obs2$GUILD)
  # 
  # CH <- read.table(file ="C://Dossier_Ervan/Guildes&co/Guildes_mobiles/IE_Gm_CHIRO_GUILDES_24?26_20200723.tsv", sep = '\t', header = TRUE)
  # names(CH) <- c("?..IMS_INS_SHORTNAMEWORLD", "CANTON", "IMS_DAT_DBIDENT", "IMS_TAXONIDCH",
  #                "IMS_ID", "IMS_SWISSCOORDINATE_X", "IMS_SWISSCOORDINATE_Y",
  #                "PR", "IMS_COORDMAXDEVIATION", "IMS_ELEVATION", "IMS_DEPTH",
  #                "IMS_YEARCOLLECTED", "IMS_MONTHCOLLECTED", "IMS_DAYCOLLECTED",
  #                "IMT_GENUS", "IMT_SUBGENUS", "IMT_SPECIES", "IMT_INFRARANK",
  #                "IMT_SUBSPECIES", "IMS_PRIORITYCH", "IMS_ORIGINCH")
  # CH$CANTON <- NULL
  # CH$PR <- NULL
  # CH$REPRO <- NA
  # CH <- merge(CH,REF,by="IMS_TAXONIDCH",all.x=TRUE)
  # obs3 <- setDT(obs2)[-which(obs2$GUILD%in%c("24-25","26") & obs2$GROUP=="CHIR"),]
  # CH <- CH[,names(obs3)]
  # names(CH) <- names(obs3)
  # obs <- rbind(obs3,CH[CH$GUILD%in%c("24-25","26"),])
  # save(obs,file="C://Dossier_Ervan/Guildes&co/Guildes_mobiles/obs_GM_06.20.Rdata")
  
  load("C://Dossier_Ervan/Guildes&co/Guildes_mobiles/obs_GM_06.20.Rdata")
  obs2 <- merge(obs,REF2,by="IMS_TAXONIDCH",all.x=TRUE)
   
  if (GUILD==23){
    obs2 <- subset(obs2,!GROUP%in%c("CHIR","AVES"))
    chiro <- read.csv2(file = "C:/Dossier_Ervan/Guildes&co/Guildes_mobiles/Extrait_gites_chiros_juv__ARE_TLMbz_20200518.csv",  header = TRUE)
    obs3 <- merge(obs,REF2[,c("IMS_TAXONIDCH","GROUP","GUILD")],by="IMS_TAXONIDCH",all.x=TRUE)
    obs3 <- subset(obs3,GROUP!="CHIR" & GUILD%in%ID)
    temp <- data.frame(matrix(NA,nrow(chiro),ncol(obs3)))
    names(temp) <- names(obs3) 
    temp$IMS_YEARCOLLECTED <- chiro$A
    temp$IMS_TAXONIDCH <- chiro$TAXON_ID
    temp$IMS_SWISSCOORDINATE_X <- chiro$CX
    temp$IMS_SWISSCOORDINATE_Y <- chiro$CY
    temp$IMS_ELEVATION <- chiro$ALT
    temp$GROUP <- NULL
    temp$GUILD <- NULL
    temp2 <- merge(temp,REF2,by="IMS_TAXONIDCH",all.x=TRUE)
    
    obs2 <- rbind(obs2,temp2)
    obs2 <- setDT(obs2)[GUILD%in%ID,c("IMS_TAXONIDCH","GENUS","SP", "NAME","IMS_SWISSCOORDINATE_X", "IMS_SWISSCOORDINATE_Y", "IMS_COORDMAXDEVIATION", "IMS_ELEVATION","IMS_YEARCOLLECTED", "IMS_MONTHCOLLECTED","REPRO","GROUP","YEAR_min", "YEAR_max", "MONTH_min", "MONTH_max", "ALT_min", "ALT_max", "PREC_min",  "GRP_WEIGHT","GUILD") ]
    
    AVE <- read.table("C:/Dossier_Ervan/Guildes&co/Guildes_mobiles/PICTIS_SOI_Gm_202009181.tsv", sep = '\t', header = TRUE) 
    AVE <- merge(AVE,REF2,by="IMS_TAXONIDCH",all.x=TRUE)
    AVE <- setDT(AVE)[!is.na(GUILD) & IMS_COORDMAXDEVIATION==71 & REPRO==3,]
    AVES <-AVE[,c("IMS_TAXONIDCH","GENUS", "SP", "NAME","IMS_SWISSCOORDINATE_X", "IMS_SWISSCOORDINATE_Y","IMS_COORDMAXDEVIATION","IMS_ELEVATION","IMS_YEARCOLLECTED","IMS_MONTHCOLLECTED","REPRO", "GROUP","YEAR_min", "YEAR_max", "MONTH_min","MONTH_max","ALT_min", "ALT_max","PREC_min","GRP_WEIGHT","GUILD")]
       
    obs <- rbind(obs2,AVES)
    names(obs)  <- c("taxonid", "genus", "species","name","x","y","coord_dev", "alt_coll","year_coll","month_coll","repro","group", "YEAR_min", "YEAR_max", "MONTH_min", "MONTH_max","ALT_min","ALT_max","PREC_min", "w","guild")
    
    } else {
  obs2 <- setDT(obs2)[GUILD%in%ID,c("IMS_TAXONIDCH","GENUS","SP", "NAME","IMS_SWISSCOORDINATE_X", "IMS_SWISSCOORDINATE_Y", "IMS_COORDMAXDEVIATION", "IMS_ELEVATION","IMS_YEARCOLLECTED", "IMS_MONTHCOLLECTED","REPRO","GROUP","YEAR_min", "YEAR_max", "MONTH_min", "MONTH_max", "ALT_min", "ALT_max", "PREC_min",  "GRP_WEIGHT","GUILD") ]
  
  summary(tapply(obs2$IMS_TAXONIDCH,obs2$NAME,function(x) length(unique(x))))  # v?rifie que les id soient unique
  
  names(obs2)  <- c("taxonid", "genus", "species","name","x","y","coord_dev", "alt_coll","year_coll","month_coll","repro","group", "YEAR_min", "YEAR_max", "MONTH_min", "MONTH_max","ALT_min","ALT_max","PREC_min", "w","guild")
  
  obs <- rbind(obs2,AVES)
    }
  
  
  # Filtre par annees, mois et precision d'obs
  library(data.table)
  obs$keep <- 1
  obs[,"YEAR_out":=ifelse(year_coll<YEAR_min,1,0),by=.SD]
  obs[,"COORD_out":=ifelse(coord_dev>PREC_min,1,0),by=.SD]
  obs[,"MONTH_out":=ifelse(month_coll<MONTH_min|month_coll>MONTH_max,1,0),by=.SD]
  obs[,"ALT_out":=ifelse(alt_coll<ALT_min,1,0),by=.SD]
  table(obs$MONTH_out)
  obs <- within(obs,keep[YEAR_out==1 | COORD_out==1 | MONTH_out==1 | ALT_out==1]<- 0)
  obs <- within(obs,keep[is.na(x)|is.na(y)]<- 0)
  table(obs$keep)
  obs <- obs[keep==1]  # filter here
  obs[,c("keep","YEAR_out","COORD_out","MONTH_out","ALT_out"):=NULL]
  
  # add hectare ID to pictis data frame
  i <- which(obs$y>=100000);length(i)
  obs[i,'CNHA'] <- paste(floor(obs[i,x]/100),floor(obs[i,y]/100),sep='')
  i <- which(obs$y<100000);length(i)
  obs[i,'CNHA'] <- paste(floor(obs[i,x]/100),'0',floor(obs[i,y]/100),sep='')
  obs$CNHA <- as.numeric(obs$CNHA)
  
  obs
  }


# # # Import data InfoFlora
# FA <- "C://Dossier_Ervan/Guildes&co/Notes/INFOFLORA_matrice_sp_guilds_final.xlsx"
# 
# SHE <- length( readxl::excel_sheets( FA ) )
# for(i in 1:SHE) {
#   if(i==1) {
#   DA <- readxl::read_xlsx(FA,sheet=i)
#   names(DA) <- c("CENTRE", "TAXON_ID", "checklist_id", "TAXON_NAME", "PRIO",
#                  "THREAT", "YEAR_min", "YEAR_max", "ALT_min", "ALT_max", "PREC_min",
#                  "Guild", "GRP_WEIGHT")
#   } else {
#   da <- readxl::read_xlsx(FA,sheet=i)
#   names(da) <- c("CENTRE", "TAXON_ID", "checklist_id", "TAXON_NAME", "PRIO",
#                  "THREAT", "YEAR_min", "YEAR_max", "ALT_min", "ALT_max", "PREC_min",
#                  "Guild", "GRP_WEIGHT")
#   DA <- rbind(DA,da)
#   }
# }
# 
# 
# data.table::setDT(DA)[,c("IMS_TAXONIDCH","GROUP","RESP","REGION","GUILD","YEAR_max","ALT_max"):=list(paste(CENTRE,TAXON_ID,sep=":"),"TRAC",NA,"CH",substr(Guild,2,4),NA,NA)]
# DA <- DA[,c("IMS_TAXONIDCH","CENTRE", "GROUP", "TAXON_ID", "TAXON_NAME",
#   "PRIO", "THREAT", "RESP", "YEAR_min", "YEAR_max", "ALT_min",
#   "ALT_max", "PREC_min", "REGION", "GUILD", "GRP_WEIGHT")]
# write.csv2(DA,file="C://Dossier_Ervan/Guildes&co/Notes/INFOFLORA_matrice_compil_06.01.20.csv",row.names=F)


#=================================================================================================================================
# Function for univariate non-parametric t-test, useful for predictor selection
# j: index of the variable
# Nrep : number of replications
# mcmc : number of resampling
# env.pres : data for the presences
# env.bck : data.for the absences
#=================================================================================================================================

var.test<-function(j,Nrep=5,mcmc=1000,env.pres,env.bck){
  p.val<-c()
  for(i in 1:Nrep){
    # x<-pull(env.pres[,..j])
    # y<-pull(env.bck[,..j])
    x<-env.pres[,j]
    y<-env.bck[,j]
    set.seed(i)
    if(length(y)>length(x)){
      yi<-y[sample(1:length(y),size=length(x))]
    }else{
      yi<-y
    }
    p.vali<-permTS(as.vector(x),as.vector(yi),alternative="two.sided", method="exact.mc",control=permControl(nmc=mcmc-1,setSEED=F))$p.value
    p.val<-c(p.val,p.vali)
  }
  p.val=round(c(mean(p.val),sd(p.val)),3)
  names(p.val)<-c("p.val","sd")
  return(p.val)
}


#=================================================================================================================================
# Functions for the model evaluations
#=================================================================================================================================
boycei<-function(interval,obs,fit){
  
  fit.bin<-fit
  obs.bin<-obs
  fit.bin[fit[]>=interval[1]&fit[]<=interval[2]]<-"i";fit.bin[fit.bin!="i"]<-0
  obs.bin[obs[]>=interval[1]&obs[]<=interval[2]]<-"i";obs.bin[obs.bin!="i"]<-0
  
  pi<-length(which(obs.bin=="i"))/length(obs)
  ei<-length(which(fit.bin=="i"))/length(fit.bin)
  fi<-pi/ei
  
  return(fi)
}


ecospat.boyce <- function (fit, obs, nclass = 0, window.w = "default", res = 100,
                           PEplot = TRUE, cor.method="spearman")
{
  if(class(fit)=="RasterLayer"){
    if(class(obs)=="data.frame"){
      obs<-extract(fit, obs)}
    fit<-getValues(fit)
    fit <- fit[!is.na(fit)]
  }
  
  if (window.w == "default") {
    window.w <- (max(fit) - min(fit))/10
  }
  interval <- c(min(fit), max(fit))
  mini <- interval[1]
  maxi <- interval[2]
  if (nclass == 0) {
    vec.mov <- seq(from = mini, to = maxi - window.w, by = (maxi -
                                                              mini - window.w)/res)
    vec.mov[res + 1] <- vec.mov[res + 1] + 1 #Trick to avoid error with closed interval in R
    interval <- cbind(vec.mov, vec.mov + window.w)
  }
  else if (length(nclass) > 1) {
    vec.mov <- c(mini, nclass)
    interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
  }
  else if (nclass > 0 & length(nclass) < 2) {
    vec.mov <- seq(from = mini, to = maxi, by = (maxi - mini)/nclass)
  }
  f <- apply(interval, 1, boycei, obs, fit)
  to.keep<-which(f!="NaN") # index to keep no NaN data
  f<-f[to.keep]
  if (length(f) < 2) {
    b <- NA #at least two points are necessary to draw a correlation
  }
  else {
    r<-c(1:length(f))[f!=c(f[-1],FALSE)] #index to remove successive duplicates
    b <- cor(f[r], vec.mov[to.keep][r], method = cor.method )# calculation of the spearman correlation (i.e. Boyce index) after removing successive duplicated values
  }
  HS <- apply(interval, 1, sum)/2 # mean habitat suitability in the moving window
  HS[length(HS)]<-HS[length(HS)]-1 #Correction of the "trick" to deal with closed interval
  HS<-HS[to.keep] #exlude the NaN
  if (PEplot == TRUE){
    plot(HS, f, xlab = "Habitat suitability",
         ylab = "Predicted/Expected ratio",col="grey",cex=0.75)
    points(HS[r],f[r],pch=19,cex=0.75)
  }
  results <- list(F.ratio = f, correlation = round(b, 3), HS = HS)
  return(results)
}


KappaRepet <-  function(Obs, Fit, TSS=FALSE)
{
  if(sum(Obs)==0) stop("\n The observed data only contains 0")
  tab <- as.data.frame(matrix(0, nrow=101, ncol=2))  ### il faut pr?ciser que: "nrow = 101" sinon le data frame se remplit de "NA" au fur et ??? mesure plut?t que de "0" (et si il y a des NA cela pose probl?me plus loin).
  
  if(length(unique(Fit))==1){
    Misc<-table(as.vector(Fit) >= as.numeric(unique(Fit)), Obs) ### Robin modified here...(avant ?a plantait l???)
    if(TSS!=TRUE) a <- KappaStat(Misc)
    else a <- TSS.Stat(Misc)
    TP <- Misc[4]
    TN <- Misc[1]
    ca0 <- (TN * 100)/sum(Misc[,1])
    ca1 <- (TP * 100)/sum(Misc[,2])
    if(is.na(ca0)) ca0<-0
    if(is.na(ca1)) ca1<-0
    if(TSS!=TRUE) invisible(list(Kappa=0, CutOff=unique(Fit), TP=TP, se=ca1, TN=TN, sp=ca0))
    else invisible(list(TSS=0, CutOff=unique(Fit), TP=TP, se=ca1, TN=TN, sp=ca0))
  }
  else{
    Quant <- quantile(Fit)
    for(j in 0:100){
      Seuil <- Quant[1] + (j*((Quant[5] - Quant[1])/100))
      Misc<-table(Fit >= Seuil, Obs)
      if(TSS!=TRUE) a <- KappaStat(Misc) else a <- TSS.Stat(Misc)
      if(!is.na(a)) if(a > 0) {tab[j+1, 1] <- Seuil; tab[j+1, 2] <- a}
      rm(Misc, Seuil)
    }
    
    t <- max(tab[,2],na.rm=TRUE)
    seuil <- tab[tab[,2]==t,1]   ### Note: Ici il se peut qu'on aie plus de 1 seuil...dans ce cas le plus bas est gard?.
    if(t > 0) {
      Misc<-table(Fit >= seuil[1], Obs)
      TP <- Misc[4]
      TN <- Misc[1]
      ca0 <- (TN * 100)/sum(Misc[,1])
      ca1 <- (TP * 100)/sum(Misc[,2])
      if(TSS!=TRUE) invisible(list(Kappa = t, CutOff = seuil[1], TP = TP, se = ca1, TN = TN, sp = ca0))
      else invisible(list(TSS = t, CutOff = seuil[1], TP = TP, se = ca1, TN = TN, sp = ca0))
    }
    else {
      if(TSS!=TRUE) invisible(list(Kappa = 0, CutOff = 0, TP = 0, se = 0, TN = 0, sp = 0))
      else invisible(list(TSS = 0, CutOff = 0, TP = 0, se = 0, TN = 0, sp = 0))
    }
  }
}

TSS.Stat<-function(Misc)
{
  if(dim(Misc)[1]==1){
    if(row.names(Misc)[1]=="FALSE") Misc<-rbind(Misc, c(0,0))
    else {
      a<-Misc
      Misc<-c(0,0)
      Misc<-rbind(Misc, a)
      
    }
  }
  n <- sum(Misc)
  a <- Misc[1,1]
  b <- Misc[1,2]
  c <- Misc[2,1]
  d <- Misc[2,2]
  sens<-a/(a+c)
  spec<-d/(b+d)
  K <- (sens + spec) - 1        #TSS
  return(K)
}

KappaStat <-  function(Misc)
{
  if(dim(Misc)[1]==1){
    if(row.names(Misc)[1]=="FALSE") Misc <- rbind(Misc, c(0,0))
    else{
      a <- Misc
      Misc <- c(0,0)
      Misc <- rbind(Misc, a)
    }
  }
  n <- sum(Misc)
  n.1 <- sum(Misc[,1])
  n.2 <- sum(Misc[,2])
  n1. <- sum(Misc[1,])
  n2. <- sum(Misc[2,])
  Po <- (1/n) * (Misc[1,1] + Misc[2,2])
  Pe <- ((1/n)^2) * ((as.numeric(n1.) * as.numeric(n.1)) + (as.numeric(n2.) * as.numeric(n.2)))
  K <- (Po - Pe)/(1 - Pe)
  #cat("\n Kappa=", K, "\n")
  return(K)
}

model.eval<-function(pred,obs){
  auc <- Hmisc::somers2(pred,obs)[1]
  SomD<-Hmisc::somers2(pred,obs)[2]
  kappa<-KappaRepet(obs,pred, TSS=TRUE)[c(1,4)]
  TSS<-kappa$TSS
  se<-kappa$se/100
  se.scale<-2*se-1
  Boyce<-ecospat.boyce(pred,obs = pred[obs==1], cor.method='kendall',PEplot=F)$correlation
  my.eval<-round(c(auc,SomD,TSS,se,se.scale,Boyce),3)
  names(my.eval)<-c('auc','SommerD','MaxTSS','MaxTSS_Sensitivity','Scaled_Sensitivity','Boyce')
  return(my.eval)
}

cv.fun <- function(this.fold, data,my.var,RF=TRUE,GBM=TRUE,GBM.param=B,GAM=TRUE,my.k='k=6',
                   ME=TRUE,me.method='maxnet',me.regmul=1,me.classes='lqh',
                   ENSEMBLE = TRUE, eval.lim=0.5){
  
  my.models<-c(RF,GBM,GAM,ME)
  names(my.models)<-c('RF','GBM','GAM','ME')
  stopifnot(length(which(my.models))>0)
  
  fmla <- as.formula(paste("Qobs ~ ", paste(my.var,collapse = "+ ")))
  
  training <- setDT(data)[my.folds != this.fold]
  validation <- setDT(data)[my.folds == this.fold]
  
  eval.rf1<-NA
  if(RF){
    rf1 <- ranger::ranger(fmla, data=training,num.trees=1000)
    pred.rf1 <- predict(rf1,validation,type="response")$predictions
    
    #Evaluation of the rf model
    eval.rf1<-model.eval(pred.rf1,validation$Qobs)
  }
  
  eval.brt<-NA
  if(GBM){  
    # Split data into training and test
    training.gbm <- training
    testing.gbm <- validation
    
    params <- list(booster = "gbtree", objective = "binary:logistic", eta=GBM.param$eta, gamma=GBM.param$gamma, 
                   max_depth=GBM.param$max_depth, min_child_weight=GBM.param$min_child_weight, subsample=GBM.param$subsample, 
                   colsample_bytree=GBM.param$colsample_bytree)
    
    # Using the inbuilt xgb.cv function, let's calculate the best nround for this model. In addition, this function also returns CV error, which is an estimate of test error.
    dtrain <- xgb.DMatrix(as.matrix(training.gbm[,..my.var]), label=training.gbm$Qobs)
    dtest <- xgb.DMatrix(as.matrix(testing.gbm[,..my.var]), label=testing.gbm$Qobs)
    
    xgbcv <- xgb.cv( params = params, data = dtrain, nrounds = 500, nfold = 5, showsd = TRUE, stratified = TRUE, print_every_n = 10, early_stopping_rounds = 20, maximize = TRUE, eval_metric = "auc")
    MIN <- which.max(xgbcv$evaluation_log$test_auc_mean)
    
    #first default - model training
    xgb1 <- xgb.train (params = params, data = dtrain, nrounds = MIN, watchlist = list(val=dtest,train=dtrain),  print_every_n= 10, early_stopping_rounds = 20, maximize = TRUE , eval_metric = "auc")
    
    TEST <-  xgb.DMatrix(as.matrix(validation[,..my.var]),label=validation$Qobs)
    pred.brt <- predict (xgb1,TEST)
    
    ### evaluation BRT
    eval.brt<-model.eval(pred.brt,validation$Qobs)
    
  }
  
  eval.gam<-NA
  if (GAM){
    # GAM model
    fmlaG <- as.formula(paste("Qobs ~ ", paste('s(',my.var,',',my.k,')',collapse = "+ ")))
    BAM <- mgcv::bam(fmlaG,data=training,family=binomial)
    pred.gam <- as.numeric(predict(BAM,newdata=validation,type="response"))
    # rm(list=c("BAM"))
    
    ### evaluation BRT
    eval.gam<-model.eval(pred.gam,validation$Qobs)
  }
  
  eval.m<-NA
  if (ME){
    # MAXENT model
    if (add.para$me.method$me.method=='dismo'){
      ME_regmul<-paste0('betamultiplier=',as.numeric(add.para$me.regmul))
      if(add.para$me.classes=='lqh'){
        my.args=c(ME_regmul,'threshold=false','product=false')
      }else{
        my.args=ME_regmul
      }
      train.m<-maxent(x=training[,..my.var],p=as.data.frame(training$Qobs),args=my.args)
      pred<-predict(train.m,as.data.frame(validation)[,my.var],type='cloglog')
    }
    if (me.method=='maxnet'){
      p=training$Qobs
      data.me=training[,..my.var]
      train.m <- maxnet(p, data.me,f=maxnet.formula(p,data.me,me.classes), regmult=me.regmul)
      pred.m<-predict(train.m,validation,type='cloglog')  
    }
    ### evaluation MAXENT
    eval.m<-model.eval(pred.m,validation$Qobs)
  }  
  ### Average model
  pred.models<-c()
  eval.average<-rep(NA,6)
  eval.w.average<-rep(NA,6)
  
  if (ENSEMBLE & length(my.models>1)){
    if(RF){pred.models<-cbind(pred.models,pred.rf1)}
    if(GBM){pred.models<-cbind(pred.models,pred.brt)}
    if(GAM){pred.models<-cbind(pred.models,pred.gam)}
    if(ME){pred.models<-cbind(pred.models,pred.m)}
    
    pred.mean <- apply(pred.models,1,mean) #average
    
    #weighted average model
    w.coef<-c()
    
    if(RF){w.coef<-c(w.coef,mean(eval.rf1[c(2,3,5,6)]))}
    if(GBM){w.coef<-c(w.coef,mean(eval.brt[c(2,3,5,6)]))}
    if(GAM){w.coef<-c(w.coef,mean(eval.gam[c(2,3,5,6)]))}
    if(ME){w.coef<-c(w.coef,mean(eval.m[c(2,3,5,6)]))} #coeff to weight the model following their performances
    
    w.coef.lim<-w.coef
    w.coef.lim[which(w.coef.lim<eval.lim)]<-0 # Model below user defined threshold have no weight !!
    
    if(sum(w.coef.lim>0)){
      pred.w.mean<-apply(pred.models,1,weighted.mean,w=w.coef.lim)
      
    }else{
      pred.w.mean<-apply(pred.models,1,weighted.mean,w=w.coef)
    }
    
    ### Evaluation of the weighted average model
    eval.w.average<-model.eval(pred.w.mean,validation$Qobs)
    
    ### Evaluation of the average model
    eval.average<-model.eval(pred.mean,validation$Qobs)
  }
  
  my.eval<-rbind(eval.rf1,eval.brt,eval.gam,eval.m,eval.average,eval.w.average)
  return(my.eval)
  
}




model.proj<-function(proj.data,selected.models,cal.data,my.var,add.para,save.dir = getwd(),
                     Npart=Npart){
  
  # # # ## Debug
  # load("G://R/G2-Eaux_dynamiques/data/2_enviro4debug.RData")
  # pacman::p_load(parallel, doParallel , ranger , xgboost , caret , data.table, mgcv , dismo , rJava  , Hmisc , dplyr , maxnet,lightgbm,methods)
  # cal.data=as.data.frame(obs.mod)
  # proj.data = env.pot.mod
  # my.var= VAR
  # Npart=list('RF'= 1, 'GBM' = 15,'GAM' =10,'MAXENT'=5)

  fmla <- as.formula(paste("Qobs ~ ", paste(my.var,collapse = "+ ")))
  if(nrow(cal.data)>5000){
    l.var<-sample(nrow(cal.data),5000)
  }else{
    l.var<-1:nrow(cal.data)
  }
  
  
  
  ### Random forest
  if(selected.models=='RF'){
    my.part<-split(sample(1:nrow(proj.data),nrow(proj.data)),as.factor(1:Npart$RF))
    
    mod <- ranger(fmla, data=cal.data,num.trees=1000) # 0.5 min
    pred<-vector(length=nrow(proj.data))
    
    for (j in 1:Npart$RF){
      print(paste0('RF_proj_partition ',j,' on ',Npart$RF))
      pred[my.part[[j]]] <- predict(mod,proj.data[my.part[[j]],],type="response")$predictions # 1,5 min
    }
    pred.lvar<-predict(mod,cal.data[l.var,],type='response')$predictions
    var.imp<-c()
    for (i in 1:length(my.var)){
      cal.datai<-cal.data[l.var,my.var]
      cal.datai[,i]<-sample(cal.datai[,i])
      var.imp<-c(var.imp,cor(pred.lvar,predict(mod,cal.datai,type='response')$predictions))
    }
  }
  
  ## Boosted regression trees
  if(selected.models=='GBM'){
    inTrain <- createDataPartition(y = cal.data$grid.id, p = 0.80, list = FALSE)  # 80% des données
    training <- cal.data[inTrain,]
    testing <- cal.data[-inTrain,]

    dtrain <- list(data=Matrix::Matrix(as.matrix(training[,my.var]),sparse=T),label=training$Qobs)
    dtest <- list(data=Matrix::Matrix(as.matrix(testing[,my.var]),sparse=T),label=testing$Qobs)
    
    mod <- lgb.load(paste(dest.folder,list.files(dest.folder)[grep("lgb.mod",list.files(dest.folder))],sep="/"))
   
    ## Predictions
    pred<-vector(length=nrow(proj.data))
    pdata <- list(data=Matrix::Matrix(as.matrix(proj.data[,VAR]),sparse=T))
    pred <- predict(mod,pdata$data)
    
    TEST.lvar<- list(data=Matrix::Matrix(as.matrix(as.data.frame(cal.data)[l.var,my.var]),sparse=T),label=l.var)
    pred.lvar<-predict(mod,TEST.lvar$data)
    cal.datai<-cal.data[l.var,my.var]
    var.imp<-c()
      for (i in 1:length(my.var)){
        cal.datai<-cal.data[l.var,my.var]
        cal.datai[,i]<-sample(cal.datai[,i])
        TEST.lvar<- list(data=Matrix::Matrix(as.matrix(cal.datai),sparse=T),label=l.var)
        var.imp<-c(var.imp,cor(pred.lvar,predict(mod,TEST.lvar$data)))
      }
     } ## end of GBM
  
  # GAM model 
  if(selected.models=='GAM'){
    my.part<-split(sample(1:nrow(proj.data),nrow(proj.data)),as.factor(1:Npart$GAM))
    my.k<-add.para$my.k  
    fmlaG <- as.formula(paste("Qobs ~ ", paste('s(',my.var,',',my.k,')',collapse = "+ ")))
    system.time(mod <- mgcv::bam(fmlaG,data=cal.data,family=binomial))
    
    pred<-vector(length=nrow(proj.data))
    
    for (j in 1:Npart$GAM){
      print(paste0('GAM_proj_partition ',j,' on ',Npart$GAM))
      pred[my.part[[j]]] <- as.numeric(predict(mod,proj.data[my.part[[j]],],type="response"))
    }
    
    pred.lvar<-predict(mod,cal.data[l.var,],type='response')
    var.imp<-c()
    for (i in 1:length(my.var)){
      cal.datai<-cal.data[l.var,my.var]
      cal.datai[,i]<-sample(cal.datai[,i])
      var.imp<-c(var.imp,cor(pred.lvar,predict(mod,cal.datai,type='response')))
    }
  } 
  
  # MAXENT model
  if(selected.models=='ME'){
    my.part<-split(sample(1:nrow(proj.data),nrow(proj.data)),as.factor(1:Npart$MAXENT))
    
    if(add.para$me.method=='dismo'){
      ME_regmul<-paste0('betamultiplier=',as.numeric(add.para$me.regmul))
      if(add.para$me.classes=='lqh'){
        my.args=c(ME_regmul,'threshold=false','product=false')
      }else{
        my.args=ME_regmul
      }
      mod<-maxent(x=cal.data[,my.var],p=as.data.frame(cal.data$Qobs),args=my.args) 
      
      pred<-vector(length=nrow(proj.data))
      for (j in 1:Npart$MAXENT){
        print(paste0('MAXENT_proj_partition ',j,' on ',Npart$MAXENT))
        pred[my.part[[j]]]<-predict(mod,as.data.frame(proj.data[my.part[[j]],] ))
      }
      
      pred.lvar<-predict(mod,as.data.frame(cal.data[l.var,]))
      var.imp<-c()
      for (i in 1:length(my.var)){
        cal.datai<-cal.data[l.var,my.var]
        cal.datai[,i]<-sample(cal.datai[,i])
        var.imp<-c(var.imp,cor(pred.lvar,predict(mod,cal.datai,type='response')))
      }
    }
    if(add.para$me.method=='maxnet'){
      p=cal.data$Qobs
      data.me=cal.data[,my.var]
      me.classes<-add.para$me.classes
      me.regmul<-add.para$me.regmul
      mod <- maxnet(p, data.me,f=maxnet.formula(p,data.me,me.classes), regmult=me.regmul)
      
      pred<-vector(length=nrow(proj.data))
      for (j in 1:Npart$MAXENT){
        print(paste0('MAXENT_proj_partition ',j,' on ',Npart$MAXENT))
        pred[my.part[[j]]]<-predict(mod,as.data.frame(proj.data[my.part[[j]],] ),type='cloglog')
      }
      pred<-predict(mod,proj.data,type='cloglog')
      
      pred.lvar<-predict(mod,as.data.frame(cal.data[l.var,]),type='cloglog')
      var.imp<-c()
      for (i in 1:length(my.var)){
        cal.datai<-cal.data[l.var,my.var]
        cal.datai[,i]<-sample(cal.datai[,i])
        var.imp<-c(var.imp,cor(pred.lvar,predict(mod,cal.datai,type='cloglog')))
      }
    }} # end of maxent

  var.imp<-1-var.imp
  full.model<-list(model=mod,varImp = var.imp,guild.pred =pred)
  save(full.model,file = paste0(save.dir,"/mod_",selected.models,'.Rdata'))
  full.model
}


#=================================================================================================================================
# ENDS of Functions for the model evaluations
#=================================================================================================================================


mod.table<-function(GAM=T,ME=T,GBM=T,cv=cv,gam.smoothing=my.k,
                    ME.method='maxnet',ME.regmul=my.regmul,ME.classes=me.classes
){
  t.table<-matrix(nrow=0,ncol=8)
  
    for (i in 1:length(cv)){
      t.tablei<-c(paste0('RF_cv_',cv[i]),'RF',cv[i],rep('',5))
      t.table<-rbind(t.table,t.tablei)
    }
  
  if(any(grepl("GBM",names(PAR)))|GBM==T){
    for (i in 1:length(cv)){
        t.tablei<-c(paste0('GBM_cv_',cv[i]),'GBM',cv[i],1,rep('',4))
        t.table<-rbind(t.table,t.tablei)
    }
  }
  
  if(any(grepl("GAM",names(PAR)))|GAM==T){
    for (i in 1:length(cv)){
      for (j in 1:length(gam.smoothing)){
        t.tablei<-c(paste0('GAM_cv_',cv[i]),'GAM',cv[i],'',gam.smoothing[j],rep('',3))
        t.table<-rbind(t.table,t.tablei)
      }
    }
  }

  if(any(grepl("ME",names(PAR)))|ME==T){
    for (i in 1:length(cv)){
      if(ME.method=='dismo'){
        for (j in 1:length(ME.classes)){
          for(z in 1:length(ME.regmul)){
            t.tablei<-c(paste0('ME_cv_',cv[i]),'ME',cv[i],rep('',2),'dismo',ME.regmul[z],ME.classes[j])
            t.table<-rbind(t.table,t.tablei)
          }
        }
      }
      if(ME.method=='maxnet')
        for (j in 1:length(ME.classes)){
          for(z in 1:length(ME.regmul)){
            t.tablei<-c(paste0('ME_cv_',cv[i]),'ME',cv[i],rep('',2),'maxnet',ME.regmul[z],ME.classes[j])
            t.table<-rbind(t.table,t.tablei)
          }
        }
    }
  }
  colnames(t.table)<-c('mod_name','model','cv','GBM_par','GAM_k','ME_method','ME_regmul',
                       'ME_classes')
  row.names(t.table)<-c()
  return(t.table) 
}

eval.mod<-function(row.number,mod_table,data=obs.mod,my.var=VAR,GBM.param){
  
  # ### ## Debug
  # mod_table=my_table
  # data=obs.mod
  # my.var=VAR
  # row.number=6
  # GBM.param = PAR$GBM_para

  mod_table<-as.list(mod_table[row.number,])
  training <- setDT(data)[my.folds != mod_table$cv]
  validation <- setDT(data)[my.folds == mod_table$cv]
  
  if(mod_table$model=='RF'){
    fmla <- as.formula(paste("Qobs ~ ", paste(my.var,collapse = "+ ")))
    rf1 <- ranger::ranger(fmla, data=training,num.trees=1000)
    pred <- predict(rf1,validation,type="response")$predictions
    
    #Evaluation of the rf model
    eval<-model.eval(pred,validation$Qobs)
  }
  
  if (mod_table$model=='GAM'){
    # GAM model
    fmlaG <- as.formula(paste("Qobs ~ ", paste('s(',my.var,',',mod_table$GAM_k,')',collapse = "+ ")))
    BAM <- mgcv::bam(fmlaG,data=training,family=binomial)
    pred<- as.numeric(predict(BAM,newdata=validation,type="response"))
    
    ### evaluation BRT
    eval<-model.eval(pred,validation$Qobs)
  }
  
  if (mod_table$model=='GBM'){
    dtrain <- lgb.Dataset(as.matrix(training[,..my.var]),label=training$Qobs)
    dtest <- lgb.Dataset(as.matrix(validation[,..my.var]),label= validation$Qobs)
    test <- list(data=Matrix::Matrix(as.matrix(validation[,..my.var]),sparse=T),label= validation$Qobs)
    
    # GBM model
    gbm_mod <- lgb.train(
      list(objective         = "binary",
           metric            = "auc",
           learning_rate     = 0.1,
           min_child_samples = 100,
           max_bin           = 100,
           subsample_freq    = 1,
           num_leaves        = GBM.param$num_leaves,
           max_depth         = GBM.param$max_depth,
           subsample         = GBM.param$subsample,
           colsample_bytree  = GBM.param$colsample_bytree,
           min_child_weight  = GBM.param$min_child_weight,
           scale_pos_weight  = GBM.param$scale_pos_weight
           ),
      dtrain,
      valids = list(validation = dtest),
      nthread = 4, 
      nrounds = 5, # increase/ decrease rounds
      verbose= 1, 
      early_stopping_rounds = 2
    )
    
    lgb.save(gbm_mod, paste0(dest.folder,"/",GUILD,"_lgb.mod.txt"))
    pred <-  predict(gbm_mod,test$data)
    
    ### evaluation BRT
    eval<- model.eval(pred,validation$Qobs)
  }
  
  if (mod_table$model=='ME'){
    # MAXENT model
    if (mod_table$ME_method=='dismo'){
      ME_regmul<-paste0('betamultiplier=',as.numeric(mod_table$ME_regmul))
      if(mod_table$ME_classes=='lqh'){
        my.args=c(ME_regmul,'threshold=false','product=false')
      }else{
        my.args=ME_regmul
      }
      train.m<-maxent(x=training[,..my.var],p=as.data.frame(training$Qobs),args=my.args)
      pred<-predict(train.m,as.data.frame(validation)[,my.var],type='cloglog')
    }
    if (mod_table$ME_method=='maxnet'){
      p=training$Qobs
      data.me=training[,..my.var]
      ME_regmul<-as.numeric(mod_table$ME_regmul)
      train.m <- maxnet(p, data.me,f=maxnet.formula(p,data.me,mod_table$ME_classes), regmult=ME_regmul)
      pred<-predict(train.m,validation,type='cloglog')  
    }
    ### evaluation MAXENT
    eval<-model.eval(pred,validation$Qobs)
  }
  model.evaluation<-list(eval=eval,mod.parameter = mod_table,prediction=pred)
}


synt.eval<-function(evali){
  mean_eval<-mean(evali$eval[c(2,3,5,6)])
  synthi<-c(evali$mod.parameter,evali$eval,mean_eval)
  return(synthi)
}


ensemble.eval<-function(data,eval.synth,model_eval,eval.lim=0.5){
  eval.synth<-as.data.table(tidyr::unnest(eval.synth,cols = colnames(eval.synth)))
  
  mod.mean_evaluation<-eval.synth[,.(mean(auc),mean(SommerD),mean(MaxTSS),mean(MaxTSS_Sensitivity),mean(Scaled_Sensitivity),
                                     mean(Boyce),mean(mean_evaluation)),
                                  by=.(model)]
  colnames(mod.mean_evaluation)[2:8]<-colnames(eval.synth)[9:15]
  
  mod.sd_evaluation<-eval.synth[,.(sd(auc),sd(SommerD),sd(MaxTSS),sd(MaxTSS_Sensitivity),sd(Scaled_Sensitivity),
                                   sd(Boyce),sd(mean_evaluation)),
                                by=.(model)]
  colnames(mod.sd_evaluation)[2:8]<-colnames(eval.synth)[9:15]
  
  mod.preds<-c()
  for(i in unique(eval.synth$model)){
    preds<-model_eval[which(eval.synth$model==i)]
    preds<-do.call(c,sapply(preds, function(x) x$prediction)) 
    mod.preds<-cbind(mod.preds,preds)  
    }
  colnames(mod.preds)<-unique(eval.synth$model)
  
  ### Average model
  eval.average<-rep(NA,6)
  eval.w.average<-rep(NA,6)
  w.coef<-mod.mean_evaluation$mean_evaluation
  w.coef.lim<-w.coef
  # w.coef.lim[which(w.coef.lim<eval.lim|is.na(w.coef))]<-0 # Model below user defined threshold have no weight !!
  names(w.coef)<-mod.mean_evaluation$model
  
  pred.mean <- apply(mod.preds,1,mean) #average
  
  # if(sum(w.coef.lim)>0){
  if(sum(w.coef.lim)>0){   ## modifi? par Ervan (04/10)
    pred.w.mean<-apply(mod.preds,1,weighted.mean,w=w.coef.lim)
    
  }else{
    pred.w.mean<-apply(mod.preds,1,weighted.mean,w=w.coef)
  }
  
  cv.order<-sort(data$my.folds)
  eval.w.average<-c()
  eval.average<-c()
  
  for (i in levels(factor(data$my.folds))){
    pred.mean.i<-pred.mean[which(data$my.folds==i)]
    pred.w.mean.i<-pred.w.mean[which(data$my.folds==i)]
    
    ### Evaluation of the weighted average model
    eval.w.averagei<-model.eval(pred.w.mean[which(cv.order==i)],data$Qobs[which(data$my.folds==i)])
    eval.w.averagei<-c(eval.w.averagei,mean(eval.w.averagei[c(2,3,5,6)]))
    eval.w.average<-rbind(eval.w.average,eval.w.averagei)
    
    ### Evaluation of the average model
    eval.averagei<-model.eval(pred.mean[which(cv.order==i)],data$Qobs[which(data$my.folds==i)])
    eval.averagei<-c(eval.averagei,mean(eval.averagei[c(2,3,5,6)]))
    eval.average<-rbind(eval.average,eval.averagei)
  }
  mod.mean_evaluation<-rbind(mod.mean_evaluation,as.data.frame(rbind(c('mean',apply(eval.average,2,mean)),
                                                                     c('w.mean',apply(eval.w.average,2,mean)))),use.names=FALSE)
  mod.sd_evaluation<-rbind(mod.mean_evaluation,as.data.frame(rbind(c('mean',apply(eval.average,2,sd)),
                                                                   c('w.mean',apply(eval.w.average,2,sd)))),use.names=FALSE)
  
  return(list('mean_evaluation'=mod.mean_evaluation,'sd_evaluation'=mod.sd_evaluation,'model_weight'=w.coef.lim))
  
}

## Benchmarking Ha to add
## Benchmarking via Qprop: on enlève les 10% des BV avec un espace prédit trop faible (= 10 ha)
pick.bench<-function(Qpred,haQprop,haSp){
  if (is.na(haSp)){
    b<-haQprop
  }else{
    if (Qpred < haSp){
      b<-Qpred
    }else{
      b<-haSp
    }
  }
  return (b)
}


BENCH<-function(PRED,OBS,th.pred,th.bench, th.qprop=1, defrag,sd.min=TRUE) {
  th.qpred<-quantile(setDT(PRED)[PRED$cal==1,]$Wmean,1-th.pred)
  PRED[,Qp:=0]
  PRED$Qp[which(PRED$Wmean >=th.qpred & is.na(PRED$ist))]<-1
  # ch.realisation<-setDT(PRED)[,.("Qobs"=length(which(Qobs==1)),"Qpred"=length(which(Qpred==1 | Qobs==1)),"EG"=length(Qobs),"BV_area"=round(unique(area))),by="BV04"]
  ch.realisation<-setDT(PRED)[,.("Qobs"=length(which(!is.na(ist))),"Qpred"=length(which(Qp==1)),"EG"=.N,"BV_area"=round(unique(area))),by="BV04"]
  ch.realisation <- ch.realisation[which(ch.realisation$Qp>10),]
  min.ch.bench <- round( 100 *median(ch.realisation[,Qobs / (Qobs + Qpred)],na.rm=TRUE),1) # median realisation rate in CH
  PRED[,"quality":=as.character(NA)]
  PRED <- within(PRED,quality[!is.na(ist)]<-'obs')
  PRED <- within(PRED,quality[which((Qp==1 & is.na(ist)))]<-'pred')
  
  res <- data.frame("BV04"=NA,"Qobs"=NA,"Qpred"=NA,"EG"=NA,"BV_area"=NA,"Qprop"=NA,"bench"=NA,"ha_to_add_Qprop"=NA,"clust"=NA,"subreg"=NA)
  for (cl in levels(factor(PRED$CLUST))){
    SUB <- data.table::setDT(PRED)[CLUST==cl,]
    TAB <- data.table::setDT(SUB)[,.("Qobs"=length(which(quality=="obs")),"Qpred"=length(which(quality=="pred")),"EG"=.N,"BV_area"=round(unique(area))),by="BV04"]
    # if(sum(TAB$Qpred)<10){next}
    ha_available <- TAB$Qpred
    TAB$Qprop = round(100*TAB$Qobs/(TAB$Qpred+TAB$Qobs),1)
    TAB <- TAB[order(TAB$Qprop,decreasing=TRUE),]
    TAB$is.bench <- 0
    TAB <- within(TAB,is.bench[(Qpred+Qobs)>10] <- 1)
    TAB$bench = round(quantile(TAB[is.bench==1,]$Qprop,th.qprop,na.rm=TRUE),1)
    if (is.na(TAB$bench[1]) | TAB$bench[1]<min.ch.bench){
      TAB$bench<-min.ch.bench
    }
    # quantile(cumsum(rep(TAB$Qprop,TAB$Qpred)/sum(TAB$Qpred)),1)
    target<- ceiling((TAB$bench/100)  * (TAB$Qobs + TAB$Qpred))
    TAB$ha_to_add_Qprop <- target-TAB$Qobs
    TAB$ha_to_add_Qprop[TAB$ha_to_add_Qprop<0]<-0
    TAB$clust=cl
    TAB$subreg <- unique(SUB$subreg)
    TAB <- TAB[,c("BV04","Qobs","Qpred","EG","BV_area","Qprop","bench","ha_to_add_Qprop","clust","subreg")]
    res <- rbind(res,TAB)
  }
  res <- res[-1,]
  
  
  ## Benchmarking via pool d'espèces
  ### Filter obs to fit with IST
  if (length(grep("subreg",names(OBS)))>0) { OBS$subreg <- NULL }
  OBS <- merge(OBS,PRED[!is.na(ist) & !duplicated(grid.id),c("CNHA","CLUST","canton","subreg","BV04")],by="CNHA")
  
  div <- data.table::setDT(OBS)[,.(N=length(unique(taxonid[which(!is.na(taxonid))])),NG=length(unique(group[which(!is.na(group))])),Qobs=ifelse(sum(w[!duplicated(taxonid)],na.rm=T)>=1,1,0)),by=.(CNHA)]
  
  ## Liste d'especes dans le cluster
  OBS <- subset(OBS,!is.na(OBS$CLUST)) #remove non-assigned pixel
  LIST <- OBS[,.("list"=unique(name[which(!is.na(name))])),by=.(CLUST,BV04,subreg)]
  CLU <- LIST[,unique(list),by=.(CLUST)]  # species by cluster
  REG <- LIST[,unique(list),by=.(subreg)]
  Y <- LIST[,.(CLUST=unique(CLUST),subreg=unique(subreg),sp=unique(list)),by=BV04]
  Nsp.ch<-length(unique(Y$sp))
  Nsp.BV<-Y[,.('Nsp'=length(unique(sp)),'subreg'=unique(subreg), 'CLUST'=unique(CLUST)),by=BV04]
  Nsp.cl<-Y[,.('Nsp'=length(unique(sp)),'subreg'=unique(subreg)),by=CLUST]
  Nsp.bioreg<-Y[,.('Nsp'=length(unique(sp))),by=subreg]
  Nsp.synthesis<-merge(
    merge(Nsp.BV,Nsp.cl[,c('Nsp','CLUST')],by='CLUST'),
    Nsp.bioreg[,c('Nsp','subreg')],by ='subreg')
  colnames(Nsp.synthesis)[4:6]<-c('Nsp_BV','Nsp_CLUST','Nsp_BIOREG')
  empty.bv<-merge(Nsp.synthesis,unique(PRED[,'BV04']),all.y=TRUE,by='BV04')
  Sp.prop<-setDT(Nsp.synthesis)[,.(BV04,CLUST,subreg,'prop_SpClust'=round(Nsp_BV/Nsp_CLUST,3), 
                                   'prop_SpBioReg'= round(Nsp_BV/Nsp_BIOREG,3))]
  
  bench.sp.prop<-Sp.prop[,.('bench_Sp_cluster'=quantile(prop_SpClust,th.bench),'bench_Sp_bioreg'=quantile(prop_SpBioReg,th.bench)),by='CLUST']
  bench.sp.prop<-bench.sp.prop[,-'bench_Sp_cluster']
  min.bench.spbioreg<-quantile(Sp.prop$prop_SpBioReg,0.5)
  bench.sp.prop$bench_Sp_bioreg[bench.sp.prop$bench_Sp_bioreg[]<min.bench.spbioreg]<-min.bench.spbioreg
  bench.sp<-merge(merge(Nsp.synthesis,Sp.prop[,-c('CLUST','subreg')],by='BV04'),bench.sp.prop,all.x=TRUE,by='CLUST')
  bench.sp[,bench_Nsp := round(bench_Sp_bioreg * Nsp_BIOREG)]
  
  clust.synth<-bench.sp[,.(unique(subreg),unique(Nsp_CLUST),unique(Nsp_BIOREG),unique(bench_Sp_bioreg),unique(bench_Nsp)),by=CLUST]
  
  BV.list<-PRED[,.(CLUST=unique(CLUST), subreg=unique(subreg),area=unique(area)),by='BV04']
  
  bench.sp<-merge(BV.list[,c('BV04','CLUST','subreg')],bench.sp[,-c('CLUST','subreg')],all.x=TRUE,by='BV04')
  bench.sp$Nsp_BV[which(is.na(bench.sp$Nsp_BV))]<-0
  bench.sp$prop_SpBioReg[which(is.na(bench.sp$prop_SpBioReg))]<-0
  bench.sp$prop_SpClust[which(is.na(bench.sp$prop_SpClust))]<-0
  bench.sp<-bench.sp[order(bench.sp$BV04),]
  
  asd<-merge(bench.sp[,c('CLUST','BV04')],clust.synth,all.x=TRUE,by='CLUST')
  asd<-asd[order(asd$BV04),]
  
  bench.sp[,c('Nsp_CLUST','Nsp_BIOREG','bench_Sp_bioreg', 'bench_Nsp')]<-asd[,c(4,5,6,7)]
  # bench.sp <- bench.sp[!is.na(CLUST),]
  L <- length(unique(bench.sp$CLUST))
  AA = 0
  ha_to_add_sp<-data.table('BV04'=NA,'ha_to_add_sp'=NA)
  ID.CLUST <- na.omit(unique(factor(bench.sp$CLUST)))
  for (i in ID.CLUST){
    AA = AA+1
    cat(paste0(L-AA,"-"))
    pop.max<-length(PRED$CLUST[PRED$CLUST==i])
    sub <- setDT(OBS)[CLUST==i,]
    tab <- as.data.frame.matrix(table(sub$CNHA,sub$name))
    if (nrow(tab)<5|ncol(tab)<2){
      ha_to_add_sp<-rbind(ha_to_add_sp,data.frame(BV04=unique(PRED[PRED$CLUST==i,]$BV04),ha_to_add_sp=NA))
    } else{
      
      curve <- vegan::specaccum(dplyr::sample_n(tab,pop.max,replace=TRUE),method="exact") # = courbe moyenne
      if (sd.min==TRUE){
        curve$richness<-round(curve$richness-curve$sd)
      }
      curve$richness[curve$richness<0]<-0
      
      delta.sites<-0
      
      
      if (max(bench.sp[CLUST==i]$bench_Nsp) > max(bench.sp[CLUST==i]$Nsp_CLUST)){
        sub.bioreg<-setDT(OBS)[subreg==sub$subreg[1],]
        tab.bioreg<-as.data.frame.matrix(table(sub.bioreg$CNHA,sub.bioreg$name))
        curve.bioreg<-vegan::specaccum(dplyr::sample_n(tab.bioreg,pop.max,replace=TRUE),method="exact") # = courbe moyenne pour la biorégion
        if (sd.min==TRUE){
          curve.bioreg$richness<-curve.bioreg$richness-curve.bioreg$sd
        }
        curve.bioreg$richness<-round(curve.bioreg$richness)
        curve.bioreg$richness[curve.bioreg$richness<0]<-0
        if (curve.bioreg$richness[1]>max(bench.sp[CLUST==i]$Nsp_CLUST)){
          curve.bioreg$richness[1]<-max(bench.sp[CLUST==i]$Nsp_CLUST)
        }
        if (curve.bioreg$richness[1]>max(bench.sp[CLUST==i]$bench_Nsp)){
          curve.bioreg$richness[1]<-max(bench.sp[CLUST==i]$bench_Nsp)
        }
        delta.sp<- c(max(bench.sp[CLUST==i]$Nsp_CLUST), max(bench.sp[CLUST==i]$bench_Nsp))
        delta.sites<-max(which(curve.bioreg$richness<=delta.sp[2]))-max(which(curve.bioreg$richness<delta.sp[1]))
      }
      
      tabi<-bench.sp[bench.sp$CLUST==i,]
      tabi$ha_to_add_sp <- as.numeric(NA)
      for (j in 1:nrow(tabi)){
        if (length(which(curve$richness>tabi[j,]$bench_Nsp))==0){ # if the max of the curve is lower than benchmark
          b<-min(which(curve$richness>=max(curve$richness)))
        }else{
          b<-min(which(curve$richness>=tabi[j,]$bench_Nsp))
        }
        if (length(which(curve$richness>tabi[j,]$Nsp_BV))==0){ # if the max of the curve is lower than Nsp in the BV
          a<-min(which(curve$richness==max(curve$richness)))
        }else{
          a<-min(which(curve$richness>=tabi[j,]$Nsp_BV))
          
        }
        delta.sites.j<-delta.sites
        if(tabi[j,]$Nsp_BV==0){delta.sites.j<-delta.sites+1}
        tabi$ha_to_add_sp[j] <- (b-a)+ delta.sites.j
        
      }
      
      ha_to_add_sp<-rbind(ha_to_add_sp,tabi[,c("BV04","ha_to_add_sp")])
    }
  }
  
  ha_to_add_sp <- ha_to_add_sp[-1,]
  ha_to_add_sp$ha_to_add_sp[which(ha_to_add_sp$ha_to_add_sp<0)]<-0
  bench.sp<-merge(bench.sp,ha_to_add_sp,all.x=TRUE,by='BV04')
  
  bench.final<-merge(bench.sp,setDT(res)[,c("subreg","clust"):=NULL],all=TRUE,by='BV04')
  bench.final[,ha_to_add_sp := pick.bench(Qpred,ha_to_add_Qprop,ha_to_add_sp),by='BV04']
  
  #### Fragmentation part -> to put into a function
  # PRED$centro<-PRED %>% 
  #   st_as_sf() %>% 
  #   st_geometry() %>% 
  #   st_centroid()

  ist.centro<-st_coordinates(st_geometry(PRED[!is.na(ist),]$centro)) # centroid for the observations
  pred.centro<-st_coordinates(st_geometry(PRED[quality=='pred',]$centro)) # centroid for the predictions
  set.seed(1)
  system.time(nn<-FNN::get.knnx(ist.centro,pred.centro,algo='kd_tree', k=1)) # 1 calculate the distance between potential sites and the closest observation
  mypts<-which(nn$nn.dist>=defrag) #remove all potential sites closer to 400 m from observations
  
  my.pts<-sp::SpatialPointsDataFrame(coords=pred.centro[mypts,], data = PRED[quality=='pred',][mypts,c('CNHA','BV04')]) # format a spatial data frame
  system.time(my.pts.dsg<-sp::remove.duplicates(my.pts,zero = defrag)) # 587 sec disaggregate potential sites with a minimal distance
  PRED.dsg<-PRED[PRED$CNHA %in% my.pts.dsg$CNHA,] # Disaggregation of the PRED table
  ## Export data for priorisation
  TEMP <- PRED.dsg[,c("CNHA","grid.id","centro")]
  save(TEMP,file=paste0(getwd(),"/data/pred_dsg",defrag,".Rdata"))
  
  ha2add.defrag<-PRED.dsg[,.('ha2add_defrag'=as.numeric(.N)),by='BV04'] # number of disagregated ha for each BV

  bench<-merge(ha2add.defrag,bench.final,by='BV04',all=TRUE) %>% 
    .[,ha2add_defrag := ifelse(is.na(ha2add_defrag),0, ha2add_defrag),by=BV04] %>% 
    .[,EB_sp_defrag:=max(c(ha2add_defrag,ha_to_add_sp),na.rm=T),by=BV04] %>% 
    .[,EB_sp_defrag_weighted := ifelse(EB_sp_defrag > ha_to_add_Qprop, ha_to_add_Qprop, EB_sp_defrag)] %>% 
    .[,EB_sp_defrag_weighted := ifelse((Qobs > 0) && (EB_sp_defrag_weighted < ha2add_defrag), ha2add_defrag, EB_sp_defrag_weighted),by=BV04] %>% 
    merge(setDT(BV)[,c('BV04','geometry')],by='BV04',all=TRUE)  %>% 
    st_as_sf()
  
  bench <- within(bench,EB_sp_defrag[is.na(EB_sp_defrag)] <- 0)
  bench <- within(bench,EB_sp_defrag_weighted[is.na(EB_sp_defrag_weighted)] <- 0)
  bench$delta<-bench$EB_sp_defrag - bench$EB_sp_defrag_weighted
  return(list(bench,PRED))
}

## Replace NA by 0 for landscape variables
NArep = function(dt) {
  ID <- match(VAR,TT$label[39:55])
  if (any(ID)){
    na.replace = function(v,value=0) { v[is.na(v)] = value; v }
    for (i in VAR[!is.na(ID)])
      eval(parse(text=paste("dt[,",i,":=na.replace(",i,")]")))
  }}


area2add <- function(DIV) {
site2area<-function(nsites,nrep=100,my.IST){
  repi<-vector(length=nrep)
  for (i in 1:nrep){
    repi[i]<-length(unique(unlist(my.IST[sample(1:length(my.IST),nsites)])))
  }
  return(list(MeanHA= mean(repi), SDha = sd(repi)))
}

my.xy<-st_centroid(DIV$geometry)
my.xy.pol<-st_buffer(my.xy,dist=142,nQuadSegs=4)

my.nn<-st_intersects(my.xy.pol,my.xy) # list all neighbouring cells with observations for each cell
# for (i in 1:length(my.nn)){           # list all neighbouring cells with observations for each cell
#   my.nn[[i]]<-my.nn[[i]][-which(my.nn[[i]]==i)]
# }

my.IST<-my.nn[which(DIV$QUAL==1)] #list all neighbouring cells with observations for each QUALITY cell
DIV$ist <- NA
DIV$ist[unique(unlist(my.nn[which(DIV$QUAL==1)]))] <- 1
dist.percentile<-quantile(1:length(my.IST),c(seq(.01,1,by =0.01))) #Percentiles of the number of IST in Switzerland

my.coef<-matrix(nrow=length(dist.percentile),ncol=2) #build the coefficient "sites/area" in regards to the saturation of the available space
for(i in 1:length(dist.percentile)){
  a<-site2area(dist.percentile[i],100,my.IST)
  my.coef[i,1]<-a$MeanHA
  my.coef[i,2]<-a$SDha
  cat(paste0(i,"-"))
}

my.coef<-data.frame(cbind(dist.percentile,my.coef))
my.coef<-cbind(my.coef,my.coef[,2]/my.coef[,1])
my.coef<-cbind(my.coef,my.coef[3]/my.coef[,1])

colnames(my.coef)<-c('Nsites','Nha_mean','Nha_sd','coef_sites_area','coef_sites_area_sd')
my.coef <- rbind(my.coef[1,],my.coef)
rownames(my.coef) <- names(quantile(1:length(my.IST),c(seq(0,1,by =0.01))))
write.table(my.coef,file=paste0(getwd(),"/data/coef_",GUILD,".txt"),sep='\t',row.names=T,quote=F)
list(DIV,my.coef)
}



#function to select quickly neighbouring cells. It returns an index corresponding to the selected points
#x : a querry point
#y : a dataset of points
#dist : neighbouring distance for selection
#k : max number of nearest neighbour (if dist = 142 -> k = 9; if dist = 284 k = 25)
st_select<-function(x,y,dist,k){
  x<-matrix(x,ncol=2)
  nn<-get.knnx(y,x,k,algo='brute')
  mypts<-nn$nn.index[which(nn$nn.dist<=dist)]
  return(mypts)
}


#function to analyse the impact of adding one cell with observation on the landscale
#x : index of the site (for apply and lapply function)
#my.data : guild space with fields 'ist' 'qual' and 'QUAL'
#'QUAL' = sites avec qualit? 1 (= ayant au-moins 2 ha avec des obs autour ou 3 ha de qualit? dans la fen?tre) (cod? 1)
#'ist' = 'ist.qual' + satellites (=cellules en p?riph?ries) (cod? 1)
#'qual' = nbr de cellules avec des obs ds la fen?tre de 300x300 (cod? de 0 ? 9)
# my.data.xy : xy coordinates of my data 

lv<-function(x,my.data,my.data.xy){
  my.sites<-st_select(my.data.xy[x,],my.data.xy,284,k=25) # delimit the size of the window for the analysis, i.e 5 by 5 ha
  my.sites.data<-my.data[my.sites,] # select the sites of the window in the dataset
  if (length(which(!is.na(my.sites.data$QUAL)))<1){
    return (c(0,1,0)) # if there is no observation, no need for analysis
  }else{
    my.new.sites<-my.sites.data
    my.new.sites$QUAL[1]<-0 # simulate that the center of the window is an observation
    
    N<-list() 
    for (i in which(!is.na(my.new.sites$QUAL))){ # select the neighbouring observations for each observation within the window
      N[[i]]<-st_select(x=my.data.xy[my.sites[i],],
                        y= matrix(my.data.xy[my.sites[which(!is.na(my.new.sites$QUAL))],],ncol=2),
                        dist =  142, k = nrow(matrix(my.data.xy[my.sites[which(!is.na(my.new.sites$QUAL))],],ncol=2)))
    }
    qual.sites<-lengths(N) # compute the quality
    sites<-which(lengths(N)>2) # count the number of sites with more than 2 neighbouring observations
    if (length(sites)>0){in.ist<-unique(unlist(N[sites]))} # detect the observation to include in the IST (core and peripheal)
    
    my.new.sites$qual[!is.na(my.new.sites$QUAL)]<-qual.sites[qual.sites[]>0] # update the values with the simulated observation
    if (length(sites)>0){
      my.new.sites$ist[which(!is.na(my.new.sites$QUAL))[in.ist]]<-1 # update the values with the simulated observation
      my.new.sites$QUAL[sites]<-1 # update the values with the simulated observation
    }
    
    area.gain<-sum(my.new.sites$ist,na.rm=T)-sum(my.sites.data$ist,na.rm=T) # measurement of the additionnal area in the IST
    qual.gain<-sum(apply(cbind(my.new.sites$qual,my.sites.data$qual)[!is.na(my.new.sites$qual),],1,max,na.rm=T))-sum(my.sites.data$qual,na.rm=T) # measurement of the additionnal quality in the IST
    site.gain<-sum(my.new.sites$QUAL,na.rm=T)-sum(my.sites.data$QUAL,na.rm=T) # measurement of the additional sites in the IST
    
    return(c(area.gain,qual.gain,site.gain))
  }
}


col.bin <- function(res2,bench,min.size=5,silence=NULL) {
  res2 <- res2[,c("BV04","ha2add_defrag", "CLUST", "subreg", "Nsp_BV", "Nsp_CLUST", "Nsp_BIOREG", 
                "prop_SpClust", "prop_SpBioReg", "bench_Sp_bioreg", "bench_Nsp", 
                "ha_to_add_sp", "Qobs", "Qpred", "EG", "BV_area", "Qprop", "bench", 
                "ha_to_add_Qprop", "EB_sp_defrag", "EB_sp_defrag_weighted","delta", "geometry")]
  
  names(res2) <- c("BV.id","EB_defrag","CLUST", "subreg", "sp.in.BV", "Nsp_CLUST", "Nsp_BIOREG", "prop_SpClust", "prop_SpBioReg", "bench_Sp_bioreg", "sp.in.bench", "EB_sp", "observed_qual", "potential_qual", "EG", "BV_area", "Qprop", "bench", "ha_to_add_Qprop", "EB", "EB_weighted", "delta", "geometry")
  
  COL <- colorRampPalette(c("#7fa046","blue"))
  COLbin <- COL(8)[seq(2,8,2)]
  TARGET <- setDT(res2)[,get(..bench)]
  BIN <- c(0,10,50,99.9,max(TARGET,na.rm=T))
  if(min(TARGET[res2$potential_qual>=min.size])==0){
    LAB <- c("0-10","10-50","50-100",">100")
  } else {
    LAB <- c("1-10","10-50","50-100",">100") 
  }
  res2$bin <- cut(TARGET,breaks=BIN,include.lowest=TRUE,labels=LAB)
  if(silence==FALSE) {  print(table(res2$bin)) }
  res2$bin <- factor(res2$bin,levels=LAB)
  res2$bin[res2$potential_qual<min.size] <- NA  # tag BV with low potential
  COLtab <- data.frame(bin=LAB,col=COLbin,lab=LAB)
  if(any(is.na(res2$bin))){
    COLtab[5,] <- c(bin=NA,col="#FF8C00",lab="no to low potential")
  }
  res3 <- plyr::join(res2,COLtab[,1:2],by="bin")
  res3$bin <- factor(res3$bin,levels=LAB)
  COLrgb <- leaflet::colorFactor(COLtab$col,res3$bin)
  st_geometry(res3) <- "geometry"
  res3$Erganzungsbedarf <- TARGET
  res3 <- res3[order(res3$BV.id),]
  return(list(res3,COLtab))
}

## TEST