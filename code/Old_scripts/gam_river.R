library(raster)
env <- stack("data/guild2/shp/env100.with.river.tif")
nlayers(env)
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
names(env) <- NAME.VAR

# Add variables one by one
library("sf")
load("C:/Dossier_Ervan/R/infraeco_git/data/guild2/buffer_river_200m_filtered_with_BV_polygons.Rdata")

a <- dim(buf)[2]
for( i in 1:(nlayers(env)-4)) {
  buf[,a+i]  <- extract(env[[i]],buf$grid.id)
  names(buf)[a+i] <- NAME.VAR[i]
}
save(buf,file="data/guild2/river_buffer200m_with_enviro_sf.Rdata")
rm(env)

# # Explore important variables by groups
# div <- setDT(obs_sf)[year_coll>=2000 & coord_dev<=300,.(N=length(unique(taxonid))),by=.(grid.id,group)] 
# div$plant <- ifelse(div$group=="PLANT",1,0)
# head(div[grid.id%in%div[duplicated(grid.id),]$grid.id,][order(grid.id)])
# BV <- left_join(riv ,as.data.frame(div),by="grid.id")
# 
# tt <- dcast(setDT(BV)[,plant,by=BV],BV ~ plant, length)
# tt$sumrow <- rowSums(tt[,2:3])
# ID.BV <- tt[which.max(tt$sumrow),]$BV
# BVm <- melt(BV %>% filter(BV==ID.BV),id.vars=c("grid.id","plant"),measure.vars=c("N"))
# BV2 <- left_join(riv %>% filter(BV==ID.BV) ,as.data.frame(BVm),by="grid.id")
# 
# library(tmap)
# tmap_mode("view")
# 
# tm_shape(BV2[!is.na(BV2$plant),]) +
#   tm_facets(by = "plant") +
#   tm_fill(
#     col = "value",
#     palette = "Reds",
#     style = "cont",
#     contrast = c(0.1, 1),
#     title = "Species richness",
#     id = "grid.id",
#     showNA = T,
#     colorNA = "gold2",
#     alpha = 0.7,
#     popup.vars = c(
#       "Species richness" = "value"
#     )) +
#   tm_layout(legend.stack="vertical",legend.outside = TRUE) +
#   tm_borders(col = "gray", lwd = 0.2)  
# +
# tm_basemap(server="Esri.WorldImagery")

# ## DMFA :
# library(FactoMineR)
# data.frame(1:ncol(BV),names(BV))
# test <- BV[!is.na(plant),c(4:41,46)]
# test <- na.omit(test)
# test[,1:38] <- data.frame(apply(test[,1:38],2,scale))
# 
# A <- DMFA(test,num.fact=39,ncp=4)
# plot(A,axes=c(1,2),choix="var",label="all")
# dimdesc(A,axes=1:3)
# res.cov <- cov(scale(test[,-39]))

## ModÃ©lisation par groupe
##########################
require(sf)
require(data.table)
require(tidyverse)
require(ggplot2)
require(reshape2)
require(qgraph)
require(vegan)
require(raster)
library(dplyr)
options("scipen"=100, "digits"=4)


load("data/guild2/obs_guild2_sf_04-09-2019.Rdata")

div <- setDT(obs_sf)[,.(N=length(unique(taxonid))),by=.(grid.id,group)]
weigh.grid <- setDT(obs_sf)[,.(W=length(unique(group))),by=.(grid.id)]
div <- merge(div,weigh.grid,by="grid.id",all.x=T)

# ## For plants only with weigh
# div_in_cells <- obs_sf %>% filter(keep= %>% group_by(grid.id) %>% summarise(N=length(unique(taxonid)),W=sum(as.numeric(tapply(w,factor(taxonid),mean,na.rm=T))))
# head(div_in_cells)
# plot(div_in_cells[,c("N","W")])
#div_in_cells$Nrel <- div_in_cells$N/max(div_in_cells$N,na.rm=T)
# st_geometry(div_in_cells) <- NULL

# Subset environmental variables
load("data/guild2/river_buffer200m_with_enviro_sf.Rdata")
# buf2 <- buf[,c("grid.id", "BV04","sinu", "slope", "debit","geometry","bio15_ps","bio12_p", "bio11_tcoldq","bio6_tminc","bio4_ts", "gdd3Y_8110_ngb5_mwconic_", "topos","ndviMEAN","alt")]
buf2 <- buf[,c("grid.id", "BV04","sinu", "slope", "debit","geometry","bio12_p", "bio6_tminc","bio4_ts", "gdd3Y_8110_ngb5_mwconic_", "topos","ndviMEAN","alt")]
buf3 <- na.omit(buf2)
SUB <- left_join(buf3,as.data.frame(div),by="grid.id")
rm(list=c("buf","buf2"))

# Do not log-transform count data!! :
# shell.exec("https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2010.00021.x")
# shell.exec("C:/Users/rutishauser/Zotero/storage/CVH63YDQ/sdm.pdf")

## RandomForest
# If the response variable is a factor (categorical), randomForest will do classification, otherwise
# it will do regression.
GROUP <- levels(factor(obs_sf$group))
LR <- as.numeric(scale(table(factor(obs_sf$group)),center=F))/10  ## learning rates by group
STACK <- list()

system.time(for (k in 1:length(GROUP)) {
  
  training <- na.omit(SUB %>% filter(group==GROUP[k]))  %>% st_drop_geometry()
  test <- left_join(buf3,div[,c("grid.id","N","group")] %>% filter(group==GROUP[k]) ,by="grid.id")
  test$obs <- 0
  test <- within(test,obs[grid.id%in%training$grid.id] <- 1)
  fmla <- as.formula(paste("N ~ ", paste(names(training)[3:12], collapse= "+")))
  PRED <- test[,c("grid.id","BV04","N")]
  
  ## Random forest unweighted
  rf1 <- ranger::ranger(fmla, data=training)
  pr <- predict(rf1,test %>% st_drop_geometry(),type="response")
  pr$predictions[test$obs==0] <- pr$predictions[test$obs==0]/max(pr$predictions[test$obs==0],na.rm=T)
  PRED[,5] <- pr$predictions
  rm(list=c("rf1","pr"))
  
  # ## Random forest weighted
  # rf1w <- ranger::ranger(fmla, data=training,case.weights = sqrt(training$W))
  # prw <- predict(rf1w,test %>% st_drop_geometry(),type="response")
  # prw$predictions[test$obs==0] <- prw$predictions[test$obs==0]/max(prw$predictions[test$obs==0],na.rm=T)
  # PRED[,6] <- prw$predictions
  # rm(list=c("rf1w","prw"))
  # 
  # # Negative binomial model
  # mod <- MASS::glm.nb(fmla, data = training,link = log,offset(log(W)))
  # prw <- as.numeric(predict(mod,test %>% st_drop_geometry()))
  # prw[test$obs==0] <- prw[test$obs==0]/max(prw[test$obs==0],na.rm=T)
  # PRED[,7] <- prw
  # rm(list=c("mod","prw"))
  
  # ## Boosted regression tree (BRT)
  ## shell.exec("https://cran.r-project.org/web/packages/dismo/vignettes/brt.pdf")
  gbm <- dismo::gbm.step(training,gbm.x=names(training)[3:12],gbm.y="N",family="poisson",tree.complexity = ,learning.rate = LR[k], bag.fraction = 0.5,max.trees=1000,silent=F)
  prw <- predict(gbm,test %>% st_drop_geometry(),n.trees=500,,type="response")
  prw[test$obs==0] <- prw[test$obs==0]/max(prw[test$obs==0],na.rm=T)
  PRED[,6] <- prw
  rm(list=c("gbm","prw"))
  
  # GAM model
  fmlaG <- as.formula(paste("N ~ ", paste("s(",names(training)[3:12],")", collapse= "+")))
  testError <- tryCatch(GAM <- mgcv::gam(fmlaG,data=training,weights=sqrt(training$W),family=poisson)
                        , error=function(e) e)
  if(!inherits(testError, "error")) {
    prw <- as.numeric(predict(GAM,test %>% st_drop_geometry(),type="response"))
    prw[test$obs==0] <- prw[test$obs==0]/max(prw[test$obs==0],na.rm=T)
    PRED[,7] <- prw} else {
      PRED[,7] <- NA
    }
  rm(list=c("GAM","prw"))
  
  ## Take average value of all models in cells without observation
  STACK[[k]] <-  PRED
}) ## 45 min to run
save(STACK,file="data/guild2/stack_models_04-09-2019.Rdata")

##################
##### Analyse #### 
##################
load("data/guild2/stack_models_04-09-2019.Rdata")

STACK <- lapply(STACK,function(x) {
  names(x) <- c("grid.id", "BV", "N", "geometry", "rf", "rfw", "nbm", "gbm", "gam")
  return(x)
})


for (i in 1:11) {
  MELT <- melt(STACK[[i]] %>% st_drop_geometry(),id.vars=c("grid.id","BV","N"))
  G <- ggplot(MELT[!is.na(MELT$N),],aes(x=N,y=value)) + geom_point() + facet_wrap(~variable) + geom_abline(slope=1,intercept=0,col=2,linetype="dashed") + labs(x="Richesse spÃ©cifique observÃ©e",y="Richesse spÃ©cifique prÃ©dite (moyenne) ") + ggpubr::stat_cor(method = "spearman",aes(label =paste(..rr.label.., "' , p='",round(..p..,4),sep = "~"))) + geom_smooth(method="lm")
  ggsave(G,file=paste("pred_richness_",GROUP[i],".pdf"))
}

STACK <- lapply(STACK,function(x) {
  x$mean <- as.numeric(rowMeans(x[,c(5:6,8)] %>% st_drop_geometry()))
  return(x)
})

MEAN <- lapply(STACK,function(x) {
  return(x$mean)
})
# grid100 = raster("data/guild2/shp/altitude100.tif")
# RAST.STACK <-  stack(lapply(STACK,function(x) fasterize::fasterize(x,grid100,field="mean")))

mlist <- unlist(MEAN)

LE <- nrow(STACK[[1]])
rm(STACK)
save(temp,file="data/guild2/stack_models_df.Rdata")

# ggplot(temp[!is.na(temp$N),],aes(x=N,y=mean)) + geom_point() + facet_wrap(~group) + geom_abline(slope=1,intercept=0,col=2,linetype="dashed") + labs(x="Richesse spÃ©cifique observÃ©e",y="Richesse spÃ©cifique prÃ©dite (moyenne) ") + ggpubr::stat_cor(method = "spearman",aes(label =paste(..rr.label.., "' , p='",round(..p..,4),sep = "~"))) + geom_smooth(method="lm")

load("data/guild2/stack_models_df.Rdata")
GEO <- temp$geometry
temp2 <- temp[,c(1:3,5)]
rm(temp)
GROUP <- c("BRYO", "COLE", "EPHE", "FUNG", "LEPI", "LICH", "ODON", "ORTH", 
           "PLEC", "PLANT", "TRIC")
temp2$group <- rep(GROUP,each= 1045787)
temp2[is.na(N),mean(mean,na.rm=T),by=group]

IDX <- which(temp2$BV==50456)


ggplot(temp2[is.na(temp2$N),],aes(mean,colour=group)) + geom_density()
sub <- st_as_sf(temp2[IDX],geometry=GEO[IDX])
MED <- temp2[,median(mean),by=group]
sub$keep <- NA
sub <- within(sub,keep[mean>=MED$V1[match(group,MED$group)]]<-1 )

library(tmap)
tmap_mode("plot") 
G1 <- tm_shape(sub) +
  tm_fill(col="cyan") +
  tm_shape(sub %>% filter(is.na(N))) +
  tm_fill(
    col = "keep",
    palette = gray.colors(2),
    id = "grid.id",
    showNA = T,
    colorNA = "cyan",
    alpha = 1,
    legend.show = FALSE) +
  tm_facets(by="group",free.coords=F) +
  tm_borders(col = "gray", lwd = 0.2) 

tmap_save(G1,file="graphs/pred_50perc_bygroup.png")

tmap_mode("plot")
G2 <- tm_shape(sub) +
  tm_fill(col="cyan") +
  tm_borders(col = "gray", lwd = 0.2) +
  tm_facets(by="group",free.coords=F,ncol=6) +
  tm_shape(sub %>% filter(!is.na(N))) +
  tm_fill(
    col = "N",
    title="Artenreichtum",
    palette = "Reds",
    id = "grid.id",
    legend.show = F)  +
  tm_layout(legend.show = F) +
  tm_borders(col = "gray", lwd = 0.2) +
  tm_facets(by="group",free.coords=F,ncol=6) 
tmap_save(G2,file="graphs/obs_bygroup.png")

## Somme diversitÃ© cumulÃ©e
tt <- setDT(sub)[,.(sum(keep,na.rm=T),sum(N,na.rm=T)),by=grid.id]
names(tt) <- c("grid.id","div","obs")
tt <- within(tt,div[div==0] <- NA)
tt <- within(tt,obs[obs==0] <- NA)
sub <- left_join(sub,tt,by="grid.id",all.x=T)

RIV <- riv[riv$BV==50456,]

library(tmap)
tmap_mode("view") 

tm_shape(sub %>% filter(is.na(N))) +
  tm_fill(
    col = "div",
    title="PrÃ¤diktionen",
    palette = "Greens",
    breaks=seq(1,12,2),
    id = "grid.id",
    showNA = F,
    colorNA = "cyan",
    alpha = 1,
    legend.show = T) +
  tm_borders(col = "gray", lwd = 0.2) +
  tm_shape(sub %>% filter(!is.na(N))) +
  tm_fill(
    col = "obs",
    style="cont",
    palette = "Reds",
    title="Artenreichtum",
    id = "grid.id")  +
  tm_layout(legend.show = T) +
  tm_borders(col = "gray", lwd = 0.2) 



#############################################
#############################################
#############################################

mod.rich <- function(X) {
  training <- na.omit(X)  %>% st_drop_geometry()
  test$obs <- 0
  test <- within(test,obs[grid.id%in%training$grid.id] <- 1)
  
  fmla <- as.formula(paste("N ~ ", paste(names(training)[3:14], collapse= "+")))
  PRED <- test[,c("grid.id","BV","N")]
  
  ## Random forest unweighted
  rf1 <- ranger::ranger(fmla, data=training)
  pr <- predict(rf1,test %>% st_drop_geometry(),type="response")
  pr$predictions[test$obs==0] <- pr$predictions[test$obs==0]/max(pr$predictions[test$obs==0],na.rm=T)
  PRED[,4] <- pr$predictions
  rm(list=c("rf1","pr"))
  
  ## Random forest weighted
  rf1w <- ranger::ranger(fmla, data=training,case.weights = sqrt(training$W))
  prw <- predict(rf1w,test %>% st_drop_geometry(),type="response")
  prw$predictions[test$obs==0] <- prw$predictions[test$obs==0]/max(prw$predictions[test$obs==0],na.rm=T)
  PRED[,5] <- prw$predictions
  rm(list=c("rf1w","prw"))
  
  # Negative binomial model
  mod <- MASS::glm.nb(fmla, data = training,link = log,offset(log(W)))
  prw <- as.numeric(predict(mod,test %>% st_drop_geometry()))
  prw[test$obs==0] <- prw[test$obs==0]/max(prw[test$obs==0],na.rm=T)
  PRED[,6] <- prw
  rm(list=c("mod","prw"))
  
  ## Boosted regression tree (BRT)
  ## shell.exec("https://cran.r-project.org/web/packages/dismo/vignettes/brt.pdf")
  LR <- as.numeric(nrow(training))/10^6  ## learning rate
  gbm <- dismo::gbm.step(training,gbm.x=names(training)[3:14],gbm.y="N",family="poisson",tree.complexity = 5,
                         learning.rate = LR, bag.fraction = 0.3,max.trees=1000,silent=T)
  prw <- predict(gbm,test %>% st_drop_geometry(),n.trees=500,,type="response")
  prw[test$obs==0] <- prw[test$obs==0]/max(prw[test$obs==0],na.rm=T)
  PRED[,7] <- prw
  rm(list=c("gbm","prw"))
  
  # GAM model
  fmlaG <- as.formula(paste("N ~ ", paste("s(",names(training)[3:14],")", collapse= "+")))
  GAM <- mgcv::gam(fmlaG,data=training,family=poisson)
  prw <- as.numeric(predict(GAM,test %>% st_drop_geometry(),type="response"))
  prw[test$obs==0] <- prw[test$obs==0]/max(prw[test$obs==0],na.rm=T)
  PRED[,8] <- prw
  rm(list=c("GAM","prw"))
  
  ## Take average value of all models in cells without observation
  return(PRED)
}
dat <- split(BV,BV$group)
test <- na.omit(buf2)

library(parallel)
ncores<-detectCores()-1
if (ncores == 0) {ncores = 1}
cl <- makeCluster(ncores)
clusterExport(cl,varlist=c("dat","mod.rich","test"),envir=.GlobalEnv)
clusterEvalQ(cl,{
  library(sf)
  library(MASS)
  library(ranger)
  library(dismo)
  library(mgcv)
  library(dplyr)
})
system.time(test<- parLapply(cl,dat, function(x) mod.rich(x)))
stopCluster(cl)



G1 <- tm_shape(PRED %>% filter(BV==50456)) +
  tm_fill(
    col = "p",
    palette = "Reds",
    id = "grid.id",
    showNA = T,
    colorNA = "cyan",
    title = "Richness pred (W)",
    alpha = 1,
    legend.show = FALSE)  +
  tm_layout(legend.show = FALSE) +
  tm_borders(col = "gray", lwd = 0.2)

G2 <- tm_shape(PRED2 %>% filter(BV==50456)) +
  tm_fill(
    col = "p",
    palette = "Reds",
    id = "grid.id",
    showNA = T,
    colorNA = "cyan",
    title = "Richness pred (N)",
    alpha = 1,
    legend.show = FALSE)  +
  tm_layout(legend.show = FALSE) +
  tm_borders(col = "gray", lwd = 0.2)

sub <- BV %>% filter(BV==50456)
sub <- within(sub,W[W<1] <- NA)

G3 <- tm_shape(sub) +
  tm_fill(
    col = "W",
    palette = "Reds",
    id = "grid.id",
    title = "Richness obs (N)",
    showNA = T,
    colorNA = "cyan",
    alpha = 1,
    legend.show = FALSE)  +
  tm_layout(legend.show = FALSE) +
  tm_borders(col = "gray", lwd = 0.2)

sub <- within(sub,N[N<5] <- NA)
G4 <- tm_shape(sub) +
  tm_fill(
    col = "N",
    palette = "Reds",
    id = "grid.id",
    title = "Richness obs (N)",
    showNA = T,
    colorNA = "cyan",
    alpha = 1,
    legend.show = FALSE)  +
  tm_layout(legend.show = FALSE) +
  tm_borders(col = "gray", lwd = 0.2)

library(grid)
grid.newpage()
page.layout <- grid.layout(nrow = 2, ncol = 2, widths=c(1,1), heights=c(1,1))
pushViewport(viewport(layout = page.layout))
print(G3, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(G4, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(G1, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(G2, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))

ggplot(PRED[PRED$mod!="obs" & !is.na(PRED$N),],aes(x=N,y=p)) + geom_point() + facet_wrap(~mod) + geom_abline(slope=1,intercept=0,col=2,linetype="dashed") + labs(x="Richesse spÃ©cifique observÃ©e",y="Richesse spÃ©cifique prÃ©dite") + ggpubr::stat_cor(method = "spearman") + geom_smooth(method="lm")


library(tmap)
tmap_mode("plot") 
tm_shape(PRED %>% filter(BV==50456)) +
  tm_fill(
    col = "p",
    palette = "Reds",
    id = "grid.id",
    showNA = T,
    colorNA = "cyan",
    alpha = 1,
    legend.show = FALSE)  +
  tm_layout(legend.show = FALSE) +
  tm_borders(col = "gray", lwd = 0.2) +
  
  tm_shape(PRED %>% filter(BV==50456)) +
  tm_fill(
    col = "p",
    palette = "Reds",
    id = "grid.id",
    showNA = T,
    colorNA = "cyan",
    alpha = 1,
    legend.show = FALSE)  +
  tm_layout(legend.show = FALSE) +
  tm_borders(col = "gray", lwd = 0.2) +
  
  
  
  rf2 <- ranger::ranger(fmla, data=training,case.weights = sqrt(training$W))
pr1 <- predict(rf2,test,type="response")
P1<- sf:::cbind.sf(BV[,c("grid.id","N")],data.frame(p=pr1$predictions,mod="rfw"))
PRED <- rbind(PRED,P1)

## Boosted regression tree (BRT)
shell.exec("https://cran.r-project.org/web/packages/dismo/vignettes/brt.pdf")
gbm <- dismo::gbm.step(training,gbm.x=names(training)[3:14],gbm.y="W",family="poisson",tree.complexity = 5,
                       learning.rate = 0.5, bag.fraction = 0.5,max.trees=1000)
pr2 <- predict(gbm,BV,n.trees=500,,type="response")
P2 <- sf:::cbind.sf(BV[,c("grid.id","N")],data.frame(p=pr2,mod="brt"))
PRED <- rbind(PRED,P2)

# GAM model
fmla <- as.formula(paste("W ~ ", paste("s(",names(training)[3:13],")", collapse= "+")))
GAM <- mgcv::gam(fmla,data=training,family=poisson)
summary(GAM)
pr3 <- predict(GAM,BV,type="response")

P3 <- sf:::cbind.sf(BV[,c("grid.id","N")],data.frame(p=as.numeric(pr2),mod="gam"))
PRED <- rbind(PRED,P3)

# Negative binomial model
fmla <- as.formula(paste("W ~ ", paste(names(training)[3:13], collapse= "+")))
library(MASS)
mod <- glm.nb(fmla, data = training,link = log)
pr4 <- predict(mod,BV,type="response")
P4 <- sf:::cbind.sf(BV[,c("grid.id","N")],data.frame(p=as.numeric(pr4),mod="nbm"))
PRED <- rbind(PRED,P4)

# Rajoute les obs comme "model"
PRED <- merge(PRED,as.data.frame(riv[,c("grid.id","BV")]),by="grid.id")
P5 <- sf:::cbind.sf(BV[,c("grid.id","N")],data.frame(p=BV$N,mod="obs"))
PRED <- rbind(PRED,P5)

ggplot(PRED[PRED$mod!="obs" & !is.na(PRED$N),],aes(x=N,y=p)) + geom_point() + facet_wrap(~mod) + geom_abline(slope=1,intercept=0,col=2,linetype="dashed") + labs(x="Richesse spÃ©cifique observÃ©e",y="Richesse spÃ©cifique prÃ©dite") + ggpubr::stat_cor(method = "spearman") + geom_smooth(method="lm")

DENS <- data.frame(dens=PRED[PRED$mod=="rfw" & !is.na(PRED$N),]$p,group="obs")
DENS <- rbind(DENS,data.frame(dens=PRED[PRED$mod=="rfw" & is.na(PRED$N),]$p,group="no.obs") )
ggplot(DENS, aes(x = dens, fill = group)) + geom_density(alpha = 0.5) + scale_x_log10() + labs(x="Richesse spÃ©cifique") 

PRED$obs <- NA
PRED <- within(PRED,obs[!is.na(N)] <- "obs")
PRED <- within(PRED,obs[is.na(N)] <- "no.obs")

# B <- ggpubr::ggdensity(PRED, x = "p",
#           add = "mean", rug = TRUE,
#           color = "obs", fill = "obs",
#           palette = c("#00AFBB", "#E7B800"))
# facet(B, facet.by = "mod")

### Plot predictions
ID.cell = setDT(riv)[BV==c(50456),]$grid.id
sub <- PRED %>% filter(grid.id%in%ID.cell)
sub$mod <- factor(sub$mod,levels=c("obs","rf", "rfw", "brt","gam","nbm"))

riv3 <- st_read("C:/Dossier_Ervan/SIG/G2-Eaux.dynamiques/Typisierung_LV03/FGT.shp")
riv3 <- st_transform(riv3,st_crs(PRED))
riv3 <- st_join(riv3,st_as_sf(riv[,"BV"],geom=riv$geometry),join=st_intersects)
RIV <- riv3 %>% filter(BV==50456)
# BORD <- BV  %>% filter(BASIS_NR%in%ID.BV)

# plot residuals
library(tmap)
tmap_mode("plot") 
tm_shape(sub) +
  tm_fill(
    col = "p",
    palette = "Reds",
    id = "grid.id",
    showNA = T,
    colorNA = "cyan",
    alpha = 1,
    legend.show = FALSE)  +
  tm_layout(legend.show = FALSE) +
  tm_borders(col = "gray", lwd = 0.2) +
  tm_facets(by="mod",ncol=3,drop.units = TRUE, free.coords = TRUE)+ 
  tm_shape(RIV) +
  tm_lines(col="blue") +
  tm_legend(show=F)+
  tm_tiles(server="Esri.WorldTopoMap")

# plot predicted richness
tm_shape(RES %>% filter(BV%in%ID.BV)) +
  tm_fill(
    col = "rich.pred",
    palette = "RdPu",
    title = "vorhergesagt Artenreichtum (# Arten)",
    id = "grid.id",
    alpha =0.95)  +
  tm_borders(col = "gray", lwd = 0.2) +
  tm_shape(RIV) +
  tm_lines(col="blue")

# plot  richness
tm_shape(RES %>% filter(BV%in%ID.BV)) +
  tm_fill(
    col = "rich",
    palette = "RdPu",
    title = "Artenreichtum (# Arten)",
    id = "grid.id",
    showNA = T,
    colorNA = "gold2",
    alpha =0.95)  +
  tm_borders(col = "gray", lwd = 0.2) +
  tm_shape(RIV) +
  tm_lines(col="blue")









# ## A. GLM
# library('MuMIn')
# fmla <- as.formula(paste("N ~ ", paste(names(training)[3:13], collapse= "+")))
# GLM <- glm(fmla,data=training,family=poisson)
# summary(GLM)
# 
# options(na.action = "na.fail")
# ms1 <- MuMIn::dredge(GLM)  # model averaging
# head(ms1)
# 
# # best model
# BM <- get.models(ms1, 1)[[1]]
# summary(BM)
# A <- sw(ms1)
# par(mar=c(4,10,0,3),las=1)
# barplot(as.numeric(sw(ms1)),names.arg=names(sw(ms1)),horiz=T)
# sw(subset(ms1, delta <= 4))

# VAR.PRED <- names(coef(BM))[-1]
# PRED <- predict(BM,newdata=BV[,VAR.PRED],type="response")
# plot(PRED,BV2$N)

# Negative binomial model
fmla <- as.formula(paste("N ~ ", paste(names(training)[1:9], collapse= "+")))
library(MASS)
mod <- glm.nb(fmla, data = training,link = log,offset(log(W)))
PRED <- predict(mod,BV2)
