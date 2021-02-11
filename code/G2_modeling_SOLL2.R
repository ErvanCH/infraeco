options("scipen"=100, "digits"=4)
require(sf)
require(data.table)
require(tidyverse)
require(ggplot2)
require(reshape2)
require(qgraph)
require(vegan)
require(raster)
source("C://Dossier_Ervan/R/functions_ER.R")

### UPDATE OUTPUT
GUILD <- 2
guild <- readxl::read_xlsx("C://Dossier_Ervan/Guildes&co/Guildes_revues_ER.xlsx",sheet=1,col_names = T)
guild.name <- guild[which(guild$No==GUILD),]$name_guild_de
DATE <- format(Sys.time(), '%d-%m-%y')
dir <- substr(getwd(),1,max(stringr::str_locate_all(getwd(),'/')[[1]]))

## IST
load("C:/Dossier_Ervan/R/Grid100/grid100_sf_with_enviro.Rdata")
dir <- substr(getwd(),1,max(stringr::str_locate_all(getwd(),'/')[[1]]))
FF <- list.files(paste0(dir,"INFOFAUNA"))[grep(paste0("_",GUILD,"_bio-idx"),list.files(paste0(dir,"INFOFAUNA")))]
ist <- read.csv(paste(dir,"INFOFAUNA", FF,sep="/"))
IST <- merge(setDT(ist),setDT(grid_sf)[,c("CNHA","grid.id","geometry")],by="CNHA",all.x=T)
st_geometry(IST) <- "geometry"


## Subset obs in river buffer
buf <- LOAD("data/buffer_river_200m_filtered_BV_polygons_sf.Rdata")
st_write(buf,"D:/SIG/SIG_BAFU/Produits_EI/2_eaux_dyn_EG.shp")
grid <- raster::raster("D://SIG/data/grid100/grid100.tif")
Q_ra <- fasterize::fasterize(buf,grid)
writeRaster(Q_ra,"D:/SIG/SIG_BAFU/Produits_EI/2_eaux_dyn_EG.tif",format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)
  


obs <- st_transform(IST,st_crs(buf))
system.time(test <- st_join(obs, buf[,c("in.riv")], join = st_intersects)) # 22 sec
table(is.na(test$in.riv)) # 3063 (8%) obs falling outside 
IST <- test[!is.na(test$in.riv) & !duplicated(test$grid.id),]

### D?sagr?gation des donn?es
system.time(IST2 <- IST.dsg(IST,300))

# 2. Compute quality per ha
# Assign spatial info & enviro variables
TT <- readxl::read_xlsx("C://Dossier_Ervan/Guildes&co/var_enviro_by_guild.xlsx")
VAR <- TT$label[!is.na(TT[,stringr::str_which(names(TT),paste0("G",GUILD,"$"))])]

if(!exists("grid_sf")) {
  load("C:/Dossier_Ervan/R/Grid100/grid100_sf_with_enviro.Rdata")
}
div2 <- setDT(grid_sf)[grid.id%in%IST2$grid.id,c("grid.id","BV04","canton","subreg","geometry","centro",..VAR)]
summary(div2)
sf::st_geometry(div2) <- "centro"

## Data for calibration
Q1 <- na.omit(data.table::setDT(div2)) 
nrow(IST2)-nrow(Q1) # 2628 ha vires
Q1$Qobs=1

# Pseudo-absence
EG <-setDT(grid_sf)[grid.id%in%IST$grid.id | grid.id%in%buf$grid.id,c("grid.id","BV04","canton","subreg","geometry","centro",..VAR)]
EG2 <- EG[!duplicated(EG$grid.id),]
NArep(EG2)
summary(EG2)

Q0 <- dplyr::sample_n(na.omit(data.table::setDT(EG2)[!grid.id%in%Q1$grid.id]),nrow(Q1),replace=F)
Q0$Qobs=0
OBS <- rbind(Q1,Q0)

## Selection des variables
library(rgdal)
library(rgeos)
library(parallel)
library(raster)
library(perm)

ncores<-detectCores()-1
if (ncores == 0) {ncores = 1}

env.pres<-as.data.frame(OBS[OBS$Qobs==1,..VAR])#env data for presences
env.bck<-as.data.frame(OBS[OBS$Qobs==0,..VAR])#env data for background
to.do<-1:ncol(env.pres) #index for parallelization

cl <- makeCluster(ncores)
clusterExport(cl,varlist=c('to.do','env.bck','env.pres','var.test'),envir=.GlobalEnv)
clusterEvalQ(cl,{
  library(perm)
})

#apply non-parametric t-test between background and presences for each variable
system.time(VAR.test<-parLapply(cl,X=to.do,fun=var.test,Nrep=5,mcmc=1000,env.pres,env.bck))  # 10 sec
stopCluster(cl)

#results formatting
select.var <- data.frame(var=VAR,do.call(rbind,VAR.test))
write.table(select.var, file = "data/selection_var.txt",sep='\t',quote=F,row.names = F)

## Generate report for variable selection
rmarkdown::render(input = paste0(getwd(),"/code/var_select.Rmd"),
                  output_format = "html_document",
                  output_file = paste0(GUILD,"_var_select2.html"),
                  output_dir = paste0(getwd(),"/report"))


## cREATION DE L'ESPACE GUILDE
TT <- readxl::read_xlsx("C://Dossier_Ervan/Guildes&co/var_enviro_selected.xlsx")
VAR <- TT$label[!is.na(TT[,stringr::str_which(names(TT),paste0("G",GUILD,"$"))])]

## Projection dans l'espace-guilde (to lower computational time)
summary(EG2[,..VAR])
env.pot <- na.omit(setDT(EG2)[!duplicated(grid.id),c("grid.id","BV04","subreg",..VAR)]) ## 45754 ha
nrow(EG2) - nrow(env.pot)  # 355 ha perdus sur la fronti?re

# st_geometry(EG2) <- "geometry"
# mapview::mapview(EG2[!EG2$grid.id%in%env.pot$grid.id,],legend=F)

## Projection dans l'espace-guilde (to lower computational time)
env.pot <- merge(env.pot,setDT(IST)[,c("grid.id","BIOIDX_TXG")],by="grid.id",all.x=T)
env.pot$ist <- env.pot$BIOIDX_TXG
env.pot$BIOIDX_TXG <- NULL

env.pot <- merge(env.pot,setDT(OBS)[,c("grid.id","Qobs")],by="grid.id",all.x=T)
env.pot$cal <- env.pot$Qobs
env.pot$Qobs <- NULL
table(duplicated(env.pot$grid.id))
env.pot <- env.pot[!duplicated(env.pot$grid.id),]
save(env.pot,file=paste0("data/",GUILD,"_envpot.RData"))

rm(list=setdiff(ls(), c("env.pot","OBS","VAR","GUILD","guild.name","dir")))
save.image(paste0("data/",GUILD,"_enviro4model.RData"))
.rs.restartR()

##########################################################
##########################################################
####   MODELING PART
##########################################################
##########################################################
GUILD <- 2
# load(paste0("data/",GUILD,"_enviro4model.RData"))
DATE <- format(Sys.time(), '%d-%m-%y')
TT <- readxl::read_xlsx("C://Dossier_Ervan/Guildes&co/var_enviro_selected.xlsx")
VAR <- TT$label[!is.na(TT[,stringr::str_which(names(TT),paste0("G",GUILD,"$"))])]

require(parallel)
require(doParallel)
require(ranger)
require(xgboost)
require(caret)
require(data.table)
require(mgcv)
require(dismo) ### Don't forget to install the program maxent.jar in the directory of the dismo library : dismo/java/maxent.jar This package is not necessary if you use maxnet (recommended)
Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_251')
require(rJava) ### This package is not necessary if you use maxnet (recommended)
require(Hmisc)
require(dplyr)
require(maxnet)
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m")) ### allocate more memory for maxent. This is not necessary if you use maxnet (recommended)
source("C://Dossier_Ervan/R/functions_ER.R")


if (!dir.exists(paste0(getwd(),"/model"))) {
  dir.create(paste0(getwd(),"/model"), showWarnings = FALSE)
}
dest.folder <- paste0(getwd(),"/model")

ncores<-detectCores()-1
cl <- makeCluster(ncores,outfiles = 'dest.folder' )

##### scaling
obs.mod<-OBS[,..VAR]

# I recommend to always store the scaling coefficients:
scale_attrib <- attributes(scale(obs.mod))[3:4]

# Make new dataset with scaled predictors
obs.mod <- as.data.table(scale(obs.mod))

obs.mod<-mutate(obs.mod,Qobs=OBS$Qobs, grid.id=OBS$grid.id)

PAR <- list()
PAR$scale <- scale_attrib
#############################################################
##### MODELING tuning (without carret)
#############################################################
k<-5 # NFold CV : five to ten is commonly accepted !
obs.mod <- mutate(obs.mod,
                  my.folds = sample(1:k,
                                    size = nrow(obs.mod),
                                    replace = TRUE))

#Maxent parameters to be tuned
my.regmul<-c(0.3,0.6,1,2,4)
me.classes= c('lqh')

#GAM parameters to be tuned
my.k<-c('k=6','k=9','k=12','k=18')

#CV used for the evaluation
cv<-2


#Generate parameters' grid
tune.table<-mod.table(GAM=F, ME=T,cv=cv,gam.smoothing=my.k,ME.method='maxnet',
                      ME.regmul=my.regmul,ME.classes=me.classes)

clusterExport(cl,varlist=c('obs.mod',
                           'tune.table',
                           'cv.fun',
                           'VAR',
                           'KappaRepet', 
                           'model.eval',
                           'boycei',
                           'ecospat.boyce',
                           'KappaStat',
                           'TSS.Stat',
                           'eval.mod'),
              envir=.GlobalEnv)

clusterEvalQ(cl,{
  library(ranger)
  library(xgboost)
  library(caret)
  library(data.table)
  library(mgcv)
  library(dismo) # only if you use maxent with dismo
  library(rJava)# only if you use maxent with dismo
  library(Hmisc)
  library(dplyr)
  library(maxnet)
})

system.time(model_tuning<-parLapply(cl,1:nrow(tune.table),eval.mod,mod_table=tune.table,data=obs.mod,my.var=VAR))#81 sec
tuning.synth<-as.data.table(t(sapply(model_tuning,synt.eval)))
colnames(tuning.synth)[ncol(tuning.synth)]<-'mean_evaluation'

if (any(grep("GAM",tune.table[,2]))){
  PAR$GAM_k<-tuning.synth[model=='GAM'& mean_evaluation==max(unlist(tuning.synth[model=='GAM']$mean_evaluation))]$GAM_k[[1]]
}
PAR$ME_regmul<-as.numeric(tuning.synth[model=='ME'& mean_evaluation==max(unlist(tuning.synth[model=='ME']$mean_evaluation))]$ME_regmul[[1]])
PAR$ME_classes<-tuning.synth[model=='ME'& mean_evaluation==max(unlist(tuning.synth[model=='ME']$mean_evaluation))]$ME_classes[[1]]


#############################################################
##### MODELING GBM tuning
#############################################################
## Select best parameters for the GBM (time & RAM consuming!!!)
doParallel::registerDoParallel(cl) ### This is to turn on parallel computing in Caret
xgb_trcontrol = trainControl(method = "cv", number = 5, allowParallel = TRUE, verboseIter = FALSE, returnData = FALSE)
inTrain <- caret::createDataPartition(y = obs.mod$grid.id, p = 0.8, list = FALSE)  # 80% des données
training <- obs.mod[inTrain,]
testing <- obs.mod[-inTrain,]


xgbGrid <- expand.grid(nrounds = c(200,500),
                       max_depth = c(3,5,7),
                       colsample_bytree = seq(0.5, 0.9, length.out = 3),
                       eta = c(0.02,0.1,0.2),
                       gamma=c(0.2, 0.5,1),
                       min_child_weight = 1, #valeur par d?faut
                       subsample = 1 #valeur par d?faut
) ### 162 parameters combinations

system.time(xgb_model <- train(training[,..VAR], as.factor(training$Qobs), trControl = xgb_trcontrol, tuneGrid = xgbGrid, method = "xgbTree",objectiv="binary:logistic",eval_metric='auc'))  # 11 min
PAR$GBM.param<- xgb_model$bestTune

save(PAR,file=paste0(dest.folder,'/model_parameter.Rdata'))

#############################################################
##### MODELING ENSEMBLE evaluation
#############################################################
load(paste0(dest.folder,'/model_parameter.Rdata'))
cv<-1:5

my_table<-mod.table(cv=cv,gam.smoothing=PAR$GAM_k,ME.method='dismo',GBM.para = nrow(PAR$GBM.param),ME.regmul=PAR$ME_regmul,ME.classes=PAR$ME_classes)

clusterExport(cl,varlist=c('PAR',
                           'my_table'),
              envir=.GlobalEnv)

#Ensemble Evaluation
system.time(model_eval<-parLapply(cl,1:nrow(my_table),eval.mod,mod_table=my_table,data=obs.mod,my.var=VAR,GBM.param=PAR$GBM.param)) #163 sec
eval.synth<-as.data.table(t(sapply(model_eval,synt.eval)))
colnames(eval.synth)[ncol(eval.synth)]<-'mean_evaluation'

eval.lim<-0.5

system.time(MyModelEvaluation<-ensemble.eval(obs.mod,eval.synth,model_eval,eval.lim = eval.lim)) # 12 secondes
MyModelEvaluation$model_weight
save(MyModelEvaluation,file = paste0(dest.folder,'/',GUILD,"_eval_mod.Rdata"))


#############################################################
##### MODELING ENSEMBLE projection
#############################################################
MyModelEvaluation <- LOAD("eval_mod","model")
load(paste0(dest.folder,'/model_parameter.Rdata'))
eval.lim<-0.5

selected.models<-MyModelEvaluation[[1]]$model[which(MyModelEvaluation$model_weight>eval.lim)]
env.pot.mod<-env.pot[,..VAR]
env.pot.mod<-as.data.table(scale(env.pot.mod,center= PAR$scale$`scaled:center`[names(env.pot.mod)], 
                                 scale=PAR$scale$`scaled:scale`[names(env.pot.mod)]))
env.pot.mod<-mutate(env.pot.mod,env.pot$Qobs)
add.para<-list(GBM.param = PAR$GBM.param, me.method = 'maxnet', me.regmul = PAR$ME_regmul,my.k=PAR$GAM_k,me.classes=PAR$ME_classes)

### individual model version useful if you want to run only one specific model
for (m in selected.models){  
  M <- model.proj(selected.model=m,cal.data=as.data.frame(obs.mod),proj.data = env.pot.mod, my.var= VAR,add.para = add.para,save.dir=dest.folder,Npart=list('RF'= 1, 'GBM' = 10,'GAM' =10,'MAXENT'=10 ))
}

### Average model
PRED <- env.pot[,c("grid.id","BV04","ist","cal")] ## Define receiving df 

var.import<-c()
mod.pred<-c()

for ( i in selected.models){
  load(paste0(dest.folder,'/mod_',i,'.Rdata'))
  mod.pred<-cbind(mod.pred,full.model$guild.pred)
  var.import<-cbind(var.import,as.vector(full.model$varImp))
}
colnames(mod.pred)<-selected.models
colnames(var.import)<-selected.models
row.names(var.import)<-VAR
PRED<-cbind(PRED,mod.pred)

#Draw variable contributions
VarIm <- data.table::data.table("var"=VAR,"mean"=apply(var.import,1,mean),"sd"=apply(var.import,1,sd),"wmean"=apply(var.import,1,weighted.mean,w=MyModelEvaluation[[3]][which(MyModelEvaluation$model_weight>eval.lim)]))
save(VarIm,file=paste0("model/",GUILD,"-importance_variable.Rdata"))

## Find best threshold for predictions
PRED[,"mean" := rowMeans(.SD), .SDcols = selected.models]
PRED$Wmean <-apply(PRED[,..selected.models],1,weighted.mean,w=MyModelEvaluation[[3]][which(MyModelEvaluation$model_weight>eval.lim)])

library(PresenceAbsence)
DD <- PRED[PRED$grid.id%in%OBS$grid.id,c("grid.id","cal","mean","Wmean")]
OPTI <- PresenceAbsence::optimal.thresholds(DD[!is.na(DD$cal),])
M <- mean(OPTI$Wmean[c(2,3,4,9)])

## Set predictions
PRED$Qpred <- 0
PRED <- within(PRED,Qpred[mean>=M]<-1)
caret::confusionMatrix(as.factor(PRED[!is.na(PRED$cal),]$Qpred),as.factor(PRED[!is.na(PRED$cal),]$cal))


# ### Add spatial information and define potential habitat
load("C:/Dossier_Ervan/R/Grid100/grid100_sf.Rdata")
PRED <- merge(PRED,grid_sf[,c("grid.id","CNHA","canton","subreg","geostat","forest","urban","lac","geometry","centro")],by="grid.id",all.x=TRUE)
PRED <- PRED[!duplicated(PRED$grid.id),]
save(PRED,file=paste0(getwd(),"/data/",GUILD,"_quality_predicted.Rdata"))

rm(list=setdiff(ls(),c("VAR","GUILD","DATE")))
.rs.restartR()

################################################
# Comparaison par cluster de bassins versants
################################################
source("C://Dossier_Ervan/R/functions_ER.R",echo=F)
library(sf)
library("tidyverse")
library(data.table)
DATE <- format(Sys.time(), '%d-%m-%y')
GUILD <- 2
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



# Create leaflets for plausibilisation
library(leafgl)
library(leaflet)
library(leafem)
library(leafpop)
library(data.table)
library(raster)
library(tidyverse)
library(sf)
GUILD=2
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
  addLegend(colors=V1[[2]]$col,labels=V1[[2]]$lab,group = "Erganzungsbedarf", position = "topright",opacity=0.8,title=paste0("Ergänzungsbedarf [ha] (max:",max(V1[[1]]$Erganzungsbedarf,na.rm=T),")"))


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
  addLayersControl(baseGroups= c("Positron","Esri"),overlayGroups = c("Erganzungsbedarf","Observed qual","Potential qual"),position="topleft") %>% 
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
  addLegend(colors=V1[[2]]$col,labels=V1[[2]]$lab,group = "EB", position = "topright",opacity=0.8,title=paste0("Ergänzungsbedarf [ha] (max:",max(V1[[1]]$Erganzungsbedarf,na.rm=T),")"))


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

