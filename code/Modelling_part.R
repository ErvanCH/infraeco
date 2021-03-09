##########################################################
##########################################################
####   MODELING PART
##########################################################
##########################################################
DATE <- format(Sys.time(), '%d-%m-%y')
Time1 <- Sys.time()
pacman::p_load(parallel, doParallel , ranger , xgboost , caret , data.table, mgcv , dismo , rJava  , Hmisc , dplyr , maxnet, lightgbm, methods)  
source("code/functions_ER.R")

if (!dir.exists(paste0(dir,"/model"))) {
  dir.create(paste0(dir,"/model"), showWarnings = FALSE)
}
dest.folder <- paste0(dir,"/model")


obs.mod<-OBS[,..VAR]

##### scaling
PAR <- list() ## object to store parameters
PAR$scale <- attributes(scale(obs.mod))[3:4]
# Make new dataset with scaled predictors
obs.mod <- as.data.table(scale(obs.mod))
obs.mod <- mutate(obs.mod,Qobs=OBS$Qobs, grid.id=OBS$grid.id)
k<-5 # NFold CV : five to ten is commonly accepted !
obs.mod <- mutate(obs.mod,
                  my.folds = sample(1:k,
                                    size = nrow(obs.mod),
                                    replace = TRUE))

### Check if parameters were estimated
if (any(grep("model_parameter",list.files(dest.folder)))) {
  rep <- menu(c("Yes", "No"), title="Model parameters are existing. Do you want to estimate new ones?")
} 
if (rep==1) {
  
#############################################################
##### MODELING tuning (without carret)
#############################################################
ncores<-detectCores()-1
cl <- makeCluster(ncores,outfiles = 'dest.folder' )


#Maxent parameters to be tuned
my.regmul<-c(0.3,0.6,1,2,4)
me.classes= c('lqh')

#GAM parameters to be tuned
my.k<-c('k=6','k=9','k=12','k=18')

#CV used for the evaluation
cv<-2


#Generate parameters' grid
tune.table<- mod.table(GAM=guild.info$GAM, ME=guild.info$ME,GBM=F,GBM_xg=F,cv=cv,gam.smoothing=my.k,ME.method=guild.info$ME_pack,ME.regmul=my.regmul,ME.classes=me.classes)

clusterExport(cl,varlist=c('obs.mod',
                           'tune.table',
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
  library(lightgbm)
})

system.time(model_tuning <- parLapply(cl,1:nrow(tune.table),eval.mod,mod_table=tune.table,data=obs.mod,my.var=VAR,GBM.param=PAR$GBM_para))#81 sec
tuning.synth<-as.data.table(t(sapply(model_tuning,synt.eval)))
colnames(tuning.synth)[ncol(tuning.synth)]<-'mean_evaluation'

if (guild.info$GAM==T) {
PAR$GAM_k<-tuning.synth[model=='GAM'& mean_evaluation==max(unlist(tuning.synth[model=='GAM']$mean_evaluation))]$GAM_k[[1]]
}
if (guild.info$ME==T) {
PAR$ME_regmul<-as.numeric(tuning.synth[model=='ME'& mean_evaluation==max(unlist(tuning.synth[model=='ME']$mean_evaluation))]$ME_regmul[[1]])
PAR$ME_classes<-tuning.synth[model=='ME'& mean_evaluation==max(unlist(tuning.synth[model=='ME']$mean_evaluation))]$ME_classes[[1]]
} 
cat(paste0("Part1: GAM/ME done in ",round(difftime(Sys.time(),Time1,units="mins"),2)," min."), "\n")

#############################################################
##### MODELING GBM tuning
#############################################################
## Select best parameters for the GBM (time & RAM consuming!!!)
## From this post: 
# https://www.kaggle.com/theoverfitter3139/lightgbm-with-gridsearch-starter-1
if (guild.info$GBM){
stopImplicitCluster()
inTrain <- createDataPartition(y = obs.mod$grid.id, p = 0.80, list = FALSE)  # 80% des données
dtrain <- lgb.Dataset(as.matrix(obs.mod[inTrain,VAR]),label=obs.mod[inTrain,]$Qobs)
dtest <- lgb.Dataset(rep(as.matrix( obs.mod[-inTrain,VAR]),10),label= rep(obs.mod[-inTrain,]$Qobs,10))

grid_search <- expand.grid(
  num_leaves        = c(3,5,7),
  max_depth         = c(4,6,8),
  subsample         = c(0.7,0.9),
  colsample_bytree  = c(0.7,0.9),
  min_child_weight  = c(0,0.01),
  scale_pos_weight  = c(100,300)
)

model <- list()
perf <- numeric(nrow(grid_search))
cat("GBM run")
for (i in 1:nrow(grid_search)) {
  if(i%in%seq(5,nrow(grid_search),10)) {cat(paste0(i,"-"))}
  model[[i]] <- lgb.train(
    list(objective         = "binary",
         metric            = "auc",
         learning_rate     = 0.1,
         min_child_samples = 100,
         max_bin           = 100,
         subsample_freq    = 1,
         num_leaves        = grid_search[i, "num_leaves"],
         max_depth         = grid_search[i, "max_depth"],
         subsample         = grid_search[i, "subsample"],
         colsample_bytree  = grid_search[i, "colsample_bytree"],
         min_child_weight  = grid_search[i, "min_child_weight"],
         scale_pos_weight  = grid_search[i, "scale_pos_weight"]),
    dtrain,
    valids = list(validation = dtest),
    nthread = 4, 
    nrounds = 5, # increase/ decrease rounds
    early_stopping_rounds = 2,
    verbose=-1  #avoid output
  )
  perf[i] <- max(unlist(model[[i]]$record_evals[["validation"]][["auc"]][["eval"]]))
  invisible(gc()) # free up memory after each model run
}


PAR$GBM_para = grid_search[which.max(perf), ]
} else if (guild.info$GBM_xg){
  ## XGBoost version
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
  
  system.time(xgb_model <- train(training[,VAR], as.factor(training$Qobs), trControl = xgb_trcontrol, tuneGrid = xgbGrid, method = "xgbTree",objectiv="binary:logistic",eval_metric='auc'))  # 18 min
  PAR$GBM_para<- xgb_model$bestTune
}  # end of GBM loop 

cat(paste0("Part2: GBM done in ",round(difftime(Sys.time(),Time1,units="mins"),2)," min."), "\n")
save(PAR,file=paste0(dest.folder,"/model_parameter.",DATE,".Rdata"))


#############################################################
##### MODELING ENSEMBLE evaluation
#############################################################
# load(paste0(dest.folder,"/model_parameter.",DATE,".Rdata"))
cv<-1:5

my_table <- mod.table(GAM=guild.info$GAM,ME=guild.info$ME,GBM=guild.info$GBM,GBM_xg=guild.info$GBM_xg,cv=cv,gam.smoothing=PAR$GAM_k,ME.method=guild.info$ME_pack,ME.regmul=PAR$ME_regmul,ME.classes=PAR$ME_classes)

clusterExport(cl,varlist=c('PAR',
                           'my_table'),
              envir=.GlobalEnv)

#Ensemble Evaluation
system.time(model_eval<-parLapply(cl,1:nrow(my_table),eval.mod,mod_table=my_table,data=obs.mod,my.var=VAR,GBM.param=PAR$GBM_para)) #63 sec

eval.synth<-as.data.table(t(sapply(model_eval,synt.eval)))
colnames(eval.synth)[ncol(eval.synth)]<-'mean_evaluation'

system.time(MyModelEvaluation<-ensemble.eval(obs.mod,eval.synth,model_eval,eval.lim = guild.info$eval.lim)) # 12 secondes
MyModelEvaluation$model_weight
save(MyModelEvaluation,file = paste0(dest.folder,'/',GUILD,"_eval_mod_",DATE,".Rdata"))

}
#############################################################
##### MODELING ENSEMBLE projection
#############################################################
MyModelEvaluation <- LOAD("eval_mod","model",dir)
PAR <- LOAD("model_parameter","model",dir)

selected.models<-MyModelEvaluation[[1]]$model[which(MyModelEvaluation$model_weight>guild.info$eval.lim)]
env.pot.mod<-env.pot[,..VAR]
env.pot.mod<-as.data.table(scale(env.pot.mod,center= PAR$scale$`scaled:center`[names(env.pot.mod)], 
                                 scale=PAR$scale$`scaled:scale`[names(env.pot.mod)]))
env.pot.mod<-mutate(env.pot.mod,env.pot$Qobs)
add.para<-list(GBM.param = PAR$GBM_para, me.method = guild.info$ME_pack, me.regmul = PAR$ME_regmul,my.k=PAR$GAM_k,me.classes=PAR$ME_classes)

### individual model version useful if you want to run only one specific model
for (m in selected.models){  
  ifelse(nrow(env.pot.mod)<8*10^5,Npart<-list('RF'= 1, 'GBM' = 15,'GAM' =10,'MAXENT'=5 ),Npart<-list('RF'= 1, 'GBM' = 20,'GAM' =20,'MAXENT'=10 ))
  M <- model.proj(selected.model=m,cal.data=as.data.frame(obs.mod),proj.data = env.pot.mod, my.var= VAR,add.para = add.para,save.dir=dest.folder,Npart=Npart)
}

cat(paste0("Part3: projections done in ",round(difftime(Sys.time(),Time1,units="mins"),2)," min."), "\n")

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
VarIm <- data.table::data.table("var"=VAR,"mean"=apply(var.import,1,mean),"sd"=apply(var.import,1,sd),"wmean"=apply(var.import,1,weighted.mean,w=MyModelEvaluation[[3]][which(MyModelEvaluation$model_weight>guild.info$eval.lim)]))
save(VarIm,file=paste0(dest.folder,"/",GUILD,"-importance_variable.Rdata"))

## Find best threshold for predictions
PRED[,"mean" := rowMeans(.SD), .SDcols = selected.models]
PRED$Wmean <-apply(PRED[,..selected.models],1,weighted.mean,w=MyModelEvaluation[[3]][which(MyModelEvaluation$model_weight>guild.info$eval.lim)])

library(PresenceAbsence)
DD <- PRED[PRED$grid.id%in%OBS$grid.id,c("grid.id","cal","mean","Wmean")]
OPTI <- PresenceAbsence::optimal.thresholds(DD[!is.na(DD$cal),])
M <- mean(OPTI$Wmean[c(2,3,4,9)])

## Set predictions
TT <- PRED[!is.na(PRED$cal),c("mean","cal")]
TT$Qpred <- 0
TT <- within(TT,Qpred[mean>=M ] <- 1)
cat(caret::confusionMatrix(as.factor(TT$Qpred),as.factor(TT$cal))) 


# ### Add spatial information and define potential habitat
load("C:/Dossier_Ervan/R/Grid100/grid100_sf.Rdata")
PRED <- merge(PRED,grid_sf[,c("grid.id","CNHA","canton","subreg","geostat","forest","urban","lac","geometry","centro")],by="grid.id",all.x=TRUE)
PRED <- PRED[!duplicated(PRED$grid.id) & lac==0,]
save(PRED,file=paste0(getwd(),"/data/",GUILD,"_quality_predicted_",DATE,".Rdata"))
cat(paste0("Modeling of guilde ",GUILD," done in ",round(difftime(Sys.time(),Time1),2)," sec."), "\n")

    