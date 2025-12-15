# GLMs for Species Distribution Modelling
# Import dataset from *.xlsx as "data"

library(biomod2)
library(rgdal)
library(raster)
library(rasterVis)
library(ggplot2)
library(gridExtra)

setwd('K:/Prove_SDM')

#getting presence values from data
myRespName <- 'Toumeyella parvicornis'
myResp <- as.numeric(data[,"status"])
#getting presence coordinates fromd data
myRespXY <- data[,c("x", "y")]

#environmental layers of Italy as explanatory variables (from Windows)
bio1 <- raster("K:/Lavoro/WorldClim/Italy (30arcsecond)/BIO01.tif")
bio2 <- raster("K:/Lavoro/WorldClim/Italy (30arcsecond)/BIO02.tif")
bio3 <- raster("K:/Lavoro/WorldClim/Italy (30arcsecond)/BIO03.tif")
bio4 <- raster("K:/Lavoro/WorldClim/Italy (30arcsecond)/BIO04.tif")
bio8 <- raster("K:/Lavoro/WorldClim/Italy (30arcsecond)/BIO08.tif")
bio9 <- raster("K:/Lavoro/WorldClim/Italy (30arcsecond)/BIO09.tif")
bio12 <- raster("K:/Lavoro/WorldClim/Italy (30arcsecond)/BIO12.tif")
bio15 <- raster("K:/Lavoro/WorldClim/Italy (30arcsecond)/BIO15.tif")
bio19 <- raster("K:/Lavoro/WorldClim/Italy (30arcsecond)/BIO19.tif")
current_bios = stack(bio1, bio2, bio3, bio4, bio8, bio9, bio12, bio15, bio19)

#environmental layers of Europe as explanatory variables (from Windows)
bio1 <- raster("K:/Lavoro/WorldClim/Europe (2.5minutes)/BIO01.tif")
bio2 <- raster("K:/Lavoro/WorldClim/Europe (2.5minutes)/BIO02.tif")
bio3 <- raster("K:/Lavoro/WorldClim/Europe (2.5minutes)/BIO03.tif")
bio4 <- raster("K:/Lavoro/WorldClim/Europe (2.5minutes)/BIO04.tif")
bio8 <- raster("K:/Lavoro/WorldClim/Europe (2.5minutes)/BIO08.tif")
bio9 <- raster("K:/Lavoro/WorldClim/Europe (2.5minutes)/BIO09.tif")
bio12 <- raster("K:/Lavoro/WorldClim/Europe (2.5minutes)/BIO12.tif")
bio15 <- raster("K:/Lavoro/WorldClim/Europe (2.5minutes)/BIO15.tif")
bio19 <- raster("K:/Lavoro/WorldClim/Europe (2.5minutes)/BIO19.tif")
current_bios = stack(bio1, bio2, bio3, bio4, bio8, bio9, bio12, bio15, bio19)

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = current_bios,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 3,
                                     PA.nb.absences = 1000,
                                     PA.strategy = 'random',
                                     na.rm = TRUE)

myBiomodData
plot(myBiomodData)

# Modeling
myBiomodOption <- BIOMOD_ModelingOptions(MAXENT = list(
  path_to_maxent.jar = "K:/Lavoro/MAXENT/"))
myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('MARS'),
  models.options = myBiomodOption,
  NbRunEval=10,
  DataSplit=80,
  Prevalence=0.5,
  VarImport=3,
  models.eval.meth = c('TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = F,
  modeling.id = 'ex2')
myBiomodModelOut

# Evaluation scores from all models
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
myBiomodModelEval

dim(myBiomodModelEval)
dimnames(myBiomodModelEval)

# Plotting models evaluation scores
models_scores_graph(
  myBiomodModelOut,
  by='models',
  metrics=c('ROC','TSS'),
  xlim=c(0.5,1),
  ylim=c(0.5,1)
)

# TSS scores
myBiomodModelEval["TSS","Testing.data",,,]
mean(myBiomodModelEval["TSS","Testing.data",,,])
sd(myBiomodModelEval["TSS","Testing.data",,,])

# ROC scores of all selected models
myBiomodModelEval["ROC","Testing.data",,,]
mean(myBiomodModelEval["ROC","Testing.data",,,])
sd(myBiomodModelEval["ROC","Testing.data",,,])

# Variable importances
var_importance<-get_variables_importance(myBiomodModelOut)
apply(var_importance,c(1,2),mean)

# Information about influence of environmental variables
myBiomod<-BIOMOD_LoadModels(myBiomodModelOut,models='MARS')

eval_strip<-response.plot2(models=myBiomod,
                               Data = get_formal_data(myBiomodModelOut,'expl.var'),
                               show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
                               do.bivariate = FALSE,
                               fixed.var.metric = 'median',
                               legend=FALSE,
                               display_title=FALSE,
                               data_species=get_formal_data(myBiomodModelOut,'resp.var')
                               )

# Ensemble Modeling
myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  em.by='all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.8),
  models.eval.meth = c('KAPPA','TSS','ROC'),
  prob.mean = F,
  prob.cv = T,
  committee.averaging = T,
  prob.mean.weight = T,
  VarImport = 0,
  prob.mean.weight.decay = 'proportional' )

# Summary
myBiomodEM

# Evaluation scores
EM <- get_evaluations(myBiomodEM)
EM

# Model Projection
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  selected.models = 'all',
  new.env = current_bios,
  proj.name = 'current',
  binary.meth = 'TSS',
  output.format = '.img',
  do.stack=FALSE)

# Ensemble models - From Guisan et al. (2017)
myBiomodProjEnsemble<-BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProj,
  binary.meth = 'TSS',
  output.format = '.img',
  do.stack = FALSE)

# summary of crated object
myBiomodProjEnsemble

# files created on hard drive
list.files("Toumeyella.parvicornis/proj_current/")

# make some plot of ensemble model by str.grep argument
plot(myBiomodProjEnsemble)





