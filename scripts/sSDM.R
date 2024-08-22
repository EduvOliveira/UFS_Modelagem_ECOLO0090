#' ---
# Disciplina: Modelagem de distribuição de espécies
# https://github.com/EduvOliveira/UFS_Modelagem_ECOLO0090
# by Eduardo V. S. Oliveira
# 22/08/2024
#' ---

###Aula SSDM###

###Passos para utilizar o algoritmo MAXENT
#(1) incluir o arquivo "maxent.jar" na pasta "java" do pacote "dismo"

# maxent
download.file(url = "https://biodiversityinformatics.amnh.org/open_source/maxent/maxent.php?op=download", destfile = paste0(system.file("java", package = "dismo"), "/maxent.zip"), mode = "wb")

unzip(zipfile = paste0(system.file("java", package = "dismo"), "/maxent.zip"), exdir = system.file("java", package = "dismo"), junkpaths = TRUE)

dir(system.file("java", package = "dismo"))

#(2) Java instalado no computador 

# install java - windows - https://www.java.com/pt_BR/download/
if(!require(rJava)) install.packages("rJava")

#--------------------


## Instalando pacotes

if(!require(raster)) install.packages("raster")
if(!require(sp)) install.packages("sp")
if(!require(virtualspecies)) install.packages("virtualspecies")
if(!require(SSDM)) install.packages("SSDM")


## Entrando com as variaveis independentes

my_path<-"C:/Users/eduar/Documents/R/Disciplina_Modelagem_Ecologica"

setwd(paste(my_path,"/source", sep="")) 

f <-list.files(path = ".",pattern = '.*\\.(tif|asc)$') 
r <- lapply(f, raster) 

preds<-stack(r) 

names(preds) <- paste0(c("bio1","bio12","bio14","bio15","bio18","bio19","bio2","bio4"))

## Preparando arquivo com as variaveis

var_path<-"C:/Users/eduar/Documents/R/Disciplina_Modelagem_Ecologica/source"

var<-load_var(path = var_path, files = c("bio1.tif","bio2.tif","bio4.tif","bio12.tif","bio14.tif","bio15.tif","bio18.tif","bio19.tif"), format = ".tif")

## Preparando arquivo com ocorrencias

path.occ<-setwd("C:/Users/eduar/Documents/R/Disciplina_Modelagem_Ecologica/occ") 

occ<-load_occ(path = path.occ, var, file = "occ_ssdm.csv", sep = ',', Xcol = "lon",Ycol = "lat", Spcol = "spp")

## plotando os pontos 

pol<-shapefile(choose.files()) 

plot(pol)
points(occ$lon, occ$lat, col='red', pch=20, cex=0.75) 

plot(var$bio12)
points(occ$lon, occ$lat, col='red', pch=20, cex=0.5) 

## Estimando a riqueza 

setwd(paste(my_path)) 

dir.create("saves4")

path.model<-setwd("C:/Users/eduar/Documents/R/Disciplina_Modelagem_Ecologica/saves4") 

t1=Sys.time()

SSDM <-stack_modelling("MAXENT", occ, var, Xcol = "lon",Ycol = "lat",
                Spcol = "spp", rep = 1,name = "modelo1", save = TRUE,
                path = path.model,cv = "k-fold", cv.param = c(4, 5),
                axes.metric = "Pearson",uncertainty = FALSE, ensemble.metric = "AUC",
                ensemble.thresh = 0.7,method = "bSSDM", metric = "SES")

t2=Sys.time()
t2-t1 

## Acessando resultados

plot(SSDM)

setwd(paste(path.model,"/modelo1/Stack/Rasters", sep="")) 

div<-raster("Diversity.tif")
plot(div)

setwd(paste(path.model,"/modelo1/Species/Vochysia_elliptica/Rasters", sep="")) 

bin<-raster("Binary.tif")
plot(bin)

## Acessando slot
# Estrutura do objeto "Stacked.SDM"

#Stack
#SSDM@name

SSDM@variable.importance

SSDM@diversity.map
#SSDM@endemism.map
#SSDM@evaluation
#SSDM@variable.importance
#SSDM@parameters

# Parametros dos modelos
slot(SSDM, "parameters")
SSDM@parameters

# Avaliacao dos modelos
slot(SSDM, "evaluation")

# Importancia das variaveis
slot(SSDM, "variable.importance")

#Species
#SSDM@esdms

SSDM@esdms$Vochysia_elliptica.Ensemble.SDM

#SSDM@esdms$Vochysia_thyrsoidea.Ensemble.SDM
#SSDM@esdms$Vochysia_tucanorum.Ensemble.SDM

# Avaliacao do modelo

SSDM@esdms$Vochysia_tucanorum.Ensemble.SDM@evaluation

# Acessando um raster binario

binary.map<-SSDM@esdms$Vochysia_elliptica.Ensemble.SDM@binary

plot(binary.map)

setwd(paste(my_path,"/saves4", sep="")) 

writeRaster(binary.map, "sp1.asc")

adeq<-SSDM@esdms$Vochysia_elliptica.Ensemble.SDM@projection

plot(adeq)

###FIM###

rm(list=ls()) 

