#' ---
# Disciplina: Modelagem de distribuição de espécies
# https://github.com/EduvOliveira/UFS_Modelagem_ECOLO0090
# by Eduardo V. S. Oliveira
# 21/08/2024
#' ---

# Modelagem de distribuição de espécies----------------------------

#####Aula: Calibração e construção de SDMs / Projetando SDMs para diferentes cenários climáticos

##Passos para utilizar o algoritmo MAXENT
#(1) incluir o arquivo "maxent.jar" na pasta "java" do pacote "dismo"

# maxent
download.file(url = "https://biodiversityinformatics.amnh.org/open_source/maxent/maxent.php?op=download", destfile = paste0(system.file("java", package = "dismo"), "/maxent.zip"), mode = "wb")

unzip(zipfile = paste0(system.file("java", package = "dismo"), "/maxent.zip"), exdir = system.file("java", package = "dismo"), junkpaths = TRUE)

dir(system.file("java", package = "dismo"))

#(2) Java instalado no computador 

# install java - windows - https://www.java.com/pt_BR/download/
if(!require(rJava)) install.packages("rJava")

#--------------------

if(!require(spThin)) install.packages("spThin")
if(!require(sdm)) install.packages("sdm")
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(raster)) install.packages("raster")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(dismo)) install.packages("dismo")

##Entrando com dados de ocorrência

my_path<-"C:/Users/eduar/Documents/R/Disciplina_Modelagem_Ecologica"

setwd(paste(my_path,"/occ", sep="")) 

occ<-read.table("Occ_capybara.txt", header = TRUE)

# Tratando viés de amostragem
# Filtragem dos pontos
# Koldasbayeva et al. 2022 (thinning > 10km)

points<-occ[,1:3]

points_thinned <- thin(loc.data = points, 
                       lat.col = "decimallatitude", long.col = "decimallongitude", 
                       spec.col = "species", 
                       thin.par = 25, reps = 100, 
                       locs.thinned.list.return = TRUE, 
                       write.files = TRUE, 
                       max.files = 3, 
                       out.dir = getwd(), out.base = "points_thinned", 
                       write.log.file = T,
                       log.file = "points_thinned.txt" )

sp <-(points_thinned)
wm <- borders("world", colour="gray50", fill="gray50", xlim = c(-80, -30), ylim = c(-60, 10))
ggplot()+ coord_fixed()+ wm +
  geom_point(data = sp[[1]], aes(x = Longitude, y = Latitude),
             colour = "black", size = 1)

setwd(paste(my_path,"/saves", sep="")) 

dev.copy(device = jpeg, file = "points_thinned_25km.jpeg", width = 1800, height = 1500, res = 300)
dev.off()

write.csv(sp[[1]],"train_points.csv")

write.csv(sp[[2]],"test_points.csv")

##Entrando com as variaveis  

setwd(paste(my_path)) 

dir.create("source")

setwd(paste(my_path,"/source", sep="")) 

# download variables ----
download.file(url = "https://github.com/EduvOliveira/UFS_Modelagem_ECOLO0090/raw/main/var/wc_10.zip", destfile = "wc_10.zip")

# unzip
unzip(zipfile = "wc_10.zip")

lst <- list.files(path=".",pattern='tif$',full.names = T) 
preds<-stack(lst)
names(preds)

bio1<-preds$bio1
bio2<-preds$bio2
bio4<-preds$bio4
bio12<-preds$bio12
bio14<-preds$bio14
bio15<-preds$bio15
bio18<-preds$bio18
bio19<-preds$bio19

my_preds<-stack(bio1,bio2,bio4,bio12,bio14,bio15,bio18,bio19)

##Criando pontos de background aleatorios
#Criar 10000 pontos (Barbet-Massin et al. 2016)
  
bg <- randomPoints(my_preds$bio1, n=1000) 
bg <- as.data.frame(bg) 

#Entrando com pontos de ocorrencia

setwd(paste(my_path,"/saves", sep="")) 
  
train<-read.csv("train_points.csv",h=TRUE)
head(train)
  
train_ptos<-train[,2:3]

##Pontos de treino

var_train<-raster::extract(my_preds,train_ptos)
train<-cbind(train_ptos,var_train)
train$pre.abs<-rep(1,times=nrow(train))
head(train)
coordinates(train)<- ~ Longitude+Latitude
crs(train) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
##Pontos de teste
  
test<-read.csv("test_points.csv",h=TRUE)
head(test)
  
test_ptos2<-test[,2:3]
  
var.tst<-raster::extract(my_preds,test_ptos2)
test<-cbind(test_ptos2,var.tst)
test$pre.abs<-rep(1,times=nrow(test))
head(test)
coordinates(test)<- ~ Longitude+Latitude
crs(test) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
##Background
var.bg<-raster::extract(my_preds,bg)
bg<-cbind(bg,var.bg)

head(bg)
colnames(bg)<-c("Longitude","Latitude", "bio1", "bio2", "bio4", "bio12", "bio14","bio15","bio18","bio19")

bgplot=bg

coordinates(bgplot)<- ~ Longitude+Latitude #para visualizar 

plot(preds$bio1)
plot(train, add=TRUE)
plot(test, add=TRUE)
plot(bgplot, add=TRUE)
  
library(sdm)
  
installAll()
  
##Gerando o modelo

setwd(paste(my_path)) 
  
dir.create("saves2")

setwd(paste(my_path,"/saves2", sep="")) 
  
t1=Sys.time()
  
#bg como data.frame
#Cria um objeto 'sdmdata'
data_occ<- sdmData(formula = pre.abs ~ bio1+bio2+bio4+bio12+bio14+bio15+bio18+bio19, train = train,test = test, predictors = my_preds, bg = bg)
  
t2=Sys.time()

t2-t1 
  
write.sdm(data_occ,"sdm_result.sdd",overwrite = TRUE)
  
#read.sdm(sdm_result.sdd")
  
t1=Sys.time()
  
model_sdm<-sdm(pre.abs ~ bio1+bio2+bio4+bio12+bio14+bio15+bio18+bio19, data_occ, methods='maxent', replication='cv', n=10) 
  
t2=Sys.time()

t2-t1 
  
write.sdm(model_sdm,"model_capybara.sdm",overwrite = TRUE)
  
#read.sdm("model_sdm.sdm")

#sem dados de teste independente
#model_ESE_S<-sdm(pre.abs ~ bio7+bio9+bio11+bio17+bio18+ph+rad, data_ESE_S, methods=c("maxent","rf","bioclim"), test.percent=25,n=2, replication='cv',cv.folds=5) 
  
#Acessar modelos por algoritmo
#model_DNE@models$pre.abs$bioclim
  
getVarImp(model_sdm)
#media dos modelos
  
gui(model_sdm)
  
#Obtendo a importancia das variaveis do modelo
vi<-getVarImp(model_sdm,id=1,wtest='training')
plot(vi,'auc')

dev.copy(device = jpeg, file = "import_var.jpeg", width = 1800, height = 1500, res = 300)
dev.off()
  
#Obtendo as medidas de avaliacao de todos os modelos
ev<-getEvaluation(model_sdm,wtest='test.indep', stat=c('TSS','Kappa','AUC','threshold','sensitivity',                              'specificity'),opt=2)
  
#obtendo a media por medida
med_auc<-mean(ev$AUC)
med_tss<-mean(ev$TSS)
med_kap<-mean(ev$Kappa)
med_th<-mean(ev$threshold)
med_se<-mean(ev$sensitivity)
med_sp<-mean(ev$specificity)
  
med_auc
med_tss
med_kap
med_th
med_se
med_sp
  
#Elaborar o mapa
#Fazendo a media ponderada dos modelos
  
t1=Sys.time()
  
ens<-ensemble(model_sdm, my_preds, setting=list(method='weighted',stat='TSS',opt = 2))
  
t2=Sys.time()
t2-t1
  
setwd(paste(my_path,"/saves2", sep="")) 
  
plot(ens)
  
writeRaster(ens, filename= "sdm_capybara.tif", format="GTiff", overwrite=TRUE)
  
##Fazendo o ensemble para cenarios futuros

setwd(paste(my_path)) 

# download variables ----
download.file(url = "https://github.com/EduvOliveira/UFS_Modelagem_ECOLO0090/raw/main/var/cmip5_85.zip", destfile = "cmip5_85.zip")

# unzip
unzip(zipfile = "cmip5_85.zip")

setwd(paste(my_path,"/cmip5_85", sep="")) 

bio1.85<-raster("bio1.85.tif")
bio2.85<-raster("bio2.85.tif")
bio4.85<-raster("bio4.85.tif")
bio12.85<-raster("bio12.85.tif")
bio14.85<- raster("bio14.85.tif")
bio15.85<- raster("bio15.85.tif")
bio18.85<- raster("bio18.85.tif")
bio19.85<- raster("bio19.85.tif")

preds_85<-stack(bio1.85,bio2.85,bio4.85,bio12.85,bio14.85,bio15.85,bio18.85,bio19.85)

names(preds_85)<-c("bio1","bio2","bio4","bio12","bio14","bio15","bio18","bio19")

ens2<-ensemble(model_sdm, preds_85, setting=list(method='weighted',stat='TSS',opt = 2))

plot(ens2)

setwd(paste(my_path)) 

dir.create("saves3")

setwd(paste(my_path,"/saves3", sep="")) 

writeRaster(ens2, filename= "ssp585_capybara.tif", format="GTiff", overwrite=TRUE)

#Raster binario

current<-ens
current[current < 0.49]<-NA  #Definir valor de threshold conforme resultado
plot(current, col = heat.colors(10, rev = TRUE))

cur.bin<-current
cur.bin[cur.bin > 0]<-1
plot(cur.bin)
writeRaster(cur.bin, "bin_capybara.asc",overwrite=TRUE)


###FIM###


