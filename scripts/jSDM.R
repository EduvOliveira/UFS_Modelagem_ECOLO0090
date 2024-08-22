#' ---
# Disciplina: Modelagem de distribuição de espécies
# https://github.com/EduvOliveira/UFS_Modelagem_ECOLO0090
# by Eduardo V. S. Oliveira
# 22/08/2024
#' ---

###Aula JSDM###

## Instalando pacotes

if(!require(raster)) install.packages("raster")
if(!require(sp)) install.packages("sp")
if(!require(dplyr)) install.packages("dplyr")
if(!require(sf)) install.packages("sf")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(devtools)) install.packages("devtools")
if(!require(Hmsc)) install.packages("Hmsc")
if(!require(corrplot)) install.packages("corrplot")

library(devtools)
install_github("macroecology/letsR")
library(letsR)

## Entrando com as variaveis independentes

my_path<-"C:/Users/eduar/Documents/R/Disciplina_Modelagem_Ecologica"

setwd(paste(my_path,"/source", sep="")) 

f <-list.files(path = ".",pattern = '.*\\.(tif|asc)$') 
r <- lapply(f, raster) 

my_preds<-stack(r) 

names(my_preds) <- paste0(c("bio1","bio12","bio14","bio15","bio18","bio19","bio2","bio4"))

## Entrando com as ocorrencias

setwd(paste(my_path,"/occ", sep="")) 

occ<-read.csv("occ_jsdm.csv", h=TRUE)
head(occ)

## Criando uma matriz de presenca e ausencia

xy<-occ[,2:3]

spp<-as.vector(occ$spp)

# Plotando os pontos 

pol<-shapefile(choose.files()) 

plot(pol)
points(occ$lon, occ$lat, col='red', pch=20, cex=0.75) 

coords<-xy
coordinates(coords) <- ~lon+lat

# Limites dos pontos

xmin(coords)
xmax(coords)
ymin(coords)
ymax(coords)

PAM <- lets.presab.points(xy=xy, species=spp, xmn = -74.14, xmx = -34.97806,
                          ymn = -33.666, ymx = 10.61457, resol = 1, remove.cells=TRUE,remove.sp = TRUE)

summary(PAM)

## Adicionando valores do raster

PAM_var_mean <- lets.addvar(PAM, my_preds)

# Acessando a matriz
matriz<-as.data.frame(PAM_var_mean)

head(matriz)

nrow(matriz)

# Atribuindo nomes das celulas

cell<-as.factor(1:1010)

data<-cbind(matriz, cell)

# Organizando as celulas

# Mover "cell" para primeira coluna

data<-data %>% select(cell, everything())

# Renomear colunas das variaveis

data <- dplyr::rename(data, bio1 = bio1_mean, bio2 = bio2_mean, bio4 = bio4_mean, bio12 = bio12_mean, bio14 = bio14_mean, bio15 = bio15_mean, bio18 = bio18_mean, bio19 = bio19_mean)

head(data)

data_filter<-na.omit(data)

setwd(paste(my_path,"/occ", sep="")) 

write.csv(data_filter,"spp_inter.csv")

data<-read.csv("spp_inter.csv", head=TRUE, row.names =1)

## Configurando o modelo

# Montagem das "caixas" que serao introduzidas dentro do modelo 

#"Plots" como fator
levels(as.factor(data$cell))

#Transformando variaveis ambientais ("caixa" x)

cov = data[,18:ncol(data)]

names(cov)

X = as.data.frame(cov[,c("bio1","bio2","bio4","bio12","bio14","bio15","bio18","bio19")])


#Criando um objeto 'xformula', onde serao criados efeitos aditivos entre as variaveis

XFormula = ~ bio1 + bio2 + bio4 + bio12 + bio14 + bio15 + bio18 + bio19

head(X)

summary(X)

#"Caixa" com ocorrencias das especies

Y = as.matrix(data[,4:17])

dim(Y)

#Coordenadas das especies

xy = data[,2:3]

rownames(xy)=data$cell
colnames(xy)=c("x-coordinate","y-coordinate")

head(xy)

#"Caixa" design do estudo 
#Contem as variaveis representando a estrutura da analise

data$cell<-as.factor(data$cell)

studyDesign<-data.frame(cell = data$cell)
#studyDesign = data[,c("site","plot")]
#Se amostras aninhadas

rL = HmscRandomLevel(sData=xy)
#site = HmscRandomLevel(units = studyDesign$site)
#plot = HmscRandomLevel(units = studyDesign$plot)#Se amostras aninhadas

## Criando o modelo a partir das diferentes "caixas"
#Montar as "caixas" produzidas anteriormente, contendo os dados e a estrutrura do modelo
#Modelo "probit" para dados de presenca
#-ausencia
#(i.e., ocorrencia de especies)
#Modelos "Poisson" e "Lognormal Poisson" para dados de contagem 
#(i.e., abundancia das especies)

simul <- Hmsc(Y=Y, XData = X,XFormula = XFormula,studyDesign = studyDesign
              ,ranLevels=list(cell=rL),distr = "probit")  

## Rodando o modelo

#Aumentar valores para a analise final:

# Quantas amostras por cadeia?
#samples = 1000
#numero de amostras posteriores para salvar

# Quantas cadeias?
#nChains = 2
#numero de cadeias independentes para rodar
#longa cadeia necessaria para alcancar convergencia 

# Quanto "thinning" aplicar?
#thin = 50
#Parametro de "afinamento" da cadeia MCMC.
#Intervalo de "afinamento" da distribuicao posterior
#amostras posteriores consecutivas nao sao independentes
#por isso, nao registrar todas as amostras
#aumentar o thin reduz a correlacao entre amostras consecutivas

# Qual o comprimento de um "transitorio"?
#transient = 500*thin
#numero de etapas executadas antes de salvar as amostras #posteriores
#cadeias iniciam da iteracao 25.000
#quando thin = 1, transient = 500, com 1000 iteracoes

#Total de iteracoes por cadeia: (1000 x 50) + (500 x 50) = 75.000
#para obter 1000 amostras, vai precisar de 50.000 #iteracoes

samples = 50
nChains = 2
thin = 2 
transient = 5
#105 iteracoes por cadeia

t1=Sys.time()

mod_HMSC = sampleMcmc(simul,
                      samples = samples,
                      thin = thin,
                      transient = transient,
                      nChains = nChains, 
                      nParallel = nChains,initPar = "fixed effects")

t2=Sys.time()
t2-t1 

setwd(paste(my_path)) 

dir.create("saves5")

setwd(paste(my_path,"/saves5", sep="")) 

save(mod_HMSC,file="mod_HMSC.Rdata")

## Carregar o modelo

load("mod_HMSC.Rdata")

## Avaliacao do modelo

# Diagnostico de Gelman (devera ser proximo de 1)

mcoda <- convertToCodaObject(mod_HMSC)

gelman.diag(mcoda$Beta[,1:10])

## Principais resultados

# Interacoes entre especies

OmegaCor = computeAssociations(mod_HMSC)

#Ver resultados
OmegaCor

supportLevel = 0.95

toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
toPlot = sign(toPlot)

plotOrder = corrMatOrder(OmegaCor[[1]]$mean,order="AOE")

corrplot(toPlot[plotOrder,plotOrder], method = "color", tl.cex=1,
         col=colorRampPalette(c("blue", "white", "red"))(255))

dev.copy(device = jpeg, file = "interacoes.jpeg", width = 2000, height = 2000, res = 300)
dev.off()

## Fazendo a predicao do modelo

setwd(paste(my_path,"/pol", sep="")) 

grid<-st_read("grid_cerrado.shp")

plot(grid$geometry)

xy<-st_coordinates(st_centroid(grid))

dim(xy)

#Extraindo dados climaticos do Raster

clim.occ<-data.frame(raster::extract(my_preds,xy))

data<-cbind(xy, clim.occ)

head(data)

#Removendo NAs
data_clean<-data %>% na.omit()

#Configurando a predicao do modelo

xy.grid = as.matrix(cbind(data_clean$X,data_clean$Y))

XData.grid = data.frame(bio1=data_clean$bio1, bio2=data_clean$bio2,
                        bio4=data_clean$bio4, bio12=data_clean$bio12,
                        bio14=data_clean$bio14,
                        bio15=data_clean$bio15, bio18=data_clean$bio18,
                        bio19=data_clean$bio19, stringsAsFactors = TRUE)


Gradient = prepareGradient(mod_HMSC, XDataNew = XData.grid, sDataNew = list(cell=xy.grid))

nParallel=2

predY = predict(mod_HMSC, Gradient=Gradient, expected = TRUE, nParallel=nParallel)

#Probabilidade de ocorrencia da especie em cada ponto
head(predY[[2]])

#Espacializando as predicoes
EpredY=Reduce("+",predY)/length(predY)

Cm = EpredY[,2]
S=rowSums(EpredY)
CWM = (EpredY%*%mod_HMSC$Tr)/matrix(rep(S,mod_HMSC$nt),ncol=mod_HMSC$nt)
xy = data_clean[,1:2]
mapData=data.frame(xy,S,Cm,CWM,stringsAsFactors=TRUE)

#Predicao de uma especie

ggplot(data = mapData, aes(x=X, y=Y, color=Cm))+geom_point(size=2) + ggtitle(expression(italic("Tersina_viridis")))+ scale_color_gradient(low="blue", high="red") + coord_equal()

#Riqueza de especies
ggplot(data = mapData, aes(x=X, y=Y, color=S))+geom_point(size=2) + ggtitle("Riqueza de especies")+ scale_color_gradient(low="blue", high="red") + coord_equal()

#Obtendo a composicao por ponto (probabilidade 0-1)

cell1<-EpredY[1,]

###FIM###

rm(list=ls()) 

