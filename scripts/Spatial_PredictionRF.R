#' ---
# Disciplina: Modelagem de distribuição de espécies
# https://github.com/EduvOliveira/UFS_Modelagem_ECOLO0090
# by Eduardo V. S. Oliveira
# 22/08/2024
#' ---

###Aula Random Forest###

if(!require(terra)) install.packages("terra")
if(!require(randomForest)) install.packages("randomForest")
if(!require(raster)) install.packages("raster")
if(!require(Metrics)) install.packages("Metrics")
if(!require(rpart)) install.packages("rpart")
if(!require(rpart.plot)) install.packages("rpart.plot")

#Diretorio de trabalho
workspace <- "C:/Users/eduar/Documents/R/Disciplina_Modelagem_Ecologica/"

setwd(paste(workspace)) 

dir.create("RF")

setwd(paste(my_path,"/RF", sep="")) 

# download variables ----
download.file(url = "https://github.com/EduvOliveira/UFS_Modelagem_ECOLO0090/raw/main/RF.zip", destfile = "RF.zip")

# unzip
unzip(zipfile = "RF.zip")

###### Carregar dados #####

AF_points <- terra::vect(paste0(workspace, "RF/sitesAF.shp"))

biomass <- terra::rast(paste0(workspace, "RF/agb_Ferreira_et_al_2017.tif"))
r1 <- terra::rast(paste0(workspace, "RF/bio1.tif"))
r2 <- terra::rast(paste0(workspace, "RF/bio12.tif"))
r3 <- terra::rast(paste0(workspace, "RF/evi2.tif"))

predictorStack<-c(biomass,r1,r2,r3) 

predictorStack
names(predictorStack)

# Renomear as camadas
names(predictorStack) <- c("agb","bio1","bio12","evi2")

# Extrair valores das camadas a partir dos pontos
pointVals <- extract(predictorStack, AF_points)

# Tabela final
modelDF <- cbind(data.frame(AF_points), pointVals)

head(modelDF)
names(modelDF)

###### Explorar os relacionamentos entre as variaveis

modelDF$bio1[is.na(modelDF$bio1)] <- 21
modelDF$bio12[is.na(modelDF$bio12)] <- 1200

plot(modelDF$agb, modelDF$bio1)
plot(modelDF$agb, modelDF$bio12)
plot(modelDF$agb, modelDF$evi12)

cor(modelDF$agb, modelDF$bio1)
cor(modelDF$agb, modelDF$bio12)
cor(modelDF$agb, modelDF$evi2)

plot(biomass, 1)
plot(r1, 1)
plot(r2, 1)


##### Ajustando um modelo de regressão de árvores aleatórias para biomassa acima do solo

##Primeiro modelo usará todas as variáveis

model1 <- agb ~ bio1 + bio12 + evi2

# Ajusta um modelo RF
rf1 <- randomForest(model1, data=modelDF)

# Importancia das variáveis
varImpPlot(rf1)

# Ajusta um modelo simplificado (o melhor)
model2 <- agb ~ bio1 + bio12 

rf2 <- randomForest(model2, data=modelDF)
rf2
varImpPlot(rf2)	

########## Aplicando o modelo para criar uma camada de predicao espacial

# Carregar as variaveis para o futuro
p1 <- terra::rast(paste0(workspace, "RF/bio1.85.tif"))
p2 <- terra::rast(paste0(workspace, "RF/bio12.85.tif"))

predictorStackSubset<-c(p1,p2)

names(predictorStackSubset) <- c("bio1","bio12")

#Aplicar o modelo para criar um raster de predição contínua (predição da biomassa total acima do solo) usando as camadas preditoras.
predRaster <- predict(predictorStackSubset, rf2)

par(mfrow=c(1,2)) # Configure para uma linha e uma coluna
plot(biomass, main = "Atual")
plot(predRaster,main = "Cenário Pessimista")
par(mfrow=c(1,1))#Padrao

dev.copy(device = jpeg, file = "AGB_AtlaticForest.jpeg", width = 1800, height = 1500, res = 300)
dev.off()

# Salvando
writeRaster(predRaster,
            paste0(workspace, "RF/biomass_prediction.tif"))

### Avaliacao do modelo

pred_model <- predict(rf2)

plot(pred_model)
#Maioria dos valores entre 180-200 t/ha

# RMSE

rmse_model <- rmse(modelDF$agb, pred_model)
rmse_model#mesma unidade do output
#quanto menor o valor do RMSE, mais acuradas sao as previsoes do modelo

print(rmse_model/mean(modelDF$agb)) 

# Arvores de decisao

model_rpart <- rpart(model2, data = modelDF)
rpart.plot(model_rpart)

###FIM###
