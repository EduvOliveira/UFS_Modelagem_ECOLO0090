#' ---
# Disciplina: Introdução à Modelagem Ecológica 
# https://github.com/EduvOliveira/UFS_Modelagem_ECOLO0090
# by Eduardo V. S. Oliveira
# 20/08/2024
#' ---

# Modelagem de distribuição de espécies----------------------------

#####Aula: Entrada de dados de ocorrencia de espécies / Filtragem e limpeza de dados de ocorrência de espécies
#####

### Obtendo / manipulando dados de ocorrencias das especies

if(!require(dismo)) install.packages("dismo")
if(!require(dplyr)) install.packages("dplyr")
if(!require(rgbif)) install.packages("rgbif")
if(!require(plyr)) install.packages("plyr")
if(!require(CoordinateCleaner)) install.packages("CoordinateCleaner")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(sf)) install.packages("sf")
if(!require(raster)) install.packages("raster")
if(!require(geodata)) install.packages("geodata")
if(!require(virtualspecies)) install.packages("virtualspecies")
if(!require(mapedit)) install.packages("mapedit")
if(!require(mapview)) install.packages("mapview")
if(!require(rnaturalearth)) install.packages("rnaturalearth")

my_path<-"C:/Users/eduar/Documents/R/Disciplina_Modelagem_Ecologica"

## GBIF (gbif.org)

## Extraindo ocorrencias de especies

gbif_data<-occ_search(scientificName=c("Tersina viridis viridis", "Elaenia cristata",                          "Thraupis sayaca sayaca","Schistochlamys ruficapillus ruficapillus"),hasCoordinate = TRUE, field = c("scientificName","decimalLatitude","decimalLongitude","country"))

gbif_data$"Tersina viridis"$data

bdspp<-rbind.fill(gbif_data$"Tersina viridis viridis"$data)

occ_aves <- dplyr::rename(bdspp, spp = scientificName, lon = decimalLongitude, lat = decimalLatitude)

head(occ_aves)

dir.create("saves")

setwd(paste(my_path,"/saves", sep="")) 

write.table(occ_aves,"Occ_aves.txt")

# Limpeza de dados

# Registros duplicados

dups <- duplicated(occ_aves[, c('spp', 'lon', 'lat')]) 
sum(dups)
occ.aves_sub <- occ_aves[!dups, ] 
dim(occ.aves_sub)

# Removendo NAs

occ_aves_clean<-na.omit(occ.aves_sub)
dim(occ_aves_clean)

write.table(occ_aves_clean,"Occ_aves_clean.txt")

# Plotando os dados

setwd(paste(my_path,"/pol", sep="")) 

pol1<-st_read("america_sul.shp")

plot(pol1$geometry, xlim=c(-60,-30), ylim=c(-60,20), axes=TRUE, col="light yellow")

points(occ_aves_clean$lon, occ_aves_clean$lat, col='red', pch=20, cex=0.75)

occ_birds <- sf::st_as_sf(occ_aves_clean, 
                       coords = c("lon", "lat"),
                       crs = 4326) 

world <- rnaturalearth::ne_countries(scale = "small",
                                     returnclass = "sf") 

ggplot() +
  geom_sf(data = world,
          mapping = aes(geometry = geometry),
          fill = "white") +
  geom_sf(data = occ_birds,
          aes(geometry = geometry),
          size = 3,
          color = "red") +
  theme_bw()

mapview::mapview(occ_birds)

# Edicao manual ----
occ_data_filter_edit <- mapedit::editFeatures(occ_birds) # atencao para o Done.
occ_data_filter_edit

mapview::mapview(occ_data_filter_edit)

# tabela
occ_coords<-occ_data_filter_edit %>% sf::st_coordinates() 
occ_drop<-occ_data_filter_edit %>% sf::st_drop_geometry() 

occ_filter<-cbind(occ_coords,occ_drop[,1])

head(occ_filter)

setwd(paste(my_path,"/saves", sep="")) 

write.table(occ_filter,"Occ_aves_filter.txt")

## Opcao 2 - GBIF

spp <- gbif(genus = "Tropidurus", species = "hygomi",removeZeros = TRUE) 

dim(spp)
sp_sub<-dplyr::select(spp, lon, lat, species)
head(sp_sub)

setwd(paste(my_path,"/saves", sep="")) 

write.table(sp_sub,"Occ_tropidurus.txt")

# Removendo NAs

occ_dados<-na.omit(sp_sub)
dim(occ_dados)

# Limpeza automatica dos dados 
# Um dataframe com apenas registros limpos

occ_filter <- CoordinateCleaner::clean_coordinates(
  x = occ_dados,
  species = "species",
  lon = "lon",
  lat = "lat",
  outliers_mtp = 2,
  tests = c("capitals", # radius around capitals
            "centroids", # radius around country and province centroids
            "duplicates", # records from one species with identical coordinates
            "equal", # equal coordinates
            "gbif", # radius around GBIF headquarters
            "institutions", # radius around biodiversity institutions
            "seas", # in the sea
            "urban", # within urban area
            "validity", # outside reference coordinate system
            "zeros" # plain zeros and lat = lon
  ))

summary(occ_filter)

plot(occ_filter, lon = "lon", lat = "lat")

#Excluindo registros problematicos
dat_cl <- occ_dados[occ_filter$.summary,]
dim(dat_cl)

# Replotando os resultados da limpeza 
wm <- borders("world", colour="gray50", fill="gray50")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = as.data.frame(dat_cl), aes(x = lon, y = lat),
             colour = "darkblue", size = 2.5)+theme_bw()

setwd(paste(my_path,"/saves", sep="")) 

write.table(occ_filter,"Occ_tropidurus_clean.txt")


# Avisos====

## speciesLink (<https://specieslink.net/search/>)
## Dados de ocorrencias de plantas a partir do pacote 'plantR' e 'BIEN'  
## Baixar dados manualmente para outros grupos
## Pacote 'spfilt' e 'bdc' são outras opções para limpeza dos dados 

##Baixando camadas ambientais

setwd(paste(my_path)) 

bio<-worldclim_global("bio", res=10, path =my_path)

names(bio)

preds<-stack(bio)

## Entrando com as variaveis

setwd(paste(my_path,"/wc2.1_10m", sep="")) 

lst <- list.files(path=".",pattern='tif$',full.names = T) 
preds<-stack(lst)
names(preds)

names(preds) <- gsub("wc2.1_10m_", " " ,names(preds))

## Recortando para a America do Sul

pol<-shapefile(choose.files()) 

plot(pol)

preds_AS<-list()
for(i in 1:19) preds_AS[[i]]<-mask(preds[[i]],pol)

preds_AS<-stack(preds_AS)
preds_AS<-crop(preds_AS, c(-85, -30, -60, 15))

sr=stack(preds_AS)

class(sr)

plot(sr[[1]])
plot(pol, add=T)

## Checando multicolinearidade

var<-removeCollinearity(sr, multicollinearity.cutoff = 0.7, select.variables = TRUE, sample.points = FALSE, plot = TRUE)

setwd(paste(my_path,"/saves", sep="")) 

dev.copy(device = jpeg, file = "Multicolinearidade.jpeg", width = 1800, height = 1500, res = 300)
dev.off()

## Selecionando as variaveis

bio1<-sr$bio_1
bio2<-sr$bio_2
bio4<-sr$bio_4
bio12<-sr$bio_12
bio14<-sr$bio_14
bio15<-sr$bio_15
bio18<-sr$bio_18
bio19<-sr$bio_19

my_preds<-stack(bio1,bio2,bio4,bio12,bio14,bio15,bio18,bio19)
plot(my_preds)

setwd(paste(my_path)) 

dir.create("var")

setwd(paste(my_path,"/var", sep="")) 

writeRaster(my_preds, file.path(setwd(paste(my_path,"/var", sep="")), names(my_preds)), bylayer = T,format="GTiff")

#Uma variavel
#writeRaster(my_preds$bio1,"bio1.tif",format="GTiff")

#--------------------

## Exercicio
# Obter registros de ocorrência de uma espécie focal
# Realizar os procedimentos de filtragem e limpeza

###Fim do script###

