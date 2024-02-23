

library(terra)
library(tmap)
library(tidyverse)
library(sp)  
library(geodata)
library(treemap)
library(forcats)


#Folder
setwd('/Users/ccruzr/Library/Mobile Documents/com~apple~CloudDocs/Cristian/Documents/Trabajo/publicaciones/Data paper mammals/DP_mammal_Col')
new.crs<-'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

## basemaps
## Colombia
colombia <- (gadm('COL' , level = 1, path=tempdir()))
colo <- as(gadm('COL' , level = 0, path=tempdir()), 'Spatial')
## World
data("World")
suda <- as(World[World$continent %in% 'South America',], 'Spatial')
suda<- suda[!suda$iso_a3 %in% c('FLK'),]

#ext
extcol <- as.polygons(ext(-80, -65, -4.9, 12.9), crs="WGS84")
plot(extcol)
## Island
col_nsprov<- as(colombia[!colombia$NAME_1 %in% 'San Andrés y Providencia',], 'Spatial')
syprov<- as(colombia[colombia$NAME_1 %in% 'San Andrés y Providencia',], 'Spatial')

####
#points
dt<- data.table::fread('03processedData/Datos_verificados.txt')
unif <- dt[!is.na(dt$decimalLatitude),]
unif <- dt[!is.na(dt$decimalLongitude),]
unif$coords <- paste(unif$decimalLatitude, unif$decimalLongitude, unif$basisOfRecord)
unif$sites <- paste(unif$decimalLatitude, unif$decimalLongitude)
nrow(unif[!duplicated(unif$sites),]) ## unique sites
unif <- unif[!duplicated(unif$coords),]
unif$basisOfRecord <- gsub('HumanObservation', 'Human Observation', gsub('MachineObservation', 'Machine Observation', gsub('PreservedSpecimen', 'Preserved Specimen', unif$basisOfRecord)))
table(unif$basisOfRecord)
coordinates(unif)<- ~ decimalLongitude + decimalLatitude
crs(unif) <- '+proj=longlat +datum=WGS84'
unif_sf <- sf::st_as_sf(unif, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)


###
# Statistics

#basis of records
BoR <- data.frame(table(dt$basisOfRecord))
BoR$percentaje <- round((BoR$Freq / 97944) * 100,2)
## deparments
Depto <- data.frame(table(dt$stateProvince))
Depto$percentaje <- round((Depto$Freq / 97944) * 100,2)
## deparments
Mpio <- data.frame(table(dt$county))
Mpio$percentaje <- round((Mpio$Freq / 97944) * 100,2)
## Orders
ord_count <- data.frame(table(dt$order))
ord_count$percentaje <- round((ordn$Freq / 97944) * 100,2)
names(ord_count) <- c('Order', "Frecuency", "Percentaje")
## Families
fam_count <- dt %>% 
  count(family, order, sort = T)
fam_count$percentaje <- round((fam_count$n / 97944) * 100,2)
names(fam_count) <- c('Family', 'Order', "Frecuency", "Percentaje")
table(fam_count$Order)
#Species
spe_count <- dt %>% 
  count(scientificName, family, order, sort = T)
spe_count$percentaje <- round((spe_count$n / 97944) * 100,2)
## most common family
table(spe_count$family)
## dates
table(dt$year)

########
# Table and Figures

# Resolución 0126 de 2024 expedida por el Ministerio de Ambiente y Desarrollo Sostenible
# https://doi.org/10.15472/frowz3
download.file("https://ipt.biodiversidad.co/sib/archive.do?r=especies-amenazadas-mads-2024", "./03ProcessedData/test.zip")
unzip("03ProcessedData/test.zip", exdir = "./03ProcessedData/taxon_distribution")
file.remove("03ProcessedData/test.zip")

res <- read_tsv("03ProcessedData/taxon_distribution/taxon.txt")
re2 <- read_tsv("03ProcessedData/taxon_distribution/distribution.txt")
resCol <- res %>%
     full_join(re2, c("id")) %>%
     filter(scientificName %in% spe_count$scientificName) %>%
     select(scientificName, threatStatus) %>%
     rename(MADS = threatStatus) 

## UICN Red List
# los datos de todas las especies presentes en la IUCN
sp<-rredlist::rl_sp(all = TRUE, key = "441e9e03832696e08879d3ff84847b5422c897d9dfd87b65f63550ea4030c78b")
#unir resultados en una tabla
all_df <- do.call(rbind, lapply(sp, "[[", "result"))
thrIUCN <- all_df %>%
     filter(scientific_name %in% spe_count$scientificName) %>%
     select(scientific_name, category) %>%
     rename(IUCN = category)%>%
     rename(scientificName = scientific_name)%>%
     distinct(., scientificName, .keep_all = T)
#Tursiops truncatus LC
#Balaenoptera physalus VU

## CITES
cit <- data.table::fread('00RawData/Index_of_CITES_Species.csv')
citT <- cit %>%
     filter(FullName %in% spe_count$scientificName) %>%
     select(FullName, CurrentListing) %>%
     rename(CITES = CurrentListing)%>%
     rename(scientificName = FullName)

## Final table
tax <- data.table::fread('03ProcessedData/taxonomy_spp_datapaper_correct_Hera.csv')

tableF <- dt %>%
  group_by(order, family, scientificName, basisOfRecord) %>%
  summarize(Count = n()) %>%
  pivot_wider(names_from = basisOfRecord, values_from = Count)%>%
  as.data.frame()
tableF$PreservedSpecimen <- gsub('Invalid Number', 0, tableF$PreservedSpecimen)

#Final table
tableF <- tableF %>% 
  left_join(resCol, c("scientificName"))%>%
  left_join(thrIUCN, c("scientificName"))%>%
  left_join(citT, c("scientificName")) %>%
  left_join(tax, c("scientificName")) %>%
  mutate(scientificName = paste(scientificName, taxonomic_authority))%>%
  select(order, family, scientificName, CITES, IUCN, MADS, PreservedSpecimen, MachineObservation, HumanObservation)%>%
  mutate(Endemic = '', .before = PreservedSpecimen )

data.table::fwrite(tableF, '05OutData/Table2.csv')


#map1
tmf<- tm_shape(as(extcol, "Spatial"))+
  tm_polygons(col = "#a3a3a3")+
  tm_shape((crop(suda, c(-82, -65, -4.9, 12.9)))) +
  tm_polygons(col = "#d7d7d7")+
  tm_shape(col_nsprov) +
  tm_grid(alpha = 0.3, col = '#f7f7f7')+
  tm_polygons()+
  tm_shape(unif_sf) +
  tm_compass(position = c("left", "bottom"), size = 4) +
  tm_scale_bar(position = c("left", "bottom", size=5))+
  tm_dots(border.col = "blue",  col = 'basisOfRecord',palette = c("#4a5ef9","#11541b","#ff8000"), size=0.3, shape=20, legend.show = TRUE, "Basis of Record")+
  tm_layout(legend.bg.color="#ffffff", legend.position = c(0.75, 0.85))

# Create viewport
vp <- grid::viewport(x = 0.15, y = 0.85, width = 0.15, height = 0.15)
vp2 <- grid::viewport(x = 0.25, y = 0.85, width = 0.2, height = 0.12)
vp3 <- grid::viewport(x = 0.83, y = 0.22, width = 0.25, height = 0.3)
#vp4 <- grid::viewport(x = 0.7, y = 0.30, width = 0.05, height = 0.05)
vpa<- list(vp, vp2, vp3)

##Framework
satm<-tm_shape(syprov, bbox=tmaptools::bb(matrix(c(-81.74,12.47,-81.68,12.60),2,2))) +tm_polygons() #san Andres
optm<-tm_shape(syprov, bbox=tmaptools::bb(matrix(c(-81.41,13.31,-81.34,13.40),2,2))) +tm_polygons() #providencia
sud <- tm_shape(suda)+ tm_polygons() + tm_shape(colo)+
     tm_polygons(col = "#000000") + tm_grid(alpha = 0.4, col = '#a7a7a7', labels.size = 0.3)# sudamerica

is_tm<- list(satm, optm, sud)
tmap_save(tmf, insets_tm= is_tm, insets_vp= vpa, filename=("04Plots/recods_BoR.jpeg"), dpi = 300)

#####3
### Diagram

## Family plot
bar_plot <- fam_count %>%
  ggplot(aes(x = factor(Family, levels = Family[order(Frecuency)]), y = Frecuency)) +
  geom_bar(stat = "identity", fill = "#f68060", color = "#fff5f5", size = 0.5) +  # Add border to bars
  geom_text(aes(label = Frecuency), vjust = 0.8, hjust = -0.1, color = "#000000", size = 3) + # Add labels for total amount
  scale_y_continuous(expand = c(0, 0.05)) +  # Adjust y-axis limits
  theme_bw() +
  theme(panel.grid.major.y  = element_blank(),
    panel.grid.major.x = element_line(color = "#c2c2c2"),  # Change color of x-axis gridlines
        #panel.grid.minor = element_blank(),  # Remove minor gridlines
        axis.line.y = element_line(color = "black"),  # Change color of y-axis line
        axis.text.y = element_text(color = "black"),  # Change color of y-axis text
        axis.ticks.y = element_line(color = "black")) +  # Change color of y-axis ticks
  coord_flip(ylim = c(0, max(fam_count$Frecuency) * 1.1)) +  # Adjust y-axis limits
  labs(x = "Family", y = "Number of Records") +
  scale_size_manual("Count", values = c(5, 1.4), guide = "none") +
  scale_x_discrete(breaks = unique(fam_count$Family))  # Increase the number of lines on x-axis

# Export the modified plot
ggsave("04Plots/recods_per_Family.jpeg", plot = bar_plot, width = 8, height = 6, dpi = 300)

# Order plot
bar_plot_ord <- ord_count %>%
  ggplot(aes(x = factor(Order, levels = Order[order(Frecuency)]), y = Frecuency)) +
  geom_bar(stat = "identity", fill = "#006527", color = "#fff5f5", size = 0.5) +  # Add border to bars
  geom_text(aes(label = Frecuency), vjust = 0.8, hjust = -0.1, color = "#000000", size = 3) + # Add labels for total amount
  scale_y_continuous(expand = c(0, 0.05)) +  # Adjust y-axis limits
  theme_bw() +
  theme(panel.grid.major.y  = element_blank(),
    panel.grid.major.x = element_line(color = "#c2c2c2"),  # Change color of x-axis gridlines
        #panel.grid.minor = element_blank(),  # Remove minor gridlines
        axis.line.y = element_line(color = "black"),  # Change color of y-axis line
        axis.text.y = element_text(color = "black"),  # Change color of y-axis text
        axis.ticks.y = element_line(color = "black")) +  # Change color of y-axis ticks
  coord_flip(ylim = c(0, max(ord_count$Frecuency) * 1.1)) +  # Adjust y-axis limits
  labs(x = "Order", y = "Number of Records") +
  scale_size_manual("Frecuency", values = c(5, 1.4), guide = "none") +
  scale_x_discrete(breaks = unique(ord_count$Order))  # Increase the number of lines on x-axis


ggsave("04Plots/recods_per_order.jpeg", plot = bar_plot_ord, width = 8, height = 6, dpi = 300)
