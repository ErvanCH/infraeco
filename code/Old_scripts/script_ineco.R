library("sf")
library("data.table")
library("tidyverse")


# Import guildes
gui <- read.csv2("data/Table_récap_milieu2.txt")
setDT(gui)
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
gui[,name:=trim(taxon_orig)]
gui <- unique(gui[,-c("taxon_orig","name_fr","habitat_orig")])
View(gui[order(v_taxon_id)])

# Remove duplicates
gui[,index:=rowSums(.SD), .SDcols = 2:3]
A <- gui[duplicated(v_taxon_id)]
AA <- gui[,.I[duplicated(gui[,.(v_taxon_id,milieu_fr)][order(v_taxon_id,index)])]]
gui2 <- gui[order(v_taxon_id,index)][-AA,]
head(gui2[v_taxon_id%in%v_taxon_id[duplicated(v_taxon_id)]][order(v_taxon_id,index)])
table(gui2[,duplicated(v_taxon_id),by=milieu_fr]$V1)  # pas de duplicat par milieu

# Import obs
obs <- read.csv("data/extraction_DB_08072019.txt")
setDT(obs)

# Add milieu
obs2 <- merge(obs,gui2[,c(1,2,3,5:7)],by.x="project_taxon_id",by.y="v_taxon_id",all.x=T,allow.cartesian=TRUE)

head(obs2[obs_id%in%obs_id[duplicated(obs_id)]])

# Add subspecies to guildes
sub <- obs2[,unique(v_accepted_taxon_id),by=project_taxon_id]
names(sub) <- c("project_taxon_id","ssp")
A <- sub[,length(unique(ssp)),by=project_taxon_id]
AA <- sub[project_taxon_id%in%A[V1!=1,project_taxon_id]]
gui3 <- merge(gui2,AA,by.x="v_taxon_id",by.y="project_taxon_id",all.x=T,allow.cartesian=TRUE)

# Add spp names
SPP <- read.csv2("data/api_potential_taxa.txt")
gui3 <- merge(gui3,SPP[,c(1,6)],by.x="ssp",by.y="id",all.x=T)

# Export table d'espèces
gui4 <- gui3[,c("v_taxon_id", "name","original_name","characteristic", "dominant", "indigenous", 
          "description_fra", "milieu_fr", "milieu_de")]
names(gui4) <- c("id_taxon","name","name_spp","characteristic","dominant","indigenous","conservation","guild_fr","guild_de")
write.csv2(gui4,file="data/guilde_08.07.2019.csv")



library(ggplot2)
library(plyr)
library(arm)
library(reshape2)
library(qgraph)
library(vegan)

# Comparaison des guildes
BB <- gui3[,name,by=milieu_fr]
BBB <- dcast(BB,milieu_fr~name)
dis <- as.matrix(vegdist(decostand(BBB[,-1],"norm"), "euclid",binary=T,diag=T, upper=T))

#dis[lower.tri(dis, diag = FALSE)] <- NA
GROUP <- data.frame(id=letters[1:18],name=as.character(as.factor(BBB[,1])))
GROUP$leg <- paste(GROUP$id,GROUP$name,sep=". ")
dist_mi <- 1/dis # one over, as qgraph takes similarity matrices as input
rownames(dist_mi) <- GROUP[,1]
colnames(dist_mi) <- GROUP[,1]

# Flux graph 
png('overlap_guilds.png', width=50, height=20, res=500,unit='cm')
qgraph(dist_mi, layout='spring', vsize=3,esize=3,groups=GROUP[,3],legend = TRUE, borders = FALSE)
dev.off()

# Add bassin versant
hydro <- st_read("Z:/PROJET Ecological Infrastructure/data/sig_data/decoupage_hydrographiquedelasuisse/Hydrografische+Gliederung/Hydrografische Gliederung_LV03/basis04.shp")
names(hydro)
print(hydro[1:5], n = 3)
#plot(hydro[5])
bassin <- st_geometry(hydro)
plot(bassin)

df <- st_sf(id = 1:length(st_geometry(hydro)), area_km2 = st_area(hydro)/10^6 , geometry = st_geometry(hydro))
class(df)

# Convert obs into sf object
library(data.table)
summary(obs)
obs[,carac := as.numeric(0)]
obs <- within(obs,carac[characteristic==1] <- 1)
obs <- within(obs,carac[characteristic==0 & dominant==1] <- 1)

obs2 <- obs %>% filter(!is.na(x_swiss) & milieu_fr=="prairies humides eutrophes" & carac==1)
vars <- c("x_swiss", "y_swiss")

obs_sf <- obs2 %>% mutate_at(vars(x_swiss, y_swiss), as.numeric) %>% 
  st_as_sf(coords = c("x_swiss", "y_swiss")  )

st_crs(obs_sf) <- st_crs(hydro)

# find points within polygons
plants_in_bassin <- st_join(obs_sf, df, join = st_within)

# count plant per bassin
plant_count <- count(as_tibble(plants_in_bassin), id) %>%
  print()

plant_density <- left_join(df, plant_count) %>%
  mutate(density = as.numeric(n / area_km2))  %>%
  print()

library(tmap)
tmap_mode("view")
A <- tm_shape(plant_density) +
  tm_fill(
    col = "density",
    palette = "Greens",
    style = "cont",
    contrast = c(0.1, 1),
    title = "Densité (obs./km2) de plantes indicatrices",
    id = "id",
    showNA = FALSE,
    alpha = 0.8,
    popup.vars = c(
      "Total obs" = "n",
      "Obs./km2" = "density"
    ),
    popup.format = list(
      n = list(format = "f", digits = 0),
      density = list(format = "f", digits = 1)
    )) +
  tm_borders(col = "gray", lwd = 0.2) +
  tm_layout(title= "Prairies humides eutrophes")

tmap_save(A,filename="maps/Prairies humides eutrophes.png",width=20,height=18,unit="cm")


#create grid of 1 ha
temp <- plant_density %>%  filter(density==max(density,na.rm=T))
grid <- temp %>% st_make_grid(cellsize = 100, what = "polygons") %>% st_intersection(temp)

grid_sf <- st_sf(id = 1:length(st_geometry(grid)),geometry = st_geometry(grid))
class(grid_sf)
print(grid_sf, n = 3)
st_crs(grid_sf)

# find points within polygons
in_grid <- st_join(obs_sf, grid_sf, join = st_within)
setDT(in_grid)[id==3712,table(taxon_id)]

# Shannon index by grid cell
in_grid <- in_grid %>% group_by(id) %>% mutate(div=diversity(as.numeric(table(taxon_id))))
AA <- as.data.table(in_grid)
head(AA[!is.na(id),])

plant2 <- left_join(grid_sf, as.data.frame(in_grid %>% select(id,div))[1:2],by="id") %>%
  print()
# # count plant per bassin
# plant_count <- count(as_tibble(plants_in_grid), id) %>%
#   print()
# 
abondance <- data.frame(a=c(200,0,0,0,0,0,151,0,0,0,0,0),b=c(60,40,20,5,10,2,30,5,15,2,12,150))
diversity(t(abondance))

tmap_mode("view")
tm_shape(plant2) +
  tm_fill(
    col = "div",
    palette = "Greens",
    style = "cont",
    contrast = c(0.1, 1),
    title = "Shannon index",
    id = "id",
    showNA = FALSE,
    alpha = 0.8,
    popup.vars = c(
      "Total obs" = "div"
    ),
    popup.format = list(
      div = list(format = "f", digits = 0)
    )) +
  tm_borders(col = "gray", lwd = 0.2) +
  tm_layout(title= "Prairies humides eutrophes")


# Différence (différence d'abondance avant/après 2000)
obs1900 <- obs[!is.na(x_swiss) & milieu_fr=="prairies humides eutrophes" & obs_year<2000]
obs2000 <- obs[!is.na(x_swiss) & milieu_fr=="prairies humides eutrophes" & obs_year>=2000 ]

vars <- c("x_swiss", "y_swiss")

obs1900_sf <- obs1900 %>% mutate_at(vars(x_swiss, y_swiss), as.numeric) %>% 
  st_as_sf(coords = c("x_swiss", "y_swiss")  )
st_crs(obs1900_sf) <- st_crs(obs_sf)
in1900 <- st_join(obs1900_sf, grid_sf, join = st_within)
count1900 <- count(as_tibble(in1900), id) %>%
  print()

obs2000_sf <- obs2000 %>% mutate_at(vars(x_swiss, y_swiss), as.numeric) %>% 
  st_as_sf(coords = c("x_swiss", "y_swiss")  )
st_crs(obs2000_sf) <- st_crs(obs_sf)
in2000 <- st_join(obs2000_sf, grid_sf, join = st_within)
count2000 <- count(as_tibble(in2000), id) %>%
  print()

AA <- left_join(count1900,count2000,by = c("id")) %>% mutate(diff=n.y-n.x) %>% filter(!is.na(id))
AA %>% summarise(min=min(diff,na.rm=T),max=max(diff,na.rm=T),mean=mean(diff,na.rm=T))

tmap_mode("plot")
A <- tm_shape(plant2) +
  tm_fill(
    col = "n",
    palette = "Greens",
    style = "cont",
    contrast = c(0.1, 1),
    title = "Densité de plantes indicatrices (obs/ha) ",
    id = "id",
    showNA = FALSE,
    alpha = 0.8,
    popup.vars = c(
      "Total obs" = "n"
    ),
    popup.format = list(
      n = list(format = "f", digits = 0)
    )) +
  tm_borders(col = "gray", lwd = 0.2) +
  tm_layout(title= "Prairies humides eutrophes")



