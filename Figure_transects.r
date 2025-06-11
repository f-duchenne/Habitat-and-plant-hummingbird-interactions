###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("data.table", "dplyr", "lubridate","lme4","glmmTMB","ggeffects","car","mixedup","ggplot2","cowplot","gridExtra","ape","caper","ggtree","bipartite","ggsflabel") 

# Check if packages are already installed
inst <- pkgs %in% installed.packages()
# Install missing packages
if (any(!inst)) install.packages(pkgs[!inst])
# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)

colo=c("chartreuse3","gold3") #color vector for habitats

# Set the working directory for the "here" package
project_folder=("C:/Users/Duchenne/Documents/esteban_chapter3/")

tab=fread(paste0(project_folder,"data/","data_for_analyses.csv"))

source("C:/Users/Duchenne/Documents/evolution_pheno_morpho/data_zenodo/scripts/functions and secondary scripts/toolbox.R")

#Load interaction data
tab=fread(paste0(project_folder,"data/","data_for_analyses.csv"))
tab$habitat=factor(tab$habitat,levels=c("forest","deforested"))
tab$Y_freq=tab$Y/tab$duration_sampling_hours 

trans=vect(paste0(project_folder,"data/","GPS_tracks_transects.gpkg"))
trans=project(trans,crs(r))
#plet(r) |> lines(trans)

marge=0.3

sat_reg=read_osm(extent(as.vector(ext(trans))+c(-3,3,-3,3)), type = "bing", zoom = 7)
region=as(sat_reg, "Raster")
df_reg <- as.data.frame(region, xy = T)
boxe=st_as_sfc(st_bbox(as.vector(ext(trans))))
st_crs(boxe) = 4326
boxe=st_transform(boxe,st_crs(sat_reg))
tabi=st_as_sf(data.frame(x=-78.4945,y=-0.198207,label="Quito"),coords=c("x","y"))
st_crs(tabi) = 4326
tabi=st_transform(tabi,st_crs(sat_reg))

pl1=ggplot() +
  geom_raster(data=df_reg,aes(x = x, y = y, fill = rgb(red, green, blue, maxColorValue = 255))) +
  scale_fill_identity() +
  geom_sf(data=boxe,fill="red",alpha=0.5)+
  geom_sf(data=tabi,size=3,fill="black",shape=21,color="white")+
  geom_sf_text_repel(data=tabi,aes(label=label),size=2.5,fill=NA,color="white",
  nudge_x = 5, nudge_y = -5)+
  scale_color_manual(values=colo)+
   scale_y_continuous(n.breaks = 4)+
  scale_x_continuous(n.breaks = 4)+
  coord_sf(expand=FALSE)+
  theme_minimal()+
  theme(panel.grid.major=element_blank(),axis.title=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
  annotation_north_arrow(location = "bl", pad_x    = unit(marge, "cm"),pad_y    = unit(marge, "cm"), style    = north_arrow_fancy_orienteering(fill= c("white", "white"), line_col  = NA, text_col  = "white", text_size = 12)) +
  # Scalebar
  annotation_scale(location = "br", pad_x    = unit(marge, "cm"), pad_y    = unit(marge, "cm"),text_col = "white")+ggtitle("a")       

trans_fig=trans[grep("Alaspungo",trans$site),]
sat=read_osm(extent(as.vector(ext(trans_fig))+c(-0.01,0.01,-0.01,0.01)), type = "bing", zoom = 17)
sat <- as(sat, "Raster")
df <- as.data.frame(sat, xy = T)
head(df)
trans_fig=project(trans_fig,crs(sat))
trans_fig$habitat=factor(c("forest","deforested"),levels=c("forest","deforested"))
ybreaks <- seq(min(df$y), max(df$y), length.out = 4)
pl2=ggplot() +
  geom_raster(data=df,aes(x = x, y = y, fill = rgb(red, green, blue, maxColorValue = 255))) +
  scale_fill_identity() +
  geom_sf(data=st_as_sf(trans_fig),aes(col=habitat),size=1.2)+
  scale_color_manual(values=colo)+
  scale_y_continuous(n.breaks = 3)+
  scale_x_continuous(n.breaks = 3)+
  coord_sf(expand=FALSE)+
  theme_minimal()+
  theme(panel.grid.major=element_blank(),axis.title=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+
  annotation_north_arrow(location = "bl", pad_x    = unit(marge, "cm"),pad_y    = unit(marge, "cm"), style    = north_arrow_fancy_orienteering(fill= c("white", "white"), line_col  = NA, text_col  = "white", text_size = 12)) +
  # Scalebar
  annotation_scale(location = "br", pad_x    = unit(marge, "cm"), pad_y    = unit(marge, "cm"),text_col = "white")+ggtitle("b")       


pdf(paste0(project_folder,"Figure_1.pdf"),width=9,height=5)
grid.arrange(pl1,pl2,ncol=2,widths=c(1,1.6))
dev.off();