###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("data.table", "dplyr", "lubridate","lme4","glmmTMB","ggeffects","car","mixedup","ggplot2","cowplot","gridExtra","ape","caper","ggtree") 

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

#Load interaction data
tab=fread(paste0(project_folder,"data/","data_for_analyses.csv"))
tab$habitat=factor(tab$habitat,levels=c("forest","deforested"))
tab$Y_freq=tab$Y/tab$duration_sampling_hours 

matg=dcast(tab,site+plant_species~hummingbird_species,value.var="Y_freq",fun.aggregate=mean,na.rm=TRUE,fill=0)

obj=split(matg,f=matg$site)

for(i in 1:length(obj)){
bidon=obj[[i]]
site=names(obj)[i]
mat=as.matrix(bidon[,-c(1:2)])
mat[apply(mat,1,sum)>0,apply(mat,2,sum)>0]
exti=second.extinct(mat, participant = "lower", method = "abun", nrep = 10)
data.frame(site=site,robust=robustness(exti)
}
