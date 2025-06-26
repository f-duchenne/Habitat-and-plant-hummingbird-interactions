###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("data.table", "dplyr", "lubridate","lme4","glmmTMB","ggeffects","car","mixedup","ggplot2","cowplot","gridExtra","ape","caper","ggtree","bipartite") 

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

image <- file.path(paste0(project_folder,"data/essai.PNG"))

# r1=terra::rast(paste0(project_folder,"data/","ETH_GlobalCanopyHeight_10m_2020_S03W081_Map.tif"))
# r2=terra::rast(paste0(project_folder,"data/","ETH_GlobalCanopyHeight_10m_2020_N00W081_Map.tif"))
# r=merge(r1,r2)
# writeRaster(r,"canoyp_height.tif")
r=rast(paste0(project_folder,"data/","canoyp_height.tif"))

trans=vect(paste0(project_folder,"data/","GPS_tracks_transects.gpkg"))
trans=project(trans,crs(r))
#plet(r) |> lines(trans)
cover_trans=extract(r,trans,weights=TRUE,exact=TRUE) 
cover_trans=merge(cover_trans,trans[,c("id","site")],by.x="ID",by.y="id")
trans=cover_trans %>% group_by(site) %>% summarise(TC=sum(ETH_GlobalCanopyHeight_10m_2020_S03W081_Map*weight,na.rm=TRUE)/sum(weight,na.rm=TRUE))
trans=merge(trans,unique(tab[,c("site","habitat","block","min_transect_elev","midpoint_Latitude","midpoint_Longitude")]),by="site")
boxplot(TC~habitat,data=trans)
trans=trans %>% group_by(block) %>% mutate(TC_ref=TC[habitat=="forest"])
trans$perturb=trans$TC_ref-trans$TC


matg=dcast(tab,site+plant_species~hummingbird_species,value.var="Y_freq",fun.aggregate=mean,na.rm=TRUE,fill=0)

obj=split(matg,f=matg$site)

resf=NULL
for(i in 1:length(obj)){
bidon=obj[[i]]
site=names(obj)[i]
mat=as.matrix(bidon[,-c(1:2)])
mat=mat[apply(mat,1,sum)>0,apply(mat,2,sum)>0]

exti=second.extinct(mat, participant = "lower", method = "random", nrep = 10)
res=data.frame(site=site,robust=robustness(exti))
nulos=vaznull(mat*100,N=10)
vec=sapply(lapply(nulos,second.extinct,participant="lower",method="random",nrep=10),robustness)
res$mean_nul=mean(vec)
res$sd_nul=sd(vec)
res$np=nrow(mat)
res$nh=ncol(mat)


for(rho in c(0.0,0.05,0.1,0.15,0.2)){
bidon=res
comp_p=matrix(rho,nrow(mat),nrow(mat))
diag(comp_p)=1
comp_h=matrix(rho,ncol(mat),ncol(mat))
diag(comp_h)=1
A=rbind(cbind(-1*comp_p,mat),cbind(t(mat),-1*comp_h))

bidon$rho=rho
bidon$stab=Omega(A)
resf=rbind(resf,bidon)
}

}


resf$zscore=(resf$robust-resf$mean_nul)/resf$sd_nul
resf=merge(resf,trans,by="site")

# resf=resf %>% group_by(block) %>% mutate(robust_ref=robust[habitat=="forest"],stab_ref=stab[habitat=="forest"])
# resf$robust_delta=resf$robust-resf$robust_ref
# resf$stab_delta=resf$stab-resf$stab_ref

pl1=ggplot(data=resf,aes(color=habitat,y=stab,x=as.factor(rho)))+geom_boxplot()+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+xlab("Strength of interspecific competition")+ylab("Structural stability")+
scale_color_manual(values=colo)+ggtitle("a")

resf$dive=resf$np+resf$nh

rhos=0.0

model=glm(stab~dive*habitat,data=subset(resf,rho==rhos),family="quasibinomial")
car::Anova(model)
b=ggpredict(model,c("dive","habitat"))


pl2=ggplot()+geom_point(data=subset(resf,rho==rhos),aes(color=habitat,y=stab,x=dive))+geom_ribbon(data=b,aes(x=x,ymin=conf.low,ymax=conf.high,y=predicted,fill=group),color=NA,alpha=0.3)+
geom_line(data=b,aes(x=x,y=predicted,color=group))+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="right")+xlab("Species richness (plant+hummingbird)")+ylab("Structural stability")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ggtitle("b")+labs(color="",fill="")



pdf(paste0(project_folder,"Figure_4.pdf"),width=7,height=3)
plot_grid(pl1,pl2,align="h")
dev.off();






unique(tab[,c("site","habitat","block","min_transect_elev","midpoint_Latitude","midpoint_Longitude")])