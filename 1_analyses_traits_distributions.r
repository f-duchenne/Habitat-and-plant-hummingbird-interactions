###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("data.table", "dplyr", "lubridate","lme4","glmmTMB","ggplot2","ggpubr","cowplot","caret","quantreg","ggeffects","car") 

# Check if packages are already installed
inst <- pkgs %in% installed.packages()
# Install missing packages
if (any(!inst)) install.packages(pkgs[!inst])
# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)

colo=c("chartreuse3","gold3") #color vector for habitats
colo2=c("gold3","dodgerblue3","chartreuse3") #color vector for habitat comparisons

# Set the working directory for the "here" package
project_folder=("C:/Users/Duchenne/Documents/esteban_chapter3/")

#Load plant transect data
transects=fread(paste0(project_folder,"data/","transects_clean.csv"))
transects$habitat=factor(transects$habitat,levels=c("forest","deforested"))

#### FIGURE 2, EFFECT OF DEFORESTATION ON PLANT COMMUNITIES
bb=transects %>% group_by(site,block,habitat,min_transect_elev) %>% summarise(rich=length(unique(plant_species)),shan=vegan::diversity(abond_plant_moy, index = "shannon"),
simp=vegan::diversity(abond_plant_moy, index = "simpson"),tubemoy=mean(Tubelength,na.rm=TRUE),tubesde=sd(Tubelength,na.rm=TRUE)/sqrt(sum(!is.na(Tubelength))),flower_abond=sum(abond_flower,na.rm=TRUE),
plant_abond=sum(abond_plant_moy,na.rm=TRUE))

modelr=glmmTMB(rich~habitat*poly(min_transect_elev,2)+(1|block),data=bb,family=nbinom1)
br=ggpredict(modelr,c("min_transect_elev[1020:3249]","habitat"))
Anova(modelr)

pl1=ggplot()+geom_ribbon(data=br,aes(x=x,y=predicted,color=group,fill=group,ymin=conf.low,ymax=conf.high),alpha=0.2,color=NA)+geom_line(data=br,aes(x=x,y=predicted,color=group,fill=group))+
geom_point(data=bb,aes(x=min_transect_elev,y=rich,color=habitat))+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position=c(0.75,0.95))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
xlab("Elevation (m.a.s.l)")+ylab("Plant species richness")+ggtitle("a")+labs(fill="",color="")

#### abundance
modela=glmmTMB(plant_abond~habitat*poly(min_transect_elev,2)+(1|block),data=bb,family=nbinom1)
ba=ggpredict(modela,c("min_transect_elev[1020:3249]","habitat"))
Anova(modela)

pl2=ggplot()+geom_ribbon(data=ba,aes(x=x,y=predicted,color=group,fill=group,ymin=conf.low,ymax=conf.high),alpha=0.2,color=NA)+geom_line(data=ba,aes(x=x,y=predicted,color=group,fill=group))+
geom_point(data=bb,aes(x=min_transect_elev, y=plant_abond,color=habitat))+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
xlab("Elevation (m.a.s.l)")+ylab("Flowering plant abundance")+ggtitle("b")

#### corolla
modelt=glmmTMB(Tubelength~habitat*poly(min_transect_elev,2)+(1|block),data=transects,family=nbinom1)
bt=ggpredict(modelt,c("min_transect_elev[1020:3249]","habitat"))
Anova(modelt)

pl3=ggplot()+geom_ribbon(data=bt,aes(x=x,y=predicted,color=group,fill=group,ymin=conf.low,ymax=conf.high),alpha=0.2,color=NA)+geom_line(data=bt,aes(x=x,y=predicted,color=group,fill=group))+
geom_pointrange(data=bb,aes(x=min_transect_elev, y=tubemoy,ymin=tubemoy-1.96*tubesde,ymax=tubemoy+1.96*tubesde,color=habitat))+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
xlab("Elevation (m.a.s.l)")+ylab("Plant flower corolla length (mm)")+ggtitle("c")

pdf(paste0(project_folder,"/Fig.2.pdf"),width=9,height=4)
plot_grid(pl1,pl2,pl3,ncol=3,align="hv")
dev.off();


##########
bb=bb[order(bb$habitat,bb$block),]

pl2=ggpaired(bb, x = "habitat", y = "plant_abond",
         color = "habitat", line.color = "gray", line.size = 0.4)+
  stat_compare_means(paired = TRUE)+scale_color_manual(values=colo)+
  xlab("Habitat")+ylab("Floewring plant abundance")+ggtitle("b")+
  theme(plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none")

### FUNCTION RAO
transects2=subset(transects,!is.na(plant_species) & !is.na(Tubelength))
trait=as.data.frame(unique(transects2[,c("plant_species","Tubelength")]))
comm=as.data.frame(dcast(transects2,site~plant_species,value.var="abond_plant_moy",fun.aggregate=mean,na.rm=TRUE,fill=0))
trait=subset(trait,plant_species %in% names(comm))

rownames(comm)=comm[,1]
comm=as.matrix(comm[,-1])
trait=trait[trait$plant_species %in% names(apply(comm,2,sum)>0),]
trait2=data.frame(Tubelength=trait$Tubelength)
rownames(trait2)=trait$plant_species
comm=comm[,apply(comm,2,sum)>0]
obj=fd_raoq(trait2,comm)

obj=merge(obj,unique(transects[,c("site","habitat","min_transect_elev","block")]),by="site")
obj=obj[order(obj$habitat,obj$block),]

pl3=ggpaired(obj, x = "habitat", y = "Q",
         color = "habitat", line.color = "gray", line.size = 0.4)+
  stat_compare_means(paired = TRUE)+scale_color_manual(values=colo,label="p.format")+
  xlab("Habitat")+ylab("Functional diversity")+ggtitle("c")+
  theme(plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none")


########################################### FIGURE s1
mat=dcast(transects,plant_species~site,value.var="abond_plant",fun.aggregate=sum,na.rm=TRUE)

lili= as.data.frame(t(combn(unique(transects$site),m=2)))
names(lili)=c("Site.x","Site.y")
lili=merge(lili,unique(transects[,c("site","min_transect_elev","habitat","block")]),by.x="Site.x",by.y="site")
lili=merge(lili,unique(transects[,c("site","min_transect_elev","habitat","block")]),by.x="Site.y",by.y="site")
lili$betad=NA
for(i in 1:nrow(lili)){
mat$X=mat[,lili$Site.x[i],with=F]
mat$Y=mat[,lili$Site.y[i],with=F]
bidon=subset(mat,X>0 | Y>0)
lili$betad[i]=vegan::vegdist(t(bidon[,c("X","Y")]),method="jaccard")
}

lili$habitat_comp=paste0(lili$habitat.x,"-",lili$habitat.y)
lili$habitat_comp[lili$habitat_comp=="deforested-forest"]="forest-deforested"
lili$elev_dist=sqrt((lili$min_transect_elev.x-lili$min_transect_elev.y)^2)
lili$block_comp="across blocks"
lili$block_comp[lili$block.x==lili$block.y]="within blocks"

pl3=ggplot(data=lili,aes(x=elev_dist,y=betad,color=habitat_comp,shape=block_comp))+geom_point()+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position=c(0.8, 0.4))+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+
xlab("Distance in elevation between sites (m.a.s.l)")+ylab("Pairwise beta-diversity between sites")+ggtitle("b")+labs(shape="",color="")

pl4=ggplot(data=lili,aes(x=habitat_comp,y=betad,color=habitat_comp))+geom_boxplot()+
theme_bw()+
theme(axis.line.x = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),axis.ticks.y=element_blank(),
axis.title=element_blank(),axis.text=element_blank(),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+ylab("Pairwise beta-diversity between sites")+ggtitle("")

bottom=plot_grid(pl3,pl4,ncol=2,rel_widths=c(4,1),align="hv")

pdf(paste0(project_folder,"/Fig.S1.pdf"),width=6,height=4)
bottom
dev.off();

