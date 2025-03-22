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

bbh=subset(tab,Y>0) %>% group_by(site,block,min_transect_elev,habitat) %>% summarise(culmen_moy=mean(culmen_length,na.rm=TRUE),culmen_sde=sd(culmen_length,na.rm=TRUE)/sqrt(sum(!is.na(culmen_length))),rich=length(unique(hummingbird_species)))

bidon=unique(subset(tab,Y>0,slect=c("site","block","min_transect_elev","habitat","hummingbird_species","culmen_length","month","year")))

#richness
modelr=glmmTMB(rich~habitat*poly(min_transect_elev,2)+(1|block),data=bbh,family=poisson)
br=ggpredict(modelr,c("min_transect_elev[1020:3249]","habitat"))
Anova(modelr)

pl1=ggplot()+geom_ribbon(data=br,aes(x=x,y=predicted,color=group,fill=group,ymin=conf.low,ymax=conf.high),alpha=0.2,color=NA)+geom_line(data=br,aes(x=x,y=predicted,color=group,fill=group))+
geom_point(data=bbh,aes(x=min_transect_elev, y=rich,color=habitat))+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position=c(0.75,0.95))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
xlab("Elevation (m.a.s.l)")+ylab("Hummingbird species richness")+ggtitle("a")+
labs(color="",fill="")

#bill length
modelb=glmmTMB(culmen_length~habitat*poly(min_transect_elev,2)+(1|block),data=bidon,family=nbinom1)
b=ggpredict(modelb,c("min_transect_elev[1020:3249]","habitat"))
Anova(modelb)

pl2=ggplot()+geom_ribbon(data=b,aes(x=x,y=predicted,color=group,fill=group,ymin=conf.low,ymax=conf.high),alpha=0.2,color=NA)+geom_line(data=b,aes(x=x,y=predicted,color=group,fill=group))+
geom_pointrange(data=bbh,aes(x=min_transect_elev, y=culmen_moy,ymin=culmen_moy-1.96*culmen_sde,ymax=culmen_moy+1.96*culmen_sde,color=habitat))+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
xlab("Elevation (m.a.s.l)")+ylab("Hummingbird bill length (mm)")+ggtitle("b")+
labs(color="",fill="")

############################################ INTERACTION FREQUENCY
tab$matching=abs(tab$Tubelength/tab$culmen_length-1) #trait matching
tab$barrier=0
tab$barrier[tab$Tubelength>2*tab$culmen_length]=1
tab$barrier2=as.factor(tab$barrier)
# Scale elevation values for numerical stability
tab$elev=scale(tab$min_transect_elev)
tab$Y2=tab$Y-tab$Y_rob
tab$log_duration=log(tab$duration_sampling_hours)

#### Handle Missing Flower Abundance Data:
# For missing or zero flower abundance values, set to a small default value (1)
tab$abond_flower[is.na(tab$abond_flower) | tab$abond_flower==0]=1
# Compute the logarithm of flower abundance for the model
tab$abond_flower_log=log(tab$abond_flower)

#### Data Cleaning:
# Remove rows where mismatch data is missing
tab=subset(tab,!is.na(matching))
unique(tab$site)

pos <- numFactor(tab$midpoint_Longitude,tab$midpoint_Latitude)
tab$group <- factor(rep(1, nrow(tab)))

 # model = glmmTMB( Y2 ~ matching*habitat + offset(log_duration) + (1|plant_species) + (1|block/hummingbird_species)+
 # (habitat|hummingbird_species), data=tab, family=nbinom1, ziformula=~barrier2*habitat + (1|block))
# save(model,file=paste0(project_folder,"data/","big_model.RData"))
load(paste0(project_folder,"data/","big_model.RData"))
summary(model)

b1=ggpredict(model,c("matching[all]","habitat"), condition=c(log_duration=log(12)))

pl3=ggplot()+geom_ribbon(data=b1,aes(x=x,y=predicted,color=group,fill=group,ymin=conf.low,ymax=conf.high),alpha=0.2,color=NA)+geom_line(data=b1,aes(x=x,y=predicted,color=group,fill=group))+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
xlab("Trait mismatch\n(% of the bill length)")+ylab(expression(paste("Interaction frequency ( ",day^-1,")")))+ggtitle("c")+labs(fill="",color="")+
scale_x_continuous(labels = scales::percent,limits=c(0,10),breaks=c(0,5,10))

b2=ggpredict(model,c("barrier2","habitat"),type="zi_prob", condition=c(log_duration=log(12)))
b2$pred=1-b2$predicted
b2$conf.l=1-b2$conf.low
b2$conf.h=1-b2$conf.high

pl4=ggplot(data=b2,aes(x=x,y=pred,color=group,fill=group))+geom_pointrange(aes(ymin=conf.l,ymax=conf.h),position=position_dodge(width=0.2))+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
xlab("Trait barrier")+ylab("Probability of interaction")+ggtitle("d")+labs(fill="",color="")+
scale_x_discrete(labels =c("no barrier","with barrier"))

bottom=plot_grid(pl1,pl2,pl3,pl4,ncol=2,align="hv")

pdf(paste0(project_folder,"/Fig.3.pdf"),width=5,height=6)
plot_grid(pl1,pl2,pl3,pl4,ncol=2,align="hv")
dev.off();

########### random effects
obj=mixedup::extract_random_coefs(model,re="hummingbird_species",component ="cond",type="zi_random")
obj=subset(obj,effect=="habitatdeforested")
obj=merge(obj,unique(subset(tab,Y>0,select=c(hummingbird_species,culmen_length))),by.x="group",by.y="hummingbird_species")

ggplot(data=obj,aes(x=culmen_length,y=value))+geom_pointrange(aes(ymin=lower_2.5,ymax=upper_97.5))+geom_hline(yintercept=0)

ht=fread("C:/Users/Duchenne/Documents/EPHI_data_clean/hummingbird_traits_2024-03-18/Hummingbird_traits_avonet.txt")
obj=merge(obj,ht,by.x="group",by.y="hummingbird_species")

phy=read.tree("C:/Users/Duchenne/Documents/EPHI_data_clean/hummingbird_traits_2024-03-18/hummingbird_phylo.tre")

phy$tip.label=gsub("_"," ",phy$tip.label)
phy$tip.label[phy$tip.label=="Colibri thalassinus"]="Colibri cyanotus"
phy$tip.label[phy$tip.label=="Schistes geoffroyi"]="Schistes albogularis"
phy=drop.tip(phy,which(duplicated(phy$tip.label)))

p4d=comparative.data(phy, obj, group, vcv=TRUE)
biche=p4d$data
biche[,1]=rownames(biche)
names(biche)[1]="id"

p <-  ggtree(p4d$phy)+geom_tiplab(size=2) + xlim_tree(30)
p2=p+geom_facet(data=biche,mapping=aes(x=value,xmin=lower_2.5,xmax=upper_97.5,color=value),geom=geom_pointrange,panel="Sensitivity to deforestation")+scale_color_gradient2()+theme(strip.background=element_blank())

pdf(paste0(project_folder,"Fig.phylo.pdf"),width=7,height=7)
facet_plot(p2,panel="Sensitivity to deforestation",data=biche,geom=geom_vline,mapping=aes(xintercept=0),linetype="dashed")+labs(color="")
dev.off();

mod2 <- pgls(value ~ culmen_length+Mass+Wing.Length, data=p4d, lambda='ML')
mod2 <-phylolm(value ~ culmen_length+Habitat+Tarsus.Length, data=p4d$data,phy=p4d$phy,model ="lambda",boot=100)
model=lm(value ~ culmen_length+Mass+Habitat, data=p4d$data)


########### forbidden links
tab2=unique(subset(tab,Y>0,select=c(hummingbird_species,site,culmen_length)))
#Load plant transect data
transects=fread(paste0(project_folder,"data/","transects_clean.csv"))
transects$habitat=factor(transects$habitat,levels=c("forest","deforested"))
ttr=unique(transects[,c("site","habitat","min_transect_elev","plant_species","Tubelength")])
tabt=merge(tab2,ttr,by=c("site"),allow.cartesian=TRUE)
tabt$matching=abs(tabt$Tubelength/tabt$culmen_length-1) #trait matching
tabt$barrier=0
tabt$barrier[tabt$Tubelength>2*tabt$culmen_length]=1
tabt$barrier2=as.factor(tabt$barrier)
tabt=subset(tabt,!is.na(matching))
tabt$prob=predict(model,newdata=tabt,type="zlink",re.form=~0)

random=ranef(model)$zi$block
names(random)="link_ef_block"
random$block=rownames(random)
random=merge(random,unique(tab[,c("site","block")]),by=c("block"))
tabt=merge(tabt,random,by="site")
tabt$prob2=boot::inv.logit(tabt$prob+tabt$link_ef_block)


b=tabt %>% group_by(site,habitat,min_transect_elev) %>% summarise(nlink=length(hummingbird_species),nb=length(unique(hummingbird_species)),np=length(unique(plant_species)),prop_forbidden=mean(prob2))

ggplot(data=b,aes(x=min_transect_elev,y=prop_forbidden,color=habitat))+geom_point()
boxplot(prop_forbidden~habitat,data=b)