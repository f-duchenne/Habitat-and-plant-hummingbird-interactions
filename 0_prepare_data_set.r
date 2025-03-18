###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("data.table", "dplyr", "lubridate","lme4") 

# Check if packages are already installed
inst <- pkgs %in% installed.packages()
# Install missing packages
if (any(!inst)) install.packages(pkgs[!inst])
# Load packages
pkg_out <- lapply(pkgs, require, character.only = TRUE)

# Set the working directory for the "here" package
project_folder=("C:/Users/Duchenne/Documents/esteban_chapter3/")


#LOAD HUMMINGBIRD DATA
dat=fread(paste0(project_folder,"data/","Interactions_data_Ecuador.txt"))
dat[dat==""]=NA
dat$date=as.IDate(dat$date,"%Y/%m/%d") #be sure date columns is recognize as date
dat$piercing[is.na(dat$piercing)]="no"
#LOAD CAMERA INFORMATION
cameras=fread(paste0(project_folder,"data/","Cameras_data_Ecuador.txt"),na.strings = c("",NA))
cameras$end_date=as.IDate(cameras$end_date,"%Y/%m/%d") #be sure date columns is recognize as date
cameras$start_date=as.IDate(cameras$start_date,"%Y/%m/%d") #be sure date columns is recognize as date
cameras$month=month(cameras$start_date) #extract month from date column
cameras$year=year(cameras$start_date) #extract year from date column
cameras=subset(cameras,year<2023) # keep only EPHI data
plant_for_res=unique(cameras[,c("plant_species","month","year","site")])
#LOAD TRANSECT DATA
transects=fread(paste0(project_folder,"data/","Transect_data_Ecuador.txt"),na.strings = c("",NA))
transects$date=as.IDate(transects$date,"%Y/%m/%d") #be sure date columns is recognize as date
transects$month=month(transects$date) #extract month from date column
transects$year=year(transects$date) #extract year from date column
#calculate a total abundance per plant species per month per site
transects=subset(transects,!is.na(plant_species)) %>% group_by(site,year,month,plant_species) %>% summarise(abond_flower=sum(total_flowers,na.rm=T),
																								abond_plant=length(total_flowers)) 
# Replace zero abundance with a minimum value
transects$abond_flower[transects$abond_flower==0]=1
#crate ID date column for each site and month
transects$combi_dates=apply(transects[,c("site","month","year")],1,paste,collapse="-")
plant_for_res$combi_dates=apply(plant_for_res[,c("site","month","year")],1,paste,collapse="-")
# labell data and replace NA abundance with a minimum value
plant_for_res$duringEPHI="no"
plant_for_res$duringEPHI[plant_for_res$combi_dates %in% transects$combi_dates]="yes"
transects=merge(transects,plant_for_res,by=c("plant_species","month","year","site","combi_dates"),all=T)
transects$abond_flower[is.na(transects$abond_flower) & transects$duringEPHI=="yes"]=1
transects$abond_plant[is.na(transects$abond_plant) & transects$duringEPHI=="yes"]=1

#LOAD SITE METADATA:
sites=fread(paste0(project_folder,"data/","Site_metadata_Ecuador.txt"),na.strings = c("",NA))

#MERGE THEM:
#interactions and cameras:

dim(dat)
dat=merge(dat,cameras,by=c("waypoint","site"),all.x=F,all.y=T) #we want to keep cameras that did not detect any hummingbird
dim(dat)
unique(subset(dat,site=="LasGralarias" & year==2019)$month)
dat=subset(dat,!is.na(plant_species)) #remove data with no plant ID
dim(dat)
unique(subset(dat,site=="LasGralarias" & year==2019)$month)
unique(subset(dat,site=="LasGralarias" & year==2019)$month)
dim(subset(dat,is.na(hummingbird_species)))
#add transect data:
dim(dat)
dat=merge(dat,transects,by=c("site","month","year","plant_species"),all.x=T,all.y=F)
dim(dat)
#add site data:
dim(dat)
dat=merge(dat,sites,by=c("site"),all=FALSE)
dim(dat)
dim(dat)
transects=merge(transects,sites,by=c("site"),all=FALSE)
dim(dat)

dat=subset(dat,year(start_date)<=2022)

# Filter out non paired forest sites
tab=tab[!(tab$site %in% c("Yanacocha", "SantaLuciaUpper", "SantaLuciaLower", "Mashpi_Capuchin")),]
transects=transects[!(transects$site %in% c("Yanacocha", "SantaLuciaUpper", "SantaLuciaLower", "Mashpi_Capuchin")),]

# Add block ID variable with 7 levels, grouping pairs 
# of forest and deforested transects
dat <- dat %>% 
  mutate(block = case_when(
    endsWith(site, "decocha") ~ "Block7",
    endsWith(site, "nacocha_disturbed") ~ "Block7",
    endsWith(site, "pungo") ~ "Block6",
    endsWith(site, "pungo_disturbed") ~ "Block6",
    endsWith(site, "alarias") ~ "Block5",
    endsWith(site, "unapi") ~ "Block5",
    endsWith(site, "tamia") ~ "Block4",
    endsWith(site, "tamia_disturbed") ~ "Block4",
    endsWith(site, "pucuna") ~ "Block3",
    endsWith(site, "misitana") ~ "Block3",
    endsWith(site, "Laguna") ~ "Block2",
    endsWith(site, "agusa") ~ "Block2",
    endsWith(site, "Choco") ~ "Block1",
    endsWith(site, "Choco_disturbed") ~ "Block1",
  ))
transects <- transects %>% 
mutate(block = case_when(
endsWith(site, "decocha") ~ "Block7",
endsWith(site, "nacocha_disturbed") ~ "Block7",
endsWith(site, "pungo") ~ "Block6",
endsWith(site, "pungo_disturbed") ~ "Block6",
endsWith(site, "alarias") ~ "Block5",
endsWith(site, "unapi") ~ "Block5",
endsWith(site, "tamia") ~ "Block4",
endsWith(site, "tamia_disturbed") ~ "Block4",
endsWith(site, "pucuna") ~ "Block3",
endsWith(site, "misitana") ~ "Block3",
endsWith(site, "Laguna") ~ "Block2",
endsWith(site, "agusa") ~ "Block2",
endsWith(site, "Choco") ~ "Block1",
endsWith(site, "Choco_disturbed") ~ "Block1",
))

#Duration from picture when non-available from cameras
dat$duration_sampling_hours[is.na(dat$duration_sampling_hours)]=dat$duration_from_pics[is.na(dat$duration_sampling_hours)]
dat=subset(dat,duration_sampling_hours>=5 & camera_problem!="yes")

#### CALCULATE INTERACTION FREQUENCIES
tab1=dat %>% dplyr::group_by(hummingbird_species,plant_species,waypoint,duration_sampling_hours,site,block,month,year,abond_flower,habitat,
midpoint_Longitude,midpoint_Latitude,min_transect_elev,Country) %>% dplyr::summarise(Y_rob=length(time[piercing=="yes"]),Y=length(time)) #number of interaction detected

#### INFER ZEROS using dcast and melt functions
mat=dcast(tab1,site+block+year+month+waypoint+plant_species+duration_sampling_hours+abond_flower+midpoint_Longitude+
midpoint_Latitude+min_transect_elev+habitat+Country~hummingbird_species,value.var="Y",fill=0)
tab=melt(mat,id.vars=c("site","block","year","month","waypoint","plant_species","duration_sampling_hours",
"abond_flower","midpoint_Longitude","midpoint_Latitude","min_transect_elev","Country","habitat"),variable.name="hummingbird_species",value.name="Y")
tab=merge(tab,tab1[,c("plant_species","hummingbird_species","waypoint","site","block","year","month","Y_rob")],by=c("plant_species","hummingbird_species","waypoint","site","block","year","month"),all=TRUE)
tab$Y_rob[is.na(tab$Y_rob)]=0
# CALCULATE TOTAL ABUNDANCE OF HUMMINGBIRDS PER SITES:
tab=tab %>% dplyr::group_by(hummingbird_species,block) %>% dplyr::mutate(abond_total_h=sum(Y)) #number of interaction detected
#REMOVE HUMMINGBIRDS THAT ARE TOTALLY ABSENT FROM ONE SITE:
tab=subset(tab,abond_total_h>0)
subset(tab,block=="Block6" & hummingbird_species=="Boissonneaua flavescens" & Y>0)

#function to calculate the mode of a vector
getmode <- function(v) {
 uniqv <- unique(na.omit(v))
 uniqv[which.max(tabulate(match(v, uniqv)))]
}

#function to calculate weighted mean
wmean <- function(v,z) {
 uniqv <- v[!is.na(v) & !is.na(z)]
 uniqz <- z[!is.na(v) & !is.na(z)]
 if(length(uniqz>0)){return(sum(uniqv*uniqz)/sum(uniqz))}else{return(NA)}
}
	
#LOAD HUMMINGBIRD TRAIT DATA AND COMBINE THEM TO HAVE ONE VALUE PER SPECIES
tr=fread(paste0(project_folder,"data/","Hummingbird_traits.txt"),na.strings = c("",NA))
tr$tail_length=as.numeric(as.character(tr$tail_length))
tr1=tr %>% dplyr::group_by(hummingbird_species,hummingbird_family) %>% dplyr::summarise(nindh=length(culmen_length[!is.na(culmen_length)]),bill_length=wmean(bill_length,N),culmen_length=wmean(culmen_length,N),
tail_length=wmean(tail_length,N))
model=lm(culmen_length~bill_length,data=tr1[!is.na(tr1$bill_length) & !is.na(tr1$culmen_length),])
tr1$culmen_length[is.na(tr1$culmen_length) & !is.na(tr1$bill_length)]=
predict(model,newdata=tr1[is.na(tr1$culmen_length) & !is.na(tr1$bill_length),])

#LOAD PLANT TRAIT DATA AND COMBINE THEM TO HAVE ONE VALUE PER SPECIES
tr=fread(paste0(project_folder,"data/","Plant_traits.txt"),na.strings = c("",NA))
tr2=subset(tr,!is.na(plant_species))  %>% group_by(plant_species,plant_family,plant_genus) %>%
summarise(nindp=length(Tubelength[!is.na(Tubelength)]),Tubelength=mean(Tubelength*10,na.rm=T),StamenLength=mean(StamenLength*10,na.rm=T),
StyleLength=mean(StyleLength*10,na.rm=T),Opening_corolla_lateral=mean(Opening_corolla_lateral,na.rm=T),
CurvatureCentral=mean(CurvatureCentral,na.rm=T),FlowerType=getmode(FlowerType),sexualsys=getmode(SexualSystem))

length(unique(tr2$plant_species[tr2$plant_species %in% dat$plant_species]))
length(unique(dat$plant_species))

#MERGE TRAITS AND DATA
tab=merge(tab,tr1,by=c("hummingbird_species"),all.x=TRUE,all.y=FALSE)
tab=merge(tab,tr2,by=c("plant_species"),all.x=TRUE,all.y=FALSE)

tab=subset(tab,!is.na(hummingbird_family))

bidon=unique(transects[,c("month","year","site","block","habitat")])
bidon_s=subset(bidon,habitat=="deforested")
bidon_s=bidon_s[,c("block","year","month")]
bidon_s=rbind(bidon_s,
data.frame(block="Block2",year=c(2017,2017,2017,2018,2019),month=c(10,11,12,8,4)),
data.frame(block="Block5",year=c(2017,2017,2018,2019,2019,2019),month=c(10,11,9,1,3,5)),data.frame(block="Block1",year=c(2019),month=c(9)))

nrow(transects)
transects=merge(transects,bidon_s,by=c("block","year","month"))
transects %>% group_by(block,site) %>% summarise(length(unique(paste0(month,year))))

bidon=unique(tab[,c("month","year","site","block","habitat")])
bidon_s=subset(bidon,habitat=="deforested")
bidon_s=bidon_s[,c("block","year","month")]
bidon_s=rbind(bidon_s,
data.frame(block="Block2",year=c(2017,2017,2018,2019),month=c(11,12,8,4)),
data.frame(block="Block5",year=c(2017,2017,2018,2019,2019,2019),month=c(10,11,9,1,3,5)),data.frame(block="Block1",year=c(2019,2019),month=c(4,9)))

nrow(tab)
tab=merge(tab,bidon_s,by=c("block","year","month"))
tab %>% group_by(block,site) %>% summarise(length(unique(paste0(month,year))))

fwrite(tab,paste0(project_folder,"data/","data_for_analyses.csv"))

#plant list based on transects
plant_res=merge(transects,tr2[,c("plant_species","Tubelength")],by=c("plant_species"),all.x=T)
plant_res=subset(plant_res,!is.na(plant_species))
plant_res=plant_res %>% group_by(site,block) %>% mutate(n_visits=length(unique(paste(month,year))))
plant_res=plant_res %>% group_by(plant_species,site,block,habitat) %>% mutate(abond_flower_moy=sum(abond_flower,na.rm=T)/n_visits,abond_plant_moy=sum(abond_plant,na.rm=T)/n_visits)
fwrite(plant_res,paste0(project_folder,"data/","transects_clean.csv"))