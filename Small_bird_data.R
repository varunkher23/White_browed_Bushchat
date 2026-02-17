library(googlesheets4)
library(tidyverse)
library(googledrive)
library(lubridate)
library(sf)
library(BiodiversityR)
library(Distance)
library(Rdistance)
library(janitor)
library(spOccupancy)
library(geosphere)

get_sp_occupancy <- function(model_pred,sp){
  curr.sp <- which(spp$species == sp)
  curr.sp.psi.samples <- model_pred$psi.0.samples[, curr.sp, ]
  curr.sp.occ <- apply(curr.sp.psi.samples, 2, mean)
  curr.sp.occ.lcl <- apply(curr.sp.psi.samples, 2, quantile, 0.05)
  curr.sp.occ.ucl <- apply(curr.sp.psi.samples, 2, quantile, 0.95)
  curr.sp.occ.ci.width = curr.sp.occ.ucl - curr.sp.occ.lcl
  sp.occ = data.frame(Species = sp,occ_median = curr.sp.occ, occ_lcl = curr.sp.occ.lcl, occ_ucl = curr.sp.occ.ucl, curr.sp.ci.width = curr.sp.occ.ci.width)%>%
    cbind(occ.cov%>%select(-transect_name_2))%>%
    unique.data.frame()
  output = list(sp.occ = sp.occ, sp.psi.samples = curr.sp.psi.samples)
  return(output)
}

grid = st_read("GIS_Layers/small_bird_grids.shp",crs= st_crs("WGS84"))%>%
  rename(transect_name_2 = ID2)

check_grid = function(location){
  loc <- data.frame(lon = unlist(strsplit(location,","))[2], lat = unlist(strsplit(location,","))[1])
  point.sf <- st_as_sf(loc, coords = c("lon", "lat"), crs = st_crs("WGS84"))
  grid_id = st_join(point.sf, grid, join = st_within)$transect_name_2
  return(grid_id)
}

form_data = read_sheet("https://docs.google.com/spreadsheets/d/1EHNQSRqDqvWNox3y6RbbZ5cjY4BBy8qhw2oqTkxBTuE/edit#gid=130138455",
                     sheet = 1)%>%filter(is.na(KEY)==F)%>%select(-`data-meta-instanceID`,-`data-meta-audit`)%>%unique.data.frame()%>%
  full_join(read_sheet("https://docs.google.com/spreadsheets/d/19PtQBv03r-dNKfrdahdmIJwgy_N5TJw2DsoLawYlsCI/edit#gid=0",
                        sheet = 1)%>%filter(is.na(KEY)==F)%>%rename(`data-Details2-timestamp` = `data-Details2-start_time`)%>%
              select(-`data-meta-instanceID`)%>%unique.data.frame())%>%
  mutate(`data-Details-start_location`=ifelse(is.na(`data-Details-start_location`),`data-Details2-auto_geopoint`,`data-Details-start_location`))%>%
                                 rowwise()%>%
  mutate(grid_id = check_grid(`data-Details-start_location`))%>%
  unique.data.frame()

form_data%>%
  filter(`data-Details-Site_Name`=="dnp")%>%
  filter(`data-Details-method`=="line")%>%
  mutate(check = ifelse(grid_id == `data-Details-Transect_Name`,1,0))%>%
  View()

metadata = data.frame()
for (i in 1:nrow(form_data)){
  if (!is.na(form_data$`data-meta-audit`[i])) {
    KEY = form_data$KEY[i]
    drive_download(form_data$`data-meta-audit`[i],path = "Metadata/odk_metadata.csv", overwrite = T)
    a = read.csv("Metadata/odk_metadata.csv")
    b = a %>%
      filter(grepl("Bird_details2",node))%>%
      separate(node, sep="/",c("Prefix1","Prefix2","Bird"))%>%
      select(Bird,latitude,longitude,accuracy)%>%
      mutate(BIRD_KEY = paste(KEY, Bird, sep = "/"))
  } else {
    b = data.frame()
  }
  metadata = rbind(metadata,b)
}

bird_data=read_sheet("https://docs.google.com/spreadsheets/d/1EHNQSRqDqvWNox3y6RbbZ5cjY4BBy8qhw2oqTkxBTuE/edit#gid=130138455",sheet = 2)%>%
  unique.data.frame()%>%rename(BIRD_KEY=KEY)%>%rename(KEY=PARENT_KEY)%>%unique.data.frame()%>%
  full_join(read_sheet("https://docs.google.com/spreadsheets/d/19PtQBv03r-dNKfrdahdmIJwgy_N5TJw2DsoLawYlsCI/edit#gid=728990291",sheet = 2)%>%
  unique.data.frame()%>%unique.data.frame()%>%rename(BIRD_KEY=KEY)%>%rename(KEY=PARENT_KEY)%>%unique.data.frame())%>%
  left_join(metadata%>%filter(latitude>25)%>%filter(accuracy<50))

data=right_join(form_data,bird_data, by = "KEY")%>%unique.data.frame()%>%
  `colnames<-`(str_extract(colnames(.),"\\w+$"))%>%
  select(unique(colnames(.)))%>%
  separate(grid_id,c("Sub_site","Transect_Number","Segment_Number"),sep = "_")%>%
  mutate(Transect_Name=paste(Sub_site,Transect_Number,sep="_"))%>%
  mutate(Transect_Name_2=paste(Transect_Name,Segment_Number,sep = "_"))%>%
  rename(Radial_Distance = Distance)%>%
  mutate(distance=abs(Radial_Distance*sin(abs(Angle)*pi/180)))%>%
  mutate(Year = year(timestamp), Month = month(timestamp))%>%
  mutate(Season = ifelse(Month > 3 & Month < 7,"Summer","Other"))%>%
  mutate(Season = ifelse(Month > 7 & Month < 10,"Monsoon",Season))%>%
  mutate(Season = ifelse(Month > 10 | Month < 3,"Winter",Season))%>%
  mutate(Survey_Season = ifelse(Season != "Winter", paste(Season,Year,sep="_"),paste(Season,paste(Year,Year+1,sep="-"),sep="_")))%>%
  mutate(Survey_Season = ifelse(Season =="Winter" & Month <3, paste(Season,paste(Year-1,Year,sep="-"),sep="_"), Survey_Season))%>%
  rowwise()%>%
  mutate(projected_long = destPoint(c(longitude,latitude),d = distance,b = Angle)[1],
         projected_lat = destPoint(c(longitude,latitude),d = distance,b = Angle)[2])

#write_rds(data, "BirdData_12.09.2025")

#### Basic parameters

effort_summary=data%>%
  filter(Site_Name=="dnp")%>%
  filter(method == "line")%>%
  group_by(Survey_Season,Transect_Name)%>%
  #filter(Species%in%c("Other","Others","None")==F)%>%
  summarise(n_line=length(unique(KEY)))%>%
  ungroup()%>%
  pivot_wider(names_from = Survey_Season,values_from = n_line)%>%
  mutate(total_winter_effort = rowSums(select(.,-Transect_Name),na.rm=T))

sampled_grids = grid%>%filter(transect_name_2%in%unique(data$Transect_Name_2))
plot(sampled_grids$geometry)

basic_parameters=data%>%
  filter(Species%in%c("Other","Others","None")==F)%>%
  filter(method=="line")%>%
  filter(Site_Name=="dnp")%>%
  filter(Transect_Name %in% (effort_summary%>%filter(!is.na(total_winter_effort)))$Transect_Name)%>% # Only for grids repeated every winter
  group_by(Survey_Season)%>%
  summarise(n_birds=sum(Number),
            n_species=n_distinct(Species),
            effort=length(unique(KEY))*0.5,
            max_detection_distance=quantile(Radial_Distance,0.95),
            avg_detection_distance=mean(Radial_Distance,na.rm=T))%>%
  mutate(ER=n_birds/effort)%>%
  mutate(ESW=as.numeric(0.34))%>% #### Based on old data, ESW is 340m for DNP
  mutate(Approx_density=ER/ESW)

WBB_df = data%>%
  clean_names()%>%
  filter(method == "line")%>%
  filter(species == "Saxicola macrorhynchus")%>%
  mutate(Sample.Label = transect_name_2)%>%
  mutate(Region.Label = "Thar")%>%
  mutate(Area = as.numeric(1))%>%
  rename(size = number)%>%
  mutate(Effort = as.numeric(1))%>%
  dplyr::select(Region.Label, Area, Sample.Label, Effort, species, size, distance,projected_lat,projected_long)

plot(WBB_df$latitude,WBB_df$longitude)

effort = data%>%
  filter(Site_Name=="dnp")%>%
  filter(method == "line")%>%
  filter(Transect_Name_2 != "JSM34_31_trial")%>%
  select(Transect_Name,Transect_Name_2,timestamp,Length,Season,Survey_Season)%>%
  unique.data.frame()%>%
  group_by(Transect_Name_2)%>%
  mutate(replicate = row_number())%>%
  ungroup()%>%
  mutate(Length = ifelse(is.na(Length),0.5,Length))%>%
  clean_names()%>%
  left_join(grid%>%as.data.frame()%>%select(transect_name_2, Land_cover))%>%
  filter(survey_season == "Winter_2023-2024")

##### Single species occupancy

occu_data = data%>%
  clean_names()%>%
  mutate(Presence = 1)%>%
  right_join(effort%>%clean_names()%>%filter(season == "Winter"))%>%
  unique.data.frame()%>%
  select(-land_cover)

sp = "Saxicola macrorhynchus"

sp_occu_data = occu_data%>%
  filter(species == sp)%>%
  right_join(effort)%>%
  group_by(transect_name_2,replicate)%>%
  summarise(Presence = ifelse(sum(number,na.rm = T)>0,1,0))%>%
  pivot_wider(id_cols = transect_name_2, values_from = Presence, names_from = replicate)%>%
  ungroup()

occ.cov = left_join(select(sp_occu_data,transect_name_2), as.data.frame(grid))%>%select(transect_name_2,Land_cover)%>%
  mutate(Presence = 1)%>%
  pivot_wider(id_cols = transect_name_2, names_from = Land_cover, values_from = Presence,values_fill = 0)%>%
  select(-transect_name_2)

y = as.matrix(sp_occu_data[,-1])%>%
  `rownames<-`(sp_occu_data$transect_name_2)%>%
  `colnames<-`(colnames(sp_occu_data)[-1])

sp_data = list(y= y, occ.covs = occ.cov)

occ.formula <- ~ CL + RL + RG
det.formula <- ~ 1

occ.formula.null <- ~ 1
det.formula.null <- ~ 1

inits <- list(alpha = 0, 
              beta = 0, 
              z = apply(sp_data$y, 1, max, na.rm = TRUE))
priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
               beta.normal = list(mean = 0, var = 2.72))

n.samples <- 10000
n.burn <- 3000
n.thin <- 2
n.chains <- 3

out.null <- PGOcc(occ.formula = occ.formula.null, 
                  det.formula = det.formula.null, 
                  data = sp_data, 
                  inits = inits, 
                  n.samples = n.samples, 
                  priors = priors, 
                  n.omp.threads = n.chains, 
                  verbose = TRUE, 
                  n.report = 1000, 
                  n.burn = n.burn, 
                  n.thin = n.thin, 
                  n.chains = n.chains, 
                  k.fold = 4, k.fold.threads = 4)

out <- PGOcc(occ.formula = occ.formula, 
             det.formula = det.formula, 
             data = sp_data, 
             inits = inits, 
             n.samples = n.samples, 
             priors = priors, 
             n.omp.threads = n.chains, 
             verbose = TRUE, 
             n.report = 1000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains,
             k.fold = 4, k.fold.threads = 4)

X.0 = as.matrix(data.frame(1,occ.cov))
out.pred = predict(out, X.0)

curr.sp.occ <- apply(out.pred$psi.0.samples, 2, mean)
curr.sp.occ.lcl <- apply(out.pred$psi.0.sample, 2, quantile, 0.05)
curr.sp.occ.ucl <- apply(out.pred$psi.0.sample, 2, quantile, 0.95)
curr.sp.occ.ci.width = curr.sp.occ.ucl - curr.sp.occ.lcl
prediction.occu = data.frame(Species = sp,occ_median = curr.sp.occ, occ_lcl = curr.sp.occ.lcl, occ_ucl = curr.sp.occ.ucl, curr.sp.ci.width = curr.sp.occ.ci.width)%>%
  cbind(occ.cov)%>%unique.data.frame()

boxplot(out.pred$psi.0.samples[,1], out.pred$psi.0.samples[,20], out.pred$psi.0.samples[,4],
        names = c("Croplands","Rangelands","Restored Grasslands"))

#### Density estimation

WBB_er = occu_data%>%
  filter(species == "Saxicola macrorhynchus")%>%
  separate(auto_geopoint,sep=",",c("Latitude","Longitude"))%>%
  mutate(Latitude = as.numeric(Latitude) + cos(animal_bearing * (3.14159 / 180))*(radial_distance/11000), 
         Longitude = as.numeric(Longitude) + sin(animal_bearing * (3.14159 / 180))*(radial_distance/11000))%>%
  mutate(object = row_number())%>%
  right_join(effort%>%clean_names())%>%
  mutate(Sample.Label = paste(transect_name_2,replicate, sep = "_"))%>%
  mutate(Region.Label = "Thar")%>%
  mutate(Area = as.numeric(1))%>%
  rename(size = number)%>%
  rename(Effort = length)%>%
  mutate(land_cover = ifelse(land_cover == "RG", "Restored","Unrestored"))%>%
  dplyr::select(Region.Label, Area, Sample.Label, Effort, species, size, distance, object,land_cover)

hist(WBB_df$distance)  

conversion <- convert_units("meter", "kilometer", "square kilometer")

hn01 <- ds(WBB_df, transect = "line", key = "hn", 
           truncation = list(left=0,right=100), formula = ~1)
plot(hn01)
summary(hn01)

density_estimate = dht2(hn01, flatfile = WBB_er, 
                        convert_units = conversion, strat_formula = ~land_cover)

WBB_er%>%
  group_by(land_cover)%>%
  summarise(n = sum(size, na.rm = T), effort = sum(Effort, na.rm = T))%>%
  ungroup()%>%
  mutate(ER = n / effort)
summary(lm(data = WBB_er,ER~land_cover))

density_estimate_plot = density_estimate%>%
  select(land_cover,Abundance,Abundance_se)%>%
  filter(land_cover != "Total")%>%
  mutate(Abundance_max = Abundance + Abundance_se)

ggplot(density_estimate_plot, aes(x = factor(land_cover, levels = c("Unrestored", "Restored")), y = Abundance, fill = land_cover))+
  geom_errorbar(aes(ymax = Abundance_max, ymin = Abundance), width = 0.2, position = position_dodge(width = 0.9)) +
  geom_bar(stat = "identity", position = "dodge")+
  labs( x = "Type of Grid", 
        y = "Estimated Density (per sq.km)")+
  theme_classic(
    base_size = 24,
    base_family = "sans"
  )+
  theme(legend.position = "none")+
  scale_fill_manual(name = "Type of Area", values = c("#F6BE00","#8C0301"))+
  ylim(0,5)



