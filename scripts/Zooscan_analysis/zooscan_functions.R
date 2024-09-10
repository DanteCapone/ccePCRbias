#Zooscan functions for CCE Metabarcoding Paper

#Load required Libraries
librarian::shelf(tidyverse, stringr,here, RColorBrewer)

#Read in Zooscan Data, process and 
readEcotaxa <- function(data) {
  #data: A .tsv file exported from Ecotaxa
  #select relevant columns to calculate data
  
  #selected columns for database upload
  data = data %>% 
    mutate(sample=sample_id,
           Haul =sample_id,
           Region ="California Current",
           Detail_Location =data$sample_id,
           Comment ="",
           process_particle_pixel_size_mm =0.0106,
          #Caclulate volume filtered to get concentration
          sample_conc=acq_sub_part/sample_tot_vol*object_depth_max,
          cycle= str_extract(object_id, "^[^_-]+"))
  
  #Select relevant columns
  data_select =data %>% dplyr::select(., sample_ship, sample_program, sample_id, Haul, Region, Detail_Location, Comment, 
                                                  object_date, object_time, object_lat, object_lon, sample_bottomdepth, object_depth_min,
                                                  object_depth_max, object_annotation_category, object_annotation_hierarchy, object_annotation_person_name,
                                                  sample, object_id, sample_id, sample_tot_vol, acq_sub_part,object_feret, 
                                                  object_area, object_major, object_minor, object_area_exc, process_particle_pixel_size_mm, acq_max_mesh,
                                                  sample_conc,cycle, object_annotation_status,acq_id)
  #Filter: only living and validated
  data_select<-data_select %>% 
    filter(object_annotation_status=="validated")%>%
    filter(!str_detect(object_annotation_hierarchy, regex("not-living", ignore_case = TRUE)))
  
  
  #Convert to mm for calculating ESD
  data_select = data_select %>% 
    mutate(area_mm2=object_area * (process_particle_pixel_size_mm**2),
           major_mm  = object_major * process_particle_pixel_size_mm,
           
           minor_mm  = object_minor * process_particle_pixel_size_mm,
           
           area_exc_mm2  = object_area_exc * (process_particle_pixel_size_mm**2),
           
           area_majmin_mm2  = pi * major_mm/2 * minor_mm/2,
           
           esd_mm  = 2 * (sqrt(area_mm2/pi)),
           
           esd_exc_mm  = 2 * (sqrt(area_exc_mm2/pi)),
           
           esd_maj_min_mm = 2 * (sqrt(area_majmin_mm2/pi))
           
    )
  

  
  print(unique(data_select$object_annotation_category))
  
  #Extract size fraction
  data_final=data_select%>%
    mutate(size_fraction=as.factor(acq_max_mesh))
  
  #Add size group as a category
  data_final %>%
    mutate(size_fraction = case_when(
      esd_mm >= 0.2 & esd_mm < 0.5  ~ '0.2-0.5',
      esd_mm >= 0.5 & esd_mm < 1    ~ '0.5-1',
      esd_mm >= 1   & esd_mm < 2    ~ '1-2',
      esd_mm > 2                          ~ '>2',
      TRUE                                      ~ 'Other'
    ))-> data_final
  
  
 return(data_final)

  
}


#Function to Convert Zooscan measurement to Carbon Biomass
# USe equations from Laveniegos and Ohman which use length-carbon regressions
# Many Use total length which 
transform_by_taxa_group <- function(df, length_type) {
  
  if (length_type == "esd"){
  df %>% mutate(dryweight_C_ug=case_when(
    object_annotation_category=="Copepoda<Maxillopoda" ~ copepods(esd_mm),
    object_annotation_category=="Calanoida" ~ calanoids_pl(esd_mm),
    object_annotation_category=="Oithonidae" ~ copepods_pl(esd_mm),
    object_annotation_category=="Harpacticoida" ~ copepods_pl(esd_mm),
    object_annotation_category=="Poecilostomatoida" ~ copepods_pl(esd_mm),
    object_annotation_category=="Eucalanidae" ~ copepods_pl(esd_mm),
    object_annotation_category=="Euphausiacea" ~ euphausiids(esd_mm ),
    object_annotation_category=="Hydrozoa" ~ hydrozoans_tl(esd_mm),
    object_annotation_category=="Polychaeta" ~ polychaetes_tl(esd_mm),
    object_annotation_category=="Ostracoda" ~ ostracods_tl(esd_mm),
    object_annotation_category=="Eumalacostraca" ~ crustacea_other_tl(esd_mm),
    object_annotation_category=="tetrazoid" ~ pyrosomes(esd_mm),
    object_annotation_category=="Salpida" ~ salps(esd_mm),
    object_annotation_category=="Hyperiidea" ~ hyperiids_tl(esd_mm),
    object_annotation_category=="Pteropoda" ~ thecosomes(esd_mm),
    object_annotation_category=="Doliolida" ~ doliolids(esd_mm),
    object_annotation_category=="Chaetognatha" ~ chaetognaths_tl(esd_mm),
    TRUE ~ NA_real_))->df
    
  } else if (length_type == "feret"){
    df %>% mutate(dryweight_C_ug=case_when(
      object_annotation_category=="Copepoda<Maxillopoda" ~ copepods(object_feret),
      object_annotation_category=="Calanoida" ~ copepods(object_feret),
      object_annotation_category=="Oithonidae" ~ copepods(object_feret),
      object_annotation_category=="Harpacticoida" ~ copepods(object_feret),
      object_annotation_category=="Poecilostomatoida" ~ copepods(object_feret),
      object_annotation_category=="Eucalanidae" ~ copepods(object_feret),
      object_annotation_category=="Euphausiacea" ~ euphausiids(object_feret),
      object_annotation_category=="Hydrozoa" ~ hydrozoans(object_feret),
      object_annotation_category=="Polychaeta" ~ polychaetes(object_feret),
      object_annotation_category=="Ostracoda" ~ ostracods(object_feret),
      object_annotation_category=="Eumalacostraca" ~ crustacea_other(object_feret),
      object_annotation_category=="tetrazoid" ~ pyrosomes(object_feret),
      object_annotation_category=="Salpida" ~ salps(object_feret),
      object_annotation_category=="Hyperiidea" ~ hyperiids(object_feret),
      object_annotation_category=="Pteropoda" ~ thecosomes(object_feret),
      object_annotation_category=="Doliolida" ~ doliolids(object_feret),
      object_annotation_category=="Chaetognatha" ~ chaetognaths(object_feret),
      object_annotation_category=="Oikopleuridae" ~ appendicularians(object_feret),
      
      
      
      TRUE ~ NA_real_))->df
      
  } else {
      # Default to feret
      print("Default to using feret")
      df %>% mutate(dryweight_C_ug=case_when(
        object_annotation_category=="Copepoda<Maxillopoda" ~ copepods(object_feret),
        object_annotation_category=="Calanoida" ~ copepods(object_feret),
        object_annotation_category=="Oithonidae" ~ copepods(object_feret),
        object_annotation_category=="Harpacticoida" ~ copepods(object_feret),
        object_annotation_category=="Poecilostomatoida" ~ copepods(object_feret),
        object_annotation_category=="Eucalanidae" ~ copepods(object_feret),
        object_annotation_category=="Euphausiacea" ~ euphausiids(object_feret ),
        object_annotation_category=="Hydrozoa" ~ hydrozoans(object_feret ),
        object_annotation_category=="Polychaeta" ~ polychaetes(object_feret ),
        object_annotation_category=="Ostracoda" ~ ostracods(object_feret ),
        object_annotation_category=="Eumalacostraca" ~ crustacea_other(object_feret ),
        object_annotation_category=="tetrazoid" ~ pyrosomes(object_feret ),
        object_annotation_category=="Salpida" ~ salps(object_feret ),
        object_annotation_category=="Hyperiidea" ~ hyperiids(object_feret ),
        object_annotation_category=="Pteropoda" ~ thecosomes(object_feret ),
        object_annotation_category=="Doliolida" ~ doliolids(object_feret ),
        object_annotation_category=="Chaetognatha" ~ chaetognaths(object_feret ),
        TRUE ~ NA_real_))->df
      }
    
    
  
  
  return(df)
}

#Taxon-specific functions for biomass from Laveniegos and Ohman 2007
# '_pl' or _tl' converts to total length from ESD Cornilis et al. 2022
copepods <- function(ESD) {
  log_C_microgram <- -6.76 + 2.512 * log10(ESD*1000)
  C_ug=10^(log_C_microgram)
  return(C_ug)
}

copepods_pl <- function(ESD) {
  a=0.031
  b=1.2
  PL=ESD/b-a
  log_C_microgram <- -6.76 + 2.512 * log10(PL*1000)
  C_ug=10^(log_C_microgram)
  return(C_ug)
}

euphausiids <- function(ESD) {
  log_C_microgram <- -0.473 + 3.174 * log10(ESD)
  C_ug=10^(log_C_microgram)
  return(C_ug)
}

euphausiids_tl <- function(ESD) {
  a=0.34
  b=2
  TL=ESD/b-a
  log_C_microgram <- -0.473 + 3.174 * log10(TL)
  C_ug=10^(log_C_microgram)
  return(C_ug)
}

ostracods <- function(ESD) {
  C_ug =(17.072*ESD^2.545)*0.398
  return(C_ug)
}


ostracods_tl <- function(ESD) {
  a=0.097
  b=1.4
  TL=ESD/b-a
  C_ug =(17.072*ESD^2.545)*0.398
  return(C_ug)
}

hyperiids <- function(ESD) {
  log_C_mg = 2.314 + 2.957*log10(ESD) #mg
  C_ug=10^(log_C_mg)*1000*0.365
  return(C_ug)
}

hyperiids_tl <- function(ESD) {
  a=0.17
  b=1.5
  TL=ESD/b-a
  log_C_mg = 2.314 + 2.957*log10(TL) #mg
  C_ug=10^(log_C_mg)*1000*0.365
  return(C_ug)
}

#Use average from decapod groups and TL conversion from Crustacea
crustacea_other <- function(ESD) {
  coefs=mean(0.133,0.322,0.810)
  exps=mean(2.44,2.31,1.77)
  C_mg =coefs*(ESD)^exps
  C_ug=C_mg*1000
  
  return(C_ug)
}

#Use average from decapod groups and TL conversion from Crustacea
crustacea_other_tl <- function(ESD) {
  a=0.064
  b=1.7
  TL=ESD/b-a
  coefs=mean(0.133,0.322,0.810)
  exps=mean(2.44,2.31,1.77)
  C_mg =coefs*(TL)^exps
  C_ug=C_mg*1000
  
  return(C_ug)
}

appendicularians <- function(ESD) {
  DW_ug = 38.8*ESD^2.574
  C_ug = 0.49*DW_ug^1.12
  return(C_ug)
}

appendicularians_tl <- function(ESD) {
  a=0.074
  b=0.59
  TL=ESD/b-a
  DW_ug = 38.8*TL^2.574
  C_ug = 0.49*DW_ug^1.12
  return(C_ug)
}

doliolids<- function(ESD) {
  C_ug = 0.51*(ESD)^2.28
  return(C_ug)
}

#USe average from all salp measurements
salps<- function(ESD) {
  coefs=mean(10.91,5.10,1.00,0.47,0.20,3.00,1.40,1.01,1.62)
  exps=mean(1.54,1.75,2.26,2.22,2.60,1.81,2.05,2.06,1.93)
  C_ug=coefs*ESD^exps
  return(C_ug)

}

pyrosomes<- function(ESD) {
  DW_mg = 0.111*ESD^1.90
  C_ug=DW_mg*.113*1000
  return(C_ug)
}

thecosomes<- function(ESD) {
  log_C_ug = 1.469 + 3.102*log10(ESD)
  C_ug=10^(log_C_ug)
  return(C_ug)
}

chaetognaths <- function(ESD) {
  ESD=as.numeric(ESD)  # Add this line to print the class of ESD
  C_ug = 0.0956 * ESD^2.9093
  return(C_ug)
}

chaetognaths_tl <- function(ESD) {
  a=0.25
  b=4.2
  TL=ESD/b-a
  ESD=as.numeric(TL)  # Add this line to print the class of ESD
  C_ug = 0.0956 * TL^2.9093
  return(C_ug)
}

polychaetes<- function(ESD) {
  C_ug = 7.58*ESD^1.3848
  return(C_ug)
}

polychaetes_tl<- function(ESD) {
  a=0.22
  b=1.9
  PL=ESD/b-a
  C_ug = 7.58*ESD^1.3848
  return(C_ug)
}

hydrozoans<- function(ESD) {
  C_ug = 20.47*ESD^0.834
  return(C_ug)
}

hydrozoans_tl<- function(ESD) {
  #Use Cnidarian Values
  a=0.22
  b=0.9
  PL=ESD/b-a
  C_ug = 20.47*ESD^0.834
  return(C_ug)
}


#References
# Cornils, A. et al. (2022) ‘Testing the usefulness of optical data for zooplankton long-term monitoring: Taxonomic composition, abundance, biomass, and size spectra from ZooScan image analysis’, 
# Limnology and Oceanography: Methods, 20(7), pp. 428–450. Available at: https://doi.org/10.1002/lom3.10495.

# Lavaniegos, B.E. and Ohman, M.D. (2007) ‘Coherence of long-term variations of zooplankton in two sectors of the California Current System’, 
# Progress in Oceanography, 75(1), pp. 42–69. Available at: https://doi.org/10.1016/j.pocean.2007.07.002.

