
map_asv_norm <- function(data_in,ranck) {
  # Load California map data
  worldmap <- map_data("world")
  states <- map_data("state")
  ca_df <- subset(states, region == "california")
  
  # title = readline("Enter the plot title: ")
  
  data=data_in %>% filter(rank==ranck)%>%
    #group by cycle, size fraction and hash
    group_by(cycle,asv_code,max_size) %>% 
    #Sum each hash occurence
    mutate(sum_asv=sum(n_reads)) %>%
    ungroup() %>%
    #sum each hash at each cycle
    group_by(Latitude,max_size) %>%
    mutate(nreads_cycle=sum(n_reads)) %>%
    mutate(percent = sum_asv/nreads_cycle)
  
  # Find the maximum value of eDNA index and its corresponding location
  max_value <- max(data$nreads_cycle)
  max_lat <- data$Latitude[data$nreads_cycle == max_value]+1
  max_lon <- data$Longitude[data$nreads_cycle == max_value]+1
  
  max_values <- data %>%
    filter(nreads_cycle == max(nreads_cycle)) %>%
    group_by(max_size) %>%
    slice_max(nreads_cycle) 
  
  
  p=ggplot(worldmap) +
    geom_map(data = worldmap, map = worldmap, aes(map_id=region), col = "white", fill = "gray50") +
    geom_point(data=data, aes(x=Longitude,y=Latitude,size=nreads_cycle, color=cycle))+ 
    #Add point at location of max
    geom_label(data = max_values, aes(x = Longitude, y = Latitude, label = "*"),
               color = "red", fill = alpha("white", 0), size = 10,
               label.padding = unit(0.3, "lines"), label.size = 0) +
    
    
    
    scale_size(range = c(2,10))+
    coord_fixed(xlim = c(-134, -119.0),  ylim = c(34, 38), ratio = 1.3)+
    ggtitle(paste(data$Species[1],"(#",data$rank[1],"most prevalent ASV, Normalized Reads)"), paste("(",data$asv_code[1],")"))+ 
    scale_x_continuous(breaks = seq(-118,-132, by = -2))+
    xlab("Latitude")+
    ylab("Longitude")+
    labs(size="Relative Read Abundance")+
    theme_classic()+
    facet_wrap(~ max_size, nrow = 3)
  
  # Save the plot with a title using here()
  file_path <- here("C:/Users/Dante Capone/OneDrive/Desktop/Scripps_PhD/Decima_Matthews_Gallego_Capone/figures/Maps/")
  
  saving=0
  if(saving==1){
    
    ggsave(paste(file_path,"/",readline("Enter the plot title: ")), p, width = 8, height = 6, units = "in")
  }
  
  return(p)
}


map_asv_raw <- function(data_in,ranck) {
  # Load California map data
  worldmap <- map_data("world")
  states <- map_data("state")
  ca_df <- subset(states, region == "california")
  
  # title = readline("Enter the plot title: ")
  
  data=data_in %>% filter(rank==ranck)%>%
    #group by cycle, size fraction and hash
    group_by(cycle,asv_code,max_size) %>% 
    #Sum each hash occurence
    mutate(sum_asv=sum(n_reads)) %>%
    ungroup() %>%
    #sum each hash at each cycle
    group_by(Latitude,max_size) %>%
    mutate(nreads_cycle=sum(n_reads)) %>%
    mutate(percent = sum_asv/nreads_cycle)
  
  # Find the maximum value of eDNA index and its corresponding location
  max_value <- max(data$nreads_cycle)
  max_lat <- data$Latitude[data$nreads_cycle == max_value]+1
  max_lon <- data$Longitude[data$nreads_cycle == max_value]+1
  
  max_values <- data %>%
    filter(nreads_cycle == max(nreads_cycle)) %>%
    group_by(max_size) %>%
    slice_max(nreads_cycle) 
  
  
  p=ggplot(worldmap) +
    geom_map(data = worldmap, map = worldmap, aes(map_id=region), col = "white", fill = "gray50") +
    geom_point(data=data, aes(x=Longitude,y=Latitude,size=nreads_cycle, color=cycle))+ 
    #Add point at location of max
    geom_label(data = max_values, aes(x = Longitude, y = Latitude, label = "*"),
               color = "red", fill = alpha("white", 0), size = 10,
               label.padding = unit(0.3, "lines"), label.size = 0) +
    
    
    
    scale_size(range = c(2,10))+
    coord_fixed(xlim = c(-134, -119.0),  ylim = c(34, 38), ratio = 1.3)+
    ggtitle(paste(data$Species[1],"(#",data$rank[1],"most prevalent ASV, Raw Reads)"), paste("(",data$asv_code[1],")"))+ 
    scale_x_continuous(breaks = seq(-118,-132, by = -2))+
    xlab("Latitude")+
    ylab("Longitude")+
    theme_classic()+
    facet_wrap(~ max_size, nrow = 3)
  
  # Save the plot with a title using here()
  file_path <- here("C:/Users/Dante Capone/OneDrive/Desktop/Scripps_PhD/Decima_Matthews_Gallego_Capone/figures/Maps/")
  
  saving=0
  if(saving==1){
    
    ggsave(paste(file_path,"/",readline("Enter the plot title: ")), p, width = 8, height = 6, units = "in")
  }
  
  return(p)
  
}
  
  
  
  ### Map for unranked data
  
  map_asv_off_on <- function(data_in,shore_ranck,primer,edna_or_normalized) {
    # Load California map data
    worldmap <- map_data("world")
    states <- map_data("state")
    ca_df <- subset(states, region == "california")
    
    # title = readline("Enter the plot title: ")
    
    data=data_in %>% 
      filter(rank_off_on==shore_ranck)%>%
      #group by cycle, size fraction and hash
      group_by(cycle,asv_code,max_size) %>% 
      #Sum each hash occurence
      mutate(sum_asv=sum(n_reads)) %>%
      ungroup() %>%
      #sum each hash at each cycle
      group_by(Latitude,max_size) %>%
      mutate(nreads_cycle=sum(n_reads)) %>%
      mutate(percent = sum_asv/nreads_cycle)
    
    # Find the maximum value of eDNA index and its corresponding location
    max_value <- max(data$nreads_cycle)
    max_lat <- data$Latitude[data$nreads_cycle == max_value]
    max_lon <- data$Longitude[data$nreads_cycle == max_value]
    
    max_values <- data %>%
      filter(nreads_cycle == max(nreads_cycle)) %>%
      group_by(max_size) %>%
      slice_max(nreads_cycle) 
    
    
    p=ggplot(worldmap) +
      geom_map(data = worldmap, map = worldmap, aes(map_id=region), col = "white", fill = "gray50") +
      geom_point(data=data, aes(x=Longitude,y=Latitude,size=nreads_cycle, color=cycle))+ 
      #Add point at location of max
      geom_label(data = max_values, aes(x = Longitude, y = Latitude, label = "*"),
                 color = "red", fill = alpha("white", 0), size = 10,
                 label.padding = unit(0.3, "lines"), label.size = 0) +
      
      
      
      scale_size(range = c(2,10))+
      coord_fixed(xlim = c(-134, -119.0),  ylim = c(34, 38), ratio = 1.3)+
      ggtitle(paste(data$Species[1],"(#",data$rank[1],"most prevalent",data$offshore_onshore_taxa[1],primer,"ASV )"), paste("(",data$asv_code[1],")"))+ 
      scale_x_continuous(breaks = seq(-118,-132, by = -2))+
      xlab("Latitude")+
      ylab("Longitude")+
      labs(size="Relative Read Abundance")+
      theme_classic()+
      facet_wrap(~ max_size, nrow = 3)
    
    # Save the plot with a title using here()
    file_path <- here("C:/Users/Dante Capone/OneDrive/Desktop/Scripps_PhD/Decima_Matthews_Gallego_Capone/figures/Maps/")
    
    saving=0
    if(saving==1){
      
      ggsave(paste(file_path,"/",readline("Enter the plot title: ")), p, width = 8, height = 6, units = "in")
    }
    
    return(p)
  }
  
