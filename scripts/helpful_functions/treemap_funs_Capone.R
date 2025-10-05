#Tree map functions

library(dplyr)
#########FUNCTIONS:
## Define function - phyloseq_transform_to_long (transforming phyloseq in dataframe)
phyloseq_transform_to_long <- function(ps) {
  otu_df <- as.data.frame(ps@otu_table@.Data, stringsAsFactors = FALSE) %>%
    rownames_to_column(var = "asv_code") %>%
    pivot_longer(cols = -asv_code,
                 names_to = "file_code",
                 values_to = "n_reads",
                 values_drop_na = TRUE) %>%
    filter(n_reads != 0) %>%
    filter(!is.na(n_reads))
  
  # See https://urldefense.com/v3/__https://github.com/joey711/phyloseq/issues/983__;!!Mih3wA!DmqmyZdyBYO-tZIXCBV7LXJzHsyr-fg8Q1Y7E940uEsgmBY-V-p8b1DwGtZvj4udaRuNQ_xDeoOPVQZz$ 
  taxo_df <- as.data.frame(ps@tax_table@.Data, stringsAsFactors = FALSE) %>%
    rownames_to_column(var = "asv_code")
  
  otu_df <- left_join(taxo_df, otu_df)
  
  metadata_df <- data.frame(sample_data(ps)) %>%
    rownames_to_column(var = "file_code")
  
  otu_df <- left_join(otu_df, metadata_df, by = c("file_code"))
  
  return(otu_df)
  
}



## Define function - phyloseq_normalize_median (actual normalization)
phyloseq_normalize_median <- function (ps) {
  ps_median = median(sample_sums(ps))
  normalize_median = function(x, t=ps_median) (if(sum(x) > 0){ round(t * (x / sum(x)))} else {x})
  ps = transform_sample_counts(ps, normalize_median)
  cat(str_c("\n========== \n") )
  print(ps)
  cat(sprintf("\n==========\nThe median number of reads used for normalization is  %.0f", ps_median))
  return(ps)
}

read_taxtable <- function (taxonomy.file, sep = ",") {
  
  s.tax <- read.csv(taxonomy.file, row.names=1, check.names=FALSE, sep = sep)
  s.taxmat <- as.matrix(s.tax)
  
  tax_table(s.taxmat)
  
}


#
custom_palette <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a")


phyloseq_long_treemap <- function(df, group1, group2, title, colors=NULL, label_group1 = TRUE) {
  
  df <- df %>%
    group_by({{group1}}, {{group2}}) %>%
    summarise(n_reads=sum(n_reads, na.rm = TRUE))
  
  # # Replace blank values with NA using mutate
  # df <- df %>%
  #   mutate({{group2}} = ifelse({{group2}} == "", NA, {{group2}}))%>%
  #   mutate({{group1}} = ifelse({{group1}} == "", NA, {{group1}}))
  
  
  #Remove NA and fake NA
  df=na.omit(df)
  df=df %>%
    filter({{group1}}!='NA' | {{group2}}!='NA')
  
  gg <- ggplot(df, aes(area = (n_reads),
                       fill = {{group1}},
                       label = {{group2}},
                       subgroup = {{group1}})) +
    ggtitle(title) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_subgroup_border() +
    treemapify::geom_treemap_text(colour = "black", place = "topleft", reflow = T,
                                  padding.x =  grid::unit(2, "mm"),
                                  padding.y = grid::unit(4, "mm"),
                                  min.size=6)  +
    theme(legend.position="none", plot.title = element_text(size = 16, face = "bold"))
  
  if (label_group1){
    gg <- gg +
      treemapify::geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.8, colour =
                                               "white", fontface = "italic", min.size = 0)
  }
  
  if (is.null(colors)){
    gg <- gg + scale_fill_viridis_d()
  } else {
    gg <- gg + scale_fill_manual(values = colors)
  }
  
  print(gg)
  return(gg)
  treemap_list <- list(gg = gg, df=df)
  return(treemap_list)
}


## TOP Taxa 

phyloseq_long_treemap_top <- function(df, group1, group2, title,top, colors=NULL, label_group1 = TRUE) {
  
  df <- df %>%
    group_by({{group1}}, {{group2}}) %>%
    summarise(n_reads=sum(n_reads, na.rm = TRUE)) %>%
    arrange(desc(n_reads))
  
  # # Replace blank values with NA using mutate
  # df <- df %>%
  #   mutate({{group2}} = ifelse({{group2}} == "", NA, {{group2}}))%>%
  #   mutate({{group1}} = ifelse({{group1}} == "", NA, {{group1}}))
  
  
  #Remove NA and fake NA
  df=na.omit(df)
  df=df %>%
    filter({{group1}}!='NA' | {{group2}}!='NA')
  #Select top 20
  df1=df[1:top,]
  
  print((df1 %>% select({{group2}})))
  gg <- ggplot(df1, aes(area = (n_reads),
                        fill = {{group1}},
                        label = {{group2}},
                        subgroup = {{group1}})) +
    ggtitle(title) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_subgroup_border() +
    treemapify::geom_treemap_text(colour = "black", place = "topleft", reflow = T,
                                  padding.x =  grid::unit(2, "mm"),
                                  padding.y = grid::unit(4, "mm"),
                                  min.size=6)  +
    theme(legend.position="none", plot.title = element_text(size = 16, face = "bold"))
  
  if (label_group1){
    gg <- gg +
      treemapify::geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.8, colour =
                                               "white", fontface = "italic", min.size = 0)
  }
  
  if (is.null(colors)){
    gg <- gg + scale_fill_viridis_d()
  } else {
    gg <- gg + scale_fill_manual(values = colors)
  }
  # print(gg)
  return(gg)
  treemap_list <- list(gg = gg, df=df)
  return(treemap_list)
}



# Next Taxa -------------------------------------------------------------

## 11-20

phyloseq_long_treemap_next <- function(df, group1, group2, title,subset_taxa, colors=NULL, label_group1 = TRUE) {
  
  df <- df %>%
    group_by({{group1}}, {{group2}}) %>%
    summarise(n_reads=sum(n_reads, na.rm = TRUE)) %>%
    arrange(desc(n_reads))
  
  # # Replace blank values with NA using mutate
  # df <- df %>%
  #   mutate({{group2}} = ifelse({{group2}} == "", NA, {{group2}}))%>%
  #   mutate({{group1}} = ifelse({{group1}} == "", NA, {{group1}}))
  
  
  #Remove NA and fake NA
  df=na.omit(df)
  df=df %>%
    filter({{group1}}!='NA' | {{group2}}!='NA')
  #Select top 20
  df1=df[subset_taxa:subset_taxa+10,]
  
  print((df1 %>% select({{group2}})))
  gg <- ggplot(df1, aes(area = (n_reads),
                        fill = {{group1}},
                        label = {{group2}},
                        subgroup = {{group1}})) +
    ggtitle(title) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_subgroup_border() +
    treemapify::geom_treemap_text(colour = "black", place = "topleft", reflow = T,
                                  padding.x =  grid::unit(2, "mm"),
                                  padding.y = grid::unit(4, "mm"),
                                  min.size=6)  +
    theme(legend.position="none", plot.title = element_text(size = 16, face = "bold"))
  
  if (label_group1){
    gg <- gg +
      treemapify::geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.8, colour =
                                               "white", fontface = "italic", min.size = 0)
  }
  
  if (is.null(colors)){
    gg <- gg + scale_fill_viridis_d()
  } else {
    gg <- gg + scale_fill_manual(values = colors)
  }
  # print(gg)
  return(gg)
  treemap_list <- list(gg = gg, df=df)
  return(treemap_list)
}
