library(pacman)
p_load(data.table, dplyr, ggplot2, RColorBrewer, viridis, ggthemes, stringr, ggpubr, 
ComplexHeatmap, colorRamp2, ggseqlogo, scales, ggrastr, GenomicRanges, ggrepel, circlize)
options(scipen=999)
source("code/common.R")

# A
readRDS("input/Fig3_A.RDS") %>% 
  ggplot(aes(technology, value, fill=stage)) +
    geom_col(position="dodge", colour="black") +
    theme_bw() +
    ylab("Number of targets, M") +
    theme_Publication() +
    facet_wrap(~Var3, scales="free") +
    scale_fill_manual(values=c(brewer.pal(6, "Greys")[2], brewer.pal(6, "Greys")[4], brewer.pal(6, "Greys")[5])) +
    scale_y_continuous(labels = label_number(suffix = " M", scale = 1e-6)) +
    theme(axis.title.x = element_blank(),
          legend.key.size = unit(0.5, 'cm'),
          legend.position="top")

# B
readRDS("input/Fig3_B.RDS") %>% 
  ggplot(aes(PC1, PC2, colour=biological, shape=replicate, alpha=1)) +
    geom_point(size=2) +
    theme_Publication() + 
    scale_shape_manual(values=c(1, 4)) +
    scale_colour_brewer(palette="Dark2") +
    guides(colour = guide_legend(nrow=2,byrow=TRUE, override.aes=list(shape=21)),
          shape = guide_legend(nrow=2,byrow=TRUE)) 
 
# C
readRDS("input/Fig3_C.RDS") %>% 
  ggplot(aes(stage, V1, fill=rn)) +
    geom_col(position="stack") +
    theme_Publication() + 
    scale_fill_manual(values = cols_posterGroups2, name="") +
    scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) +
    ylab("Fraction of regions, %") +
    theme(axis.title.x = element_blank(),
          legend.key.size = unit(0.5, 'cm'),
          legend.position = "top") +
    guides(fill=guide_legend(nrow=2,byrow=FALSE))

# D
cmr <- "complete"
cdr <- "euclidean"                
col_fun <- colorRamp2(c(0, 0.3, 1), viridis(3) )
col_uCG <- colorRamp2(c(0, 0.2, 0.7), plasma(3) )
col_hmC <- colorRamp2(c(0, 0.1, 0.2), inferno(4)[2:4]) 

ht_opt$simple_anno_size <- unit(3, "mm")
modWidth <- 5.8
openWidth <- 5.8
heatmapHeight <- 10

data <- readRDS("input/Fig3_D.RDS")

ht <- Heatmap(data[, .(up_TT_S0_GCCG, up_TT_S1_GCCG, up_TT_S2_GCCG) ], 
              name = "upstream",
              heatmap_legend_param = list(direction = "horizontal", legend_width = unit(2, "cm"), at = c(0, 0.5, 1)),
              column_title_rot = 90, 
              cluster_rows=FALSE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_fun,
              use_raster=TRUE,
              width = unit(openWidth, "mm"),
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),                   
              row_title_rot = 0,
              column_title=NULL,
              show_column_names=FALSE,
              row_title=NULL) +               
      Heatmap(data[, .(pu_TT_S0_GCCG, pu_TT_S1_GCCG,pu_TT_S2_GCCG) ], 
              name = "1kbUP",
              heatmap_legend_param = list(direction = "horizontal", legend_width = unit(2, "cm"), at = c(0, 0.5, 1)),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_fun,
              show_heatmap_legend = FALSE,
              use_raster=TRUE,
              width = unit(openWidth, "mm"),
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),                           
               row_title_rot = 0,
              column_title=NULL,
              show_column_names=FALSE,
              row_title=NULL) +
      Heatmap(data[, .(pd_TT_S0_GCCG, pd_TT_S1_GCCG, pd_TT_S2_GCCG) ], 
              name = "1kDown",
              heatmap_legend_param = list(direction = "horizontal", legend_width = unit(2, "cm"), at = c(0, 0.5, 1)),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_fun,
              show_heatmap_legend = FALSE,
              use_raster=TRUE,
              width = unit(openWidth, "mm"),
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),             
               row_title_rot = 0,
              column_title=NULL,
              show_column_names=FALSE,
              row_title=NULL) +                
      Heatmap(data[, .(gb_TT_S0_GCCG, gb_TT_S1_GCCG, gb_TT_S2_GCCG) ], 
              name = "body",
              heatmap_legend_param = list(direction = "horizontal", legend_width = unit(2, "cm"), at = c(0, 0.5, 1)),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_fun,
              show_heatmap_legend = FALSE,
              use_raster=TRUE,
              width = unit(openWidth, "mm"),
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),                  
               row_title_rot = 0,
              column_title=NULL,
              show_column_names=FALSE,
              row_title=NULL) +  
###################### uCG fraction                            
      Heatmap(data[, .(up_TT_S0_GCCG_uCG, up_TT_S1_GCCG_uCG, up_TT_S2_GCCG_uCG) ], 
              name = "upstream_uCG",
              heatmap_legend_param = list(direction = "horizontal", legend_width = unit(2, "cm"), at = c(0, 0.4, 0.8)),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_uCG,
              width = unit(modWidth, "mm"),
              use_raster=TRUE,
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),                 
               row_title_rot = 0,
              column_title=NULL,
              show_column_names=FALSE,
              row_title=NULL) +
      Heatmap(data[, .(pu_TT_S0_GCCG_uCG, pu_TT_S1_GCCG_uCG, pu_TT_S2_GCCG_uCG) ], 
              name = "1kbUP_uCG",
              heatmap_legend_param = list(direction = "horizontal", legend_width = unit(2, "cm"), at = c(0, 0.4, 0.8)),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_uCG,
              show_heatmap_legend = FALSE,
              width = unit(modWidth, "mm"),
              use_raster=TRUE,
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),                             
               row_title_rot = 0,
              column_title=NULL,
              show_column_names=FALSE,
              row_title=NULL) +
      Heatmap(data[, .(pd_TT_S0_GCCG_uCG, pd_TT_S1_GCCG_uCG, pd_TT_S2_GCCG_uCG) ], 
              name = "1kDown_uCG",
              heatmap_legend_param = list(direction = "horizontal", legend_width = unit(2, "cm"), at = c(0, 0.4, 0.8)),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_uCG,
              show_heatmap_legend = FALSE,
              width = unit(modWidth, "mm"),
              use_raster=TRUE,
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),                
               row_title_rot = 0,
              column_title=NULL,
              show_column_names=FALSE,
              row_title=NULL) +                
      Heatmap(data[, .(gb_TT_S0_GCCG_uCG, gb_TT_S1_GCCG_uCG, gb_TT_S2_GCCG_uCG) ], 
              name = "body_uCG",
              heatmap_legend_param = list(direction = "horizontal", legend_width = unit(2, "cm"), at = c(0, 0.4, 0.8)),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_uCG,
              show_heatmap_legend = FALSE,
              width = unit(modWidth, "mm"),
              use_raster=TRUE,
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),                
               row_title_rot = 0,
              column_title=NULL,
              show_column_names=FALSE,
              row_title=NULL) + 
####################### hmC
      Heatmap(data[, .(up_TT_S0_GCCG_hmC, up_TT_S1_GCCG_hmC, up_TT_S2_GCCG_hmC) ], 
              name = "upstream_hmC",
              heatmap_legend_param = list(direction = "horizontal", legend_width = unit(2, "cm"), at = c(0, 0.1, 0.2)),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_hmC,
              width = unit(modWidth, "mm"),
              use_raster=TRUE,
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),        
        
               row_title_rot = 0,
              column_title=NULL,
              show_column_names=FALSE,
              row_title=NULL) +
      Heatmap(data[, .(pu_TT_S0_GCCG_hmC, pu_TT_S1_GCCG_hmC, pu_TT_S2_GCCG_hmC) ], 
              name = "1kbUP_hmC",
              heatmap_legend_param = list(direction = "horizontal", legend_width = unit(2, "cm"), at = c(0, 0.1, 0.2)),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_hmC,
              show_heatmap_legend = FALSE,
              width = unit(modWidth, "mm"),
              use_raster=TRUE,
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),                         
              row_title_rot = 0,
              column_title=NULL,
              show_column_names=FALSE,
              row_title=NULL) +
      Heatmap(data[, .(pd_TT_S0_GCCG_hmC, pd_TT_S1_GCCG_hmC, pd_TT_S2_GCCG_hmC) ], 
              name = "1kDown_hmC",
              heatmap_legend_param = list(direction = "horizontal", legend_width = unit(2, "cm"), at = c(0, 0.1, 0.2)),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_hmC,
              show_heatmap_legend = FALSE,
              width = unit(modWidth, "mm"),
              use_raster=TRUE,
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),                    
              row_title_rot = 0,
              column_title=NULL,
              show_column_names=FALSE,
              row_title=NULL) +                    
      Heatmap(data[, .(gb_TT_S0_GCCG_hmC, gb_TT_S1_GCCG_hmC, gb_TT_S2_GCCG_hmC) ], 
              name = "body_hmC",
              heatmap_legend_param = list(direction = "horizontal", legend_width = unit(2, "cm"), at = c(0, 0.1, 0.2)),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_hmC,
              show_heatmap_legend = FALSE,
              width = unit(modWidth, "mm"),
              use_raster=TRUE,
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),               
              row_title_rot = 0,
              column_title=NULL,
              show_column_names=FALSE,
              row_title=NULL)      
  draw(ht, 
      main="1kbUP", 
      merge_legend = TRUE, 
      heatmap_legend_side = "bottom", 
      row_split=data$skirtukas,
      show_heatmap_legend = FALSE,
      cluster_row_slices = FALSE,
      ht_gap = unit(c(rep(1, 3), 3, rep(1, 3), 3, rep(1, 3), 1), "mm"),
      row_gap = unit(c(1), "mm"),
      annotation_legend_side = "bottom",
      show_row_dend = FALSE) 

# E
data <- readRDS("input/Fig3_E.RDS")
col_fpkm <- colorRamp2(c(0, 5, 10),colorRampPalette(brewer.pal(9,"Reds"))(3))
col_fun <- colorRamp2(c(0, 0.5, 1), viridis(3))  
col_fun2 <- colorRamp2(c(0, 0.3, 0.6), cividis(3))  
col_RNA_groups <- c(
    "dideja" = brewer.pal(n = 3, name = "Reds")[3],
    "padidejo\nnekito" = brewer.pal(n = 3, name = "Reds")[2],
    "nekito\npadidejo" = brewer.pal(n = 3, name = "Reds")[1],
    "mazeja" = brewer.pal(n = 3, name = "Blues")[3],
    "mazejo\nnekito" = brewer.pal(n = 3, name = "Blues")[2],
    "nekito\nmazejo" = brewer.pal(n = 3, name = "Blues")[1],
    "mazejo\ndidejo" = brewer.pal(n = 3, name = "Greens")[2],
    "didejo\nmazejo" = brewer.pal(n = 3, name = "Purples")[2],
    "KITA" = "black"
)  
col_scaled <- colorRamp2(c(-1, 0, 1), viridis(3))  

cmr <- "complete"
cdr <- "euclidean"    

ht <- Heatmap((data[, .(up_S0, up_S1, up_S2) ] %>% t %>% scale %>% t), 
              name = "upstream",
              heatmap_legend_param = list(direction = "vertical"),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_scaled ,
              use_raster=TRUE,
              width = unit(modWidth, "mm"),
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),    
              row_title_rot = 0) +  
      Heatmap((data[, .(prom_S0, prom_S1, prom_S2) ] %>% t %>% scale %>% t), 
              name = "Prom1kbaround",
              heatmap_legend_param = list(direction = "vertical"),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_scaled ,
              use_raster=TRUE,
              width = unit(modWidth, "mm"),
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),                 
              row_title_rot = 0) +       
      Heatmap((data[, .(body_S0, body_S1, body_S2) ] %>%  t %>% scale %>% t ), 
              name = "GeneBody",
              heatmap_legend_param = list(direction = "vertical"),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              clustering_distance_rows = cdr,
              clustering_method_rows = cmr, 
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              col=col_scaled ,
              use_raster=TRUE,
              width = unit(modWidth, "mm"),
              height=unit(heatmapHeight, "cm"),
              top_annotation = HeatmapAnnotation(                                            
                                                 stage = c("S0", "S1", "S2"), 
                                                 col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                                                 show_annotation_name = FALSE,
                                                 show_legend = TRUE,
                                                 border=TRUE,
                                                 annotation_name_side = "left"),  
              row_title_rot = 0) +                            
       Heatmap(data[, .(RNA_group) ], 
              name = "RNA_group",
              heatmap_legend_param = list(direction = "vertical"),
              column_title_rot = 90, 
              cluster_rows=TRUE,
              cluster_column_slices = FALSE,
              cluster_column = FALSE,
              use_raster=TRUE,
              width = unit(2, "mm"),
              height=unit(heatmapHeight, "cm"),
              
              col=col_RNA_groups)

draw(ht, 
    row_split=data$prom_group, 
    merge_legend = TRUE, 
    heatmap_legend_side = "right", 
    show_heatmap_legend = TRUE, 
    ht_gap = unit(c(2), "mm"),
    row_gap = unit(c(2), "mm"),
    row_order = data$eil_NR,
    annotation_legend_side = "right",
    show_row_dend = FALSE) 