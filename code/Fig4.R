library(pacman)
p_load(data.table, dplyr, ggplot2, RColorBrewer, viridis, ggthemes, stringr, ggpubr, ComplexHeatmap, 
  colorRamp2, ggseqlogo, scales, ggrastr, GenomicRanges, ggrepel, EnrichedHeatmap)
options(scipen=999)
source("./code/common.R")

# A
lis <- readRDS(paste0("input/Fig4_A.RDS"))
S0 <- lis[["S0"]]
S1 <- lis[["S1"]]
S2 <- lis[["S2"]]
d2plot_open <- lis[["d2plot_open"]]
d2plot_closes <- lis[["d2plot_closes"]]
rowO <- readRDS("input/Fig4_ArowO.RDS")

nomM_S0 <- normalizeToMatrix(S0,
                             d2plot_open, 
                             value_column = "group_poster", extend = 3000, target_ratio = 0.3)
nomM_S1 <- normalizeToMatrix(S1,
                             d2plot_open,
                            value_column = "group_poster", extend = 3000, target_ratio = 0.3)
nomM_S2 <- normalizeToMatrix(S2, 
                             d2plot_open,
                              value_column = "group_poster", extend = 3000, target_ratio = 0.3) 
ht_list <- EnrichedHeatmap(nomM_S0, col=states_col, column_title="S0", row_title_rot=0,  use_raster = TRUE, cluster_rows = TRUE, 
                             name="Mod. group", show_heatmap_legend = TRUE, axis_name=axis_names) +
              EnrichedHeatmap(nomM_S1, col=states_col, column_title="S1", use_raster = TRUE, cluster_rows = TRUE,  
                             name="Mod. group S1", show_heatmap_legend = FALSE, axis_name=axis_names) +
              EnrichedHeatmap(nomM_S2, col=states_col, column_title="S2", use_raster = TRUE, cluster_rows = TRUE,  
                             name="Mod. group S2", show_heatmap_legend = FALSE, axis_name=axis_names) 
draw(ht_list, annotation_legend_side = "bottom",  heatmap_legend_side = "bottom", row_order=rowO[[1]])              


nomM_S0 <- normalizeToMatrix(S0,
                             d2plot_closes, 
                             value_column = "group_poster", extend = 3000, target_ratio = 0.3)
nomM_S1 <- normalizeToMatrix(S1,
                             d2plot_closes,
                             value_column = "group_poster", extend = 3000, target_ratio = 0.3)
nomM_S2 <- normalizeToMatrix(S2, 
                             d2plot_closes,
                             value_column = "group_poster", extend = 3000, target_ratio = 0.3) 
ht_list <- EnrichedHeatmap(nomM_S0, col=states_col, column_title="S0", row_title_rot=0,  use_raster = TRUE, cluster_rows = TRUE, 
                             name="Mod. group", show_heatmap_legend = TRUE, axis_name=axis_names) +
              EnrichedHeatmap(nomM_S1, col=states_col, column_title="S1", use_raster = TRUE, cluster_rows = TRUE, 
                             name="Mod. group S1", show_heatmap_legend = FALSE, axis_name=axis_names) +
              EnrichedHeatmap(nomM_S2, col=states_col, column_title="S2", use_raster = TRUE, cluster_rows = TRUE, 
                             name="Mod. group S2", show_heatmap_legend = FALSE, axis_name=axis_names) 

draw(ht_list, annotation_legend_side = "bottom",  heatmap_legend_side = "bottom", row_order=rowO[[2]])              

# B
rez <- readRDS("input/Fig4_B1.RDS")
ggplot(data=rez)  + 
    geom_line(data=rez[stage =="S0", ], aes(rn, fraction, group="1"), colour=brewer.pal(6, "Dark2")[1]) +
    geom_line(data=rez[stage =="S1", ], aes(rn, fraction, group="1"), colour=brewer.pal(6, "Dark2")[2]) +
    geom_line(data=rez[stage =="S2", ], aes(rn, fraction, group="1"), colour=brewer.pal(6, "Dark2")[3]) +
    theme_bw() +
    geom_vline(xintercept = 100, linetype="dotted", color = "grey", size=0.5) +
    geom_vline(xintercept = 187, linetype="dotted", color = "grey", size=0.5) +
    theme_Publication() +
    facet_wrap(~type) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.ticks.x= element_blank(),
          axis.title=element_blank(),
          axis.text.x = element_blank(),
          legend.position="bottom") 

tmp <- readRDS("input/Fig4_B2.RDS") # hmC
ggplot(tmp, aes(colour=stage))  + 
    geom_line(data=tmp[stage =="S0", ], aes(rn, fraction, group="1"), colour=brewer.pal(6, "Dark2")[1]) +
    geom_line(data=tmp[stage =="S1", ], aes(rn, fraction, group="1"), colour=brewer.pal(6, "Dark2")[2]) +
    geom_line(data=tmp[stage =="S2", ], aes(rn, fraction, group="1"), colour=brewer.pal(6, "Dark2")[3]) +
    theme_bw() +
    facet_wrap(~type) +
    theme_bw() +
    geom_vline(xintercept = 10, linetype="dotted", color = "grey", size=0.5) +
    geom_vline(xintercept = 20, linetype="dotted", color = "grey", size=0.5) +
    theme_Publication() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.ticks.x= element_blank(),
          axis.title=element_blank(),
          axis.text.x = element_blank()) 
  
tmp <- readRDS("input/Fig4_B3.RDS") # uCG
ggplot(tmp, aes(colour=stage))  + 
    geom_line(data=tmp[stage =="S0", ], aes(rn, fraction, group="1"), colour=brewer.pal(6, "Dark2")[1]) +
    geom_line(data=tmp[stage =="S1", ], aes(rn, fraction, group="1"), colour=brewer.pal(6, "Dark2")[2]) +
    geom_line(data=tmp[stage =="S2", ], aes(rn, fraction, group="1"), colour=brewer.pal(6, "Dark2")[3]) +
    theme_bw() +
    facet_wrap(~type) +
    theme_bw() +
    geom_vline(xintercept = 10, linetype="dotted", color = "grey", size=0.5) +
    geom_vline(xintercept = 20, linetype="dotted", color = "grey", size=0.5) +
    theme_Publication() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.ticks.x= element_blank(),
          axis.title=element_blank(),
          axis.text.x = element_blank()) 

tmp <- readRDS("input/Fig4_B4.RDS") # hmC full 
ggplot(tmp, aes(colour=stage))  + 
    geom_line(data=tmp[stage =="S0", ], aes(rn, fraction, group="1"), colour=brewer.pal(6, "Dark2")[1]) +
    geom_line(data=tmp[stage =="S1", ], aes(rn, fraction, group="1"), colour=brewer.pal(6, "Dark2")[2]) +
    geom_line(data=tmp[stage =="S2", ], aes(rn, fraction, group="1"), colour=brewer.pal(6, "Dark2")[3]) +
    theme_bw() +
    facet_wrap(~type) +
    theme_bw() +
    geom_vline(xintercept = 10, linetype="dotted", color = "grey", size=0.5) +
    geom_vline(xintercept = 20, linetype="dotted", color = "grey", size=0.5) +
    theme_Publication() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.ticks.x= element_blank(),
          axis.title=element_blank(),
          axis.text.x = element_blank()) 
# C
lis <- readRDS(paste0("input/Fig4_C.RDS"))

fig1 <- rbind(lis[["p1"]], lis[["p5"]]) %>% 
   as.data.table() %>% 
   .[, modas := factor(modas, levels=c("hmC_OCR", "hmC_FULL"))]  %>% 
   ggplot(aes(elementas, value, fill=stadija)) +
    geom_boxplot(outlier.shape=NA, lwd=0.5, fatten = 0.5 ) +
    facet_wrap(~modas, nrow=1) +
    theme_Publication() +
      scale_fill_manual(values=c(brewer.pal(6, "Dark2")[1], brewer.pal(6, "Dark2")[2], brewer.pal(6, "Dark2")[3])) +
      scale_y_continuous(limits=c(0, 0.6))

fig2 <- rbind(lis[["p2"]], lis[["p6"]]) %>% 
   as.data.table() %>% 
   .[, modas := factor(modas, levels=c("hmC_OCR", "hmC_FULL"))]  %>% 
   ggplot(aes(elementas, value, fill=stadija )) +
    geom_boxplot(outlier.shape=NA, lwd=0.5, fatten = 0.5) +
    facet_wrap(~modas, nrow=1) +
    theme_Publication() +
      scale_fill_manual(values=c(brewer.pal(6, "Dark2")[1], brewer.pal(6, "Dark2")[2], brewer.pal(6, "Dark2")[3])) +
      scale_y_continuous(limits=c(0, 0.6))

fig3 <- rbind(lis[["p3"]], lis[["p7"]]) %>% 
   as.data.table() %>% 
   .[, modas := factor(modas, levels=c("hmC_OCR", "hmC_FULL"))]  %>% 
   ggplot(aes(elementas, value, fill=stadija)) +
    geom_boxplot(outlier.shape=NA, lwd=0.5, fatten = 0.5 ) +
    facet_wrap(~modas, nrow=1) +
    theme_Publication() +
      scale_fill_manual(values=c(brewer.pal(6, "Dark2")[1], brewer.pal(6, "Dark2")[2], brewer.pal(6, "Dark2")[3])) +
      scale_y_continuous(limits=c(0, 0.6))

fig4 <- rbind(lis[["p4"]], lis[["p8"]]) %>% 
   as.data.table() %>% 
   .[, modas := factor(modas, levels=c("hmC_OCR", "hmC_FULL"))]  %>% 
   ggplot(aes(elementas, value, fill=stadija)) +
    geom_boxplot(outlier.shape=NA, lwd=0.5, fatten = 0.5 ) +
    facet_wrap(~modas, nrow=1) +
    theme_Publication() +
      scale_fill_manual(values=c(brewer.pal(6, "Dark2")[1], brewer.pal(6, "Dark2")[2], brewer.pal(6, "Dark2")[3])) +
      scale_y_continuous(limits=c(0, 0.6))

p <- (fig1 |  fig2 | fig3 | fig4)
p <- p  + plot_layout(guides = "collect")
