library(pacman)
p_load(data.table, dplyr, ggplot2, RColorBrewer, viridis, ggthemes, stringr, ggpubr, ComplexHeatmap, 
  colorRamp2, ggseqlogo, scales, ggrastr, GenomicRanges, ggrepel, ggvenn, circlize, ggvenn)
options(scipen=999)
source("code/common.R")
# A
cols_posterGroups4  <- c("GC_only"="lightskyblue", 
                        "uCGhmC_minus"="plum3", 
                        "hmC_plius"=viridis(12)[5],
                        "uCG_plius"=viridis(12)[11],
                        "uCGhmC_plius"=viridis(12)[9]
                        )
readRDS("input/Fig6_A.RDS") %>% 
  .[, V1 := factor(V1, levels=c("GC_only", "uCG-_hmC-", "hmC+", "uCG+", "uCG+_hmC+"))] %>% 
  ggplot(aes(stage, N, fill=V1)) +
    geom_col() +
    theme_bw() +
    scale_fill_manual(values=cols_posterGroups) +
    scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) +
    theme_Publication()

# B
d <- readRDS("input/Fig6_B.RDS")
total <- d[, lapply(.SD, sum), by="stage", .SDcols="N"] %>% .[]

merge(d, total, by="stage") %>% 
  .[, frakcijos := N.x/N.y] %>% 
  .[, V1 := factor(V1, levels=c("On_gene", "by_gene", "proximal_intergenic", "distal_intergenic", "KITA"))] %>% 
  ggplot(aes(stage, frakcijos, fill=V1)) +
    geom_col(position="stack") +
    theme_bw() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 2)) +
    scale_fill_manual(values=alpha(c("On_gene"="grey90", "by_gene"="cadetblue4", "proximal_intergenic"="darkseagreen", "distal_intergenic"="grey50", "KITA"="grey10"), 0.5)) +
    theme_Publication()

# C
d <- readRDS("input/Fig6_C.RDS")
ggvenn(d[["S1denovo"]], text_size=2, 
      set_name_size=2, show_percentage = FALSE, 
      stroke_size=0.25, fill_color = c("blue", "yellow", "green", "red"), fill_alpha = 0.0) 

ggvenn(d[["S2denovo"]], text_size=2, 
      set_name_size=2, show_percentage = FALSE, 
      stroke_size=0.25, fill_color = c("blue", "yellow", "green", "red"), fill_alpha = 0.0) 

ggvenn(d[["S1S2denovo"]], text_size=2, 
      set_name_size=2, show_percentage = FALSE, 
      stroke_size=0.25, fill_color = c("blue", "yellow", "green", "red"), fill_alpha = 0.0) 

# D
readRDS("input/Fig6_D1.RDS") %>% 
  as.data.table() %>% 
  .[, RNA := factor(RNA, levels=c("RNA_inc", "RNA_dec"))] %>% 
  as.data.table() %>% 
  .[zymuo != "hmC"] %>% 
  ggplot(aes(variable, value, fill=zymuo, colour=zymuo)) +
    geom_linerange(aes(x=variable, ymin=0, ymax=value, colour=zymuo, fill=zymuo), position = position_dodge(width = 0.75), size = 0.5) +
    geom_point(position = position_dodge(width = 0.75), size=1) +
    theme_bw() +
    geom_hline(yintercept=0, linetype="dashed", colour="red", size=0.24) +
    facet_grid(vars(RNA), vars(tipas_anotacija)) +
    ylab("Fold change") +
    theme(axis.title.x=element_blank()) +
    scale_colour_manual(name = "type", values = c("open"=brewer.pal(5, "PuRd")[3], "uCG"=brewer.pal(5, "Greens")[3], "hmC"=brewer.pal(5, "Oranges")[3])) +
    scale_fill_manual(name = "type", values = c("open"=brewer.pal(5, "PuRd")[3], "uCG"=brewer.pal(5, "Greens")[3], "hmC"=brewer.pal(5, "Oranges")[3])) +
    coord_cartesian(ylim=c(-0.65, 0.5))

readRDS("input/Fig6_D2.RDS") %>% 
  as.data.table() %>% 
  .[, RNA := factor(RNA, levels=c("RNA_inc", "RNA_dec"))] %>% 
  as.data.table() %>% 
  .[zymuo == "hmC"] %>%
  .[variable == "LFC_S0S2", ] %>%  
  ggplot(aes(variable, value, fill=zymuo, colour=zymuo)) +
    geom_linerange(aes(x=variable, ymin=0, ymax=value, colour=zymuo, fill=zymuo), position = position_dodge(width = 0.75), size = 0.5) +
    geom_point(position = position_dodge(width = 0.75), size=1) +
    theme_bw() +
    geom_hline(yintercept=0, linetype="dashed", colour="red", size=0.24) +
    facet_grid(vars(RNA), vars(tipas_anotacija)) +
    ylab("Fold change") +
    theme(axis.title.x=element_blank()) +
    scale_colour_manual(name = "type", values = c("open"=brewer.pal(5, "PuRd")[3], "uCG"=brewer.pal(5, "Greens")[3], "hmC"=brewer.pal(5, "Oranges")[3])) +
    scale_fill_manual(name = "type", values = c("open"=brewer.pal(5, "PuRd")[3], "uCG"=brewer.pal(5, "Greens")[3], "hmC"=brewer.pal(5, "Oranges")[3])) +
    coord_cartesian(ylim=c(-1, 0.7))

a <- readRDS("../Individual_reps/output/code/Data_6D/Data_6D_background.RDS")
saveRDS(a, "input/Fig6_D3.RDS")
readRDS("input/Fig6_D3.RDS") %>% 
  as.data.table() %>% 
  .[, RNA := factor(RNA, levels=c("RNA_inc", "RNA_dec"))] %>% 
  as.data.table() %>% 
  .[zymuo != "hmC"] %>% 
  ggplot(aes(variable, value, fill=zymuo, colour=zymuo)) +
  geom_linerange(aes(x=variable, ymin=0, ymax=value, colour=zymuo, fill=zymuo), position = position_dodge(width = 0.75), size = 0.5) +
  geom_point(position = position_dodge(width = 0.75), size=1) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", colour="red", size=0.24) +
  facet_grid(vars(RNA), vars(tipas_anotacija)) +
  ylab("Fold change") +
  theme(axis.title.x=element_blank()) +
  scale_colour_manual(name = "type", values = c("open"=brewer.pal(5, "PuRd")[3], "uCG"=brewer.pal(5, "Greens")[3], "hmC"=brewer.pal(5, "Oranges")[3])) +
  scale_fill_manual(name = "type", values = c("open"=brewer.pal(5, "PuRd")[3], "uCG"=brewer.pal(5, "Greens")[3], "hmC"=brewer.pal(5, "Oranges")[3])) +
  coord_cartesian(ylim=c(-0.65, 0.5))

readRDS("input/Fig6_D4.RDS") %>% 
  as.data.table() %>% 
  .[, RNA := factor(RNA, levels=c("RNA_inc", "RNA_dec"))] %>% 
  as.data.table() %>% 
  .[zymuo == "hmC"] %>% 
  .[variable == "LFC_S0S2", ] %>%  
  ggplot(aes(variable, value, fill=zymuo, colour=zymuo)) +
  geom_linerange(aes(x=variable, ymin=0, ymax=value, colour=zymuo, fill=zymuo), position = position_dodge(width = 0.75), size = 0.5) +
  geom_point(position = position_dodge(width = 0.75), size=1) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", colour="red", size=0.24) +
  facet_grid(vars(RNA), vars(tipas_anotacija)) +
  ylab("Fold change") +
  theme(axis.title.x=element_blank()) +
  scale_colour_manual(name = "type", values = c("open"=brewer.pal(5, "PuRd")[3], "uCG"=brewer.pal(5, "Greens")[3], "hmC"=brewer.pal(5, "Oranges")[3])) +
  scale_fill_manual(name = "type", values = c("open"=brewer.pal(5, "PuRd")[3], "uCG"=brewer.pal(5, "Greens")[3], "hmC"=brewer.pal(5, "Oranges")[3])) +
  coord_cartesian(ylim=c(-1, 0.7))

# E
d <- readRDS("input/Fig6_E.RDS")
tmp_scaled_uCG <- d[["tmp_scaled_uCG"]] %>% .[tipas != "RNA_lateDec", ]
tmp_scaled_hmC <- d[["tmp_scaled_hmC"]] %>% .[tipas != "RNA_lateDec", ]
tmp_scaled_uCG_body <- d[["tmp_scaled_uCG_body"]] %>% .[tipas != "RNA_lateDec", ]
tmp_scaled_uCG_bodyFull <- d[["tmp_scaled_uCG_bodyFull"]] %>% .[tipas != "RNA_lateDec", ] 
tmp_scaled_hmC_body  <- d[["tmp_scaled_hmC_body"]] %>% .[tipas != "RNA_lateDec", ]
tmp_scaled_hmC_bodyFull <- d[["tmp_scaled_hmC_bodyFull"]] %>% .[tipas != "RNA_lateDec", ]
tmp_scaled_uCG$tipas <- factor(tmp_scaled_uCG$tipas, levels=c("RNA_inc", "RNA_earlyInc", "RNA_lateInc",  "RNA_dec", "RNA_earlyDec"))
heatmapWidth <- 18
heatmapHeight <- 4
colors <- structure(c(brewer.pal(3, "Reds")[3], brewer.pal(3, "Reds")[2], brewer.pal(3, "Reds")[1], brewer.pal(3, "Blues")[3], brewer.pal(3, "Blues")[2], brewer.pal(3, "Blues")[1]), names = c("RNA_inc", "RNA_earlyInc", "RNA_lateInc", "RNA_dec", "RNA_earlyDec", "RNA_lateDec"))
col_fun = colorRamp2(c(-1.5, 0, 1.5), viridis(3))

top_A = HeatmapAnnotation(                                            
                          stage = c("S0", "S1", "S2"), 
                          col = list(stage = c("S0" = brewer.pal(6, "Greys")[2], "S1" = brewer.pal(6, "Greys")[4], "S2" =brewer.pal(6, "Greys")[5] ) ), 
                          show_annotation_name = FALSE,
                          show_legend = TRUE,
                          border=TRUE,
                          simple_anno_size = unit(3, "mm"),
                          annotation_name_side = "left")

ht <- Heatmap(tmp_scaled_uCG[, -(1:2)], name="z-score", clustering_distance_rows = "spearman", column_title=NULL, cluster_column = FALSE, cluster_column_slices = FALSE, column_title_rot=90, row_title_rot=0, show_heatmap_legend = TRUE, width = unit(heatmapWidth, "mm"), height=unit(heatmapHeight, "cm"), top_annotation = top_A, use_raster=TRUE, col=col_fun) +       
      Heatmap(tmp_scaled_hmC[, -(1:2)], cluster_rows = TRUE, name="a", column_title=NULL, cluster_column = FALSE, cluster_column_slices = FALSE, column_title_rot=90, show_heatmap_legend = FALSE, width = unit(heatmapWidth, "mm"), height=unit(heatmapHeight, "cm"), top_annotation = top_A, use_raster=TRUE, col=col_fun) +
      Heatmap(tmp_scaled_uCG_body[, -(1:2)], name="b", column_title=NULL, cluster_column = FALSE, cluster_column_slices = FALSE, column_title_rot=90, show_heatmap_legend = FALSE, width = unit(heatmapWidth, "mm"), height=unit(heatmapHeight, "cm"), top_annotation = top_A, use_raster=TRUE, col=col_fun) + 
      Heatmap(tmp_scaled_uCG_bodyFull[, -(1:2)], name=".", column_title=NULL, cluster_column = FALSE, cluster_column_slices = FALSE, column_title_rot=90, show_heatmap_legend = FALSE, width = unit(heatmapWidth, "mm"), height=unit(heatmapHeight, "cm"), top_annotation = top_A, use_raster=TRUE, col=col_fun) + 
      Heatmap(tmp_scaled_hmC_body[, -(1:2)], name="d", column_title=NULL, cluster_column = FALSE, cluster_column_slices = FALSE, column_title_rot=90, show_heatmap_legend = FALSE, width = unit(heatmapWidth, "mm"), height=unit(heatmapHeight, "cm"), top_annotation = top_A, use_raster=TRUE, col=col_fun) + 
      Heatmap(tmp_scaled_hmC_bodyFull[, -(1:2)], name="e", column_title=NULL, cluster_column = FALSE, cluster_column_slices = FALSE, column_title_rot=90, show_heatmap_legend = FALSE, width = unit(heatmapWidth, "mm"), height=unit(heatmapHeight, "cm"), top_annotation = top_A, use_raster=TRUE, col=col_fun) +
      Heatmap(tmp_scaled_uCG$tipas, name=NULL, col=colors, width=unit(3, "mm")) 

draw(ht, row_split = tmp_scaled_uCG$tipas, cluster_column_slices = NULL, cluster_row_slices = FALSE, cluster_column = NULL, column_title_gp=gpar(fontsize = 7))