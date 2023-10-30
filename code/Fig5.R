library(pacman)
p_load(data.table, dplyr, ggplot2, RColorBrewer, viridis, ggthemes, stringr, ggpubr, ComplexHeatmap, 
  colorRamp2, ggseqlogo, scales, ggrastr, GenomicRanges, ggrepel, scales)
options(scipen=999)

# A
readRDS("input/Fig5_A.RDS") %>% 
    ggplot(aes(nr, value, fill=elementas)) +
    geom_col(, position="dodge") +
    theme_bw() +
    theme_Publication() +
    scale_fill_manual(values=c(brewer.pal(6, "Greys")[3], brewer.pal(6, "Greys")[5])) 

# B
# Created using R 4.0.1
tmp <- readRDS("input/Fig5_B.RDS")

tmp[[1]] + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey", size=0.5) + 
  geom_vline(xintercept = 0, linetype="dashed", color = "grey", size=0.5) +
  scale_fill_distiller(palette=4, direction=1, limits=c(0, 0.16)) 

tmp[[2]] + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey", size=0.5) + 
  geom_vline(xintercept = 0, linetype="dashed", color = "grey", size=0.5) + 
  scale_fill_distiller(palette=4, direction=1, limits=c(0, 0.16)) 

tmp[[3]] + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey", size=0.5) + 
  geom_vline(xintercept = 0, linetype="dashed", color = "grey", size=0.5) + 
  scale_fill_distiller(palette=4, direction=1, limits=c(0, 0.16)) 

tmp[[4]] +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "grey", size=0.5) + 
  geom_vline(xintercept = 0, linetype="dashed", color = "grey", size=0.5) 

# C
readRDS("input/Fig5_C.RDS") %>% 
  ggplot(., aes(variable, value, colour=variable.1, group=variable.1)) +
    geom_line(key_glyph = draw_key_label) +
    facet_wrap(~RNA_group, nrow=1) +
    theme_bw() +
    theme_Publication() +
    scale_colour_brewer(palette="Dark2")

# D
readRDS("input/Fig5_D1.RDS") %>% 
 ggplot(aes(quant, value)) +
  geom_boxplot(outlier.shape=NA, lwd=0.5, fatten=1) +
  facet_wrap(~subtype, scales="free", nrow=1) +
  coord_cartesian(ylim =c(0, 0.5)) +
  theme_Publication()

readRDS("input/Fig5_D2.RDS") %>% 
 ggplot(aes(quant, value)) +
  geom_boxplot(outlier.shape=NA, lwd=0.5, fatten=1) +
  facet_wrap(~subtype, scales="free", nrow=1) +
  coord_cartesian(ylim =c(0, 0.5)) +
  theme_Publication()

# E
readRDS("input/Fig5_E.RDS") %>% 
  ggplot(aes(bin, fr, colour=signal, group=signal)) +
    geom_line(key_glyph = draw_key_label) +
    geom_vline(xintercept = 10.5, linetype="dashed", color = "grey", size=0.5) + 
    facet_wrap(~exp+ stage, nrow=1) +
    theme_bw() +
    scale_colour_brewer(palette="Set1", direction=1)