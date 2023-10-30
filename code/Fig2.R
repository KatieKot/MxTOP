library(pacman)
p_load(data.table, dplyr, ggplot2, RColorBrewer, viridis, ggthemes, stringr, ggpubr, ComplexHeatmap, 
  colorRamp2, ggseqlogo, scales, ggrastr, GenomicRanges, ggrepel)
options(scipen=999)

# B 
# This figure is a part of Fig S3 F

# C
d <- readRDS("input/Fig2_C.RDS")
d %>% 
  as.data.table() %>% 
  .[, tipas := factor(tipas, levels=c("shallow_GC", "shallow_GCCG", "TT_GC", "TT_GCCG"), labels=c("idGC\nshallow", "idT\nshallow", "idGC", "idT"))] %>% 
  ggplot(aes(refas, new_svoris, fill=group_poster)) +
    geom_col(position="stack", width=.65) +
    facet_wrap(~tipas, ncol=1, strip.position="right") +
    coord_flip() +
    theme_Publication() +
    scale_fill_manual(values=cols_posterGroups3) +
    theme(legend.position="top", 
          panel.spacing.y = unit(0.1, "lines"),
          axis.title.y=element_blank()) +
    ylim(c(0, 0.91))

# D
d <- readRDS("input/Fig2_D.RDS")
rbind(
  d[["ATAC"]] %>% .[, Tech := "ATAC-seq"],
  d[["DNA"]] %>% .[, Tech := "DNase-seq"])  %>% 
  .[samplas %in% c("TT_GCCG", "TT_GC", "TT_shallow_GCCG", "TT_shallow_GC", "gDNA_GCCG", "gDNA_GC"), ] %>% 
  .[, .(samplas, Tech, p.value, estimate)] %>% 
  .[, samplas := factor(samplas, levels=c("TT_GCCG", "TT_GC", "Ts_GCCG", "Ts_GC", "gDNA_GCCG", "gDNA_GC"), labels=c("idT", "idGC", "idT\nshallow", "idGC\nshallow", "idT\ngDNA", "idGC\ngDNA" ))] %>% 
  .[, val := log2(estimate)] %>% 
  ggplot(aes(samplas, val, fill=Tech)) +
    geom_point(colour="black", shape=21) +
    theme_Publication() + 
    ylab("Odds ratio, log2") +
    geom_hline(yintercept=0, color = "red", size=0.25, linetype = "dashed") +
    scale_fill_brewer(palette="Dark2", name="") +
    coord_flip()    

# E
readRDS("input/Fig2_E.RDS") %>% 
  ggplot(aes(group_poster, met, fill=tipas)) +
    geom_boxplot(outlier.shape=NA) +
    theme_bw() +
    scale_fill_brewer(palette="Greys", name="Regions source") +
    facet_grid(~elementas, scales="free", space="free") +
    theme_Publication() +
    ylab("Methylation") +
    theme(legend.position="top",
          legend.key.size = unit(1, 'cm'),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank())     

# F
colours  <- c("hmCm_uCGm"="plum3", 
                "hmCp"=viridis(12)[5],
                "uCGp"=viridis(12)[11],
                "hmCp_uCGp"=viridis(12)[9]
                )   

readRDS("input/Fig2_F.RDS") %>% 
  ggplot(aes(uCG_fraction, metilinimas, fill=region_type, shape=enhancer_type, label=regionName)) +
    geom_smooth(aes(uCG_fraction, metilinimas, group="1"), method = "lm", color = "grey", se=FALSE, linetype="dashed", lwd=0.5) +
    geom_point() +
    theme_bw() +
    scale_fill_manual(values=colours) +
    geom_text_repel(size=1.5) + 
    stat_cor(aes(uCG_fraction, metilinimas, group="1"), label.x = 0.5, size=2) +
    scale_shape_manual(values = c(21, 25)) +
    theme_Publication()  

# G
readRDS("./input/Fig2_G.RDS") %>%      
  ggplot(aes(group, value)) +
    geom_boxplot(outlier.shape=NA) +
    facet_wrap(~variable, nrow=2) +
    xlab("Expression group") +
    ylab("Open fraction, %") +
    theme_Publication()