library(pacman)
p_load(data.table, dplyr, ggplot2, RColorBrewer, viridis, ggthemes, stringr, ggpubr, ComplexHeatmap, 
  colorRamp2, ggseqlogo, scales, ggrastr, GenomicRanges, ggrepel)
options(scipen=999)

# C
readRDS("input/Fig1_C.RDS") %>% 
  ggplot(aes(x="", y=value, fill=variable)) +
    geom_bar(stat="identity", width=1, colour="black") +
    geom_text(data=subset(d_targets, value>0.9), 
                    aes(y=value-0.5, label=paste0(round(value, 4)*100, "%")),
                    arrow = arrow(length = unit(0.015, "npc")),
                    segment.linetype = 1,
                    size=7
                  ) +
    coord_polar("y", start=0.1) +
    facet_wrap(~target, nrow=1) + 
    scale_fill_brewer(palette="Dark2", name=NULL) +
    theme_void() +
    geom_col(aes(x=0,y=0), width=0.25) + 
    theme(strip.text.x = element_text(size=rel(1.2), color="black"),
          legend.key.size = unit(1, 'cm'),
          legend.position="top"
    )

# D
readRDS("input/Fig1_D.RDS") %>% 
  .[, fragment := "Fragment #1"] %>% 
  .[, sample := factor(sample, levels=c("hmC_0", "hmC_10", "hmC_20", "hmC_50", "hmC_100"), labels=c("uCG", "1:10", "1:5", "1:1", "5hmC"))] %>% 
  .[, CG_NR := factor(CG_NR, levels=c("CCGG", "CG", "GC"))] %>%  
    .[is.na(T), T := 0] %>% 
    .[is.na(A), A := 0] %>%
    .[, A := A + 1] %>% 
    .[, T := T + 1] %>% 
    .[, Fraction := log(T/A)] %>% 
  ggplot(aes(strand, Fraction, fill=sample)) + 
      geom_point(aes(fill=sample),  colour="black", shape=21, size=1.5, position = position_dodge2(width = 0.5)) +
      facet_wrap(~fragment+CG_NR) +
      theme_bw() +
      ylab("5hmC/5uCG, log") +
      ylim(c(-5.5, 5)) +
      scale_fill_brewer(palette = "RdPu", direction=1) +
      theme_Publication() 