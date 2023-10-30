library(pacman)
p_load(data.table, dplyr, ggplot2, RColorBrewer, viridis, ggthemes, stringr, ggpubr, ComplexHeatmap, 
  colorRamp2, ggseqlogo, scales, ggrastr, GenomicRanges, ggrepel)
options(scipen=999)

# A
# B
w <- 30
tmp <-  readRDS(paste0("output/EnriHeatOpenGroupsVSRNAgroups/NormMatrix_open_regions", w, ".RDS"))
rez <- rbind(
colMeans(tmp[["IncIncS0"]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S0"] %>% 
    .[, type := "open_inc"] %>% 
    .[, mod := "atvirumas"], 
colMeans(tmp[["IncIncS1"]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S1"] %>% 
    .[, type := "open_inc"] %>% 
    .[, mod := "atvirumas"], 

colMeans(tmp[["IncIncS2"]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S2"] %>% 
    .[, type := "open_inc"] %>% 
    .[, mod := "atvirumas"], 

colMeans(tmp[["IncdecS0"]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S0"] %>% 
    .[, type := "open_dec"] %>% 
    .[, mod := "atvirumas"], 

colMeans(tmp[["IncdecS1"]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S1"] %>% 
    .[, type := "open_dec"] %>% 
    .[, mod := "atvirumas"], 

colMeans(tmp[["IncdecS2"]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S2"] %>% 
    .[, type := "open_dec"] %>% 
    .[, mod := "atvirumas"], 

colMeans(tmp[["DecIncS0"]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S0"] %>% 
    .[, type := "closes_inc"] %>% 
    .[, mod := "atvirumas"] ,

colMeans(tmp[["DecIncS1"]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S1"] %>% 
    .[, type := "closes_inc"] %>% 
    .[, mod := "atvirumas"], 

colMeans(tmp[["DecIncS2"]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S2"] %>% 
    .[, type := "closes_inc"] %>% 
    .[, mod := "atvirumas"] ,

colMeans(tmp[["DecdecS0"]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S0"] %>% 
    .[, type := "closes_dec"] %>% 
    .[, mod := "atvirumas"] ,

colMeans(tmp[["DecdecS1"]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S1"] %>% 
    .[, type := "closes_dec"] %>% 
    .[, mod := "atvirumas"] ,

colMeans(tmp[["DecdecS2"]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S2"] %>% 
    .[, type := "closes_dec"] %>% 
    .[, mod := "atvirumas"] ) %>% 
     as.data.table() %>% 
   .[, rn := factor(rn, levels=colnames(tmp[[1]][]))]  %>% 
   .[, type := factor(type, levels=c("open_inc", "open_dec", "closes_inc", "closes_dec"), labels=c("open_inc", "open_dec", "closes_inc", "closes_dec"))] 



prof <- ggplot(data=rez)  + 
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

prof
svg(paste0(outdatadir, "fig_profiliai_atvirumai.svg"), width=1.5, heigh=1.5)
   prof  + theme(legend.key.size = unit(0.5, 'cm'),
        strip.text.y = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        strip.text.x = element_blank(),
        legend.position="bottom")
dev.off() 




w <- 300
opens_increase <- readRDS(paste0("output/EnriHeatOpenGroupsVSRNAgroups/open_increase_frac_Regions_", w, "_normMatrix.RDS"))
opens_decrease <- readRDS(paste0("output/EnriHeatOpenGroupsVSRNAgroups/open_decrease_frac_Regions_", w, "_normMatrix.RDS"))
closes_increase <- readRDS(paste0("output/EnriHeatOpenGroupsVSRNAgroups/closes_increase_frac_Regions_", w, "_normMatrix.RDS"))
closes_decrease <- readRDS(paste0("output/EnriHeatOpenGroupsVSRNAgroups/closes_decrease_frac_Regions_", w, "_normMatrix.RDS"))

get_prof <- function(x, k, l, m) {
tmp <- rbind(
  colMeans(opens_increase[[k]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S0"] %>% 
    .[, type := "open_inc"] %>% 
    .[, mod := x],
  colMeans(opens_increase[[l]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S1"] %>% 
    .[, type := "open_inc"] %>% 
    .[, mod := x],
  colMeans(opens_increase[[m]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S2"] %>% 
    .[, type := "open_inc"] %>% 
    .[, mod := x], 
  colMeans(opens_decrease[[k]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S0"] %>% 
    .[, type := "open_dec"] %>% 
    .[, mod := x],
  colMeans(opens_decrease[[l]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S1"] %>% 
    .[, type := "open_dec"] %>% 
    .[, mod := x],
  colMeans(opens_decrease[[m]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S2"] %>% 
    .[, type := "open_dec"] %>% 
    .[, mod := x],    
  colMeans(closes_increase[[k]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S0"] %>% 
    .[, type := "closes_inc"] %>% 
    .[, mod := x],
  colMeans(closes_increase[[l]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S1"] %>% 
    .[, type := "closes_inc"] %>% 
    .[, mod := x],
  colMeans(closes_increase[[m]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S2"] %>% 
    .[, type := "closes_inc"] %>% 
    .[, mod := x], 
  colMeans(closes_decrease[[k]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S0"] %>% 
    .[, type := "closes_dec"] %>% 
    .[, mod := x],
  colMeans(closes_decrease[[l]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S1"] %>% 
    .[, type := "closes_dec"] %>% 
    .[, mod := x],
  colMeans(closes_decrease[[m]][]) %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    setnames(".", "fraction") %>% 
    .[, stage := "S2"] %>% 
    .[, type := "closes_dec"] %>% 
    .[, mod := x]       ) %>% 
  as.data.table() %>% 
  .[, rn := factor(rn, levels=colnames(opens_increase[[1]][]))] 

prof <- tmp %>% 
  .[, type := factor(type, levels=c("open_inc", "open_dec", "closes_inc", "closes_dec"))] %>% 
    ggplot(aes(colour=stage))  + 
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
  
return(prof) } 

p1 <- get_prof("hmC", 1, 2, 3)
p2 <- get_prof("uCG", 4, 5, 6)
#p3 <- get_prof("idT", 7, 8, 9)

p1
p2

#p3 <- p3 + theme(legend.key.size = unit(0.5, 'cm'),
#        strip.text.y = element_blank(),
#        axis.title=element_blank(),
#        axis.text=element_blank(),
#        strip.text.x = element_blank(),
#        legend.position="bottom")

svg(paste0(outdatadir, "hmC_profile.svg"), width=1.5, heigh=1.5)
   p1 + theme(legend.key.size = unit(0.5, 'cm'),
        strip.text.y = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        strip.text.x = element_blank(),
        legend.position="bottom")
dev.off()  


svg(paste0(outdatadir, "uCG_profile.svg"), width=1.5, heigh=1.5)
   p2 + theme(legend.key.size = unit(0.5, 'cm'),
        strip.text.y = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        strip.text.x = element_blank(),
        legend.position="bottom")
dev.off()  


opens_increase <- readRDS(paste0("output/EnriHeatOpenGroupsVSRNAgroups/open_increase_frac_total_", w, "_normMatrix.RDS"))
opens_decrease <- readRDS(paste0("output/EnriHeatOpenGroupsVSRNAgroups/open_decrease_frac_total_", w, "_normMatrix.RDS"))
closes_increase <- readRDS(paste0("output/EnriHeatOpenGroupsVSRNAgroups/closes_increase_frac_total_", w, "_normMatrix.RDS"))
closes_decrease <- readRDS(paste0("output/EnriHeatOpenGroupsVSRNAgroups/closes_decrease_frac_total_", w, "_normMatrix.RDS"))
p1_full <- get_prof("hmC", 1, 2, 3)
p1_full

svg(paste0(outdatadir, "hmC_full.svg"), width=1.5, heigh=1.5)
   p1_full + theme(legend.key.size = unit(0.5, 'cm'),
        strip.text.y = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        strip.text.x = element_blank(),
        legend.position="bottom")
dev.off() 





# C

d <- readRDS("output/ComplexHeat_part_fullSignal/ComplexHeatmap_data.RDS")

p5 <- d[gene_id %in% open_inc, ] %>% 
  as.data.table() %>% 
  .[, .(gb_TT_S0_GCCG_hmC, gb_TT_S1_GCCG_hmC, gb_TT_S2_GCCG_hmC)] %>% 
  melt() %>% 
  .[, variable := as.character(variable)] %>% 
  .[, stadija := str_extract(variable, "S.")] %>% 
  .[, elementas := str_extract(variable, "^..") %>% factor(., levels=c("pr", "gb"))] %>% 
  .[, modifikacija := str_extract(variable, "_...$") %>% gsub("_", "", .) %>% factor(., levels=c("idT", "uCG", "hmC"))] %>%
  .[, potipis := "FULL"] %>% 
  .[, modas := paste0(modifikacija, "_", potipis)]

p6 <- d[gene_id %in% open_dec, ] %>% 
  as.data.table() %>% 
  .[, .(gb_TT_S0_GCCG_hmC, gb_TT_S1_GCCG_hmC, gb_TT_S2_GCCG_hmC)] %>% 
  melt() %>% 
  .[, variable := as.character(variable)] %>% 
  .[, stadija := str_extract(variable, "S.")] %>% 
  .[, elementas := str_extract(variable, "^..") %>% factor(., levels=c("pr", "gb"))] %>% 
  .[, modifikacija := str_extract(variable, "_...$") %>% gsub("_", "", .) %>% factor(., levels=c("idT", "uCG", "hmC"))] %>%
  .[, potipis := "FULL"] %>% 
  .[, modas := paste0(modifikacija, "_", potipis)]

p7 <- d[gene_id %in% closes_inc, ] %>% 
  as.data.table() %>% 
  .[, .(gb_TT_S0_GCCG_hmC, gb_TT_S1_GCCG_hmC, gb_TT_S2_GCCG_hmC)] %>% 
  melt() %>% 
  .[, variable := as.character(variable)] %>% 
  .[, stadija := str_extract(variable, "S.")] %>% 
  .[, elementas := str_extract(variable, "^..") %>% factor(., levels=c("pr", "gb"))] %>% 
  .[, modifikacija := str_extract(variable, "_...$") %>% gsub("_", "", .) %>% factor(., levels=c("idT", "uCG", "hmC"))] %>%
  .[, potipis := "FULL"] %>% 
  .[, modas := paste0(modifikacija, "_", potipis)]

p8 <- d[gene_id %in% closes_dec, ] %>% 
  as.data.table() %>% 
  .[, .(gb_TT_S0_GCCG_hmC, gb_TT_S1_GCCG_hmC, gb_TT_S2_GCCG_hmC)] %>% 
  melt() %>% 
  .[, variable := as.character(variable)] %>% 
  .[, stadija := str_extract(variable, "S.")] %>% 
  .[, elementas := str_extract(variable, "^..") %>% factor(., levels=c("pr", "gb"))] %>% 
  .[, modifikacija := str_extract(variable, "_...$") %>% gsub("_", "", .) %>% factor(., levels=c("idT", "uCG", "hmC"))] %>%
  .[, potipis := "FULL"] %>% 
  .[, modas := paste0(modifikacija, "_", potipis)]


fig1 <- rbind(p1, p5) %>% 
   as.data.table() %>% 
   .[, modas := factor(modas, levels=c("hmC_OCR", "hmC_FULL"))]  %>% 
   ggplot(aes(elementas, value, fill=stadija)) +
    geom_boxplot(outlier.shape=NA, lwd=0.5, fatten = 0.5 ) +
    facet_wrap(~modas, nrow=1) +
    theme_Publication() +
      scale_fill_manual(values=c(brewer.pal(6, "Dark2")[1], brewer.pal(6, "Dark2")[2], brewer.pal(6, "Dark2")[3])) +
      scale_y_continuous(limits=c(0, 0.6))

fig2 <- rbind(p2, p6) %>% 
   as.data.table() %>% 
   .[, modas := factor(modas, levels=c("hmC_OCR", "hmC_FULL"))]  %>% 
   ggplot(aes(elementas, value, fill=stadija )) +
    geom_boxplot(outlier.shape=NA, lwd=0.5, fatten = 0.5) +
    facet_wrap(~modas, nrow=1) +
    theme_Publication() +
      scale_fill_manual(values=c(brewer.pal(6, "Dark2")[1], brewer.pal(6, "Dark2")[2], brewer.pal(6, "Dark2")[3])) +
      scale_y_continuous(limits=c(0, 0.6))

fig3 <- rbind(p3, p7) %>% 
   as.data.table() %>% 
   .[, modas := factor(modas, levels=c("hmC_OCR", "hmC_FULL"))]  %>% 
   ggplot(aes(elementas, value, fill=stadija)) +
    geom_boxplot(outlier.shape=NA, lwd=0.5, fatten = 0.5 ) +
    facet_wrap(~modas, nrow=1) +
    theme_Publication() +
      scale_fill_manual(values=c(brewer.pal(6, "Dark2")[1], brewer.pal(6, "Dark2")[2], brewer.pal(6, "Dark2")[3])) +
      scale_y_continuous(limits=c(0, 0.6))

fig4 <- rbind(p4, p8) %>% 
   as.data.table() %>% 
   .[, modas := factor(modas, levels=c("hmC_OCR", "hmC_FULL"))]  %>% 
   ggplot(aes(elementas, value, fill=stadija)) +
    geom_boxplot(outlier.shape=NA, lwd=0.5, fatten = 0.5 ) +
    facet_wrap(~modas, nrow=1) +
    theme_Publication() +
      scale_fill_manual(values=c(brewer.pal(6, "Dark2")[1], brewer.pal(6, "Dark2")[2], brewer.pal(6, "Dark2")[3])) +
      scale_y_continuous(limits=c(0, 0.6))


fig1
fig2
fig3
fig4

p <- (fig1 |  fig2 | fig3 | fig4)
p <- p  + plot_layout(guides = "collect")

svg(paste0(outdatadir, "fig_boxplot_Combinuotas_hmCBody_suskaiciais.svg"), width=6, heigh=2)
   p
dev.off() 

fig1 <- fig1 + theme(legend.key.size = unit(1, 'cm'),
        strip.text.y = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        strip.text.x = element_blank(),
        legend.position="none")

fig2 <- fig2 + theme(legend.key.size = unit(1, 'cm'),
        strip.text.y = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        strip.text.x = element_blank(),
        legend.position="none")

fig3 <- fig3 + theme(legend.key.size = unit(1, 'cm'),
        strip.text.y = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        strip.text.x = element_blank(),
        legend.position="none")


fig4 <- fig4 + theme(legend.key.size = unit(1, 'cm'),
        strip.text.y = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        strip.text.x = element_blank(),
        legend.position="none")

p <- (fig1 | fig2 | fig3 | fig4)

p <- p  + plot_layout(guides = "collect")

svg(paste0(outdatadir, "fig_boxplot_Combinuotas_hmCBody.svg"), width=6, heigh=1)
   p
dev.off() 


# D