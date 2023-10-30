library(pacman)
p_load(data.table, dplyr, ggplot2, RColorBrewer, viridis, ggthemes, stringr, ggpubr, ComplexHeatmap, 
  colorRamp2, ggseqlogo, scales, ggrastr, GenomicRanges, ggrepel)
options(scipen=999)

###############################################################################
###############################################################################
############################ Sup. Figure 3
###############################################################################
###############################################################################
# A
d <- rbind(
  (d1 <- readRDS(paste0("input/FigS3_A1.RDS")) ),
  (d2 <- readRDS(paste0("input/FigS3_A1.RDS")) %>% 
  .[sample %in% c("mESC_R1", "mESC_R2"), ] %>% 
  .[sample == "mESC_R1", sample := "dESC_R1" ] %>% 
  .[sample == "mESC_R2", sample := "dESC_R2" ] )) %>% 
  .[, sample := factor(sample, levels=c("gDNA_R1", "gDNA_R2", "mESC_R1", "mESC_R2", "dESC_R1", "dESC_R2"))] 

ggplot(d, aes(sample, value, fill=Adapter, label=labels)) +
  geom_col(position = position_stack(reverse = TRUE), colour="black") +
  theme_bw() +
  facet_wrap(~Target, nrow=1) +
  theme_Publication() +
  scale_fill_manual(values=c("#ff5555ff", "#0067ffff")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.title.x = element_blank(),
        legend.position="top",
        legend.key.size = unit(1.5,"line") ,
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
        axis.title.y=element_blank())      

# B
d1 <- readRDS(paste0("input/FigS3_B1.RDS"))
d2 <- readRDS(paste0("input/FigS3_B2.RDS"))
d3 <- readRDS(paste0("input/FigS3_B3.RDS"))

listas <- list(nuclei=d3, gDNA=d2, control=d1)

ggseqlogo(listas, met="prob", colour="black", seq_type='dna') +
  scale_x_continuous(breaks=1:3, labels=-1:1) +
  theme_Publication() +
  facet_wrap(~seq_group) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(size = rel(1.5), color="black")) +
  scale_y_continuous(breaks=pretty_breaks(n=3))

# C
d <- readRDS("input/FigS3_C.RDS") %>% 
  .[samplas %in% c("GCCG_75", "GC_75")] %>% 
  .[]

ggplot(d, aes(samplas, plociai)) +
  geom_boxplot(outlier.shape=NA) +
  coord_cartesian(ylim=c(0, 750)) +
  theme_Publication()

# D
d <- readRDS("input/FigS3_D.RDS")
d <- reshape(d, idvar="i", timevar="j", direction="wide")
colnames(d) <- gsub("rez.", "", colnames(d))
rownames(d) <- d$i
d <- d[, -1]

d <- d[, c("TD_CB_GCCG", "TD_CB_GC", "TD_CS_GCCG", "TD_CS_GC", "shallow_GCCG",  "shallow_GC", "deep_GCCG", "deep_GC")]
d <- d[c("TD_CB_GCCG", "TD_CB_GC", "TD_CS_GCCG", "TD_CS_GC", "shallow_GCCG",  "shallow_GC", "deep_GCCG", "deep_GC"), ]

gyliai <- str_split(colnames(d), "_") %>% sapply(., `[`, 1) %>% factor(., levels=c("TD", "deep", "shallow"), labels=c("shallow", "deep", "shallow")) %>% as.character.factor()
type <- str_split(colnames(d), "_") %>% sapply(., length) %>% factor(., levels=c(3, 2), labels=c("2X", "3X")) %>% as.character.factor()
targets <- str_extract(colnames(d), "_[GC]+$") %>% gsub("_", "", .) %>% factor(., levels=c("GC", "GCCG"), labels=c("idGC", "idT")) %>% as.character.factor()
col_fun <- colorRamp2(c(0, 0.5, 1), viridis(3))

column_ha = HeatmapAnnotation(targets = targets, 
                              depth = gyliai,
                              type = type, 
                              gp = gpar(col = "black"),
                              simple_anno_size = unit(0.2, "cm"),
                              show_legend = FALSE,
                              col = list(type=c("2X"=viridis(20)[19], "3X"=viridis(20)[17]), 
                                         depth=c("shallow"=viridis(20)[13], "deep"=viridis(20)[9]),
                                         targets=c("idT"=viridis(20)[6], "idGC"=viridis(20)[2])
                                         )
                              )
row_ha = rowAnnotation(type = type, 
                            depth = gyliai,
                            targets = targets, 
                            gp = gpar(col = "black"),
                            show_legend = FALSE, 
                            simple_anno_size = unit(0.2, "cm"),
                              col = list(type=c("2X"=viridis(20)[19], "3X"=viridis(20)[17]), 
                                         depth=c("shallow"=viridis(20)[13], "deep"=viridis(20)[9]),
                                         targets=c("idT"=viridis(20)[6], "idGC"=viridis(20)[2])
                                        
                                         )
                              )     
Heatmap(d, col = col_fun,  
      cluster_rows = FALSE,
      rect_gp = gpar(col = "black", lwd = 0.5), 
      cluster_columns = FALSE,
      top_annotation = column_ha,
      right_annotation = row_ha,
      show_column_dend = FALSE,
      show_row_dend = FALSE,
      show_row_names = FALSE, 
      show_column_names = FALSE,
      show_heatmap_legend = FALSE,
      name="Overlap",
      cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
         grid.text((round(d[i, j], 2)*1), x, y,  gp=gpar(fontsize = 6))  
    })


# E
d <- readRDS("input/FigS3_E.RDS")
d <- d$raw_counts %>% 
  as.data.table() %>% 
  .[, facete := "A"] %>% 
  .[element %in% c("Gene", "Intergenic"), facete := "B"] %>% 
  .[]

ggplot(d, aes(element, amount, fill=group)) +
    geom_col(position="stack", width=.75) +
    theme_bw() +  
    scale_fill_manual(values=cols_posterGroups) + 
    facet_wrap(~facete, scales="free") +
    theme_Publication()

# F
regions <- readRDS(paste0("input/FigS3_F1.RDS")) %>% 
  .[["TT_S0"]] %>% 
  as.data.table() %>% 
  .[, .(nCG, nGC, group_poster, uGC, uCG)] %>% 
  .[nGC > 300, nGC := 300] %>% 
  .[nCG > 200, nCG := 200] %>% 
  .[, uCG_grupe := "uCG 0"] %>% 
  .[uCG>0, uCG_grupe := "uCG 1-2"] %>% 
  .[uCG>2, uCG_grupe := "uCG 3-5"] %>% 
  .[uCG>5, uCG_grupe := "uCG 6-9"] %>% 
  .[uCG>9, uCG_grupe := "uCG 10.."] %>% 
  .[, uCG_grupe := factor(uCG_grupe, levels=c("uCG 0", "uCG 1-2", "uCG 3-5", "uCG 6-9", "uCG 10.."))] %>%  
  .[uGC>=0, idGC_grupe := "idGC 0-9"] %>% 
  .[uGC>9, idGC_grupe := "idGC 10-49"] %>% 
  .[uGC>49, idGC_grupe := "idGC 50-99"] %>% 
  .[uGC>99, idGC_grupe := "idGC 100.."] %>% 
  .[, idGC_grupe := factor(idGC_grupe, levels=c("idGC 0-9", "idGC 10-49", "idGC 50-99", "idGC 100.."))] %>% 
  .[, group_poster := factor(group_poster, 
          levels=c("GC_only", "uCG-_hmC-", "hmC+", "uCG+", "uCG+_hmC+"),
          labels=c("Only GC", "No 5uCG/5hmC", "5hmC", "5uCG", "5hmC&5uCG"))] %>% 
  .[, alpha := 0] %>% 
  .[uCG_grupe == "uCG 0", alpha := 0.25] %>%                   
  .[uCG_grupe == "uCG 1-2", alpha := 0.05] %>%                   
  .[uCG_grupe == "uCG 3-5", alpha := 0.1] %>%                   
  .[uCG_grupe == "uCG 6-9", alpha := 0.1] %>%                   
  .[uCG_grupe == "uCG 10..", alpha := 0.02] 

cols_scatter  <- c("uCG 0"=brewer.pal(3, "Greys")[2], 
                        "uCG 1-2"=brewer.pal(7, "Dark2")[1], 
                        "uCG 3-5"=brewer.pal(7, "Dark2")[2],
                        "uCG 6-9"=brewer.pal(7, "Dark2")[3],
                        "uCG 10.."=brewer.pal(7, "Dark2")[6]
                        )    

ggplot(regions, aes(nGC, nCG, colour=uCG_grupe, alpha=alpha)) +
  geom_point(size=0.05) +
  theme_Publication() +
  xlab("Genomic GC number") +
  ylab("Genomic CG number") +
  facet_grid(vars(idGC_grupe), vars(group_poster), scales="free") +
  theme(legend.key.size = unit(1, 'cm'),
        strip.text.y = element_text(angle = 90),
          legend.position="top") +
  scale_colour_manual(values=cols_scatter, name="Number of uCG") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size=3)), alpha=FALSE)  

regions <- readRDS(paste0("input/FigS3_F1.RDS")) %>% 
  .[["TT_S0"]] %>% 
  as.data.table() %>% 
  .[, .(nCG, nGC, group_poster, uGC, hmC)] %>% 
  .[nGC > 300, nGC := 300] %>% 
  .[nCG > 200, nCG := 200] %>% 
  .[, hmC_grupe := "hmC 0"] %>% 
  .[hmC>0, hmC_grupe := "hmC 1-2"] %>% 
  .[hmC>2, hmC_grupe := "hmC 3-5"] %>% 
  .[hmC>5, hmC_grupe := "hmC 6-9"] %>% 
  .[hmC>9, hmC_grupe := "hmC 10.."] %>% 
  .[, hmC_grupe := factor(hmC_grupe, levels=c("hmC 0", "hmC 1-2", "hmC 3-5", "hmC 6-9", "hmC 10.."))] %>%  
 #  .[, idGC_grupe := "idGC 0"] %>% 
  .[uGC>=0, idGC_grupe := "idGC 0-9"] %>% 
  .[uGC>9, idGC_grupe := "idGC 10-49"] %>% 
  .[uGC>49, idGC_grupe := "idGC 50-99"] %>% 
  .[uGC>99, idGC_grupe := "idGC 100.."] %>% 
  .[, idGC_grupe := factor(idGC_grupe, levels=c("idGC 0-9", "idGC 10-49", "idGC 50-99", "idGC 100.."))] %>% 
  .[, group_poster := factor(group_poster, 
          levels=c("GC_only", "uCG-_hmC-", "hmC+", "uCG+", "uCG+_hmC+"),
          labels=c("Only GC", "No 5uCG/5hmC", "5hmC", "5uCG", "5hmC&5uCG"))] %>% 
  .[, alpha := 0] %>% 
  .[hmC_grupe == "hmC 0", alpha := 0.25] %>%                   
  .[hmC_grupe == "hmC 1-2", alpha := 0.05] %>%                   
  .[hmC_grupe == "hmC 3-5", alpha := 0.1] %>%                   
  .[hmC_grupe == "hmC 6-9", alpha := 0.1] %>%                   
  .[hmC_grupe == "hmC 10..", alpha := 0.02] 

cols_scatter  <- c("hmC 0"=brewer.pal(3, "Greys")[2], 
                        "hmC 1-2"=brewer.pal(7, "Dark2")[1], 
                        "hmC 3-5"=brewer.pal(7, "Dark2")[2],
                        "hmC 6-9"=brewer.pal(7, "Dark2")[3],
                        "hmC 10.."=brewer.pal(7, "Dark2")[6]
                        )    

ggplot(regions, aes(nGC, nCG, colour=hmC_grupe, alpha=alpha)) +
  geom_point(size=0.05) +
  theme_Publication() +
  xlab("Genomic GC number") +
  ylab("Genomic CG number") +
  facet_grid(vars(idGC_grupe), vars(group_poster), scales="free") +
  theme(legend.key.size = unit(1, 'cm'),
        strip.text.y = element_text(angle = 90),
          legend.position="top") +
  scale_colour_manual(values=cols_scatter, name="Number of hmC") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size=3)), alpha=FALSE)  


###############################################################################
###############################################################################
############################ Sup. Figure 4
###############################################################################
###############################################################################
#A 
d_reg <- rbind(readRDS("input/FigS4_A1.RDS") %>% .[, elementas := gsub("serum_", "", elementas)] %>% .[!(elementas %in% c("Promoter-nonactive", "Promoter-active")), ], 
               readRDS("input/FigS4_A2.RDS") )
tvarka <- d_reg[signalas == "uCG+", ] %>% 
  setorder(., -estimate) %>% 
  .[, elementas] 


d_reg %>% 
   .[pval <= 0.05, ] %>% 
  .[, elementas := factor(elementas, levels=tvarka)] %>% 
  .[, estimate := log2(estimate)] %>% 
  ggplot(aes(elementas, estimate, fill=signalas)) +
    geom_jitter(shape=21, size=1.5, width=0.1) +
    theme_bw() +
    theme_Publication() + 
    geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
    theme(legend.key.size = unit(1, 'cm')) +
    scale_fill_manual(values = cols_posterGroups) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))        

rbind(readRDS("input/FigS4_A3.RDS") %>% .[, elementas := gsub("serum_", "", elementas)]  %>% .[!(elementas %in% c("Promoter-nonactive", "Promoter-active")), ],
      readRDS("input/FigS4_A4.RDS")) %>% 
  .[pval <= 0.05, ] %>% 
  .[, elementas := factor(elementas, levels=tvarka)] %>% 
  .[, estimate := log2(estimate)] %>% 
  ggplot(aes(elementas, estimate, fill=signalas)) +
    geom_jitter(shape=21, size=1.5, width=0.1) +
    theme_bw() +
    theme_Publication() + 
    geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
    theme(legend.key.size = unit(1, 'cm'), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) )

# B
nbins=15
d <- readRDS("input/FigS4_B1.RDS")
d1 <- readRDS("input/FigS4_B2.RDS")
viskas_viskas_GC <- d[tipas == "viskas" & targets == "viskas", ]
viskas_idCGplius_GC <- d[tipas == "viskas" & targets == "idCG_plius", ] 
viskas_viskas_GCCG <- d1[tipas == "viskas" & targets == "viskas", ]
viskas_idCGplius_GCCG <- d1[tipas == "viskas" & targets == "idCG_plius", ] 

rbind(viskas_viskas_GC %>% .[, facete := "A"], 
      viskas_idCGplius_GC %>% .[, facete := "B"],
      viskas_viskas_GCCG %>% .[, facete := "C"], 
      viskas_idCGplius_GCCG %>% .[, facete := "D"]) %>% 
  ggplot(aes(sqrt(mod_density2), met)) +
        geom_point(size=0.1, alpha=0.1, colour="grey") +
        geom_density_2d(aes(color = ..level..), bins=nbins, linewidth=0.25) +
        scale_color_viridis_c(option="C") +
        theme_Publication() +
        coord_cartesian(xlim=c(0, 3)) +
        facet_wrap(~facete, nrow=1) +
        theme(legend.position="bottom")

# C
d <- readRDS("input/FigS4_C.RDS")
ggplot(d, aes(methylation_group, uCG_fr)) +
   geom_violin(position=position_dodge(0.25), lwd=0.5)  +
    geom_boxplot(outlier.shape=NA, width=0.1, position = position_dodge(width = 0.9), lwd=0.1) +
    theme_bw() +
    facet_wrap(~elem, nrow=1) +
    theme_Publication() +
    stat_summary(fun=mean, geom='point', shape=20, size=1, colour="black")
    
# D
d <- readRDS("input/FigS4_D.RDS")

d %>%
  .[tipas != "visi", ] %>%  
  .[nCG > 3, ] %>% 
  ggplot(aes(x=met, colour=tipas)) +
    geom_density(key_glyph = "timeseries") +
    theme_bw() +
    scale_colour_brewer(palette="Dark2") +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    theme_Publication()

# E
dTAB <- readRDS("input/FigS4_E1.RDS") %>% 
    as.data.table %>%
    .[, frakcija := TAB_seq/nCG] %>% 
    .[, tipas := "TAB-seq"] %>% 
    .[, .(grupe, frakcija, nCG, tipas)] %>% 
    .[]

dMxTOP <- readRDS("input/FigS4_E2.RDS") %>% 
    .[, frakcija := hmC_fraction] %>% 
    .[, grupe := group_poster] %>% 
    .[, tipas := "MxTOP"] %>% 
    .[, .(grupe, frakcija, nCG, tipas)] %>% 
    .[grupe == "uCG-_hmC-", grupe := "hmC-uCG-"] %>% 
    .[grupe == "uCG+_hmC+", grupe := "hmC+uCG+"] %>% 
    .[]

rbind(dTAB, dMxTOP) %>% 
  as.data.table() %>% 
  .[nCG > 4, ] %>% 
  ggplot(aes(grupe, frakcija, fill=tipas)) +
    geom_boxplot(outlier.shape=NA) +
    theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme_Publication() +
    ylab("hmC fraction") +
    xlab("group")

###############################################################################
###############################################################################
############################ Sup. Figure 5
###############################################################################
###############################################################################

# D
d <- readRDS("input/FigS5_D.RDS") %>%  
   .[!is.na(type), ] %>% 
   .[type == "Cluster 5", type := "Cluster X"] %>% 
   .[type == "Cluster 4", type := "Cluster 5"] %>% 
   .[type == "Cluster X", type := "Cluster 4"] %>% 
   .[type != "Cluster X"] %>% 
   .[, type := droplevels(type)]

ggplot(d, aes(stage, cts, color="grey")) +
         geom_line(aes(group = gene_id), alpha=0.05, colour="grey") +
         geom_line(stat = "summary", fun = "mean", colour = "darkgrey", alpha=0.5, linewidth = 1, 
                  aes(group = 1)) +
         geom_line(aes(stage, cts,color=gene, group = gene_id), linewidth=0.75, alpha=1, data=d[gene != "",]) +          
         ylab("Expression change") +
         xlab("Stages") +           
         theme_bw() +
         theme_Publication() +
         scale_x_discrete(expand = c(0.01,0)) + 
         facet_wrap(~type, ncol=4) +
         scale_color_brewer(palette="Set1", name = "Gene", direction=-1) +
         guides(colour = guide_legend(override.aes = list(alpha = 1, size=2))) +
         theme(legend.position="top", 
               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
               panel.spacing = unit(1.5, "lines"),
               legend.key.size = unit(1, 'cm')) +
         guides(colour=guide_legend(nrow=1, 
                                  byrow=FALSE, 
                                 alpha=1, 
                                 override.aes = list(linewidth = 4),  
                                 keywidth = unit(10, units="mm"), 
                                 keyheight = unit(10, units="mm") ))
     
# E 
sox2 <- "ENSMUSG00000074637" 
nanog <- "ENSMUSG00000012396"
POU5F1 <- "ENSMUSG00000024406" 
SOX11 <- "ENSMUSG00000063632"
TET1 <- "ENSMUSG00000047146"
pax6 <- "ENSMUSG00000027168"

gene_ids <- c("ENSMUSG00000074637", "ENSMUSG00000024406", "ENSMUSG00000012396", "ENSMUSG00000047146", "ENSMUSG00000027168", "ENSMUSG00000063632")
gene_symbol <- c("sox2", "POU5F1", "nanog",  "TET1", "pax6", "SOX11"  )

d <- readRDS("input/FigS5_E.RDS")
d[gene_id %in% c(sox2, nanog, POU5F1, SOX11, TET1, pax6), ] %>% 
   .[, .(pu_TT_S0_GCCG, pu_TT_S1_GCCG, pu_TT_S2_GCCG, pd_TT_S0_GCCG, pd_TT_S1_GCCG, pd_TT_S2_GCCG, pr_TT_S0_GCCG, pr_TT_S1_GCCG, pr_TT_S2_GCCG, gene_id)] %>% 
   melt() %>% 
   .[, tipas := "upstream"] %>% 
   .[grepl("pd", variable), tipas := "downstream"] %>% 
   .[grepl("pr", variable), tipas := "promoter"] %>% 
   .[, stage := str_extract(variable, "S.") ] %>% 
   .[, gene_id := factor(gene_id, 
                         levels=gene_ids,
                         labels=gene_symbol) ]%>% 
   .[tipas == "promoter",] %>%                          
   ggplot(aes(stage, value, group="1")) +
      geom_line() +
      geom_point(size=2) +
      facet_wrap(~gene_id, nrow=1) +
      theme_Publication() +
      ylab("Open fraction") +
      theme(axis.title.x = element_blank()) 

###############################################################################
###############################################################################
############################ Sup. Figure 6
###############################################################################
###############################################################################

# A
d <- readRDS("input/FigS6_A.RDS")

d %>% as.data.table() %>% 
  .[, estimate := log2(estimate)] %>% 
  .[p.value <= 0.05, ] %>% 
  ggplot(aes(samplas, estimate, fill=regionas)) +
    geom_jitter(shape=21, size=1.5, width=0.1, key_glyph = draw_key_rect) +
    theme_Publication()  + 
    coord_cartesian(ylim=c(1, 2.5)) +
    scale_fill_manual(values=cols_posterGroups) +
    theme_Publication() + 
    theme(legend.key.size = unit(1, 'cm')) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3))
     
# B
d <- readRDS("input/FigS6_B.RDS")

d[["S1"]] %>% as.data.table() %>% 
    .[, pval := p.adjust] %>% 
    .[p.adjust <= 0.1] %>% 
    .[order(p.adjust), ] %>% 
    head(10) %>%  
    .[, pval := (-1)*log10(p.adjust)] %>%   
    ggplot(aes(reorder(Description, pval), pval)) +
      geom_segment(aes(xend=Description, yend=0)) +
      geom_point(size=2, fill="black", shape=21) +
      coord_flip() +
      theme_bw() +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
      theme(legend.position="bottom") +
      theme_Publication()

d[["S2"]] %>% as.data.table() %>% 
    .[, pval := p.adjust] %>% 
    .[p.adjust <= 0.1] %>% 
    .[order(p.adjust), ] %>% 
    head(10) %>%  
    .[, pval := (-1)*log10(p.adjust)] %>%   
    ggplot(aes(reorder(Description, pval), pval)) +
      geom_segment(aes(xend=Description, yend=0)) +
      geom_point(size=2, fill="black", shape=21) +
      coord_flip() +
      theme_bw() +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
      theme(legend.position="bottom") +
      theme_Publication()

d[["S1S2"]] %>% as.data.table() %>% 
    .[, pval := p.adjust] %>% 
    .[p.adjust <= 0.1] %>% 
    .[order(p.adjust), ] %>% 
    head(10) %>%  
    .[, pval := (-1)*log10(p.adjust)] %>%   
    ggplot(aes(reorder(Description, pval), pval)) +
      geom_segment(aes(xend=Description, yend=0)) +
      geom_point(size=2, fill="black", shape=21) +
      coord_flip() +
      theme_bw() +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
      theme(legend.position="bottom") +
      theme_Publication()

# C
middle <- readRDS("input/FigS6_C1.RDS")
lower <- readRDS("input/FigS6_C2.RDS")
upper <- readRDS("input/FigS6_C3.RDS")
d <- readRDS("input/FigS6_C4.RDS")
ID2DO <- d[[1]][, 1] 
RNA2DO <- c("RNA_dideja", "RNA_fastInc", "RNA_lateInc", "RNA_mazeja", "RNA_fastDec")

midle_dt <- middle %>% as.data.table() %>% 
  .[, LFC := log2((NPC+1)/(ESC+1))] %>% 
  .[anotacija %in% ID2DO$anotacija, ] %>% 
  .[RNA %in% RNA2DO, ] %>% 
  .[, potipis := "De novo OCR on gene body"] 

lower_dt <- lower %>% as.data.table() %>% 
  .[, LFC := log2((NPC+1)/(ESC+1))] %>% 
  .[anotacija %in% ID2DO$anotacija, ] %>% 
  .[RNA %in% RNA2DO, ] %>% 
  .[, potipis := "OCR gene body"] 

upper_dt <- upper %>% as.data.table() %>% 
  .[, LFC := log2((NPC+1)/(ESC+1))] %>% 
  .[anotacija %in% ID2DO$anotacija, ] %>% 
  .[RNA %in% RNA2DO, ] %>% 
  .[, potipis := "OCR distal"] 

dt <- rbind(midle_dt, lower_dt, upper_dt) %>% 
  as.data.table() %>% 
  .[, RNA := factor(RNA, levels=RNA2DO)]

ggplot(dt, aes(x=LFC, colour=potipis)) +
      geom_density(key_glyph = "timeseries") +
      theme(plot.title = element_text(size=8), strip.text.x = element_text(size = 15)) +
      facet_wrap(~RNA, nrow=1) +
      scale_colour_manual(values=c(brewer.pal(8, "Dark2")[1], brewer.pal(8, "Dark2")[2], brewer.pal(8, "Dark2")[8])) +
      theme_Publication()  +
      guides(timeseries = guide_legend(override.aes = list(size = 1))) + 
      theme(strip.text.x = element_text(size = 15))
    
