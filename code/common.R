cols_posterGroups  <- c("GC_only"="lightskyblue", 
                        "uCG-_hmC-"="plum3", 
                        "hmC+"=viridis(12)[5],
                        "uCG+"=viridis(12)[11],
                        "uCG+_hmC+"=viridis(12)[9]
                        )

cols_posterGroups2  <- c("Only GC"="lightskyblue", 
                        "No 5uCG/5hmC"="plum3", 
                        "5hmC"=viridis(12)[5],
                        "5uCG"=viridis(12)[11],
                        "5hmC&5uCG"=viridis(12)[9]
                        )      


cols_posterGroups3  <- c("CG-"="lightskyblue", 
                        "idCG-"="plum3", 
                        "5hmCG+"=viridis(12)[5],
                        "uCG+"=viridis(12)[11],
                        "idCG+"=viridis(12)[9]
                        )                                                    

cols_posterGroups4  <- c("GC_only"="lightskyblue", 
                        "uCGhmC_minus"="plum3", 
                        "hmC_plius"=viridis(12)[5],
                        "uCG_plius"=viridis(12)[11],
                        "uCGhmC_plius"=viridis(12)[9]
                        )


theme_Publication <- function(base_size=10, base_family="Arial") {

      (theme_foundation(base_size=base_size, base_family=base_family) + 
        theme(
              #plot.title = element_text(face = "bold",
              #                           size = rel(1.5), 
              #                           hjust = 0.5),
               text = element_text(family = base_family, size = base_size, color = "black"),

               panel.background = element_rect(fill = "transparent", colour = NA),
               #panel.border = element_rect(colour = NA),
               panel.border = element_blank(),

               plot.background = element_rect(fill = "transparent", colour = NA),
               
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               
               #axis.title.y = element_text(angle=90, vjust = 2),
               #axis.title.x = element_text(vjust = -0.2),
               axis.line = element_line(colour="black"),
               axis.text = element_text(size=base_size, color="black"),
               #axis.text.x =element_text(size = rel(1), angle = 90, vjust = 0.5, hjust=1),
               #axis.text.y =element_text(size = rel(1)),
               axis.ticks = element_line(size=0.3, color="black"),
               axis.title = element_text(size = rel(1.2), color="black"),
               
               strip.background=element_blank(),
               #strip.text = element_text(face="bold"),
               strip.text.x = element_text(size=base_size, color="black"),
               strip.text.y = element_text(size=base_size, color="black", angle = 0),

              legend.background = element_rect(
                    fill = alpha("white", 0),
                    color = alpha("white", 0)
                  ),
              legend.key = element_rect(color = NA, fill = NA),
              plot.margin=unit(c(2, 2, 2, 2),"mm"),
              legend.key.size = unit(1, "line"),
              legend.text = element_text(size = base_size),

              plot.tag=element_text(size=14, face="bold")

#               legend.key = element_rect(colour = NA),
#               legend.position = "bottom",
#               legend.direction = "horizontal",
#               legend.key.size= unit(0.2, "cm"),
#               legend.margin = margin(t = 0, unit = "cm"),
#               legend.title = element_text(face="italic"),
#               plot.margin=unit(c(5, 5, 5, 5),"mm")
               
          ))      
}