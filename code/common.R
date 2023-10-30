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
              text = element_text(family = base_family, size = base_size, color = "black"),

               panel.background = element_rect(fill = "transparent", colour = NA),
               panel.border = element_blank(),

               plot.background = element_rect(fill = "transparent", colour = NA),
               
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               
               axis.line = element_line(colour="black"),
               axis.text = element_text(size=base_size, color="black"),
               axis.ticks = element_line(size=0.3, color="black"),
               axis.title = element_text(size = rel(1.2), color="black"),
               
               strip.background=element_blank(),
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
              
          ))      
}

###############################################################################
# Function used to perform MFUZZ clustering 
# d - data to be clustered. Should be a data.table 
# x -  columns in d that should be used for clustering. 
# cl - number of clusters 
# nam - output file name 
###############################################################################

do_mfuzz <- function(d, x, cl, nam) {
  sel_counts <- d[, x, with=FALSE] 
  dd <- rbind(c(1, 2, 3), as.matrix(sel_counts))
  rownames(dd) <- c("time", d$gene_id)
  colnames(dd) <- c("S0", "S1", "S2")
  tmp <- tempfile()
  write.table(dd,file=tmp, sep='\t', quote = F, col.names=NA)
  data <- table2eset(file=tmp)
  data.s <- standardise(data)
  data.s <- data.s[rowSums(is.na(exprs(data.s)))==0] 
  saveRDS(data.s, paste0(outdatadir, nam, "_", cl, "_dataS.RDS"))
  m1 <- mestimate(data.s)
  pdf(paste0(outdatadir, "Dmin_", nam, ".pdf"), width = 9, height = 6)
  Dmin(data.s, m=m1, crange=seq(2,21,1), repeats=3, visu=TRUE)
  dev.off()
  set.seed(1987)
  c <- mfuzz(data.s,c=cl,m=m1)
  pdf(paste0(outdatadir, "Profiles_", nam, ".pdf"), width = 9, height = 6)
  mfuzz.plot(data.s, cl=c, mfrow=c(3, 3), time.labels=c("S0","S1", "S2"), new.window=FALSE)
  dev.off()
  acore <- acore(data.s,c,min.acore=0)
  acore_list <- do.call(rbind, lapply(seq_along(acore), function(i){ data.frame(CLUSTER=i, acore[[i]])}))
  saveRDS(acore_list, paste0(outdatadir, nam, "_mfuzz.RDS"))
}