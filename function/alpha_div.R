alpha_div <- function(ps, metadata, group_levels=NULL,
                      colorRGB, xlab_title="Group",
                      outdir="./", tabidx=1, figidx=1){
  library(agricolae)
  library(dplyr)
  library(ggplot2)
  
  print(sample_sums(ps))
  alpha <- estimate_richness(ps, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson"))
  
  metadata_alpha <- cbind(metadata, alpha[rownames(metadata),])
  write.csv(metadata_alpha, paste0(outdir, "Table/Table",tabidx,".metadata_alpha.csv"))
  
  # Initialize an empty data frame to store stats if it doesn't exist
  
  colorRGB <- c("#263859", "#3b8686", "#f06966", "#1F6ED4", "#79bd9a", "#f1ac9d", "#39BAE8")
  method <- c("Observed", "Shannon", "Simpson", "Chao1")
  sub_design <- metadata_alpha
  
  # Use user-defined group levels if provided
  if(!is.null(group_levels)) {
    sub_design$Group <- factor(sub_design$Group, levels = group_levels)
  }
  
  sub_design <- sub_design %>%
    filter(rownames(sub_design) %in% sample_names(ps))
  
  for(m in method){
    model <- aov(sub_design[[m]] ~ Group, data=sub_design)
    Tukey_HSD <- TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
    Tukey_HSD_table <- as.data.frame(Tukey_HSD$Group) 
    write.table(paste(m, "\n\t", sep=""), file=paste(outdir, "Table/Table",tabidx,".alpha_",m,".txt",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)
    suppressWarnings(write.table(Tukey_HSD_table, file=paste(outdir, "Table/Table",figidx,".alpha_",m,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))
    
    # LSD test for stat label
    out <- LSD.test(model, "Group", p.adj="fdr") # alternative fdr
    stat <- out$groups
    sub_design$stat <- stat[as.character(sub_design$Group),]$groups
    max_val <- max(sub_design[,c(m)])
    min_val <- min(sub_design[,c(m)])
    x <- sub_design[,c("Group", m)]
    y <- x %>% group_by(Group) %>% summarise_(Max=paste('max(',m,')',sep=""))
    y <- as.data.frame(y)
    rownames(y) <- y$Group
    sub_design$y <- y[as.character(sub_design$Group),]$Max + (max_val-min_val)*0.1
    
    p <- ggplot(sub_design, aes(x=Group, y=sub_design[[m]])) +
      geom_boxplot(aes(color=Group), fill="white", outlier.shape  = NA, linewidth=.3) +
      geom_jitter(aes(fill=Group,color=Group), position=position_jitter(0.17), size=.5, alpha=0.7) +
      scale_y_continuous(limits = c(0, 1.2 * max_val)) +
      labs(title="", x=xlab_title, y=m) +
      geom_text(data=sub_design, aes(x=Group, y=y, color=Group, label= stat)) +
      theme(
        legend.position="none",
        axis.title.y=element_text(hjust=.5),
        plot.title=element_text(size=11, hjust=.5)
      ) +
      theme_bw() +
      scale_fill_manual(values = colorRGB) +
      scale_color_manual(values = colorRGB) +
      theme(
        axis.text.x = element_text(size=12, angle = 45, hjust = 1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size = 15, color="black")
      )
    
    ggsave(paste(outdir, "Figure/Figure", figidx, ".", m, ".pdf", sep=""), p, width = 10, height = 9, units = "cm")
    print(p)
  }
} # function END
