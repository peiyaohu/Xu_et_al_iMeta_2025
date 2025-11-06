
beta_div<- function(ps, method="bray",sample_data_ps=metadata, 
                    colname="sampleid", group_levels=NULL, perm_group= NULL,
                    colorRGB, outdir="./", tabidx=1, figidx=1){
  library(vegan)
  library(dplyr)
  library(tibble)
  
  beta<- vegdist(t(otu_table(ps)),method = "bray")

  beta_df<- beta %>%  
    cmdscale(k=4) %>% as.data.frame()%>% rownames_to_column(var="ID")
  eig<- beta %>%  
    cmdscale(k=4, eig=T) %>% .$eig
  
  colnames(beta_df)<- c("ID",paste0("PC",1:4))
  beta_df<- inner_join(beta_df,sample_data_ps, by=c("ID" = colname)) %>% as_tibble()
  
  beta_df$Group<- factor(beta_df$Group, levels=group_levels)
  centroid<-beta_df %>%
    select(starts_with("PC"),Group) %>%
    group_by(Group) %>%
    summarise(PC1=mean(PC1),
              PC2=mean(PC2))
  # PCoA
  ggplot(beta_df, aes(x=PC1, y=PC2, color=Group))+
    geom_point(alpha=.7, size=2) +
    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
         title="PCoA")+
    geom_point(data=centroid,aes(fill=Group),color="black",shape=21,size=2,show.legend = F)+
    scale_fill_manual(values = colorRGB)+
    scale_color_manual(values = colorRGB)+
    theme_light()+
    theme(axis.title = element_text(size=15,colour = "black"),
          axis.text = element_text(size=12,colour = "black"),
          legend.text= element_text(size=6,colour = "black"),
          legend.title = element_text(size=7),
          legend.key.size = unit(.3, 'cm'),
          legend.margin = margin(t=1,r=1,b=1,l=1,unit="pt"))+stat_ellipse(show.legend = F)
  ggsave(paste(outdir, "Figure/Figure", figidx, ".", "PCoA", ".pdf", sep=""),  width = 10, height = 9, units = "cm")
  
  
  # Create the formula dynamically based on perm_group
  perm_formula <- as.formula(paste("beta ~", paste(perm_group, collapse = " + ")))
  
  # Perform the adonis2 test
  test <- adonis2(perm_formula, data = sample_data_ps, permutations = 1e4)
  test
  data.frame(id=c(perm_group, "Residual","Total"),
             Df=test$Df,
             SumOfSqs = test$SumOfSqs,
             R2 = test$R2,
             F= test$F,
             Pr = test$`Pr(>F)`) %>% write.csv(paste(outdir, "Table/Table",tabidx,".PERMANOVA_1e4.csv",sep=""))
  print(test)
  
}
