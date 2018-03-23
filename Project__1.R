library("tidyverse")
library("phyloseq")
library("magrittr")
load("qiime2_phyloseq.RData")
set.seed(9376)
q.norm = rarefy_even_depth(qiime2, sample.size=100000)
q.perc = transform_sample_counts(q.norm, function(x) 100 * x/sum(x))

q.alpha = estimate_richness(q.norm, measures = c("Chao1", "Shannon"))

q.meta.alpha = full_join(rownames_to_column(q.alpha), rownames_to_column(data.frame(q.perc@sam_data)), by = "rowname")

q.meta.alpha %>% 
  
  ggplot() +
  geom_point(aes(x=Depth_m, y=Shannon)) +
  geom_smooth(method='auto', aes(x=as.numeric(Depth_m), y=Shannon)) +
  labs(title="Alpha-diversity across depth", y="Shannon's diversity index", x="Depth (m)")  +
  geom_line(aes(x=Depth_m, y=O2_uM/15, colour="O2_uM")) +
  geom_point(aes(x=Depth_m, y=O2_uM/15, colour="O2_uM"))+
  scale_y_continuous(sec.axis = sec_axis(~.*(15), name = "O2 (uM)")) +
  scale_colour_manual(values = c("blue", "red"))+
  labs(title="Alpha-diversity and Oxygen across depth", y = "Shannon’s diversity index" , x = "Depth (m)" , colour = "Parameter") +
  theme(legend.position = c(0.8, 0.9))

q.meta.alpha %>% 
  
  ggplot() +
  geom_point(aes(x=O2_uM, y=Shannon)) +
  labs(title="Alpha-diversity across oxygen", y="Shannon's diversity index", x="Oxygen (uM)")

q.meta.alpha %>% 
  mutate(O2_group = ifelse(O2_uM == 0, "anoxic", "oxic")) %>% 
  
  ggplot() +
  geom_boxplot(aes(x=O2_group, y=Shannon)) +
  labs(title="Alpha-diversity by oxic/anoxic", y="Shannon's diversity index", x="Oxygen")

#Taxa presence and abundance

q.perc %>% 
  
  plot_bar(fill="Domain") + 
  geom_bar(aes(fill=Domain), stat="identity") +
  labs(title="Domains across samples")

q.perc %>% 
  
  plot_bar() + 
  geom_bar(aes(fill=Phylum), stat="identity") +
  facet_wrap(~Phylum, scales="free_y")+
  labs(title="Clases across samples; by Phylum")

q.norm %>%
  subset_taxa(Genus=="Candidatus_Scalindua") %>%
  tax_glom(taxrank = 'Genus') %>%
  psmelt() %>%
  
  lm(Abundance)

library(magrittr)
q.perc %>%
  subset_taxa(Genus=="Candidatus_Scalindua") %>%
  psmelt() %>%
  group_by(Sample) %>%
  summarize(Abundance_sum=sum(Abundance), Depth_m=mean(Depth_m)) %>%
  
  ggplot() +
  geom_point(aes(x=Depth_m, y=Abundance_sum)) +
  geom_smooth(method='lm', aes(x=as.numeric(Depth_m), y=Abundance_sum)) +
  labs(title="Abundance Candidatus Scalindua across depth")

       
  #Within your taxon, what is the richness (number of OTUs/ASVs)?
  
  q.norm %>% 
    subset_taxa(Genus=="Candidatus_Scanlindua") 
  
  q.norm %>% 
    subset_taxa(Genus=="Candidatus_Scanlindua") %>%
    estimate_richness(measures = c("Observed"))
  
  #Do the abundances of OTUs/ASVs within your taxon of interest change significantly with depth and/or oxygen concentration?
  #####################
  q.norm %>% 
    psmelt() %>% 
    filter(OTU=="Otu0242") %>% 
    
    lm(Abundance ~ Depth_m, .) %>% 
    summary()
-----
subset_taxa(Genus== "Candidatus_Scalindua") %>%
    otu_table()
  
    
    
    
    
    
    
    
    
    
    
  
  
  
  
  
  
  
  
  
  p.adjust(c(0.501, 0.031, 0.005, 0.324), method="fdr")
  
  q.perc %>% 
    subset_taxa(Domain=="unknown") %>% 
    psmelt() %>% 
    
    ggplot() +
    geom_point(aes(x=Depth_m, y=Abundance)) +
    geom_smooth(method='lm', aes(x=Depth_m, y=Abundance)) +
    facet_wrap(~OTU, scales="free_y") +
    labs(title="Abundance of OTUs within unclassified domain across depth")
  
  q.perc %>% 
    subset_taxa(Domain=="unknown") %>%
    psmelt() %>% 
    
    ggplot() +
    geom_point(aes(x=Sample, y=OTU, size=Abundance, color=OTU)) + 
    scale_size_continuous(range = c(0,5)) +
    labs(title="Abundance of OTUs within unclassified domain across depth")
  
  
  q.meta.alpha %>% 
    
    ggplot() +
    geom_point(aes(x=Depth_m, y=Shannon), scale="free") +
    geom_smooth(method='auto', aes(x=as.numeric(Depth_m), y=Shannon)) +
    labs(title="Alpha-diversity across depth", y="Shannon's diversity index", x="Depth (m)") + 
    geom_point(aes(x=Depth_m, y=O2_uM), scale="free") +
    scale_y_continuous(position = "left" ) +
    geom_smooth(method='auto', aes(x=as.numeric(Depth_m), y=Shannon)) + 
    scale_y_continuous(position = "right")
  