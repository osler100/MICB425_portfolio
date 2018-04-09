library(tidyverse)

qiime2_taxa=data.frame(qiime2@tax_table)
qiime2_otu=data.frame(qiime2@otu_table)
qiime2_otu_t=t(qiime2_otu)
qiime2_taxa_t=t(qiime2_taxa)
qiime2_total=rbind(qiime2_otu_t, qiime2_taxa_t)
qiime2_total_t=t(qiime2_total)
qiime2_total_t2=data.frame(qiime2_total_t)
qiime2_total_tax=qiime2_total_t2%>%filter(Genus== "D_5__Candidatus Scalindua")

gp2 = subset_taxa(qiime2, Genus== "D_5__Candidatus Scalindua")
plot_bar(gp2, fill="Genus")

library("tidyverse")
library("phyloseq")
library("magrittr")
load("qiime2_phyloseq.RData")

#Setting seed for our group
set.seed(9376)

#Rarefying
q.norm = rarefy_even_depth(qiime2, sample.size=100000)

#Transforming to Abundance Perc.
q.perc = transform_sample_counts(q.norm, function(x) 100 * x/sum(x))

#estimating richness
q.alpha = estimate_richness(q.norm, measures = c("Chao1", "Shannon"))

#combining alpha diversity data with biogeochemical data
q.meta.alpha = full_join(rownames_to_column(q.alpha), rownames_to_column(data.frame(m.perc@sam_data)), by = "rowname")

# to generate Alpha-diversity and Oxygen across depth plot:
q.meta.alpha %>%  
  ggplot() +
  geom_point(aes(x=Depth_m, y=Shannon, colour= "Shannon Diversity")) +
  geom_smooth(method='auto', aes(x=as.numeric(Depth_m), y=Shannon)) +
  labs(title="Alpha-diversity across depth", y="Shannon's diversity index", x="Depth (m)") +
  geom_line(aes(x=Depth_m, y=O2_uM/15, colour="O2_uM")) +
  geom_point(aes(x=Depth_m, y=O2_uM/15, colour="O2_uM"))+
  scale_y_continuous(sec.axis = sec_axis(~.*(15), name = "O2 (uM)")) +
  scale_colour_manual(values = c("blue", "red"))+
  labs(title="Alpha-diversity and Oxygen across depth", y = "Shannon’s diversity index" , x = "Depth (m)" , colour = "Parameter") +
  theme(legend.position = c(0.8, 0.9))

#Generating Oxygen vs Alpha diversity plot
q.meta.alpha %>%
  
  ggplot() +
  geom_point(aes(x=O2_uM, y=Shannon)) +
  labs(title="Alpha-diversity across oxygen", y="Shannon's diversity index", x="Oxygen (uM)")

#Grouping depths based on oxic/anoxic conditions
q.meta.alpha %>%
  mutate(O2_group = ifelse(O2_uM == 0, "anoxic", "oxic")) %>%
  
  #plot to compare alpha diversity of oxic vs anoxic regions of the ocean 
  ggplot() +
  geom_boxplot(aes(x=O2_group, y=Shannon)) +
  labs(title= "Alpha-diversity by oxic/anoxic", y="Shannon's diversity index", x="Oxygen")

#Taxa Presence and Abundance:
q.perc %>%
  plot_bar(fill="Phylum") +
  geom_bar(aes(fill=Phylum), stat="identity") +
  labs(title="Phyla across samples")

#class across samples
q.perc %>%
  plot_bar(fill="Class") +
  geom_bar(aes(fill=Class), stat="identity") +
  labs(title="Class across samples")

#Phyla across Samples; by domain
q.perc %>%
  plot_bar() +
  geom_bar(aes(fill=Domain), stat="identity") +
  facet_wrap(~Phylum, scales="free_y")+
  labs(title="Phyla across samples")

#Classes across Samples; by Phylum
q.perc %>%
  plot_bar() +
  geom_bar(aes(fill=Phylum), stat="identity") +
  facet_wrap(~Class, scales="free_y", nrow=10)+
  labs(title="Clases across samples; by Phylum")

#Fitting data to linear model
q.norm %>%
  subset_taxa(Genus=="D_5__Candidatus Scalindua") %>%
  tax_glom(taxrank = 'Genus') %>%
  psmelt() %>%
  lm(abundance)


library(magrittr)
q.perc %>%
  subset_taxa(Genus=="D_5__Candidatus Scalindua") %>%
  psmelt() %>%
  group_by(Sample) %>%
  summarize(Abundance_sum=sum(Abundance), Depth_m=mean(Depth_m)) %>%
  
  ggplot() +
  geom_point(aes(x=Depth_m, y=Abundance_sum)) +
  geom_smooth(method='lm', aes(x=as.numeric(Depth_m), y=Abundance_sum)) +
  labs(title="Abundance Candidatus Scalindua across depth")

#Linear model: Abundance Against Oxygen
q.perc %>%
  subset_taxa(Genus=="D_5__Candidatus Scalindua") %>%
  psmelt() %>%
  group_by(Sample) %>%
  summarize(Abundance_sum=sum(Abundance), O2_uM=mean(O2_uM)) %>%
  
  ggplot() +
  geom_point(aes(x=O2_uM, y=Abundance_sum)) +
  geom_smooth(method='lm', aes(x=as.numeric(O2_uM), y=Abundance_sum)) +
  labs(title="Abundance Candidatus Scalindua across Oxygen conc.")

#Richness across all samples:

q.norm %>% 
  subset_taxa(Genus="D_5__Candidatus Scalindua")

#Richness across taxon:
q.norm %>%
  subset_taxa(Genus== "D_5__Candidatus Scalindua")

#Richness of taxon across samples
q.norm %>%
  subset_taxa(Genus=="D_5__Candidatus Scalindua") %>%
  estimate_richness(measures = c("Observed"))


#Generate list of OTUs from that genus:
q.norm %>%
  subset_taxa(Genus== "D_5__Candidatus Scalindua") %>%
  otu_table()
