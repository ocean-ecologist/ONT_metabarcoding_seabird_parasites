##SEABIRD BLOOD PARASITE ANALYSIS 2025 - script put together by Sofia Olivier 

#Set working directory 
setwd (""<path/to/your/wd>"")

##PART 1: Map of seabird colonies to show host breeding range and distribution and relative colony size (Methods) 
install.packages("ctv")
ctv::install.views("Spatial", coreOnly = TRUE)
install.packages(c("tidyverse", "sp", "rgdal", "sf", "lwgeom", "stars", "exactextractr"))
install.packages("rosm")
install.packages("ggspatial")
install.packages("prettymapr")
install.packages("ggrepel")
install.packages("ggmap")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")

library(tidyverse)
library(sf)
library(sp)
library(rosm)
library(ggspatial)
library(ggplot2)
library(readr)
library(readxl)
library(prettymapr)
library(ggrepel)

#Import seabird distribution data (in this case, the data consists of the co-ordinates of breeding locations and corresponding breeding pair number)
sb_distr <- read_xlsx("seabirdcolonies.xlsx")

#Reshape the data frame from wide format into long format 
sb_distr_tidy <- sb_distr %>%
  pivot_longer(
    cols = c(`African penguins`, `Cape gannets`, `Cape cormorants`),
    names_to = "Species",
    values_to = "Value"
  )

#Turn the co-ordinates into a spatial object 
sb_spatial <- st_as_sf(sb_distr_tidy, coords = c("Longitude", "Latitude"), crs = 4326)

#Remove all the 0 colony values from the spatial data as these colonies are "empty"
sb_spatial <- sb_spatial[sb_spatial$Value != 0, ]

#Change the order of the species so that maps appear in this order
sb_spatial$Species <- factor(sb_spatial$Species, 
                             levels = c("African penguins", "Cape gannets", "Cape cormorants"))

##Create three maps adjacent to each other of the breeding colony distribution of each host species 
ggplot() +
  annotation_map_tile(type = "osm", progress = "none", zoom = 7) +
  geom_sf(data=sb_spatial, aes(size = Value, fill = Species),
          color = "black", shape = 21, alpha = 1) +
  scale_fill_manual(values = c("African penguins" = "#DA78B7",
                               "Cape gannets" = "#95CD90",
                               "Cape cormorants" = "#8AC9E5")) +
  scale_size_area(max_size = 30) +
  facet_wrap(~Species, ncol = 1) +
  geom_text_repel(
    data = sb_spatial,
    aes(label = Location, geometry = geometry),
    stat = "sf_coordinates",
    color = "black",
    size = 3,
    segment.color = NA) +
  theme_minimal() +
  theme(strip.text = element_blank(),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(
    size = "Relative Colony size",
    fill = "Species"
  )


#The same map with no legends 
ggplot() +
  annotation_map_tile(type = "osm", progress = "none", zoom = 7) +
  geom_sf(data=sb_spatial, aes(size = Value, fill = Species),
          color = "black", shape = 21, alpha = 1) +
  scale_fill_manual(values = c("African penguins" = "#DA78B7",
                               "Cape gannets" = "#95CD90",
                               "Cape cormorants" = "#8AC9E5")) +
  scale_size_area(max_size = 30) +
  facet_wrap(~Species, ncol = 1) +
  theme_minimal() +
  theme(legend.position = "none")

#Map of just Southern Africa - no points plotted 
ggplot() +
  annotation_map_tile(type = "osm", zoom = 5) +
  coord_sf(xlim = c(10, 34), ylim = c(-36, -22), expand = FALSE, 
           crs = 4326, default_crs = NULL) +
  theme_void()

##PART 2: Analysis of blood parasite taxa sequenced from the blood of seabird hosts (Results)
#For importing and cleaning data
library(readr)
library(readxl)
library(tidyverse)
library(dplyr)
library(phyloseq)
library(pheatmap)
library(RColorBrewer)
library(vegan)
library(ggplot2)
library(pairwiseAdonis)
library(ecole)
library(writexl)
library(indicspecies)

##Short fragment 18S Sequencing Results 
#Import Operational Taxonomic Unit (OTU) tables (for three hosts) made using VSEARCH 
sr_pen_otu <- read_tsv("SR_pen_otutable.0.97.tsv")
sr_gan_otu <- read_tsv("SR_gan_otutable.0.97.tsv")
sr_com_otu <- read_tsv("SR_com_otutable.0.97.tsv")
  
#Import taxonomy tables (for three hosts) made using a SINTAX algorithm and curated 18S rRNA Apicomplexa database 
sr_pen_tax <- read_xlsx("SR_pen_tax_class.xlsx")
sr_gan_tax <- read_xlsx("SR_gan_tax_class.xlsx")
sr_com_tax <- read_xlsx("SR_com_tax_class.xlsx")

#Pivot OTU tables into long format 
sr_pen_otu_long <- sr_pen_otu %>%
  pivot_longer(-OTU_ID, names_to = "Sample", values_to = "Abundance")

sr_com_otu_long <- sr_com_otu %>%
  pivot_longer(-OTU_ID, names_to = "Sample", values_to = "Abundance")

sr_gan_otu_long <- sr_gan_otu %>%
  pivot_longer(-OTU_ID, names_to = "Sample", values_to = "Abundance")

#Import metadata for host species
pen_meta <- read_xlsx("Penguin_metadata.xlsx")
com_meta <- read_xlsx("Cormorant_metadata.xlsx")
gan_meta <- read_xlsx("Gannet_metadata.xlsx")

#Rename OTUs by host species in OTU tables (i.e. add prefix to distinguish host species)
sr_pen_otu <- sr_pen_otu %>%
  mutate(OTU_ID = paste0("PEN_", OTU_ID))

sr_com_otu <- sr_com_otu %>%
  mutate(OTU_ID = paste0("COM_", OTU_ID))

sr_gan_otu <- sr_gan_otu %>%
  mutate(OTU_ID = paste0("GAN_", OTU_ID))

#Rename OTUs by host species in taxonomy tables 
sr_pen_tax <- sr_pen_tax %>%
  mutate(OTU_ID = paste0("PEN_", OTU_ID))

sr_com_tax <- sr_com_tax %>%
  mutate(OTU_ID = paste0("COM_", OTU_ID))

sr_gan_tax <- sr_gan_tax %>%
  mutate(OTU_ID = paste0("GAN_", OTU_ID))

#Convert OTU Tables into long format 
sr_pen_otu_long <- sr_pen_otu %>%
  pivot_longer(cols = -OTU_ID, names_to = "Sample_ID", values_to = "Abundance") %>%
  left_join(sr_pen_tax, by = "OTU_ID") %>%
  left_join(pen_meta, by = "Sample_ID") 

sr_gan_otu_long <- sr_gan_otu %>%
  pivot_longer(cols = -OTU_ID, names_to = "Sample_ID", values_to = "Abundance") %>%
  left_join(sr_gan_tax, by = "OTU_ID") %>%
  left_join(gan_meta, by = "Sample_ID") 

sr_com_otu_long <- sr_com_otu %>%
  pivot_longer(cols = -OTU_ID, names_to = "Sample_ID", values_to = "Abundance") %>%
  left_join(sr_com_tax, by = "OTU_ID") %>%
  left_join(com_meta, by = "Sample_ID")

#Filter for OTU presence in host samples 
sr_pen_otu_present <- sr_pen_otu_long %>%
  filter(Abundance > 0)

sr_com_otu_present <- sr_com_otu_long %>%
  filter(Abundance > 0)

sr_gan_otu_present <- sr_gan_otu_long %>%
  filter(Abundance > 0)

##Presence of OTUs in host species samples classified to family level with ≥0.75 SINTAX bootstrap confidence
#Penguins
sr_pen_fam_75 <- sr_pen_otu_present %>%
  filter(F_conf >= 0.75)

#Gannets 
sr_gan_fam_75 <- sr_gan_otu_present %>%
  filter(F_conf >= 0.75)

#Cormorants
sr_com_fam_75 <- sr_com_otu_present %>%
  filter(F_conf >= 0.75)

#Combine all three datasets
seabird_tidy_sr <- bind_rows(sr_pen_fam_75, sr_gan_fam_75, sr_com_fam_75)

#Family-level parasite abundances across species with same size bar graphs 
#Step 1 — Sum OTU abundances per host species and family-level parasite classification
sb_host_family_sr <- seabird_tidy_sr %>%
  group_by(Host_Species, Family) %>%
  summarise(Total_Abund = sum(Abundance), .groups = "drop")

#Step 2 — Convert to relative abundance *within each host species*
sb_host_rel_sr <- sb_host_family_sr %>%
  group_by(Host_Species) %>%
  mutate(RelAbund = Total_Abund / sum(Total_Abund)) %>%
  ungroup()

#Step 3 — Change order of appearance 
sb_host_rel_sr$Host_Species <- factor(sb_host_rel_sr$Host_Species, 
                                          levels = c("African penguin", "Cape gannet", "Cape cormorant"))

#Step 4 - Plot the relative family abundances in three adjacent stacked barplots 
ggplot(sb_host_rel_sr, aes(x = Host_Species, y = RelAbund, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Babesiidae" = "darkred",
                               "Cryptosporidiidae" = "lemonchiffon2",
                               "Eimeriidae" = "cadetblue4",
                               "Sarcocystidae" = "darkgreen",
                               "Theileriidae" = "thistle")) +
  labs(
    x = "Host Species",
    y = "Relative Abundance (%)",
    fill = "Family"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
    axis.text.y = element_text(size = 18),               
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)), 
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)), 
    legend.title = element_text(size = 20, face = "bold"), 
    legend.text = element_text(size = 18),                 
    panel.grid.major.x = element_blank()
  )

##Presence of OTUs classified to genus level with ≥0.75 SINTAX bootstrap confidence
#Penguins
sr_pen_gen_75 <- sr_pen_otu_present %>%
  filter(G_conf >= 0.75)

#Gannets 
sr_gan_gen_75 <- sr_gan_otu_present %>%
  filter(G_conf >= 0.75)

#Cormorants
sr_com_gen_75 <- sr_com_otu_present %>%
  filter(G_conf >= 0.75)

#Combine all three datasets
seabird_tidy_gen_sr <- bind_rows(sr_pen_gen_75, sr_gan_gen_75, sr_com_gen_75)

#Family-level parasite abundances across species with same size bar graphs 
#Step 1 — Sum OTU abundances per host species and family-level parasite classification
sb_host_family_gen_sr <- seabird_tidy_gen_sr %>%
  group_by(Host_Species, Family) %>%
  summarise(Total_Abund = sum(Abundance), .groups = "drop")

#Step 2 — Convert to relative abundance *within each host species*
sb_host_rel_gen_sr <- sb_host_family_gen_sr %>%
  group_by(Host_Species) %>%
  mutate(RelAbund = Total_Abund / sum(Total_Abund)) %>%
  ungroup()

#Step 3 — Change order of appearance 
sb_host_rel_gen_sr$Host_Species <- factor(sb_host_rel_gen_sr$Host_Species, 
                                      levels = c("African penguin", "Cape gannet", "Cape cormorant"))

#Step 4 - Plot the relative family abundances in three adjacent stacked barplots
ggplot(sb_host_rel_gen_sr, aes(x = Host_Species, y = RelAbund, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Babesiidae" = "darkred",
                               "Cryptosporidiidae" = "lemonchiffon2",
                               "Eimeriidae" = "cadetblue4",
                               "Sarcocystidae" = "darkgreen",
                               "Theileriidae" = "thistle")) +
  labs(
    x = "Host Species",
    y = "Relative Abundance (%)",
    fill = "Family"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
    axis.text.y = element_text(size = 18),               
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)), 
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)), 
    legend.title = element_text(size = 20, face = "bold"), 
    legend.text = element_text(size = 18),                 
    panel.grid.major.x = element_blank()
  )

##Long fragment 18S Sequencing Results 
#Import Operational Taxonomic Unit (OTU) tables (for three hosts) made using VSEARCH 
lr_pen_otu <- read_tsv("LR_pen_otutable.0.97.tsv")
lr_gan_otu <- read_tsv("LR_gan_otutable.0.97.tsv")
lr_com_otu <- read_tsv("LR_com_otutable.0.97.tsv") 

#Import taxonomy tables (for three hosts) made using a SINTAX algorithm and curated 18S rRNA Apicomplexa database 
lr_pen_tax <- read_xlsx("LR_pen_tax_class.xlsx")
lr_gan_tax <- read_xlsx("LR_gan_tax_class.xlsx")
lr_com_tax <- read_xlsx("LR_com_tax_class.xlsx")

#Pivot OTU table into long format 
lr_pen_otu_long <- lr_pen_otu %>%
  pivot_longer(-OTU_ID, names_to = "Sample", values_to = "Abundance")

lr_com_otu_long <- lr_com_otu %>%
  pivot_longer(-OTU_ID, names_to = "Sample", values_to = "Abundance")

lr_gan_otu_long <- lr_gan_otu %>%
  pivot_longer(-OTU_ID, names_to = "Sample", values_to = "Abundance")

#Rename OTUs by host species in OTU tables 
lr_pen_otu <- lr_pen_otu %>%
  mutate(OTU_ID = paste0("PEN_", OTU_ID))

lr_com_otu <- lr_com_otu %>%
  mutate(OTU_ID = paste0("COM_", OTU_ID))

lr_gan_otu <- lr_gan_otu %>%
  mutate(OTU_ID = paste0("GAN_", OTU_ID))

#Rename OTUs by host species in tax tables 
lr_pen_tax <- lr_pen_tax %>%
  mutate(OTU_ID = paste0("PEN_", OTU_ID))

lr_com_tax <- lr_com_tax %>%
  mutate(OTU_ID = paste0("COM_", OTU_ID))

lr_gan_tax <- lr_gan_tax %>%
  mutate(OTU_ID = paste0("GAN_", OTU_ID))

#Convert OTU Tables into long format 
lr_pen_otu_long <- lr_pen_otu %>%
  pivot_longer(cols = -OTU_ID, names_to = "Sample_ID", values_to = "Abundance") %>%
  left_join(lr_pen_tax, by = "OTU_ID") %>%
  left_join(pen_meta, by = "Sample_ID") 

lr_gan_otu_long <- lr_gan_otu %>%
  pivot_longer(cols = -OTU_ID, names_to = "Sample_ID", values_to = "Abundance") %>%
  left_join(lr_gan_tax, by = "OTU_ID") %>%
  left_join(gan_meta, by = "Sample_ID") 

lr_com_otu_long <- lr_com_otu %>%
  pivot_longer(cols = -OTU_ID, names_to = "Sample_ID", values_to = "Abundance") %>%
  left_join(lr_com_tax, by = "OTU_ID") %>%
  left_join(com_meta, by = "Sample_ID")

#Filter for OTU presence in host samples 
lr_pen_otu_present <- lr_pen_otu_long %>%
  filter(Abundance > 0)

lr_com_otu_present <- lr_com_otu_long %>%
  filter(Abundance > 0)

lr_gan_otu_present <- lr_gan_otu_long %>%
  filter(Abundance > 0)

##Presence of OTUs classified to family level with ≥0.75 SINTAX bootstrap confidence
#Penguins
lr_pen_fam_75 <- lr_pen_otu_present %>%
  filter(F_conf >= 0.75)

#Gannets 
lr_gan_fam_75 <- lr_gan_otu_present %>%
  filter(F_conf >= 0.75)

#Cormorants
lr_com_fam_75 <- lr_com_otu_present %>%
  filter(F_conf >= 0.75)

#Combine all three datasets
seabird_tidy_lr <- bind_rows(lr_pen_fam_75, lr_gan_fam_75, lr_com_fam_75)

#Family-level parasite abundances across species with same size bar graphs 
#Step 1 — Sum OTU abundances per host species and family-level parasite classification
sb_host_family_lr <- seabird_tidy_lr %>%
  group_by(Host_Species, Family) %>%
  summarise(Total_Abund = sum(Abundance), .groups = "drop")

#Step 2 — Convert to relative abundance *within each host species*
sb_host_rel_lr <- sb_host_family_lr %>%
  group_by(Host_Species) %>%
  mutate(RelAbund = Total_Abund / sum(Total_Abund)) %>%
  ungroup()

#Step 3 — Change order of appearance 
sb_host_rel_lr$Host_Species <- factor(sb_host_rel_lr$Host_Species, 
                                      levels = c("African penguin", "Cape gannet", "Cape cormorant"))

#Step 4 - Plot the relative family abundances in three adjacent stacked barplots
ggplot(sb_host_rel_lr, aes(x = Host_Species, y = RelAbund, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Babesiidae" = "darkred",
                               "Cryptosporidiidae" = "lemonchiffon2",
                               "Eimeriidae" = "cadetblue4",
                               "Sarcocystidae" = "darkgreen",
                               "Theileriidae" = "thistle")) +
  labs(
    x = "Host Species",
    y = "Relative Abundance (%)",
    fill = "Family"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
    axis.text.y = element_text(size = 18),               
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)), 
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)), 
    legend.title = element_text(size = 20, face = "bold"), 
    legend.text = element_text(size = 18),                 
    panel.grid.major.x = element_blank()
  )

##Presence of OTUs classified to genus level with ≥0.75 SINTAX bootstrap confidence
#Penguins
lr_pen_gen_75 <- lr_pen_otu_present %>%
  filter(G_conf >= 0.75)

#Gannets 
lr_gan_gen_75 <- lr_gan_otu_present %>%
  filter(G_conf >= 0.75)

#Cormorants
lr_com_gen_75 <- lr_com_otu_present %>%
  filter(G_conf >= 0.75)

#Combine all three datasets
seabird_tidy_gen_lr <- bind_rows(lr_pen_gen_75, lr_gan_gen_75, lr_com_gen_75)

#Family-level parasite abundances across species with same size bar graphs 
#Step 1 — Sum OTU abundances per host species and family-level parasite classification
sb_host_family_gen_lr <- seabird_tidy_gen_lr %>%
  group_by(Host_Species, Family) %>%
  summarise(Total_Abund = sum(Abundance), .groups = "drop")

#Step 2 — Convert to relative abundance *within each host species*
sb_host_rel_gen_lr <- sb_host_family_gen_lr %>%
  group_by(Host_Species) %>%
  mutate(RelAbund = Total_Abund / sum(Total_Abund)) %>%
  ungroup()

#Step 3 — Change order of appearance 
sb_host_rel_gen_lr$Host_Species <- factor(sb_host_rel_gen_lr$Host_Species, 
                                      levels = c("African penguin", "Cape gannet", "Cape cormorant"))

#Step 4 - Plot the relative family abundances in three adjacent stacked barplots
ggplot(sb_host_rel_gen_lr, aes(x = Host_Species, y = RelAbund, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Babesiidae" = "darkred",
                               "Cryptosporidiidae" = "lemonchiffon2",
                               "Eimeriidae" = "cadetblue4",
                               "Sarcocystidae" = "darkgreen",
                               "Theileriidae" = "thistle")) +
  labs(
    x = "Host Species",
    y = "Relative Abundance (%)",
    fill = "Family"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
    axis.text.y = element_text(size = 18),               
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)), 
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)), 
    legend.title = element_text(size = 20, face = "bold"), 
    legend.text = element_text(size = 18),                 
    panel.grid.major.x = element_blank()
  )


##Counts of samples with different parasite *families* at ≥0.75 SINTAX bootstrap confidence using *short* fragment reads
#Penguins
sr_fam75_pen_count <- sr_pen_fam_75 %>%
  distinct(Sample_ID) %>%
  count()

sr_fam75_pen_count

#Gannets 
sr_fam75_gan_count <- sr_gan_fam_75 %>%
  distinct(Sample_ID) %>%
  count()

sr_fam75_gan_count

#Cormorants
sr_fam75_com_count <- sr_com_fam_75 %>%
  distinct(Sample_ID) %>%
  count()

sr_fam75_com_count

##Counts of samples with different parasite *genera* at ≥0.75 SINTAX bootstrap confidence using *short* fragment reads
#Penguins
sr_gen75_pen_count <- sr_pen_gen_75 %>%
  distinct(Sample_ID) %>%
  count()

sr_gen75_pen_count

#Gannets 
sr_gen75_gan_count <- sr_gan_gen_75 %>%
  distinct(Sample_ID) %>%
  count()

sr_gen75_gan_count

#Cormorants
sr_gen75_com_count <- sr_com_gen_75 %>%
  distinct(Sample_ID) %>%
  count()

sr_gen75_com_count

##Count of individuals with genus *Babesia* at ≥0.75 SINTAX bootstrap confidence using *short* fragment reads
#Penguins
sr_pen_bab_75 <- sr_pen_otu_present %>%
  filter(G_conf >= 0.75, Genus == "Babesia")

sr_pen_bab_75_count <- sr_pen_bab_75 %>%
  distinct(Sample_ID) %>%
  count()

sr_pen_bab_75_count

#Gannets 
sr_gan_bab_75 <- sr_gan_otu_present %>%
  filter(G_conf >= 0.75, Genus == "Babesia")

sr_gan_bab_75_count <- sr_gan_bab_75 %>%
  distinct(Sample_ID) %>%
  count()

sr_gan_bab_75_count

#Cormorants
sr_com_bab_75 <- sr_com_otu_present %>%
  filter(G_conf >= 0.75, Genus == "Babesia")

sr_com_bab_75_count <- sr_com_bab_75 %>%
  distinct(Sample_ID) %>%
  count()

sr_com_bab_75_count

##Count of individuals with *non-Babesia* genera at at ≥0.75 SINTAX bootstrap confidence using *short* fragment reads
#Penguins
sr_pen_nonbab_75 <- sr_pen_otu_present %>%
  filter(G_conf >= 0.75, Genus != "Babesia")

sr_pen_nonbab_75_count <- sr_pen_nonbab_75 %>%
  distinct(Sample_ID) %>%
  count()

sr_pen_nonbab_75_count

#Gannets 
sr_gan_nonbab_75 <- sr_gan_otu_present %>%
  filter(G_conf >= 0.75, Genus != "Babesia")

sr_gan_nonbab_75_count <- sr_gan_nonbab_75 %>%
  distinct(Sample_ID) %>%
  count()

sr_gan_nonbab_75_count

#Cormorants
sr_com_nonbab_75 <- sr_com_otu_present %>%
  filter(G_conf >= 0.75, Genus != "Babesia")

sr_com_nonbab_75_count <- sr_com_nonbab_75 %>%
  distinct(Sample_ID) %>%
  count()

sr_com_nonbab_75_count

##Counts of samples with different parasite *families* at ≥0.75 SINTAX bootstrap confidence using *long* fragment reads
#Penguins
lr_fam75_pen_count <- lr_pen_fam_75 %>%
  distinct(Sample_ID) %>%
  count()

lr_fam75_pen_count

#Gannets 
lr_fam75_gan_count <- lr_gan_fam_75 %>%
  distinct(Sample_ID) %>%
  count()

lr_fam75_gan_count

#Cormorants
lr_fam75_com_count <- lr_com_fam_75 %>%
  distinct(Sample_ID) %>%
  count()

lr_fam75_com_count

##Getting counts of individuals with different parasite *genera* at ≥0.75 SINTAX bootstrap confidence using *long* fragment reads
#Penguins
lr_gen75_pen_count <- lr_pen_gen_75 %>%
  distinct(Sample_ID) %>%
  count()

lr_gen75_pen_count

#Gannets 
lr_gen75_gan_count <- lr_gan_gen_75 %>%
  distinct(Sample_ID) %>%
  count()

lr_gen75_gan_count

#Cormorants
lr_gen75_com_count <- lr_com_gen_75 %>%
  distinct(Sample_ID) %>%
  count()

lr_gen75_com_count

##Count of individuals with genus *Babesia* at ≥0.75 SINTAX bootstrap confidence using *long* fragment reads
#Penguins
lr_pen_bab_75 <- lr_pen_otu_present %>%
  filter(G_conf >= 0.75, Genus == "Babesia")

lr_pen_bab_75_count <- lr_pen_bab_75 %>%
  distinct(Sample_ID) %>%
  count()

lr_pen_bab_75_count

#Gannets 
lr_gan_bab_75 <- lr_gan_otu_present %>%
  filter(G_conf >= 0.75, Genus == "Babesia")

lr_gan_bab_75_count <- lr_gan_bab_75 %>%
  distinct(Sample_ID) %>%
  count()

lr_gan_bab_75_count

#Cormorants
lr_com_bab_75 <- lr_com_otu_present %>%
  filter(G_conf >= 0.75, Genus == "Babesia")

lr_com_bab_75_count <- lr_com_bab_75 %>%
  distinct(Sample_ID) %>%
  count()

lr_com_bab_75_count

##Count of individuals with *non-Babesia* genera at ≥0.75 SINTAX bootstrap confidence using *long* fragment reads
#Penguins
lr_pen_nonbab_75 <- lr_pen_otu_present %>%
  filter(G_conf >= 0.75, Genus != "Babesia")

lr_pen_nonbab_75_count <- lr_pen_nonbab_75 %>%
  distinct(Sample_ID) %>%
  count()

lr_pen_nonbab_75_count

#Gannets 
lr_gan_nonbab_75 <- lr_gan_otu_present %>%
  filter(G_conf >= 0.75, Genus != "Babesia")

lr_gan_nonbab_75_count <- lr_gan_nonbab_75 %>%
  distinct(Sample_ID) %>%
  count()

lr_gan_nonbab_75_count

#Cormorants
lr_com_nonbab_75 <- lr_com_otu_present %>%
  filter(G_conf >= 0.75, Genus != "Babesia")

lr_com_nonbab_75_count <- lr_com_nonbab_75 %>%
  distinct(Sample_ID) %>%
  count()

lr_com_nonbab_75_count

##Non-Babesia merged OTUs table 
seabird_nonbab_sr_OTUs <- bind_rows(sr_pen_nonbab_75, sr_gan_nonbab_75, sr_com_nonbab_75)
writexl::write_xlsx(seabird_nonbab_sr_OTUs, "SB_nonbab_sr_OTUs.xlsx")

seabird_nonbab_lr_OTUs <- bind_rows(lr_pen_nonbab_75, lr_gan_nonbab_75, lr_com_nonbab_75)
writexl::write_xlsx(seabird_nonbab_lr_OTUs, "SB_nonbab_lr_OTUs.xlsx")

##Babesia presence at ≥0.75 SINTAX bootstrap confidence (both short and long fragment reads)
sr_pen_babesia <- sr_pen_otu_present %>%
  filter(Genus == "Babesia", G_conf >= 0.75)
writexl::write_xlsx(sr_pen_babesia, "SR_pen_Babesia.xlsx")

lr_pen_babesia <- lr_pen_otu_present %>%
  filter(Genus == "Babesia", G_conf >= 0.75)
writexl::write_xlsx(lr_pen_babesia, "LR_pen_Babesia.xlsx")

sr_com_babesia <- sr_com_otu_present %>%
  filter(Genus == "Babesia", G_conf >= 0.75)
writexl::write_xlsx(sr_com_babesia, "SR_com_Babesia.xlsx")

lr_com_babesia <- lr_com_otu_present %>%
  filter(Genus == "Babesia", G_conf >= 0.75)
writexl::write_xlsx(lr_com_babesia, "LR_com_Babesia.xlsx")

sr_gan_babesia <- sr_gan_otu_present %>%
  filter(Genus == "Babesia", G_conf >= 0.75)
writexl::write_xlsx(sr_com_babesia, "SR_com_Babesia.xlsx")

lr_gan_babesia <- lr_gan_otu_present %>%
  filter(Genus == "Babesia", G_conf >= 0.75)
writexl::write_xlsx(lr_com_babesia, "LR_com_Babesia.xlsx")

##Venn diagrams comparing samples with Babesia vs with different genera at ≥0.75 SINTAX bootstrap confidence
install.packages("ggvenn")
library("ggvenn")

#Venn diagrams per host species using short fragment reads 
#Penguins
#Step 1 - Separate samples with Babesia and those without Babesia
sr_pen_babesia_samples <- sr_pen_gen_75 %>%
  filter(Genus == "Babesia") %>%
  pull(Sample_ID) %>%
  unique()

sr_pen_othertaxa_samples <- sr_pen_gen_75 %>%
  filter(Genus != "Babesia") %>%
  pull(Sample_ID) %>%
  unique()

#Step 2 - Define the two groups to be compared in the Venn diagram 
sr_pen_venn_groups <- list(
  Babesia = sr_pen_babesia_samples,
  Other_taxa = sr_pen_othertaxa_samples
)

#Step 3 - Plot Venn diagram of two groups including overlap (i.e. samples with both Babesia and other taxa)
ggvenn(
  sr_pen_venn_groups,
  fill_color = c("darkred", "#CB7CB3"),
  stroke_size = 0.8,
  set_name_size = 4,
  text_size = 14
)

#Gannets
#Step 1 - Separate samples with Babesia and those without Babesia
sr_gan_babesia_samples <- sr_gan_gen_75 %>%
  filter(Genus == "Babesia") %>%
  pull(Sample_ID) %>%
  unique()

sr_gan_othertaxa_samples <- sr_gan_gen_75 %>%
  filter(Genus != "Babesia") %>%
  pull(Sample_ID) %>%
  unique()

#Step 2 - Define the two groups to be compared in the Venn diagram 
sr_gan_venn_groups <- list(
  Babesia = sr_gan_babesia_samples,
  Other_taxa = sr_gan_othertaxa_samples
)

#Step 3 - Plot Venn diagram of two groups including overlap (i.e. samples with both Babesia and other taxa)
ggvenn(
  sr_gan_venn_groups,
  fill_color = c("darkgreen", "#95CD90"),
  stroke_size = 0.8,
  set_name_size = 4,
  text_size = 14
)

#Cormorants
#Step 1 - Separate samples with Babesia and those without Babesia
sr_com_babesia_samples <- sr_com_gen_75 %>%
  filter(Genus == "Babesia") %>%
  pull(Sample_ID) %>%
  unique()

sr_com_othertaxa_samples <- sr_com_gen_75 %>%
  filter(Genus != "Babesia") %>%
  pull(Sample_ID) %>%
  unique()

#Step 2 - Define the two groups to be compared in the Venn diagram 
sr_com_venn_groups <- list(
  Babesia = sr_com_babesia_samples,
  Other_taxa = sr_com_othertaxa_samples
)

#Step 3 - Plot Venn diagram of two groups including overlap (i.e. samples with both Babesia and other taxa)
ggvenn(
  sr_com_venn_groups,
  fill_color = c("cadetblue", "#8AC9E5"),
  stroke_size = 0.8,
  set_name_size = 4,
  text_size = 14
)

#Venn diagrams per host species using short fragment reads 
#Penguins
#Step 1 - Separate samples with Babesia and those without Babesia
lr_pen_babesia_samples <- lr_pen_gen_75 %>%
  filter(Genus == "Babesia") %>%
  pull(Sample_ID) %>%
  unique()

lr_pen_othertaxa_samples <- lr_pen_gen_75 %>%
  filter(Genus != "Babesia") %>%
  pull(Sample_ID) %>%
  unique()

#Step 2 - Define the two groups to be compared in the Venn diagram 
lr_pen_venn_groups <- list(
  Babesia = lr_pen_babesia_samples,
  Other_taxa = lr_pen_othertaxa_samples
)

#Step 3 - Plot Venn diagram of two groups including overlap (i.e. samples with both Babesia and other taxa)
ggvenn(
  lr_pen_venn_groups,
  fill_color = c("darkred", "#CB7CB3"),
  stroke_size = 0.8,
  set_name_size = 4,
  text_size = 14
)

#Gannets
#Step 1 - Separate samples with Babesia and those without Babesia
lr_gan_babesia_samples <- lr_gan_gen_75 %>%
  filter(Genus == "Babesia") %>%
  pull(Sample_ID) %>%
  unique()

lr_gan_othertaxa_samples <- lr_gan_gen_75 %>%
  filter(Genus != "Babesia") %>%
  pull(Sample_ID) %>%
  unique()

#Step 2 - Define the two groups to be compared in the Venn diagram 
lr_gan_venn_groups <- list(
  Babesia = lr_gan_babesia_samples,
  Other_taxa = lr_gan_othertaxa_samples
)

#Step 3 - Plot Venn diagram of two groups including overlap (i.e. samples with both Babesia and other taxa)
ggvenn(
  lr_gan_venn_groups,
  fill_color = c("darkgreen", "#95CD90"),
  stroke_size = 0.8,
  set_name_size = 4,
  text_size = 14
)

#Cormorants
#Step 1 - Separate samples with Babesia and those without Babesia
lr_com_babesia_samples <- lr_com_gen_75 %>%
  filter(Genus == "Babesia") %>%
  pull(Sample_ID) %>%
  unique()

lr_com_othertaxa_samples <- lr_com_gen_75 %>%
  filter(Genus != "Babesia") %>%
  pull(Sample_ID) %>%
  unique()

#Step 2 - Define the two groups to be compared in the Venn diagram 
lr_com_venn_groups <- list(
  Babesia = lr_com_babesia_samples,
  Other_taxa = lr_com_othertaxa_samples
)

#Step 3 - Plot Venn diagram of two groups including overlap (i.e. samples with both Babesia and other taxa)
ggvenn(
  lr_com_venn_groups,
  fill_color = c("cadetblue", "#8AC9E5"),
  stroke_size = 0.8,
  set_name_size = 4,
  text_size = 14
)

##END OF CODE