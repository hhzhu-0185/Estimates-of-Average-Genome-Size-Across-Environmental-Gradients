library(vegan)
library(ggplot2)
library(scales)
library(agricolae)
library(ggtext)
library(readxl)
library(data.table)
library(ggpubr)
library(tidyverse)

rm(list = ls())

setwd("D:/data_others/Saline-alkali-NE/results/")

alkaline_env <- as.data.frame(fread("DA_factor.csv"))
alkaline_metag_gs_kraken <- as.data.frame(fread("genomesize/prok_ags.txt"))
alkaline_metag_gs_eukrep_vibrant <- as.data.frame(fread("genomesize/prok_eukrep_vibrant.txt"))

alkaline_asv_gs <- readRDS("amplicon/asv/asv_genomesize.RDS")

alkaline_env$metagenome_AGS_kraken <- alkaline_metag_gs_kraken$Average_Genome_Size/1000000
alkaline_env$metagenome_AGS_eukrep_vibrant <- alkaline_metag_gs_eukrep_vibrant$Average_Genome_Size/1000000
alkaline_env <- left_join(alkaline_env, alkaline_asv_gs, by="samplename")

# 350 450
# metag_euk_vibrant
p_metagenome_AGS_eukrep_vibrant <- ggplot(alkaline_env, aes(x=pH, y=metagenome_AGS_eukrep_vibrant, color = pH, size = metagenome_AGS_eukrep_vibrant)) + 
  geom_point() +
  stat_cor(method = "spearman", size = 5) + 
  scale_size_continuous(range = c(2,6),guide = "none") +
  geom_smooth(method = "lm", color = "white", se = TRUE, aes(group = 1)) +  # 线性拟合
  scale_color_gradientn(colors = "#2e6cae", values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(title = "Metagenomes",x = "pH", y = "Average genome size (Mbp)") +
  scale_y_continuous(labels = label_number(scale = 1)) + 
  theme_bw()+
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black",face = "bold"),
        axis.text.y = element_text(size = 12, color = "black",face = "bold"),
        legend.title = element_text(size = 14, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 14, face = "bold", color = "black"),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
        legend.position = "none")
p_metagenome_AGS_eukrep_vibrant
ggsave("../figure/metagenome_AGS_eukrep_vibrant.svg", p_metagenome_AGS_eukrep_vibrant, width = 3.5, height = 4.5)

# metag_kraken
p_metagenome_AGS_kraken <- ggplot(alkaline_env, aes(x=pH, y=metagenome_AGS_kraken, color = pH, size = metagenome_AGS_kraken)) + 
  geom_point() +
  stat_cor(method = "spearman", size = 5) + 
  scale_size_continuous(range = c(2,6),guide = "none") +
  geom_smooth(method = "lm", color = "white", se = TRUE, aes(group = 1)) +  # 线性拟合
  scale_color_gradientn(colors = "#2e6cae", values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(title = "Metagenomes",x = "pH", y = "Average genome size (Mbp)") +
  scale_y_continuous(labels = label_number(scale = 1)) + 
  theme_bw()+
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black",face = "bold"),
        axis.text.y = element_text(size = 12, color = "black",face = "bold"),
        legend.title = element_text(size = 14, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 14, face = "bold", color = "black"),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
        legend.position = "none")
p_metagenome_AGS_kraken
ggsave("../figure/metagenome_AGS_kraken.svg", p_metagenome_AGS_kraken, width = 3.5, height = 4.5)


# metag_asv

p_amplicon_ASV_AGS <- ggplot(alkaline_env, aes(x=pH, y=amplicon_ASV_AGS, color = pH, size = amplicon_ASV_AGS)) + 
  geom_point() +
  stat_cor(method = "spearman", size = 5, hjust = -0.2) + 
  scale_size_continuous(range = c(2,6),guide = "none") +
  geom_smooth(method = "lm", color = "white", se = TRUE, aes(group = 1)) +  # 线性拟合
  scale_color_gradientn(colors = "#2e6cae", values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(title = "16S metabarcoding",x = "pH", y = "Average genome size (Mbp)") +
  scale_y_continuous(labels = label_number(scale = 1)) + 
  theme_bw()+
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black",face = "bold"),
        axis.text.y = element_text(size = 12, color = "black",face = "bold"),
        legend.title = element_text(size = 14, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 14, face = "bold", color = "black"),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
        legend.position = "none")
p_amplicon_ASV_AGS
ggsave("../figure/amplicon_ASV_AGS.svg", p_amplicon_ASV_AGS, width = 3.5, height = 4.5)



# unmated_proportion

p_ASV_unmatched_proportion <- ggplot(alkaline_env, aes(x=pH, y=ASV_unmatched)) + 
  geom_point(size = 2, color = "purple") +
  stat_cor(method = "spearman", size = 5) + 
  geom_smooth(method = "lm", color = "white", se = TRUE, aes(group = 1)) +  # 线性拟合
  labs(title = NULL, x = "pH", y = "Unmatched proportion (%)") +
  scale_y_continuous(labels = label_number(scale = 1)) + 
  theme_bw()+
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black",face = "bold"),
        axis.text.y = element_text(size = 12, color = "black",face = "bold"),
        legend.title = element_text(size = 14, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 14, face = "bold", color = "black"),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
        legend.position = "none")
p_ASV_unmatched_proportion
ggsave("../figure/ASV_unmatched_proportion.svg", p_ASV_unmatched_proportion, width = 4, height = 4)



saveRDS(alkaline_env,"alkaline_env.RDS")


##################
# NC 

rm(list=ls())
setwd("D:/data_others/Saline-alkali-NC/result/")

nc_env <- as.data.frame(fread("../metadata/Ec.csv"))
names_ampicon <- as.data.frame(fread("../metadata/Name_amplicon.csv"))
colnames(names_ampicon)[1] = "amplicon_Run"

names_metag <- as.data.frame(fread("../metadata/Name_metag.csv"))
colnames(names_metag)[1] = "metag_Run"


nc_env <- left_join(nc_env, names_ampicon, by=c("site" = "Library Name"))
nc_env <- left_join(nc_env, names_metag, by=c("site" = "Library Name"))

nc_metag_gs_kraken <- as.data.frame(fread("genomesize/prok_kraken_ags.txt"))[,1:2]
colnames(nc_metag_gs_kraken) <- c("samplename","metagenome_AGS_kraken")

nc_metag_gs_eukrep_vibrant <- as.data.frame(fread("genomesize/prok_eukrep_vibrant_ags.txt"))[,1:2]
colnames(nc_metag_gs_eukrep_vibrant) <- c("samplename","metagenome_AGS_eukrep_vibrant")

nc_asv_gs <- readRDS("amplicon/asv/asv_genomesize.RDS")

nc_env <- left_join(nc_env, nc_metag_gs_kraken, by = c("metag_Run" = "samplename"))
nc_env <- left_join(nc_env, nc_metag_gs_eukrep_vibrant, by = c("metag_Run" = "samplename"))
nc_env <- left_join(nc_env, nc_asv_gs, by = c("amplicon_Run" = "samplename"))

nc_env$metagenome_AGS_kraken <- nc_env$metagenome_AGS_kraken/1000000
nc_env$metagenome_AGS_eukrep_vibrant <- nc_env$metagenome_AGS_eukrep_vibrant/1000000

# 350 450
# metag_euk_vibrant
p_metagenome_AGS_eukrep_vibrant <- ggplot(nc_env, aes(x = log(`EC  (dS/m)`), y=metagenome_AGS_eukrep_vibrant, color = log(`EC  (dS/m)`), size = metagenome_AGS_eukrep_vibrant)) + 
  geom_point() +
  stat_cor(method = "spearman", size = 5) + 
  scale_size_continuous(range = c(2,6),guide = "none") +
  geom_smooth(method = "lm", color = "white", se = TRUE, aes(group = 1)) +  # 线性拟合
  scale_color_gradientn(colors = "#4a9e48", values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(title = "Metagenomes",x = "Log2(EC)", y = "Average genome size (Mbp)") +
  scale_y_continuous(labels = label_number(scale = 1)) + 
  theme_bw()+
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black",face = "bold"),
        axis.text.y = element_text(size = 12, color = "black",face = "bold"),
        legend.title = element_text(size = 14, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 14, face = "bold", color = "black"),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
        legend.position = "none")
p_metagenome_AGS_eukrep_vibrant
ggsave("../figure/metagenome_AGS_eukrep_vibrant.svg", p_metagenome_AGS_eukrep_vibrant, width = 3.5, height = 4.5)

# metag_kraken
p_metagenome_AGS_kraken <- ggplot(nc_env, aes(x=log(`EC  (dS/m)`), y=metagenome_AGS_kraken, color = log(`EC  (dS/m)`), size = metagenome_AGS_kraken)) + 
  geom_point() +
  stat_cor(method = "spearman", size = 5) + 
  scale_size_continuous(range = c(2,6),guide = "none") +
  geom_smooth(method = "lm", color = "white", se = TRUE, aes(group = 1)) +  # 线性拟合
  scale_color_gradientn(colors = "#4a9e48", values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(title = "Metagenomes",x = "Log2(EC)", y = "Average genome size (Mbp)") +
  scale_y_continuous(labels = label_number(scale = 1)) + 
  theme_bw()+
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black",face = "bold"),
        axis.text.y = element_text(size = 12, color = "black",face = "bold"),
        legend.title = element_text(size = 14, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 14, face = "bold", color = "black"),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
        legend.position = "none")
p_metagenome_AGS_kraken
ggsave("../figure/metagenome_AGS_kraken.svg", p_metagenome_AGS_kraken, width = 3.5, height = 4.5)


# metag_asv

p_amplicon_ASV_AGS <- ggplot(nc_env, aes(x=log(`EC  (dS/m)`), y=amplicon_ASV_AGS, color = log(`EC  (dS/m)`), size = amplicon_ASV_AGS)) + 
  geom_point() +
  stat_cor(method = "spearman", size = 5, hjust = -0.2) + 
  scale_size_continuous(range = c(2,6),guide = "none") +
  geom_smooth(method = "lm", color = "white", se = TRUE, aes(group = 1)) +  # 线性拟合
  scale_color_gradientn(colors = "#4a9e48", values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(title = "16S metabarcoding",x = "Log2(EC)", y = "Average genome size (Mbp)") +
  scale_y_continuous(labels = label_number(scale = 1)) + 
  theme_bw()+
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black",face = "bold"),
        axis.text.y = element_text(size = 12, color = "black",face = "bold"),
        legend.title = element_text(size = 14, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 14, face = "bold", color = "black"),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
        legend.position = "none")
p_amplicon_ASV_AGS
ggsave("../figure/amplicon_ASV_AGS.svg", p_amplicon_ASV_AGS, width = 3.5, height = 4.5)



# unmated_proportion

p_ASV_unmatched_proportion <- ggplot(nc_env, aes(x=log(`EC  (dS/m)`), y=ASV_unmatched)) + 
  geom_point(size = 2, color = "purple") +
  stat_cor(method = "spearman", size = 5) + 
  geom_smooth(method = "lm", color = "white", se = TRUE, aes(group = 1)) +  # 线性拟合
  labs(title = NULL, x = "Log2(EC)", y = "Unmatched proportion (%)") +
  scale_y_continuous(labels = label_number(scale = 1)) + 
  theme_bw()+
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black",face = "bold"),
        axis.text.y = element_text(size = 12, color = "black",face = "bold"),
        legend.title = element_text(size = 14, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 14, face = "bold", color = "black"),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
        legend.position = "none")
p_ASV_unmatched_proportion
ggsave("../figure/ASV_unmatched_proportion.svg", p_ASV_unmatched_proportion, width = 4, height = 4)



saveRDS(nc_env,"nc_env.RDS")



