
rm(list=ls())
library(data.table)
library(ggplot2)
library(ggpubr)
library(vegan)
library(scales)
library(ggtext)
library(ggnewscale)
library(tidyverse)
library(ggpmisc)

setwd("D:/data_others/pH/")
env <- as.data.frame(fread("Fig.2.csv"))
name <- as.data.frame(fread("SraRunTable (7).csv"))
env <- left_join(env, name, by=c("sample_name" = "sample.name"))

ags_eukrep_vibrant <- as.data.frame(fread("ags_eukrep_vibrant.txt"))

ags_kraken <- as.data.frame(fread("ags_kraken.txt"))

table(ags_kraken$Sample_ID == ags_eukrep_vibrant$Sample_ID)

ags_all <- data.frame(sample_name = ags_kraken$Sample_ID, 
                      metag_AGS_kraken = ags_kraken$Average_Genome_Size/1000000, 
                      metag_AGS_eukrep_vibrant = ags_eukrep_vibrant$Average_Genome_Size/1000000,
                      Study = rep("Wang et al.",36))

ags_all <- left_join(ags_all, env[,c("sample_name","pH")], by=c("sample_name" = "sample_name"))




setwd("D:/data_others/Saline-alkali-NE/results/")
env <- readRDS("alkaline_env.RDS")

ags_all2 <- data.frame(sample_name = env$samplename, 
                      metag_AGS_kraken = env$metagenome_AGS_kraken, 
                      metag_AGS_eukrep_vibrant = env$metagenome_AGS_eukrep_vibrant,
                      Study = rep("This study",40),
                      pH = env$pH)



ags_all3 <- rbind(ags_all, ags_all2)

ags_all3$Study <- factor(ags_all3$Study, levels=c("Wang et al.", "This study"))

plotcol <- c("#3286BC","red")
p1 <- ggplot(ags_all3, aes(x=pH, y=metag_AGS_eukrep_vibrant, group = Study, colour = Study)) + 
  geom_point(aes(shape = Study), size = 5, alpha = 0.7) + 
  stat_poly_eq(formula = y ~ poly(x, 1, raw = TRUE),
               aes(label = paste(after_stat(rr.label), 
                                 gsub("P", "p", after_stat(p.value.label)),
                                 sep = "*\", \"*")),
               parse = TRUE,
               label.x = "right", label.y = "top", size = 5, hjust=1) +
  geom_smooth(method = "lm", color = "white", se = TRUE,aes(group = Study)) +  # 线性拟合
  scale_color_manual(name = "Type", values = plotcol)+
  scale_shape_manual(name = "Type", values=c(15,16)) +
  labs(x = "pH", y = "Average genome size (Mbp)") +
  scale_x_continuous(limits = c(3,11), breaks = seq(3,11,1)) +
  guides(color = guide_legend(order = 1),  # Habitat 图例的顺序为 1
         shape = "none")+ # Location 图例的顺序为 2
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
        axis.title = element_text(size = 17, face = "bold"),
        axis.text = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 17, face = "bold"))

p1
ggsave("../figure/eukrep_vibrant_combined_AGS.svg", p1, height = 5, width = 6.8)

p2 <- ggplot(ags_all3, aes(x=pH, y=metag_AGS_kraken, group = Study, colour = Study)) + 
  geom_point(aes(shape = Study), size = 5) + 
  stat_poly_eq(formula = y ~ poly(x, 1, raw = TRUE),
               aes(label = paste(after_stat(rr.label), 
                                 gsub("P", "p", after_stat(p.value.label)),
                                 sep = "*\", \"*")),
               parse = TRUE,
               label.x = "right", label.y = "bottom", size = 5, hjust=1, vjust=-1) +
  geom_smooth(method = "lm", color = "white", se = TRUE,aes(group = Study)) +  # 线性拟合
  scale_color_manual(name = "Type", values = plotcol)+
  scale_shape_manual(name = "Type", values=c(15,16)) +
  labs(x = "pH", y = "Average genome size (Mbp)") +
  scale_x_continuous(limits = c(3,11), breaks = seq(3,11,1)) +
  guides(color = guide_legend(order = 1),  # Habitat 图例的顺序为 1
         shape = "none")+ # Location 图例的顺序为 2
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
        axis.title = element_text(size = 17, face = "bold"),
        axis.text = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 17, face = "bold"))


p2
ggsave("../figure/kraken_combined_AGS.svg", p2, height = 5, width = 6.8)
