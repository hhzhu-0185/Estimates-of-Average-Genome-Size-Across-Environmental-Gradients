library(data.table)
library(ggplot2)
library(ggpubr)
library(vegan)
library(scales)
library(ggtext)
library(ggnewscale)
library(tidyverse)
library(ggpmisc)
library(scales)

rm(list = ls())

setwd("D:/data_others/Saline-alkali-NE/results/")

all_alpha <- readRDS("all_alpha.RDS")
alkaline_env <- readRDS("alkaline_env.RDS") 
all_alpha <- left_join(all_alpha, alkaline_env, by="samplename")


# kraken
####Fig3a####

all_alpha$kraken_Shannon_norm <- rescale(all_alpha$kraken_Shannon, to = c(0, 1))
all_alpha$amplic_shannon_norm <- rescale(all_alpha$amplicon_Shannon, to = c(0, 1))
norm <- all_alpha[,c("samplename", "pH", "kraken_Shannon_norm", "amplic_shannon_norm")]


norm1 <- norm %>% mutate(
  grp = if_else(kraken_Shannon_norm > amplic_shannon_norm, "FA", "AF")
  
)
df <- norm1 %>%
  pivot_longer(
    cols = c(kraken_Shannon_norm, amplic_shannon_norm), 
    names_to = "types",
    values_to = "values"
  )

p1 <- ggplot(df, aes(x = pH, y = values, color = types)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      "kraken_Shannon_norm" = "#E41A1C",         
      "amplic_shannon_norm" = "#377EB8" ),
    labels = c("kraken_Shannon_norm" = "H'.fun.scaled", "amplic_shannon_norm" = "H'.tax.scaled"))+
  geom_smooth(se = F, method = "lm",formula = y ~ poly(x, 2, raw = TRUE),size = 2) +  
  stat_poly_eq(formula = y ~ poly(x, 2, raw = TRUE),
               aes(label = paste(after_stat(rr.label), 
                                 gsub("P", "p", after_stat(p.value.label)),
                                 sep = "*\", \"*")),
               parse = TRUE,
               label.x = "right", label.y = "bottom", size = 5, hjust=1.5) +
  new_scale_color()+
  geom_path(aes(group = samplename, color = grp),linewidth = 0.8) +
  scale_color_manual(values = c("#377EB8","#E41A1C" )) +
  guides(color = "none")  +
  labs(x = "pH", y = "Shannon diversity (scaled)") +
  theme_bw()+
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5),
        legend.position = "top",
        legend.margin = margin(t = 0, r = 0, b = -5, l = 0))
p1

ggsave("../figure/kraken_ph_diversity.svg", p1, height = 5, width = 5)

####Fig3b####
col1 <- c("blue","#336dff","#1E90FF","#80B1DA","#ffcc33","#ffb133","lightcoral","red")
all_alpha$kraken_dis <- all_alpha$kraken_Shannon_norm - all_alpha$amplic_shannon_norm

p2 <- ggplot(all_alpha, aes(x = pH, y = kraken_dis, color = pH)) +
  geom_point(size = 5) +
  stat_cor(method = "spearman", size = 5) +
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(color = "white", se = TRUE, aes(group = 1), method = "lm") +
  scale_color_gradientn(colors = col1, values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(x = "pH", y = "H'.fun.scaled - H'.tax.scaled (ΔH')") +
  #ggtitle("KO richness") +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5))
p2
ggsave("../figure/kraken_ph_delta.svg", p2, height = 4.8, width = 5.8)


####Fig3c####
p3 <- ggplot(all_alpha, aes(x = kraken_dis, y = metagenome_AGS_kraken, color = pH)) +
  geom_point(size = 5) +
  stat_cor(method = "spearman", size = 5) +
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(color = "white", se = TRUE, aes(group = 1), method = "lm") +  # 线性拟合
  scale_color_gradientn(colors = col1, values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(x = "ΔH'", y = "Metagenomic AGS") +
  #ggtitle("KO richness") +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5))

p3
ggsave("../figure/kraken_metagAGS_delta.svg", p3, height = 4.8, width = 5.8)



####Fig3d####
p4 <- ggplot(all_alpha, aes(x = kraken_dis, y = amplicon_ASV_AGS, color = pH)) +
  geom_point(size = 5) +
  stat_cor(method = "spearman", size = 5, hjust= -0.8) +
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(color = "white", se = TRUE, aes(group = 1), method = "lm") +  # 线性拟合
  scale_color_gradientn(colors = col1, values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(x = "ΔH'", y = "16S metabacording based AGS") +
  #ggtitle("KO richness") +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5))

p4
ggsave("../figure/kraken_ampliconAGS_delta.svg", p4, height = 4.8, width = 5.8)






## eukrep_vibrant
####Fig3a####

all_alpha$eukrep_vibrant_Shannon_norm <- rescale(all_alpha$eukrep_vibrant_Shannon, to = c(0, 1))
norm <- all_alpha[,c("samplename", "pH", "eukrep_vibrant_Shannon_norm", "amplic_shannon_norm")]


norm1 <- norm %>% mutate(
  grp = if_else(eukrep_vibrant_Shannon_norm > amplic_shannon_norm, "FA", "AF")
  
)
df <- norm1 %>%
  pivot_longer(
    cols = c(eukrep_vibrant_Shannon_norm, amplic_shannon_norm), 
    names_to = "types",
    values_to = "values"
  )

p5 <- ggplot(df, aes(x = pH, y = values, color = types)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      "eukrep_vibrant_Shannon_norm" = "#E41A1C",         
      "amplic_shannon_norm" = "#377EB8" ),
    labels = c("eukrep_vibrant_Shannon_norm" = "H'.fun.scaled", "amplic_shannon_norm" = "H'.tax.scaled"))+
  geom_smooth(se = F, method = "lm",formula = y ~ poly(x, 2, raw = TRUE),size = 2) +  
  stat_poly_eq(formula = y ~ poly(x, 2, raw = TRUE),
               aes(label = paste(after_stat(rr.label), 
                                 gsub("P", "p", after_stat(p.value.label)),
                                 sep = "*\", \"*")),
               parse = TRUE,
               label.x = "right", label.y = "bottom", size = 5, hjust=1.5) +
  new_scale_color()+
  geom_path(aes(group = samplename, color = grp),linewidth = 0.8) +
  scale_color_manual(values = c("#377EB8","#E41A1C" )) +
  guides(color = "none")  +
  labs(x = "pH", y = "Shannon diversity (scaled)") +
  theme_bw()+
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5),
        legend.position = "top",
        legend.margin = margin(t = 0, r = 0, b = -5, l = 0))
p5

ggsave("../figure/eukrep_vibrant_ph_diversity.svg", p5, height = 5, width = 5)

####Fig3b####
col1 <- c("blue","#336dff","#1E90FF","#80B1DA","#ffcc33","#ffb133","lightcoral","red")
all_alpha$eukrep_vibrant_dis <- all_alpha$eukrep_vibrant_Shannon_norm - all_alpha$amplic_shannon_norm

p6 <- ggplot(all_alpha, aes(x = pH, y = eukrep_vibrant_dis, color = pH)) +
  geom_point(size = 5) +
  stat_cor(method = "spearman", size = 5) +
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(color = "white", se = TRUE, aes(group = 1), method = "lm") +
  scale_color_gradientn(colors = col1, values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(x = "pH", y = "H'.fun.scaled - H'.tax.scaled (ΔH')") +
  #ggtitle("KO richness") +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5))


p6
ggsave("../figure/eukrep_vibrant_ph_delta.svg", p6, height = 4.8, width = 5.8)

####Fig3c####
p7 <- ggplot(all_alpha, aes(x = eukrep_vibrant_dis, y = metagenome_AGS_eukrep_vibrant, color = pH)) +
  geom_point(size = 5) +
  stat_cor(method = "spearman", size = 5) +
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(color = "white", se = TRUE, aes(group = 1), method = "lm") +  # 线性拟合
  scale_color_gradientn(colors = col1, values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(x = "ΔH'", y = "Metagenomic AGS") +
  #ggtitle("KO richness") +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5))

p7
ggsave("../figure/eukrep_vibrant_metagAGS_delta.svg", p7, height = 4.8, width = 5.8)



####Fig3d####
p8 <- ggplot(all_alpha, aes(x = eukrep_vibrant_dis, y = amplicon_ASV_AGS, color = pH)) +
  geom_point(size = 5) +
  stat_cor(method = "spearman", size = 5, hjust= -0.8) +
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(color = "white", se = TRUE, aes(group = 1), method = "lm") +  # 线性拟合
  scale_color_gradientn(colors = col1, values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(x = "ΔH'", y = "16S metabacording based AGS") +
  #ggtitle("KO richness") +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5))


p8
ggsave("../figure/eukrep_vibrant_ampliconAGS_delta.svg", p8, height = 4.8, width = 5.8)


