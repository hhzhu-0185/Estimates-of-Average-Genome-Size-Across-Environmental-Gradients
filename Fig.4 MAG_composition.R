library(ggplot2)
library(ggpubr)
library(tidyverse)
library(data.table)

rm(list=ls())


setwd("D:/data_others/Saline-alkali-NE/results/MAG/")

env <- readRDS("../alkaline_env.RDS")

arc <- as.data.frame(fread("gtdbtk.ar53.summary.tsv"))
bac <- as.data.frame(fread("gtdbtk.bac120.summary.tsv"))
all_taxonomy <- rbind(arc, bac) 

abundance <- as.data.frame(fread("mag_coverage.tsv"))
checkm <- as.data.frame(fread("checkm_results.txt"))

genomesize <- as.data.frame(fread("genomesize_stat.txt"))
genomesize$V1 <- gsub(".fna", "", genomesize$V1)
genomesize <- left_join(genomesize, checkm, by=c("V1" = "Bin Id"))
#genomesize$V2 <- genomesize$V2 *((100 - 0)/genomesize$Completeness)


genomesize <- genomesize %>% 
  mutate(group = ifelse(V2 < 2, "< 2 Mb",
                        ifelse(V2 >=2 & V2 < 3, "2-3 Mb",
                               ifelse(V2 >= 3 & V2 < 4, "3-4 Mb",
                                      ifelse(V2 >= 4 & V2 < 5, "4-5 Mb","> 5Mb")))))



abundance1 <- left_join(abundance, genomesize[,c("V1", "group")], by=c("Genome" = "V1"))

abundance2 <- aggregate(.~group, abundance1[2:42], FUN = "sum")

# 将 2:41 列转换为相对丰度（百分比）
abundance2[, 2:41] <- sweep(abundance2[, 2:41], 2, colSums(abundance2[, 2:41]), FUN = "/")

abundance3 <- reshape2::melt(abundance2, id="group")
abundance3$variable <- substr(abundance3$variable, 1,4)
abundance3 <- left_join(abundance3, env, by=c("variable" = "samplename"))



abundance4 <- aggregate(.~variable, abundance3[,c("variable", "pH")], FUN = "mean")
order_new <- abundance4$variable[order(abundance4$pH, decreasing = F)]
text_new <- sprintf("%.2f", abundance4$pH[order(abundance4$pH, decreasing = F)])

p1 <- ggplot(abundance3, aes(x=factor(variable, levels = order_new), y=value, fill = factor(group, levels = c("< 2 Mb", "2-3 Mb", "3-4 Mb", "4-5 Mb", "> 5Mb")))) + 
  geom_col(position = position_stack()) +
  scale_x_discrete(labels = text_new) +
  scale_fill_manual(name = "Genome size", values = c("#ff6b35", "#ffa726", "#42a5f5", "#1565c0", "#0d47a1")) +
  labs(x= "pH", y="Relative abundance") +
  theme_bw() + 
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, face = "bold", size = 7, color = "black"),
        axis.text.y = element_text(face = "bold", size = 8, color = "black"))

p1
ggsave("../../figure/MAG_abundance.svg", p1, height = 3.5, width = 5.5)


p2 <- ggplot(abundance3, aes(x = pH, y=value, color = factor(group, levels = c("< 2 Mb", "2-3 Mb", "3-4 Mb", "4-5 Mb", "> 5Mb")))) + 
  geom_point(size=1) + 
  stat_cor(method = "spearman",  size=3.5, label.x = 8.75, label.y.npc = 1) +
  scale_color_manual(name = "Genome size", values = c("#ff6b35", "#ffa726", "#42a5f5", "#1565c0", "#0d47a1")) +
  geom_smooth(method = "lm", se = F, linewidth=1.5) +
  labs(x= "pH", y="Relative abundance") +
  ylim(0,0.65) +
  theme_bw() + 
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold", size = 10, color = "black"))

p2
ggsave("../../figure/MAG_abundance_correlation.svg", p2, height = 3.5, width = 4.4)


abundance5 <- left_join(abundance, genomesize[,c("V1", "V2")], by=c("Genome" = "V1"))
abundance5[, 2:41] <- sweep(abundance5[, 2:41], 2, colSums(abundance5[, 2:41]), FUN = "/")


for (i in 2:41) {
  abundance5[,i] <- abundance5[,i] * abundance5[,"V2"]
} 


ags_MAG <- as.data.frame(colSums(abundance5[,2:41]))
names(ags_MAG) <- "ags"

ags_MAG$samplename <- substr(rownames(ags_MAG), 1,4)

ags_MAG <- left_join(ags_MAG, env, by="samplename")

p3 <- ggplot(ags_MAG, aes(x = pH, y=ags)) + 
  geom_point(size=3, alpha = 0.5) + 
  stat_cor(method = "spearman",  label.x = 8.75, label.y.npc = 1, color = "#d35400") +
  geom_smooth(method = "lm", se = F, 	linewidth=1.5, color = "#d35400", alpha = 0.5) +
  labs(x= "pH", y="Relative abundance") +
  theme_bw() + 
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold", size = 10, color = "black"))

p3
ggsave("../../figure/MAG_AGS.svg", p3, height = 3.5, width = 3.3)





###  nc





rm(list=ls())


setwd("D:/data_others/Saline-alkali-NC/result/MAG/")

env <- readRDS("../nc_env.RDS")

arc <- as.data.frame(fread("gtdbtk.ar53.summary.tsv"))
bac <- as.data.frame(fread("gtdbtk.bac120.summary.tsv"))
all_taxonomy <- rbind(arc, bac) 

abundance <- as.data.frame(fread("all_bin_relative.tsv"))
abundance <- abundance[-1,]
checkm <- as.data.frame(fread("checkm_results.txt"))

genomesize <- as.data.frame(fread("genomesize_mag.txt"))
genomesize$V1 <- gsub(".fa", "", genomesize$V1)
genomesize <- left_join(genomesize, checkm, by=c("V1" = "Bin Id"))
#genomesize$V2 <- genomesize$V2 *((100 - genomesize$Contamination)/genomesize$Completeness)


genomesize <- genomesize %>% 
  mutate(group = ifelse(V2 < 2, "< 2 Mb",
                        ifelse(V2 >=2 & V2 < 3, "2-3 Mb",
                               ifelse(V2 >= 3 & V2 < 4, "3-4 Mb",
                                      ifelse(V2 >= 4 & V2 < 5, "4-5 Mb","> 5Mb")))))



abundance1 <- left_join(abundance, genomesize[,c("V1", "group")], by=c("Genome" = "V1"))

abundance2 <- aggregate(.~group, abundance1[2:39], FUN = "sum")

# 将 2:41 列转换为相对丰度（百分比）
abundance2[, 2:38] <- sweep(abundance2[, 2:38], 2, colSums(abundance2[, 2:38]), FUN = "/")

abundance3 <- reshape2::melt(abundance2, id="group")
abundance3$variable <- substr(abundance3$variable, 1,11)
abundance3 <- left_join(abundance3, env, by=c("variable" = "metag_Run"))



abundance4 <- aggregate(.~variable, abundance3[,c("variable", "EC  (dS/m)")], FUN = "mean")
order_new <- abundance4$variable[order(abundance4$`EC  (dS/m)`, decreasing = F)]
text_new <- sprintf("%.2f", abundance4$`EC  (dS/m)`[order(abundance4$`EC  (dS/m)`, decreasing = F)])

p1 <- ggplot(abundance3, aes(x=factor(variable, levels = order_new), y=value, fill = factor(group, levels = c("< 2 Mb", "2-3 Mb", "3-4 Mb", "4-5 Mb", "> 5Mb")))) + 
  geom_col(position = position_stack()) +
  scale_x_discrete(labels = text_new) +
  scale_fill_manual(name = "Genome size", values = c("#ff6b35", "#ffa726", "#42a5f5", "#1565c0", "#0d47a1")) +
  labs(x= "EC (dS/m)", y="Relative abundance") +
  theme_bw() + 
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, face = "bold", size = 7, color = "black"),
        axis.text.y = element_text(face = "bold", size = 8, color = "black"))

p1
ggsave("../../figure/MAG_abundance.svg", p1, height = 3.5, width = 5.5)


p2 <- ggplot(abundance3, aes(x = `EC  (dS/m)`, y=value, color = factor(group, levels = c("< 2 Mb", "2-3 Mb", "3-4 Mb", "4-5 Mb", "> 5Mb")))) + 
  geom_point(size=1) + 
  stat_cor(method = "spearman",  size=3.5, label.x = 3.5, label.y.npc = 1) +
  scale_color_manual(name = "Genome size", values = c("#ff6b35", "#ffa726", "#42a5f5", "#1565c0", "#0d47a1")) +
  geom_smooth(method = "lm", se = F, linewidth=1.5) +
  labs(x= "EC (dS/m)", y="Relative abundance") +
  ylim(0,0.7) +
  theme_bw() + 
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold", size = 10, color = "black"))

p2
ggsave("../../figure/MAG_abundance_correlation.svg", p2, height = 3.5, width = 4.4)


abundance5 <- left_join(abundance, genomesize[,c("V1", "V2")], by=c("Genome" = "V1"))
abundance5[, 2:38] <- sweep(abundance5[, 2:38], 2, colSums(abundance5[, 2:38]), FUN = "/")


for (i in 2:38) {
  abundance5[,i] <- abundance5[,i] * abundance5[,"V2"]
} 


ags_MAG <- as.data.frame(colSums(abundance5[,2:38]))
names(ags_MAG) <- "ags"

ags_MAG$samplename <- substr(rownames(ags_MAG), 1,11)

ags_MAG <- left_join(ags_MAG, env, by=c("samplename" = "metag_Run"))

p3 <- ggplot(ags_MAG, aes(x = log(`EC  (dS/m)`), y=ags)) + 
  geom_point(size=3, alpha = 0.5) + 
  stat_cor(method = "spearman",  label.x = -1, label.y.npc = 1, color = "#d35400") +
  geom_smooth(method = "lm", se = F, 	linewidth=1.5, color = "#d35400", alpha = 0.5) +
  labs(x= "Log2(EC)", y="Relative abundance") +
  theme_bw() + 
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold", size = 10, color = "black"))

p3
ggsave("../../figure/MAG_AGS.svg", p3, height = 3.5, width = 3.3)



