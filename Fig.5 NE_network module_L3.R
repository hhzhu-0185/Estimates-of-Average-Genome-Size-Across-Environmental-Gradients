# NE

rm(list = ls())

library(ggplot2)
library(MuMIn)
library(nlme)
library(reshape2)
library(edgeR)
library(ggrepel)

setwd("D:/data_others/Saline-alkali-NE/results/")

ko_eukrep_vibrant <- as.data.frame(fread("./emapper/eukrep_vibrant/KEGG_ko.sum.raw.txt"))
ko_eukrep_vibrant <- column_to_rownames(ko_eukrep_vibrant,"KEGG_ko")
ko_eukrep_vibrant <- ko_eukrep_vibrant[-1,]
ko_eukrep_vibrant <- round(ko_eukrep_vibrant)
ko_eukrep_vibrant <- ko_eukrep_vibrant[rowSums(ko_eukrep_vibrant > 0) >= 20,]
ko_eukrep_vibrant$ko <- as.factor(rownames(ko_eukrep_vibrant))

ko_eukrep_vibrant1 <- reshape2::melt(ko_eukrep_vibrant)
module <- readRDS("ko.modu.memb_eukrep_vibrant1.RDS")
module$ko <- rownames(module)
table(module$c.membership.fun.fc..)

ko_module1 <- ko_eukrep_vibrant1[ko_eukrep_vibrant1$ko %in% module[module$c.membership.fun.fc.. == "1", "ko"], ]

env <- readRDS("alkaline_env.RDS")
env$site <- paste0(1:40, ".quant")
ko_module1 <- left_join(ko_module1, env, by=c("variable" = "site"))


kegg <- as.data.frame(fread("emapper/kegg htext.txt"))

ko_eukrep_vibrant2 <- left_join(ko_eukrep_vibrant, module, by="ko")
ko_eukrep_vibrant2$ko <- substr(ko_eukrep_vibrant2$ko, 4, 9)
ko_eukrep_vibrant2 <- left_join(ko_eukrep_vibrant2, kegg, by=c("ko" = "KO"))
ko_eukrep_vibrant2 <- ko_eukrep_vibrant2[ko_eukrep_vibrant2$c.membership.fun.fc.. %in% 1:2, ]
ko_eukrep_vibrant2 <- ko_eukrep_vibrant2[!duplicated(paste0(ko_eukrep_vibrant2$ko,ko_eukrep_vibrant2$level3, ko_eukrep_vibrant2$c.membership.fun.fc..)),]

ko_eukrep_vibrant2 <- aggregate(ko~level3 + c.membership.fun.fc.., data = ko_eukrep_vibrant2[,c(41,42,45)], FUN = function(x) length(unique(x)))



ko_eukrep_vibrant2_1 <- ko_eukrep_vibrant2[ko_eukrep_vibrant2$c.membership.fun.fc.. == "1", c(1,3)]
ko_eukrep_vibrant2_2 <- ko_eukrep_vibrant2[ko_eukrep_vibrant2$c.membership.fun.fc.. == "2", c(1,3)]

ko_eukrep_vibrant2_combined <- left_join(ko_eukrep_vibrant2_1, ko_eukrep_vibrant2_2, by="level3")
ko_eukrep_vibrant2_combined <- column_to_rownames(ko_eukrep_vibrant2_combined, "level3")
#ko_eukrep_vibrant2_combined[is.na(ko_eukrep_vibrant2_combined)] <- 0
ko_eukrep_vibrant2_combined <- na.omit(ko_eukrep_vibrant2_combined)
#ko_eukrep_vibrant2_combined <- ko_eukrep_vibrant2_combined[!(rownames(ko_eukrep_vibrant2_combined) %in% remove), ]

group <- factor(rep(c("Control", "Treat"), each = 1))


dgelist <- DGEList(counts = ko_eukrep_vibrant2_combined, group = group)

bcv = 0.01
lrt <- exactTest(dgelist, dispersion=bcv^2)

results <- as.data.frame(topTags(lrt, n = nrow(dgelist$counts)))
remove <- c(kegg[kegg$level1 %in% c("Human Diseases", "Organismal Systems"), "level3"], "Cellular community - eukaryotes", "Poorly characterized")


results <- results[!(rownames(results) %in% remove), ]
results <- results[-(grep("animal", rownames(results))),]
results <- results[-(grep("yeast", rownames(results))),]


results$level3 <- rownames(results)

results_sorted <- results %>%
  mutate(level3 = reorder(level3, logFC)) %>%
  mutate(sig = ifelse(logFC > 0.5 & FDR < 0.05 , "Down", 
                      ifelse(logFC < -0.5 & FDR < 0.05 , "Up", "Stable")))

set.seed(123)
p3 <- ggplot(results_sorted, aes(x=-logFC, y=-log10(FDR), color = sig)) + 
  geom_point(alpha=0.6, size=1) +
  xlim(-6,6)+
  scale_color_manual(values=c("blue", "darkgrey", "red"))+
  geom_vline(xintercept=c(-0.5,0.5),lty=2,col="grey",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="grey",lwd=0.5) +
  scale_y_continuous(breaks=seq(0,20,by=4)) +
  geom_text_repel(
    data = subset(results_sorted, sig != "Stable"),
    aes(label = level3, color = sig, size = -log(FDR)/2),
    fontface = "bold",
    force = 15,            # 增加斥力
    max.overlaps = Inf,    # 允许处理无限重叠（强制计算出位置）
    box.padding = 0.1,     # 增加文字与周围空间的距离
    point.padding = 0.1,   # 增加文字与数据点的距离
    segment.size = 0.2,    # 连接线细一点，减少视觉干扰
    seed = 123             # 固定随机种子，这样每次跑出来的图是一样的
  ) +
  scale_size_continuous(range = c(1.2, 4)) +
  guides(size = "none") +
  labs(x="log2(fold change)", y="-log10 (FDR)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        legend.title = element_blank(),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"))

p3
ggsave("../figure/edger_L3_volcano_eukrep_vibrant.svg", p3, height = 6, width = 6)


# kraken

rm(list = ls())

library(ggplot2)
library(MuMIn)
library(nlme)
library(reshape2)
library(edgeR)

setwd("D:/data_others/Saline-alkali-NE/results/")

ko_kraken <- as.data.frame(fread("./emapper/kraken/KEGG_ko.sum.raw.txt"))
ko_kraken <- column_to_rownames(ko_kraken,"KEGG_ko")
ko_kraken <- ko_kraken[-1,]
ko_kraken <- round(ko_kraken)
ko_kraken <- ko_kraken[rowSums(ko_kraken > 0) >= 20,]
ko_kraken$ko <- as.factor(rownames(ko_kraken))

ko_kraken1 <- reshape2::melt(ko_kraken)
module <- readRDS("ko.modu.memb_kraken1.RDS")
module$ko <- rownames(module)
table(module$c.membership.fun.fc..)

ko_module1 <- ko_kraken1[ko_kraken1$ko %in% module[module$c.membership.fun.fc.. == "1", "ko"], ]

env <- readRDS("alkaline_env.RDS")
env$site <- paste0(1:40, ".quant")
ko_module1 <- left_join(ko_module1, env, by=c("variable" = "site"))


kegg <- as.data.frame(fread("emapper/kegg htext.txt"))

ko_kraken2 <- left_join(ko_kraken, module, by="ko")
ko_kraken2$ko <- substr(ko_kraken2$ko, 4, 9)
ko_kraken2 <- left_join(ko_kraken2, kegg, by=c("ko" = "KO"))
ko_kraken2 <- ko_kraken2[ko_kraken2$c.membership.fun.fc.. %in% 1:2, ]
ko_kraken2 <- ko_kraken2[!duplicated(paste0(ko_kraken2$ko,ko_kraken2$level2, ko_kraken2$c.membership.fun.fc..)),]

ko_kraken2 <- aggregate(ko~level3 + c.membership.fun.fc.., data = ko_kraken2[,c(41,42,45)], FUN = function(x) length(unique(x)))



ko_kraken2_1 <- ko_kraken2[ko_kraken2$c.membership.fun.fc.. == "1", c(1,3)]
ko_kraken2_2 <- ko_kraken2[ko_kraken2$c.membership.fun.fc.. == "2", c(1,3)]

ko_kraken2_combined <- left_join(ko_kraken2_1, ko_kraken2_2, by="level3")
ko_kraken2_combined <- column_to_rownames(ko_kraken2_combined, "level3")
#ko_kraken2_combined[is.na(ko_kraken2_combined)] <- 0
ko_kraken2_combined <- na.omit(ko_kraken2_combined)
#ko_kraken2_combined <- ko_kraken2_combined[!(rownames(ko_kraken2_combined) %in% remove), ]

group <- factor(rep(c("Control", "Treat"), each = 1))


dgelist <- DGEList(counts = ko_kraken2_combined, group = group)

bcv = 0.01
lrt <- exactTest(dgelist, dispersion=bcv^2)

results <- as.data.frame(topTags(lrt, n = nrow(dgelist$counts)))
remove <- c(kegg[kegg$level1 %in% c("Human Diseases", "Organismal Systems"), "level3"], "Cellular community - eukaryotes", "Poorly characterized")



results <- results[!(rownames(results) %in% remove), ]
results <- results[-(grep("animal", rownames(results))),]
results <- results[-(grep("yeast", rownames(results))),]


results$level3 <- rownames(results)

results_sorted <- results %>%
  mutate(level3 = reorder(level3, logFC)) %>%
  mutate(sig = ifelse(logFC > 0.5 & FDR < 0.05 , "Up", 
                      ifelse(logFC < -0.5 & FDR < 0.05 , "Down", "Stable")))

set.seed(123)
p3 <- ggplot(results_sorted, aes(x=logFC, y=-log10(FDR), color = sig)) + 
  geom_point(alpha=0.6, size=1) +
  xlim(-6,6)+
  scale_color_manual(values=c("blue", "darkgrey", "red"))+
  geom_vline(xintercept=c(-0.5,0.5),lty=2,col="grey",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="grey",lwd=0.5) +
  scale_y_continuous(breaks=seq(0,20,by=4)) +
  geom_text_repel(
    data = subset(results_sorted, sig != "Stable"),
    aes(label = level3, color = sig, size = -log(FDR)/2),
    fontface = "bold",
    force = 20,            # 增加斥力
    max.overlaps = Inf,    # 允许处理无限重叠（强制计算出位置）
    box.padding = 0.1,     # 增加文字与周围空间的距离
    point.padding = 0.1,   # 增加文字与数据点的距离
    segment.size = 0.2,    # 连接线细一点，减少视觉干扰
    seed = 123             # 固定随机种子，这样每次跑出来的图是一样的
  ) +
  scale_size_continuous(range = c(1.5, 4)) +
  guides(size = "none") +
  labs(x="log2(fold change)", y="-log10 (FDR)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        legend.title = element_blank(),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"))

p3

ggsave("../figure/edger_L3_volcano_kraken.svg", p3, height = 6, width = 6)
