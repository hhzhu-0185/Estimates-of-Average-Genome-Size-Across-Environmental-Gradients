# NC

rm(list = ls())

library(ggplot2)
library(MuMIn)
library(nlme)
library(reshape2)
library(edgeR)
library(ggrepel)

setwd("D:/data_others/Saline-alkali-NC/result/")

ko_eukrep_vibrant <- as.data.frame(fread("./KEGG/eukrep_vibrant/output3.KEGG_ko.sum.raw.txt"))
ko_eukrep_vibrant <- column_to_rownames(ko_eukrep_vibrant,"KEGG_ko")
ko_eukrep_vibrant <- ko_eukrep_vibrant[-1,]
ko_eukrep_vibrant <- round(ko_eukrep_vibrant)
ko_eukrep_vibrant <- ko_eukrep_vibrant[rowSums(ko_eukrep_vibrant > 0) >= 19,]
ko_eukrep_vibrant$ko <- rownames(ko_eukrep_vibrant)

ko_eukrep_vibrant1 <- reshape2::melt(ko_eukrep_vibrant)
module <- readRDS("ko.modu.memb_eukrep_vibrant1.RDS")
module$ko <- rownames(module)
table(module$c.membership.fun.fc..)

ko_module1 <- ko_eukrep_vibrant1[ko_eukrep_vibrant1$ko %in% module[module$c.membership.fun.fc.. == "1", "ko"], ]

env <- readRDS("nc_env.RDS")
env$metag_Run <- paste0(env$metag_Run, ".quant")
names(env)[2] <- "EC"
ko_module1 <- left_join(ko_module1, env, by=c("variable" = "metag_Run"))

ko_module1_mean <- aggregate(cbind(EC, value) ~ variable, data = ko_module1, FUN = mean)

p1 <- ggplot() + 
  geom_smooth(data=ko_module1, aes(x= log(EC), y= log(value+1), group = ko), col= "grey", method = "lm" , size = 0.05,  fill=NA)+
  geom_smooth(data=ko_module1_mean, aes(x= log(EC), y= log(value+1) ), col= "red", method = "lm", size = 5)+
  labs(x = "Log2(EC)",y = "Log (abundance + 1)")+ylim(-1,10)+
  theme_bw()+
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"))

ggsave("../figure/module1_eukrep_vibrant.svg", p1, height = 4, width = 4.3)
ggsave("../figure/module1_eukrep_vibrant.png", p1, height = 4, width = 4.3, dpi = 1000)
ggsave("../figure/module1_eukrep_vibrant2.png", p1, height = 3, width = 3.2, dpi = 1000)

lme0<-lme(value~log(EC),random=~1|ko,data=ko_module1, control=lmeControl(opt = "optim"))
summary(lme0)

anova(lme0)
r.squaredGLMM(lme0)



ko_module2 <- ko_eukrep_vibrant1[ko_eukrep_vibrant1$ko %in% module[module$c.membership.fun.fc.. == "3", "ko"], ]
ko_module2 <- left_join(ko_module2, env, by=c("variable" = "metag_Run"))

ko_module2_mean <- aggregate(cbind(EC, value) ~ variable, data = ko_module2, FUN = mean)

p2 <- ggplot() + 
  geom_smooth(data=ko_module2, aes(x= log(EC), y= log(value+1), group = ko), col= "grey", method = "lm" , size = 0.05,  fill=NA)+
  geom_smooth(data=ko_module2_mean, aes(x= log(EC), y= log(value+1) ), col= "blue", method = "lm", size = 5)+
  labs(x = "Log2(EC)",y = "Log (abundance + 1)")+ylim(-1,10)+
  theme_bw()+
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"))

ggsave("../figure/module2_eukrep_vibrant.svg", p2, height = 4, width = 4.3)
ggsave("../figure/module2_eukrep_vibrant.png", p2, height = 4, width = 4.3, dpi = 1000)
ggsave("../figure/module2_eukrep_vibrant2.png", p2, height = 3, width = 3.2, dpi = 1000)

lme0<-lme(value~log(EC),random=~1|ko,data=ko_module2, control=lmeControl(opt = "optim"))
summary(lme0)

anova(lme0)
r.squaredGLMM(lme0)

kegg <- as.data.frame(fread("kegg/kegg htext.txt"))

ko_eukrep_vibrant2 <- left_join(ko_eukrep_vibrant, module, by="ko")
ko_eukrep_vibrant2$ko <- substr(ko_eukrep_vibrant2$ko, 4, 9)
ko_eukrep_vibrant2 <- left_join(ko_eukrep_vibrant2, kegg, by=c("ko" = "KO"))
ko_eukrep_vibrant2 <- ko_eukrep_vibrant2[ko_eukrep_vibrant2$c.membership.fun.fc.. %in% c(1,3), ]
ko_eukrep_vibrant2 <- ko_eukrep_vibrant2[!duplicated(paste0(ko_eukrep_vibrant2$ko,ko_eukrep_vibrant2$level2, ko_eukrep_vibrant2$c.membership.fun.fc..)),]

ko_eukrep_vibrant2 <- aggregate(ko~level2 + c.membership.fun.fc.., data = ko_eukrep_vibrant2[,c(38,39,41)],  FUN = function(x) length(unique(x)))



ko_eukrep_vibrant2_1 <- ko_eukrep_vibrant2[ko_eukrep_vibrant2$c.membership.fun.fc.. == "1", c(1,3)]
ko_eukrep_vibrant2_2 <- ko_eukrep_vibrant2[ko_eukrep_vibrant2$c.membership.fun.fc.. == "3", c(1,3)]

ko_eukrep_vibrant2_combined <- left_join(ko_eukrep_vibrant2_1, ko_eukrep_vibrant2_2, by="level2")
ko_eukrep_vibrant2_combined <- column_to_rownames(ko_eukrep_vibrant2_combined, "level2")
#ko_eukrep_vibrant2_combined[is.na(ko_eukrep_vibrant2_combined)] <- 0
ko_eukrep_vibrant2_combined <- na.omit(ko_eukrep_vibrant2_combined)
#ko_eukrep_vibrant2_combined <- ko_eukrep_vibrant2_combined[!(rownames(ko_eukrep_vibrant2_combined) %in% remove), ]

group <- factor(rep(c("Control", "Treat"), each = 1))


dgelist <- DGEList(counts = ko_eukrep_vibrant2_combined, group = group)

bcv = 0.01
lrt <- exactTest(dgelist, dispersion=bcv^2)

results <- as.data.frame(topTags(lrt, n = nrow(dgelist$counts)))
remove <- c(kegg[kegg$level1 %in% c("Human Diseases", "Organismal Systems"), "level2"], "Cellular community - eukaryotes", "Poorly characterized")


results <- results[!(rownames(results) %in% remove), ]

results$level2 <- rownames(results)

results_sorted <- results %>%
  mutate(level2 = reorder(level2, logFC)) %>%
  mutate(sig = ifelse(logFC > 0.5, "Down", ifelse(logFC > 0 & FDR < 0.05, "Down", 
                                                  ifelse(logFC < - 0.5, "Up", ifelse(logFC < 0 & FDR < 0.05, "Up", 
                                                                                     "Stable")))))

p3 <- ggplot(results_sorted, aes(x=-logFC, y=-log10(FDR), color = sig)) + 
  geom_point(alpha=0.6, size=3) +
  xlim(-6,6)+
  scale_color_manual(values=c("blue", "darkgrey", "red"))+
  geom_vline(xintercept=c(-0.5,0.5),lty=2,col="grey",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="grey",lwd=0.5) +
  scale_y_continuous(breaks=seq(0,10,by=2)) +
  geom_text_repel(
    data = subset(results_sorted, sig != "Stable"),
    aes(label = level2, color = sig),
    size = 2,              # 稍微调小，给避让留空间
    fontface = "bold",
    force = 20,            # 增加斥力
    max.overlaps = Inf,    # 允许处理无限重叠（强制计算出位置）
    box.padding = 0.2,     # 增加文字与周围空间的距离
    point.padding = 0.2,   # 增加文字与数据点的距离
    segment.size = 0.2,    # 连接线细一点，减少视觉干扰
    seed = 123             # 固定随机种子，这样每次跑出来的图是一样的
  ) +
  labs(x="log2(fold change)", y="-log10 (FDR)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        legend.title = element_blank(),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"))

p3
ggsave("../figure/edger_volcano_eukrep_vibrant.svg", p3, height = 4, width = 4)


# kraken

rm(list = ls())


ko_kraken <- as.data.frame(fread("KEGG/kraken/output3.KEGG_ko.sum.raw.txt"))
ko_kraken <- column_to_rownames(ko_kraken,"KEGG_ko")
ko_kraken <- ko_kraken[-1,]
ko_kraken <- round(ko_kraken)
ko_kraken <- ko_kraken[rowSums(ko_kraken > 0) >= 19,]
ko_kraken$ko <- rownames(ko_kraken)

ko_kraken1 <- reshape2::melt(ko_kraken)
module <- readRDS("ko.modu.memb_kraken1.RDS")
module$ko <- rownames(module)
table(module$c.membership.fun.fc..)

ko_module1 <- ko_kraken1[ko_kraken1$ko %in% module[module$c.membership.fun.fc.. == "2", "ko"], ]

env <- readRDS("nc_env.RDS")
env$metag_Run <- paste0(env$metag_Run, ".quant")
names(env)[2] <- "EC"
ko_module1 <- left_join(ko_module1, env, by=c("variable" = "metag_Run"))

ko_module1_mean <- aggregate(cbind(EC, value) ~ variable, data = ko_module1, FUN = mean)

p1 <- ggplot() + 
  geom_smooth(data=ko_module1, aes(x= log(EC), y= log(value+1), group = ko), col= "grey", method = "lm" , size = 0.05,  fill=NA)+
  geom_smooth(data=ko_module1_mean, aes(x= log(EC), y= log(value+1) ), col= "red", method = "lm", size = 5)+
  labs(x = "Log2(EC)",y = "Log (abundance + 1)")+ylim(-1,10)+
  theme_bw()+
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"))

ggsave("../figure/module1_kraken.svg", p1, height = 4, width = 4.3)
ggsave("../figure/module1_kraken.png", p1, height = 4, width = 4.3, dpi = 1000)


lme0<-lme(value~log(EC),random=~1|ko,data=ko_module1, control=lmeControl(opt = "optim"))
summary(lme0)

anova(lme0)
r.squaredGLMM(lme0)



ko_module2 <- ko_kraken1[ko_kraken1$ko %in% module[module$c.membership.fun.fc.. == "1", "ko"], ]
ko_module2 <- left_join(ko_module2, env, by=c("variable" = "metag_Run"))

ko_module1_mean <- aggregate(cbind(EC, value) ~ variable, data = ko_module2, FUN = mean)

p2 <- ggplot() + 
  geom_smooth(data=ko_module2, aes(x= log(EC), y= log(value+1), group = ko), col= "grey", method = "lm" , size = 0.05,  fill=NA)+
  geom_smooth(data=ko_module1_mean, aes(x= log(EC), y= log(value+1) ), col= "blue", method = "lm", size = 5)+
  labs(x = "Log2(EC)",y = "Log (abundance + 1)")+ylim(-1,10)+
  theme_bw()+
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, unit = "cm"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"))

ggsave("../figure/module2_kraken.svg", p2, height = 4, width = 4.3)
ggsave("../figure/module2_kraken.png", p2, height = 4, width = 4.3, dpi = 1000)


lme0<-lme(value~log(EC),random=~1|ko,data=ko_module2, control=lmeControl(opt = "optim"))
summary(lme0)

anova(lme0)
r.squaredGLMM(lme0)

kegg <- as.data.frame(fread("KEGG/kegg htext.txt"))

ko_kraken2 <- left_join(ko_kraken, module, by="ko")
ko_kraken2$ko <- substr(ko_kraken2$ko, 4, 9)
ko_kraken2 <- left_join(ko_kraken2, kegg, by=c("ko" = "KO"))
ko_kraken2 <- ko_kraken2[ko_kraken2$c.membership.fun.fc.. %in% 1:2, ]
ko_kraken2 <- ko_kraken2[!duplicated(paste0(ko_kraken2$ko,ko_kraken2$level2, ko_kraken2$c.membership.fun.fc..)),]

ko_kraken2 <- aggregate(ko~level2 + c.membership.fun.fc.., data = ko_kraken2[,c(38,39,41)], FUN = function(x) length(unique(x)))



ko_kraken2_1 <- ko_kraken2[ko_kraken2$c.membership.fun.fc.. == "1", c(1,3)]
ko_kraken2_2 <- ko_kraken2[ko_kraken2$c.membership.fun.fc.. == "2", c(1,3)]

ko_kraken2_combined <- left_join(ko_kraken2_1, ko_kraken2_2, by="level2")
ko_kraken2_combined <- column_to_rownames(ko_kraken2_combined, "level2")
#ko_kraken2_combined[is.na(ko_kraken2_combined)] <- 0
ko_kraken2_combined <- na.omit(ko_kraken2_combined)
#ko_kraken2_combined <- ko_kraken2_combined[!(rownames(ko_kraken2_combined) %in% remove), ]

group <- factor(rep(c("Control", "Treat"), each = 1))


dgelist <- DGEList(counts = ko_kraken2_combined, group = group)

bcv = 0.01
lrt <- exactTest(dgelist, dispersion=bcv^2)


results <- as.data.frame(topTags(lrt, n = nrow(dgelist$counts)))
remove <- c(kegg[kegg$level1 %in% c("Human Diseases", "Organismal Systems"), "level2"], "Cellular community - eukaryotes", "Poorly characterized")


results <- results[!(rownames(results) %in% remove), ]

results$level2 <- rownames(results)

results_sorted <- results %>%
  mutate(level2 = reorder(level2, logFC)) %>%
  mutate(sig = ifelse(logFC > 0.5, "Up", ifelse(logFC > 0 & FDR < 0.05, "Up", 
                                                ifelse(logFC < - 0.5, "Down", ifelse(logFC < 0 & FDR < 0.05, "Down", 
                                                                                     "Stable")))))



p3 <- ggplot(results_sorted, aes(x=logFC, y=-log10(FDR), color = sig)) + 
  geom_point(alpha=0.6, size=3) +
  xlim(-6,6)+
  scale_color_manual(values=c("blue", "darkgrey", "red"))+
  geom_vline(xintercept=c(-0.5,0.5),lty=2,col="grey",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="grey",lwd=0.5) +
  scale_y_continuous(breaks=seq(0,10,by=2)) +
  geom_text_repel(
    data = subset(results_sorted, sig != "Stable"),
    aes(label = level2, color = sig),
    size = 2.5,              # 稍微调小，给避让留空间
    fontface = "bold",
    force = 20,            # 增加斥力
    max.overlaps = Inf,    # 允许处理无限重叠（强制计算出位置）
    box.padding = 0.3,     # 增加文字与周围空间的距离
    point.padding = 0.3,   # 增加文字与数据点的距离
    segment.size = 0.3,    # 连接线细一点，减少视觉干扰
    seed = 123             # 固定随机种子，这样每次跑出来的图是一样的
  ) +
  labs(x="log2(fold change)", y="-log10 (FDR)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        legend.title = element_blank(),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"))


p3

ggsave("../figure/edger_volcano_kraken.svg", p3, height = 4, width = 4)

