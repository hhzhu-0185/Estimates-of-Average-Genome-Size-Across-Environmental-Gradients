library(data.table)
library(Hmisc)
library(igraph)


rm(list = ls())
setwd("D:/data_others/Saline-alkali-NE/results/"  ) 
# ko_kraken <- as.data.frame(fread("kraken/KEGG_ko.sum.raw.txt"))
# ko_kraken <- column_to_rownames(ko_kraken,"KEGG_ko")
# ko_kraken <- ko_kraken[-1,]
# ko_kraken <- round(ko_kraken)
# ko_kraken <- ko_kraken[rowSums(ko_kraken > 0) >= 20,]
# 
# 
# ko_eukrep_vibrant <- as.data.frame(fread("eukrep_vibrant/KEGG_ko.sum.raw.txt"))
# ko_eukrep_vibrant <- column_to_rownames(ko_eukrep_vibrant,"KEGG_ko")
# ko_eukrep_vibrant <- ko_eukrep_vibrant[-1,]
# ko_eukrep_vibrant <- round(ko_eukrep_vibrant)
# ko_eukrep_vibrant <- ko_eukrep_vibrant[rowSums(ko_eukrep_vibrant > 0) >= 20,]
# 
# cor_matrix_kraken <- rcorr(as.matrix(t(ko_kraken)), type = "spearman")
# saveRDS(cor_matrix_kraken, "cor_matrix_kraken.RDS")
# 
# cor_matrix_eukrep_vibrant <- rcorr(as.matrix(t(ko_eukrep_vibrant)), type = "spearman")
# saveRDS(cor_matrix_eukrep_vibrant, "cor_matrix_eukrep_vibrant.RDS")


cor_matrix_eukrep_vibrant <- readRDS("cor_matrix_eukrep_vibrant.RDS")
cor_matrix_kraken <- readRDS("cor_matrix_kraken.RDS")

r_mat <- cor_matrix_eukrep_vibrant$r
p_mat <- cor_matrix_eukrep_vibrant$P

# 1. 提取上三角矩阵，k=1 表示跳过对角线 (即不保留自相关)
# 这一步会把对角线和下三角设为 FALSE
upper_idx <- upper.tri(r_mat, diag = FALSE)

# 2. 根据该索引提取 R 值和 P 值
r_upper <- r_mat[upper_idx]
p_upper <- p_mat[upper_idx]

# 3. 获取对应的 OTU 名称
# 获取行列名
nms <- rownames(r_mat)
# 使用 expand.grid 创建行列组合索引，并只保留上三角部分
grid <- expand.grid(Var1 = nms, Var2 = nms)
grid_upper <- grid[upper_idx, ]

# 4. 构建数据框
sig_links_upper <- data.frame(
  KO1 = grid_upper$Var1,
  KO2 = grid_upper$Var2,
  r = r_upper,
  p = p_upper
)

sig_links_upper$p <- p.adjust(sig_links_upper$p, method = "fdr")
sig_links_filtered <- subset(sig_links_upper, p < 0.05 & r > 0.6)


g <- graph_from_data_frame(sig_links_filtered, directed=FALSE)

fun.fc<-cluster_fast_greedy(g)
print(modularity(fun.fc))
modularity(fun.fc,membership(fun.fc))

ko.modu.memb <- data.frame(c(membership(fun.fc) ))
saveRDS(ko.modu.memb, "ko.modu.memb_eukrep_vibrant1.RDS")

print(sizes(fun.fc))###


fun.comps <- membership(fun.fc)
node_colors <- ifelse(fun.comps == 2, "blue",
                      ifelse(fun.comps == 1, "red", "transparent"))

V(g)$color <- node_colors


set.seed(123)
pdf("Network_eukrep_vibrant1.pdf", width = 4, height = 4)
plot(g, layout = layout_with_kk, edge.width=0.07,edge.color="grey", vertex.frame.color=NA,vertex.label=NA,edge.lty=1,
     edge.curved=T,vertex.size=1,margin=c(0, 0,0,0))
dev.off()

# set.seed(123)
# png("Network_eukrep_vibrant.png", width = 4, height = 4, units = "in", res = 1000)
# 
# plot(g, layout = layout_with_kk, edge.width = 0.07, edge.color = "grey", 
#      vertex.frame.color = NA, vertex.label = NA, edge.lty = 1,
#      edge.curved = TRUE, vertex.size = 1, margin = 0) # margin设为0或较小值
# 
# dev.off()


###################

r_mat <- cor_matrix_kraken$r
p_mat <- cor_matrix_kraken$P

# 1. 提取上三角矩阵，k=1 表示跳过对角线 (即不保留自相关)
# 这一步会把对角线和下三角设为 FALSE
upper_idx <- upper.tri(r_mat, diag = FALSE)

# 2. 根据该索引提取 R 值和 P 值
r_upper <- r_mat[upper_idx]
p_upper <- p_mat[upper_idx]

# 3. 获取对应的 OTU 名称
# 获取行列名
nms <- rownames(r_mat)
# 使用 expand.grid 创建行列组合索引，并只保留上三角部分
grid <- expand.grid(Var1 = nms, Var2 = nms)
grid_upper <- grid[upper_idx, ]

# 4. 构建数据框
sig_links_upper <- data.frame(
  KO1 = grid_upper$Var1,
  KO2 = grid_upper$Var2,
  r = r_upper,
  p = p_upper
)

sig_links_upper$p <- p.adjust(sig_links_upper$p, method = "fdr")
sig_links_filtered <- subset(sig_links_upper, p < 0.05 & r > 0.6)


g <- graph_from_data_frame(sig_links_filtered, directed=FALSE)

fun.fc<-cluster_fast_greedy(g)
print(modularity(fun.fc))
modularity(fun.fc,membership(fun.fc))

ko.modu.memb <- data.frame(c(membership(fun.fc) ))
saveRDS(ko.modu.memb, "ko.modu.memb_kraken1.RDS")

print(sizes(fun.fc))###


fun.comps <- membership(fun.fc)
node_colors <- ifelse(fun.comps == 1, "blue",
                      ifelse(fun.comps == 2, "red", "transparent"))


V(g)$color <- node_colors


set.seed(123)
pdf("Network_kraken1.pdf", width = 4, height = 4)
plot(g, layout = layout_with_kk, edge.width=0.07,edge.color="grey", vertex.frame.color=NA,vertex.label=NA,edge.lty=1,
     edge.curved=T,vertex.size=1,margin=c(0, 0,0,0))
dev.off()

# set.seed(123)
# png("Network_kraken.png", width = 4, height = 4, units = "in", res = 1000)
# 
# plot(g, layout = layout_with_kk, edge.width = 0.07, edge.color = "grey", 
#      vertex.frame.color = NA, vertex.label = NA, edge.lty = 1,
#      edge.curved = TRUE, vertex.size = 1, margin = 0) # margin设为0或较小值
# 
# dev.off()


# nc




rm(list = ls())
setwd("D:/data_others/Saline-alkali-NC/result/") 
ko_kraken <- as.data.frame(fread("kegg/kraken/output3.KEGG_ko.sum.raw.txt"))
ko_kraken <- column_to_rownames(ko_kraken,"KEGG_ko")
ko_kraken <- ko_kraken[-1,]
ko_kraken <- round(ko_kraken)
ko_kraken <- ko_kraken[rowSums(ko_kraken > 0) >= 19,]


ko_eukrep_vibrant <- as.data.frame(fread("kegg/eukrep_vibrant/output3.KEGG_ko.sum.raw.txt"))
ko_eukrep_vibrant <- column_to_rownames(ko_eukrep_vibrant,"KEGG_ko")
ko_eukrep_vibrant <- ko_eukrep_vibrant[-1,]
ko_eukrep_vibrant <- round(ko_eukrep_vibrant)
ko_eukrep_vibrant <- ko_eukrep_vibrant[rowSums(ko_eukrep_vibrant > 0) >= 19,]

cor_matrix_kraken <- rcorr(as.matrix(t(ko_kraken)), type = "spearman")
saveRDS(cor_matrix_kraken, "cor_matrix_kraken.RDS")

cor_matrix_eukrep_vibrant <- rcorr(as.matrix(t(ko_eukrep_vibrant)), type = "spearman")
saveRDS(cor_matrix_eukrep_vibrant, "cor_matrix_eukrep_vibrant.RDS")


cor_matrix_eukrep_vibrant <- readRDS("cor_matrix_eukrep_vibrant.RDS")
cor_matrix_kraken <- readRDS("cor_matrix_kraken.RDS")

r_mat <- cor_matrix_eukrep_vibrant$r
p_mat <- cor_matrix_eukrep_vibrant$P

# 1. 提取上三角矩阵，k=1 表示跳过对角线 (即不保留自相关)
# 这一步会把对角线和下三角设为 FALSE
upper_idx <- upper.tri(r_mat, diag = FALSE)

# 2. 根据该索引提取 R 值和 P 值
r_upper <- r_mat[upper_idx]
p_upper <- p_mat[upper_idx]

# 3. 获取对应的 OTU 名称
# 获取行列名
nms <- rownames(r_mat)
# 使用 expand.grid 创建行列组合索引，并只保留上三角部分
grid <- expand.grid(Var1 = nms, Var2 = nms)
grid_upper <- grid[upper_idx, ]

# 4. 构建数据框
sig_links_upper <- data.frame(
  KO1 = grid_upper$Var1,
  KO2 = grid_upper$Var2,
  r = r_upper,
  p = p_upper
)

sig_links_upper$p <- p.adjust(sig_links_upper$p, method = "fdr")
sig_links_filtered <- subset(sig_links_upper, p < 0.05 & r > 0.6)


g <- graph_from_data_frame(sig_links_filtered, directed=FALSE)

fun.fc<-cluster_fast_greedy(g)
print(modularity(fun.fc))
modularity(fun.fc,membership(fun.fc))

ko.modu.memb <- data.frame(c(membership(fun.fc) ))
saveRDS(ko.modu.memb, "ko.modu.memb_eukrep_vibrant1.RDS")

print(sizes(fun.fc))###


fun.comps <- membership(fun.fc)
node_colors <- ifelse(fun.comps == 3, "blue",
                      ifelse(fun.comps == 1, "red", "transparent"))

V(g)$color <- node_colors


set.seed(123)
pdf("Network_eukrep_vibrant1.pdf", width = 4, height = 4)
plot(g, layout = layout_with_kk, edge.width=0.07,edge.color="grey", vertex.frame.color=NA,vertex.label=NA,edge.lty=1,
     edge.curved=T,vertex.size=1,margin=c(0, 0,0,0))
dev.off()

# set.seed(123)
# png("Network_eukrep_vibrant.png", width = 4, height = 4, units = "in", res = 1000)
# 
# plot(g, layout = layout_with_kk, edge.width = 0.07, edge.color = "grey", 
#      vertex.frame.color = NA, vertex.label = NA, edge.lty = 1,
#      edge.curved = TRUE, vertex.size = 1, margin = 0) # margin设为0或较小值
# 
# dev.off()


###################

r_mat <- cor_matrix_kraken$r
p_mat <- cor_matrix_kraken$P

# 1. 提取上三角矩阵，k=1 表示跳过对角线 (即不保留自相关)
# 这一步会把对角线和下三角设为 FALSE
upper_idx <- upper.tri(r_mat, diag = FALSE)

# 2. 根据该索引提取 R 值和 P 值
r_upper <- r_mat[upper_idx]
p_upper <- p_mat[upper_idx]

# 3. 获取对应的 OTU 名称
# 获取行列名
nms <- rownames(r_mat)
# 使用 expand.grid 创建行列组合索引，并只保留上三角部分
grid <- expand.grid(Var1 = nms, Var2 = nms)
grid_upper <- grid[upper_idx, ]

# 4. 构建数据框
sig_links_upper <- data.frame(
  KO1 = grid_upper$Var1,
  KO2 = grid_upper$Var2,
  r = r_upper,
  p = p_upper
)

sig_links_upper$p <- p.adjust(sig_links_upper$p, method = "fdr")
sig_links_filtered <- subset(sig_links_upper, p < 0.05 & r > 0.6)


g <- graph_from_data_frame(sig_links_filtered, directed=FALSE)

fun.fc<-cluster_fast_greedy(g)
print(modularity(fun.fc))
modularity(fun.fc,membership(fun.fc))

ko.modu.memb <- data.frame(c(membership(fun.fc) ))
saveRDS(ko.modu.memb, "ko.modu.memb_kraken1.RDS")

print(sizes(fun.fc))###


fun.comps <- membership(fun.fc)
node_colors <- ifelse(fun.comps == 1, "blue",
                      ifelse(fun.comps == 2, "red", "transparent"))


V(g)$color <- node_colors


set.seed(123)
pdf("Network_kraken1.pdf", width = 4, height = 4)
plot(g, layout = layout_with_kk, edge.width=0.07,edge.color="grey", vertex.frame.color=NA,vertex.label=NA,edge.lty=1,
     edge.curved=T,vertex.size=1,margin=c(0, 0,0,0))
dev.off()

# set.seed(123)
# png("Network_kraken.png", width = 4, height = 4, units = "in", res = 1000)
# 
# plot(g, layout = layout_with_kk, edge.width = 0.07, edge.color = "grey", 
#      vertex.frame.color = NA, vertex.label = NA, edge.lty = 1,
#      edge.curved = TRUE, vertex.size = 1, margin = 0) # margin设为0或较小值
# 
# dev.off()



