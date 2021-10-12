# load packages
library(readxl) # read .xlsx files
library(dplyr) # data
library(tidyverse) # data
library(ggplot2) # viz
library(hrbrthemes) # viz themes
library(viridisLite) # color viz
library(viridis) # color viz
library(plotly) # interactive viz
library(heatmaply) # heatmap
library(ggpubr) # stats in ggplot


#load data
data <- read_excel("./data/data.xlsx") # all variables
scores <- read_excel("./data/scores.xlsx")
liver <- read_excel("./data/liver.xlsx")


# Figure 3B | Spleen score at 2TP

# Obtain spleen score data
spleenscore <- scores |>
  dplyr::select(id, TP, Spleen)

# Stats
# mean spleen score x TP
aggregate(spleenscore[,3], list(spleenscore$TP), mean)
# Mann-Whitney x TP
wilcox.test(spleenscore$Spleen ~ spleenscore$TP, paired = FALSE, conf.int = TRUE, conf.level = 0.95, exact = TRUE, correct = TRUE)

# Viz
spleenscore |>
  ggplot(aes(x=TP, y=Spleen, fill=TP)) +
  geom_boxplot(show.legend = FALSE) +
  expand_limits(y=0) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="E") +
  geom_jitter(show.legend = FALSE, width = 0.15, height = 0.10, shape = 21, color = "black") +
  stat_compare_means(label = "p.format") +
  ggtitle("") +
  xlab("") +
  ylab("Pathology score") +
  theme(plot.title = element_text(size = 11)) +
  theme_classic()


# Figure 3C | qPCR at 2TP

# Obtain spleen score data in log scale
pcr <- data |>
  select(id, TP, qPCR_spleen) |>
  mutate(logpcr = log10(qPCR_spleen))

# Stats
wilcox.test(pcr$qPCR_spleen ~ pcr$TP, paired = FALSE, conf.int = TRUE, conf.level = 0.95, exact = TRUE, correct = TRUE)

# Viz
pcr |>
  ggplot(aes(x=TP, y=logpcr, fill=TP)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="E") +
  geom_jitter(show.legend = FALSE, width = 0.15, height = 0.10, shape = 21, color = "black") +
  stat_compare_means(label = "p.format") +
  ggtitle("") +
  xlab("")+
  ylab("CFU equivalents (qPCR)") +
  theme(plot.title = element_text(size = 11)) +
  theme_classic()


## Figure 3D | Correlation spleen score qPCR

# Data
sc <- pcr[-1,] |>
  mutate(spleen = spleenscore$Spleen)

# Viz
sc |>
  ggscatter(x = "spleen", y = "logpcr",
            add = "reg.line",
            add.params = list(color = "red", fill = "lightgray"),
            conf.int = TRUE,
  ) +
  stat_cor(method = "spearman", label.x = 1,
           p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = "rho") +
  scale_x_continuous(name = "Pathology score") +
  scale_y_continuous(name = "CFU equivalents (qPCR)") +
  labs(title = "") +
  theme_classic()


## Figure 4 | Granuloma area at 2TP

# data
liver
liver$area <- as.numeric(liver$area)

# stats
wilcox.test(liver$area ~ liver$TP, paired = FALSE, conf.int = TRUE, conf.level = 0.95, exact = TRUE, correct = TRUE)

# viz
liver |>
  ggplot(aes(x=TP, y=area, fill=TP)) +
  geom_violin(show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6, option = "E") +
  geom_jitter(show.legend = FALSE, width = 0.15, height = 0.10, shape = 21, color = "black") +
  stat_compare_means(label = "p.format") +
  ggtitle("Granuloma area") +
  xlab("")+
  ylab("Granuloma area (um2)") +
  theme(plot.title = element_text(size = 11)) +
  scale_color_viridis_d(option = "E", guide = "none") +
  theme_classic()


## Liver score

# data
liver
liver$area <- as.numeric(liver$area)

# 1. reorder
liver$id = with(liver, reorder(id, area, median))
liver

# 2. quantile distribution | median
quantile(liver$area, p = c(.25, .5, .75))
median_total <- median(liver$area) # median total liver granuloma area
median_area <- aggregate(liver$area, list(id = liver$id), median) # median total liver granuloma area x animal

# 3. liver score
liver_score <- median_area |>
  select(id) |>
  mutate(score = ifelse(median_area$x > median_total, 1, 0))
liver_score

# 4. viz
liver |>
  ggplot(aes(x = id, y = area, fill = TP)) +
  geom_boxplot() +
  rotate() +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6, option = "E") +
  geom_jitter(show.legend = FALSE, width = 0.15, height = 0.10, shape = 21, color = "black") +
  geom_hline(yintercept = quantile(liver$area, p = c(.25, .5, .75)), linetype = 2) +
  ggtitle("Liver index") +
  xlab("")+
  ylab("Granuloma area (um2)") +
  theme(plot.title = element_text(size = 11)) +
  theme_classic()


## Patindex

# 1. select columns from scores
patindex <- scores |>
  mutate(Liver = rowSums(scores[,8:10], na.rm = TRUE)) |>
  select(id, Thoracic, cLN, Peritoneum, ingLN, Spleen, Liver)

# 2. obtain index
patindex <- patindex |>
  mutate(Index = rowSums(patindex[,2:7], na.rm = TRUE))

# 3. add TP col
patindex <- patindex |>
  mutate(TP = ifelse(scores$TP == "Early infection", 0, 1))


## Figure 5A | Pathology Summary

# data
mat <- patindex
mat <- mat %>%
  remove_rownames |>
  column_to_rownames(var = "id") |>
  mutate(TP = ifelse(mat$TP == 1, "Late infection", "Early infection"))
mat

# viz | ggheatmap() to obtain a static plot
heatmaply(mat,
          ylab = "", xlab = "",
          main = "Pathology index",
          dendogram = "row",
          show_dendrogram = F,
          label_names = c("Id", "Feature", "Value"),
          grid_color = "white",
          grid_size = 0.00000001,
          margins = c(60, 100, 40, 20),
          heatmap_layers = theme(axis.line = element_blank()),
          file = "./plots/pathology_index.html")

# open plot
browseURL("./plots/pathology_index.html")


## Figure 5B | Pathology index at 2TP

# Stats
wilcox.test(patindex$Index ~ patindex$TP, 
            alternative="two.sided", # alternative hypothesis
            paired = FALSE,
            conf.int = TRUE,
            conf.level = 0.95)

# TP as Early/Late
patindex_tp <- patindex |>
  select(Index, TP) |>
  mutate(TP = ifelse(TP == 0, "Early infection", "Late infection"))
patindex_tp

# Viz
patindex_tp |>
  ggplot(aes(x=TP, y=Index, fill=TP)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="E") +
  geom_jitter(show.legend = FALSE, width = 0.15, height = 0.10, shape = 21, color = "black") +
  stat_compare_means(label = "p.format") +
  ggtitle("Pathology index") +
  xlab("")+
  ylab("Pathology index") +
  theme(plot.title = element_text(size = 11)) +
  theme_classic()