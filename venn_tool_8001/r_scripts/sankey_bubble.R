# 桑基气泡图生成脚本
library(ggplot2)
library(ggalluvial)
library(patchwork)
library(dplyr)
library(tidyr)
library(jsonlite)

# 从命令行参数获取数据文件路径
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript sankey_bubble.R <input_json> <output_png>")
}

input_file <- args[1]
output_file <- args[2]

# 读取JSON数据
data <- fromJSON(input_file)

# 确保数据有效
if (length(data$terms) == 0) {
  stop("No terms in data")
}

# ========== 前端可控排序和数量 ==========
# 支持 sort_by: pvalue / adj_pvalue / combined_score / intersection_size / ratio
# 支持 sort_order: asc / desc
# 支持 top_n: 1-100
adj_pvalue <- if (!is.null(data$adj_pvalue)) as.numeric(data$adj_pvalue) else as.numeric(data$pvalue)
combined_score <- if (!is.null(data$combined_score)) as.numeric(data$combined_score) else rep(0, length(data$terms))
intersection_size <- if (!is.null(data$intersection_size)) as.numeric(data$intersection_size) else sapply(data$genes, length)
sort_by <- if (!is.null(data$sort_by)) as.character(data$sort_by)[1] else "pvalue"
sort_order <- if (!is.null(data$sort_order)) as.character(data$sort_order)[1] else "asc"
top_n <- if (!is.null(data$top_n)) suppressWarnings(as.integer(data$top_n)[1]) else 25
if (is.na(top_n) || top_n < 1) top_n <- 25
top_n <- min(100, top_n)

temp_df <- data.frame(
  idx = seq_along(data$terms),
  Term = as.character(data$terms),
  fg_term_num = sapply(data$genes, length),
  Ratio = as.numeric(data$ratio),
  PValue = as.numeric(data$pvalue),
  AdjPValue = adj_pvalue,
  CombinedScore = combined_score,
  IntersectionSize = intersection_size,
  stringsAsFactors = FALSE
)

if (sort_by == "adj_pvalue") {
  sort_vec <- temp_df$AdjPValue
} else if (sort_by == "combined_score") {
  sort_vec <- temp_df$CombinedScore
} else if (sort_by == "intersection_size") {
  sort_vec <- temp_df$IntersectionSize
} else if (sort_by == "ratio") {
  sort_vec <- temp_df$Ratio
} else {
  sort_vec <- temp_df$PValue
}

if (sort_order == "desc") {
  ordered_indices <- temp_df$idx[order(sort_vec, decreasing = TRUE, na.last = TRUE)]
} else {
  ordered_indices <- temp_df$idx[order(sort_vec, decreasing = FALSE, na.last = TRUE)]
}

n_terms <- min(top_n, length(ordered_indices))
top_indices <- ordered_indices[1:n_terms]

df <- data.frame(
  TermID = 1:n_terms,
  Term = as.character(data$terms[top_indices]),
  fg_term_num = sapply(data$genes[top_indices], length),
  Ratio = as.numeric(data$ratio[top_indices]),
  PValue = as.numeric(data$pvalue[top_indices]),
  stringsAsFactors = FALSE
)

# 添加基因ID字符串
df$IDs <- sapply(data$genes[top_indices], function(g) paste(g, collapse = ","))

cat("Selected", n_terms, "terms by", sort_by, sort_order, "\n")
cat("PValue range:", min(df$PValue), "-", max(df$PValue), "\n")
# ========== 排序修改结束 ==========

# 整理数据为长格式
dfLong <- df %>%
  separate_rows(IDs, sep = ",") %>%
  mutate(Term = factor(Term, levels = unique(Term))) %>%
  mutate(IDs = as.character(IDs))

# 准备桑基图数据
dfSankey <- to_lodes_form(dfLong %>% select(c("IDs", "Term")),
                           key = "x",
                           axes = c(1, 2)) %>%
  mutate(flowColor = rep(as.character(dfLong$Term), 2))

# 定义颜色（根据实际数量动态扩展）
n_colors <- length(unique(dfSankey$stratum))
colors <- c("#3b82f6", "#10b981", "#f59e0b", "#ef4444", "#8b5cf6",
            "#ec4899", "#06b6d4", "#84cc16", "#f97316", "#14b8a6",
            "#6366f1", "#22c55e", "#eab308", "#f43f5e", "#a855f7",
            "#0ea5e9", "#84cc16", "#f97316", "#14b8a6", "#6366f1")
if (n_colors > length(colors)) {
  colors <- rep(colors, ceiling(n_colors / length(colors)))
}

# 根据数据量动态调整图片尺寸
img_width <- max(14, min(30, n_terms * 0.4))
img_height <- max(8, min(20, n_terms * 0.25))

# 绘制桑基图
sankeyPlot <- ggplot(data = dfSankey,
       aes(x = x,
           stratum = factor(stratum, levels = unique(stratum)),
           alluvium = alluvium,
           y = 1,
           label = stratum,
           fill = stratum
           )) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_flow(aes(fill = flowColor), alpha = 0.3, width = 0, knot.pos = 0.1) +
  geom_stratum(width = 0.05, color = "white") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2,
            hjust = 1, nudge_x = -0.03) +
  guides(fill = "none", color = "none") +
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), units = "cm")
  ) +
  scale_x_discrete(expand = c(0.2, 0, 0, 0)) +
  scale_fill_manual(values = colors[1:n_colors])

# 准备气泡图数据
bubbleDf <- df %>%
  mutate(Term = factor(Term, levels = rev(Term))) %>%
  arrange(Term) %>%
  mutate(Term_num = cumsum(fg_term_num) - fg_term_num / 2)

# 绘制气泡图
dot_plot <- ggplot(bubbleDf, aes(x = Ratio, y = Term_num, color = -log10(PValue))) +
  geom_point(aes(size = fg_term_num)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, sum(bubbleDf$fg_term_num, na.rm = TRUE))) +
  scale_color_gradient(low = "#fbbf24", high = "#dc2626") +
  scale_radius(range = c(1, 4), name = "Size") +
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  ) +
  theme_bw() +
  labs(size = "Count", color = "-log10(P)", y = "", x = "Ratio") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "inches"),
    panel.grid = element_blank()
  )

# 合并两个图
final_plot <- sankeyPlot + dot_plot + plot_layout(widths = c(2, 1))

# 保存图片（动态尺寸）
ggsave(output_file, final_plot, width = img_width, height = img_height, dpi = 150)

cat("Plot saved to:", output_file, "\n")
cat("Terms:", n_terms, "\n")
