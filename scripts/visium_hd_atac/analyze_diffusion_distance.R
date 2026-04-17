library(data.table)
library(argparse)
library(tidyverse)


parser <- ArgumentParser(description = "Analyze Spatial Diffusion Distance from 6-Column Fragments")
parser$add_argument("--fragments_file", type="character", required=TRUE, 
                    help="[Required] 6-column (chr, start, end, barcode, count, strand) fragments.tsv.gz file")
parser$add_argument("--outdir", type="character", required=TRUE, 
                    help="[Required] Output directory for saving plot and tsv file")
parser$add_argument("--bin_width_um", type="integer", default=16, 
                    help="Width (in um) for binning distances. Should match resolution.")
parser$add_argument("--max_dist_um", type="integer", default=1000, 
                    help="Maximum distance (in um) to analyze.")

Args <- parser$parse_args()

cat("--- R script starting: Spatial Diffusion Analysis ---\n")
dir.create(Args$outdir, recursive = TRUE, showWarnings = FALSE)

output_tsv <- file.path(Args$outdir, "diffusion_by_distance.tsv")
output_pdf <- file.path(Args$outdir, "diffusion_plot.pdf")

cat(paste("Loading fragments file:", Args$fragments_file, "\n"))
tryCatch({
  dt <- fread(
    cmd = paste("cat", Args$fragments_file), 
    col.names = c("chr", "start", "end", "barcode", "count", "strand")
  )
}, error = function(e) {
  stop(paste("Failed to read fragments file. Is it a 6-column .tsv.gz file?", e))
})

cat("Fragments loaded. Total reads:", dt[, sum(count)], "\n")

cat("Parsing barcodes to extract coordinates...\n")
# 示例: s_016um_00365_00044
dt[, c("s", "res_str", "x_str", "y_str") := tstrsplit(barcode, "_", fixed=TRUE)]
dt[, x := as.integer(x_str)]
dt[, y := as.integer(y_str)]
dt[, resolution := as.integer(gsub("um", "", res_str))]

# 检查解析是否成功
if (anyNA(dt$x) || anyNA(dt$y)) {
  stop("Failed to parse X/Y coordinates from barcodes. Check barcode format.")
}
# 获取此文件的分辨率 (假设所有 bins 分辨率相同)
file_resolution <- dt$resolution[1]
cat(paste("Detected file resolution:", file_resolution, "um\n"))

cat("Aggregating reads by fragment key and barcode...\n")
dt[, fragment_key := paste(chr, start, strand, sep = ":")] # 5'prime 已在 start/end

# (注意: 我们假设 start/end 已经是 1bp 的 5'prime 坐标)
# 聚合 reads, 得到每个 (key, barcode) 组合的总数
key_summary <- dt[, .(read_count = sum(count)), by = .(fragment_key, barcode, x, y)]

cat("Identifying dominant barcode for each fragment key...\n")
# 1. 找出每个 key 的最大 read 数
key_summary[, max_count := max(read_count), by = fragment_key]

# 2. 标记所有“获胜者”
dominant_barcodes <- key_summary[read_count == max_count]

# 3. 处理平局：如果多个 barcode "获胜"，我们只选第一个
# (这是 v2 脚本的逻辑, 但在这里我们只用它来*定义*源头)
dominant_map <- dominant_barcodes[, .SD[1], by = fragment_key]
dominant_map <- dominant_map[, .(fragment_key, dominant_x = x, dominant_y = y)]

cat("Dominant barcode map created.\n")

cat("Merging maps and calculating distances for non-dominant reads...\n")
# 将 "获胜者" 的坐标合并回 "所有 reads" 的表
all_reads_with_dominant <- merge(key_summary, dominant_map, by = "fragment_key")

# 找出 *所有* 非优势 reads (即 X 或 Y 坐标不匹配的)
non_dominant_reads <- all_reads_with_dominant[x != dominant_x | y != dominant_y]

cat(paste("Found", non_dominant_reads[, sum(read_count)], "non-dominant reads.\n"))

# 计算欧氏距离 (单位: bins)
non_dominant_reads[, dist_bins := sqrt((x - dominant_x)^2 + (y - dominant_y)^2)]
# 转换为物理距离 (um)
non_dominant_reads[, dist_um := dist_bins * file_resolution]

cat("Generating diffusion histogram...\n")
# 对距离进行分箱 (例如, 0-16um, 16-32um, ...)
bins <- seq(0, Args$max_dist_um, by = Args$bin_width_um)
non_dominant_reads[, dist_bin_label := cut(dist_um, breaks = bins, right = FALSE)]

# 过滤掉NA (即 > max_dist_um 的点)
non_dominant_reads <- non_dominant_reads[!is.na(dist_bin_label)]

# 按距离分箱汇总
diffusion_histogram <- non_dominant_reads[, 
  .(total_diffused_reads = sum(read_count)), 
  by = .(dist_bin_label)
]

total_reads <- dt[, sum(count)]
diffusion_histogram[, fraction_of_total := total_diffused_reads / total_reads]

setkey(diffusion_histogram, dist_bin_label)

cat(paste("Saving results to:", output_tsv, "\n"))
fwrite(diffusion_histogram, file = output_tsv, sep = "\t")

cat(paste("Saving plot to:", output_pdf, "\n"))

# (为了绘图，我们需要将 factor 转回数值，取区间的起始点)
diffusion_histogram_plot <- diffusion_histogram %>%
  mutate(
    dist_start = as.numeric(str_extract(gsub("\\[|\\)", "", dist_bin_label), "^[0-9\\.]+"))
  )

g <- ggplot(diffusion_histogram_plot, aes(x = dist_start, y = total_diffused_reads)) +
  geom_col(fill = "#B3CDE3", color = "black") + # 柱状图
  scale_x_continuous(
    breaks = seq(0, Args$max_dist_um, by = Args$bin_width_um * 4), # 每 4 个 bin 显示一个标签
    name = "Distance from Dominant Bin (μm)"
  ) +
  scale_y_log10( # 强烈推荐使用 log 刻度来查看长尾
    name = "Total Diffused Reads (Log Scale)"
  ) +
  labs(
    title = "Spatial Diffusion Profile",
    subtitle = paste("Total non-dominant reads (NDRR) vs. Distance from origin (", file_resolution, "um bins)", sep = "")
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(output_pdf, g, width = 10, height = 6)

cat("--- R script execution completed ---\n")
