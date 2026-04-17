library(tidyverse)
library(patchwork)
args = commandArgs(T)

stat = args[1]
outdir = args[2] # ./prefix

df <- read_delim(stat, delim = "\t", col_names = F)

p_density <- ggplot(df, aes(x=X1)) + 
  geom_density(color="firebrick3", linewidth = 1.2 ) + 
    theme_bw() +
    scale_x_continuous(limits = c(-10, 160), breaks = seq(0, 160, 20), labels = c("0", "20", "40", "60", "80", "100", "120", "140", "160")) + 
    labs(x = "PolyA trimmed R2 length", y = "Density") +
  theme(axis.text.y = element_text(size=8),
        axis.title.x =  element_text(size=14),
        axis.title.y =  element_text(size=14),
        axis.text.x = element_text(size=8))

result = c()
for(cutoff in seq(0, 150, 10)){
    result = c(result, round(sum(df$X1 >= cutoff) / nrow(df) * 100, 2))
}
cumulative_df = data.frame(cutoff = seq(0, 150, 10), fraction = result)

p_cumsum <- ggplot(cumulative_df, aes(x=cutoff, y = fraction)) + 
  geom_point(color="black", size = 1.5 ) + 
  geom_line(color="firebrick3", linewidth = 1.2 ) + 
    theme_bw() +
    scale_x_continuous(limits = c(-10, 160), breaks = seq(0, 160, 20), labels = c("0", "20", "40", "60", "80", "100", "120", "140", "160")) + 
    labs(x = ">= certain length", y = "Fraction (%)") +
  theme(axis.text.y = element_text(size=12),
        axis.title.x =  element_text(size=16),
        axis.title.y =  element_text(size=16),
        axis.text.x = element_text(size=12))
        

ggsave(p_density + p_cumsum, 
        filename = paste0(outdir, ".abortive_stat.pdf"),
        width = 240,
        height = 120,
        units = "mm",
        dpi = 600)
