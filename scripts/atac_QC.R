library(tidyverse)
library(patchwork)

args = commandArgs(T)
min_contacts = as.numeric(args[1])
tss_file = args[2]
out_pdf = args[3]
in_tissue = as.logical(args[4]) # true or false
barcode_file = args[5]

if(is.na(min_contacts)){
     min_contacts = 1000
}
if(is.na(in_tissue)){
     in_tissue = FALSE
}
if(is.na(in_tissue)){
     in_tissue = FALSE
}
df <- read_delim(tss_file, delim = "\t")
colnames(df)[1] = "barcode"

print(paste0("min_contacts: ", min_contacts))

df %>%
     filter(n_fragment > min_contacts) -> df

print(paste0("Row of df: ", nrow(df)))

plot_stat <- function(df){
     frag_count <- ggplot(df, aes(x = "", y=log10(n_fragment))) + 
     geom_violin( scale = 'width', bw = 0.5, fill = "#219ebc") + 
     geom_boxplot(outlier.size = 0, width = 0.12, show.legend = F, fill = "#219ebc") +  
     labs(y="", 
          x="Unique fragment counts",
          subtitle = paste0("Median: ", median(df$n_fragment)),
          fill="") +
     theme_bw() +
     theme(legend.position = "none",
               axis.title.x =element_text(size=18), 
               axis.title.y =element_text(size=18),
               axis.text.x = element_text(size=10,colour = "black"),
               axis.text.y = element_text(size=10,colour = "black"),
               panel.grid = element_line(colour = 'white')) 

     frac_dup <- ggplot(df, aes(x = "", y=frac_dup * 100)) + 
     geom_violin( scale = 'width', bw = 0.5,fill = "#ffb703") + 
     geom_boxplot(outlier.size = 0, width = 0.12, show.legend = F,fill = "#ffb703") +  
     labs(y="%", 
          x="Duplication rate",
               subtitle = paste0("Median: ", round(median(df$frac_dup),4) * 100),
          fill="") +
     theme_bw() +
     theme(legend.position = "none",
               axis.title.x =element_text(size=18), 
               axis.title.y =element_text(size=18),
               axis.text.x = element_text(size=10,colour = "black"),
               axis.text.y = element_text(size=10,colour = "black"),
               panel.grid = element_line(colour = 'white')) 

     frac_mito <- ggplot(df, aes(x = "", y=frac_mito * 100)) + 
     geom_violin( scale = 'width', bw = 0.5,fill = "#8ecae6") + 
     geom_boxplot(outlier.size = 0, width = 0.12, show.legend = F,fill = "#8ecae6") +  
     labs(y="%", 
          x="Mitocondria fragments",
          fill="") +
     theme_bw() +
     theme(legend.position = "none",
               axis.title.x =element_text(size=18), 
               axis.title.y =element_text(size=18),
               axis.text.x = element_text(size=10,colour = "black"),
               axis.text.y = element_text(size=10,colour = "black"),
               panel.grid = element_line(colour = 'white')) 

     frac_tss <- ggplot(df, aes(x = "", y=tsse )) + 
     geom_violin( scale = 'width', bw = 0.5,fill = "#fb8500") + 
     geom_boxplot(outlier.size = 0, width = 0.12, show.legend = F,fill = "#fb8500") +  
     labs(y="", 
          x="TSS enrichment score",
               subtitle = paste0("Median: ", round(median(df$tsse),2)),
          fill="") +
     theme_bw() +
     theme(legend.position = "none",
               axis.title.x =element_text(size=18), 
               axis.title.y =element_text(size=18),
               axis.text.x = element_text(size=10,colour = "black"),
               axis.text.y = element_text(size=10,colour = "black"),
               panel.grid = element_line(colour = 'white')) 

     ggsave(frac_tss + frac_mito + frag_count  + frac_dup, filename = out_pdf, width = 200, height = 200, units = "mm", dpi = 600)
}

if(in_tissue){
     barcode_in_tissue = read_delim(barcode_file, delim = "\t", col_names = F)
     colnames(barcode_in_tissue) <- c("barcode")
     df %>%
          inner_join(barcode_in_tissue) -> df.in_tissue

     plot_stat(df.in_tissue)

}else{
     plot_stat(df)
}
