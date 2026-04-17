library(tidyverse)
library(patchwork)

args = commandArgs(T)
file_raw = args[1]
file_UMI = args[2]
file_pdf = args[3]

process_file <- function(file_name, color, title_name){
    raw_reads <- read_delim(file_name, delim = "\t", col_names = F)
    colnames(raw_reads) <- c("chr", "start", "end", "cell_barcode", "count", "strand")
    raw_reads %>%
        mutate(group_bin = case_when(
            count == 1 ~ "1",
            count == 2 ~ "2",
            count == 3 ~ "3",
            count == 4 ~ "4",
            count > 4 & count <= 10 ~ "[5, 10]",
            count > 10 & count <= 20 ~ "[11, 20]",
            count > 20 & count <= 30 ~ "[21, 30]",
            count > 30 & count <= 40 ~ "[31, 40]",
            count > 40 & count <= 50 ~ "[41, 50]",
            count > 50 & count <= 100 ~ "[51, 100]",
            count > 100 ~ "> 100"    )) -> raw_reads
    
    raw_reads %>%
        group_by(group_bin) %>%
        summarise(frequency_count =  sum(count)) %>%
        ungroup() %>%
        mutate(percentage = frequency_count / sum(frequency_count) * 100) -> raw_reads.df

    raw_reads.df$group_bin <- factor(raw_reads.df$group_bin,
         levels = c("1", "2", "3", "4",  
                   "[5, 10]", "[11, 20]",
                    "[21, 30]", "[31, 40]",
                    "[41, 50]", "[51, 100]", "> 100"))

    raw_reads.qc <- ggplot(raw_reads.df, aes(x=group_bin, y=percentage)) +
        geom_col(color='black',position = position_dodge(width = 1), width = 0.7, fill = color) +
        geom_text(aes(label=paste0(round(percentage,2),"%"), y=percentage+2), position=position_dodge(0.9), vjust=0.1) +
        xlab('')+ylab("%")+ labs(title = title_name)+
        scale_y_continuous(limits = c(0, 100)) + 
    theme_bw()+theme(panel.grid = element_line(colour = 'white'),
                    legend.position = 'none',
                        axis.title.x =element_text(size=16, vjust = -0.2), 
                    axis.title.y=element_text(size=14),
                    axis.text.x=element_text(size=14,colour = "black", angle = 20, vjust = -0.02))
    
    return(raw_reads.qc)

}

raw.plot <- process_file(file_raw, "orange", "Raw read count per fragment per spot")
UMI.plot <- process_file(file_UMI, "firebrick3", "UMI count per fragment per spot")

ggsave(raw.plot / UMI.plot, filename = file_pdf, width = 150, height = 320, units = "mm", dpi = 600)