library(tidyverse)
library(patchwork)

args = commandArgs(T)
indir = args[1] # directory stores.metadata.txt
sampleid = args[2]
file_pdf = args[3] # svg
in_tissue = as.logical(args[4]) # true or false
barcode_dir = args[5]
if(is.na(in_tissue)){
    in_tissue = FALSE
}

df = data.frame()
files = list.files(indir, pattern= "*in_tissue.*metadata.txt")
files_all = c(files,
    list.files(indir, pattern = "*SC.metadata.txt")
)
for(file_ in files_all){
    metadata = read_delim(paste0(indir, "/", file_), )
    colnames(metadata)[1] = "barcode"

    print(head(df))
    replce_string = paste0(sampleid, "_")
    label = gsub(".metadata.txt", "", gsub(replce_string, "", file_))
    label = gsub("in_tissue_", "", label)
    metadata$res = rep(label, nrow(metadata))

    if(in_tissue){
        level_num = gsub("L", "", label)
        if(level_num != "SC"){ 
            if(level_num == "B"){
                level_num = "13"
            }
            barcode_path = paste0(barcode_dir, "/level_", level_num, "/barcodes_in_tissue.tsv.gz")
            print(paste0("Extract barcodes in tissue: ", barcode_path))
            barcode_in_tissue = read_delim(barcode_path, delim = "\t", col_names = F)
            colnames(barcode_in_tissue) <- c("barcode")
            metadata %>%
                inner_join(barcode_in_tissue) -> metadata
        }
    }
    df = bind_rows(metadata, df)
}

df %>%
    mutate(resolution = case_when(
        res == "SC" ~ "SC",
        res == "LB" ~ "100μm",
        res == "L7" ~ "50μm",
        res == "L6" ~ "42μm",
        res == "L5" ~ "35μm",
        res == "L4" ~ "27μm",
        res == "L3" ~ "20μm",
        res == "L2" ~ "10μm",
        res == "L1" ~ "5μm")) -> df

df$resolution <- factor(df$resolution, 
                        levels = c("SC","5μm", "10μm", "20μm",
                                  "27μm", "35μm",
                                  "42μm", "50μm",
                                  "100μm"))

df$unique_frag <- log10(df$n_fragment)

df %>%
    group_by(resolution) %>%
    summarise(mean_labe = median(n_fragment)) %>%
    ungroup() -> df.mean
    
df.mean$median_unique_frag <- log10(df.mean$mean_labe)

n_frag.plot <-  ggplot(df, aes(x=resolution, y=unique_frag, fill=resolution)) +
    geom_violin(position = position_dodge(width=1.2), scale = 'width', bw = 0.4) + 
    geom_boxplot(position = position_dodge(width=1.2),
               outlier.size = 0, width = 0.3, show.legend = F)+
    geom_text(inherit.aes = FALSE, 
                data = df.mean,
            aes(label=round(mean_labe,0), x = resolution, y=median_unique_frag+1.4), vjust=0.1) +
    scale_fill_manual(values=c(
        "SC" = "#80989b",
        "5μm"="#eae2b7",
        "10μm"="#fcbf49",
        "20μm"="#f77f00",
        "27μm"="#d62828",
        "35μm"="#e26d5c",
        "42μm"="#ecf39e",
        "50μm"="#90a955",
        "100μm"="#4f772d"
        )) +
    scale_x_discrete(
        labels = c(
            "SC" = "SC",
            "5μm" = expression(paste("5",mu, "m")),
            "10μm" = expression(paste("10",mu, "m")),
            "20μm" = expression(paste("20",mu, "m")),
            "27μm" = expression(paste("27",mu, "m")),
            "35μm" = expression(paste("35",mu, "m")),
            "42μm" = expression(paste("42",mu, "m")),
            "50μm" = expression(paste("50",mu, "m")),
            "100μm" = expression(paste("100",mu, "m"))))+
    xlab('')+ylab("log10 Unique Fragments")+ labs("")+
  scale_y_continuous(limits = c(1, 6.2), breaks = c(1, 2, 3, 4, 5), labels = c("1", "2", "3", "4", "5"))+
theme_bw()+theme(panel.grid = element_line(colour = 'white'),
                   legend.position = "None",
                   axis.title.y=element_text(size=16),
                 axis.title.x=element_blank(),
                   axis.text.x= element_text(angle = 0, vjust = 1))


ggsave(filename = file_pdf,
       n_frag.plot,
       width=140,height = 100, units = "mm",dpi = 600)
