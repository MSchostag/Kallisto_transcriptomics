rm(list=ls()); gc()

####   DESeq2 analysis   ####

#### Installing needed R-packeages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rhdf5")
BiocManager::install("tximportData")
BiocManager::install("DESeq2")

if (!require("pacman")) install.packages("pacman")

# https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto

####   Loading all need packages #### 

pacman::p_load(here,stringr,vctrs,vegan,pheatmap,readxl,ggplot2,tximport,tidyverse,DESeq2,dplyr, install = T)


#####Creating folders one would need####

here() # check that you are the right place
getwd()
dir.create(here("figs"),showWarnings = F)
dir.create(here("data"),showWarnings = F)

#### loading the data ####

#the folder were all the files are located 
dir <- here()
list.files(dir) #make sure the right folders appears

samples <- read.table(file.path(dir, "S26samples.txt"), header = TRUE)
samples
meta <- read.table(file.path(dir, "meta.txt"), header = TRUE)
meta
rownames(meta) <- meta$sample

files <- file.path(dir, "kallisto", samples$run, "abundance.h5")
files

# here you need to fill in the all your samples names - You could copy it all from the sample file 

names(files) <- meta$sample


files

txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)

# check that the input files looks okay
head(txi.kallisto$counts)

# making data frame for DeSeq2
dds_data <- DESeqDataSetFromTximport(txi.kallisto, 
                                     colData=meta,
                                     ~condition)
# creating DeSeq2 analysis and object
dds <- DESeq(dds_data)



      ####SizeFactor####

# https://support.bioconductor.org/p/97676/

nm <- assays(dds)[["avgTxLength"]]
sf <- estimateSizeFactorsForMatrix(counts(dds) / nm)
sf
      ####Results from DeSeq2####

res <- results(dds)
res <- res[order(res$padj),]
summary(res)

plotDispEsts(dds) 

number_hits <-  as.data.frame(txi.kallisto$counts)

colSums(number_hits)

      ####PCA plot####


vsdata <- vst(dds, blind=FALSE)
plotPCA(object = vsdata, intgroup="condition") + 
  theme_classic() + 
  geom_point(size=4)
ggsave(here("figs", "PCA_plot.png"))

      ####Making contrasts####

# THIS WHOLE SECTION NEEDS TO ADJUSTED FOR your own ANALYSIS

#Mutant 24 h vs Mutant 72 h 
res_m_24_m_72 <- results(dds, contrast=c("condition","m_24","m_72"), tidy = T)
#Mutant 24 h vs Wild type 24 h
res_m_24_wt_24 <- results(dds, contrast=c("condition","m_24","wt_24"), tidy = T)
#Wild type 24 h vs Wild type 72 h
res_wt_24_wt_72 <- results(dds, contrast=c("condition","wt_24","wt_72"), tidy = T)
#Mutant 24 h vs Wild type 72 h
res_m_72_wt_72 <- results(dds, contrast=c("condition","m_72","wt_72"), tidy = T)


# Again you would need to fill in the names of the contrast data.frames and their names in the following two lines
contrast_list <- list(res_m_24_m_72,res_m_24_wt_24,res_wt_24_wt_72,res_m_72_wt_72)
contrast_names <- c("Mutant_24 vs Mutant_72","Mutant_24 vs Wildtype_24","Wildtype_24 vs Wildtype_72","Mutant_72 vs Wildtype_72")

all_contrast <- map2(contrast_list, contrast_names, function(x, name) {
  x %>% 
    mutate(contrast = name) %>% as_tibble()
})



all_contrast_filtered <- do.call("rbind", all_contrast) %>% 
  mutate(sig = case_when(padj > 0.01 ~ "not_rel",
                         log2FoldChange > 2 & padj < 0.01 ~ "up_sig",
                         log2FoldChange < -2 & padj < 0.01 ~ "down_sig",
                         log2FoldChange < 2 & padj < 0.01 & log2FoldChange > -2 ~ "not_rel")) %>% 
  filter(padj!="") %>% 
  as_tibble()


map2(all_contrast, contrast_names, function(contrast, names){
  write.csv(contrast, here("data", paste0(names, ".csv" )))
  
})
    

     ####  Volcano plots  ####

volcano_plots <- all_contrast_filtered %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col = sig)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values=c("Red", "black", "blue")) +
  facet_wrap(.~contrast)
# save the figure in figure folder
ggsave(here("figs", paste0("volcano_plots", ".png")))


     #### Plot the gene of interest ####


plot_gene <- function(genename){
  plotCounts(dds,
             gene = genename,
             intgroup = "condition", 
             returnData = T) %>% 
    ggplot(aes(x=condition, y=count)) + 
    geom_point(size=4, position = position_jitter()) + 
    ggtitle(genename) + 
    theme_classic() + 
    ylab("Normalized Count") + 
    xlab("Group") + 
    theme(plot.title = element_text(hjust = 0.5, face="bold"), axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"))  
}

plot_gene("gene-OL67_RS18380") 

      ####  most upregulated  ####

up_sig_list <- all_contrast_filtered %>% 
  group_by(contrast) %>% 
  filter(sig =="up_sig") %>% 
  arrange(desc(baseMean), .by_group = T, decrease = T) %>% 
  dplyr::slice(1:9)

down_sig_list <- all_contrast_filtered %>% 
  group_by(contrast) %>% 
  filter(sig =="down_sig") %>% 
  arrange(desc(baseMean), .by_group = T, decrease = T) %>% 
  dplyr::slice(1:9)


count_plots_up <- map(pull(distinct(up_sig_list, contrast)), function(x){
  genes <- up_sig_list %>% 
    filter(contrast == x) %>% 
    pull(row)
  map(genes, function(y){
    plotCounts(dds,
               gene = y,
               intgroup = "condition",
               returnData = T) %>% 
      mutate(gene = y)
  })
})

map2(count_plots_up, pull(distinct(up_sig_list, contrast)), function(x, name) {
  do.call("rbind", x) %>% 
    rownames_to_column(var = "group") %>% 
    ggplot(aes(x = condition, y = count)) +
    geom_point(position = position_jitter()) +
    theme_bw()+
    facet_wrap(.~gene, scales = "free") +
    labs(title = paste0("Top 9 Up regulated genes in ", name))  +
    xlab("Group") + 
    ylab("Normalized Count")
    ggsave(here("figs", paste0("up_regulated_",name, ".png")))
})

count_plots_down <- map(pull(distinct(down_sig_list, contrast)), function(x){
  genes <- up_sig_list %>% 
    filter(contrast == x) %>% 
    pull(row)
  map(genes, function(y){
    plotCounts(dds,
               gene = y,
               intgroup = "condition",
               returnData = T) %>% 
      mutate(gene = y)
  })
})

map2(count_plots_down, pull(distinct(down_sig_list, contrast)), function(x, name) {
  do.call("rbind", x) %>% 
    rownames_to_column(var = "group") %>% 
    ggplot(aes(x = condition, y = count)) +
    geom_point(position = position_jitter()) +
    theme_bw()+
    facet_wrap(.~gene, scales = "free") +
    ylab("Normalized Count") + 
    xlab("Group") + 
    labs(title = paste0("Top 9 Down regulated genes in ", name))  
  ggsave(here("figs", paste0("down_regulated_",name, ".png")))
})



####Ipath_analysis####

# you would need to type in the filename of the eggnog mapper output file 

eggnog_data <- read_tsv(file.path(here(), "MM_dtjj1h_k.emapper.annotations.tsv")) %>% # type between ""
  dplyr::rename(row = query)

all_contrast_filtered_ipath <- do.call("rbind", all_contrast) %>% 
  mutate(sig = case_when(padj > 0.01 ~ "not_rel",
                         log2FoldChange > 0.5 & padj < 0.01 ~ "up_sig",
                         log2FoldChange < -0.5 & padj < 0.01 ~ "down_sig",
                         log2FoldChange < 0.5 & padj < 0.01 & log2FoldChange > -0.5 ~ "not_rel")) %>% 
  filter(padj!="") %>% 
  as_tibble()



all_contrast_filtered_ipath_f <- all_contrast_filtered_ipath %>%
  filter(sig != "not_rel") %>% 
  group_by(contrast) %>% 
  group_split() %>%
  setNames(unique(all_contrast_filtered_ipath$contrast))


all_contrast_filtered_ipath_ready <- map(all_contrast_filtered_ipath_f, function(contrast){
  contrast %>% 
    left_join(eggnog_data, by = "row")
})


# Write iPath-ready table

genes_to_iPath <- function(gene_list, file_name_append = NULL) {
  map2(gene_list, names(gene_list), function(contrast, name){
    int_list <- contrast %>% 
      dplyr::select(KEGG_ko, sig) %>% 
      filter(!is.na(KEGG_ko), !grepl("-", KEGG_ko)) %>% 
      mutate(KEGG_ko = strsplit(KEGG_ko, ",")) %>%
      unnest(KEGG_ko) %>% 
      mutate(KEGG_ko = str_remove(KEGG_ko, "ko:")) %>%
      mutate(sig = ifelse(sig == "up_sig", "#038BE7", "#EC3804"),
             width = "W10") # make column named change that gives you up and down reg genes
    
    write_tsv(int_list, here("data", 
                             paste0("iPath_", name, "_", file_name_append, ".tsv")), 
              col_names = F)
  })
}

genes_to_iPath(all_contrast_filtered_ipath_ready)



####COG categories plot####

# create COG table for each contrast with significant genes + saving table with all data

cog_tally <- map2(all_contrast_filtered_ipath_ready, contrast_names, function(x, names){
  cog_tally <- x %>% 
    dplyr::filter(sig =="up_sig" | sig =="down_sig") %>% 
    filter(!is.na(COG_category), 
           !COG_category=='-') %>%
    mutate(COG_category = strsplit(COG_category, "")) %>%
    unnest(COG_category)
  
  write_tsv(cog_tally, here("data", paste0("COG_data_", names, ".tsv" )))
  
  cog_tally <- cog_tally %>% 
    dplyr::select(COG_category,sig) %>%   
    group_by(COG_category, sig) %>% 
    tally(name = "count")
  
  
  eggnog <- eggnog_data %>% 
    select(COG_category) %>%
    filter(!is.na(COG_category), 
           !COG_category=='-') %>% 
    mutate(COG_category = strsplit(COG_category, "")) %>%
    unnest(COG_category) %>%
    group_by(COG_category) %>% 
    tally(name = "total_COG")
  
  cog_tally <- cog_tally %>% 
    left_join(eggnog, by = "COG_category") %>% 
    mutate(relative = count/total_COG * 100,
           count = ifelse(sig == "down_sig", count*-1, count),
           relative = ifelse(sig == "down_sig", relative*-1, relative))
    })



map2(cog_tally, contrast_names, function(x, names){
    x %>%   
    ggplot(aes(fill=sig, x=COG_category, y=relative)) +
    geom_col() + 
    theme_classic() +
    coord_flip() +
    xlab("COG categories") +
    scale_x_discrete(limits=rev) + 
    scale_fill_manual(values=c("#EC3804" , "#038BE7")) +
    scale_y_continuous(limits = c(-100, 100)) +
    xlab("COG categories") +
    ggtitle(paste0("DE COG categories in", " ", names)) +
    theme(plot.title = element_text(hjust = 0.5, face="bold"))
  
  ggsave(here("figs", paste0("Relative DE COG_categories in",names, ".png")))
  
  x %>%   
    ggplot(aes(fill=sig, x=COG_category, y=count)) +
    geom_col() + 
    theme_classic() +
    coord_flip() +
    xlab("COG categories") +
    scale_x_discrete(limits=rev) + 
    scale_fill_manual(values=c("#EC3804" , "#038BE7")) +
    xlab("COG categories") +
    ggtitle(paste0("DE COG categories in", " ", names)) +
    theme(plot.title = element_text(hjust = 0.5, face="bold"))
  
  ggsave(here("figs", paste0("DE COG_categories in",names, ".png")))
})
























