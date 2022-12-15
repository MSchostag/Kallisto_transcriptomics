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

here() # checke that you are the right place
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

head(txi.kallisto$counts)


dds_data <- DESeqDataSetFromTximport(txi.kallisto, 
                                     colData=meta,
                                     ~condition)

dds <- DESeq(dds_data)


      ####SizeFactor####

# https://support.bioconductor.org/p/97676/

nm <- assays(dds)[["avgTxLength"]]
sf <- estimateSizeFactorsForMatrix(counts(dds) / nm)

      ####Results from DeSeq2####

res <- results(dds)
res <- res[order(res$padj),]
summary(res)

plotDispEsts(dds) 

      ####PCA plot####


vsdata <- vst(dds, blind=FALSE)
plotPCA(object = vsdata, intgroup="condition") + 
  theme_classic() + 
  geom_point(size=4)
ggsave(here("figs", "PCA_plot.png"))

      ####Making contrasts####


#Mutant 24 h vs Mutant 72 h 
res_m_24_m_72 <- results(dds, contrast=c("condition","m_24","m_72"), tidy = T)
#Mutant 24 h vs Wild type 24 h
res_m_24_wt_24 <- results(dds, contrast=c("condition","m_24","wt_24"), tidy = T)
#Wild type 24 h vs Wild type 72 h
res_wt_24_wt_72 <- results(dds, contrast=c("condition","wt_24","wt_72"), tidy = T)
#Mutant 24 h vs Wild type 72 h
res_m_72_wt_72 <- results(dds, contrast=c("condition","m_72","wt_72"), tidy = T)


# 

contrats_list <- list(res_m_24_m_72,res_m_24_wt_24,res_wt_24_wt_72,res_m_72_wt_72)
contrats_names <- c("Mutant_24 vs Mutant_72","Mutant_24 vs Wildtype_24","Wildtype_24 vs Wildtype_72","Mutant_72 vs Wildtype_72")

all_contrast <- map2(contrats_list, contrats_names, function(x, name) {
  x %>% 
    mutate(contrast = name) %>% as_tibble()
})

all_contrats_filtered <- do.call("rbind", all_contrast) %>% 
  mutate(sig = case_when(padj > 0.01 ~ "not_rel",
                         log2FoldChange > 2 & padj < 0.01 ~ "up_sig",
                         log2FoldChange < -2 & padj < 0.01 ~ "down_sig",
                         log2FoldChange < 2 & padj < 0.01 & log2FoldChange > -2 ~ "not_rel")) %>% 
  filter(padj!="") %>% 
  as_tibble()



map2(all_contrast, contrats_names, function(contrast, names){
  write.csv(contrast, here("data", paste0(names, ".csv" )))
  
})
    

     ####  Volcano plots  ####

volcano_plots <- all_contrats_filtered %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col = sig)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values=c("Red", "black", "blue")) +
  facet_wrap(.~contrast)
# save the figure in figure folder
ggsave(here("figs", paste0("volcano_plots", ".png")))

      ####  most upregulated  ####

up_sig_list <- all_contrats_filtered %>% 
  group_by(contrast) %>% 
  filter(sig =="up_sig") %>% 
  arrange(desc(baseMean), .by_group = T, decrease = T) %>% 
  dplyr::slice(1:9)

down_sig_list <- all_contrats_filtered %>% 
  group_by(contrast) %>% 
  filter(sig =="down_sig") %>% 
  arrange(desc(baseMean), .by_group = T, decrease = T) %>% 
  dplyr::slice(1:9)




map(pull(distinct(up_sig_list, contrast)), function(x) {
  up_sig_list %>% 
    filter(contrast == x)
})
 #######http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html####


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
    labs(title = name)  
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
    labs(title = name)  
  ggsave(here("figs", paste0("down_regulated_",name, ".png")))
})








#_----------------

####Ipath_analysis####

setwd("C:/Users/mdesc/OneDrive - Danmarks Tekniske Universitet/1-Cemist/1-projects/12-transcriptomic analysis/eggnog mapper output/")
dir <- "C:/Users/mdesc/OneDrive - Danmarks Tekniske Universitet/1-Cemist/1-projects/12-transcriptomic analysis/eggnog mapper output/"
list.files(dir)
eggnog_data <- read.table(file.path(dir, "all_genes.txt"), header = TRUE) %>% 
  dplyr::rename(row = query)

all_contrats_filtered_ipath <- do.call("rbind", all_contrast) %>% 
  mutate(sig = case_when(padj > 0.01 ~ "not_rel",
                         log2FoldChange > 0.5 & padj < 0.01 ~ "up_sig",
                         log2FoldChange < -0.5 & padj < 0.01 ~ "down_sig",
                         log2FoldChange < 0.5 & padj < 0.01 & log2FoldChange > -0.5 ~ "not_rel")) %>% 
  filter(padj!="") %>% 
  as_tibble()



all_contrats_filtered_ipath_f <- all_contrats_filtered_ipath %>%
  filter(sig != "not_rel") %>% 
  group_by(contrast) %>% 
  group_split() %>%
  setNames(unique(all_contrats_filtered_ipath$contrast))


all_contrats_filtered_ipath_ready <- map(all_contrats_filtered_ipath_f, function(contrast){
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

genes_to_iPath(all_contrats_filtered_ipath_ready)





colSums(counts(dds, normalized=T))



library(pheatmap)

pheatmap(counts(dds), labels_col = meta$condition)

pheatmap(log10(counts(dds, normalized=T)+1), labels_col = meta$condition)




