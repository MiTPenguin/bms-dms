### TYK2 IL-23 RNA-seq

1. [Dataset Properties](#part1)
2. [Unsupervised Profiles](#part2)
3. [Statistical Modeling](#part3)
4. [Differentially Expressed Genes](#part4)
5. [Gene Visualizations](#part5)


```R
library(ggpubr)
library(DESeq2)
library(ggcorrplot)
library(furrr)
library(ggbiplot)
library(patchwork)
library(ComplexHeatmap)
library(tidyverse)
```

#### Dataset Properties <a name="part1"></a>


```R
alignments <- read_tsv("../pipeline/diffexp/multiqc_data/multiqc_star.txt") %>%
    mutate(id = as.numeric(gsub("_.*", "", Sample)))

count_summary <- read_tsv("../pipeline/diffexp/multiqc_data/multiqc_featureCounts.txt")  %>%
    mutate(id = as.numeric(gsub("_.*", "", Sample)))
```


```R
options(repr.plot.width = 8, repr.plot.height = 18)
alignments %>%
    ggplot() +
        geom_bar(aes(x = total_reads, y = reorder(Sample, -id)), stat = "identity") +
        theme_pubr(base_size = 16, x.text.angle = 45) +
        xlab("Total Read Depth") + ylab("") +
        ggtitle("Total Read Depth")
```


    
![png](IL23-RNAseq_files/IL23-RNAseq_4_0.png)
    



```R
alignments %>%
    ggplot() +
        geom_bar(aes(x = uniquely_mapped_percent/100, y = reorder(Sample, -1*as.numeric(id))), stat = "identity") +
        theme_pubr(base_size = 16, x.text.angle = 45) +
        xlab("Unique Alignment Rate") + ylab("") +
        ggtitle("Unique Alignment Rate")
```


    
![png](IL23-RNAseq_files/IL23-RNAseq_5_0.png)
    



```R
count_summary %>%
    ggplot() +
        geom_bar(aes(x = percent_assigned/100, y = reorder(Sample, -1*as.numeric(id))), stat = "identity") +
        theme_pubr(base_size = 16, x.text.angle = 45) +
        xlab("Gene Assignment Rate") + ylab("") +
        ggtitle("Gene Assignment Rate")
```


    
![png](IL23-RNAseq_files/IL23-RNAseq_6_0.png)
    


#### Unsupervised Profiles <a name="part2"></a>


```R
paths <- str_c("../pipeline/diffexp/",
                    alignments$Sample,
                    ".counts.tsv")
names(paths) <- paths
df <- map_dfr(paths,
              ~read_tsv(.x,
                        col_names = NULL,
                        skip = 2,
                        col_select = c(1,7)),
              .id = "sample") %>%
    rename("gene" = "X1",
           "count" = "X7") %>%
    mutate(sample = gsub("../pipeline/diffexp/|.counts.tsv", "", sample))

samp_prop <- read_tsv("../sumstats/RNASEQ/run2/sample-properties-il23.tsv") %>%
    arrange(sample_id) %>%
    mutate(rep = rep(c(1,2), 40),
           covariate = as.factor(str_c(group, "_", cytokine, "_", dosage, "_", time)))
rownames(samp_prop) <- samp_prop$sample_id

df_id <- df %>%
    mutate(id = as.numeric(gsub("_.*", "", sample))) %>%
    select(-sample) %>%
    arrange(id)
```


```R
cor_data_wide <- inner_join(df_id, samp_prop, by = c("id" = "sample_id")) %>%
    filter(time == "t6", group == "WT",  count != 0) %>%
    arrange(cytokine, group, dosage, rep) %>%
    select(gene, covariate, rep, count) %>%
    mutate(count = log2(count)) %>%
    pivot_wider(names_from = covariate:rep, values_from = count)

cor_data_wide$v <- rowMeans(as.matrix(cor_data_wide[,-1]), na.rm = TRUE)

cor_data_filt <- cor_data_wide %>%
    arrange(-v) %>%
    head(1000) %>%
    select(-gene, -v)

cor_mat <- cor(cor_data_filt, use = "pairwise.complete.obs")

options(repr.plot.width = 10, repr.plot.height = 10)
ggcorrplot(cor_mat, lab = TRUE, lab_size = 4, hc.order = FALSE, type = "upper", show.diag = TRUE) +
    scale_fill_gradient2(limit = c(0.5,1), low = "blue", high =  "red", mid = "white", midpoint = 0.75)
```

    [1m[22mScale for [32mfill[39m is already present.
    Adding another scale for [32mfill[39m, which will replace the existing scale.



    
![png](IL23-RNAseq_files/IL23-RNAseq_9_1.png)
    



```R
pca_data <- inner_join(df_id, samp_prop, by = c("id" = "sample_id")) %>%
    filter(!is.na(covariate)) %>%
    arrange(cytokine, group, dosage, rep) %>%
    select(gene, covariate, rep, count) %>%
    mutate(count = log2(count + 1)) %>%
    pivot_wider(names_from = covariate:rep, values_from = count)

pca_data <- t(as.matrix(pca_data[,-1]))
```


```R
pca1 <- ggbiplot(prcomp(pca_data[,colVars(pca_data) > 0]),labels.size = 3,
         var.axes = FALSE,
         labels = rownames(pca_data),
         groups = gsub("_1|_2", "",  rownames(pca_data)),
         choices = 1:2,) +
        ggtitle("PC1 // PC2")

pca2 <- ggbiplot(prcomp(pca_data[,colVars(pca_data) > 0]),labels.size = 3,
         var.axes = FALSE,
         labels = rownames(pca_data),
         groups = gsub("_1|_2", "",  rownames(pca_data)),
         choices = 2:3,) +
        ggtitle("PC2 // PC3")

pca3 <- ggbiplot(prcomp(pca_data[,colVars(pca_data) > 0]),labels.size = 3,
         var.axes = FALSE,
         labels = rownames(pca_data),
         groups = gsub("_1|_2", "",  rownames(pca_data)),
         choices = 3:4,) +
        ggtitle("PC3 // PC4")

options(repr.plot.width = 15, repr.plot.height = 5, warn = -1)
pca1 + pca2 + pca3  &
  theme_minimal(base_size = 14) &
  theme(legend.position = "none") &
  coord_equal(ratio = 0.5)
```

    [1m[22mCoordinate system already present. Adding new coordinate system, which will
    replace the existing one.
    [1m[22mCoordinate system already present. Adding new coordinate system, which will
    replace the existing one.
    [1m[22mCoordinate system already present. Adding new coordinate system, which will
    replace the existing one.



    
![png](IL23-RNAseq_files/IL23-RNAseq_11_1.png)
    


#### Statistical Modeling <a name="part3"></a>


```R
de_data <- df_id %>%
    filter(id %in% 1:56) %>%
    pivot_wider(names_from = id, values_from = count)

de_prop <- samp_prop %>% filter(sample_id %in% 1:56)

# de_prop[21,1] <- 49
# de_prop[49,1] <- 21
# de_prop <- de_prop %>% arrange(sample_id)

rownames(de_prop) <- de_prop$sample_id
```


```R
raw_counts <- de_data
names(raw_counts)[-1] <- str_c(de_prop$covariate, "_", de_prop$rep)
write_tsv(raw_counts, "../sumstats/RNASEQ/run2/raw-counts.tsv")
```


```R
deobj_counts <- DESeqDataSetFromMatrix(countData = de_data %>% select(-gene),
    colData = de_prop,
    design = ~covariate)

deresult <- DESeq(deobj_counts)
```


```R
cov_group <- de_prop %>%
    filter(!grepl("none", covariate)) %>%
    distinct(covariate) %>%
    pull(covariate) %>%
    as.character()

none_group <- gsub("IFNa|IL10|IL23|low|high", "none", cov_group) %>%
    as.character()

plan(multicore, workers = 25)
norm_result <- future_map2(.x = cov_group,
                    .y = none_group,
                    ~results(deresult,
                             contrast = c("covariate", .x, .y),
                             independentFiltering = FALSE))

sumstats <- map2_dfr(norm_result,
                     cov_group,
                     ~bind_cols("gene" = de_data$gene,
                                            as_tibble(.x),
                                            "condition" = .y))

sumstats %>%
    separate(condition, c("background", "cytokine", "dosage", "time"), "_") %>%
    write_tsv("../sumstats/RNASEQ/run2/deseq2-sumstats-vs-none.tsv")
```


```R
sumstats <- read_tsv( "../sumstats/RNASEQ/run2/deseq2-sumstats-vs-none.tsv")
split_sumstats <- sumstats %>%
    mutate(group = case_when(padj == 1 ~ "FDR = 1",
                             padj < 0.01 ~ "FDR < 0.01",
                             TRUE ~ "NS")) %>%
    mutate(time = relevel(as.factor(time), ref = "t6"),
           dosage = relevel(as.factor(dosage), ref = "low"))

ma_grid <- split_sumstats %>%
    ggplot() +
        geom_point(aes(x = log2(baseMean),
                       y = log2FoldChange,
                       color = group)) +
        theme_pubr(base_size = 15) +
        facet_grid(rows = vars(dosage, time),
                   cols = vars(cytokine, background)) +
        scale_color_manual(values = c("FDR = 1" = "gray",
                                      "NS" = "black",
                                      "FDR < 0.01" = "red"))
```


```R
options(repr.plot.width = 20, repr.plot.height = 20, warn = -1)
ma_grid
```


    
![png](IL23-RNAseq_files/IL23-RNAseq_18_0.png)
    



```R
volcano_grid <- split_sumstats %>%
    ggplot() +
        geom_point(aes(x = log2FoldChange,
                       y = -log10(pvalue),
                       color = group)) +
        theme_pubr(base_size = 15) +
        facet_grid(rows = vars(dosage, time),
                   cols = vars(cytokine, background)) +
        scale_color_manual(values = c("FDR = 1" = "gray",
                                      "NS" = "black",
                                      "FDR < 0.01" = "red")) +
        coord_cartesian(ylim = c(0, 50))

options(repr.plot.width = 20, repr.plot.height = 20, warn = -1)
volcano_grid
```


    
![png](IL23-RNAseq_files/IL23-RNAseq_19_0.png)
    



```R
options(repr.plot.width = 10, repr.plot.height = 5)
split_sumstats  %>%
    filter(group == "FDR < 0.01") %>%
    group_by(background, cytokine, dosage, time) %>%
    count(group, .drop = FALSE) %>%
    ggplot() +
        geom_bar(aes(x = dosage,
                     y = n,
                     fill = time),
                 position = position_dodge(0.5),
                 width = 0.5,
                 stat = "identity")  +
        facet_grid(cols = vars(cytokine, background)) +
        theme_pubr(base_size = 16)

split_sumstats  %>%
    filter(group == "FDR < 0.01") %>%
    group_by(background, cytokine, dosage, time) %>%
    count(group, .drop = FALSE) %>%
    select(-group) %>%
    arrange(time, cytokine, rev(background), dosage) %>%
    rename("FDR < 0.01 vs Untreated" = "n")
```


<table class="dataframe">
<caption>A grouped_df: 24 × 5</caption>
<thead>
	<tr><th scope=col>background</th><th scope=col>cytokine</th><th scope=col>dosage</th><th scope=col>time</th><th scope=col>FDR &lt; 0.01 vs Untreated</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>WT    </td><td>IFNa</td><td>low </td><td>t6 </td><td>1532</td></tr>
	<tr><td>WT    </td><td>IFNa</td><td>high</td><td>t6 </td><td>2058</td></tr>
	<tr><td>P1104A</td><td>IFNa</td><td>low </td><td>t6 </td><td>2095</td></tr>
	<tr><td>P1104A</td><td>IFNa</td><td>high</td><td>t6 </td><td>2900</td></tr>
	<tr><td>WT    </td><td>IL10</td><td>low </td><td>t6 </td><td>1098</td></tr>
	<tr><td>WT    </td><td>IL10</td><td>high</td><td>t6 </td><td>   2</td></tr>
	<tr><td>P1104A</td><td>IL10</td><td>low </td><td>t6 </td><td>2082</td></tr>
	<tr><td>P1104A</td><td>IL10</td><td>high</td><td>t6 </td><td> 302</td></tr>
	<tr><td>WT    </td><td>IL23</td><td>low </td><td>t6 </td><td>2212</td></tr>
	<tr><td>WT    </td><td>IL23</td><td>high</td><td>t6 </td><td>2129</td></tr>
	<tr><td>P1104A</td><td>IL23</td><td>low </td><td>t6 </td><td>2499</td></tr>
	<tr><td>P1104A</td><td>IL23</td><td>high</td><td>t6 </td><td>2975</td></tr>
	<tr><td>WT    </td><td>IFNa</td><td>low </td><td>t24</td><td>  32</td></tr>
	<tr><td>WT    </td><td>IFNa</td><td>high</td><td>t24</td><td>  65</td></tr>
	<tr><td>P1104A</td><td>IFNa</td><td>low </td><td>t24</td><td>   8</td></tr>
	<tr><td>P1104A</td><td>IFNa</td><td>high</td><td>t24</td><td>  57</td></tr>
	<tr><td>WT    </td><td>IL10</td><td>low </td><td>t24</td><td>   0</td></tr>
	<tr><td>WT    </td><td>IL10</td><td>high</td><td>t24</td><td>  17</td></tr>
	<tr><td>P1104A</td><td>IL10</td><td>low </td><td>t24</td><td>   0</td></tr>
	<tr><td>P1104A</td><td>IL10</td><td>high</td><td>t24</td><td> 778</td></tr>
	<tr><td>WT    </td><td>IL23</td><td>low </td><td>t24</td><td> 603</td></tr>
	<tr><td>WT    </td><td>IL23</td><td>high</td><td>t24</td><td> 607</td></tr>
	<tr><td>P1104A</td><td>IL23</td><td>low </td><td>t24</td><td> 436</td></tr>
	<tr><td>P1104A</td><td>IL23</td><td>high</td><td>t24</td><td> 464</td></tr>
</tbody>
</table>




    
![png](IL23-RNAseq_files/IL23-RNAseq_20_1.png)
    


#### Differentially Expressed Genes <a name="part4"></a>


```R
vsd <- assay(vst(deresult, blind = FALSE))
vsd <- cbind(vsd, "row_var" = rowVars(vsd, na.rm = TRUE))

vsd_gene <- bind_cols("gene" = de_data$gene,
                      vsd) %>%
    pivot_longer(names_to = "sample", values_to = "norm_value", `1`:`56`) %>%
    inner_join(de_prop %>% mutate(sample_id = as.character(sample_id)),
               by = c("sample" = "sample_id")) %>%
    arrange(cytokine, group, dosage, time, rep)

var_genes <- vsd_gene %>%
    select(gene, row_var) %>%
    distinct() %>%
    arrange(-row_var) %>%
    head(500) %>%
    pull(gene)

vsd_gene_wide_sig <- vsd_gene %>%
    filter(gene %in% var_genes) %>%
    mutate(sample_name = str_c(cytokine, " ", dosage, " ", group, " ", time, " ", rep)) %>%
    select(gene, sample_name, norm_value) %>%
    pivot_wider(names_from = sample_name, values_from = norm_value)

vsd_gene_wide_all <- vsd_gene %>%
    mutate(sample_name = str_c(cytokine, " ", dosage, " ", group, " ", time, " ", rep)) %>%
    select(gene, sample_name, norm_value) %>%
    pivot_wider(names_from = sample_name, values_from = norm_value)
```


```R
options(repr.plot.width = 13.5, repr.plot.height = 12)
Heatmap(t(as.matrix(vsd_gene_wide_sig[,-1])),
        column_title = "Top 500 Most Variable Genes",
        cluster_columns = TRUE,
        cluster_rows = FALSE,
        row_split = gsub(" .*", "", names(vsd_gene_wide_sig[,-1])),
        use_raster = TRUE,
        raster_by_magick = TRUE,
        name = "Variance\nStabilized\nLog-Mean")
```


    
![png](IL23-RNAseq_files/IL23-RNAseq_23_0.png)
    



```R
de_genes_sig <- split_sumstats %>%
    filter(log2FoldChange > 1, padj < 0.01, !grepl("ENSG", gene)) %>%
    pull(gene) %>%
    unique()

heat_genes_wide <- vsd_gene_wide_all[vsd_gene_wide_all$gene %in% de_genes_sig,]
heat_data <- t(as.matrix(heat_genes_wide[,-1]))
colnames(heat_data) <- unlist(heat_genes_wide[,1])

options(repr.plot.width = 10, repr.plot.height = 15)
Heatmap(t(heat_data),show_row_names = FALSE,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        column_split = gsub("low |high |t6|t24| 1| 2| none", "", colnames(heat_genes_wide[,-1])),
        use_raster = TRUE,
        column_title_rot = 45,
        name = "Variance\nStabilized\nLog-Mean",
        row_title = "Genes DE vs None in any condition\nwith Log2FoldChange > 1 and FDR < 0.01")
```


    
![png](IL23-RNAseq_files/IL23-RNAseq_24_0.png)
    


#### Gene Visualizations <a name="part5"></a>


```R
plot_gene <- function(gene_id) {

    the_plot <- split_sumstats %>%
        filter(gene == gene_id) %>%
        mutate(dosage = relevel(as.factor(dosage), ref = "low")) %>%
        ggplot() +
            geom_pointrange(aes(x = dosage,
                                y = log2FoldChange,
                                ymin = log2FoldChange - 2*lfcSE,
                                ymax = log2FoldChange + 2*lfcSE,
                                color = time), position = position_dodge(width = 0.4)) +
            theme_pubr(base_size = 16,
                       x.text.angle = 45) +
            ggtitle(gene_id) +
            geom_hline(yintercept = 0) +
            facet_grid(cols = vars(cytokine, background))

    return(the_plot)
    
}


stat1 <- plot_gene("STAT1")
stat2 <- plot_gene("STAT2")
stat3 <- plot_gene("STAT3")

jak1 <- plot_gene("JAK1")
jak2 <- plot_gene("JAK2")
jak3 <- plot_gene("JAK3")
```


```R
options(repr.plot.width = 20, repr.plot.height = 6)
stat1 + stat2 + stat3
jak1 + jak2 + jak3
```


    
![png](IL23-RNAseq_files/IL23-RNAseq_27_0.png)
    



    
![png](IL23-RNAseq_files/IL23-RNAseq_27_1.png)
    



```R
options(repr.plot.width = 20, repr.plot.height = 6)
plot_gene("SOCS3") + 
    plot_gene("BCL3") +
    plot_gene("NFIL3")
```


    
![png](IL23-RNAseq_files/IL23-RNAseq_28_0.png)
    

