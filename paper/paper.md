---
title: "Octant-BMS TYK2 Manuscript Data Visuals"
author: "Nathan Abell and Conor Howard"
date: 'May 07, 2025'
output: github_document
---



# Common Data Processing



# Figure 1



```
## Error in file(filename, "r", encoding = encoding): cannot open the connection
```

### Main Heatmaps

![plot of chunk main-heatmap-activity](./fig-1/main-heatmap-activity-1.png)

![plot of chunk main-heatmap-stability](./fig-1/main-heatmap-stability-1.png)



# Figure S1























# I. SETUP

# Chunk options



### Packages


``` r
pacman::p_load(
  colorspace,
  ggbeeswarm,
  ggnewscale,
  ggh4x,
  ggpubr,
  ggrastr,
  ggsci,
  ggsignif,
  magrittr,
  readxl,
  paletteer,
  patchwork,
  scico,
  tidyverse
)
```

```
## Error in loadNamespace(x): there is no package called 'pacman'
```

### Variables


``` r
cbPalette <- c("#DC5E65", "#E69F00", "#56B4E9", "#0072B2", "darkgrey", "#009E73", "#F0E442", "pink", "#CC79A7", "lightgrey", "grey")

aa_order <- c(
  "*", "P", "G", "A", "M", "V", "L", "I", "T", "S",
  "C", "Q", "N", "Y", "W", "F", "E", "D", "K", "H", "R"
)

setwd("~/Analyses/bms-dms/paper")
```

```
## Error in setwd("~/Analyses/bms-dms/paper"): cannot change working directory
```

``` r
source("../../dms/src/model/dms-analysis-utils.R")
```

```
## Error in file(filename, "r", encoding = encoding): cannot open the connection
```

# Functions

## compute_difference


``` r
compute_difference <- function(test, control, sumstats) {
  df1 <- sumstats %>%
    dplyr::filter(condition == test) %>%
    select(-condition)
  df2 <- sumstats %>%
    dplyr::filter(condition == control) %>%
    select(-condition)

  df <- inner_join(df1, df2,
    by = c("pos", "chunk", "aa", "version")
  )

  new_stats <- df %>%
    mutate(
      log2FoldChange = log2FoldChange.x - log2FoldChange.y,
      log2StdError = sqrt(log2StdError.x^2 + log2StdError.y^2)
    ) %>%
    select(pos, chunk, aa, log2FoldChange, log2StdError, version) %>%
    ungroup()

  new_stats$condition <- paste0(test, " - ", control)

  return(new_stats)
}
```

## theme_pub


``` r
theme_pub <- function(base_size = 11, base_family = "") {
  # based on https://github.com/noamross/noamtools/blob/master/R/theme_nr.R
  # start with theme_bw and modify from there!
  theme_bw(base_size = base_size, base_family = base_family) + # %+replace%
    theme(
      # grid lines
      panel.grid.major.x = element_line(colour = "#ECECEC", linewidth = 0.5, linetype = 1),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(colour = "#ECECEC", linewidth = 0.5, linetype = 1),
      panel.background = element_blank(),

      # axis options
      # axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = rel(1), vjust = 0.25),
      axis.title.y = element_text(size = rel(1), vjust = 0.35),
      # axis.text      = element_text(color="black", size=rel(1)),
      axis.text.y = element_text(size = rel(1), color = "black"),
      axis.text.x = element_text(angle = 45, color = "black", hjust = 1, size = rel(1)),

      # legend options
      legend.title = element_text(size = rel(1)),
      legend.key = element_rect(fill = "white"),
      legend.key.size = unit(0.02, "npc"),
      legend.text = element_text(size = rel(1)),

      # facet options
      strip.text = element_text(size = rel(1)),
      strip.background = element_blank(),

      # title options
      plot.title = element_text(size = rel(1), hjust = 0.5),
      plot.subtitle = element_text(size = rel(0.8), hjust = 0),
    )
}
```

# Load and format data

After discussing with Conor, will use 1 and 10 U/mL IFNa data from run3, and 100 U/mL IFNa from run7


``` r
# load IFNa sumstats
ifna_sumstats_run3 <- read_tsv("sumstats/TYK2-run3-combined-cleaned.sumstats.tsv", show_col_types = FALSE) %>%
  # dplyr::filter(condition == "IFNalpha100_0") %>%
  mutate(
    aa = if_else(aa %in% c("X", "Stop", "*"), "*", aa),
    aa = factor(aa, levels = aa_order)
  )
```

```
## Error: 'sumstats/TYK2-run3-combined-cleaned.sumstats.tsv' does not exist in current working directory ('/bms-dms/paper').
```

``` r
# compute contrasts for IFNa data
conditions <- unique(ifna_sumstats_run3$condition)
```

```
## Error: object 'ifna_sumstats_run3' not found
```

``` r
noneIdx <- which(conditions == "None_0")
```

```
## Error: object 'conditions' not found
```

``` r
ifna_run3 <- map(
  conditions[-noneIdx],
  ~ compute_difference(
    test = .,
    control = "None_0",
    sumstats = ifna_sumstats_run3
  )
) %>%
  list_rbind() %>%
  mutate(condition = factor(condition,
    levels = c(
      "IFNbeta_100 - None_0",
      "IFNalpha_1 - None_0",
      "IFNalpha_10 - None_0",
      "IFNalpha_100 - None_0",
      "IFNalphaWithDrug_100 - None_0"
    )
  )) %>%
  dplyr::filter(condition %in% c("IFNalpha_1 - None_0", "IFNalpha_10 - None_0", "IFNalpha_100 - None_0")) %>%
  mutate(
    z_statistic = log2FoldChange / log2StdError,
    p.value = pmin(
      pnorm(z_statistic, mean = 0, sd = 1) * 2,
      (1 - pnorm(z_statistic, sd = 1)) * 2
    )
  ) %>%
  mutate(run = "run3")
```

```
## Error: object 'conditions' not found
```

``` r
# load run7 IFN100 data
ifna_run7 <- read_tsv("sumstats/TYK2-run7-combined-cleaned.sumstats.tsv", show_col_types = FALSE) %>% 
  dplyr::filter(condition == "IFNalpha100_0") %>%
  mutate(condition = "IFNalpha_100 - None_0") %>% 
  mutate(
    aa = if_else(aa %in% c("X", "Stop", "*"), "*", aa),
    aa = factor(aa, levels = aa_order)
  ) %>%
  mutate(run = "run7") %>%
  rename(z_statistic = statistic)
```

```
## Error: 'sumstats/TYK2-run7-combined-cleaned.sumstats.tsv' does not exist in current working directory ('/bms-dms/paper').
```

``` r
# merge experiments
ifna_merge <- bind_rows(ifna_run3, ifna_run7) %>%
  mutate(run = factor(run, levels = c("run7", "run3"))) %>% 
  arrange(run) %>% 
  # for IFNalpha100_0, keep sumstats from run7 only where possible, default to run3 when not found in run7
  group_by(pos, aa, condition) %>% 
  slice_head(n = 1) %>%
  ungroup() %>% 
  pivot_longer(c("log2FoldChange"), names_to = "param", values_to = "param_estimate") %>%
  select(pos, aa, chunk, param, run, param_estimate, param_se = log2StdError, condition, z_statistic, p.value) %>%
  mutate(
    assay = "IFNa",
    p.adj = p.adjust(p.value, method = "BH")
  )
```

```
## Error: object 'ifna_run3' not found
```

``` r
stability <- read_tsv("sumstats/TYK2-FLOW-flow-cleaned.midpoints.tsv", show_col_types = FALSE) %>%
  mutate(
    aa = if_else(aa %in% c("X", "Stop", "*"), "*", aa),
    aa = factor(aa, levels = aa_order)
  ) %>%
  pivot_longer(c("midpoint_shift"), names_to = "param", values_to = "param_estimate") %>%
  select(pos, aa, chunk, param, param_estimate, param_se = midpoint_shift_se, z_statistic = statistic, p.value, p.adj) %>%
  mutate(assay = "Flow")
```

```
## Error: 'sumstats/TYK2-FLOW-flow-cleaned.midpoints.tsv' does not exist in current working directory ('/bms-dms/paper').
```

``` r
merge_sumstats <- bind_rows(ifna_merge, stability) %>%
  # add annotation labels for plotting
  mutate(facet_label = case_when(
    assay == "Flow" ~ "Abundance",
    assay == "IFNa" ~ paste0(
      "IFN\u03B1 ", str_extract(condition,
        pattern = "IFNalpha_?(\\d+)",
        group = 1
      ),
      " U/mL"
    )
  )) %>%
  mutate(facet_label = fct_relevel(facet_label, "Abundance", after = Inf))
```

```
## Error: object 'ifna_merge' not found
```

``` r
#annotate data with 3 letter aa
aalet <- read_tsv("data/aa-letters.tsv", col_names = c("full", "thr", "aa")) %>%
    select(-full) 

wtseq <- read_tsv("data/TYK2-sequence.tsv", col_names = c("pos", "thr"), skip = 1) %>%
    left_join(aalet, by = "thr") %>%
    rename("wt_aa" = "aa",
           "wt_thr" = "thr")
  
merge_sumstats <- merge_sumstats %>%
  left_join(wtseq, by = "pos") %>% 
  left_join(aalet %>% rename("mut_thr" = "thr"), by = c("aa")) %>% 
  mutate(aa = factor(aa, levels = aa_order)) %>%
  replace_na(replace = list(mut_thr = "Stop")) %>% 
  mutate(mut = paste0(wt_thr, pos, mut_thr)) %>% 
  #remove amino acid 1188
  dplyr::filter(pos != 1188)
```

```
## Error: object 'merge_sumstats' not found
```

# Load metadata


``` r
#load variant annotations / scores
tyk2_clinvar <- read_xlsx("data/clinvar_TYK2_20250109.xlsx")  %>% 
    select(Name,
           gene = `Gene(s)`,
           clinvar_class = `Germline classification`,
           mut = `Protein Change (header)`) %>% 
  mutate(clinvar_class = gsub("\\(.*","",clinvar_class)) %>%
  #remove Gly512Arg variant with Uncertain significance
  dplyr::filter(!(mut == "Gly512Arg" & clinvar_class == "Uncertain significance")) %>% 
  #aggregate benign and likely benign into single category
  mutate(clinvar_class = case_when(clinvar_class %in% c("Benign", "Benign/Likely benign", "Likely benign") ~ "Benign/Likely benign",
                                   .default = clinvar_class)) %>% 
  distinct()
```

```
## Error in read_xlsx("data/clinvar_TYK2_20250109.xlsx"): could not find function "read_xlsx"
```

``` r
eve_scores_b05 <- read_csv("data/TYK2_HUMAN_b05_20000_samples.csv") %>%
    mutate(mutations = gsub("^.","",mutations)) %>%
    select(-protein_name) %>%
    dplyr::filter(mutations != "t") %>%
    separate(mutations, c("pos", "aa"), sep = -1) %>%
    rename("EVE Score 1" = "evol_indices") %>%
    mutate(pos = as.numeric(pos))

eve_scores_b03 <- read_csv("data/TYK2_HUMAN_b03_20000_samples.csv") %>%
    mutate(mutations = gsub("^.","",mutations)) %>%
    select(-protein_name) %>%
    dplyr::filter(mutations != "t") %>%
    separate(mutations, c("pos", "aa"), sep = -1) %>%
    rename("EVE Score 2" = "evol_indices") %>%
    mutate(pos = as.numeric(pos))

eve_scores <- inner_join(eve_scores_b05, eve_scores_b03, by = c("pos", "aa"))

tyk2_sift <- read_delim("data/TYK2-sift.tsv.gz", delim = " ") %>%
    pivot_longer(names_to = "aa", values_to = "sift", A:Y) %>%
    mutate(mut_sift = paste0(wt, aa)) %>%
    select(mut_sift, sift)

tyk2_polyp <- read_tsv("data/TYK2-polyphen2.tsv.gz") %>%
    select(o_pos, o_aa2, pph2_prob) %>%
    rename("PolyPhen2" = "pph2_prob")

alphamissense <- read_tsv("data/TYK2.Uniprot.P29597.AlphaMissense_aa.substitutions.tsv",
                          col_names = c("uniprot_id" ,"mut_id", "alphamissense_score", "AlphaMissense Class")) 

#annotate sumstats
anno_sumstats <- merge_sumstats %>%
  left_join(tyk2_clinvar, by = "mut", relationship = "many-to-many") %>% 
  replace_na(replace = list(clinvar_class = "not in ClinVar")) %>% 
  mutate(mut_sift = paste0(pos, wt_aa, aa),
         mut_alphaM = paste0(wt_aa, pos, aa)) %>%
  left_join(tyk2_sift, by = "mut_sift") %>% 
  left_join(tyk2_polyp, by = c("pos" = "o_pos",
                                 "aa" = "o_aa2")) %>%
  left_join(alphamissense, by = c("mut_alphaM" = "mut_id")) %>% 
  select(!c(mut_sift, mut_alphaM))
```

```
## Error: object 'merge_sumstats' not found
```

``` r
anno_df <- anno_sumstats %>% 
  dplyr::filter(mut_thr != "Stop") %>% 
  dplyr::filter(assay == "Flow" | run == "run7") %>%
  pivot_longer(c("clinvar_class", "AlphaMissense Class"),
               names_pattern = "(clinvar|AlphaMissense).*",
               names_to = "classifier",
               values_to = "classifier_value") %>%
  mutate(classifier_value = str_to_title(classifier_value)) %>% 
  mutate(classifier_value = factor(classifier_value, levels = c("Pathogenic",
                                                               "Conflicting Classifications Of Pathogenicity",
                                                               "Uncertain Significance",
                                                               "Ambiguous",
                                                               "Likely Benign",
                                                               "Benign/Likely Benign",
                                                               "Benign",
                                                               "Not In Clinvar")))
```

```
## Error: object 'anno_sumstats' not found
```

# Figures

## Fig S1x1 - Flow + IFNa-1,10,100 heatmaps


``` r
S1x1 <- merge_sumstats %>%
  mutate(facet_label = str_replace_all(facet_label, "([0-1]+)", "\n\\1")) %>% 
  mutate(facet_label = fct_relevel(facet_label, "Abundance", after = Inf)) %>% 
  mutate(aa = fct_recode(aa, "Stop" = "*")) %>% 
  ggplot() +
  facet_grid2(
    rows = "facet_label",
    drop = FALSE,
    scales = "fixed",
    axes = "all"
  ) + 
  geom_tile(aes(x = pos, y = aa, fill = z_statistic)) +
  scale_fill_scico_mid("Z statistic", limits = c(-10, 2.5),
                       palette = "vik", mid = 0,
                       oob = scales::squish) +
  scale_x_continuous(
    breaks = seq(0, 1000, by = 250),
    expand = c(0, 0)
  ) +
  ylab("") +
  xlab("TYK2 Amino Acid Position") +
  theme_pub(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5,
                               size = unit(8, "points")),
    legend.key.width = unit(0.05, "npc"),
    axis.text.y = element_text(size = unit(8, "points")),
    legend.position = "top",
    legend.text = element_text(size = unit(8, "points")),
    legend.title = element_text(margin = margin(0,10,0,0)),
    legend.box.spacing = unit(0.001, "npc"),
    strip.text.y = element_text(hjust = 0.5, angle = 0, size = unit(8, "points")),
  ) 
```

```
## Error: object 'merge_sumstats' not found
```

``` r
S1x1
```

```
## Error: object 'S1x1' not found
```

## Fig S1x2

### Fig S1x2-A - Barplots with LoF/GoF counts

I'm going to just use a 1% FDR to determine GoF and LoF variants with no minimum log2fold or midpoint shift threshold


``` r
A <- merge_sumstats %>%
  # dplyr::filter(assay == "Flow" |) %>%
  mutate(facet_label = str_replace_all(facet_label, "([0-1]+)", "\n\\1")) %>% 
  mutate(facet_label = fct_relevel(facet_label, "Abundance", after = Inf)) %>%
  mutate(status = case_when(param_estimate < 0 & p.adj < 0.01 ~ "LoF",
    param_estimate > 0 & p.adj < 0.01 ~ "GoF",
    .default = "no effect"
  )) %>%
  mutate(status = factor(status, levels = c("GoF", "no effect", "LoF"))) %>%
  ggplot(aes(x = facet_label, fill = status)) +
  geom_bar(position = "fill") +
  geom_label(aes(label = after_stat(count)),
    stat = "count",
    position = position_fill(vjust = 0.5),
    show.legend = FALSE,
    size = rel(2),
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual("", values = c("#DC5E65", "lightgrey", "#0072B2")) +
  theme_pub() +
  labs(
    
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "top",
        panel.grid.major.x = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        legend.margin=margin(0,0,-25,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.text = element_text(size = rel(0.5)),
        legend.key.size = unit(0.01, "npc"),
        axis.text= element_text(size = unit(6, "points")),
        plot.tag.position = c(0.01, 1.01)
        )
```

```
## Error: object 'merge_sumstats' not found
```

### Fig S1x2-B - stop effect density plots


``` r
B <- merge_sumstats %>%
  # dplyr::filter(grepl("- None0", condition)) %>%
  mutate(group = factor(if_else(aa == "*", "stop variants", "other variants"),
                        levels = c("stop variants", "other variants"))) %>%
  ggplot() +
  geom_density(
    aes(
      x = z_statistic, group = group,
      fill = group
    ),
    color = NA,
    alpha = 0.5
  ) +
  geom_text(aes(label = facet_label), x = -15, y = 0.4, size = rel(2), hjust = 0) +
  facet_grid2(facet_label ~ ., scales = "fixed") +
  scale_fill_manual("", values = c(
    "stop variants" = "red",
    "other variants" = "black")) +
  labs(
    x = "Z-Statistic",
    y = "Density",
    # title = "Fig S1D",
    # subtitle = "IFN\u03B1 data normalized 0 cytokine condition"
  ) +
  coord_cartesian(xlim = c(-15, 7)) +
  scale_y_continuous(breaks = 0.5) +
  # scale_x_continuous(expand = c(0, 0)) +
  theme_pub(base_size = 12) +
  theme(
    axis.title = element_text(size = rel(0.6)),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.box = element_blank(),
    legend.position = "top",
    legend.key = element_rect(color = NA),
    legend.margin = margin(0,0,-25,0),
    legend.box.spacing = unit(0.001, "npc"),
    legend.title = element_blank(),
    legend.text = element_text(size = rel(0.5)),
    legend.key.height = unit(0.01, "npc"),
    strip.text.y = element_blank(),# element_text(angle = 0, size = rel(0.7)),
    axis.text= element_text(size = unit(6, "points")),
    plot.tag.position = c(0.01, 0.9)
  )
```

```
## Error: object 'merge_sumstats' not found
```

### Fig S1x2-C - AlphaMissense


``` r
am_df <-  anno_df %>% 
  dplyr::filter(classifier == "AlphaMissense")
```

```
## Error: object 'anno_df' not found
```

``` r
C <- am_df %>%
    ggplot(aes(x = classifier_value, y = z_statistic,
             color = classifier_value)) +
  geom_quasirandom_rast(size = 0.5,
                        show.legend = FALSE) +
  # geom_violin(
  #   aes(fill = classifier_value),
  #   scale = "width",
  #   show.legend = FALSE
  #   ) +
  geom_boxplot(show.legend = FALSE,
               col = "black",
               outliers = FALSE,
               width = 0.09) +
  scale_color_manual(values = c(
    "Pathogenic" = "#DC5E65",
      "Ambiguous" =  "#E69F00",
      "Benign" =  "#0072B2")
  ) +
    scale_fill_manual(values = c(
    "Pathogenic" = "#DC5E65",
      "Ambiguous" =  "#E69F00",
      "Benign" =  "#0072B2")
  ) +
  # add significance calculation using KS test
  geom_signif(
    color = "black",
    comparisons = list(
      c("Pathogenic", "Ambiguous"),
      c("Pathogenic", "Benign"),
      c("Ambiguous", "Benign")
    ),
    map_signif_level = function(p) {
      if(p < 0.001) return("***")
      else if(p < 0.01) return("**")
      else if(p < 0.05) return("*")
      else return("ns")
    },
    test = "ks.test",
    margin_top = 0.01,
    step_increase = 0.025,
    tip_length = 0.005,
    vjust = 0.65,
    size = 0.5,
    textsize = 3
  ) +
  facet_grid2(. ~ facet_label,
              scales = "free_y",
              space = "free") +
  labs(
    x = "AlphaMissense annotation",
    y = "Z statistic",
  ) +
  ylim(-15,7) +
  theme_pub(base_size = 11) +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(size = rel(0.8)),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = rel(0.6)),
    strip.text.x = element_text(size = rel(0.7)),
    plot.tag.position = c(0.01, 0.98)
    )
```

```
## Error: object 'am_df' not found
```

### Fig S1x2-D - Fiducial variant forest plot


``` r
spikeins <- read_tsv("../paper/data/tyk2-spikeins.tsv",
  col_names = c("pos", "mut_id", "class")
) %>%
  mutate(
    class = if_else(class == "-", "Unknown", class),
    class = factor(class, levels = c("LOF", "Unknown", "GOF"))
  ) %>%
  arrange(class, pos) %>%
  separate_wider_regex(mut_id,
    patterns = c(".*\\d+", "aa" = ".*"),
    cols_remove = FALSE
  ) %>%
  mutate(aa = ifelse(aa == "Stop", "*", aa)) 


spikeins$mut_id <- factor(spikeins$mut_id, levels = spikeins$mut_id)

spikeins_join <- spikeins %>%
  inner_join(merge_sumstats, by = c("pos", "aa")) %>%
  arrange(class, pos) %>%
  mutate(class = factor(class, levels = c("GOF", "Unknown", "LOF"))) %>%
  mutate(facet_label = str_replace_all(facet_label, "([0-1]+)", "\n\\1")) %>%
  mutate(facet_label = fct_relevel(facet_label, "Abundance", after = Inf)) %>% 
  mutate(significance = factor(if_else(p.adj < 0.01, "FDR < 1%", "not significant")))
```

```
## Error: object 'merge_sumstats' not found
```

``` r
D1 <- spikeins_join %>% 
  filter(facet_label != "Abundance") %>% 
  ggplot() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 26.5,
             color = "#DC5E65",
             linetype = "dashed") +
  geom_text(
    data = tibble(facet_label = "IFNα \n1 U/mL",
                  label = "GOF"),
    aes(label = label),
    color = "#DC5E65",
    y = 28.5,
    x = -1.5
  ) +
  geom_hline(yintercept = 15.5,
             color = "#0072B2",
             linetype = "dashed") +
    geom_text(
    data = tibble(facet_label = "IFNα \n1 U/mL",
                  label = "LOF"),
    aes(label = label),
    color = "#0072B2",
    y = 12.5,
    x = -1.5
  ) +
  geom_pointrange(
    aes(
      x = param_estimate,
      y = factor(mut_id),
      xmin = param_estimate - 2 * param_se,
      xmax = param_estimate + 2 * param_se,
      color = significance
    ),
    linewidth = 0.25,
    size = 0.1
  ) +
  theme_pubr(base_size = 13) +
  labs(
    x = expression(Log[2]~fold~change~"\u00B1"~2~standard~errors),
    y = "") +
    #," ± 2 standard errors"),
  #  subtitle = "IFNa data normalized to 0 cytokine control\nEffect size = log2FoldChange for IFNa data\nEffect size = midpoint shift for Abundance DMS (see methods)\nPoints in black (padj < 0.01)"
  scale_color_manual("", values = c("black", "lightgrey")) +
  xlim(c(-2.5, 1)) +
  theme_pub(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    strip.text.x = element_text(size = unit(6, "points")),
    axis.title = element_text(size = rel(0.6)),
    legend.position = "top",
    legend.margin=margin(0,0,-10,0),
    axis.text.y = element_text(size = unit(6, "points")),
    axis.text.x = element_text(angle = 0, hjust = 0.6,
                               size = unit(6, "points")),
    panel.grid.major.y = element_blank(),
    plot.tag.position = c(0.01, 0.98)
  ) +
  facet_grid2(. ~ facet_label, scales = "fixed")
```

```
## Error: object 'spikeins_join' not found
```

``` r
#panel for abundance DMS
D2 <- spikeins_join %>% 
  filter(facet_label == "Abundance") %>% 
  ggplot() +
  geom_vline(xintercept = 0) +
   geom_hline(yintercept = 26.5,
             color = "#DC5E65",
             linetype = "dashed") +
  geom_hline(yintercept = 15.5,
             color = "#0072B2",
             linetype = "dashed") +
  geom_pointrange(
    aes(
      x = param_estimate,
      y = mut_id,
      xmin = param_estimate - 2 * param_se,
      xmax = param_estimate + 2 * param_se,
      color = significance
    ),
    linewidth = 0.25,
    size = 0.1,
    show.legend = FALSE
  ) +
  theme_pubr(base_size = 13) +
  labs(
    x = "Midpoint shift\n± 2 standard errors",
    y = "",
    # title = "Fig S1E",
  #  subtitle = "IFNa data normalized to 0 cytokine control\nEffect size = log2FoldChange for IFNa data\nEffect size = midpoint shift for Abundance DMS (see methods)\nPoints in black (padj < 0.01)"
  ) +
  scale_color_manual("", values = c("black", "lightgrey")) +
  scale_x_continuous(limits = c(-1,1),
                     breaks = c(-1,0,1)) +
  theme_pub(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    strip.text.x = element_text(size = unit(6, "points")),
    axis.title = element_text(size = rel(0.6)),
    axis.text.y = element_blank(),
    plot.margin = margin(0,0,0,10),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = unit(4, "points")),
    axis.text.x = element_text(angle = 0, hjust = 0.6,
                               size = unit(6, "points")),
    panel.grid.major.y = element_blank(),
    plot.tag.position = c(0.01, 0.98)
  ) +
  facet_grid2(. ~ facet_label)
```

```
## Error: object 'spikeins_join' not found
```

``` r
D <- D1 + (D2 + plot_layout(tag_level = "new")) +
  plot_layout(widths = c(3,1))
```

```
## Error: object 'D1' not found
```

### Fig S1x2-E - ClinVar annotation


``` r
anno_df <- anno_df %>% 
  mutate(classifier_value = fct_recode(classifier_value, "Conflicting\nClassifications\nOf Pathogenicity" = "Conflicting Classifications Of Pathogenicity",
                                       "Uncertain\nSignificance"  = "Uncertain Significance",
                                       "Benign/\nLikely Benign" = "Benign/Likely Benign"))
```

```
## Error: object 'anno_df' not found
```

``` r
counts <- anno_df %>% 
  dplyr::filter(classifier == "clinvar") %>% 
  select(mut, classifier_value) %>% 
  distinct() %>% 
  group_by(classifier_value) %>% 
  summarise(count = n())
```

```
## Error: object 'anno_df' not found
```

``` r
E <- anno_df %>% 
  dplyr::filter(classifier == "clinvar") %>% 
  left_join(counts) %>% 
  arrange(classifier_value) %>%
  # mutate(classifier_value_count = paste("count", count)) %>% 
  mutate(classifier_value_count =  fct_inorder(paste0(classifier_value, " (", format(count, big.mark = ",", trim = TRUE),")"))) %>%
  ggplot(aes(
    x = z_statistic,
    y = classifier_value_count,
    color = classifier_value_count,
             fill = classifier_value_count)) +
  geom_vline(xintercept = 0, color = "lightgrey", width = 0.5) +
  geom_quasirandom_rast(orientation = "y",
                   size = 0.5,
                   width = 0.2,
                   varwidth = TRUE,
                   show.legend = FALSE) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  facet_grid2(. ~ facet_label,
              scales = "free_x",
              space = "free") +
  labs(
    x = "Z statistic",
    y = "Clinvar annotation",
  ) +
  theme_pub(base_size = 11) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.6),
        axis.text.y = element_text(angle = 0, hjust = 0.5),
        axis.text= element_text(size = unit(6, "points")),
        strip.text.y = element_text(angle = 0, size = unit(8, "points")),
        axis.title = element_text(size = unit(10, "points")),
        plot.tag.position = c(0.01, 0.98)
        )
```

```
## Error: object 'anno_df' not found
```

## Merged Figure S1x2


``` r
layout <- c(
  area(1, 1, 1, 1),
  area(1, 2, 1, 2),
  area(1, 3, 1, 3),
  area(2, 1, 2, 2),
  area(2, 3, 2, 3)
)
A + free(B, type = "label") + free(C, side = "l") + free(D) + E +
  plot_layout(heights = c(1, 1),
              guides = "keep",
              design = layout) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.key = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title = element_text(size = unit(8, "points")),
    axis.text = element_text(size = unit(6, "points")),
    plot.tag = element_text(size = unit(14, "points")),
    # plot.title = element_text(size = unit(14, "points")),
    legend.text = element_text(size = unit(8, "points")),
    # strip.text.x = element_blank(),
    strip.text.x = element_text(size = unit(6, "points")),
    plot.tag.position = c(0.01, 1.01)
  )
```

```
## Error: object 'A' not found
```

# Figure 2



``` r
knitr::opts_chunk$set(
  fig.path = "./fig-2/"
)

library(tidyverse)

ifna_signaling <- read_tsv("../sumstats/TYK2-run3.sumstats.tsv", show_col_types = FALSE) %>%
    filter(condition == "IFNalpha100_0") %>%
    mutate(p.adj_dms = p.adjust(p.value, method = "BH")) %>%
    rename("DMS Statistic" = "statistic")
stability <- read_tsv("../sumstats/TYK2-vamp2.midpoints.tsv", show_col_types = FALSE)  %>%
    rename("p.adj_flow" = "p.adj",
           "FlowDMS Statistic" = "statistic")

signal_stab <- inner_join(ifna_signaling,
                          stability,
                          by = c("chunk", "pos", "aa")) %>%
  mutate(class = if_else(log2FoldChange < 0 &
                         p.adj_dms < 0.05 &
                         (p.adj_flow > 0.1 | midpoint_shift > 0), "signaling-only LOF", "other"))
```

### IFN-alpha Signaling vs Stability


``` r
ggplot(signal_stab) +
  geom_point(aes(x = log2FoldChange,
                 y = midpoint_shift,
                 color = class,
                 alpha = if_else(class == "signaling-only LOF", 1, 0.3)), stroke = 0) +
  xlab("DMS Mutant vs WT Log2 Fold Change") +
  ylab("FlowDMS Midpoint Shift") +
  scale_color_manual(values = c("signaling-only LOF" = "lightblue", "other" = "gray")) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_alpha_identity()
```

![plot of chunk signal-vs-stability](./fig-2/signal-vs-stability-1.png)


``` r
pos_counts <- signal_stab %>%
  group_by(pos) %>%
  count(class) %>%
  filter(class == "signaling-only LOF")

ggplot(pos_counts) +
  geom_bar(aes(x = pos, y = n), stat = "identity") +
  coord_flip() +
  scale_x_reverse() +
  xlab("TYK2 Position") +
  ylab("Number of\nSignaling-Only\nLOF Variants")
```

![plot of chunk signal-only-hist](./fig-2/signal-only-hist-1.png)


# Figure 3



``` r
knitr::opts_chunk$set(
  fig.path = "./fig-3/"
)

colors <- c("Neither" = "gray",
            "Both" = "purple",
            "BMS-986202" = "red",
            "Zasocitinib" = "blue")

aa_order <- c("G","A","P","K","H",
              "R","T","S","C","E",
              "D","M","V","I","L",
              "Y","W","F","Q","N")
```

### Drug Resistance


``` r
bms_resist_all <- contrast_sumstats %>%
    filter(grepl("BMS", condition)) %>%
    mutate(assay = factor(assay, levels = c("assay1", "assay2", "assay3")))

ggplot(bms_resist_all) +
    geom_hline(yintercept = -log10(0.01), color = "red") +
    geom_point(aes(x = log2Contrast,
                   y = -log10(fdr_contrast),
                   color = if_else(-log10(fdr_contrast) > -log10(0.01), "red", "black")),
               size = 2) +
    xlab("Log2 Fold Change") + ylab("-log10(FDR)") +
    scale_color_identity() +
    ggtitle("IFN-alpha + BMS-986202, Normalized to Untreated") +
    facet_wrap(~assay + condition, nrow = 1)
```

![plot of chunk bms-drig-resist-1](./fig-3/bms-drig-resist-1-1.png)


``` r
bms_resist <- contrast_sumstats %>%
    filter(assay == "assay1",
           condition == "IFNalpha100+BMS-986202_100")

ggplot(bms_resist) +
    geom_hline(yintercept = -log10(0.01), color = "red") +
    geom_point(aes(x = log2Contrast,
                   y = -log10(fdr_contrast),
                   color = if_else(fdr_contrast < 0.01, "red", "black")),
               size = 2) +
    geom_label_repel(data = bms_resist %>% filter(fdr_contrast < 0.01),
                     aes(x = log2Contrast,
                         y = -log10(fdr_contrast),
                         color = if_else(fdr_contrast < 0.01, "red", "black"),
                         label = str_c(pos, aa)),
                     size = 2, nudge_x = 1) +
    xlab("Log2 Fold Change") + ylab("-log10(FDR)") +
    scale_color_identity() +
    ggtitle("IFN-alpha + BMS-986202, Normalized to Untreated")
```

![plot of chunk bms-drig-resist-2](./fig-3/bms-drig-resist-2-1.png)

### Binding Site (GoF) Comparisons


``` r
pos_counts <- contrast_sumstats_bms %>%
  filter(fdr_contrast < 0.01, log2Contrast > 0) %>%
  count(pos, condition, assay)
```

```
## Error: object 'contrast_sumstats_bms' not found
```

``` r
ggplot(pos_counts) +
  geom_bar(aes(x = pos, y = n), stat = "identity") +
  coord_cartesian(xlim = c(1, 1147)) +
  facet_wrap(~condition+assay, ncol = 1) +
  xlab("TYK2 Position") + ylab("Number of Significant Positive Variants")
```

```
## Error in `combine_vars()`:
## ! At least one layer must contain all faceting variables: `condition` and `assay`
## x Plot is missing `c("condition", "assay")`
## x Layer 1 is missing `c("condition", "assay")`
```

### Chemical Footprints

Group assignment for Zasocitinib and BMS-986202 drug resistance profiles:

``` r
sumstats_resist <- contrast_sumstats %>%
    filter(assay == "assay4", condition %in% c("IFNalpha+Zasocitinib_1e-06", "IFNalpha+BMS-986202_2e-08")) %>%
    select(pos, condition, aa, statistic, fdr_contrast) %>%
    pivot_wider(names_from = condition, values_from = c(statistic, fdr_contrast)) %>%
    mutate(`FDR < 0.01` = case_when(`fdr_contrast_IFNalpha+BMS-986202_2e-08` < 0.01 & `fdr_contrast_IFNalpha+Zasocitinib_1e-06` < 0.01 ~ "Both",
                           `fdr_contrast_IFNalpha+BMS-986202_2e-08` < 0.01 & `fdr_contrast_IFNalpha+Zasocitinib_1e-06` > 0.01 ~ "BMS-986202",
                           `fdr_contrast_IFNalpha+BMS-986202_2e-08` > 0.01 & `fdr_contrast_IFNalpha+Zasocitinib_1e-06` < 0.01 ~ "Zasocitinib",
                           TRUE ~ "Neither"),
          `FDR < 0.01` = if_else(`statistic_IFNalpha+Zasocitinib_1e-06` < 0 & `statistic_IFNalpha+BMS-986202_2e-08` < 0,
                                 "Neither",
                                 `FDR < 0.01`),
          aa = factor(aa, levels = aa_order))
```

```
## Error in `mutate()`:
## i In argument: `FDR < 0.01 = case_when(...)`.
## Caused by error in `case_when()`:
## ! Failed to evaluate the left-hand side of formula 1.
## Caused by error:
## ! object 'fdr_contrast_IFNalpha+BMS-986202_2e-08' not found
```

``` r
sig_vars_resist <- sumstats_resist %>% filter(`FDR < 0.01` != "Neither") %>% pull(pos) %>% unique()
```

```
## Error: object 'sumstats_resist' not found
```

``` r
prepped_data_resist <- sumstats_resist %>% filter(pos %in% sig_vars_resist)
```

```
## Error: object 'sumstats_resist' not found
```

``` r
names(prepped_data_resist) <- c("pos", "aa", "Z_BMS", "FDR_BMS")
```

```
## Error: object 'prepped_data_resist' not found
```

Group assignment for Zasocitinib and BMS-986202 drug potentiation profiles:

``` r
sumstats_potentiate <- contrast_sumstats %>%
    filter(assay == "assay7", condition %in% c("IFNalpha100+Zasocitinib_7e-09", "IFNalpha100+BMS-986202_2e-08")) %>%
    select(pos, condition, aa, statistic, fdr_contrast) %>%
    pivot_wider(names_from = condition, values_from = c(statistic, fdr_contrast)) %>%
    mutate(`FDR < 0.01` = case_when(`fdr_contrast_IFNalpha100+BMS-986202_2e-08` < 0.01 & `fdr_contrast_IFNalpha100+Zasocitinib_7e-09` < 0.01 ~ "Both",
                           `fdr_contrast_IFNalpha100+BMS-986202_2e-08` < 0.01 & `fdr_contrast_IFNalpha100+Zasocitinib_7e-09` > 0.01 ~ "BMS-986202",
                           `fdr_contrast_IFNalpha100+BMS-986202_2e-08` > 0.01 & `fdr_contrast_IFNalpha100+Zasocitinib_7e-09` < 0.01 ~ "Zasocitinib",
                           TRUE ~ "Neither"),
          `FDR < 0.01` = if_else(`fdr_contrast_IFNalpha100+Zasocitinib_7e-09` < 0 & `fdr_contrast_IFNalpha100+BMS-986202_2e-08` < 0,
                                 "Neither",
                                 `FDR < 0.01`),
          aa = factor(aa, levels = aa_order))
```

```
## Error in `mutate()`:
## i In argument: `FDR < 0.01 = case_when(...)`.
## Caused by error in `case_when()`:
## ! Failed to evaluate the left-hand side of formula 1.
## Caused by error:
## ! object 'fdr_contrast_IFNalpha100+BMS-986202_2e-08' not found
```

``` r
sig_vars_potentiate <- sumstats_potentiate %>% filter(`FDR < 0.01` != "Neither") %>% pull(pos) %>% unique()
```

```
## Error: object 'sumstats_potentiate' not found
```

``` r
prepped_data_potentiate <- sumstats_potentiate %>% filter(pos %in% sig_vars_potentiate)
```

```
## Error: object 'sumstats_potentiate' not found
```
