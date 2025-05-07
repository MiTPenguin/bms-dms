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



# Figure 2




### IFN-alpha Signaling vs Stability

![plot of chunk signal-vs-stability](./fig-2/signal-vs-stability-1.png)

![plot of chunk signal-only-hist](./fig-2/signal-only-hist-1.png)


# Figure 3




### Drug Resistance

![plot of chunk bms-drig-resist-1](./fig-3/bms-drig-resist-1-1.png)

![plot of chunk bms-drig-resist-2](./fig-3/bms-drig-resist-2-1.png)

### Binding Site (GoF) Comparisons


```
## Error: object 'contrast_sumstats_bms' not found
```

```
## Error in `combine_vars()`:
## ! At least one layer must contain all faceting variables: `condition` and `assay`
## x Plot is missing `c("condition", "assay")`
## x Layer 1 is missing `c("condition", "assay")`
```

### Chemical Footprints

Group assignment for Zasocitinib and BMS-986202 drug resistance profiles:

```
## Error in `mutate()`:
## i In argument: `FDR < 0.01 = case_when(...)`.
## Caused by error in `case_when()`:
## ! Failed to evaluate the left-hand side of formula 1.
## Caused by error:
## ! object 'fdr_contrast_IFNalpha+BMS-986202_2e-08' not found
```

```
## Error: object 'sumstats_resist' not found
```

```
## Error: object 'sumstats_resist' not found
```

```
## Error: object 'prepped_data_resist' not found
```

Group assignment for Zasocitinib and BMS-986202 drug potentiation profiles:

```
## Error in `mutate()`:
## i In argument: `FDR < 0.01 = case_when(...)`.
## Caused by error in `case_when()`:
## ! Failed to evaluate the left-hand side of formula 1.
## Caused by error:
## ! object 'fdr_contrast_IFNalpha100+BMS-986202_2e-08' not found
```

```
## Error: object 'sumstats_potentiate' not found
```

```
## Error: object 'sumstats_potentiate' not found
```
