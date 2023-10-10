## Differential Expression by TYK2 Genotype

- [Data Structure](#part1)
- [Bioinformatics Processing](#part2)
- [Differential Expression](#part3)
- [Primary Results: Cytokine Effects Within Each Genotype](#part4)
  - [Luciferase](#part4a)
  - [Pathway Components](#part4b)
  - [Volcano Plots](#part4c)
- [Secondary Results: Comparing Cytokine Effects Between Genotypes](#part5)
  - [Heatmaps](#part5a)

### Data Structure <a name="part1"></a>

We want to describe how gene expression changes based on the genotype of TYK2 and based on cytokine treatment, while taking into account doxycycline concentration. Accordingly, we have RNA-Seq libraries with the following sample structure:

| Genotype | Doxycycline (nM) | Cytokine (U/mL) | Replicate |
| :-- | :-- | :-- | :-- |
| WT | 100 | 0 | 1 |
| WT | 100 | 0 | 2 |
| EFF | 100 | 0 | 1 |
| EFF | 100 | 0 | 2 |
| P1104A | 100 | 0 | 1 |
| P1104A | 100 | 0 | 2 |
| WT | 100 | 100 | 1 |
| WT | 100 | 100 | 2 |
| EFF | 100 | 100 | 1 |
| EFF | 100 | 100 | 2 |
| P1104A | 100 | 100 | 1 |
| P1104A | 100 | 100 | 2 |
| WT | 1000 | 0 | 1 |
| WT | 1000 | 0 | 2 |
| EFF | 1000 | 0 | 1 |
| EFF | 1000 | 0 | 2 |
| P1104A | 1000 | 0 | 1 |
| P1104A | 1000 | 0 | 2 |
| WT | 1000 | 100 | 1 |
| WT | 1000 | 100 | 2 |
| EFF | 1000 | 100 | 1 |
| EFF | 1000 | 100 | 2 |
| P1104A | 1000 | 100 | 1 |
| P1104A | 1000 | 100 | 2 |

### Bioinformatics Processing <a name="part2"></a>

After raw data generation and demultiplexing, adapters were removed with cutadapt and then processed through two quantification procedures. Using STAR, we generated traditional alignments with very good results (80-90% unique across samples) followed by gene-level read counting with featureCounts and import into DESeq2. Separately, we analyzed the trimmed FASTQ's with kallisto to obtain pseudoalignments and imported them into DESeq2 with tximport.

The references used had the following modifications:

- hg38 contains an additional contig with the entire sequence of our synthetic construct except for the TYK2 gene itself
- gencode v41 transcripts contain an additional transcript corresponding to our luciferase sequence
- gencode v41 annotations contain additional gene, transcript, and exon entries corresponding to our luciferase sequence


    
![png](TYK2-Differential-Expression_files/TYK2-Differential-Expression_5_0.png)
    


### Differential Expression <a name="part3"></a>

Since it is more compatible with visual inspection and interpretation of results, the rest of this notebook uses the output of the "traditional" procedure using alignment and feature counting. To identify differentially expressed genes, we separate samples into groups by doxycycline concentration. Then, we can apply two possible modeling approaches to capture different kinds of effects. First, we can use a nested model structure which tests for significant cytokine effects in each genotype separately, while controlling for cytokine-independent differences in genotype groups. This model structure is shown below:


    
![png](TYK2-Differential-Expression_files/TYK2-Differential-Expression_7_0.png)
    


There are two types of comparisons we can make between sample groups:

1. First, we can test whether each gene is differentially expressed upon cytokine treatment in each genotype individually. We do this in a single model, and that model returns three sets of summary statistics that contain the fold changes upon cytokine treatment in each genotype, and the test against the null hypothesis that cytokine treated and untreated are equal.
2. Second, we can test whether the cytokine effect itself is different (a difference-of-differences) between genotypes. This secondary contrast returns two sets of summary statistics that compare the difference in cytokine-induced fold change between two genotypes, either WT vs P1104A or WT vs EFF.


    
![png](TYK2-Differential-Expression_files/TYK2-Differential-Expression_9_0.png)
    


### Primary Results: Cytokine Effects Within Each Genotype <a name="part4"></a>

There are a lot of different questions we can ask of these summary statistics, but the "primary" comparison is the cytokine response of each gene within each genotype. If we control the false discovery rate to 1%, we obtain the following numbers of significantly differentially expressed genes upon cytokine treatment in each condition:


<table class="dataframe">
<caption>A tibble: 6 × 3</caption>
<thead>
	<tr><th scope=col>doxy</th><th scope=col>DEG_count</th><th scope=col>explanation</th></tr>
	<tr><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 100</td><td> 16</td><td>Cytokine Effect in EFF   </td></tr>
	<tr><td>1000</td><td>  0</td><td>Cytokine Effect in EFF   </td></tr>
	<tr><td> 100</td><td> 11</td><td>Cytokine Effect in P1104A</td></tr>
	<tr><td>1000</td><td> 46</td><td>Cytokine Effect in P1104A</td></tr>
	<tr><td> 100</td><td> 11</td><td>Cytokine Effect in WT    </td></tr>
	<tr><td>1000</td><td>518</td><td>Cytokine Effect in WT    </td></tr>
</tbody>
</table>



#### Luciferase <a name="part4a"></a>

We can visualize the expression of the luciferase reporter in two ways. First, we can just plot the summary statistics, showing the effect size and confidence interval from each comparison. Second, we can extract and plot the normalized expression values directly.


    
![png](TYK2-Differential-Expression_files/TYK2-Differential-Expression_21_0.png)
    


The cytokine effect results in an approximately 10-16 fold change in luciferase in WT and P1104A, but not in EFF. There is no significant difference between genotypes on average, separate from the cytokine effect. While we could formally test whether there is a significant difference between the WT and P1104A effects, there is no real need to do so since the intervals are nearly identical.

We can also extract and plot the normalized expression of luciferase directly in each sample. These are the log2(counts per million):


    
![png](TYK2-Differential-Expression_files/TYK2-Differential-Expression_25_0.png)
    


#### Pathway Components <a name="part4b"></a>

Beyond luciferase, we are particularly interested in the relative expression of TYK2, JAK1, STAT1, STAT2, IFNAR1, IFNAR2, and IRF9. We extract the same normalized counts as above as `log2(counts per million)` and plot them as a grid from each doxycycline concentration. Note: these values are plotted regardless of whether the gene is expressed sufficiently for differential testing. For example, `IRF9` has far fewer than the required counts per sample for testing but is plotted below for comparison:


    
![png](TYK2-Differential-Expression_files/TYK2-Differential-Expression_27_0.png)
    



    
![png](TYK2-Differential-Expression_files/TYK2-Differential-Expression_27_1.png)
    


#### Volcano Plots <a name="part4c"></a>

In volcano plots, the unadjusted p-value is transformed (negative log) such that higher means more significant, and then plotted relative to the estimated log2 fold change. Luciferase clearly jumps out, but we also see some elevation of pathway components namely STAT1. EFF is entirely flat, while P1104A and WT have a similar structure:


    
![png](TYK2-Differential-Expression_files/TYK2-Differential-Expression_30_0.png)
    


### Secondary Results: Comparing Cytokine Effects Between Genotypes <a name="part5"></a>

In the above plots, we quantified the effects of cytokine in each genotype indidually. However, we are most interested in the genes whose cytokine response is different between genotypes. So we computed the secondary comparisons of EFF vs WT and P1104A vs WT, and tested for whether the cytokine effect itself in each genotype was equal. For both 100 nM and 1000 nM doxycycline, there were zero significantly different genes in the P1104A vs WT comparison, while there were 6 and 70 for EFF respectively.

#### Heatmaps <a name="part5a"></a>

The easiest way to view all of these genes together is with a heatmap. These plots show the cytokine ("primary") log2(Fold Change) for each gene in each genotype, but subset to only show the genes with significant differences between either EFF vs WT or P1104A vs WT ("secondary").


    
![png](TYK2-Differential-Expression_files/TYK2-Differential-Expression_34_0.png)
    



    
![png](TYK2-Differential-Expression_files/TYK2-Differential-Expression_35_0.png)
    

