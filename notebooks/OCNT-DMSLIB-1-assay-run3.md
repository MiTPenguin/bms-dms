## OCNT-DMSLIB-1 DMS Assay Run #3: TYK2 IFN-alpha Gradient

This data represents the first sequencing run of what will be two total runs. It contains 24 samples in the following experimental design, with replicates in each sample group labeled `A` through `D`:

| ID | Condition | Dosage | Replicates
| --- | ----------- | --- | ----------- |
| 1 | IFN-beta | 100 | 4 |
| 2 | None | 0 | 4 |
| 3 | IFN-alpha | 1 | 4 |
| 4 | IFN-alpha | 10 | 4 |
| 5 | IFN-alpha | 100 | 4 |
| 6 | IFN-alpha + Drug | 100 + Drug | 4 |

Let's check some effects we expect to observe (mainly QC parameters and stop codons), then examine the global distribution of mutant vs WT effects across conditions.

1. [Barcode Sequencing Distributions](#part1)
2. [Inference and Stop Codon Effects](#part2)
3. [Visualizations](#part3)
4. [Gain-of-Function](#part4)
5. [Spike-Ins](#part5)
6. [Positionwise Aggregation](#part6)


### Barcode Sequencing Distributions <a name="part1"></a>

In most samples, we just miss the desired 30 barcodes per sample (though the medians are generally pretty close).There is also substantial dropout in chunk 16, with only about 31% of those variants making it through. Otherwise, most amino acid level variants have a median 24 unique barcode per residue.


    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_5_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_5_1.png)
    


To get a sense of positional distribution, we can show the same data as lineplots across the length of TYK2. Below is an example using sample `1A` only; the remaining plots can be found [here](./coverage-plots):


    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_9_0.png)
    


Finally, we can count up the number of unique barcodes per chunk per sample, and observe that the relative proportions are highly consistent across replicate groups:


    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_11_0.png)
    


Or, if we add all chunks together and just get the total number of unique barcodes per sample, we generally hover around 650,000 or so:


    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_13_0.png)
    


### Inference and Stop Codon Effects <a name="part2"></a>

We'll normalize each condition to either IFN-beta or None using the summary statistics, and show results for both. Our usual starting point is asking whether we observe the expected, significantly negative stop codon effects. The easiest way to show this, while allowing for some chunk specificity, is to plot the density of Z-statistics for all non-stop variants in all chunks in each condition, but separate the stop variant distributions by chunk. 


    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_19_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_19_1.png)
    


We can also examine directly the number of variants across all comparisons that are significant at a 1% FDR:


    
    
    |condition                        | Non-Significant| Significant (FDR < 0.01)|
    |:--------------------------------|---------------:|------------------------:|
    |IFNbeta100 - None0               |           23677|                       46|
    |IFNalpha1 - None0                |           23706|                       17|
    |IFNalpha10 - None0               |           21524|                     2199|
    |IFNalpha100 - None0              |           21076|                     2647|
    |IFNalphaWithDrug100 - None0      |           23713|                       10|
    |None0 - IFNbeta100               |           23677|                       46|
    |IFNalpha1 - IFNbeta100           |           23620|                      103|
    |IFNalpha10 - IFNbeta100          |           21129|                     2594|
    |IFNalpha100 - IFNbeta100         |           20904|                     2819|
    |IFNalphaWithDrug100 - IFNbeta100 |           23680|                       43|


### Visualizations <a name="part4"></a>

Heatmaps and volcano plots are very efficient ways to take in the landscape of these data. In these heatmaps, the color scale is constant across all facets and the plotted values include the log2FoldChange and the Z-Statistic. The scale is unrestricted - in both plot sets, it was allowed to take on the most extreme values which existed in the data. Note how some regions have clearly a noisier baseline, and so stick out in the log2FoldChange heatmap, but which is correctly accounted for and thus flat in the Z-Statistic heatmap.


    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_24_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_24_1.png)
    



    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_25_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_25_1.png)
    



    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_27_0.png)
    


### Gain-of-Function <a name="part4"></a>

We clearly can identify Loss-of-Function variants very well, but can we detect any significant Gain-of-Function variants? Below, we normalize by subtracting `None` from `IFNalpha10` and sorting by the most significant _positive_ effects. We use `IFNalpha10` reasoning that a lighter stimulus might allow a larger dynamic range for positive, expression-increasing effects:


    
    
    |  pos|clone |chunk |aa | log2FoldChange| std.error| dispersion|version |condition          | statistic|   p.value|     p.adj|
    |----:|:-----|:-----|:--|--------------:|---------:|----------:|:-------|:------------------|---------:|---------:|---------:|
    |  220|c8    |4     |S  |      1.4234249| 0.3456665|  0.8625597|v2.1.0  |IFNalpha10 - None0 |  4.117914| 0.0000382| 0.0011271|
    |  544|c8    |8     |W  |      0.8458308| 0.2185625|  0.8692214|v2.1.0  |IFNalpha10 - None0 |  3.869972| 0.0001088| 0.0028774|
    |   96|c8    |2     |W  |      1.2677983| 0.3493724|  0.8751519|v2.1.0  |IFNalpha10 - None0 |  3.628788| 0.0002848| 0.0067691|
    |  262|c8    |4     |I  |      1.2390584| 0.3476599|  0.8625785|v2.1.0  |IFNalpha10 - None0 |  3.563996| 0.0003653| 0.0084223|
    | 1142|c8    |17    |Q  |      1.0683177| 0.3001314|  0.8645981|v2.1.0  |IFNalpha10 - None0 |  3.559500| 0.0003716| 0.0085459|
    |  678|c8    |10    |F  |      0.2880072| 0.0813095|  0.9019933|v2.1.0  |IFNalpha10 - None0 |  3.542110| 0.0003969| 0.0090539|
    |  370|c8    |6     |K  |      1.4636668| 0.4159884|  0.8564798|v2.1.0  |IFNalpha10 - None0 |  3.518528| 0.0004339| 0.0098034|
    |  908|c8    |13    |W  |      1.4129279| 0.4063805|  0.8789082|v2.1.0  |IFNalpha10 - None0 |  3.476859| 0.0005073| 0.0112437|
    |  998|c8    |15    |C  |      1.6703148| 0.4825362|  0.8769124|v2.1.0  |IFNalpha10 - None0 |  3.461533| 0.0005371| 0.0118276|
    |  966|c8    |14    |Q  |      0.9093408| 0.2641394|  0.8713134|v2.1.0  |IFNalpha10 - None0 |  3.442655| 0.0005760| 0.0125900|


And the same, normalized to `IFNbeta100`:


    
    
    |  pos|clone |chunk |aa | log2FoldChange| std.error| dispersion|version |condition               | statistic|   p.value|     p.adj|
    |----:|:-----|:-----|:--|--------------:|---------:|----------:|:-------|:-----------------------|---------:|---------:|---------:|
    |  678|c8    |10    |F  |      0.5511683| 0.0802311|  0.9019933|v2.1.0  |IFNalpha10 - IFNbeta100 |  6.869756| 0.0000000| 0.0000000|
    |  644|c8    |10    |G  |      1.0318508| 0.2562107|  0.9031483|v2.1.0  |IFNalpha10 - IFNbeta100 |  4.027352| 0.0000564| 0.0015976|
    |  994|c8    |15    |Q  |      1.0903653| 0.2915067|  0.8772489|v2.1.0  |IFNalpha10 - IFNbeta100 |  3.740447| 0.0001837| 0.0045847|
    |  865|c8    |13    |E  |      1.0837661| 0.2953856|  0.8828775|v2.1.0  |IFNalpha10 - IFNbeta100 |  3.668987| 0.0002435| 0.0059020|
    | 1176|c8    |17    |T  |      1.1696930| 0.3333135|  0.8651627|v2.1.0  |IFNalpha10 - IFNbeta100 |  3.509288| 0.0004493| 0.0101138|
    |  129|c8    |2     |E  |      0.9372831| 0.2747876|  0.8733427|v2.1.0  |IFNalpha10 - IFNbeta100 |  3.410936| 0.0006474| 0.0139634|
    |  687|c8    |10    |F  |      0.6529143| 0.1961317|  0.9070758|v2.1.0  |IFNalpha10 - IFNbeta100 |  3.328959| 0.0008717| 0.0181336|
    |  286|c8    |5     |S  |      0.9007315| 0.2769944|  0.8819481|v2.1.0  |IFNalpha10 - IFNbeta100 |  3.251804| 0.0011468| 0.0230234|
    |   23|c8    |1     |F  |      0.9441966| 0.2911887|  0.8671667|v2.1.0  |IFNalpha10 - IFNbeta100 |  3.242559| 0.0011846| 0.0236589|
    |  229|c8    |4     |W  |      1.0475061| 0.3231032|  0.8648528|v2.1.0  |IFNalpha10 - IFNbeta100 |  3.242017| 0.0011869| 0.0236984|


Examining the heatmaps, there seem to be some very specific gain of function mutations that jump out of the combined treatment with IFN-alpha and Drug. Let's subtract None and IFN-beta from IFN-alpha with Drug and see:


    
    
    | pos|clone |chunk |aa | log2FoldChange| std.error| dispersion|version |condition                   | statistic|  p.value|     p.adj|
    |---:|:-----|:-----|:--|--------------:|---------:|----------:|:-------|:---------------------------|---------:|--------:|---------:|
    | 671|c8    |10    |Y  |      1.8762327| 0.2528983|  0.9051413|v2.1.0  |IFNalphaWithDrug100 - None0 |  7.418921| 0.00e+00| 0.0000000|
    | 671|c8    |10    |Q  |      1.6281940| 0.2250351|  0.9051413|v2.1.0  |IFNalphaWithDrug100 - None0 |  7.235291| 0.00e+00| 0.0000000|
    | 671|c8    |10    |E  |      1.6032344| 0.2870867|  0.9051413|v2.1.0  |IFNalphaWithDrug100 - None0 |  5.584496| 0.00e+00| 0.0000015|
    | 690|c8    |10    |P  |      1.2785861| 0.2790964|  0.9042129|v2.1.0  |IFNalphaWithDrug100 - None0 |  4.581163| 4.60e-06| 0.0001722|
    | 803|c8    |12rc  |W  |      1.5363432| 0.3629702|  0.8927898|v2.1.0  |IFNalphaWithDrug100 - None0 |  4.232698| 2.31e-05| 0.0007182|
    | 758|c8    |11    |H  |      1.6713900| 0.3956623|  0.9110504|v2.1.0  |IFNalphaWithDrug100 - None0 |  4.224284| 2.40e-05| 0.0007422|
    | 671|c8    |10    |L  |      1.1474319| 0.2842464|  0.9051413|v2.1.0  |IFNalphaWithDrug100 - None0 |  4.036751| 5.42e-05| 0.0015405|
    | 687|c8    |10    |R  |      0.9494082| 0.2379089|  0.9070758|v2.1.0  |IFNalphaWithDrug100 - None0 |  3.990638| 6.59e-05| 0.0018320|
    | 758|c8    |11    |D  |      1.6511111| 0.4172947|  0.9110504|v2.1.0  |IFNalphaWithDrug100 - None0 |  3.956703| 7.60e-05| 0.0020846|
    | 671|c8    |10    |K  |      1.1590278| 0.3095617|  0.9051413|v2.1.0  |IFNalphaWithDrug100 - None0 |  3.744093| 1.81e-04| 0.0045244|



    
    
    |  pos|clone |chunk |aa | log2FoldChange| std.error| dispersion|version |condition                        | statistic|  p.value|     p.adj|
    |----:|:-----|:-----|:--|--------------:|---------:|----------:|:-------|:--------------------------------|---------:|--------:|---------:|
    |  671|c8    |10    |Q  |       1.802742| 0.2221113|  0.9051413|v2.1.0  |IFNalphaWithDrug100 - IFNbeta100 |  8.116392| 0.00e+00| 0.0000000|
    |  758|c8    |11    |D  |       2.754826| 0.4134491|  0.9110504|v2.1.0  |IFNalphaWithDrug100 - IFNbeta100 |  6.663036| 0.00e+00| 0.0000000|
    |  671|c8    |10    |Y  |       1.558015| 0.2478449|  0.9051413|v2.1.0  |IFNalphaWithDrug100 - IFNbeta100 |  6.286251| 0.00e+00| 0.0000000|
    |  671|c8    |10    |L  |       1.417474| 0.2769616|  0.9051413|v2.1.0  |IFNalphaWithDrug100 - IFNbeta100 |  5.117946| 3.00e-07| 0.0000153|
    |  671|c8    |10    |E  |       1.383624| 0.2823585|  0.9051413|v2.1.0  |IFNalphaWithDrug100 - IFNbeta100 |  4.900239| 1.00e-06| 0.0000421|
    |  952|c8    |14    |G  |       1.352558| 0.2893725|  0.8710758|v2.1.0  |IFNalphaWithDrug100 - IFNbeta100 |  4.674107| 3.00e-06| 0.0001153|
    | 1102|c8    |16    |N  |       3.048775| 0.6733919|  0.9133294|v2.1.0  |IFNalphaWithDrug100 - IFNbeta100 |  4.527490| 6.00e-06| 0.0002157|
    |  780|c8    |12rc  |Q  |       1.594299| 0.3522103|  0.8940226|v2.1.0  |IFNalphaWithDrug100 - IFNbeta100 |  4.526556| 6.00e-06| 0.0002166|
    |  942|c8    |14    |E  |       2.759433| 0.6155687|  0.8693487|v2.1.0  |IFNalphaWithDrug100 - IFNbeta100 |  4.482738| 7.40e-06| 0.0002605|
    |  687|c8    |10    |R  |       1.025473| 0.2346518|  0.9070758|v2.1.0  |IFNalphaWithDrug100 - IFNbeta100 |  4.370191| 1.24e-05| 0.0004148|


Since position 671 jumps out, let's extract the summary statistics for that position in all conditions. Significant effects (FDR < 1%) are highlighted in magenta:


    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_35_0.png)
    


There are some _very_ strong (compared to what we usually see) gain-of-function effects in chunk 10 associated with drug treatment (there is visually no trace of these in any condition other than IFN-alpha 100 + Drug). We obviously detect far, far more LoF variation but we do appear to have sufficient power to begin to curate a set of putative GoF variants in addition to LoF.

### Spike-Ins <a name="part5"></a>

Similar to Alanine671, we can make the same plots for all our spike-in variants:


    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_39_0.png)
    


### Positionwise Aggregation <a name="part6"></a>

For painting on a structure, we need a single number per position instead of 20 numbers (and 20 standard errors on those numbers). We could subset particular amino acids to collapse, but the simplest starting point is to just collapse everything except for stop codons - that is, take the average across all mutation effects in all non-stop amino acids. This dampens the effect of most mutations, since in most cases this will average inactive variants at the same position with active variants, but returns tracks where each TYK2 position is assigned a single number:


    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_43_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_43_1.png)
    


To prevent the dampening, we could instead take the average variant effect at each position for only significant (FDR < 1%) variants, and assign zero to all other non-significant variants. Now, many positions shrink to zero and significant positions are very spiky:


    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_46_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_46_1.png)
    


Another approach to aggregation is to take the maximum significant effect at each position, where the "maximum effect" means the largest absolute z-statistic per position:


    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_48_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run3_files/OCNT-DMSLIB-1-assay-run3_48_1.png)
    

