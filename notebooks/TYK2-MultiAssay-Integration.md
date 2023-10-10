### TYK2 Multi-Assay Integration

For the TYK2 variant pool, we now have data across many variants from two independent experimental methods: IFN-alpha DMS and VAMP-seq. Due to VAMP-seq dropout, we have three chunks (3, 10, and 14) with meaningful overlap between the two experiments. We also have multiple IFN-alpha concentrations in the DMS only, and multiple activity bins in VAMP-seq only.

1. [Primary LoF Signal Comparison](#part1)
2. [Contrast and Pairwise Correlations](#part2)
3. [Dimensional Projections](#part3)
4. [Biochemical Aggregation](#part4)
5. [Comparison with Midpoints](#part5)

#### Primary LoF Signal Comparison <a name="part1"></a>

We know from each assay individually that the clearest observed variant effect patterns are the strong LoF under either high IFN-alpha treatment (IFNalpha100) in DMS, and also that the slopes in conditions C and D were the strongest in VAMP-seq. We can extract just these summary statistics for the overlapping chunks, and plot them. Both sets of plots below are identical except the scale is zoomed in on the lower set:


    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_4_0.png)
    



    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_4_1.png)
    


We can combine the chunks and break out the stop effects as side densities, to avoid overloading the color scheme:


    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_6_0.png)
    



    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_6_1.png)
    


This is consistent with what we observe visually in the heatmaps - the VAMP-seq signal (particularly in chunks 10 and 14) is much broader than the DMS signal. There are also several variants that strongly deviate between the two assays, with the strongest such cases being those that are neutral in DMS but highly negative in VAMP-seq. Let's pull out that extreme variant, which is I684S:


    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_8_0.png)
    


There are also several variants which are signifiacant at a 1% FDR in both assays, but switch direction. Let's extract and profile those next:


    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_11_0.png)
    


#### Contrast and Pairwise Correlations <a name="part2"></a>




    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_15_0.png)
    



<table class="dataframe">
<caption>A tibble: 10 × 9</caption>
<thead>
	<tr><th scope=col>pos</th><th scope=col>chunk</th><th scope=col>aa</th><th scope=col>estimate</th><th scope=col>std.error</th><th scope=col>condition</th><th scope=col>statistic</th><th scope=col>p.value</th><th scope=col>p.adj</th></tr>
	<tr><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>638</td><td>10</td><td>*</td><td>1.433331</td><td>0.1212419</td><td>VAMP - DMS</td><td>11.822082</td><td>0.000000e+00</td><td>0.000000e+00</td></tr>
	<tr><td>694</td><td>10</td><td>K</td><td>3.077613</td><td>0.4067985</td><td>VAMP - DMS</td><td> 7.565448</td><td>3.863576e-14</td><td>1.802573e-11</td></tr>
	<tr><td>660</td><td>10</td><td>E</td><td>2.842491</td><td>0.4156396</td><td>VAMP - DMS</td><td> 6.838836</td><td>7.983836e-12</td><td>1.676206e-09</td></tr>
	<tr><td>696</td><td>10</td><td>L</td><td>3.149479</td><td>0.4601242</td><td>VAMP - DMS</td><td> 6.844846</td><td>7.655876e-12</td><td>1.676206e-09</td></tr>
	<tr><td>696</td><td>10</td><td>H</td><td>2.062465</td><td>0.3058917</td><td>VAMP - DMS</td><td> 6.742469</td><td>1.557177e-11</td><td>2.842863e-09</td></tr>
	<tr><td>652</td><td>10</td><td>E</td><td>1.748647</td><td>0.3053662</td><td>VAMP - DMS</td><td> 5.726395</td><td>1.025873e-08</td><td>7.977113e-07</td></tr>
	<tr><td>660</td><td>10</td><td>Q</td><td>2.204317</td><td>0.3868810</td><td>VAMP - DMS</td><td> 5.697663</td><td>1.214607e-08</td><td>8.947603e-07</td></tr>
	<tr><td>950</td><td>14</td><td>E</td><td>2.989064</td><td>0.5415613</td><td>VAMP - DMS</td><td> 5.519346</td><td>3.402631e-08</td><td>1.984396e-06</td></tr>
	<tr><td>696</td><td>10</td><td>G</td><td>2.291883</td><td>0.4238188</td><td>VAMP - DMS</td><td> 5.407695</td><td>6.384108e-08</td><td>3.309490e-06</td></tr>
	<tr><td>660</td><td>10</td><td>D</td><td>1.535333</td><td>0.3237177</td><td>VAMP - DMS</td><td> 4.742815</td><td>2.107685e-06</td><td>5.754518e-05</td></tr>
</tbody>
</table>




    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_16_1.png)
    



    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_18_0.png)
    


#### Dimensional Projections <a name="part3"></a>

Two very simple projections we can apply are t-SNE (non-linear) and PCA (traditional and linear), shown below in that order and with stops shown in red in all plots:


```R
mat <- as.matrix(dms_vamp[,-1:-3])
tsneobj <- Rtsne(mat)
projections <- as_tibble(tsneobj$Y) %>% mutate(pos = dms_vamp$pos, aa = dms_vamp$aa)
```


    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_21_0.png)
    



    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_22_0.png)
    


In this particular PCA, the first two components explain ~90% of the variance, and the first four get up to ~99%. We can show those four components across the heatmap grid as intensities:


    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_24_0.png)
    



    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_24_1.png)
    



    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_24_2.png)
    



    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_24_3.png)
    


We can also invert this view, and instead of considering each variant as a vector across all our conditions and assays, we can instead consider each sample as a vector across all tested variants. Doing so puts each sample into a reduced space with an informative pattern:


    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_26_0.png)
    


#### Biochemical Aggregation <a name="part4"></a>

Chemically, the structural range of encoded amino acids suggests several ways we can reduce the complexity of our data by combining groups of residues based on biochemical properties. Below are three previously used residue groupings, the first two applied to MC4R and the third applied to ADRB2:


    
    
    |aa            |group               |
    |:-------------|:-------------------|
    |N,Q           |Amide               |
    |V,I,L,M,F,Y,W |AromaticHydrophobic |
    |C             |Cysteine            |
    |D,E           |Negative            |
    |K,R,H         |Positive            |
    |P             |Proline             |
    |G,A,S,T       |SmallSerThr         |
    |*             |Stop                |



    
    
    |aa      |group            |
    |:-------|:----------------|
    |N,Q     |Amide            |
    |F,Y,W,H |Aromatic         |
    |C       |Cysteine         |
    |V,I     |Hydrophobic1     |
    |L,M     |Hydrophobic2     |
    |D,E     |Negative         |
    |K,R     |Positive         |
    |P       |Proline          |
    |S,T     |Serine/Threonine |
    |G,A     |Small            |
    |*       |Stop             |



    
    
    |aa      |group         |
    |:-------|:-------------|
    |N,Q     |Amide         |
    |F,W,Y   |Aromatic      |
    |I,L,V,M |Hydrophobic   |
    |D,E     |Negative      |
    |C,S,T   |Nucleophillic |
    |R,H,K   |Positive      |
    |P       |Proline       |
    |A,G     |Small         |
    |*       |Stop          |


Joining these and pivoting highlights where they diverge:


    
    
    |aa |set1                |set2             |set3          |
    |:--|:-------------------|:----------------|:-------------|
    |G  |SmallSerThr         |Small            |Small         |
    |A  |SmallSerThr         |Small            |Small         |
    |S  |SmallSerThr         |Serine/Threonine |Nucleophillic |
    |T  |SmallSerThr         |Serine/Threonine |Nucleophillic |
    |C  |Cysteine            |Cysteine         |Nucleophillic |
    |P  |Proline             |Proline          |Proline       |
    |V  |AromaticHydrophobic |Hydrophobic1     |Hydrophobic   |
    |I  |AromaticHydrophobic |Hydrophobic1     |Hydrophobic   |
    |L  |AromaticHydrophobic |Hydrophobic2     |Hydrophobic   |
    |M  |AromaticHydrophobic |Hydrophobic2     |Hydrophobic   |
    |F  |AromaticHydrophobic |Aromatic         |Aromatic      |
    |Y  |AromaticHydrophobic |Aromatic         |Aromatic      |
    |W  |AromaticHydrophobic |Aromatic         |Aromatic      |
    |N  |Amide               |Amide            |Amide         |
    |Q  |Amide               |Amide            |Amide         |
    |D  |Negative            |Negative         |Negative      |
    |E  |Negative            |Negative         |Negative      |
    |K  |Positive            |Positive         |Positive      |
    |R  |Positive            |Positive         |Positive      |
    |H  |Positive            |Aromatic         |Positive      |
    |*  |Stop                |Stop             |Stop          |


To combine summary statistics at each position, we can compute the average effect for each group. For some groups (like Proline) this is identical to the original. We can visualize these in all the usual ways:


    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_36_0.png)
    


As always, these make a little more sense as heatmaps:


    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_38_0.png)
    


The main LoF effects remain what we have observed in the past (stops, prolines, etc). The main interesting GoF effects we observed in the past were in the inhibitor-treated samples, and while these positions and effects seemed very specific, we can see if there are any significant GoF effects in this analysis:


    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_41_0.png)
    


This standout position is 671, which is strongly expected given our previous observations with these data. Since this number is not too large, we can extract all significant effects in this framework and compare. The main initial takeaways are:

  - position 671: resistance induced by almost anything
  - position 690: resistance induced by proline (consistent with initial analysis)
  - position 758: resistance induced by negatively charged residues


    
    
    |agg  |  pos|group               |condition                   |  estimate|     p.adj|
    |:----|----:|:-------------------|:---------------------------|---------:|---------:|
    |set2 |  434|Hydrophobic2        |IFNalphaWithDrug100 - None0 | 0.8887010| 0.0071012|
    |set2 |  544|Aromatic            |IFNalphaWithDrug100 - None0 | 0.3931423| 0.0040784|
    |set3 |  544|Aromatic            |IFNalphaWithDrug100 - None0 | 0.4413305| 0.0026723|
    |set1 |  550|Positive            |IFNalphaWithDrug100 - None0 | 0.4494357| 0.0095595|
    |set1 |  644|AromaticHydrophobic |IFNalphaWithDrug100 - None0 | 0.4047989| 0.0035017|
    |set1 |  671|Amide               |IFNalphaWithDrug100 - None0 | 0.7728190| 0.0000630|
    |set1 |  671|AromaticHydrophobic |IFNalphaWithDrug100 - None0 | 0.6030824| 0.0000001|
    |set1 |  671|Negative            |IFNalphaWithDrug100 - None0 | 0.7687693| 0.0001947|
    |set1 |  671|Positive            |IFNalphaWithDrug100 - None0 | 0.6238538| 0.0003423|
    |set2 |  671|Amide               |IFNalphaWithDrug100 - None0 | 0.7728190| 0.0000835|
    |set2 |  671|Aromatic            |IFNalphaWithDrug100 - None0 | 0.6575425| 0.0000275|
    |set2 |  671|Hydrophobic2        |IFNalphaWithDrug100 - None0 | 0.8157731| 0.0000835|
    |set2 |  671|Negative            |IFNalphaWithDrug100 - None0 | 0.7687693| 0.0002008|
    |set2 |  671|Positive            |IFNalphaWithDrug100 - None0 | 0.7381155| 0.0008986|
    |set3 |  671|Amide               |IFNalphaWithDrug100 - None0 | 0.7728190| 0.0000708|
    |set3 |  671|Aromatic            |IFNalphaWithDrug100 - None0 | 0.7449465| 0.0000426|
    |set3 |  671|Hydrophobic         |IFNalphaWithDrug100 - None0 | 0.4966844| 0.0004074|
    |set3 |  671|Negative            |IFNalphaWithDrug100 - None0 | 0.7687693| 0.0002190|
    |set3 |  671|Positive            |IFNalphaWithDrug100 - None0 | 0.6238538| 0.0003850|
    |set2 |  687|Positive            |IFNalphaWithDrug100 - None0 | 0.6205051| 0.0106421|
    |set1 |  690|Proline             |IFNalphaWithDrug100 - None0 | 1.2785861| 0.0000370|
    |set2 |  690|Proline             |IFNalphaWithDrug100 - None0 | 1.2785861| 0.0000509|
    |set3 |  690|Proline             |IFNalphaWithDrug100 - None0 | 1.2785861| 0.0000416|
    |set1 |  728|AromaticHydrophobic |IFNalphaWithDrug100 - None0 | 0.4052645| 0.0059830|
    |set1 |  758|Negative            |IFNalphaWithDrug100 - None0 | 1.4330368| 0.0000312|
    |set2 |  758|Negative            |IFNalphaWithDrug100 - None0 | 1.4330368| 0.0000428|
    |set3 |  758|Negative            |IFNalphaWithDrug100 - None0 | 1.4330368| 0.0000351|
    |set1 |  934|Cysteine            |IFNalphaWithDrug100 - None0 | 1.1701179| 0.0046257|
    |set2 |  934|Cysteine            |IFNalphaWithDrug100 - None0 | 1.1701179| 0.0063604|
    |set2 |  950|Small               |IFNalphaWithDrug100 - None0 | 0.6990701| 0.0072861|
    |set3 |  950|Small               |IFNalphaWithDrug100 - None0 | 0.6990701| 0.0059613|
    |set1 | 1028|AromaticHydrophobic |IFNalphaWithDrug100 - None0 | 0.3621301| 0.0082273|
    |set3 | 1038|Aromatic            |IFNalphaWithDrug100 - None0 | 0.6256203| 0.0073756|
    |set1 | 1141|AromaticHydrophobic |IFNalphaWithDrug100 - None0 | 0.3573828| 0.0070193|
    |set1 | 1186|Positive            |IFNalphaWithDrug100 - None0 | 0.7050422| 0.0018178|
    |set3 | 1186|Positive            |IFNalphaWithDrug100 - None0 | 0.7050422| 0.0020450|


#### Comparison with Midpoints <a name="part5"></a>

For chunk 10, we have the complete midpoint estimates with standard errors and Z-statistics, as well as the DMS profile across several assays. We can contrast these explicitly, since inspection (see the midpoint notebook) indicates the Z-statistics are similarly scaled across both assays, including with similar spike-in patterns.


    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_47_0.png)
    



    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_47_1.png)
    



    
![png](TYK2-MultiAssay-Integration_files/TYK2-MultiAssay-Integration_50_0.png)
    

