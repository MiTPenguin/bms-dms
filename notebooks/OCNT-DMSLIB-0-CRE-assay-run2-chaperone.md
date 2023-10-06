## OCNT-DMSLIB-0-CRE DMS Assay Run #2: MC4R, Gs, Ipsen (DMS11)

| ID | Condition | Concentration | Replicates
| --- | ----------- | --- | ----------- |
| 1 | Forskolin | 2.5e-5 | 4 |
| 2 | None | 0 | 4 |
| 3 | aMSH | 1e-6 | 4 |
| 4 | Forskolin+Ipsen | 2.5e-5 | 4 |
| 5 | None+Ipsen | 4e-10 | 4 |
| 6 | aMSH+Ipsen | 1e-6 | 4 |

This report generates similar results to the full-length assay (DMS5) with aMSH and THIQ, but focused on identifying variants whose effect is mitigated or rescued by Ipsen treatment.

1. [Barcode Sequencing Distributions](#part1)
2. [Statistical Models of Chaperone Rescue](#part2)
3. [Visualizations](#part3)
4. [The "None" Condition](#part4)

### Barcode Sequencing Distributions <a name="part1"></a>


    
![png](OCNT-DMSLIB-0-CRE-assay-run2-chaperone_files/OCNT-DMSLIB-0-CRE-assay-run2-chaperone_5_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run2-chaperone_files/OCNT-DMSLIB-0-CRE-assay-run2-chaperone_5_1.png)
    


### Statistical Models of Chaperone Rescue <a name="part2"></a>

We could extract many different comparisons from these data. A good starting point is to normalize each of None and aMSH to Forskolin within each Ipsen group (so normalize aMSH-Ipsen to Forskolin-Ipsen, and so on). This corresponds to our previous normalization comparisons (though the aMSH concentration is very different), and we can plot the distribution of stop effects per chunk as before:


    
![png](OCNT-DMSLIB-0-CRE-assay-run2-chaperone_files/OCNT-DMSLIB-0-CRE-assay-run2-chaperone_9_0.png)
    


We can directly plot the Z-statistics of aMSH with Ipsen (normalized to Forskolin) versus the Z-statistic of aMSH without Ipsen (normalized to Forskolin without Ipsen) from the above densities, separating out stops and non-stops:


    
![png](OCNT-DMSLIB-0-CRE-assay-run2-chaperone_files/OCNT-DMSLIB-0-CRE-assay-run2-chaperone_11_0.png)
    


This is easier to interpret if we color each point by the direction of the difference between Ipsen and No Ipsen effects. In this plot, the gold points are those where the mutant is higher activity than WT in aMSH+Ipsen, but lower than WT in aMSH alone:


    
![png](OCNT-DMSLIB-0-CRE-assay-run2-chaperone_files/OCNT-DMSLIB-0-CRE-assay-run2-chaperone_13_0.png)
    


Stops are generally a bit _more_ negative in Ipsen compared to NoIpsen. However, we will be mostly interested in identifying variants which are on the other end of the distribution. We can summarize the number that are significant in either direction:


    
    
    |group    |direction              | NonSignificant| Significant|
    |:--------|:----------------------|--------------:|-----------:|
    |Stop     |More Active With Ipsen |             47|           1|
    |Stop     |Less Active With Ipsen |            256|          27|
    |Non-Stop |More Active With Ipsen |           2805|         149|
    |Non-Stop |Less Active With Ipsen |           2836|         512|


### Visualizations <a name="part3"></a>

To better understand what these comparisons are doing, we can generate both the heatmaps for each aMSH condition separately, and then show the heatmap of the difference which corresponds to the comparison shown above. First, we show the aMSH and aMSH + Ipsen heatmaps (both normalized to their respective Forskolin):


    
![png](OCNT-DMSLIB-0-CRE-assay-run2-chaperone_files/OCNT-DMSLIB-0-CRE-assay-run2-chaperone_19_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run2-chaperone_files/OCNT-DMSLIB-0-CRE-assay-run2-chaperone_19_1.png)
    


Then, we can compute the heatmap representing the difference between the two aMSH heatmaps above:


    
![png](OCNT-DMSLIB-0-CRE-assay-run2-chaperone_files/OCNT-DMSLIB-0-CRE-assay-run2-chaperone_21_0.png)
    


And equivalently, represent the results as a volcano plot using the log2FoldChange and the unadjusted p-values:


    
![png](OCNT-DMSLIB-0-CRE-assay-run2-chaperone_files/OCNT-DMSLIB-0-CRE-assay-run2-chaperone_23_0.png)
    


The resulting significant, positive variant set excluding stops contains 149 variants, which we can summarize:


    
![png](OCNT-DMSLIB-0-CRE-assay-run2-chaperone_files/OCNT-DMSLIB-0-CRE-assay-run2-chaperone_25_0.png)
    



    
![png](OCNT-DMSLIB-0-CRE-assay-run2-chaperone_files/OCNT-DMSLIB-0-CRE-assay-run2-chaperone_25_1.png)
    



    
    
    | pos|  n|
    |---:|--:|
    |  69|  7|
    | 161|  6|
    | 162|  5|
    | 280|  5|
    | 284|  5|
    | 135|  4|



    
    
    |aa |  n|
    |:--|--:|
    |S  | 19|
    |L  | 15|
    |R  | 15|
    |P  | 12|
    |A  | 10|
    |G  | 10|


### The "None" Condition <a name="part4"></a>

The None condition is a bit strange, since for the main question described in the previous section we do not really need it. However, it is worth mentioning that there was a very unusual and striking pattern in the None data when normalized to Forskolin:


    
![png](OCNT-DMSLIB-0-CRE-assay-run2-chaperone_files/OCNT-DMSLIB-0-CRE-assay-run2-chaperone_28_0.png)
    


There appears to be very little basal activity, but that is only because the overwhelming signal from the Ipsen-treated sample is pushing the scale up. We can regenerate these with a more constrained scale:


    
![png](OCNT-DMSLIB-0-CRE-assay-run2-chaperone_files/OCNT-DMSLIB-0-CRE-assay-run2-chaperone_30_0.png)
    


Even with this adjustment, there is a very clear variant set which is strongly activated (their activity becomes significantly higher than WT) upon Ipsen treatment, in the absence of any aMSH. Separately, basal activity in the lower panel is somewhat visible but much weaker than in the Gs/DMS5 assay.
