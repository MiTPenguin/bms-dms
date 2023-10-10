## OCNT-DMSLIB-1 DMS Assay Run #8: Protocol Comparisons

| ID | Protocol | Condition | Concentration (M) | Replicates
| --- | --- | -------- | --- | ----------- |
| 1 | A | None | 0 | 4 |
| 2 | A | IFN-alpha (100 U/mL) | 0 | 4 |
| 3 | A | BMS-986202 + IFN-alpha (100U/mL) | 2e-8 | 4 |
| 4 | B | None | 0 | 4 |
| 5 | B | IFN-alpha (100 U/mL) | 0 | 4 |
| 6 | B | BMS-986202 + IFN-alpha (100U/mL) | 2e-8 | 4 |

1. [Barcode Sequencing Distributions](#part1)
2. [Inference and Stop Codon Effects](#part2)
3. [Visualizations](#part3)
4. [Protocol Comparisons](#part4)
5. [Potentiation Comparisons](#part5)

### Barcode Sequencing Distributions <a name="part1"></a>


    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_4_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_4_1.png)
    


To get a sense of positional distribution, we can show the same data as lineplots across the length of TYK2. Below is an example using sample `1A` only; the remaining plots can be found [here](./coverage-plots):


    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_7_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_8_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_9_0.png)
    


### Inference and Stop Codon Effects <a name="part2"></a>


    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_16_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_16_1.png)
    


We can also examine directly the number of variants across all comparisons that are significant at a 1% FDR:


    
    
    |protocol  |drug                   |  conc|direction | Non-Significant| Significant (FDR < 0.01)|
    |:---------|:----------------------|-----:|:---------|---------------:|------------------------:|
    |protocol1 |IFNalpha100            | 0e+00|GoF       |            8749|                        0|
    |protocol1 |IFNalpha100            | 0e+00|LoF       |           12596|                     2375|
    |protocol2 |IFNalpha100            | 0e+00|GoF       |            8822|                        3|
    |protocol2 |IFNalpha100            | 0e+00|LoF       |           12721|                     2172|
    |protocol1 |IFNalpha100+BMS-986202 | 2e-08|GoF       |            8155|                       28|
    |protocol1 |IFNalpha100+BMS-986202 | 2e-08|LoF       |           13181|                     2356|
    |protocol2 |IFNalpha100+BMS-986202 | 2e-08|GoF       |            8511|                       29|
    |protocol2 |IFNalpha100+BMS-986202 | 2e-08|LoF       |           13243|                     1935|


### Visualizations <a name="part3"></a>


    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_21_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_21_1.png)
    



    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_23_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_23_1.png)
    


### Protocol Comparisons <a name="part4"></a>

The simplest comparison is to plot the same summary statistic from both protocols, across all variants. We can use either the normalized or unnormalized results, and can compare the log2 fold change, z-statistic, and the standard errors:


    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_25_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_25_1.png)
    



    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_26_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_26_1.png)
    



    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_27_0.png)
    


This is a little hard to see, so below is the same plot but zoomed in a bit and with just the boxplots:


    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_29_0.png)
    


Using the normalized summary statistics, we can count and compare the number of significant effects in either direction:


    
    
    |protocol  |drug                   |  LoF|    NS| GoF|
    |:---------|:----------------------|----:|-----:|---:|
    |protocol1 |IFNalpha100            | 2375| 21345|   0|
    |protocol2 |IFNalpha100            | 2172| 21543|   3|
    |protocol1 |IFNalpha100+BMS-986202 | 2356| 21336|  28|
    |protocol2 |IFNalpha100+BMS-986202 | 1935| 21754|  29|


For the three sets where both protocols have some significant variants (IFNalpha100 LoF, IFNalpha100+BMS-986202 LoF, and IFNalpha100+BMS-986202 GoF) we can count how many variants are significant in one, the other, or both:


    
    
    |group     | IFNalpha100 LoF| IFNalpha100+BMS-986202 LoF| IFNalpha100+BMS-986202 GoF|
    |:---------|---------------:|--------------------------:|--------------------------:|
    |Both      |            1989|                       1683|                         16|
    |protocol1 |             386|                        673|                         12|
    |protocol2 |             183|                        252|                         13|


For the IFNalpha100+BMS-986202 GoF variants in particular, we can examine the underlying summary statistics more closely for those that are significant in either or both protocols:


    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_36_0.png)
    


### Potentiation Comparisons <a name="part5"></a>

In this dataset, we care mainly about evaluating potentiation effects. To do so, we don't normalize to "None", but instead subtract the IFN-alpha condition from the IFN-alpha + BMS condition within each protocol separately. Comparing the results will let us evaluate whether we identify similar potentiators in both datasets. The simple heatmap below shows which variants are significant in both or one protocol: 


    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_39_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run8_files/OCNT-DMSLIB-1-assay-run8_40_0.png)
    

