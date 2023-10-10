### TYK2 Replicate Subsampling

For a full-length TYK2 DMS dataset, re-sampling replicates is difficult due to the computational requirements of re-fitting the mixed model for each subsample in sufficient throughputt to enable inference. While we may do this regardless (or do so for an individual position or subset of positions), a much faster approach is to simply remove the last sample ("D") from each replicate group and regenerate the summary statistics as before. Having done so, we can directly compare the results for each variant:


    
![png](TYK2-Rep-Subsample_files/TYK2-Rep-Subsample_5_0.png)
    


More interestingly, we can compare the Z-statistics to see what the differences in the final significant effect set look like:


    
![png](TYK2-Rep-Subsample_files/TYK2-Rep-Subsample_7_0.png)
    



    
![png](TYK2-Rep-Subsample_files/TYK2-Rep-Subsample_8_0.png)
    


We can count up the number of variants that fall into each of these categories (significant in both, neither, or only one) in each of the comparisons. These values, shown as counts and percentages of total variants tested, are in the below tables:


    
    
    |condition           |normalization | Both| N=3 Only| N=4 Only| Neither|
    |:-------------------|:-------------|----:|--------:|--------:|-------:|
    |IFNalpha1           |IFNbeta100    |   43|       16|       57|   23498|
    |IFNalpha10          |IFNbeta100    | 1882|       49|      573|   21110|
    |IFNalpha100         |IFNbeta100    | 2257|       43|      345|   20969|
    |IFNalphaWithDrug100 |IFNbeta100    |   25|       13|       18|   23558|
    |IFNalpha1           |None0         |    5|        9|       10|   23590|
    |IFNalpha10          |None0         | 1438|       71|      632|   21473|
    |IFNalpha100         |None0         | 2030|       58|      391|   21135|
    |IFNalphaWithDrug100 |None0         |    7|        4|        3|   23600|



    
    
    |condition           |normalization |      Both|  N=3 Only|  N=4 Only|  Neither|
    |:-------------------|:-------------|---------:|---------:|---------:|--------:|
    |IFNalpha1           |IFNbeta100    | 0.1820954| 0.0677564| 0.2413822| 99.50877|
    |IFNalpha10          |IFNbeta100    | 7.9698484| 0.2075040| 2.4265266| 89.39612|
    |IFNalpha100         |IFNbeta100    | 9.5578894| 0.1820954| 1.4609977| 88.79902|
    |IFNalphaWithDrug100 |IFNbeta100    | 0.1058694| 0.0550521| 0.0762260| 99.76285|
    |IFNalpha1           |None0         | 0.0211739| 0.0381130| 0.0423478| 99.89837|
    |IFNalpha10          |None0         | 6.0896079| 0.3006691| 2.6763784| 90.93334|
    |IFNalpha100         |None0         | 8.5965952| 0.2456170| 1.6557974| 89.50199|
    |IFNalphaWithDrug100 |None0         | 0.0296434| 0.0169391| 0.0127043| 99.94071|

