### TYK2 Partial Loss-of-Function

#### Full Length Partial LoF

We regenerated summary statistics using the stop aggregation model (where all stops are considered within each chunk as their own "position", and compared to all WT). Plotting all variants along the x-axis is very visually cluttered, but scatter plots of the Z-statsistics between mutant effects relative to WT and relative to aggregated Stop provide a way to quickly see the trends:


    
![png](TYK2-Partial-LoF_files/TYK2-Partial-LoF_5_0.png)
    


Counting up the number of partial LoF (red points) in each condition, we find:


    
    
    |condition                | partial LoF|
    |:------------------------|-----------:|
    |IFNalpha10 - None0       |           4|
    |IFNalpha100 - None0      |          37|
    |IFNalpha10 - IFNbeta100  |          16|
    |IFNalpha100 - IFNbeta100 |          74|


In contrast, we have the following total counts for LoF variants per comparison:


    
    
    |condition                        | total LoF|
    |:--------------------------------|---------:|
    |IFNalpha1 - IFNbeta100           |        53|
    |IFNalpha1 - None0                |         9|
    |IFNalpha10 - IFNbeta100          |      1665|
    |IFNalpha10 - None0               |      1439|
    |IFNalpha100 - IFNbeta100         |      1686|
    |IFNalpha100 - None0              |      1569|
    |IFNalphaWithDrug100 - IFNbeta100 |         2|
    |IFNbeta100 - None0               |        33|
    |None0 - IFNbeta100               |         8|


Let's look at the overlap between normalization conditions:


    
![png](TYK2-Partial-LoF_files/TYK2-Partial-LoF_11_0.png)
    


We can examine the relationship between partial LoF called in either of the IFNalpha100 conditions directly:


    
![png](TYK2-Partial-LoF_files/TYK2-Partial-LoF_14_0.png)
    


There are thresholding effects - all variants which are significantly partial LoF in one of the above two conditions are also in the correct directions in the other condition, just not necessairly strongly enough to be significant.

#### Comparison to mini-DMS

We performed the same analysis on the first mini-DMS, so we can make the same comparison we did above for the two conditions. We extract all variants that were significantly partial LoF in one or the other (or both), and then plot their Z-scores for WT and Stop effects:


    
![png](TYK2-Partial-LoF_files/TYK2-Partial-LoF_18_0.png)
    


The above plot shows that all partial LoF variants identified in the full length assay (in chunk 10) are directionally consistent with the original mini-DMS, though not always significant in the mini-DMS. We can check the reverse too - that significant partial LoF variants in the mini-DMS are directionally consistent in the new data, and it looks like they indeed are (with a few exceptions):


    
![png](TYK2-Partial-LoF_files/TYK2-Partial-LoF_20_0.png)
    

