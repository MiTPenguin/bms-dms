#### TYK2 Library Diversity QC

We have two chunks (15 and 16), each with two types of samples (plasmid and gDNA) and two replicates of each group (A and  B) for 2\*2\*2 = 8 samples total. We demultiplexed, counted reads per unique barcode sequence, and mapped those barcode counts back to their parent oligos. Considering the set of barcodes succesfully associated to an oligo, we have the following counts of unique barcodes per sample:


    
![png](TYK2-Lib-Div-QC_files/TYK2-Lib-Div-QC_4_0.png)
    



    
![png](TYK2-Lib-Div-QC_files/TYK2-Lib-Div-QC_4_1.png)
    


Next, we can plot the number of reads for each "unique barcode" detected. We see that in addition to there being fewer mapped barcodes per variant in chunk 16, there are also far more barcodes with very low counts in chunk 16 in gDNA (but not in plasmid):


    
![png](TYK2-Lib-Div-QC_files/TYK2-Lib-Div-QC_6_0.png)
    


If we check the total number of unique barcodes in each sample, we find similar numbers in both chunk 15 and chunk 16 gDNA. The difference (based on the above plot) is that most of these barcodes are very low abundance in chunk 16, while they have the expected distribution in chunk 15:


    
    
    |sample   | number of mapped barcodes detected|
    |:--------|----------------------------------:|
    |15A_gdna |                              70342|
    |15A_maxi |                            1350060|
    |15B_gdna |                              90503|
    |15B_maxi |                            1424653|
    |16A_gdna |                              74227|
    |16A_maxi |                             496077|
    |16B_gdna |                              53241|
    |16B_maxi |                             384350|



    
    
    |sample   | median read count per mapped barcode|
    |:--------|------------------------------------:|
    |15A_gdna |                                   69|
    |15A_maxi |                                   53|
    |15B_gdna |                                   61|
    |15B_maxi |                                   56|
    |16A_gdna |                                    1|
    |16A_maxi |                                   20|
    |16B_gdna |                                    1|
    |16B_maxi |                                   16|

